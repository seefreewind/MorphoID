#!/usr/bin/env python3
"""Freeze the Vannan cohort and run non-expression Phase 0D metadata audits.

This script performs dependency/alias, donor-fold, within-disease, spatial
dependence, cell-composition coupling, and target-coverage audits. It never
reads expression values and never trains a gene-expression model.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import chi2_contingency
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, normalized_mutual_info_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / "results/phase0/phase0c_v"
VR = ROOT / "results/phase0/phase0c_vr"
OUT = ROOT / "results/phase0/phase0d"
FIG = ROOT / "figures/phase0d"
CONFIG = ROOT / "configs"
SEED = 20260927
GRID_UM = 150
LINEAGES = ["epithelial", "endothelial", "fibroblast", "smooth_muscle_pericyte",
            "macrophage", "monocyte", "neutrophil", "CD4_T", "CD8_T", "B",
            "plasma", "mast", "other_immune", "other_stromal", "other"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def cramers_v(a: pd.Series, b: pd.Series) -> float:
    tab = pd.crosstab(a.astype(str), b.astype(str))
    if min(tab.shape) < 2 or tab.to_numpy().sum() == 0:
        return float("nan")
    chi2 = chi2_contingency(tab, correction=False)[0]
    n = tab.to_numpy().sum()
    r, k = tab.shape
    phi2 = chi2 / n
    phi2corr = max(0.0, phi2 - ((k - 1) * (r - 1)) / max(n - 1, 1))
    rcorr = r - ((r - 1) ** 2) / max(n - 1, 1)
    kcorr = k - ((k - 1) ** 2) / max(n - 1, 1)
    denom = min(kcorr - 1, rcorr - 1)
    return float(math.sqrt(phi2corr / denom)) if denom > 0 else float("nan")


def assoc_class(v: float, exact: bool, informative: bool = True) -> str:
    if not informative:
        return "NOT_INFORMATIVE_NO_VARIATION"
    if exact:
        return "EXACT_ALIAS"
    if not np.isfinite(v):
        return "UNAVAILABLE"
    if v >= 0.70:
        return "HIGH_ASSOCIATION"
    if v >= 0.40:
        return "MODERATE"
    return "LOW"


def deterministic_directions(a: pd.Series, b: pd.Series) -> tuple[bool, bool]:
    frame = pd.DataFrame({"a": a.astype(str), "b": b.astype(str)}).dropna()
    return (bool(frame.groupby("a").b.nunique().le(1).all()),
            bool(frame.groupby("b").a.nunique().le(1).all()))


def make_fold_assignment(donors: pd.DataFrame, seed: int = SEED, k: int = 5) -> pd.DataFrame:
    """Greedy deterministic donor assignment, disease-stratified then balanced.

    Donors are first distributed round-robin within disease strata. Same-disease
    donor swaps then minimize imbalance in section counts, microregion counts,
    TMA counts, and run counts without changing the disease allocation.
    """
    rng = np.random.default_rng(seed)
    d = donors.sort_values("donor_id").reset_index(drop=True).copy()
    strata = []
    for disease, g in d.groupby("disease", sort=True):
        ids = g.donor_id.to_numpy().copy()
        rng.shuffle(ids)
        strata.append((disease, ids))
    assignment: dict[str, int] = {}
    for _, ids in strata:
        # Rotate start position between strata to avoid deterministic alignment.
        offset = len(assignment) % k
        for i, donor in enumerate(ids):
            assignment[str(donor)] = (i + offset) % k

    tmas = sorted({f"TMA{x}" for vals in d.TMA_set for x in vals.split(";") if x})
    runs = sorted({str(x) for vals in d.run_set for x in vals.split(";") if x})
    features = []
    for row in d.itertuples(index=False):
        values = [float(row.n_sections), float(row.n_microregions)]
        values += [float(f"TMA{x}" in row.TMA_set.split(";")) for x in range(1, 5)]
        values += [float(run in row.run_set.split(";")) for run in runs]
        features.append(values)
    X = np.asarray(features, dtype=float)
    target = X.sum(axis=0) / k
    scales = np.maximum(target, 1.0)
    donor_idx = {str(x): i for i, x in enumerate(d.donor_id)}

    def cost(mapping: dict[str, int]) -> float:
        fold_x = np.zeros((k, X.shape[1]), dtype=float)
        for donor, fold in mapping.items():
            fold_x[fold] += X[donor_idx[donor]]
        imbalance = np.square((fold_x - target) / scales).sum()
        counts = np.bincount(list(mapping.values()), minlength=k)
        imbalance += 0.2 * np.square(counts - len(d) / k).sum()
        return float(imbalance)

    # Deterministic best-improving swaps within disease strata.
    for _ in range(100):
        improved = False
        current = cost(assignment)
        for disease, _ids in strata:
            group_ids = d.loc[d.disease.eq(disease), "donor_id"].astype(str).tolist()
            pairs = list(itertools.combinations(group_ids, 2))
            rng.shuffle(pairs)
            for a, b in pairs:
                if assignment[a] == assignment[b]:
                    continue
                trial = assignment.copy()
                trial[a], trial[b] = trial[b], trial[a]
                new_cost = cost(trial)
                if new_cost + 1e-12 < current:
                    assignment = trial
                    current = new_cost
                    improved = True
        if not improved:
            break
    out = d[["donor_id", "disease", "TMA_set", "n_sections", "n_microregions"]].copy()
    out["fold"] = out.donor_id.map(lambda x: f"fold_{assignment[str(x)] + 1}")
    out["grouping_unit"] = "donor_id"
    out["seed"] = seed
    return out.sort_values(["fold", "donor_id"]).reset_index(drop=True)


def build_association_audit(sections: pd.DataFrame, donors: pd.DataFrame) -> pd.DataFrame:
    records = []
    pairs = [
        ("donor", "TMA", donors, "donor_id", "TMA_set"),
        ("donor", "run", donors, "donor_id", "run_set"),
        ("donor", "disease", donors, "donor_id", "disease"),
        ("disease", "TMA", sections, "disease_status", "TMA"),
        ("disease", "run", sections, "disease_status", "run"),
        ("disease", "collection_center", sections, "disease_status", "collection_center"),
        ("disease", "affected_status", sections, "disease_status", "affected_status"),
        ("TMA", "run", sections, "TMA", "run"),
        ("TMA", "collection_center", sections, "TMA", "collection_center"),
        ("TMA", "scanner", sections, "TMA", "he_scanner"),
        ("donor", "section_count", donors, "donor_id", "n_sections"),
    ]
    for left_name, right_name, frame, left_col, right_col in pairs:
        temp = frame[[left_col, right_col]].dropna().copy()
        if right_name == "TMA" and left_name == "donor":
            temp[right_col] = temp[right_col].astype(str).str.split(";")
            temp = temp.explode(right_col)
        elif right_name == "run" and left_name == "donor":
            temp[right_col] = temp[right_col].astype(str).str.split(";")
            temp = temp.explode(right_col)
        v = cramers_v(temp[left_col], temp[right_col])
        l2r, r2l = deterministic_directions(temp[left_col], temp[right_col])
        exact = l2r and r2l
        informative = temp[left_col].nunique() > 1 and temp[right_col].nunique() > 1
        if left_name == "donor":
            category = "DONOR_UNIT_DETERMINISTIC" if l2r else "DONOR_CROSS_CLASS_STRUCTURE"
            if informative and exact:
                category = "EXACT_ALIAS"
        else:
            category = assoc_class(v, exact, informative)
        records.append({
            "pair": f"{left_name}__by__{right_name}", "unit": "unique_donor" if frame is donors else "section",
            "n_rows": len(temp), "n_left_levels": temp[left_col].nunique(),
            "n_right_levels": temp[right_col].nunique(), "cramers_v_bias_corrected": v,
            "normalized_mutual_information": float(normalized_mutual_info_score(
                temp[left_col].astype(str), temp[right_col].astype(str), average_method="arithmetic")),
            "left_to_right_deterministic": l2r, "right_to_left_deterministic": r2l,
            "exact_deterministic_alias": exact,
            "association_class": category,
            "notes": ("Donor identity is the grouping unit; one row per donor makes Cramer's V uninformative. Directional determinism/NMI are retained." if left_name == "donor"
                     else "Disease is donor-constant; outer folds group all sections by donor." if "disease" in {left_name, right_name}
                     else "Bias-corrected Cramer's V on frozen admitted sections/donors."),
        })
    records.append({
        "pair": "TMA__by__acquisition_date", "unit": "section", "n_rows": len(sections),
        "n_left_levels": sections.TMA.nunique(), "n_right_levels": 0,
        "cramers_v_bias_corrected": np.nan, "normalized_mutual_information": np.nan,
        "left_to_right_deterministic": False, "right_to_left_deterministic": False,
        "exact_deterministic_alias": False, "association_class": "UNAVAILABLE_NOT_IN_SOURCE_METADATA",
        "notes": "Per-sample acquisition date was not found; GEO upload timestamps are not acquisition dates.",
    })
    return pd.DataFrame(records)


def make_graph(sections: pd.DataFrame, donors: pd.DataFrame, edges: pd.DataFrame) -> None:
    graph = nx.DiGraph()
    nodes_by_type: dict[str, list[str]] = {}
    for edge in edges.itertuples(index=False):
        src = f"{edge.source_type}:{edge.source_id}"
        dst = f"{edge.target_type}:{edge.target_id}"
        graph.add_edge(src, dst, edge_type=edge.edge_type)
        nodes_by_type.setdefault(edge.source_type, []).append(src)
        nodes_by_type.setdefault(edge.target_type, []).append(dst)
    nodes_by_type = {k: sorted(set(v)) for k, v in nodes_by_type.items()}
    layers = ["donor", "core", "section", "registration_status", "primary_inclusion",
              "TMA", "run", "registered_HE", "disease", "affected_status", "institution", "scanner", "panel", "instrument",
              "analysis_version", "technical_batch"]
    x = {name: float(i) for i, name in enumerate(layers)}
    pos = {}
    for layer in layers:
        nodes = nodes_by_type.get(layer, [])
        if not nodes:
            continue
        ys = np.linspace(1.0, -1.0, len(nodes)) if len(nodes) > 1 else [0.0]
        pos.update({node: (x[layer], float(y)) for node, y in zip(nodes, ys)})
    type_by_node = {node: typ for typ, nodes in nodes_by_type.items() for node in nodes}
    colors = []
    for node in graph.nodes:
        typ = type_by_node.get(node, "")
        if typ == "donor":
            disease = donors.set_index("donor_id").loc[node.split(":", 1)[1], "disease"]
            colors.append("#C94C4C" if disease == "pulmonary_fibrosis" else "#4776A8")
        elif typ == "section":
            disease = sections.set_index("sample_id").loc[node.split(":", 1)[1], "disease_status"]
            colors.append("#C94C4C" if str(disease).lower() in {"pf", "disease", "pulmonary_fibrosis"} else "#4776A8")
        elif typ == "donor":
            colors.append("#666666")
        else:
            colors.append("#E7A83E")
    fig, ax = plt.subplots(figsize=(18, 9))
    nx.draw_networkx_edges(graph, pos, ax=ax, alpha=0.13, width=0.5, arrows=False)
    nx.draw_networkx_nodes(graph, pos, node_color=colors, node_size=50, ax=ax, linewidths=0.15,
                           edgecolors="#333333")
    labels = {n: n.split(":", 1)[-1] for n in graph.nodes}
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=4.4, ax=ax)
    for layer in layers:
        ax.text(x[layer], 1.12, layer.replace("_", " "), ha="center", va="bottom", fontsize=8,
                fontweight="bold", rotation=0)
    ax.set_title("Frozen Vannan cohort dependency graph · outer split group = donor", fontsize=13)
    ax.set_axis_off()
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "vannan_dependency_graph.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_freeze_and_static_audits() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    CONFIG.mkdir(parents=True, exist_ok=True)
    score = pd.read_csv(VR / "registration_evidence_scorecard.tsv", sep="\t")
    cross = pd.read_csv(VR / "author_sample_run_crosswalk.tsv", sep="\t")
    source = pd.read_csv(V / "vannan_sample_crosswalk.tsv", sep="\t")
    source = source[source.requested_28_19_primary_subset.astype(str).eq("YES")].copy()
    primary = cross.drop(columns=["TMA"], errors="ignore").merge(
        score[["sample_id", "final_registration_status"]], on="sample_id", validate="one_to_one")
    primary = primary.merge(source[["sample_id", "disease_status", "affected_status", "TMA", "run",
                                   "collection_center", "he_scanner", "instrument_id",
                                   "xenium_analysis_version", "core_id", "replicate_id", "panel_version"]],
                            on="sample_id", validate="one_to_one")
    admitted = primary[primary.final_registration_status.eq("REGISTRATION_PASS_WITH_QUALIFICATION")].copy()
    excluded = primary[~primary.sample_id.isin(admitted.sample_id)].copy()
    if len(primary) != 28 or len(admitted) != 26 or admitted.donor_id.nunique() != 19:
        raise RuntimeError(f"Frozen Phase0CVR cohort mismatch: sections={len(admitted)}, donors={admitted.donor_id.nunique()}")
    cohort = primary[["sample_id", "donor_id", "TMA", "run", "disease_status", "affected_status",
                      "final_registration_status"]].copy()
    cohort["primary_inclusion"] = np.where(
        cohort.final_registration_status.eq("REGISTRATION_PASS_WITH_QUALIFICATION"), "YES", "NO")
    cohort["exclusion_reason"] = np.where(cohort.primary_inclusion.eq("YES"), "",
                                          "FIXED_COORDINATE_SUPPORT_FAILURE")
    cohort = cohort.rename(columns={"final_registration_status": "registration_status"})
    cohort = cohort.sort_values(["donor_id", "sample_id"])
    cohort_path = CONFIG / "frozen_vannan_primary_cohort.tsv"
    cohort.to_csv(cohort_path, sep="\t", index=False)
    (CONFIG / "frozen_vannan_primary_cohort.sha256").write_text(
        f"{sha256(cohort_path)}  {cohort_path.name}\n", encoding="utf-8")

    admitted = admitted.sort_values(["donor_id", "sample_id"])
    micro = pd.read_csv(V / "microregion_feasibility.tsv", sep="\t", usecols=[
        "sample_id", "donor_id", "disease_status", "affected_status", "region_x", "region_y", "grid_um",
        *[f"n_{x}" for x in LINEAGES], "total_cells", "n_adjacent_nonempty_regions_8", "has_nonempty_neighbor"])
    # The Phase 0C-V feasibility source carries parallel 100/150/200-um grids.
    # Every Phase 0D primary audit is frozen to 150 um; never pool scales.
    micro = micro[micro.sample_id.isin(admitted.sample_id) & micro.grid_um.eq(GRID_UM)].copy()
    if "Control" in micro.disease_status.astype(str).unique():
        micro["disease"] = micro.disease_status.astype(str).str.lower().replace({"control": "control", "disease": "pulmonary_fibrosis"})
    else:
        micro["disease"] = micro.disease_status
    section_regions = micro.groupby("sample_id", as_index=False).size().rename(columns={"size": "n_microregions"})
    section_regions["sample_id"] = section_regions.sample_id.astype(str)
    admitted = admitted.merge(section_regions, on="sample_id", how="left", validate="one_to_one")

    # Donor summary and fixed 5-fold candidate allocation.
    donor_rows = []
    for donor, g in admitted.groupby("donor_id", sort=True):
        diseases = sorted(g.disease_status.astype(str).str.lower().unique())
        disease = "control" if all(x in {"control", "unaffected"} for x in diseases) else "pulmonary_fibrosis"
        tma_set = ";".join(sorted(set("TMA" + g.TMA.astype(str).str.replace("TMA", "", regex=False))))
        run_set = ";".join(sorted(set("Run" + g.run.astype(str).str.replace("Run", "", regex=False))))
        donor_rows.append({"donor_id": donor, "disease": disease, "TMA_set": tma_set,
                           "run_set": run_set, "n_sections": len(g),
                           "n_microregions": int(g.n_microregions.sum())})
    donors = pd.DataFrame(donor_rows)
    folds = make_fold_assignment(donors)
    folds = folds.rename(columns={"TMA_set": "TMA_set", "disease": "disease"})
    folds.to_csv(CONFIG / "phase1_donor_folds.tsv", sep="\t", index=False)

    # Compare the frozen custom allocation with ordinary GroupKFold on the
    # same donor rows. Disease remains an evaluation attribute, never a split unit.
    gkf_rows = []
    splitter = GroupKFold(n_splits=5)
    for i, (_train, test) in enumerate(splitter.split(
            np.zeros((len(donors), 1)), donors.disease, groups=donors.donor_id), 1):
        held = donors.iloc[test]
        sec = admitted[admitted.donor_id.isin(held.donor_id)]
        gkf_rows.append({"fold": f"fold_{i}", "heldout_donors": len(held),
                         "heldout_control_donors": int(held.disease.eq("control").sum()),
                         "heldout_PF_donors": int(held.disease.eq("pulmonary_fibrosis").sum()),
                         "heldout_sections": len(sec), "heldout_microregions": int(sec.n_microregions.sum()),
                         "TMA_sections": ";".join(f"{k}:{v}" for k,v in sec.TMA.astype(str).value_counts().sort_index().items()),
                         "run_sections": ";".join(f"{k}:{v}" for k,v in sec.run.astype(str).value_counts().sort_index().items())})
    pd.DataFrame(gkf_rows).to_csv(OUT / "ordinary_groupkfold_donor_balance.tsv", sep="\t", index=False)

    fold_metrics = []
    for fold, g in folds.groupby("fold", sort=True):
        donor_ids = set(g.donor_id)
        sec = admitted[admitted.donor_id.isin(donor_ids)]
        fold_metrics.append({"fold": fold, "heldout_donors": len(g),
                             "heldout_control_donors": int(g.disease.eq("control").sum()),
                             "heldout_PF_donors": int(g.disease.eq("pulmonary_fibrosis").sum()),
                             "heldout_sections": len(sec), "heldout_microregions": int(sec.n_microregions.sum()),
                             "TMA_sections": ";".join(f"{k}:{v}" for k,v in sec.TMA.astype(str).value_counts().sort_index().items()),
                             "run_sections": ";".join(f"{k}:{v}" for k,v in sec.run.astype(str).value_counts().sort_index().items())})
    pd.DataFrame(fold_metrics).to_csv(OUT / "donor_fold_balance.tsv", sep="\t", index=False)

    # Dependency edges retain missing technical fields explicitly.
    edges = []
    def edge(src_type, src_id, dst_type, dst_id, evidence):
        edges.append({"source_type": src_type, "source_id": str(src_id), "target_type": dst_type,
                      "target_id": str(dst_id), "edge_type": f"{src_type}_to_{dst_type}", "evidence": evidence})
    # Keep all 28 audited sections in the dependency graph so that the two
    # permanent fixed-coordinate exclusions remain visible as provenance.
    for r in primary.itertuples(index=False):
        edge("donor", r.donor_id, "core", r.core_id, "GEO/author crosswalk")
        edge("core", r.core_id, "section", r.sample_id, "frozen sample crosswalk")
        edge("section", r.sample_id, "registration_status", r.final_registration_status,
             "Phase 0C-VR frozen scorecard")
        edge("section", r.sample_id, "primary_inclusion",
             "YES" if r.final_registration_status == "REGISTRATION_PASS_WITH_QUALIFICATION" else "NO",
             "frozen cohort inclusion rule")
        edge("section", r.sample_id, "TMA", f"TMA{str(r.TMA).replace('TMA','')}", "GEO sample metadata")
        edge("section", r.sample_id, "run", f"Run{str(r.run).replace('Run','')}", "GEO sample metadata + author code audit")
        edge("section", r.sample_id, "registered_HE", r.registered_HE_file, "GSM-specific registered H&E")
        edge("section", r.sample_id, "disease", "PF" if str(r.disease_status).lower() not in {"control", "unaffected"} else "Control", "GEO sample metadata")
        edge("section", r.sample_id, "affected_status", r.affected_status, "GEO sample metadata")
        edge("section", r.sample_id, "institution", r.collection_center, "GEO sample metadata")
        edge("section", r.sample_id, "scanner", r.he_scanner, "publication-level scanner; no per-image EXIF")
        edge("section", r.sample_id, "panel", r.panel_version, "author-corrected 343-gene panel")
        edge("section", r.sample_id, "instrument", r.instrument_id, "GEO sample metadata")
        edge("section", r.sample_id, "analysis_version", r.xenium_analysis_version, "GEO sample metadata")
        edge("section", r.sample_id, "technical_batch", f"{r.TMA}|{r.run}|{r.instrument_id}", "TMA/run/instrument proxy; no separate batch ID reported")
    edge_df = pd.DataFrame(edges)
    edge_df.to_csv(OUT / "dependency_edges.tsv", sep="\t", index=False)
    # Add unreported attributes as explicit audit records; do not invent their values.
    missing_tech = pd.DataFrame([
        {"variable": "acquisition_date", "status": "UNAVAILABLE_PER_SAMPLE", "basis": "GEO upload time is not image acquisition time"},
        {"variable": "scanner_per_image", "status": "UNAVAILABLE_PER_IMAGE", "basis": "scanner stated only at publication-level methods; TIFF EXIF does not establish per-image device"},
        {"variable": "technical_batch_id", "status": "UNREPORTED", "basis": "TMA/run/instrument retained as documented proxies"},
    ])
    missing_tech.to_csv(OUT / "technical_metadata_availability.tsv", sep="\t", index=False)
    make_graph(primary, donors, edge_df)

    # Donor/core structure flags.
    donor_detail = []
    for donor, g in admitted.groupby("donor_id", sort=True):
        statuses = set(g.affected_status.astype(str).str.lower())
        donor_detail.append({"donor_id": donor, "n_admitted_sections": len(g),
                             "section_ids": ";".join(sorted(g.sample_id.astype(str))),
                             "affected_levels": ";".join(sorted(statuses)),
                             "has_less_and_more_affected_pair": {"less_affected", "more_affected"}.issubset(statuses),
                             "TMA_set": ";".join(sorted(set("TMA" + g.TMA.astype(str).str.replace("TMA", "", regex=False)))),
                             "run_set": ";".join(sorted(set("Run" + g.run.astype(str).str.replace("Run", "", regex=False))))})
    pd.DataFrame(donor_detail).to_csv(OUT / "donor_core_dependency_summary.tsv", sep="\t", index=False)

    # Contingency/alias statistics. Use the frozen sections and donor-level summaries.
    sections = admitted.copy()
    sections["disease_status"] = sections.disease_status.astype(str).str.lower().replace({"control": "control"})
    sections["TMA"] = sections.TMA.astype(str).map(lambda x: x if x.startswith("TMA") else f"TMA{x}")
    sections["run"] = sections.run.astype(str).map(lambda x: x if x.startswith("Run") else f"Run{x}")
    assoc = build_association_audit(sections, donors)
    assoc.to_csv(OUT / "contingency_alias_audit.tsv", sep="\t", index=False)
    tables = OUT / "contingency_tables"
    tables.mkdir(exist_ok=True)
    pair_specs = [
        ("donor", "TMA", donors.assign(TMA=donors.TMA_set.str.split(";" )).explode("TMA"), "donor_id", "TMA"),
        ("donor", "run", donors.assign(run=donors.run_set.str.split(";")).explode("run"), "donor_id", "run"),
        ("donor", "disease", donors, "donor_id", "disease"),
        ("disease", "TMA", sections, "disease_status", "TMA"),
        ("disease", "run", sections, "disease_status", "run"),
        ("disease", "collection_center", sections, "disease_status", "collection_center"),
        ("disease", "affected_status", sections, "disease_status", "affected_status"),
        ("TMA", "run", sections, "TMA", "run"),
        ("TMA", "collection_center", sections, "TMA", "collection_center"),
        ("TMA", "scanner", sections, "TMA", "he_scanner"),
    ]
    for a_name, b_name, frame, a_col, b_col in pair_specs:
        tab = pd.crosstab(frame[a_col], frame[b_col])
        tab.to_csv(tables / f"{a_name}_by_{b_name}.tsv", sep="\t")
    donor_section_counts = donors[["donor_id", "n_sections"]].copy()
    donor_section_counts.to_csv(tables / "donor_by_section_count.tsv", sep="\t", index=False)
    # No per-sample acquisition date is available; record the requested pair
    # explicitly rather than substituting GEO submission dates.
    (tables / "TMA_by_acquisition_date.tsv").write_text(
        "availability\tUNAVAILABLE_PER_SAMPLE\n"
        "basis\tGEO upload timestamps are not image acquisition dates\n",
        encoding="utf-8")

    # PF-only and control-only feasibility: regions and lineage cell counts.
    feasibility = []
    count_cols = [f"n_{x}" for x in LINEAGES]
    for disease, g in micro.groupby("disease", sort=True):
        dg = g.groupby("donor_id").agg(n_regions=("sample_id", "size"), n_sections=("sample_id", "nunique"))
        n_donors = len(dg)
        state = "FEASIBLE" if n_donors >= 5 else "WEAK" if n_donors >= 3 else "INFEASIBLE"
        feasibility.append({"disease": disease, "independent_donors": n_donors,
                            "sections": g.sample_id.nunique(), "150um_regions": len(g),
                            "cells": int(g.total_cells.sum()), "median_cells_per_region": float(g.total_cells.median()),
                            "classification": state})
        for lineage in LINEAGES:
            eligible = g[g[f"n_{lineage}"] >= 20]
            donors_eligible = eligible.groupby("donor_id").size()
            feasibility.append({"disease": disease, "lineage": lineage,
                                "independent_donors": int((donors_eligible >= 10).sum()),
                                "sections": eligible.sample_id.nunique(), "150um_regions": len(eligible),
                                "cells": int(eligible[f"n_{lineage}"].sum()),
                                "median_cells_per_region": float(eligible[f"n_{lineage}"].median()) if len(eligible) else 0,
                                "classification": "CELL_FEASIBLE_10_REGIONS_PER_DONOR" if (donors_eligible >= 10).sum() >= 5 else "WEAK_CELL_SUPPORT"})
    pd.DataFrame(feasibility).to_csv(OUT / "within_disease_feasibility.tsv", sep="\t", index=False)

    # Microregion clustering and spatial dependence (counts only; no expression).
    region_rows = []
    per_section = []
    for sample, g in micro.groupby("sample_id", sort=True):
        coords = g[["region_x", "region_y"]].to_numpy(float) * GRID_UM
        tree = cKDTree(coords)
        d2, _ = tree.query(coords, k=2)
        nn = d2[:, 1]
        xx = np.log1p(g.total_cells.to_numpy(float))
        coords_int = g[["region_x", "region_y"]].to_numpy(int)
        lookup = {tuple(p): i for i, p in enumerate(coords_int)}
        num = 0.0
        W = 0
        for i, (gx, gy) in enumerate(coords_int):
            for dx, dy in ((-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)):
                j = lookup.get((gx+dx, gy+dy))
                if j is not None:
                    num += (xx[i]-xx.mean()) * (xx[j]-xx.mean())
                    W += 1
        denom = float(np.square(xx-xx.mean()).sum())
        moran = (len(g)/W) * num / denom if W and denom > 0 else np.nan
        rec = {"sample_id": sample, "donor_id": g.donor_id.iloc[0], "disease": g.disease.iloc[0],
               "n_regions": len(g), "regions_with_neighbor": int(g.has_nonempty_neighbor.sum()),
               "fraction_regions_with_8_neighbor": float(g.has_nonempty_neighbor.mean()),
               "median_nearest_nonempty_region_distance_um": float(np.median(nn)),
               "mean_nearest_nonempty_region_distance_um": float(nn.mean()),
               "moran_I_log1p_cell_count_8_neighbor": float(moran)}
        per_section.append(rec)
        region_rows.append(rec)
    pd.DataFrame(per_section).to_csv(OUT / "spatial_dependence_by_section.tsv", sep="\t", index=False)
    donor_stats = micro.groupby("donor_id").agg(n_regions=("sample_id", "size"), n_sections=("sample_id", "nunique"),
                                                mean_cells=("total_cells", "mean"),
                                                var_cells=("total_cells", "var")).reset_index()
    overall_mean = micro.total_cells.mean()
    within = sum(((g.total_cells - g.total_cells.mean())**2).sum() for _,g in micro.groupby("donor_id"))
    between = sum(len(g)*(g.total_cells.mean()-overall_mean)**2 for _,g in micro.groupby("donor_id"))
    msw = within / max(len(micro)-len(donor_stats), 1)
    msb = between / max(len(donor_stats)-1,1)
    n0 = (len(micro) - (donor_stats.n_regions.pow(2).sum()/len(micro))) / max(len(donor_stats)-1,1)
    icc = max(0.0, (msb-msw)/(msb+(n0-1)*msw)) if (msb+(n0-1)*msw)>0 else np.nan
    spatial_summary = pd.DataFrame([{
        "admitted_donors": len(donor_stats), "admitted_sections": micro.sample_id.nunique(),
        "nonempty_150um_regions": len(micro), "mean_regions_per_donor": donor_stats.n_regions.mean(),
        "median_regions_per_donor": donor_stats.n_regions.median(),
        "mean_sections_per_donor": donor_stats.n_sections.mean(),
        "fraction_regions_with_nonempty_8_neighbor": micro.has_nonempty_neighbor.mean(),
        "median_section_moran_I_log1p_cell_count": pd.Series([x["moran_I_log1p_cell_count_8_neighbor"] for x in per_section]).median(),
        "one_way_donor_ICC_total_cell_count_approx": icc,
        "independent_inference_unit": "donor_id",
        "primary_microregion_bootstrap_allowed": "NO",
    }])
    spatial_summary.to_csv(OUT / "microregion_independence_summary.tsv", sep="\t", index=False)

    # Composition-only, donor-held-out disease coupling. C is ground-truth
    # coarse composition; no expression or morphology features enter.
    donor_comp = micro.groupby(["donor_id", "disease"], as_index=False)[count_cols].sum()
    for col in count_cols:
        donor_comp["p_" + col[2:]] = donor_comp[col] / donor_comp[count_cols].sum(axis=1).clip(lower=1)
    comp_features = ["p_" + x for x in LINEAGES]
    donor_comp = donor_comp.merge(folds[["donor_id", "fold"]], on="donor_id", validate="one_to_one")
    y = donor_comp.disease.astype(str).to_numpy()
    Xc = donor_comp[comp_features].to_numpy(float)
    pred = np.empty(len(donor_comp), dtype=object)
    prob = np.zeros((len(donor_comp), 2), dtype=float)
    classes = np.array(["control", "pulmonary_fibrosis"])
    for fold in sorted(donor_comp.fold.unique()):
        te = donor_comp.fold.eq(fold).to_numpy()
        tr = ~te
        scale = StandardScaler().fit(Xc[tr])
        model = LogisticRegression(C=0.1, max_iter=1000, random_state=SEED)
        if len(np.unique(y[tr])) < 2:
            pred[te] = y[tr][0]
            prob[te, classes == y[tr][0]] = 1.0
        else:
            model.fit(scale.transform(Xc[tr]), y[tr])
            p = model.predict_proba(scale.transform(Xc[te]))
            pred[te] = model.classes_[np.argmax(p, axis=1)]
            for j, cls in enumerate(model.classes_):
                prob[np.ix_(te, classes == cls)] = p[:, [j]]
    observed_ba = balanced_accuracy_score(y, pred)
    rng = np.random.default_rng(SEED + 5)
    null = []
    for _ in range(500):
        yp = rng.permutation(y)
        pp = np.empty(len(yp), dtype=object)
        for fold in sorted(donor_comp.fold.unique()):
            te = donor_comp.fold.eq(fold).to_numpy(); tr = ~te
            scale = StandardScaler().fit(Xc[tr])
            model = LogisticRegression(C=0.1, max_iter=1000, random_state=SEED)
            if len(np.unique(yp[tr])) < 2:
                pp[te] = yp[tr][0]
            else:
                model.fit(scale.transform(Xc[tr]), yp[tr])
                pp[te] = model.predict(scale.transform(Xc[te]))
        null.append(balanced_accuracy_score(yp, pp))
    p_perm = (1 + sum(x >= observed_ba for x in null)) / (len(null)+1)
    boot_rng = np.random.default_rng(SEED + 6)
    boot = []
    for _ in range(2000):
        ix = boot_rng.integers(0, len(y), len(y))
        if len(set(y[ix])) < 2: continue
        boot.append(balanced_accuracy_score(y[ix], pred[ix]))
    comp_result = pd.DataFrame([{
        "target": "disease", "predictors": "ground-truth 15-lineage donor aggregate composition C only",
        "n_donors": len(y), "folds": 5, "balanced_accuracy": observed_ba,
        "macro_auroc": np.nan, "bootstrap_95ci_low": float(np.quantile(boot, .025)),
        "bootstrap_95ci_high": float(np.quantile(boot, .975)), "permutation_replicates": len(null),
        "permutation_baseline_mean": float(np.mean(null)), "permutation_baseline_95pct": float(np.quantile(null, .95)),
        "permutation_p": p_perm, "interpretation": "composition-disease coupling diagnostic; not expression prediction or causal attribution",
    }])
    comp_result.to_csv(OUT / "composition_disease_coupling.tsv", sep="\t", index=False)
    donor_comp.assign(oof_pred=pred, oof_p_control=prob[:,0], oof_p_PF=prob[:,1]).to_csv(
        OUT / "composition_disease_oof_by_donor.tsv", sep="\t", index=False)

    # Test whether coarse section-level composition also carries technical or
    # local-severity identity. The evaluation unit is a section summary, while
    # folds and bootstrap blocks remain donors. Donor identity itself is not
    # predictable in a donor-held-out split and is recorded as not estimable.
    section_comp = micro.groupby("sample_id", as_index=False)[count_cols].sum()
    section_comp = section_comp.merge(
        admitted[["sample_id", "donor_id", "TMA", "run", "disease_status", "affected_status"]],
        on="sample_id", validate="one_to_one")
    section_comp_features = []
    denom = section_comp[count_cols].sum(axis=1).clip(lower=1)
    for col in count_cols:
        name = "p_" + col[2:]
        section_comp[name] = section_comp[col] / denom
        section_comp_features.append(name)
    section_comp = section_comp.merge(folds[["donor_id", "fold"]], on="donor_id", validate="many_to_one")
    composition_shortcut_rows = []
    for target, source_col in [("TMA", "TMA"), ("run", "run"),
                               ("disease", "disease_status"), ("affected_status", "affected_status")]:
        frame = section_comp.copy()
        frame["label"] = frame[source_col].astype(str)
        classes_local = np.array(sorted(frame.label.unique()))
        y_local = frame.label.to_numpy()
        x_local = frame[section_comp_features].to_numpy(float)
        pred_local = np.empty(len(frame), dtype=object)
        prob_local = np.zeros((len(frame), len(classes_local)), dtype=float)
        unsupported_local = set()
        train_w = frame.donor_id.map(lambda d: 1.0 / int(frame.donor_id.eq(d).sum())).to_numpy(float)
        for fold_name in sorted(frame.fold.unique()):
            te = frame.fold.eq(fold_name).to_numpy()
            tr = ~te
            train_classes = set(y_local[tr]); unsupported_local.update(set(y_local[te]) - train_classes)
            if len(train_classes) < 2:
                only = sorted(train_classes)[0]
                pred_local[te] = only
                prob_local[np.ix_(te, np.flatnonzero(classes_local == only))] = 1.0
                continue
            scaler = StandardScaler().fit(x_local[tr])
            model = LogisticRegression(C=0.1, max_iter=1000, class_weight="balanced", random_state=SEED)
            model.fit(scaler.transform(x_local[tr]), y_local[tr], sample_weight=train_w[tr])
            pp = model.predict_proba(scaler.transform(x_local[te]))
            pred_local[te] = model.classes_[np.argmax(pp, axis=1)]
            for j, cls in enumerate(model.classes_):
                prob_local[np.ix_(te, np.flatnonzero(classes_local == cls))] = pp[:, [j]]
        donor_nsections = frame.groupby("donor_id").sample_id.nunique()
        score_w = frame.donor_id.map(lambda d: 1.0 / donor_nsections[d]).to_numpy(float)
        ba_local = balanced_accuracy_score(y_local, pred_local, sample_weight=score_w)
        try:
            if unsupported_local:
                raise ValueError("class absent from one or more training folds")
            auc_local = roc_auc_score(y_local, prob_local, labels=classes_local,
                                      multi_class="ovr", average="macro", sample_weight=score_w)
        except Exception:
            auc_local = np.nan
        rng_local = np.random.default_rng(SEED + 31 + len(composition_shortcut_rows))
        boot_local = []
        donor_ids_local = frame.donor_id.unique()
        for _ in range(2000):
            chosen = rng_local.choice(donor_ids_local, len(donor_ids_local), replace=True)
            blocks = [np.flatnonzero(frame.donor_id.to_numpy() == d) for d in chosen]
            ix = np.concatenate(blocks)
            if len(set(y_local[ix])) < 2:
                continue
            ww = np.concatenate([np.full(len(block), 1.0 / len(block)) for block in blocks])
            boot_local.append(balanced_accuracy_score(y_local[ix], pred_local[ix], sample_weight=ww))
        composition_shortcut_rows.append({
            "target": target, "predictors": "ground-truth 15-lineage section-aggregate composition C only",
            "n_sections": len(frame), "n_donors": frame.donor_id.nunique(), "n_classes": len(classes_local),
            "class_counts_sections": ";".join(f"{k}:{v}" for k,v in frame.label.value_counts().sort_index().items()),
            "balanced_accuracy": float(ba_local), "macro_auroc": float(auc_local),
            "bootstrap_95ci_low": float(np.quantile(boot_local, .025)) if boot_local else np.nan,
            "bootstrap_95ci_high": float(np.quantile(boot_local, .975)) if boot_local else np.nan,
            "unseen_test_classes_across_folds": ";".join(sorted(unsupported_local)),
            "outer_group": "donor_id", "inference_unit": "donor_id",
        })
    composition_shortcut_rows.append({
        "target": "donor_id", "predictors": "ground-truth composition C only",
        "n_sections": len(section_comp), "n_donors": section_comp.donor_id.nunique(),
        "n_classes": section_comp.donor_id.nunique(), "balanced_accuracy": np.nan, "macro_auroc": np.nan,
        "bootstrap_95ci_low": np.nan, "bootstrap_95ci_high": np.nan,
        "unseen_test_classes_across_folds": "all held-out donors",
        "outer_group": "donor_id", "inference_unit": "donor_id",
        "interpretation": "NOT_ESTIMABLE: donor identity is held out by design",
    })
    pd.DataFrame(composition_shortcut_rows).to_csv(
        OUT / "composition_metadata_shortcut.tsv", sep="\t", index=False)

    # Core eight target coverage: thresholds are fixed before any expression
    # outcome analysis. Coverage uses signature crosswalk generated in Phase 0C-V.
    coverage = pd.read_csv(V / "signature_coverage.tsv", sep="\t")
    spec = [
        ("proliferation", "core", "proliferation"), ("hypoxia", "core", "hypoxia"),
        ("EMT / transitional epithelial", "core", "EMT"),
        ("ECM remodeling", "core", "ECM remodeling"),
        ("IFN response", "core", "IFN-gamma response"),
        ("inflammatory response", "core", "inflammatory response"),
        ("cytotoxicity", "core", "cytotoxicity"),
        ("fibrosis / myofibroblast", "lung_specific", "fibroblast activation"),
    ]
    target_rows = []
    for name, domain, program in spec:
        row = coverage[(coverage.domain.eq(domain)) & coverage.program.eq(program)]
        if row.empty:
            raise RuntimeError(f"Target signature missing: {domain}/{program}")
        row = row.iloc[0]
        n = int(row.n_panel_genes); frac = float(row.coverage_fraction)
        status = "PRIMARY_USABLE" if n >= 10 and frac >= .20 else (
            "SECONDARY_USABLE" if n >= 5 and frac >= .10 else "LOW_COVERAGE")
        target_rows.append({"candidate_program": name, "reference_signature": program,
                            "n_reference_genes": int(row.n_total_unique_genes), "n_panel_genes": n,
                            "coverage_fraction": frac, "coverage_status_phase0cv": row.coverage_status,
                            "phase1_target_class": status, "available_panel_genes": row.retained_genes,
                            "classification_rule": "primary: >=10 genes and >=20%; secondary: >=5 and >=10%; otherwise low",
                            "selection_basis": "predefined signature-to-343-gene-panel overlap; no expression values or correlation used"})
    target_df = pd.DataFrame(target_rows)
    target_df.to_csv(CONFIG / "frozen_phase1_candidate_programs.tsv", sep="\t", index=False)

    # Within-lineage program coverage and donor-level cell support at >=20 cells
    # in >=10 distinct 150-um regions per donor.
    lin_specs = [
        ("epithelial", "proliferation", "core", "proliferation"),
        ("epithelial", "injury", "lung_specific", "alveolar epithelial injury"),
        ("epithelial", "transitional epithelial", "lung_specific", "transitional epithelial state"),
        ("fibroblast", "ECM remodeling", "lung_specific", "ECM deposition"),
        ("fibroblast", "fibroblast activation", "lung_specific", "fibroblast activation"),
        ("fibroblast", "myofibroblast", "lung_specific", "myofibroblast"),
        ("macrophage", "inflammatory", "lung_specific", "macrophage inflammatory"),
        ("macrophage", "profibrotic", "lung_specific", "profibrotic macrophage"),
        ("T_cell", "IFN / cytotoxic state", "core", "IFN-gamma response"),
    ]
    lineage_rows = []
    for lineage, target, domain, program in lin_specs:
        row = coverage[(coverage.domain.eq(domain)) & coverage.program.eq(program)].iloc[0]
        n_genes, frac = int(row.n_panel_genes), float(row.coverage_fraction)
        cls = "PRIMARY_USABLE" if n_genes >= 5 and frac >= .50 else (
            "SECONDARY_USABLE" if n_genes >= 3 and frac >= .25 else "LOW_COVERAGE")
        cell_col = {"T_cell": None}.get(lineage, f"n_{lineage}")
        if cell_col:
            elig = micro[micro[cell_col] >= 20].groupby("donor_id").size()
            eligible_donors = int((elig >= 10).sum())
            eligible_regions = int(micro[cell_col].ge(20).sum())
            eligible_sections = int(micro.loc[micro[cell_col].ge(20), "sample_id"].nunique())
        else:
            eligible_donors=eligible_regions=eligible_sections=0
        cell_support = "FEASIBLE" if eligible_donors >= 5 else "WEAK" if eligible_donors >= 3 else "INFEASIBLE"
        lineage_rows.append({"lineage": lineage, "target": target, "reference_signature": program,
                             "n_reference_genes": int(row.n_total_unique_genes), "n_panel_genes": n_genes,
                             "coverage_fraction": frac, "coverage_class": cls,
                             "eligible_donors_ge_10_regions_ge_20_cells": eligible_donors,
                             "eligible_sections": eligible_sections, "eligible_regions": eligible_regions,
                             "cell_support_class": cell_support,
                             "target_feasibility": "FEASIBLE" if cls in {"PRIMARY_USABLE", "SECONDARY_USABLE"} and cell_support == "FEASIBLE" else "LIMITED_OR_INFEASIBLE",
                             "available_panel_genes": row.retained_genes,
                             "classification_rule": "primary: >=5 genes and >=50%; secondary: >=3 and >=25%; otherwise low; cell support: >=20 cells in >=10 regions for >=5 donors",
                             "T_cell_note": "No lineage-specific T-cell state signature with adequate coverage; broad IFN/cytotoxicity sets remain low-coverage." if lineage == "T_cell" else ""})
    lineage_df = pd.DataFrame(lineage_rows)
    lineage_df.to_csv(CONFIG / "frozen_within_lineage_targets.tsv", sep="\t", index=False)

    # Freeze Z covariates, image preprocessing and unit/section rules.
    z = pd.DataFrame([
        {"variable": "normalized_x", "known_at_inference": "YES", "biological_or_technical": "anatomical_position", "risk_of_absorbing_signal": "MODERATE", "included_primary": "YES", "included_sensitivity": "YES", "justification": "fixed image coordinates are available at inference; included to control systematic within-section position"},
        {"variable": "normalized_y", "known_at_inference": "YES", "biological_or_technical": "anatomical_position", "risk_of_absorbing_signal": "MODERATE", "included_primary": "YES", "included_sensitivity": "YES", "justification": "paired with normalized x as a minimal positional nuisance adjustment"},
        {"variable": "disease_status", "known_at_inference": "YES_IF_TASK_IS_CONDITIONED_ON_COHORT_LABEL", "biological_or_technical": "biological", "risk_of_absorbing_signal": "HIGH", "included_primary": "YES_FOR_DISEASE_ADJUSTED_ESTIMAND", "included_sensitivity": "YES", "justification": "controls pooled PF/control differences; mandatory PF-only analysis remains separate"},
        {"variable": "affected_status", "known_at_inference": "YES_FOR_SECTION_METADATA", "biological_or_technical": "biological_local_severity", "risk_of_absorbing_signal": "VERY_HIGH", "included_primary": "NO", "included_sensitivity": "YES", "justification": "local severity may be part of the morphology signal; report as sensitivity only"},
        {"variable": "TMA", "known_at_inference": "YES", "biological_or_technical": "technical_batch_proxy", "risk_of_absorbing_signal": "HIGH", "included_primary": "NO", "included_sensitivity": "YES", "justification": "audit indicates disease/TMA association; adjust only as sensitivity if overlap supports within-donor-held-out comparison"},
        {"variable": "run", "known_at_inference": "YES", "biological_or_technical": "technical_batch_proxy", "risk_of_absorbing_signal": "HIGH", "included_primary": "NO", "included_sensitivity": "YES", "justification": "run is not exactly aliased with TMA/disease; retain for sensitivity rather than loading primary Z"},
        {"variable": "collection_center", "known_at_inference": "YES", "biological_or_technical": "source", "risk_of_absorbing_signal": "HIGH", "included_primary": "NO", "included_sensitivity": "NO", "justification": "strongly coupled to disease and duplicates donor/site identity; report descriptive shortcut only"},
        {"variable": "scanner", "known_at_inference": "NO_PER_IMAGE", "biological_or_technical": "technical", "risk_of_absorbing_signal": "UNRESOLVED", "included_primary": "NO", "included_sensitivity": "NO", "justification": "only publication-level scanner statement found, with no variation or verified per-image metadata"},
    ])
    z.to_csv(CONFIG / "frozen_covariates_Z.tsv", sep="\t", index=False)
    (CONFIG / "frozen_image_preprocessing.yaml").write_text(
        "primary:\n  image: author_registered_HE\n  color_space: raw_RGB\n  normalization: none\n  patch_unit_um: 150\n  grid_origin: fixed_image_origin_zero\n  orientation: identity_as_author_registered\n  fit_from_heldout_donor: false\nsensitivity:\n  - name: Macenko\n    parameters: {alpha_percentile: 1, beta_od_threshold: 0.15}\n    fit_reference: training_donors_only_within_each_outer_fold\n    target_stain_matrix: donor_equal_weighted_median_of_training_section_stain_vectors\n    target_concentration_percentiles: donor_equal_weighted_median_of_training_section_99th_percentiles\n    apply_to: heldout_donors_after_fit\n    heldout_statistics_used_for_fit: false\n  - name: grayscale\n    conversion: fixed_luminance_Y_0.2126R_0.7152G_0.0722B\nshortcut_controls:\n  - raw_color_only\n  - grayscale_only\nselection_rule: freeze_before_any_Phase1_performance_comparison\n", encoding="utf-8")
    (CONFIG / "frozen_morphology_scale.yaml").write_text(
        "primary_morphology_unit: 150um_microregion\nsensitivity_units_um: [100, 200]\nprohibited_primary_units:\n  - single_nucleus_crop\n  - single_cell_crop\n  - nucleus_shape_prediction\n  - cell_level_morphology_embedding\nrationale: single_cell_or_nucleus_level_registration_precision_not_validated\n", encoding="utf-8")
    (CONFIG / "frozen_phase1_estimand.yaml").write_text(
        "primary_estimand: DeltaR2_morph\nformula: R2(C + N + Z + X) - R2(C + N + Z)\noutcome_Y: gene_expression_or_predefined_molecular_program\nC:\n  definition: ground_truth_cell_composition\nN:\n  definition: ground_truth_neighborhood_composition\n  primary_radius_um: 150\n  sensitivity_radii_um: [100, 200]\nZ:\n  primary_spatial_covariates: [normalized_x, normalized_y]\n  disease_adjusted_estimand: include_disease_status\n  sensitivity_only: [affected_status, TMA, run]\nX: frozen_pathology_morphology_embedding\nmorphology_unit_um: 150\nouter_validation: donor_held_out\nsame_donor_crosses_folds: false\ninference_unit: donor_id\nprimary_uncertainty: donor_level_bootstrap\nmicroregion_level_bootstrap_allowed: false\nresidual_analysis:\n  required: true\n  procedure: outer_fold_cross_fitted\n  full_data_residualization_allowed: false\npreprocessing:\n  primary: raw_author_registered_HE\n  sensitivity: training_donor_fit_Macenko\n  color_only_and_grayscale_controls: mandatory\nselection_and_hyperparameter_tuning: nested_donor_level_CV\nphase1_execution: not_started_requires_separate_human_authorization\n",
        encoding="utf-8")
    (CONFIG / "frozen_section_selection.yaml").write_text(
        "default: include_all_26_admitted_sections\nouter_group: donor_id\npermanent_exclusions:\n  VUILD105MA1: FIXED_COORDINATE_SUPPORT_FAILURE\n  VUILD48LA1: FIXED_COORDINATE_SUPPORT_FAILURE\nfuture_technical_exclusion_rule:\n  - section_file_missing_after_retry\n  - SHA256_mismatch_against_frozen_manifest\n  - image_decode_failure_after_retry\n  - reproducible_failure_of_frozen_fixed_coordinate_QC\n  - never_exclude_for_image_quality_rank_or_model_performance\npaired_sections: retain_all_and_keep_within_same_outer_donor_fold\n", encoding="utf-8")

    # TMA group counts, donor same-donor LF/MF pair flags, and all major outputs.
    sample_summary = pd.DataFrame([{
        "frozen_candidates": len(primary), "admitted_sections": len(admitted),
        "permanently_excluded_sections": len(excluded), "independent_donors": donors.donor_id.nunique(),
        "control_donors": int(donors.disease.eq("control").sum()),
        "PF_donors": int(donors.disease.eq("pulmonary_fibrosis").sum()),
        "admitted_sections_within_donor_multiple_cores": int(donors.n_sections.gt(1).sum()),
        "donors_with_LF_MF_pair": int(pd.DataFrame(donor_detail).has_less_and_more_affected_pair.sum()),
        "all_outer_splits_grouped_by": "donor_id",
    }])
    sample_summary.to_csv(OUT / "frozen_cohort_summary.tsv", sep="\t", index=False)
    print(f"Frozen cohort: {len(admitted)} sections / {donors.donor_id.nunique()} donors; 5 donor folds; outputs under {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    build_freeze_and_static_audits()
