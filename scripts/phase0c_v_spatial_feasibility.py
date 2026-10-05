#!/usr/bin/env python3
"""Outcome-blind cell ontology and microregion feasibility for GSE250346.

Uses only the corrected author cell metadata and published primary TMA1-4.
No expression scores, prediction, clustering, transform fitting, or cell
exclusions are introduced here.
"""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/phase0/phase0c_v"
CELL_FILE = BASE / "vannan_primary_cell_metadata.tsv.gz"
CROSSWALK = BASE / "vannan_sample_crosswalk.tsv"
TIFF_META = BASE / "registered_he_tiff_metadata.tsv"

COARSE = {
    "AT1": "epithelial", "AT2": "epithelial", "Basal": "epithelial",
    "Goblet": "epithelial", "KRT5-/KRT17+": "epithelial",
    "Multiciliated": "epithelial", "PNEC": "epithelial",
    "Proliferating Airway": "epithelial", "Proliferating AT2": "epithelial",
    "RASC": "epithelial", "Secretory": "epithelial",
    "Transitional AT2": "epithelial",
    "Arteriole": "endothelial", "Capillary": "endothelial",
    "Lymphatic": "endothelial", "Venous": "endothelial",
    "Activated Fibrotic FBs": "fibroblast", "Adventitial FBs": "fibroblast",
    "Alveolar FBs": "fibroblast", "Inflammatory FBs": "fibroblast",
    "Myofibroblasts": "fibroblast", "Proliferating FBs": "fibroblast",
    "Subpleural FBs": "fibroblast",
    "SMCs/Pericytes": "smooth_muscle_pericyte",
    "Alveolar Macrophages": "macrophage", "Interstitial Macrophages": "macrophage",
    "Macrophages - IFN-activated": "macrophage", "SPP1+ Macrophages": "macrophage",
    "Monocytes/MDMs": "monocyte", "Neutrophils": "neutrophil",
    "CD4+ T-cells": "CD4_T", "Tregs": "CD4_T", "CD8+ T-cells": "CD8_T",
    "B cells": "B", "Proliferating B cells": "B", "Plasma": "plasma",
    "Mast": "mast", "Mesothelial": "other_stromal",
    "Basophils": "other_immune", "cDCs": "other_immune",
    "Langerhans cells": "other_immune", "Migratory DCs": "other_immune",
    "NK/NKT": "other_immune", "pDCs": "other_immune",
    "Proliferating Myeloid": "other_immune", "Proliferating NK/NKT": "other_immune",
    "Proliferating T-cells": "other_immune",
}
ORDER = ["epithelial", "endothelial", "fibroblast", "smooth_muscle_pericyte",
         "macrophage", "monocyte", "neutrophil", "CD4_T", "CD8_T", "B",
         "plasma", "mast", "other_immune", "other_stromal", "other"]

PROGRAMS = {
    "epithelial": "Core proliferation POOR (12/200); injury markers 5/10 coverage-only",
    "fibroblast": "Core ECM POOR (34/321); fibroblast activation markers 7/13 and myofibroblast markers 3/10 are coverage-only",
    "macrophage": "Macrophage inflammatory markers 6/10 and profibrotic markers 4/10 are coverage-only",
    "CD4_T": "Core IFN-gamma POOR (24/200); cytotoxicity POOR (18/155); no validated adequate T-cell state module",
    "CD8_T": "Core exhaustion POOR (13/200; mouse-origin source signature); cytotoxicity POOR (18/155)",
}


def point_in_polygon(x: np.ndarray, y: np.ndarray, poly: np.ndarray) -> np.ndarray:
    inside = np.zeros(len(x), dtype=bool)
    j = len(poly) - 1
    for i in range(len(poly)):
        xi, yi = poly[i]
        xj, yj = poly[j]
        crossing = ((yi > y) != (yj > y)) & (
            x < (xj - xi) * (y - yi) / ((yj - yi) if yj != yi else np.finfo(float).eps) + xi
        )
        inside ^= crossing
        j = i
    return inside


def correction_audit(cross: pd.DataFrame) -> None:
    audit = pd.read_csv(BASE / "seurat_object_audit.tsv", sep="\t").set_index("key")["value"]
    inputs = pd.read_csv(BASE / "phase0c_v_input_manifest.tsv", sep="\t")
    source = {Path(r.file_path).name: r for r in inputs.itertuples(index=False)}
    ids = source["GSE250346_remove_nuclei_VUILD105MA1.csv.gz"]
    poly = source["GSE250346_remove_coordinates_VUILD105MA1.csv.gz"]
    rds = source["GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds"]
    rows = [
        {"correction_type": "author_removed_nuclei_IDs", "sample_id": "VUILD105MA1",
         "listed_items": int(audit["n_exclusion_nuclei_ids"]),
         "matched_corrected_object_items": int(audit["matched_removed_nuclei"]),
         "excluded_from_joined_primary_table": 0,
         "interpretation": "All listed IDs were already absent from the latest explicitly corrected author RDS; no unrelated cells were removed.",
         "source_file": ids.file_path, "source_sha256": ids.sha256, "source_size_bytes": ids.size_bytes,
         "replacement_object": rds.file_path, "replacement_sha256": rds.sha256},
        {"correction_type": "author_coordinate_exclusion_polygon", "sample_id": "VUILD105MA1",
         "listed_items": int(audit["n_exclusion_coordinate_vertices"]),
         "matched_corrected_object_items": int(audit["matched_removed_polygon_cells"]),
         "excluded_from_joined_primary_table": int(audit["matched_removed_polygon_cells"]),
         "interpretation": "The five-vertex author polygon overlapped 22 cells in the corrected RDS; those cells were removed before the joined primary table was created.",
         "source_file": poly.file_path, "source_sha256": poly.sha256, "source_size_bytes": poly.size_bytes,
         "replacement_object": rds.file_path, "replacement_sha256": rds.sha256},
    ]
    pd.DataFrame(rows).to_csv(BASE / "correction_application_audit.tsv", sep="\t", index=False)


def mapping_table(cells: pd.DataFrame) -> None:
    unknown = sorted(set(cells.original_cell_type.dropna()) - set(COARSE))
    if unknown:
        raise ValueError(f"Outcome-blind coarse ontology missing author labels: {unknown}")
    cells["morphoid_coarse_ontology"] = cells.original_cell_type.map(COARSE)
    map_df = (cells.groupby(["original_cell_type", "lineage", "morphoid_coarse_ontology"],
                            dropna=False).size().rename("n_cells_primary").reset_index())
    map_df["mapping_basis"] = "author cell type + broad lineage; outcome-blind frozen grouping"
    map_df.to_csv(BASE / "vannan_celltype_mapping.tsv", sep="\t", index=False)
    with gzip.open(CELL_FILE, "wt", encoding="utf-8", newline="") as f:
        cells.to_csv(f, sep="\t", index=False)


def image_bounds(cells: pd.DataFrame) -> pd.DataFrame:
    cross = pd.read_csv(CROSSWALK, sep="\t")
    tm = pd.read_csv(TIFF_META, sep="\t")
    primary = cross[cross.published_primary == "YES"]
    rows = []
    for _, sm in primary.iterrows():
        sub = cells[cells.sample_id == sm.sample_id]
        img = tm[tm.sample_id == sm.sample_id]
        if len(img) != 1:
            raise ValueError(f"Expected exactly one registered H&E row for {sm.sample_id}")
        img = img.iloc[0]
        has_unit = img.ImageJ_spatial_unit in {"micron", "microns", "um"}
        px_x = float(img.physical_pixel_size_x_um) if pd.notna(img.physical_pixel_size_x_um) else np.nan
        px_y = float(img.physical_pixel_size_y_um) if pd.notna(img.physical_pixel_size_y_um) else np.nan
        if has_unit and np.isfinite(px_x) and np.isfinite(px_y):
            use_x, use_y, scale_basis = px_x, px_y, "per-image TIFF ImageJ micron resolution"
        else:
            # Diagnostic only. Do not use this shared-pipeline scale to pass R1/R2.
            use_x = use_y = 0.2125
            scale_basis = "diagnostic consensus scale 0.2125 um/px; TIFF unit missing; not admissible for R1"
        xp = sub.x_um.to_numpy(float) / use_x
        yp = sub.y_um.to_numpy(float) / use_y
        inb = (xp >= 0) & (xp < int(img.width_px)) & (yp >= 0) & (yp < int(img.height_px))
        valid = has_unit and np.isfinite(px_x) and np.isfinite(px_y)
        r2 = float(inb.mean()) if valid else np.nan
        rows.append({
            "sample_id": sm.sample_id, "donor_id": sm.donor_id, "GSM": sm.GSM,
            "TMA": sm.TMA, "disease_status": sm.disease_status,
            "registered_HE": img.file, "registered_HE_sha256": img.sha256,
            "image_width_px": int(img.width_px), "image_height_px": int(img.height_px),
            "tiff_scale_unit": img.ImageJ_spatial_unit,
            "pixel_size_x_um": px_x, "pixel_size_y_um": px_y,
            "scale_basis": scale_basis, "n_cells": len(sub),
            "centroids_inside_image_fraction_R2": r2,
            "centroids_inside_image_fraction_consensus_diagnostic": float(inb.mean()),
            "R2_threshold": 0.95,
            "R2_status": ("PASS" if r2 >= 0.95 else "FAIL") if valid else "UNRESOLVED",
            "x_um_min": float(sub.x_um.min()), "x_um_max": float(sub.x_um.max()),
            "y_um_min": float(sub.y_um.min()), "y_um_max": float(sub.y_um.max()),
            "transform": "identity orientation, zero translation; author H&E scale only; no fitted/refit transform",
        })
    out = pd.DataFrame(rows)
    out.to_csv(BASE / "spatial_bounds_diagnostic.tsv", sep="\t", index=False)
    reg = out[["sample_id", "donor_id", "registered_HE", "n_cells",
               "centroids_inside_image_fraction_R2", "R2_status"]].copy()
    reg = reg.rename(columns={"centroids_inside_image_fraction_R2": "centroids_inside_image"})
    reg["tissue_coverage"] = np.nan
    reg["dapi_rho"] = np.nan
    reg["visual_pass"] = "NOT_REVIEWED"
    reg["n_visual_patches"] = 0
    reg["qc_deviation"] = "No transform fitting; R3 image mask, R4 landmarks, R5 DAPI and R6 visual patch review pending/unavailable in bounded inputs."
    reg["registration_verdict"] = "UNRESOLVED"
    reg["reason"] = out.R2_status + "; frozen R4/R5/R6 gates incomplete"
    reg.to_csv(BASE / "registration_qc.tsv", sep="\t", index=False)
    return out


def build_regions(cells: pd.DataFrame, grid_um: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = cells[["cell_id", "sample_id", "donor_id", "disease_status", "affected_status",
               "x_um", "y_um", "morphoid_coarse_ontology"]].copy()
    d["region_x"] = np.floor(d.x_um / grid_um).astype(np.int32)
    d["region_y"] = np.floor(d.y_um / grid_um).astype(np.int32)
    grouped = d.groupby(["sample_id", "donor_id", "disease_status", "affected_status",
                         "region_x", "region_y", "morphoid_coarse_ontology"], observed=True).size()
    wide = grouped.unstack("morphoid_coarse_ontology", fill_value=0)
    for typ in ORDER:
        if typ not in wide.columns:
            wide[typ] = 0
    wide = wide[ORDER]
    wide.index.names = ["sample_id", "donor_id", "disease_status", "affected_status", "region_x", "region_y"]
    regions = wide.reset_index()
    count_cols = [f"n_{typ}" for typ in ORDER]
    regions.columns = list(regions.columns[:6]) + count_cols
    regions["grid_um"] = grid_um
    regions["region_id"] = (regions.sample_id.astype(str) + f"_{grid_um}_" +
                            regions.region_x.astype(str) + "_" + regions.region_y.astype(str))
    regions["total_cells"] = regions[count_cols].sum(axis=1).astype(int)
    # Tissue-mask area is left NA until registered H&E tissue support is fully audited.
    regions["tissue_area_um2"] = np.nan
    props = regions[count_cols].to_numpy(float)
    totals = props.sum(axis=1)
    p = np.divide(props, totals[:, None], out=np.zeros_like(props), where=totals[:, None] > 0)
    regions["composition_entropy_bits"] = -np.sum(np.where(p > 0, p * np.log2(np.maximum(p, 1e-300)), 0), axis=1)
    for j, typ in enumerate(ORDER):
        regions[f"p_{typ}"] = p[:, j]
    occupied = {(r.sample_id, int(r.region_x), int(r.region_y)) for r in regions.itertuples()}
    neighbor_counts = []
    for r in regions.itertuples():
        neighbor_counts.append(sum((r.sample_id, int(r.region_x) + dx, int(r.region_y) + dy) in occupied
                                   for dx in (-1, 0, 1) for dy in (-1, 0, 1)
                                   if not (dx == 0 and dy == 0)))
    regions["n_adjacent_nonempty_regions_8"] = neighbor_counts
    regions["has_nonempty_neighbor"] = regions.n_adjacent_nonempty_regions_8 > 0

    bydonor = regions.groupby(["donor_id", "disease_status", "grid_um"], observed=True)
    summary = bydonor.agg(
        n_regions=("region_id", "size"), n_sections=("sample_id", "nunique"),
        median_cells=("total_cells", "median"), q1_cells=("total_cells", lambda x: x.quantile(.25)),
        q3_cells=("total_cells", lambda x: x.quantile(.75)),
        median_composition_entropy_bits=("composition_entropy_bits", "median"),
        regions_with_neighbor=("has_nonempty_neighbor", "sum"),
        regions_without_neighbor=("has_nonempty_neighbor", lambda x: (~x).sum()),
    ).reset_index()
    for typ in ORDER:
        summary[f"regions_with_{typ}"] = bydonor[f"n_{typ}"].apply(lambda x: int((x > 0).sum())).to_numpy()
    summary["IQR_cells"] = summary.q3_cells - summary.q1_cells
    summary["neighbor_feasible_fraction"] = summary.regions_with_neighbor / summary.n_regions
    return regions, summary


def lineage_feasibility(regions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    primary = regions[regions.grid_um == 150]
    for typ in ORDER:
        col = f"n_{typ}"
        for threshold in (20, 10):
            eligible = primary[primary[col] >= threshold]
            donor_counts = eligible.groupby("donor_id").size()
            rows.append({
                "lineage": typ, "target_cell_threshold_per_150um_region": threshold,
                "eligible_regions": len(eligible), "eligible_donors": eligible.donor_id.nunique(),
                "eligible_sections": eligible.sample_id.nunique(),
                "median_target_cells_per_eligible_region": float(eligible[col].median()) if len(eligible) else np.nan,
                "min_regions_per_eligible_donor": int(donor_counts.min()) if len(donor_counts) else 0,
                "median_regions_per_eligible_donor": float(donor_counts.median()) if len(donor_counts) else np.nan,
                "program_coverage_status": PROGRAMS.get(typ, "No prespecified within-lineage molecular program in this audit"),
                "donor_threshold_5_met": "YES" if eligible.donor_id.nunique() >= 5 else "NO",
                "registration_gate_status": "UNRESOLVED; see registration_qc.tsv",
                "interpretation": "Cell-count feasibility only; does not fit scores or expression models.",
            })
    return pd.DataFrame(rows)


def main() -> None:
    BASE.mkdir(parents=True, exist_ok=True)
    cross = pd.read_csv(CROSSWALK, sep="\t")
    cells = pd.read_csv(CELL_FILE, sep="\t")
    if cells.cell_id.duplicated().any():
        raise ValueError("Duplicate cell_id in corrected primary cell metadata")
    if cells.donor_id.isna().any() or cells.sample_id.isna().any():
        raise ValueError("Missing sample/donor identifiers in primary cell table")
    if not (cells.x_um.map(np.isfinite) & cells.y_um.map(np.isfinite)).all():
        raise ValueError("Non-finite cell coordinates in corrected primary table")
    sample_donor = cross.set_index("sample_id").donor_id
    mismatch = cells.donor_id != cells.sample_id.map(sample_donor)
    if mismatch.any():
        raise ValueError(f"Sample/donor crosswalk mismatches: {int(mismatch.sum())}")
    correction_audit(cross)
    mapping_table(cells)
    bounds = image_bounds(cells)
    all_regions, all_summaries = [], []
    for grid in (150, 100, 200):
        region, summary = build_regions(cells, grid)
        all_regions.append(region)
        all_summaries.append(summary)
    regions = pd.concat(all_regions, ignore_index=True)
    summaries = pd.concat(all_summaries, ignore_index=True)
    regions.to_csv(BASE / "microregion_feasibility.tsv", sep="\t", index=False, na_rep="NA")
    summaries.to_csv(BASE / "donor_microregion_summary.tsv", sep="\t", index=False, na_rep="NA")
    lineage_feasibility(regions).to_csv(BASE / "within_lineage_feasibility.tsv", sep="\t", index=False, na_rep="NA")
    summary = {
        "primary_cells": int(len(cells)), "primary_sections": int(cells.sample_id.nunique()),
        "primary_donors": int(cells.donor_id.nunique()),
        "coarse_ontology_labels": int(cells.morphoid_coarse_ontology.nunique()),
        "grid_region_counts": {str(g): int((regions.grid_um == g).sum()) for g in (150, 100, 200)},
        "grid_median_cells": {str(g): float(regions.loc[regions.grid_um == g, "total_cells"].median()) for g in (150, 100, 200)},
        "grid_fraction_with_nonempty_8_neighbor": {str(g): float(regions.loc[regions.grid_um == g, "has_nonempty_neighbor"].mean()) for g in (150, 100, 200)},
        "tissue_area": "NA pending section-level H&E tissue mask and registration audit",
        "composition_entropy_between_donor_sd": float(summaries[summaries.grid_um == 150].median_composition_entropy_bits.std(ddof=1)),
        "r2_pass_sections": int((bounds.R2_status == "PASS").sum()),
        "r2_fail_sections": int((bounds.R2_status == "FAIL").sum()),
        "r2_unresolved_sections": int((bounds.R2_status == "UNRESOLVED").sum()),
    }
    (BASE / "spatial_feasibility_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
