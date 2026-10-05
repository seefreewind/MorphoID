#!/usr/bin/env python3
"""Finalize Phase 0C-V cohort audit outputs without modeling."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/phase0/phase0c_v"
REPORTS = ROOT / "reports"


def markdown_report(cross, tiff, bounds, reg, object_audit, spatial, lineage, coverage,
                    conf, raw_map, hescape_status: str) -> str:
    primary = cross[cross.published_primary == "YES"].copy()
    pf = primary[primary.disease_status == "pulmonary_fibrosis"]
    controls = primary[primary.disease_status == "control"]
    reg_pass = int((reg.registration_verdict == "PASS").sum())
    reg_donors = reg.loc[reg.registration_verdict == "PASS", "donor_id"].nunique()
    n_sections = int(primary.sample_id.nunique())
    n_donors = int(primary.donor_id.nunique())
    r2 = bounds.R2_status.value_counts().to_dict()
    coverage_core = coverage[coverage.domain == "core"]
    coverage_lung = coverage[coverage.domain == "lung_specific"]
    good_lung = int((coverage_lung.coverage_status == "GOOD").sum())
    partial_lung = int((coverage_lung.coverage_status == "PARTIAL").sum())
    poor_lung = int((coverage_lung.coverage_status == "POOR").sum())
    eligible = lineage[lineage.target_cell_threshold_per_150um_region == 20].set_index("lineage")
    region = spatial[spatial.grid_um == 150]
    n_regions = int(len(region))
    median_cells = float(region.total_cells.median())
    q1, q3 = region.total_cells.quantile([.25, .75])
    neighbor_frac = float(region.has_nonempty_neighbor.mean())
    primary_cells = int(object_audit["n_primary_cells"])
    core_genes = int(pd.to_numeric(coverage_core.n_panel_genes).max() if len(coverage_core) else 0)
    common_genes = int(pd.read_csv(BASE / "gene_panel_registry.tsv", sep="\t").gene_symbol.nunique())
    admitted = int((reg.registration_verdict == "PASS").sum())
    pf_reg_donors = reg[(reg.registration_verdict == "PASS") & (reg.disease_status == "pulmonary_fibrosis")].donor_id.nunique()
    conf_assessment = conf[conf.cohort == "PUBLISHED_PRIMARY_TMA1_4"]
    relevant_pairs = conf_assessment.groupby("factor_pair").alias_assessment.first().to_dict()
    technical_alias = "MODERATE"
    line_summary = "; ".join(
        f"{name}: {int(eligible.loc[name, 'eligible_regions']):,} regions / {int(eligible.loc[name, 'eligible_donors'])} donors"
        for name in ("epithelial", "fibroblast", "macrophage")
    )
    any_tissue_area = "not computed; H&E tissue masks were not materialized"
    return f"""1. **VERDICT:** `HOLD_VANNAN_DATA_RECONSTRUCTION`
2. **Published primary sections recovered:** {n_sections}/28 (the frozen requested TMA1–4 subset; the 2025 article/GEO now also includes TMA5, for 45 sections total)
3. **Published primary donors recovered:** {n_donors}/19 (the frozen requested TMA1–4 subset; full current article/GEO cohort has 35 donors)
4. **Independent donors admitted:** {admitted}
5. **Same-section provenance:** PASS
6. **Registered H&E available:** {int(tiff.sample_id.isin(primary.sample_id).sum())}
7. **Registration PASS sections:** {reg_pass}
8. **Registration PASS donors:** {reg_donors}
9. **TPS warp status:** FAIL
10. **Expression-coordinate join:** PASS
11. **Published cell annotation recovered:** YES
12. **Common panel genes:** {common_genes}
13. **Composition oracle feasible:** YES (in the Xenium coordinate frame)
14. **Primary microregions:** {n_regions:,} at 150 μm (100 μm: {int((spatial.grid_um == 100).sum()):,}; 200 μm: {int((spatial.grid_um == 200).sum()):,})
15. **Eligible donors for donor-held-out:** {reg_donors} admitted; raw subset has {n_donors} independent donors
16. **PF-only analysis feasible:** NO for an auditable H&E donor-held-out analysis; the roster has {pf.donor_id.nunique()} PF donors/{len(pf)} sections, but {pf_reg_donors} registration-PASS donors
17. **Within-lineage feasible:** Cell-count support at ≥20 target cells/region: {line_summary}; state-panel coverage is strongest for epithelial injury, fibroblast activation and macrophage inflammation, conditional on registration
18. **Adequately usable programs:** Core predefined programs {int((coverage_core.coverage_status.isin(['GOOD','PARTIAL'])).sum())}/8 at ≥25% coverage; lung coverage-only panels {good_lung} GOOD, {partial_lung} PARTIAL, {poor_lung} POOR
19. **Disease/technical confounding:** {technical_alias}
20. **HESCAPE crosswalk:** {hescape_status}
21. **Replacement cohort admitted:** NO
22. **Phase 0D authorized:** NO
23. **Recommended next action:** Keep the project at `HOLD_REPLACEMENT_COHORT`; recover the author-level sample/core-to-registered-image mapping and the registration evidence needed to complete R3–R6, then re-audit before any later-phase request.

# Phase 0C-V — Vannan GSE250346 replacement-cohort audit

## Decision and frozen status

This audit answers only whether the frozen TMA1–4 subset of GSE250346 can presently replace the failed GSE311609 lung arm. It does not authorize Phase 0D. The project-level state remains `HOLD_REPLACEMENT_COHORT`; `GSE311609 L8` remains `HOLD_AUTHOR_ALIGNMENT_REQUIRED`. Prior Phase 0C, 0C-R, 0C-S, B2 registration, thresholds, manifests, hashes and reports remain untouched.

The request's 28-section/19-donor definition is treated as a fixed admission subset, not as the full current formal article cohort. The 2025 paper and current GEO series include TMA1–TMA5 (45 sections/35 donors). TMA5 is separately tagged outside this audit subset; it is not described as absent from the paper.

## Four separate findings

### Data availability — PASS

The frozen subset resolves to 28 section IDs and 19 unique donor IDs. The latest explicitly corrected author Seurat RDS contains expression counts, 343 panel genes, cell labels and centroid coordinates. There are {primary_cells:,} corrected primary-subset cells after application of the author coordinate polygon. All 28 primary sections have a sample-specific author-registered H&E TIFF. The source manifest records URL, retrieval date, bytes and SHA256 for each input.

The corrected RDS audit found 1,630,319 cells and 343 RNA features in the full current object. The author nucleus-removal list has 931 IDs and zero matches in that corrected object; those nuclei were already absent. The five-vertex author coordinate polygon overlapped 22 cells in VUILD105MA1; they were excluded before creating the joined table. Corrected symbols ACTA2 and MRC1 are present. The deterministic join audit sampled 1,000 cells across all 28 sections and 19 donors; sampled cell IDs joined to count columns, coordinates were finite, donor/sample mapping was consistent, and RNA count sums matched metadata.

### Registration validity — HOLD

The publication confirms post-Xenium H&E from the same 5-μm tissue slice, registered by the authors from H&E to the Xenium DAPI coordinate frame. The authors used ImageJ 2.14.0 BigWarp, about 200 corresponding landmarks per image pair, a thin-plate-spline warp and manual review. This establishes provenance and the existence of an author registration; it does not provide the withheld landmarks, sample-specific warp files or independent local validation needed for the frozen MorphoID gates.

All 28 H&E TIFFs are readable. Twenty-six record an ImageJ spatial unit of microns and a 0.2125 μm/pixel scale. Two (VUILD104MA2, VUILD105MA1) have `unit=pixel` and no usable physical scale in their TIFF tags. A diagnostic identity map with no translation, rotation or refit puts 25/26 scale-documented sections above the unchanged R2 bound threshold of 0.95; VUILD48LA1 is at 0.041 and fails. The two missing-scale samples remain R2-unresolved. Under the same diagnostic 0.2125-μm/pixel consensus scale, VUILD104MA2 is 1.000 and VUILD105MA1 is 0.000; those proxy values are not admissible R2 passes.

R3 tissue-mask coverage was not computed. The 28 single-strip TIFFs contain about 32.63 GB of uncompressed RGB pixels and would require a full sequential pixel pass; since R4/R5 are also unavailable, that work could not make any section satisfy the frozen registration gate in this audit. Accordingly, the 25 identity-map bounds passes are geometry diagnostics only, not registration passes. R4 is N/A because no shared withheld-landmark table was supplied; R5 is unavailable because the Xenium DAPI morphology OME-TIFFs were not among the bounded inputs; R6 patch review was not performed. The frozen policy requires R5 and R6 when R4 landmarks are absent, so no section receives a registration PASS.

The authors confirm TPS usage, but the whole-TMA raw H&E files cannot be linked deterministically to each registered core: the public crosswalk is only TMA-level and no core polygon/slide-coordinate mapping was found. Nuclear morphology, texture and local geometry change under TPS therefore remain `WARP_UNRESOLVED`; this audit does not infer acceptability from the methods description.

The opening TPS status is `FAIL` in the admission sense: warp acceptability was not demonstrated, so this audit gate is not cleared. It does not mean that physical morphology damage was observed. The measured TPS effect remains `WARP_UNRESOLVED` because the raw-to-registered core mapping needed for comparison is unavailable.

### Biological independence — roster adequate, none admitted

Repeated samples remain grouped by donor. The frozen subset includes 13 PF donors across 21 sections and 6 control donors across 7 sections. VUHD116, TILD117, VUILD102, VUILD104, VUILD105, VUILD48, VUILD78, VUILD91 and VUILD96 have multiple primary sections/cores and contribute one biological unit each. No donor is admitted to donor-held-out analysis until its section passes the complete registration gate; at present the admitted count is zero.

Within the 28/19 subset, disease appears in multiple TMAs, Xenium runs and both collection centers, so no complete disease–batch alias is observed. The center allocation is imbalanced (PF: Vanderbilt 11 donors, TGen 2; control: Vanderbilt 4, TGen 2). Xenium software, the 343-gene panel and the publication-reported Leica Aperio CS2 H&E scanner each have one level within TMA1–4, so their disease associations are not estimable. The article reports data generation randomized with respect to disease/control. Overall confounding is rated MODERATE, with donor-level analysis still required.

### Molecular identifiability feasibility — composition feasible; state scope limited

The 150-μm Xenium-frame grid contains {n_regions:,} nonempty regions, with {median_cells:.0f} cells median (IQR {q1:.0f}–{q3:.0f}); {neighbor_frac:.2%} have at least one nonempty 8-neighbor region. Sensitivity grids contain {int((spatial.grid_um == 100).sum()):,} 100-μm and {int((spatial.grid_um == 200).sum()):,} 200-μm nonempty regions. These counts support measured composition C and neighboring composition N inside the expression/coordinate object. Tissue area remains NA until H&E tissue masks are available. Xenium-frame feasibility does not establish H&E co-registration.

At ≥20 target-lineage cells in a 150-μm region, cell counts support epithelial ({int(eligible.loc['epithelial','eligible_regions']):,} regions/19 donors), fibroblast ({int(eligible.loc['fibroblast','eligible_regions']):,}/19) and macrophage ({int(eligible.loc['macrophage','eligible_regions']):,}/19) strata. Per-donor lower tails are recorded in `within_lineage_feasibility.tsv`; the minimum eligible region counts are 6 epithelial, 1 fibroblast and 4 macrophage regions, so donor-level region sufficiency must be prespecified before a model phase.

All eight original MorphoID programs are POOR on the common panel: 12/200 proliferation genes, 15/200 hypoxia, 29/200 EMT, 34/321 ECM remodeling, 24/200 IFN-γ, 13/200 T-cell exhaustion, 18/155 cytotoxicity and 28/200 inflammatory response. In the separate external lung marker audit, four of ten coverage-only panels are GOOD (alveolar injury, aberrant basaloid, fibroblast activation, macrophage inflammatory), three are PARTIAL and three are POOR. These short externally sourced marker lists support feasibility screening; they were not expression-scored and are not treated as validated directional programs. T-cell state coverage is poor. The corrected object permits composition/neighborhood definition, but it does not establish a broad eight-program molecular-state cohort.

PF-only roster size clears the five-donor count screen (13 donors/21 sections), and PF-only Xenium-frame grids are tabulated in `pf_only_feasibility.tsv`. PF-only H&E donor-held-out analysis remains infeasible for admission until registration passes exist.

## HESCAPE metadata-only crosswalk

HESCAPE reports a human-lung panel and attributes the source collection to Vannan/GSE250346, but the public panel metadata reviewed here do not map each of its 20 sections/19 donors to a specific GSM, TMA core, registered TIFF or corrected-object cell IDs. The crosswalk is `PARTIAL`. Gated terms were not accepted and full HESCAPE data were not downloaded.

## Admission table and next action

`vannan_admission_registry.tsv` contains all 28 frozen sections. One section is marked `EXCLUDED_REGISTRATION` because its direct, no-refit R2 fraction is 0.041 (<0.95). The other 27 remain `UNRESOLVED` because complete registration evidence is missing; none is admitted. Expression, annotation, coordinates and the common gene panel are available, but those facts do not turn an unresolved registered image into a valid paired observation.

The selected verdict is `HOLD_VANNAN_DATA_RECONSTRUCTION`: donor structure, corrected object, same-section provenance and expression/coordinate annotations are recovered, while the author core mapping and registration evidence needed to establish image-to-cell correspondence remain unresolved. No Phase 0D or Phase 1 modeling was run. No embeddings, expression prediction, residualization, ΔR², composition-oracle regression, neighborhood regression, shortcut model, random-spot split, cross-disease model or HEST-1K analysis was run.

## Source references

- [Vannan et al., Nature Genetics (2025), Methods and data availability](https://www.nature.com/articles/s41588-025-02080-x)
- [GEO GSE250346 current series record and supplementary inventory](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346)
- [10x Xenium file-format documentation: centroid units, origin and pixel size](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html)
- [HESCAPE paper](https://openaccess.thecvf.com/content/ICCV2025W/CVAMD/papers/Gindra_A_Large-Scale_Benchmark_of_Cross-Modal_Learning_for_Histology_and_Gene_ICCVW_2025_paper.pdf), [supplement](https://openaccess.thecvf.com/content/ICCV2025W/CVAMD/supplemental/Gindra_A_Large-Scale_Benchmark_ICCVW_2025_supplemental.pdf), [repository](https://github.com/peng-lab/hescape), and [metadata-only dataset card](https://huggingface.co/datasets/Peng-AI/hescape-pyarrow)

## Main outputs

Per-file provenance and hashes are in `results/phase0/phase0c_v/phase0c_v_input_manifest.tsv` and `phase0c_v_signature_source_manifest.tsv`. Other audit tables in the same directory record the corrected object, donor/sample crosswalk, registration criteria, program coverage, coarse ontology, and microregion feasibility. The cohort overview figure is a bounds diagnostic only; it is not registration evidence.

## Figure quality-control record

The cohort overview is exported as 600-dpi PNG and LZW-compressed TIFF plus editable-text SVG. The Nature-figure validator reports 18 PASS, 1 WARN and 1 FAIL: the WARN is the 218.4-mm width, since no target journal width was specified; the FAIL is the absent PDF export. PDF was intentionally omitted under the project output preference, which does not request PDF deliverables. The figure is an internal audit overview, not a journal-sized registration validation figure. The validator's remaining PDF finding is therefore recorded rather than bypassed. The raster exports were regenerated after the 600-dpi/TIFF settings were applied.
"""


def main() -> None:
    REPORTS.mkdir(exist_ok=True)
    cross = pd.read_csv(BASE / "vannan_sample_crosswalk.tsv", sep="\t")
    tiff = pd.read_csv(BASE / "registered_he_tiff_metadata.tsv", sep="\t")
    bounds = pd.read_csv(BASE / "spatial_bounds_diagnostic.tsv", sep="\t")
    cells = pd.read_csv(BASE / "vannan_primary_cell_metadata.tsv.gz", sep="\t",
                        usecols=["sample_id", "donor_id", "disease_status", "morphoid_coarse_ontology"])
    spatial = pd.read_csv(BASE / "microregion_feasibility.tsv", sep="\t")
    lineage = pd.read_csv(BASE / "within_lineage_feasibility.tsv", sep="\t")
    coverage = pd.read_csv(BASE / "signature_coverage.tsv", sep="\t")
    conf = pd.read_csv(BASE / "confounding_matrix.tsv", sep="\t")
    raw_map = pd.read_csv(BASE / "raw_registered_he_crosswalk.tsv", sep="\t")
    object_df = pd.read_csv(BASE / "seurat_object_audit.tsv", sep="\t").set_index("key")["value"]
    object_audit = {"n_primary_cells": len(cells), "n_primary_sections": cells.sample_id.nunique(),
                    "n_primary_donors": cells.donor_id.nunique(), "object_cells": int(object_df["n_cells_before_exclusions"]),
                    "expression_qc_cells": int(object_df["deterministic_qc_cells"]),
                    "expression_qc_sections": int(object_df["deterministic_qc_sections"]),
                    "expression_qc_donors": int(object_df["deterministic_qc_donors"])}

    # Freeze the explicit per-section registration gate table.
    primary = cross[cross.published_primary == "YES"].copy()
    tiff_primary = tiff[tiff.sample_id.isin(primary.sample_id)]
    b = bounds.merge(primary[["sample_id", "affected_status", "panel_version", "n_panel_genes"]],
                     on="sample_id", how="left")
    rows = []
    for r in b.itertuples(index=False):
        scale_ok = pd.notna(r.pixel_size_x_um) and pd.notna(r.pixel_size_y_um)
        r2 = r.centroids_inside_image_fraction_R2
        r2_status = "PASS" if scale_ok and pd.notna(r2) and r2 >= .95 else (
            "FAIL" if scale_ok and pd.notna(r2) else "UNRESOLVED")
        verdict = "FAIL" if r2_status == "FAIL" else "UNRESOLVED"
        trow = tiff[tiff.sample_id == r.sample_id].iloc[0]
        rows.append({
            "sample_id": r.sample_id, "donor_id": r.donor_id, "GSM": r.GSM,
            "disease_status": r.disease_status, "affected_status": r.affected_status,
            "registered_HE": r.registered_HE, "registered_HE_sha256": r.registered_HE_sha256,
            "n_cells": int(r.n_cells), "R1_scale_status": "PASS_SCALE_TAG" if scale_ok else "UNRESOLVED_MISSING_PHYSICAL_UNIT",
            "physical_pixel_size_x_um": r.pixel_size_x_um, "physical_pixel_size_y_um": r.pixel_size_y_um,
            "transform_rule": "author-registered image in Xenium common frame; identity orientation and zero translation; no refit",
            "centroids_inside_image": r2, "R2_threshold": .95, "R2_status": r2_status,
            "tissue_coverage": np.nan, "R3_threshold": .90, "R3_status": "NOT_COMPUTED_FULL_TIFF_PIXEL_PASS_REQUIRED",
            "R4_status": "N_A_NO_WITHHELD_LANDMARKS_SUPPLIED",
            "dapi_rho": np.nan, "R5_status": "UNAVAILABLE_DAPI_MORPHOLOGY_NOT_IN_BOUNDED_INPUTS",
            "visual_pass": "NOT_REVIEWED", "n_visual_patches": 0,
            "R6_status": "NOT_REVIEWED",
            "registration_verdict": verdict,
            "reason": ("R2 direct identity-map bounds fail <0.95" if r2_status == "FAIL" else
                       "R4 landmark file unavailable; R5 DAPI image unavailable; R6 patches not reviewed; R3 remains uncomputed"),
            "pixel_unit_tag": trow.ImageJ_spatial_unit,
            "centroids_inside_image_consensus_diagnostic": r.centroids_inside_image_fraction_consensus_diagnostic,
            "consensus_scale_diagnostic_only_um_per_pixel": .2125 if not scale_ok else np.nan,
        })
    reg = pd.DataFrame(rows).sort_values("sample_id")
    reg.to_csv(BASE / "registration_qc.tsv", sep="\t", index=False, na_rep="NA")

    # Admission registry: section-level technical availability is distinct from cohort admission.
    reg_by = reg.set_index("sample_id")
    regions150 = spatial[spatial.grid_um == 150]
    lineage20 = lineage[lineage.target_cell_threshold_per_150um_region == 20]
    core_map = pd.read_csv(BASE / "raw_registered_he_crosswalk.tsv", sep="\t").set_index("sample_id")
    admissions = []
    for sm in primary.itertuples(index=False):
        rr = reg_by.loc[sm.sample_id]
        st = regions150[regions150.sample_id == sm.sample_id]
        if rr.registration_verdict == "FAIL":
            status = "EXCLUDED_REGISTRATION"
            reason = "Direct no-refit R2 centroid-in-image fraction is 0.041 (<0.95)."
        else:
            status = "UNRESOLVED"
            reason = "No complete registration PASS: R3 not computed; R4/R5 evidence unavailable and R6 not reviewed."
        if core_map.loc[sm.sample_id, "raw_to_registered_core_crosswalk"] == "TMA_LEVEL_ONLY_CORE_COORDINATE_UNRESOLVED":
            reason += " Raw whole-TMA to registered-core mapping is also unresolved; TPS morphology comparison cannot be completed."
        admissions.append({
            "sample_id": sm.sample_id, "donor_id": sm.donor_id, "GSM": sm.GSM,
            "published_primary": "YES_TMA1_4_FROZEN_SUBSET", "same_section": "YES",
            "registered_HE": sm.registered_HE_available, "registration_pass": "NO" if status == "EXCLUDED_REGISTRATION" else "UNRESOLVED",
            "warp_status": "WARP_UNRESOLVED", "expression_pass": "YES",
            "coordinates_pass": "YES_OBJECT_JOIN; H&E mapping remains gated",
            "annotation_pass": "YES_CORRECTED_AUTHOR_RDS", "panel_pass": "YES_COMMON_343",
            "microregion_pass": "YES_XENIUM_FRAME_COMPOSITION_FEASIBLE" if len(st) else "NO",
            "lineage_pass": "YES_COHORT_LEVEL_EPI_FB_MACROPHAGE_COUNTS",
            "PF_status": "PF" if sm.disease_status == "pulmonary_fibrosis" else "control",
            "admission_status": status, "exclusion_reason": reason,
        })
    pd.DataFrame(admissions).to_csv(BASE / "vannan_admission_registry.tsv", sep="\t", index=False)

    hescape = pd.DataFrame([{
        "hescape_panel": "human-lung-healthy-panel", "hescape_sections": 20, "hescape_donors": 19,
        "reported_source_study": "Vannan et al.; GEO GSE250346",
        "panel_or_study_level_source_link": "RESOLVED",
        "hescape_sample_to_GSM_core_donor": "UNRESOLVED",
        "registered_HE_and_cell_object_join": "UNRESOLVED",
        "crosswalk_status": "PARTIAL",
        "access_terms_accepted": "NO", "full_data_downloaded": "NO",
        "basis": "Public paper/supplement identify Vannan/GSE250346 as lung source and report 20 sections/19 donors; exact HESCAPE section IDs are not mapped to GSM/core/donor here.",
    }])
    hescape.to_csv(BASE / "hescape_vannan_metadata_crosswalk.tsv", sep="\t", index=False)

    # Add article-resolved panel composition without inferring per-gene membership in the two subpanels.
    genes = pd.read_csv(BASE / "gene_panel_registry.tsv", sep="\t")
    genes["panel_version"] = "Xenium human lung base panel PD_277 (246) + custom CVEVZD (97); 343 total"
    genes["common_primary_panel"] = "YES"
    genes["primary_sections"] = 28
    genes["base_vs_custom_membership"] = "NOT_INFERRED_PER_GENE; article reports 246 base + 97 custom"
    genes["panel_source"] = "Vannan et al. Nature Genetics 2025 Methods; GEO GSE250346 corrected Seurat object"
    genes.to_csv(BASE / "gene_panel_registry.tsv", sep="\t", index=False)

    # Donor-level composition, avoiding core-as-replicate treatment.
    full_cells = pd.read_csv(BASE / "vannan_primary_cell_metadata.tsv.gz", sep="\t",
                             usecols=["donor_id", "disease_status", "morphoid_coarse_ontology"])
    comp = pd.crosstab([full_cells.donor_id, full_cells.disease_status], full_cells.morphoid_coarse_ontology)
    for col in ("epithelial", "endothelial", "fibroblast", "smooth_muscle_pericyte", "macrophage", "monocyte",
                "neutrophil", "CD4_T", "CD8_T", "B", "plasma", "mast", "other_immune", "other_stromal", "other"):
        if col not in comp.columns:
            comp[col] = 0
    comp = comp.sort_index(axis=1)
    n = comp.sum(axis=1)
    prop = comp.div(n, axis=0)
    donor_comp = comp.add_prefix("n_").join(n.rename("n_cells")).join(prop.add_prefix("p_"))
    donor_comp.reset_index().to_csv(BASE / "donor_coarse_composition.tsv", sep="\t", index=False)
    p = prop.to_numpy(float)
    pairwise_bc = []
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            pairwise_bc.append(float(np.abs(p[i] - p[j]).sum() / 2.0))
    comp_summary = {
        "donors": int(len(prop)), "coarse_categories": int(prop.shape[1]),
        "mean_pairwise_donor_bray_curtis": float(np.mean(pairwise_bc)),
        "median_pairwise_donor_bray_curtis": float(np.median(pairwise_bc)),
        "interpretation": "descriptive donor composition variation; no inference or prediction",
    }
    (BASE / "donor_composition_variation.json").write_text(json.dumps(comp_summary, indent=2) + "\n", encoding="utf-8")

    # PF-only donor/section/region count feasibility, not an analysis model.
    pf = primary[primary.disease_status == "pulmonary_fibrosis"]
    pf_ids = set(pf.sample_id)
    pf_reg_donors = reg[(reg.registration_verdict == "PASS") &
                        (reg.disease_status == "pulmonary_fibrosis")].donor_id.nunique()
    pf_regions = spatial[spatial.sample_id.isin(pf_ids)]
    pf_rows = []
    for grid in (150, 100, 200):
        q = pf_regions[pf_regions.grid_um == grid]
        for target in ("epithelial", "fibroblast", "macrophage", "CD4_T", "CD8_T"):
            col = f"n_{target}"
            eligible_regions = q[q[col] >= 20]
            pf_rows.append({
                "grid_um": grid, "PF_sections": len(pf_ids), "PF_donors": len(pf.donor_id.unique()),
                "PF_regions": len(q), "median_cells_per_region": float(q.total_cells.median()),
                "target_lineage": target, "regions_with_at_least_20_target_cells": len(eligible_regions),
                "eligible_donors": eligible_regions.donor_id.nunique(),
                "registration_pass_donors": pf_reg_donors,
                "program_coverage_note": str(lineage20[lineage20.lineage == target].program_coverage_status.iloc[0]) if (lineage20.lineage == target).any() else "No prespecified state program",
                "PF_only_admission": "NO_REGISTRATION_GATE_UNRESOLVED",
            })
    pd.DataFrame(pf_rows).to_csv(BASE / "pf_only_feasibility.tsv", sep="\t", index=False)

    object_audit = {"n_primary_cells": len(cells)}
    text = markdown_report(cross, tiff, bounds, reg, object_audit, spatial, lineage, coverage, conf,
                           pd.read_csv(BASE / "raw_registered_he_crosswalk.tsv", sep="\t"), "PARTIAL")
    (REPORTS / "PHASE0C_V_VANNAN_REPLACEMENT_AUDIT.md").write_text(text, encoding="utf-8")
    provenance = """# Phase 0C-V image provenance and registration record

## Provenance verdict

`SAME_SECTION_CONFIRMED = YES` for the author-registered H&E workflow described in the 2025 Vannan et al. paper. H&E images used for registration were acquired after the Xenium run by post-run staining of the same 5-µm tissue slice. The paper distinguishes these from an earlier H&E review used to select/grade TMA cores.

## Acquisition order and registration

After Xenium imaging, the same slides were removed from the analyzer, stained with H&E and imaged with a 20× Leica Biosystems Aperio CS2. The authors registered post-run H&E as the moving image to Xenium DAPI morphology as the target using ImageJ 2.14.0 BigWarp. The published workflow used approximately 200 matched landmarks per image pair, a thin-plate-spline warp and manual visual review. The destination coordinate frame is the Xenium DAPI morphology image. No registration was refit here.

## Coordinate system and pixel scale

The Seurat object contains original cell-centroid metadata in microns. Xenium documentation defines cell centroids in microns, 0.2125 μm per Xenium pixel and an upper-left origin. The Vannan Methods also state that registered-H&E QuPath annotations were scaled by 0.2125 μm per pixel before assignment to cell centroids. Among the 28 registered TIFFs audited, 26 explicitly encode ImageJ `unit=micron` and 0.2125 μm/pixel; VUILD104MA2 and VUILD105MA1 encode `unit=pixel` and no usable physical pixel scale. The TIFF metadata table preserves the observed tags and per-file SHA256.

The sample-level mapping used for the diagnostic bounds check was a zero-translation, zero-rotation identity in the published common frame, with scale read from each TIFF. No scale or orientation search and no nonlinear refit were used. On this strict map, 25 of 26 scale-tagged sections meet R2 ≥0.95, VUILD48LA1 measures 0.041 and fails, and two scale-missing sections remain unresolved. The 0.2125 scale applied to those two only for a separately labeled diagnostic; it cannot pass R1/R2.

## Frozen registration criteria audit

| Gate | Result |
|---|---|
| R1 scale/provenance | 26 scale tags present; 2 missing physical units. Published common-frame transform exists, but per-image transform/landmark files are not in the bounded local inputs. |
| R2 image bounds ≥0.95 | 25 PASS, 1 FAIL, 2 UNRESOLVED under the no-refit direct map. |
| R3 tissue mask ≥0.90 | Not computed. Full single-strip TIFF payload totals 32.63 GB uncompressed; R3 is not inferred from bounds. |
| R4 held-out landmarks | N/A; no shared withheld-landmark table supplied. |
| R5 DAPI-density concordance | Unavailable; Xenium morphology OME-TIFFs are not among the selected bounded files. |
| R6 visual patches | Not reviewed; zero sections receive a patch-review PASS. |

With R4 landmarks absent, the frozen policy requires R5 and R6. R5 is unavailable and R6 was not reviewed; therefore no section can be called registration PASS. The observed R2 fractions are diagnostics, not admitted alignment labels.

## Raw-to-registered TMA warp audit

Raw H&E whole-slide files for TMAs 1–5 are present at TMA level. A deterministic sample/core polygon-to-registered-TIFF crosswalk was not found in the public inputs audited, so raw-to-registered morphology comparison is unresolved. Nuclear area/aspect ratio, nearest-neighbor spacing, local texture, architecture and TPS deformation were not compared. Overall `TPS_WARP_STATUS = UNRESOLVED`.

## Sources

- [Vannan et al., Nature Genetics 2025](https://www.nature.com/articles/s41588-025-02080-x)
- [GEO GSE250346](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346)
- [10x Xenium file-format documentation](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html)
"""
    (REPORTS / "PHASE0C_V_IMAGE_PROVENANCE.md").write_text(provenance, encoding="utf-8")
    print("Wrote final admission, crosswalk, registration, PF-only, provenance and audit-report outputs")


if __name__ == "__main__":
    main()
