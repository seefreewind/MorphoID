# Reproducibility routes

## Inspect results

Use results/final/master_result_registry.tsv and results/phase1/tables/ for full precision estimates. Figures 1–5 use the corresponding frozen TSVs, not manuscript-rounded values. Fibroblast adjusted residual CI upper endpoint is 0.06847325464224278; nearest-three-decimal display is 0.068. Preserve all negative estimates and invalid BCa flags.

## Refit from derived features

The original scripts retain their repository-relative paths. Keep scripts/, configs/ and results/ in the repository root. The 26 CPU embedding NPZs, region_analysis_matrix.tsv, original D5 features, common universe and input hash receipts are included. Historical hash gates and chronological status checks may also require predecessor manifests/raw inputs; inspect each script's prerequisites before execution. Refit in a separate copy so the archived frozen outputs remain intact.

Observed analysis runtime: Python 3.9.6; NumPy 1.26.4; pandas 2.3.3; SciPy 1.13.1; torch 2.8.0; transformers 4.48.1; scikit-learn 1.6.1; matplotlib 3.9.4; Pillow 11.3.0; R 4.4.3. This receipt is distinct from the original environment.yml. No portable environment or numerical parity was tested during release.

Original execution order: phase1a build matrix/crop audit; phase1d extraction; phase1d_r2 common-universe amendment; phase1e_g_fit_models; phase1h_audit; phase1i_audit; phase1j precision/feasibility. These commands are historical entry points, not instructions to reinterpret the frozen study or perform an unauthorized rescue analysis.

## Source extraction

Acquire GEO GSE250346, original registration/molecular objects and official Phikon-v2 weights. Historical manifests retain upstream hashes and local paths as provenance; those paths need to be mapped to your own source installation. Missing raw images/weights are intentional. Do not assume that aggregate/model-derived tables replace their source assay.

## Figure redraw

Create the original output directories and run BIB_FIGURE_TABLE_OPTIMIZATION/figures_final/render_optimized_figures.py in a separate checkout. It reads frozen TSVs and produces presentation files; it performs no model fit. Final published vector copies are provided directly under BIB_PRESUBMISSION_QC/figures_final/.
