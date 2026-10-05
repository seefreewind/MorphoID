# MorphoID

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23150000.svg)](https://doi.org/10.5281/zenodo.23150000)

Frozen v1.0.0 archive: https://doi.org/10.5281/zenodo.23150000

MorphoID audits the conditional predictive identifiability of molecular state from registered histology. This frozen research release contains analysis code, derived de-identified regional features, full-precision results and figure source data for the H&E–Xenium pulmonary fibrosis pilot.

## Frozen scientific scope

26 admitted sections, 19 independent donors (6 control, 13 pulmonary fibrosis); 12,144 common eligible microregions before target-specific cell-support gates. The morphology unit is 150 μm. The three targets are epithelial injury, fibroblast activation and macrophage inflammatory programs. Single-cell molecular inference and external validation have not been established.

The model ladder compares M0=Z, M1=C+Z, M2=C+N+Z, M3=X+Z and M4=C+N+Z+X. Donor-held-out prediction, training-side donor-cross-fitted residuals, independently retrained PF-only analysis and stain/color/coordinate controls are preserved. Negative R² and invalid finite-donor inference entries must remain visible. Direction recurrence across dependent partitions does not establish donor-level precision or causal identifiability.

## Repository layout

- scripts/: frozen analysis, inference, reporting and original figure code.
- configs/: cohort, targets, donor folds, estimands, preprocessing and registered seeds.
- results/phase1/: processed matrix, color/crop audits, 26 CPU-provenance NPZ embeddings, OOF predictions and primary/residual/control tables.
- results/phase1h/, phase1i/, phase1j/: stability, omission sensitivity, finite-donor inference and planning/feasibility outputs.
- results/phase0/phase0d/: nuisance diagnostics and original D5 color features.
- results/final/master_result_registry.tsv: integrated frozen estimate registry.
- BIB_FIGURE_TABLE_OPTIMIZATION/: original plotting source and full-precision TSV copies (historical directory name; the repository itself is journal-neutral).
- BIB_PRESUBMISSION_QC/figures_final/: final vector figure files and provenance.
- environment/: original declared environment; observed runtime is results/final/observed_runtime.json.
- manifests/release_files.tsv: original file hashes and byte sizes.
- docs/: release scope, upstream access and reproducibility instructions.

## Upstream data and encoder

Source spatial data are public at GEO GSE250346: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346. Download raw data from the original source under its terms. No raw H&E TIFF/NDPI/CZI, upstream RDS or proprietary image copies are included. The official Owkin Phikon-v2 revision is 2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd; checkpoint SHA256 is 261ae680fa699b3b951597fd57aa19c02ef735805acb104b93af69b36d928569. Obtain weights from https://huggingface.co/owkin/phikon-v2 under the provider's terms; weights and RGB crops are not redistributed.

## Reproduction

Inspect frozen results and hashes before any computation. See docs/REPRODUCIBILITY.md for the distinction between result inspection, model reruns from derived features and full source extraction. This release was assembled without running models, resampling or tests. The original declared environment is historical; the observed analysis runtime is the primary execution receipt. Platform-dependent reproduction has not been tested in this release.

## Authors

Da Lin; Ying Chen; Yue Liu; Yu Zhang. Correspondence: Yu Zhang (zhangyu1@wzhealth.com). Affiliations are supplied in CITATION.cff and .zenodo.json. Author order follows the manuscript.

## Reuse

Original project code is licensed under MIT (LICENSE). Phikon-v2-derived embeddings and result outputs are shared under the Owkin non-commercial research license, including its share-alike conditions (licenses/OWKIN_PHIKON_V2_LICENSE.pdf). These outputs are not granted under CC BY or MIT. Public availability does not remove the upstream restrictions. Raw images and model weights are not redistributed.
