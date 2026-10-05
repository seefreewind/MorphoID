1. **PHASE 1J VERDICT:** HOLD_TRANSPORTABILITY
2. **Current independent donors:** 19
3. **Primary cohort changed:** NO
4. **Primary point estimates changed:** NO
5. **Significance used as planning target:** NO
6. **Current epithelial CI width:** DeltaR2_morph: reference 0.312902, jackknife 0.583165, BCa 0.316116; Residual_R2: reference 0.380162, jackknife 1.535035, BCa 0.421494
7. **Current fibroblast CI width:** DeltaR2_morph: reference 2.214818, jackknife 5.284424, BCa 2.890219; Residual_R2: reference 0.207540, jackknife 0.584243, BCa UNRELIABLE / NOT USED
8. **Current macrophage CI width:** DeltaR2_morph: reference 0.453734, jackknife 1.475788, BCa 0.478073; Residual_R2: reference 0.427151, jackknife 0.784201, BCa 0.363556
9. **Added donors for 33% precision gain — epithelial:** 29 at the requested ≥80% probability threshold (total 48); 32 with the conservative Monte Carlo safeguard (total 51)
10. **Added donors for 33% precision gain — fibroblast:** PRECISION_RESCUE_NOT_FEASIBLE_AT_PRACTICAL_N (evaluated through 40 added / 59 total)
11. **Added donors for 33% precision gain — macrophage:** 29 at the requested ≥80% probability threshold (total 48); 30 with the conservative Monte Carlo safeguard (total 49)
12. **Vannan extension new independent donors:** 16
13. **Vannan extension technically admissible donors:** 0 verified under frozen admission gates; 13 have author-registered per-core H&E metadata
14. **Independent external cohorts identified:** 13 study sources screened; 0 admitted independent validation cohorts
15. **Best external candidate:** NONE ADMITTED. Same-section image-access follow-up priority: Columbia CosMx COVID22; suitable only for supportive lung replication unless unchanged estimand transport is established
16. **Same-section registration feasible:** NO — complete frozen admission not established for a new cohort; same-slide acquisition is documented for some candidates
17. **Frozen targets transportable:** NO — no external candidate has cleared all unchanged gene, measurement-unit and lineage gates
18. **Composition oracle transportable:** NO — no candidate has cleared the complete frozen oracle/unit/ontology gate
19. **150 μm morphology transportable:** NO — new-donor numerical coordinate, tissue and full-crop admission not established
20. **Augmentation feasible:** NO — strict admission remains pending; a 13-donor metadata-supported source exists
21. **Independent validation feasible:** NO
22. **External validation authorized:** NO
23. **Augmentation authorized:** NO
24. **Recommended next phase:** Await explicit authorization for a bounded Vannan TMA5 admission audit of the 13 author-registered new donors; original discovery remains frozen. Do not start augmentation or validation models.

# Phase 1J precision and additional-donor feasibility

Audit date: 2026-10-03. Scope: planning and public metadata only. The frozen Phase 1I verdict remains `HOLD_FINITE_DONOR_PRECISION`. This audit neither changes that inference nor creates an external effect estimate.

## Decision

There are plausible additional donor sources, so project closure for lack of any source is premature. Their availability does not establish admission. The closest augmentation source is Vannan TMA5: 16 donors outside the original 19, including 13 with author-level per-core registered H&E. Coordinate validity, tissue coverage, 150 μm crop support and the nuclear transcript unit have not been cleared for those new donors. No independent cohort currently passes the complete frozen estimand gate.

The empirical projections indicate that 13–16 potential extension donors would not meet the conservative 33% width-reduction target for any of the three targets: the minima are 32 and 30 for epithelial and macrophage, and unattained for fibroblast through 40 added donors. This does not mean the extension would add no information; it means the prespecified practical precision target is not established by that source alone.

## Frozen design and planning contract

The primary analysis retains 19 donors, 26 admitted sections, the common universe, ≥20 lineage cells per eligible region, the three signatures, 150 μm morphology, C/N/Z, Phikon-v2, Ridge model family, original donor folds, OOF predictions and all Phase 1H/1I results. New observations would be labeled AUGMENTATION or INDEPENDENT_VALIDATION and could never replace the original discovery results.

The planning contract was written to `configs/phase1j_precision_plan.yaml` before simulation and captured with SHA256 in `results/phase1j/pre_run_contract.json`. Both start and end checks compare the frozen scientific inputs, Phase 1I outputs and reports; the final gate is recorded in `input_hash_gate.tsv`. Model training and candidate image/count processing were not run.

### Current precision and donor dispersion

| target | effect | original_estimate | reference_CI_width | jackknife_CI_width | BCa_CI_width | pseudovalue_SD | top2_influence_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| epithelial_injury | DeltaR2_morph | 0.04744 | 0.3129 | 0.5832 | 0.3161 | 0.605 | 0.3302 |
| epithelial_injury | Residual_R2 | -0.163 | 0.3802 | 1.535 | 0.4215 | 1.592 | 0.4214 |
| fibroblast_activation | DeltaR2_morph | 0.8816 | 2.215 | 5.284 | 2.89 | 5.482 | 0.3947 |
| fibroblast_activation | Residual_R2 | 0.03637 | 0.2075 | 0.5842 | NA | 0.6061 | 0.243 |
| macrophage_inflammatory | DeltaR2_morph | -0.2366 | 0.4537 | 1.476 | 0.4781 | 1.531 | 0.2957 |
| macrophage_inflammatory | Residual_R2 | -0.1195 | 0.4272 | 0.7842 | 0.3636 | 0.8135 | 0.3545 |


Fibroblast has the greatest relative precision need by the preregistered maximum jackknife/reference-width ratio (approximately 2.82 for residual R²); its Delta R² also has the largest absolute width. The original target family is unchanged. The unreliable fibroblast residual BCa interval remains unavailable. The top-two influence fraction and pseudovalue dispersion are descriptive summaries; no donor was removed.

### Heterogeneity-preserving projections

Each hypothetical cohort retains all 19 original donors exactly once and appends 0–40 draws from their empirical joint donor distribution. The same donor identities are used across all targets and both estimands. Repeated historical blocks are a planning approximation to future independent observations; they are not additional real donors and do not increase the actual observed donor count.

For the OOF projection, 1,000 hypothetical augmentation compositions are evaluated at every integer total n=19–59. Each composition receives 500 whole-donor bootstrap draws; donor sums reconstruct the original region-weighted Delta R² and residual R² algebra. These intervals are conditional on frozen predictions and do not include future model-training uncertainty. A 10,000-draw n=19 baseline supplies a multiplicative calibration to the existing frozen reference width. Calibration factors (0.955–1.007) and simulated baseline widths are retained in `current_precision.tsv`; no observed effect is recentered.

The second method appends donor-paired Phase 1I refit jackknife pseudovalues and projects width as 2 × t(0.975,n−1) × sample SD / √n. This exactly reproduces the original jackknife width at n=19. It is a dispersion proxy for refitted-estimator uncertainty, not a new jackknife refit or a prediction of future training stability. The dotted 1/√n scaling is a reference only.

P1, P2 and P3 require width ≤80%, ≤67% and ≤50% of each method’s own frozen current width. Probabilities use all 1,000 outer projections; invalid projections cannot count as successes. Nonpositive target variance invalidates a bootstrap draw and >5% invalid inner draws invalidates its projection. No unreliable projections were observed. The two methods answer different precision questions and their widths are not pooled.

The unmodified requested probability threshold gives a minimum of 29 added donors (48 total) for both epithelial injury and macrophage inflammatory across both effects/methods; fibroblast remains unattained. Method/effect minima are retained in `requested_probability_threshold.tsv`. The stricter preregistered Wilson safeguard gives 32 and 30. These are separate thresholds, not interchangeable definitions. Visible jumps in pseudovalue median widths are retained: a discrete influential-donor multiplicity changes which mode crosses the median as n grows; curves were not smoothed or forced to decrease.

The requested threshold of ≥80% probability for P2 is reported separately by effect and method. The conservative target threshold also requires the Wilson 95% Monte Carlo lower bound ≥0.80 for both effects and both methods, with the criterion retained at every larger evaluated n. This additional safeguard was specified before simulation. It measures numerical planning reliability, not biological certainty. No criterion uses p values, FDR, effect direction or whether an interval crosses zero.

| target | criterion | minimum_total_donors | MINIMUM_ADDITIONAL_DONORS | status |
| --- | --- | --- | --- | --- |
| epithelial_injury | CONSERVATIVE_BOTH_EFFECTS_BOTH_METHODS | 51 | 32 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| epithelial_injury | OOF_BLOCK_EMPIRICAL_DeltaR2_morph | 45 | 26 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| epithelial_injury | REFIT_PSEUDOVALUE_PROXY_DeltaR2_morph | 45 | 26 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| epithelial_injury | OOF_BLOCK_EMPIRICAL_Residual_R2 | 51 | 32 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| epithelial_injury | REFIT_PSEUDOVALUE_PROXY_Residual_R2 | 50 | 31 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| fibroblast_activation | CONSERVATIVE_BOTH_EFFECTS_BOTH_METHODS | NA | NA | PRECISION_RESCUE_NOT_FEASIBLE_AT_PRACTICAL_N |
| fibroblast_activation | OOF_BLOCK_EMPIRICAL_DeltaR2_morph | NA | NA | PRECISION_RESCUE_NOT_FEASIBLE_AT_PRACTICAL_N |
| fibroblast_activation | REFIT_PSEUDOVALUE_PROXY_DeltaR2_morph | 45 | 26 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| fibroblast_activation | OOF_BLOCK_EMPIRICAL_Residual_R2 | 44 | 25 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| fibroblast_activation | REFIT_PSEUDOVALUE_PROXY_Residual_R2 | 43 | 24 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| macrophage_inflammatory | CONSERVATIVE_BOTH_EFFECTS_BOTH_METHODS | 49 | 30 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| macrophage_inflammatory | OOF_BLOCK_EMPIRICAL_DeltaR2_morph | 45 | 26 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| macrophage_inflammatory | REFIT_PSEUDOVALUE_PROXY_DeltaR2_morph | 43 | 24 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| macrophage_inflammatory | OOF_BLOCK_EMPIRICAL_Residual_R2 | 49 | 30 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |
| macrophage_inflammatory | REFIT_PSEUDOVALUE_PROXY_Residual_R2 | 43 | 24 | FEASIBLE_UNDER_EMPIRICAL_PROJECTION |


At 59 total donors, fibroblast Delta R² meets P2 in only 63.4% of OOF projections (Wilson lower bound 60.4%), despite a median relative width of 0.639. The simpler square-root rule would conceal this tail risk. `PRECISION_LIMIT_NOT_PRACTICALLY_RESCUABLE` applies to its conservative family-level rescue under this empirical planning contract, not as proof that all future data could never improve precision. No n>59 was simulated.

## Table J1. Precision projections

Every paired cell is OOF / pseudovalue proxy. Widths are medians across hypothetical cohorts; probabilities are the proportions meeting each reduction target. The complete continuous-n table includes planning percentiles, Wilson bounds and invalid-draw diagnostics in `results/phase1j/precision_projection.tsv`.

| Target | Effect | Total / added | Width OOF / proxy | Relative width | Pr ≤80% | Pr ≤67% | Pr ≤50% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| epithelial_injury | DeltaR2_morph | 19 / 0 | 0.313 / 0.583 | 1.000 / 1.000 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 |
| epithelial_injury | DeltaR2_morph | 24 / 5 | 0.271 / 0.495 | 0.867 / 0.849 | 0.145 / 0.156 | 0.001 / 0.000 | 0.000 / 0.000 |
| epithelial_injury | DeltaR2_morph | 29 / 10 | 0.244 / 0.452 | 0.780 / 0.775 | 0.609 / 0.630 | 0.040 / 0.046 | 0.000 / 0.000 |
| epithelial_injury | DeltaR2_morph | 34 / 15 | 0.224 / 0.414 | 0.717 / 0.710 | 0.898 / 0.888 | 0.250 / 0.290 | 0.000 / 0.000 |
| epithelial_injury | DeltaR2_morph | 39 / 20 | 0.208 / 0.383 | 0.664 / 0.656 | 0.977 / 0.984 | 0.548 / 0.571 | 0.001 / 0.000 |
| epithelial_injury | DeltaR2_morph | 49 / 30 | 0.183 / 0.339 | 0.585 / 0.581 | 1.000 / 1.000 | 0.941 / 0.938 | 0.038 / 0.076 |
| epithelial_injury | DeltaR2_morph | 59 / 40 | 0.166 / 0.307 | 0.530 / 0.527 | 1.000 / 1.000 | 0.997 / 0.998 | 0.247 / 0.300 |
| epithelial_injury | Residual_R2 | 19 / 0 | 0.380 / 1.535 | 1.000 / 1.000 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 |
| epithelial_injury | Residual_R2 | 24 / 5 | 0.321 / 1.241 | 0.845 / 0.808 | 0.324 / 0.337 | 0.008 / 0.000 | 0.000 / 0.000 |
| epithelial_injury | Residual_R2 | 29 / 10 | 0.287 / 1.070 | 0.756 / 0.697 | 0.633 / 0.570 | 0.167 / 0.198 | 0.000 / 0.000 |
| epithelial_injury | Residual_R2 | 34 / 15 | 0.263 / 1.118 | 0.693 / 0.728 | 0.769 / 0.815 | 0.406 / 0.434 | 0.008 / 0.000 |
| epithelial_injury | Residual_R2 | 39 / 20 | 0.245 / 1.003 | 0.646 / 0.653 | 0.874 / 0.919 | 0.568 / 0.685 | 0.028 / 0.009 |
| epithelial_injury | Residual_R2 | 49 / 30 | 0.218 / 0.848 | 0.574 / 0.553 | 0.963 / 0.991 | 0.816 / 0.808 | 0.171 / 0.187 |
| epithelial_injury | Residual_R2 | 59 / 40 | 0.199 / 0.802 | 0.523 / 0.523 | 0.995 / 1.000 | 0.924 / 0.969 | 0.375 / 0.357 |
| fibroblast_activation | DeltaR2_morph | 19 / 0 | 2.215 / 5.284 | 1.000 / 1.000 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 |
| fibroblast_activation | DeltaR2_morph | 24 / 5 | 1.887 / 4.589 | 0.852 / 0.868 | 0.212 / 0.228 | 0.000 / 0.000 | 0.000 / 0.000 |
| fibroblast_activation | DeltaR2_morph | 29 / 10 | 1.717 / 4.096 | 0.775 / 0.775 | 0.543 / 0.637 | 0.084 / 0.088 | 0.000 / 0.000 |
| fibroblast_activation | DeltaR2_morph | 34 / 15 | 1.770 / 3.713 | 0.799 / 0.703 | 0.502 / 0.907 | 0.363 / 0.297 | 0.000 / 0.000 |
| fibroblast_activation | DeltaR2_morph | 39 / 20 | 1.661 / 3.448 | 0.750 / 0.653 | 0.652 / 0.987 | 0.353 / 0.605 | 0.014 / 0.003 |
| fibroblast_activation | DeltaR2_morph | 49 / 30 | 1.500 / 3.068 | 0.677 / 0.580 | 0.842 / 1.000 | 0.476 / 0.944 | 0.157 / 0.087 |
| fibroblast_activation | DeltaR2_morph | 59 / 40 | 1.414 / 2.774 | 0.639 / 0.525 | 0.984 / 1.000 | 0.634 / 0.999 | 0.123 / 0.313 |
| fibroblast_activation | Residual_R2 | 19 / 0 | 0.208 / 0.584 | 1.000 / 1.000 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 |
| fibroblast_activation | Residual_R2 | 24 / 5 | 0.162 / 0.500 | 0.782 / 0.856 | 0.532 / 0.027 | 0.311 / 0.000 | 0.109 / 0.000 |
| fibroblast_activation | Residual_R2 | 29 / 10 | 0.133 / 0.449 | 0.640 / 0.769 | 0.711 / 0.724 | 0.551 / 0.002 | 0.272 / 0.000 |
| fibroblast_activation | Residual_R2 | 34 / 15 | 0.113 / 0.412 | 0.545 / 0.705 | 0.806 / 0.978 | 0.685 / 0.230 | 0.414 / 0.000 |
| fibroblast_activation | Residual_R2 | 39 / 20 | 0.100 / 0.381 | 0.480 / 0.653 | 0.860 / 1.000 | 0.750 / 0.641 | 0.533 / 0.000 |
| fibroblast_activation | Residual_R2 | 49 / 30 | 0.079 / 0.337 | 0.382 / 0.576 | 0.950 / 1.000 | 0.881 / 0.990 | 0.715 / 0.013 |
| fibroblast_activation | Residual_R2 | 59 / 40 | 0.066 / 0.306 | 0.320 / 0.523 | 0.984 / 1.000 | 0.953 / 1.000 | 0.824 / 0.243 |
| macrophage_inflammatory | DeltaR2_morph | 19 / 0 | 0.454 / 1.476 | 1.000 / 1.000 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 |
| macrophage_inflammatory | DeltaR2_morph | 24 / 5 | 0.381 / 1.266 | 0.840 / 0.858 | 0.343 / 0.088 | 0.017 / 0.000 | 0.000 / 0.000 |
| macrophage_inflammatory | DeltaR2_morph | 29 / 10 | 0.343 / 1.141 | 0.756 / 0.773 | 0.678 / 0.693 | 0.175 / 0.017 | 0.001 / 0.000 |
| macrophage_inflammatory | DeltaR2_morph | 34 / 15 | 0.313 / 1.043 | 0.690 / 0.706 | 0.856 / 0.958 | 0.419 / 0.261 | 0.007 / 0.000 |
| macrophage_inflammatory | DeltaR2_morph | 39 / 20 | 0.290 / 0.963 | 0.639 / 0.653 | 0.956 / 0.998 | 0.632 / 0.631 | 0.027 / 0.000 |
| macrophage_inflammatory | DeltaR2_morph | 49 / 30 | 0.257 / 0.854 | 0.566 / 0.579 | 0.998 / 1.000 | 0.910 / 0.983 | 0.171 / 0.032 |
| macrophage_inflammatory | DeltaR2_morph | 59 / 40 | 0.231 / 0.776 | 0.509 / 0.526 | 1.000 / 1.000 | 0.982 / 1.000 | 0.440 / 0.255 |
| macrophage_inflammatory | Residual_R2 | 19 / 0 | 0.427 / 0.784 | 1.000 / 1.000 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 |
| macrophage_inflammatory | Residual_R2 | 24 / 5 | 0.350 / 0.677 | 0.820 / 0.864 | 0.442 / 0.140 | 0.081 / 0.000 | 0.000 / 0.000 |
| macrophage_inflammatory | Residual_R2 | 29 / 10 | 0.309 / 0.608 | 0.722 / 0.775 | 0.648 / 0.660 | 0.376 / 0.033 | 0.008 / 0.000 |
| macrophage_inflammatory | Residual_R2 | 34 / 15 | 0.281 / 0.556 | 0.657 / 0.710 | 0.762 / 0.936 | 0.521 / 0.275 | 0.099 / 0.000 |
| macrophage_inflammatory | Residual_R2 | 39 / 20 | 0.259 / 0.515 | 0.606 / 0.656 | 0.844 / 0.996 | 0.653 / 0.574 | 0.212 / 0.000 |
| macrophage_inflammatory | Residual_R2 | 49 / 30 | 0.225 / 0.456 | 0.527 / 0.582 | 0.943 / 1.000 | 0.831 / 0.974 | 0.414 / 0.054 |
| macrophage_inflammatory | Residual_R2 | 59 / 40 | 0.200 / 0.414 | 0.469 / 0.528 | 0.988 / 1.000 | 0.918 / 1.000 | 0.617 / 0.270 |


## Figure J1. CI width versus hypothetical total donor count

![J1 precision planning](../figures/phase1j/J1_precision_projection.png)

**Legend.** Columns represent the three retained targets; upper and lower rows show Delta R² and residual R² interval widths. Blue curves are calibrated frozen OOF block projections; orange curves are refit-pseudovalue dispersion proxies. Shading is the 10th–90th planning percentile across 1,000 hypothetical augmentation compositions, not a sampling confidence band. Dotted curves show the corresponding 1/√n reference. The original 19 donors remain fixed; appended donors are sampled from their joint empirical distribution. No biological model was fitted and no statistical significance test was used. Source: `precision_projection.tsv`.

## Vannan extension provenance

The current GSE250346 record contains 45 sample records and 35 donor IDs, matching the historical complete crosswalk. The frozen discovery subset was drawn from the original 28 cores/19 donors, of which 26 sections passed admission. The 17 TMA5 records contain 16 new donors plus another core of the original donor TILD117. TMA5 is already part of the final 2025 publication; it is an extension of the frozen preprint-era subset, not a newly discovered post-final-publication cohort. No later additional GSM was found in this audit. [Current GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346), [final author article](https://www.nature.com/articles/s41588-025-02080-x).

The shared author panel is 343 genes and contains all 18 frozen signature genes. TMA5 uses Xenium 2.0.0.10 and multi-tissue staining, versus 1.1.0.2 in the original TMAs. Author annotations for all 17 TMA5 sample IDs already exist locally. The author analysis restricts cell-aware expression to nuclear transcripts; extension admission must preserve that unit instead of silently substituting updated full-cell counts. Public registered per-core H&E is absent for TILD028LA, TILD167LA and VUHD049. Whole-TMA or alternate-modality images are not counted as an admitted per-core replacement.

Historical corrections remain provenance events: the August 2024 codeword correction affects ACTA2/MRC1; the corrected annotation/expression source is retained. The May 2025 TILD117 rename separates original MA1 (GSM7990534) from TMA5 MA2 (GSM8505453); it does not create a donor. The June 2025 processed-data update and current February 2026 GEO update are tracked without changing the primary inputs. Same-section author provenance is supported; numerical registration PASS has not been inferred from filenames.

### Table J2. TMA5 donor audit

| sample_id | donor_id | GSM | new_independent_donor | registered_HE | coordinate_validity | morphology_150um |
| --- | --- | --- | --- | --- | --- | --- |
| TILD028LA | TILD028 | GSM8505448 | YES | NO_PER_CORE_REGISTERED_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| TILD049MA | TILD049 | GSM8505449 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| TILD080LA | TILD080 | GSM8505450 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| TILD111LA | TILD111 | GSM8505451 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| TILD113LA | TILD113 | GSM8505452 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| TILD117MA2 | TILD117 | GSM8505453 | NO | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| TILD130LA | TILD130 | GSM8505454 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| TILD167LA | TILD167 | GSM8505455 | YES | NO_PER_CORE_REGISTERED_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| TILD299MA | TILD299 | GSM8505456 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| TILD315MA | TILD315 | GSM8505457 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| VUHD038 | VUHD038 | GSM8505458 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| VUHD049 | VUHD049 | GSM8505459 | YES | NO_PER_CORE_REGISTERED_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| VUHD090 | VUHD090 | GSM8505460 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| VUILD141MA | VUILD141 | GSM8505461 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| VUILD142MA | VUILD142 | GSM8505462 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| VUILD49LA | VUILD49 | GSM8505463 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |
| VUILD58MA | VUILD58 | GSM8505464 | YES | YES_AUTHOR_FILE | PENDING_FROZEN_ADMISSION_QA | PENDING_TISSUE_CROP_QA |


`VANNAN_EXTENSION_INSUFFICIENT` under the complete technical-admission rule: 16 new donor IDs, 13 metadata-supported registered partners, 0 newly technically admitted. This verdict records unresolved admission evidence, not a finding that 16 donors are numerically below five. The 45-row crosswalk retains all requested fields, correction notes, source URLs and original admission membership.

## Public-cohort search and independence

A bounded public search through 2026-10-03 screened Xenium first, followed by CosMx, MERFISH, Visium HD and conventional Visium. Thirteen independent study sources were examined, alongside Vannan extension, derived repositories and a single-donor platform contingency. Unknown specimen duplication, donor crosswalk, gene coverage or registration is kept unknown. Repeated cores, modalities, visits and panel releases are never counted as independent donors. The registry is a screening record, not a declaration that all listed studies are compatible. No exhaustive-absence claim is made.

The most direct independent same-slide acquisition lead is the 22-donor Columbia CosMx study. Its protocol explicitly stains the imaged slides with H&E after CosMx. Public raw spatial files, registration transforms, exact frozen gene coverage and FOV support remain unverified. Its GEO superseries points to sequence experiments; GSE287111 is scRNA-seq and must not be mislabeled as the CosMx donor dataset. COVID/control is not the frozen PF/control adjustment. [Author spatial study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12463091/).

The COPD custom panel was checked directly from all four small public panel JSON files: only 10 of 18 frozen genes are present. Missing genes are COL1A1, COL1A2, COL3A1, ITGB6, KRT18, S100A8, S100A9, SOX4. All three target definitions lose genes. They are `TARGET_NOT_TRANSPORTABLE`; no subset signature or proxy gene can rescue admission. Gene-level evidence is in `frozen_gene_coverage.tsv`. Linked participant metadata are controlled through BioLINCC; fluorescence morphology files are not H&E. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313006), [author article](https://www.nature.com/articles/s41588-025-02480-z).

### Table J3. Candidate registry

| candidate_id | accession | n_donors | same_section | candidate_role | admission_status | reason | source_URL |
| --- | --- | --- | --- | --- | --- | --- | --- |
| VANNAN_TMA5 | GSE250346 | 16 | YES | AUGMENTATION_CANDIDATE | PENDING | 16 independent discovery-external donors; same study. 13 have per-core registered H&E; 0 newly admitted under frozen numerical gates. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346 |
| CHUV_XENIUM | GSE311609 | 10 | YES_AUTHOR_PROTOCOL | SUPPORTIVE_ONLY | NOT_ADMITTED | 22 panel/section records are 10 donors; 5K has six distinct donors. Existing registration audit admitted none; PF/control disease structure absent. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311609 |
| COPD_XENIUM | GSE313006 | 38 | HE_PAIR_UNKNOWN | SUPPORTIVE_ONLY | NOT_ADMITTED | Four TMA records represent 38 participants, not four donors. All four panels lack eight frozen genes. No signature substitution allowed. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313006 |
| COSMX_COVID22 | LungSpatialDB;GSE287114_SEQUENCE_ONLY | 22 | YES_POST_COSMX_HE_PROTOCOL | SUPPORTIVE_ONLY | PENDING | 22 donors, 400 FOVs. GEO GSE287111 is organoid scRNA-seq, not the 22-donor CosMx matrix. Sparse FOV sampling and PF/control absence block unchanged estimand. | https://pmc.ncbi.nlm.nih.gov/articles/PMC12463091/ |
| MAYR_IPF | Zenodo10012934;10015169 | 7 | XENIUM_ADJACENT_VISIUM_HE | SUPPORTIVE_ONLY | NOT_ADMITTED | Seven donors across Visium; do not infer seven distinct Xenium donors. Adjacent-section pairing fails strict same-section Xenium requirement. | https://pmc.ncbi.nlm.nih.gov/articles/PMC11313858/ |
| PF_VISIUM8 | S-BSST1410 | 8 | VISIUM_SAME_SECTION | SUPPORTIVE_ONLY | NOT_ADMITTED | Four IPF and four control donors. Spot mixtures cannot supply unchanged single-cell within-lineage targets and true composition oracle. | https://www.nature.com/articles/s41588-024-01819-2 |
| ADULT_MERFISH6 | HuBMAP/LungMAP publication a10041ad9ebae0b42d3c7f602ba37b82 | 6 | EXACT_HE_MERFISH_PAIR_UNVERIFIED | SUPPORTIVE_ONLY | PENDING | Six spatial donors; eleven overall atlas donors are not six additional spatial donors. Same-section H&E not established; healthy-only cannot reproduce PF-adjusted comparison. | https://pmc.ncbi.nlm.nih.gov/articles/PMC12140004/ |
| NSCLC_MERFISH46 | SCP2510;Zenodo11198494 | 46 | EXACT_HE_PAIR_UNKNOWN | SUPPORTIVE_ONLY | PENDING | 46 processed patient datasets documented by author deposit. RNA/protein coregistration does not establish a registered H&E partner or PF/control transport. | https://zenodo.org/records/11198494 |
| NEUMAP_HD8 | GSE266680;GSE298732;GSM9021740 | 8 | YES_VISIUM_HD_HE_WORKFLOW | SUPPORTIVE_ONLY | PENDING | Eight patients, three cores per patient, not 24 donors. Core-to-patient mapping and cell-state recovery need verification; PF/control Z absent. | https://www.nature.com/articles/s41586-025-09807-0 |
| EGFR_HD4 | GSE288758 | 4 | YES_VISIUM_HD_WORKFLOW | SUPPORTIVE_ONLY | PENDING | Four named patients ER, SR, LR1, LR2; pre/post specimens are repeated donors. Below preferred five, no frozen PF/control comparison. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288758 |
| ASTHMA_XENIUM | GSE269354 | 8_CONFIRMED_COHORT1;16_SPECIMEN_IDS_UPPER_BOUND | HE_ADDITIONAL_SECTION | SUPPORTIVE_ONLY | NOT_ADMITTED | Metadata lists eight distinct biopsy IDs in each of two historical cohorts; patient cross-cohort duplication cannot be independently ruled out. Sixteen specimen IDs are a donor upper bound, not a certified independent count. | https://www.nature.com/articles/s41590-025-02161-3 |
| AML_XENIUM | GSE319763 | UNKNOWN | UNKNOWN_HE_PAIR | SUPPORTIVE_ONLY | PENDING | Three GSM records do not certify three independent donors. Disease context differs; unknown donor crosswalk and registered H&E block admission. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319763 |
| COSMX_COVID6 | GSE327879 | 6 | UNKNOWN_HE_PAIR | SUPPORTIVE_ONLY | PENDING | GEO design identifies six patients on two slides; no strict H&E pairing/registration or unchanged PF/control Z established. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327879 |
| VENDOR_COSMX3 | Bruker CosMx NSCLC FFPE | 3 | UNKNOWN_HE_PAIR | SUPPORTIVE_ONLY | NOT_ADMITTED | Three patients; repeated tissue samples are not independent donors. Below preferred five and no unchanged disease structure. | https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/ |
| VANNAN_DERIVED | CAMEO-Lung | 19 | DERIVED | SUPPORTIVE_ONLY | EXCLUDED_DUPLICATE_DISCOVERY | Derived current Vannan donors add zero independent information; do not add repository counts across releases. | https://huggingface.co/datasets/theislab/CAMEO-Lung |
| SINGLE_DONOR_DEMO | 10x public lung examples | 1 | VERSION_DEPENDENT | SUPPORTIVE_ONLY | EXCLUDED_SINGLE_DONOR | Single-donor demonstration cannot meet the independent donor gate; no count aggregation across panel demos without donor provenance. | https://www.10xgenomics.com/datasets |


Roles are prospective suitability classes. Candidates missing a required gate have no permission for modeling. SUPPORTIVE_ONLY includes partial tissue/platform replication that cannot apply the complete frozen estimand. No independent study currently qualifies for EXTERNAL_VALIDATION_CANDIDATE admission. Existing CHUV registration failures are reused from the frozen Phase 0C-S evidence, without rerunning transform optimization. Asthma has eight confirmed cohort-1 donors and 16 specimen IDs across historical cohorts; independent cross-cohort duplication is not certified. Its H&E uses additional sections, which fails strict same-section admission.

## Table J4. Frozen-estimand transportability

| candidate_id | target_definitions | C_ontology | N_structure | Z_structure | morphology_150um | encoder_preprocessing | admission_status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| VANNAN_TMA5 | PASS_18_OF_18 | AUTHOR_MAPPING_AVAILABLE_PENDING_UNIT_CHECK | PENDING_COORDINATE_QA | PF_CONTROL_METADATA_AVAILABLE | PENDING_CROP_QA | RGB_HE_AVAILABLE_PENDING_QA | PENDING |
| CHUV_XENIUM | 5K_EXACT_GENES_NOT_VERIFIED;OTHER_PANELS_INCOMPLETE | PENDING | PENDING_REGISTERED_COORDINATES | FAIL_PF_CONTROL_ABSENT | FAIL_OR_UNRESOLVED_REGISTRATION | HE_AVAILABLE_GEOMETRY_UNRESOLVED | NOT_ADMITTED |
| COPD_XENIUM | FAIL_10_OF_18 | PENDING | PENDING | FAIL_COPD_NOT_FROZEN_PF_STATUS | UNKNOWN_HE_SCALE | NO_PUBLIC_REGISTERED_RGB_HE | NOT_ADMITTED |
| COSMX_COVID22 | UNKNOWN_EXACT_18_GENES | PENDING | PENDING_FOV_SAMPLING | FAIL_COVID_NOT_PF | PENDING_FOV_BOUNDARIES | POST_STAIN_HE_EXISTS_PUBLIC_MAPPING_UNKNOWN | PENDING |
| MAYR_IPF | 289_PANEL_EXACT_18_NOT_VERIFIED | PENDING | VISIUM_SPOT_NEIGHBORHOODS_DIFFER | PF_CONTROL_AVAILABLE | FAIL_XENIUM_ADJACENT_HE | VISIUM_HE_SUPPORTIVE_ONLY | NOT_ADMITTED |
| PF_VISIUM8 | WHOLE_TRANSCRIPTOME;EXACT_ASSAY_COVERAGE_UNVERIFIED | FAIL_INFERRED_COMPOSITION | FAIL_SPOT_UNIT | PF_CONTROL_AVAILABLE | PENDING_NUMERICAL_QA | HE_AVAILABLE | NOT_ADMITTED |
| ADULT_MERFISH6 | 503_GENES_EXACT_18_UNKNOWN | PENDING | PENDING | FAIL_HEALTHY_ONLY | UNKNOWN_SAME_HE_SECTION | UNKNOWN_REGISTERED_RGB_HE | PENDING |
| NSCLC_MERFISH46 | 484_GENES_EXACT_18_UNKNOWN | PENDING | PENDING | FAIL_PF_CONTROL_ABSENT | UNKNOWN_HE_PAIR | UNKNOWN_REGISTERED_RGB_HE | PENDING |
| NEUMAP_HD8 | WHOLE_TRANSCRIPTOME_PROBE_ASSAY;EXACT_18_UNVERIFIED | PENDING_BIN2CELL_ONTOLOGY | PENDING_RECOVERED_CELLS | FAIL_PF_CONTROL_ABSENT | AUTHOR_SCALE_SUPPORTED_QA_PENDING | HE_AVAILABLE | PENDING |
| EGFR_HD4 | WHOLE_TRANSCRIPTOME_PROBE_ASSAY;EXACT_18_UNVERIFIED | PENDING | PENDING | FAIL_PF_CONTROL_ABSENT | PENDING | HE_METADATA_AVAILABLE | PENDING |
| ASTHMA_XENIUM | 339_GENES_EXACT_18_UNVERIFIED | PENDING | PENDING | FAIL_ASTHMA_NOT_PF | FAIL_ADJACENT_HE | FAIL_SAME_SECTION_PAIR | NOT_ADMITTED |
| AML_XENIUM | UNKNOWN_EXACT_18 | PENDING | PENDING | FAIL_AML_NOT_PF | UNKNOWN | UNKNOWN | PENDING |
| COSMX_COVID6 | 1000_GENES_EXACT_18_UNKNOWN | PENDING | PENDING | FAIL_COVID_NOT_PF | UNKNOWN | UNKNOWN | PENDING |
| VENDOR_COSMX3 | UNKNOWN_EXACT_18 | PENDING | PENDING | FAIL_PF_CONTROL_ABSENT | UNKNOWN | UNKNOWN | NOT_ADMITTED |
| VANNAN_DERIVED | SAME_PANEL | DERIVED | DERIVED | DERIVED | DERIVED | DERIVED | EXCLUDED_DUPLICATE_DISCOVERY |
| SINGLE_DONOR_DEMO | VERSION_DEPENDENT | PENDING | PENDING | UNKNOWN | PENDING | PENDING | EXCLUDED_SINGLE_DONOR |


“NO” in the opening transport fields means no complete new-donor PASS has been established. Pending evidence is distinguished from biological incompatibility in Table J4. Commercial availability of high-plex genes is not exact gene coverage. Cell labels must map to the frozen coarse C ontology and lineage units; deconvolved spot mixtures or image-inferred labels are not automatically true composition oracles. N requires the original physical neighborhood, Z retains normalized x/y and PF/control status, and H&E support must permit full 150 μm crops with unchanged encoder preprocessing. DAPI registration or RNA/protein panel coregistration cannot substitute for H&E admission.

## Figure J2. Candidate feasibility gate map

![J2 candidate feasibility](../figures/phase1j/J2_candidate_feasibility.png)

**Legend.** Rows represent screened study sources, the Vannan augmentation source and excluded duplicate/single-donor contingencies. Columns summarize independent donor evidence, same H&E section, frozen genes, C/N, frozen Z, 150 μm admission and overall admission. “+” denotes metadata support, “?” unresolved evidence, and “x” incompatibility or exclusion. Metadata support does not establish a numerical registration PASS. No biological effects or predicted outcomes enter the map. Asthma’s eight confirmed cohort-1 donors and 16 historical specimen IDs are distinguished. Source: `J2_figure_source.tsv` and the full candidate registry.

## Next-stage options — defined, not executed

**Option A: independent validation.** Admission would require ≥5 independently identified donors from an external cohort, unchanged signatures and measurement units, supported C/N/Z, same-section H&E registration PASS and valid 150 μm morphology. Keep the same encoder, scale, Ridge family and donor-held-out logic. Define an external cohort’s own frozen donor split without rewriting the original split. Fit/evaluate within that separate cohort only after explicit authorization; report positive, null and negative effects equally. No candidate is selected for its anticipated fibroblast result.

**Option B: augmentation.** If the 13 author-registered TMA5 donors clear admission, retain the original 19 as the primary discovery cohort and identify new donors as a prospective augmentation set. Prespecify: original discovery fixed, extension-only estimation, then combined sensitivity. Any combined output remains a sensitivity result and cannot replace the frozen discovery result. Annotation/unit/geometry admission comes before model execution. The 16-donor count cannot be treated as an assured solution to the 30–32 donor precision needs.

**Option C: no further compatible data.** If the bounded admission audit cannot establish new compatible donors, retain FINITE_DONOR_PRECISION_LIMITATION and move to manuscript interpretation only with authorization. Closure is not triggered here because a plausible augmentation source remains. No larger simulation, new target or encoder comparison is proposed.

### Future cohort-level synthesis

Only after independently authorized compatible cohorts exist: estimate Delta R² and residual R² separately per cohort, then perform fixed-effect descriptive synthesis with explicit estimator and variance compatibility. Add random-effects sensitivity only with ≥3 cohorts; with two, do not give heterogeneity variance a strong interpretation. Cohort directions must all be reported. No meta-analysis, concatenation, batch correction or embedding harmonization was performed in Phase 1J.

## Evidence, QA and limitations

Source metadata are retained under `results/phase1j/source_metadata/` with URL, bytes, SHA256 and response status in `source_metadata_manifest.tsv`. These are article XML, GEO descriptive/platform records, small panel files and one participant/specimen metadata spreadsheet, not full expression matrices or H&E datasets. Large GEO platform SOFT records are metadata; their size does not imply a downloaded candidate count matrix. Unsupported GEO `targ=samples` responses are explicitly flagged and not used as evidence. Public portal access failure is an access limitation, not proof that no files exist.

The projections assume future donor contributions resemble the present 19-donor empirical distribution. They cannot identify unseen disease/site heterogeneity, future refit instability, measurement shifts or the information gain of genuinely new donors outside that support. OOF and pseudovalue methods preserve the observed heterogeneity in different ways and must not be read as guaranteed future CIs. The strict target minimum combines two effects and two methods and is more conservative than a single-method threshold.

Figure QA uses the nature-figure Python workflow: source preflight, editable SVG, 300-dpi PNG and complete six-panel/matrix visual inspection at the declared 183-mm size. No PDF was requested. The raw static validator remains NOT READY under its default submission bundle: it requires PDF/TIFF and its 600-dpi default, and misreads the 183/25.4 width expression. These are documented planning-export exceptions in figure_qa.json; actual editable SVG width (183 mm), raster size and 300-dpi resolution are independently verified. Visual and numerical checks are separate. The final completion record is written only after visual inspection and the final frozen-input hash gate.

## Required assertions

```text
PRIMARY_19_DONOR_ANALYSIS_CHANGED = NO
DONORS_REMOVED = NO
TARGETS_CHANGED = NO
ENCODER_CHANGED = NO
MODEL_CHANGED = NO
SIGNIFICANCE_TARGETED = NO
EXTERNAL_MODELING_RUN = NO
AUGMENTATION_MODELING_RUN = NO
EXTERNAL_VALIDATION_AUTHORIZED = NO
AUGMENTATION_AUTHORIZED = NO
```

Phase 1J stops here. No candidate full dataset download, new embedding extraction, predictive refit, external Delta R² calculation or biological signal validation is authorized by this report.
