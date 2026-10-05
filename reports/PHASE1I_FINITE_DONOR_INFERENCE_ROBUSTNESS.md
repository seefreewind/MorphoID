1. **PHASE 1I VERDICT:** HOLD_FINITE_DONOR_PRECISION
2. **Independent donors:** 19
3. **Primary targets:** 3
4. **Original Phase 1H verdict:** HOLD_INFERENCE_INSTABILITY
5. **Primary point estimates changed:** NO
6. **Original folds changed:** NO
7. **Original models refit for primary analysis:** NO
8. **Reference donor bootstrap retained:** YES
9. **Delete-1 jackknife completed:** YES
10. **BCa valid:** NO; 5/6 effect CIs pass feasibility
11. **Donor-level sign tests completed:** YES
12. **Epithelial direction consensus:** DeltaR2_morph: DIRECTION_CONSISTENT; Residual_R2: DIRECTION_MIXED
13. **Fibroblast direction consensus:** DeltaR2_morph: DIRECTION_CONSISTENT; Residual_R2: DIRECTION_MIXED
14. **Macrophage direction consensus:** DeltaR2_morph: DIRECTION_CONSISTENT; Residual_R2: DIRECTION_CONSISTENT
15. **Epithelial inference classification:** DIRECTION_UNSTABLE
16. **Fibroblast inference classification:** DIRECTION_UNSTABLE
17. **Macrophage inference classification:** NEGATIVE_BUT_INFERENCE_SENSITIVE
18. **Methods with conflicting direction:** [{"target": "epithelial_injury", "effect": "Residual_R2"}]
19. **Methods with support/uncertain switching:** [{"target": "epithelial_injury", "effect": "Residual_R2"}, {"target": "fibroblast_activation", "effect": "DeltaR2_morph"}, {"target": "macrophage_inflammatory", "effect": "DeltaR2_morph"}, {"target": "macrophage_inflammatory", "effect": "Residual_R2"}]
20. **Largest jackknife influence target:** epithelial_injury/Residual_R2; donor VUILD115; max abs studentized=7.22
21. **Reference CI width:** epithelial_injury/DeltaR2_morph: 0.3129; epithelial_injury/Residual_R2: 0.3802; fibroblast_activation/DeltaR2_morph: 2.215; fibroblast_activation/Residual_R2: 0.2075; macrophage_inflammatory/DeltaR2_morph: 0.4537; macrophage_inflammatory/Residual_R2: 0.4272
22. **Jackknife CI width:** epithelial_injury/DeltaR2_morph: 0.5832; epithelial_injury/Residual_R2: 1.535; fibroblast_activation/DeltaR2_morph: 5.284; fibroblast_activation/Residual_R2: 0.5842; macrophage_inflammatory/DeltaR2_morph: 1.476; macrophage_inflammatory/Residual_R2: 0.7842
23. **BCa CI width:** epithelial_injury/DeltaR2_morph: 0.3161; epithelial_injury/Residual_R2: 0.4215; fibroblast_activation/DeltaR2_morph: 2.89; fibroblast_activation/Residual_R2: NA; macrophage_inflammatory/DeltaR2_morph: 0.4781; macrophage_inflammatory/Residual_R2: 0.3636
24. **Finite-donor precision limitation:** HIGH
25. **External validation authorized:** NO
26. **Recommended next phase:** Separate pre-registered independent-donor precision/transportability design, subject to human authorization. No external analysis performed.

## Scope, frozen quantities and uncertainty layers

The audit uses the existing disease-adjusted primary analysis: 19 donors, 26 admitted sections, 12,144 common regions, and the original three targets. Target-eligible region counts remain 3,219, 4,869 and 2,661. The original five donor folds, primary and training-side cross-fitted residual predictions, Ridge family, alpha grid, 150 µm scale, lineage threshold of 20, C/N/Z and CPU Phikon-v2 embeddings are unchanged. No model was trained in Phase1I. All 19 donor-deletion refits and all 50 POOLED/PF_ONLY partitions remain frozen; only their existing estimates were read.

Three uncertainty layers are reported separately. REFERENCE_BOOTSTRAP is the original 2,000-replicate donor-cluster OOF bootstrap, including paired M2/M4 donor draws and its original three-target FDR. BALANCED_DONOR_BOOTSTRAP uses 10,000 frozen OOF donor-block resamples and remains conditional on the fitted predictions. The refitted delete-1 jackknife measures donor sensitivity of the fitted pipeline. The 50 partition distributions measure PARTITION_STABILITY; their 5–95% ranges are never presented as sampling confidence intervals. No omnibus interval combines these layers.

## Predeclared inference constructions

Seeds and all operational interpretation rules were saved before calculations in [phase1i_inference_seeds.yaml](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/configs/phase1i_inference_seeds.yaml>). Mean and median are both reported; no trimmed mean or optional donor-count precision curve was run. No method or interval was selected according to significance.

For n=19, pseudovalue_i = 19 × theta_full − 18 × theta_without_i. Jackknife bias = 18 × (mean theta_without_i − theta_full). The pseudovalue mean is the bias-corrected estimate; SE is pseudovalue sample SD / sqrt(19). JACKKNIFE_NORMAL_CI uses that estimate ± t_18(0.975) × SE, despite the historical NORMAL label. The original full-cohort estimator remains the primary point estimate. Its approximate two-sided t p-value is reported as an approximation, not an exact finite-sample guarantee.

Standardized influence is the centered pseudovalue divided by its sample SD. Externally studentized influence compares one pseudovalue with the other 18: (pseudo_i − mean_others)/(SD_others × sqrt(19/18)). Absolute delete-1 influence, median, p90, top donor and top-two combined magnitude/share remain descriptive; they do not authorize exclusion.

Balanced resampling shuffles a pool containing exactly 10,000 copies of each donor and reshapes it into 10,000 blocks of 19 draws. Every replicate retains all regions/sections of each drawn donor, including multiplicity; the same donor draw is paired across M2/M4 and residual losses. It samples prediction blocks only, so duplicate IDs never enter model training. Point and block statistics retain the original region-count R² algebra; donor-equal loss summaries are explicitly separate sensitivities.

BCa is INFERENCE_SENSITIVITY only. Its acceleration must describe the same statistic as its bootstrap distribution: frozen OOF-block delete-1 estimates supply acceleration for frozen OOF-block bootstrap. The refitted Phase1H acceleration is also recorded for comparison, but mixing refit acceleration with conditional prediction-block bootstrap would mix uncertainty layers and is not used. BCa is unreliable with undefined acceleration/bias, degeneracy, >5% invalid draws, nonpositive adjustment denominator, an adjusted tail with fewer than ten expected draws, or deterministic four-chunk endpoint deviations exceeding 25% of CI width. These checks were fixed before computation. No BCa p-value was manufactured.

Donor SSE gain is SSE_M2 − SSE_M4; residual SSE gain is SSE_zero − SSE_prediction. Donor MAE gain uses the corresponding mean absolute errors. Every donor has equal weight in the sign test and the mean/median summaries; raw SSE totals naturally also reflect that donor’s number of regions and are not a donor-size normalization. Exact two-sided binomial sign tests exclude exact-zero ties; Wilcoxon sensitivities use exact sign enumeration of average ranks, conditional on magnitudes and the signed-rank symmetry assumption. SSE and MAE sensitivity families are separate.

For each method/effect/loss, BH uses exactly the three original targets when a p-value exists; the original reference FDR is copied unchanged. Bootstrap sign-tail p-values and jackknife t p-values are approximate, and their multiplicity outputs retain that qualification. BCa has CI classification only. No smallest-p combination is used. Agreement labels use valid CI exclusion where available; CI-free sign/Wilcoxon sensitivities use p<0.05 and the paired-loss direction. FDR classification is a separate column, so unadjusted directional interval support is distinguishable from three-target inference.

Direction consensus requires the original full estimate, refit jackknife mean, LOO majority, POOLED partition median, donor-equal SSE mean and SSE median to agree, with ≥80% of LOO and partition signs in that direction. MAE directions are shown separately. Direction consistency is not significance. Opposing point/donor-summary directions are recorded as INFERENCE_DIRECTION_CONFLICT; support/uncertain or FDR switching is INFERENCE_METHOD_SENSITIVE. A valid method’s uncertainty interval spanning zero is not itself evidence of the opposite direction.

ROBUST_POSITIVE/NEGATIVE requires stable direction, at least two valid supported CI constructions and no inferential/FDR switching; a positive target additionally needs compatible supported positive residual evidence. Direction-mixed delta or residual yields DIRECTION_UNSTABLE. ROBUST_NULL requires all valid primary CIs inside ±0.05 and absolute estimates≤0.01, with no conflict; non-significance alone cannot establish a null. Relative precision is CI width / abs(original point estimate), reported NA when abs(point)<0.01. HIGH precision limitation means at least two of six target/effect combinations have reference or jackknife width at least twice the original effect (or a near-zero effect with width>0.1). These descriptive rules do not estimate a required future donor count.

## Stage-level interpretation

All three delta effects satisfy the saved direction-consensus rule. The epithelial residual has a donor-equal mean/median direction conflict, and fibroblast residual signs agree in 15/19 deletions (78.9%), narrowly below the fixed80% requirement; no threshold was rounded or relaxed. Fibroblast and macrophage delta intervals switch between support and uncertainty across conditional bootstrap and refitted jackknife constructions. Jackknife interval widths are approximately0.58,5.28 and1.48 for delta, versus original effects0.047,0.882 and−0.237. The broad intervals and refit donor sensitivity support HOLD_FINITE_DONOR_PRECISION. No two valid primary uncertainty constructions support opposite scientific directions, so the stronger NO_GO criterion is not established. This verdict concerns finite precision in the current cohort and does not prove that every reasonable estimator or future donor cohort must fail.

## Target-specific results

### epithelial injury

| effect | method | estimate | CI_low | CI_high | p_value | FDR | classification | FDR_classification | valid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeltaR2_morph | REFERENCE_BOOTSTRAP | 0.04744 | -0.1092 | 0.2037 | 0.5187 | 0.5187 | UNCERTAIN | UNCERTAIN | True |
| DeltaR2_morph | JACKKNIFE_NORMAL_CI | 0.06857 | -0.223 | 0.3602 | 0.6272 | 0.6272 | UNCERTAIN | UNCERTAIN | True |
| DeltaR2_morph | BALANCED_DONOR_BOOTSTRAP | 0.04744 | -0.1045 | 0.213 | 0.5247 | 0.5247 | UNCERTAIN | UNCERTAIN | True |
| DeltaR2_morph | BCA_INFERENCE_SENSITIVITY | 0.04744 | -0.1153 | 0.2009 | NA | NA | UNCERTAIN | NOT_DEFINED | True |
| Residual_R2 | REFERENCE_BOOTSTRAP | -0.163 | -0.3763 | 0.003854 | 0.06197 | 0.09295 | UNCERTAIN | UNCERTAIN | True |
| Residual_R2 | JACKKNIFE_NORMAL_CI | -1.031 | -1.799 | -0.2637 | 0.01127 | 0.03382 | NEGATIVE_SUPPORT | NEGATIVE_SUPPORT | True |
| Residual_R2 | BALANCED_DONOR_BOOTSTRAP | -0.163 | -0.3831 | -0.001929 | 0.0458 | 0.06869 | NEGATIVE_SUPPORT | UNCERTAIN | True |
| Residual_R2 | BCA_INFERENCE_SENSITIVITY | -0.163 | -0.433 | -0.01149 | NA | NA | NEGATIVE_SUPPORT | NOT_DEFINED | True |

| effect | loss | donors_improved | donors_worsened | ties | exact_sign_p | exact_Wilcoxon_p |
| --- | --- | --- | --- | --- | --- | --- |
| DeltaR2_morph | SSE_GAIN | 13 | 6 | 0 | 0.1671 | 0.3955 |
| DeltaR2_morph | MAE_GAIN | 13 | 6 | 0 | 0.1671 | 0.3736 |
| Residual_R2 | SSE_GAIN | 10 | 9 | 0 | 1 | 0.5678 |
| Residual_R2 | MAE_GAIN | 10 | 9 | 0 | 1 | 0.7983 |

| effect | LOO_same_full_fraction | partition_same_full_fraction | partition_median | donor_equal_SSE_mean | donor_equal_SSE_median | donor_equal_MAE_mean | donor_equal_MAE_median | direction_consensus |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeltaR2_morph | 0.9474 | 0.98 | 0.09297 | 3.514 | 4.14 | 0.02202 | 0.04589 | DIRECTION_CONSISTENT |
| Residual_R2 | 0.9474 | 0.84 | -0.04941 | -9.083 | 0.03656 | -0.01692 | 0.001987 | DIRECTION_MIXED |

| target | target_classification | delta_direction_consensus | residual_direction_consensus | delta_inference | residual_inference | interpretation_flip_donors | external_positive_eligible |
| --- | --- | --- | --- | --- | --- | --- | --- |
| epithelial_injury | DIRECTION_UNSTABLE | DIRECTION_CONSISTENT | DIRECTION_MIXED | INFERENCE_ROBUST | INFERENCE_DIRECTION_CONFLICT | 5 | False |

### fibroblast activation

| effect | method | estimate | CI_low | CI_high | p_value | FDR | classification | FDR_classification | valid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeltaR2_morph | REFERENCE_BOOTSTRAP | 0.8816 | 0.00608 | 2.221 | 0.04498 | 0.06747 | POSITIVE_SUPPORT | UNCERTAIN | True |
| DeltaR2_morph | JACKKNIFE_NORMAL_CI | 2.01 | -0.6325 | 4.652 | 0.1275 | 0.3824 | UNCERTAIN | UNCERTAIN | True |
| DeltaR2_morph | BALANCED_DONOR_BOOTSTRAP | 0.8816 | -0.001928 | 2.226 | 0.05359 | 0.08039 | UNCERTAIN | UNCERTAIN | True |
| DeltaR2_morph | BCA_INFERENCE_SENSITIVITY | 0.8816 | 0.05841 | 2.949 | NA | NA | POSITIVE_SUPPORT | NOT_DEFINED | True |
| Residual_R2 | REFERENCE_BOOTSTRAP | 0.03637 | -0.1391 | 0.06847 | 0.4698 | 0.4698 | UNCERTAIN | UNCERTAIN | True |
| Residual_R2 | JACKKNIFE_NORMAL_CI | 0.1306 | -0.1615 | 0.4228 | 0.3598 | 0.3598 | UNCERTAIN | UNCERTAIN | True |
| Residual_R2 | BALANCED_DONOR_BOOTSTRAP | 0.03637 | -0.1347 | 0.06347 | 0.4464 | 0.4464 | UNCERTAIN | UNCERTAIN | True |
| Residual_R2 | BCA_INFERENCE_SENSITIVITY | 0.03637 | NA | NA | NA | NA | UNCERTAIN | NOT_DEFINED | False |

| effect | loss | donors_improved | donors_worsened | ties | exact_sign_p | exact_Wilcoxon_p |
| --- | --- | --- | --- | --- | --- | --- |
| DeltaR2_morph | SSE_GAIN | 12 | 7 | 0 | 0.3593 | 0.1336 |
| DeltaR2_morph | MAE_GAIN | 11 | 8 | 0 | 0.6476 | 0.953 |
| Residual_R2 | SSE_GAIN | 12 | 7 | 0 | 0.3593 | 0.3955 |
| Residual_R2 | MAE_GAIN | 11 | 8 | 0 | 0.6476 | 0.3321 |

| effect | LOO_same_full_fraction | partition_same_full_fraction | partition_median | donor_equal_SSE_mean | donor_equal_SSE_median | donor_equal_MAE_mean | donor_equal_MAE_median | direction_consensus |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeltaR2_morph | 1 | 1 | 0.6535 | 105.8 | 11.04 | -0.01134 | 0.04714 | DIRECTION_CONSISTENT |
| Residual_R2 | 0.7895 | 0.98 | 0.07923 | 6.304 | 0.3515 | 0.008425 | 0.004885 | DIRECTION_MIXED |

| target | target_classification | delta_direction_consensus | residual_direction_consensus | delta_inference | residual_inference | interpretation_flip_donors | external_positive_eligible |
| --- | --- | --- | --- | --- | --- | --- | --- |
| fibroblast_activation | DIRECTION_UNSTABLE | DIRECTION_CONSISTENT | DIRECTION_MIXED | INFERENCE_METHOD_SENSITIVE | INFERENCE_METHOD_SENSITIVE | 13 | False |

### macrophage inflammatory

| effect | method | estimate | CI_low | CI_high | p_value | FDR | classification | FDR_classification | valid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeltaR2_morph | REFERENCE_BOOTSTRAP | -0.2366 | -0.4603 | -0.006614 | 0.04398 | 0.06747 | NEGATIVE_SUPPORT | UNCERTAIN | True |
| DeltaR2_morph | JACKKNIFE_NORMAL_CI | -0.3357 | -1.074 | 0.4022 | 0.3518 | 0.5277 | UNCERTAIN | UNCERTAIN | True |
| DeltaR2_morph | BALANCED_DONOR_BOOTSTRAP | -0.2366 | -0.4615 | 0.002264 | 0.05239 | 0.08039 | UNCERTAIN | UNCERTAIN | True |
| DeltaR2_morph | BCA_INFERENCE_SENSITIVITY | -0.2366 | -0.5036 | -0.02553 | NA | NA | NEGATIVE_SUPPORT | NOT_DEFINED | True |
| Residual_R2 | REFERENCE_BOOTSTRAP | -0.1195 | -0.4428 | -0.01565 | 0.02299 | 0.06897 | NEGATIVE_SUPPORT | UNCERTAIN | True |
| Residual_R2 | JACKKNIFE_NORMAL_CI | -0.3708 | -0.7629 | 0.02131 | 0.06239 | 0.09358 | UNCERTAIN | UNCERTAIN | True |
| Residual_R2 | BALANCED_DONOR_BOOTSTRAP | -0.1195 | -0.4364 | -0.0126 | 0.0204 | 0.06119 | NEGATIVE_SUPPORT | UNCERTAIN | True |
| Residual_R2 | BCA_INFERENCE_SENSITIVITY | -0.1195 | -0.3688 | -0.005273 | NA | NA | NEGATIVE_SUPPORT | NOT_DEFINED | True |

| effect | loss | donors_improved | donors_worsened | ties | exact_sign_p | exact_Wilcoxon_p |
| --- | --- | --- | --- | --- | --- | --- |
| DeltaR2_morph | SSE_GAIN | 8 | 11 | 0 | 0.6476 | 0.1447 |
| DeltaR2_morph | MAE_GAIN | 9 | 10 | 0 | 1 | 0.1956 |
| Residual_R2 | SSE_GAIN | 5 | 14 | 0 | 0.06357 | 0.1447 |
| Residual_R2 | MAE_GAIN | 5 | 14 | 0 | 0.06357 | 0.04013 |

| effect | LOO_same_full_fraction | partition_same_full_fraction | partition_median | donor_equal_SSE_mean | donor_equal_SSE_median | donor_equal_MAE_mean | donor_equal_MAE_median | direction_consensus |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeltaR2_morph | 1 | 0.92 | -0.1456 | -6.331 | -0.4777 | -0.03926 | -0.006693 | DIRECTION_CONSISTENT |
| Residual_R2 | 1 | 0.88 | -0.03393 | -2.701 | -0.3934 | -0.02598 | -0.01093 | DIRECTION_CONSISTENT |

| target | target_classification | delta_direction_consensus | residual_direction_consensus | delta_inference | residual_inference | interpretation_flip_donors | external_positive_eligible |
| --- | --- | --- | --- | --- | --- | --- | --- |
| macrophage_inflammatory | NEGATIVE_BUT_INFERENCE_SENSITIVE | DIRECTION_CONSISTENT | DIRECTION_CONSISTENT | INFERENCE_METHOD_SENSITIVE | INFERENCE_METHOD_SENSITIVE | 18 | False |

## Table I1 — Delta inference comparison

| target | method | original_estimate | estimate | CI_low | CI_high | CI_width | relative_precision_original | p_value | FDR | classification |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| epithelial_injury | REFERENCE_BOOTSTRAP | 0.04744 | 0.04744 | -0.1092 | 0.2037 | 0.3129 | 6.596 | 0.5187 | 0.5187 | UNCERTAIN |
| epithelial_injury | JACKKNIFE_NORMAL_CI | 0.04744 | 0.06857 | -0.223 | 0.3602 | 0.5832 | 12.29 | 0.6272 | 0.6272 | UNCERTAIN |
| epithelial_injury | BALANCED_DONOR_BOOTSTRAP | 0.04744 | 0.04744 | -0.1045 | 0.213 | 0.3175 | 6.694 | 0.5247 | 0.5247 | UNCERTAIN |
| epithelial_injury | BCA_INFERENCE_SENSITIVITY | 0.04744 | 0.04744 | -0.1153 | 0.2009 | 0.3161 | 6.664 | NA | NA | UNCERTAIN |
| fibroblast_activation | REFERENCE_BOOTSTRAP | 0.8816 | 0.8816 | 0.00608 | 2.221 | 2.215 | 2.512 | 0.04498 | 0.06747 | POSITIVE_SUPPORT |
| fibroblast_activation | JACKKNIFE_NORMAL_CI | 0.8816 | 2.01 | -0.6325 | 4.652 | 5.284 | 5.994 | 0.1275 | 0.3824 | UNCERTAIN |
| fibroblast_activation | BALANCED_DONOR_BOOTSTRAP | 0.8816 | 0.8816 | -0.001928 | 2.226 | 2.228 | 2.527 | 0.05359 | 0.08039 | UNCERTAIN |
| fibroblast_activation | BCA_INFERENCE_SENSITIVITY | 0.8816 | 0.8816 | 0.05841 | 2.949 | 2.89 | 3.279 | NA | NA | POSITIVE_SUPPORT |
| macrophage_inflammatory | REFERENCE_BOOTSTRAP | -0.2366 | -0.2366 | -0.4603 | -0.006614 | 0.4537 | 1.918 | 0.04398 | 0.06747 | NEGATIVE_SUPPORT |
| macrophage_inflammatory | JACKKNIFE_NORMAL_CI | -0.2366 | -0.3357 | -1.074 | 0.4022 | 1.476 | 6.238 | 0.3518 | 0.5277 | UNCERTAIN |
| macrophage_inflammatory | BALANCED_DONOR_BOOTSTRAP | -0.2366 | -0.2366 | -0.4615 | 0.002264 | 0.4638 | 1.96 | 0.05239 | 0.08039 | UNCERTAIN |
| macrophage_inflammatory | BCA_INFERENCE_SENSITIVITY | -0.2366 | -0.2366 | -0.5036 | -0.02553 | 0.4781 | 2.021 | NA | NA | NEGATIVE_SUPPORT |

[tableI1_delta_inference.tsv](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/tableI1_delta_inference.tsv>)

## Table I2 — Residual inference comparison

| target | method | original_estimate | estimate | CI_low | CI_high | CI_width | relative_precision_original | p_value | FDR | classification |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| epithelial_injury | REFERENCE_BOOTSTRAP | -0.163 | -0.163 | -0.3763 | 0.003854 | 0.3802 | 2.332 | 0.06197 | 0.09295 | UNCERTAIN |
| epithelial_injury | JACKKNIFE_NORMAL_CI | -0.163 | -1.031 | -1.799 | -0.2637 | 1.535 | 9.418 | 0.01127 | 0.03382 | NEGATIVE_SUPPORT |
| epithelial_injury | BALANCED_DONOR_BOOTSTRAP | -0.163 | -0.163 | -0.3831 | -0.001929 | 0.3811 | 2.338 | 0.0458 | 0.06869 | NEGATIVE_SUPPORT |
| epithelial_injury | BCA_INFERENCE_SENSITIVITY | -0.163 | -0.163 | -0.433 | -0.01149 | 0.4215 | 2.586 | NA | NA | NEGATIVE_SUPPORT |
| fibroblast_activation | REFERENCE_BOOTSTRAP | 0.03637 | 0.03637 | -0.1391 | 0.06847 | 0.2075 | 5.706 | 0.4698 | 0.4698 | UNCERTAIN |
| fibroblast_activation | JACKKNIFE_NORMAL_CI | 0.03637 | 0.1306 | -0.1615 | 0.4228 | 0.5842 | 16.06 | 0.3598 | 0.3598 | UNCERTAIN |
| fibroblast_activation | BALANCED_DONOR_BOOTSTRAP | 0.03637 | 0.03637 | -0.1347 | 0.06347 | 0.1981 | 5.448 | 0.4464 | 0.4464 | UNCERTAIN |
| fibroblast_activation | BCA_INFERENCE_SENSITIVITY | 0.03637 | 0.03637 | NA | NA | NA | NA | NA | NA | UNCERTAIN |
| macrophage_inflammatory | REFERENCE_BOOTSTRAP | -0.1195 | -0.1195 | -0.4428 | -0.01565 | 0.4272 | 3.576 | 0.02299 | 0.06897 | NEGATIVE_SUPPORT |
| macrophage_inflammatory | JACKKNIFE_NORMAL_CI | -0.1195 | -0.3708 | -0.7629 | 0.02131 | 0.7842 | 6.565 | 0.06239 | 0.09358 | UNCERTAIN |
| macrophage_inflammatory | BALANCED_DONOR_BOOTSTRAP | -0.1195 | -0.1195 | -0.4364 | -0.0126 | 0.4238 | 3.548 | 0.0204 | 0.06119 | NEGATIVE_SUPPORT |
| macrophage_inflammatory | BCA_INFERENCE_SENSITIVITY | -0.1195 | -0.1195 | -0.3688 | -0.005273 | 0.3636 | 3.043 | NA | NA | NEGATIVE_SUPPORT |

[tableI2_residual_inference.tsv](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/tableI2_residual_inference.tsv>)

## Table I3 — Donor-level sign / Wilcoxon sensitivity

| target | effect | loss | analysis_role | n_donors | donors_improved | donors_worsened | ties | exact_sign_p | exact_Wilcoxon_p | Wilcoxon_statistic | mean | median |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| epithelial_injury | DeltaR2_morph | SSE_GAIN | SENSITIVITY | 19 | 13 | 6 | 0 | 0.1671 | 0.3955 | 73 | 3.514 | 4.14 |
| epithelial_injury | DeltaR2_morph | MAE_GAIN | SENSITIVITY | 19 | 13 | 6 | 0 | 0.1671 | 0.3736 | 72 | 0.02202 | 0.04589 |
| epithelial_injury | Residual_R2 | SSE_GAIN | SENSITIVITY | 19 | 10 | 9 | 0 | 1 | 0.5678 | 80 | -9.083 | 0.03656 |
| epithelial_injury | Residual_R2 | MAE_GAIN | SENSITIVITY | 19 | 10 | 9 | 0 | 1 | 0.7983 | 88 | -0.01692 | 0.001987 |
| fibroblast_activation | DeltaR2_morph | SSE_GAIN | SENSITIVITY | 19 | 12 | 7 | 0 | 0.3593 | 0.1336 | 57 | 105.8 | 11.04 |
| fibroblast_activation | DeltaR2_morph | MAE_GAIN | SENSITIVITY | 19 | 11 | 8 | 0 | 0.6476 | 0.953 | 93 | -0.01134 | 0.04714 |
| fibroblast_activation | Residual_R2 | SSE_GAIN | SENSITIVITY | 19 | 12 | 7 | 0 | 0.3593 | 0.3955 | 73 | 6.304 | 0.3515 |
| fibroblast_activation | Residual_R2 | MAE_GAIN | SENSITIVITY | 19 | 11 | 8 | 0 | 0.6476 | 0.3321 | 70 | 0.008425 | 0.004885 |
| macrophage_inflammatory | DeltaR2_morph | SSE_GAIN | SENSITIVITY | 19 | 8 | 11 | 0 | 0.6476 | 0.1447 | 58 | -6.331 | -0.4777 |
| macrophage_inflammatory | DeltaR2_morph | MAE_GAIN | SENSITIVITY | 19 | 9 | 10 | 0 | 1 | 0.1956 | 62 | -0.03926 | -0.006693 |
| macrophage_inflammatory | Residual_R2 | SSE_GAIN | SENSITIVITY | 19 | 5 | 14 | 0 | 0.06357 | 0.1447 | 58 | -2.701 | -0.3934 |
| macrophage_inflammatory | Residual_R2 | MAE_GAIN | SENSITIVITY | 19 | 5 | 14 | 0 | 0.06357 | 0.04013 | 44 | -0.02598 | -0.01093 |

[tableI3_donor_sensitivity.tsv](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/tableI3_donor_sensitivity.tsv>)

## Table I4 — Inference consensus

| target | target_classification | delta_direction_consensus | residual_direction_consensus | delta_inference | residual_inference | interpretation_flip_donors | external_positive_eligible |
| --- | --- | --- | --- | --- | --- | --- | --- |
| epithelial_injury | DIRECTION_UNSTABLE | DIRECTION_CONSISTENT | DIRECTION_MIXED | INFERENCE_ROBUST | INFERENCE_DIRECTION_CONFLICT | 5 | False |
| fibroblast_activation | DIRECTION_UNSTABLE | DIRECTION_CONSISTENT | DIRECTION_MIXED | INFERENCE_METHOD_SENSITIVE | INFERENCE_METHOD_SENSITIVE | 13 | False |
| macrophage_inflammatory | NEGATIVE_BUT_INFERENCE_SENSITIVE | DIRECTION_CONSISTENT | DIRECTION_CONSISTENT | INFERENCE_METHOD_SENSITIVE | INFERENCE_METHOD_SENSITIVE | 18 | False |

| target | effect | direction_consensus | inference_consensus | support_uncertain_switching | FDR_switching | opposing_supported_CI_methods | method_point_direction_conflict | primary_supported_CI_methods | valid_BCa |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| epithelial_injury | DeltaR2_morph | DIRECTION_CONSISTENT | INFERENCE_ROBUST | False | False | False | False | 0 | True |
| epithelial_injury | Residual_R2 | DIRECTION_MIXED | INFERENCE_DIRECTION_CONFLICT | True | True | False | True | 3 | True |
| fibroblast_activation | DeltaR2_morph | DIRECTION_CONSISTENT | INFERENCE_METHOD_SENSITIVE | True | False | False | False | 2 | True |
| fibroblast_activation | Residual_R2 | DIRECTION_MIXED | INFERENCE_METHOD_SENSITIVE | False | False | False | False | 0 | False |
| macrophage_inflammatory | DeltaR2_morph | DIRECTION_CONSISTENT | INFERENCE_METHOD_SENSITIVE | True | False | False | False | 2 | True |
| macrophage_inflammatory | Residual_R2 | DIRECTION_CONSISTENT | INFERENCE_METHOD_SENSITIVE | True | False | False | False | 3 | True |

[tableI4_inference_consensus.tsv](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/tableI4_inference_consensus.tsv>)

## Influence and BCa feasibility

| target | effect | original_estimate | JACKKNIFE_BIAS_CORRECTED_ESTIMATE | jackknife_bias | jackknife_SE | CI_low | CI_high | p_value | df | max_absolute_standardized_influence | max_absolute_studentized_influence | median_absolute_influence | p90_absolute_influence | top_donor | top2_donors | top2_combined_abs_influence | top2_fraction_abs_influence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| epithelial_injury | DeltaR2_morph | 0.04744 | 0.06857 | -0.02114 | 0.1388 | -0.223 | 0.3602 | 0.6272 | 18 | 2.654 | 3.459 | 0.01818 | 0.05473 | VUILD115 | VUILD115;VUILD110 | 0.1548 | 0.3302 |
| epithelial_injury | Residual_R2 | -0.163 | -1.031 | 0.8682 | 0.3653 | -1.799 | -0.2637 | 0.01127 | 18 | 3.586 | 7.219 | 0.03555 | 0.09695 | VUILD115 | VUILD115;TILD175 | 0.4836 | 0.4214 |
| fibroblast_activation | DeltaR2_morph | 0.8816 | 2.01 | -1.128 | 1.258 | -0.6325 | 4.652 | 0.1275 | 18 | 2.373 | 2.895 | 0.09342 | 0.5673 | VUILD106 | VUILD106;THD0008 | 1.512 | 0.3947 |
| fibroblast_activation | Residual_R2 | 0.03637 | 0.1306 | -0.09428 | 0.139 | -0.1615 | 0.4228 | 0.3598 | 18 | 2.212 | 2.616 | 0.0269 | 0.04727 | VUHD095 | VUHD095;VUILD106 | 0.1279 | 0.243 |
| macrophage_inflammatory | DeltaR2_morph | -0.2366 | -0.3357 | 0.09914 | 0.3512 | -1.074 | 0.4022 | 0.3518 | 18 | 2.123 | 2.472 | 0.04814 | 0.1406 | VUILD115 | VUILD115;VUILD106 | 0.3534 | 0.2957 |
| macrophage_inflammatory | Residual_R2 | -0.1195 | -0.3708 | 0.2513 | 0.1866 | -0.7629 | 0.02131 | 0.06239 | 18 | 2.152 | 2.518 | 0.01669 | 0.07652 | VUILD115 | VUILD115;TILD175 | 0.2167 | 0.3545 |

| BCa_valid | BCa_status | BCa_reason | CI_low | CI_high | acceleration | bias_correction_z0 | adjusted_quantile_low | adjusted_quantile_high | chunk_endpoint_max_deviation | invalid_fraction | target | effect | refit_jackknife_acceleration | acceleration_source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| True | BCA_VALID | NA | -0.1153 | 0.2009 | -0.01356 | -0.03284 | 0.01872 | 0.9675 | 0.005639 | 0 | epithelial_injury | DeltaR2_morph | -0.04755 | Frozen OOF-block delete1, matched to conditional bootstrap |
| True | BCA_VALID | NA | -0.433 | -0.01149 | -0.08677 | 0.03309 | 0.01129 | 0.9584 | 0.02918 | 0 | epithelial_injury | Residual_R2 | -0.09515 | Frozen OOF-block delete1, matched to conditional bootstrap |
| True | BCA_VALID | NA | 0.05841 | 2.949 | 0.1475 | 0.0813 | 0.0823 | 0.9987 | 0.1214 | 0 | fibroblast_activation | DeltaR2_morph | 0.03325 | Frozen OOF-block delete1, matched to conditional bootstrap |
| False | BCA_UNRELIABLE | adjusted_tail_under_10_replicates | -0.01531 | 0.1177 | 0.08699 | 0.3951 | 0.163 | 0.9996 | 0.015 | 0 | fibroblast_activation | Residual_R2 | -0.01704 | Frozen OOF-block delete1, matched to conditional bootstrap |
| True | BCA_VALID | NA | -0.5036 | -0.02553 | -0.06297 | -0.01454 | 0.01162 | 0.9572 | 0.01459 | 0 | macrophage_inflammatory | DeltaR2_morph | -0.01839 | Frozen OOF-block delete1, matched to conditional bootstrap |
| True | BCA_VALID | NA | -0.3688 | -0.005273 | -0.07891 | 0.2801 | 0.04881 | 0.9855 | 0.009639 | 0 | macrophage_inflammatory | Residual_R2 | -0.02542 | Frozen OOF-block delete1, matched to conditional bootstrap |

All 114 refit pseudovalues: [jackknife_pseudovalues.tsv](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/jackknife_pseudovalues.tsv>). All standardized/studentized donor influences: [studentized_donor_influence.tsv](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/studentized_donor_influence.tsv>). Donor-equal means/medians: [donor_equal_effect_summary.tsv](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/donor_equal_effect_summary.tsv>). Donor-paired observations: [donor_paired_losses.tsv](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/donor_paired_losses.tsv>).

## Figures

![I1_delta_inference_comparison](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/figures/I1_delta_inference_comparison.png>)

I1 | Delta inference comparison. Panels a–c show the three targets; circles are original full estimates except the bias-corrected jackknife center. Lines are 95% intervals. Reference and balanced/BCa intervals condition on frozen OOF predictions. Jackknife intervals use 19 existing refits, df18. Unreliable BCa endpoints are omitted with an explicit label.

![I2_residual_inference_comparison](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/figures/I2_residual_inference_comparison.png>)

I2 | Residual inference comparison. Layout and methods match I1; outcome is the original training-side cross-fitted residual R². Negative estimates and zero-spanning intervals remain visible.

![I3_jackknife_donor_pseudovalues](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/figures/I3_jackknife_donor_pseudovalues.png>)

I3 | All donor pseudovalues. Top row is delta, bottom row residual R²; columns are targets. Each point represents one of 19 refit deletion pseudovalues. Dashed lines mark the unchanged full estimate, not a null threshold.

![I4_donor_paired_loss_improvement](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/figures/I4_donor_paired_loss_improvement.png>)

I4 | Paired donor losses. Rows show M2-minus-M4 predictive SSE gain and MAE gain, followed by zero-baseline-minus-residual-prediction SSE gain and MAE gain. Each panel contains all19 donor points, one per donor without independent-unit weighting by region counts. Positive means improved prediction; zero is a reference line. Raw paired observations are deterministic conditional on the frozen predictions, so per-donor error bars are not supplied.

![I5_inference_agreement](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/figures/I5_inference_agreement.png>)

I5 | Qualitative inference agreement. Top/bottom show delta/residual evidence. +/blue indicates positive support; −/brown negative support; ?/gray uncertain; NA/white unreliable BCa. Primary CI methods and donor-loss sensitivity tests are labeled separately; this figure never equates their estimands. FDR status is separately reported in the agreement source table.

Donor-index mapping for I3/I4: 1=THD0008; 2=THD0011; 3=TILD117; 4=TILD175; 5=VUHD069; 6=VUHD095; 7=VUHD113; 8=VUHD116; 9=VUILD102; 10=VUILD104; 11=VUILD105; 12=VUILD106; 13=VUILD107; 14=VUILD110; 15=VUILD115; 16=VUILD48; 17=VUILD78; 18=VUILD91; 19=VUILD96.

## Verification and interpretation

| check | pass |
| --- | --- |
| 19_donors_and_3_targets | True |
| 114_refit_pseudovalues | True |
| no_primary_model_refit | True |
| 60000_balanced_effect_estimates | True |
| reference_bootstrap_retained | True |
| only_three_targets | True |
| source_hashes_unchanged | True |
| donor_losses_unit_weight | True |
| epithelial_injuryDeltaR2_morph_jackknife_algebra | True |
| epithelial_injuryResidual_R2_jackknife_algebra | True |
| fibroblast_activationDeltaR2_morph_jackknife_algebra | True |
| fibroblast_activationResidual_R2_jackknife_algebra | True |
| macrophage_inflammatoryDeltaR2_morph_jackknife_algebra | True |
| macrophage_inflammatoryResidual_R2_jackknife_algebra | True |
| epithelial_injuryDeltaR2_morphMAE_GAIN_loss_recomputed | True |
| epithelial_injuryDeltaR2_morphSSE_GAIN_loss_recomputed | True |
| epithelial_injuryResidual_R2MAE_GAIN_loss_recomputed | True |
| epithelial_injuryResidual_R2SSE_GAIN_loss_recomputed | True |
| fibroblast_activationDeltaR2_morphMAE_GAIN_loss_recomputed | True |
| fibroblast_activationDeltaR2_morphSSE_GAIN_loss_recomputed | True |
| fibroblast_activationResidual_R2MAE_GAIN_loss_recomputed | True |
| fibroblast_activationResidual_R2SSE_GAIN_loss_recomputed | True |
| macrophage_inflammatoryDeltaR2_morphMAE_GAIN_loss_recomputed | True |
| macrophage_inflammatoryDeltaR2_morphSSE_GAIN_loss_recomputed | True |
| macrophage_inflammatoryResidual_R2MAE_GAIN_loss_recomputed | True |
| macrophage_inflammatoryResidual_R2SSE_GAIN_loss_recomputed | True |
| BALANCED_DONOR_BOOTSTRAPDeltaR2_morphPRIMARY_R2_separate_BH_family | True |
| BALANCED_DONOR_BOOTSTRAPResidual_R2PRIMARY_R2_separate_BH_family | True |
| BCA_INFERENCE_SENSITIVITYDeltaR2_morphPRIMARY_R2_separate_BH_family | True |
| BCA_INFERENCE_SENSITIVITYResidual_R2PRIMARY_R2_separate_BH_family | True |
| EXACT_SIGN_TESTDeltaR2_morphMAE_GAIN_separate_BH_family | True |
| EXACT_SIGN_TESTDeltaR2_morphSSE_GAIN_separate_BH_family | True |
| EXACT_SIGN_TESTResidual_R2MAE_GAIN_separate_BH_family | True |
| EXACT_SIGN_TESTResidual_R2SSE_GAIN_separate_BH_family | True |
| JACKKNIFE_NORMAL_CIDeltaR2_morphPRIMARY_R2_separate_BH_family | True |
| JACKKNIFE_NORMAL_CIResidual_R2PRIMARY_R2_separate_BH_family | True |
| WILCOXON_SENSITIVITYDeltaR2_morphMAE_GAIN_separate_BH_family | True |
| WILCOXON_SENSITIVITYDeltaR2_morphSSE_GAIN_separate_BH_family | True |
| WILCOXON_SENSITIVITYResidual_R2MAE_GAIN_separate_BH_family | True |
| WILCOXON_SENSITIVITYResidual_R2SSE_GAIN_separate_BH_family | True |
| epithelial_injuryDeltaR2_morph_reference_values_unchanged | True |
| epithelial_injuryResidual_R2_reference_values_unchanged | True |
| fibroblast_activationDeltaR2_morph_reference_values_unchanged | True |
| fibroblast_activationResidual_R2_reference_values_unchanged | True |
| macrophage_inflammatoryDeltaR2_morph_reference_values_unchanged | True |
| macrophage_inflammatoryResidual_R2_reference_values_unchanged | True |
| epithelial_injuryDeltaR2_morphSSE_GAIN_exact_sign_verified | True |
| epithelial_injuryDeltaR2_morphSSE_GAIN_signed_rank_verified | True |
| epithelial_injuryDeltaR2_morphMAE_GAIN_exact_sign_verified | True |
| epithelial_injuryDeltaR2_morphMAE_GAIN_signed_rank_verified | True |
| epithelial_injuryResidual_R2SSE_GAIN_exact_sign_verified | True |
| epithelial_injuryResidual_R2SSE_GAIN_signed_rank_verified | True |
| epithelial_injuryResidual_R2MAE_GAIN_exact_sign_verified | True |
| epithelial_injuryResidual_R2MAE_GAIN_signed_rank_verified | True |
| fibroblast_activationDeltaR2_morphSSE_GAIN_exact_sign_verified | True |
| fibroblast_activationDeltaR2_morphSSE_GAIN_signed_rank_verified | True |
| fibroblast_activationDeltaR2_morphMAE_GAIN_exact_sign_verified | True |
| fibroblast_activationDeltaR2_morphMAE_GAIN_signed_rank_verified | True |
| fibroblast_activationResidual_R2SSE_GAIN_exact_sign_verified | True |
| fibroblast_activationResidual_R2SSE_GAIN_signed_rank_verified | True |
| fibroblast_activationResidual_R2MAE_GAIN_exact_sign_verified | True |
| fibroblast_activationResidual_R2MAE_GAIN_signed_rank_verified | True |
| macrophage_inflammatoryDeltaR2_morphSSE_GAIN_exact_sign_verified | True |
| macrophage_inflammatoryDeltaR2_morphSSE_GAIN_signed_rank_verified | True |
| macrophage_inflammatoryDeltaR2_morphMAE_GAIN_exact_sign_verified | True |
| macrophage_inflammatoryDeltaR2_morphMAE_GAIN_signed_rank_verified | True |
| macrophage_inflammatoryResidual_R2SSE_GAIN_exact_sign_verified | True |
| macrophage_inflammatoryResidual_R2SSE_GAIN_signed_rank_verified | True |
| macrophage_inflammatoryResidual_R2MAE_GAIN_exact_sign_verified | True |
| macrophage_inflammatoryResidual_R2MAE_GAIN_signed_rank_verified | True |
| all19_donor_draws_balanced | True |
| five_editable_svg_and_png | True |
| I1_delta_inference_comparison_export | True |
| I2_residual_inference_comparison_export | True |
| I3_jackknife_donor_pseudovalues_export | True |
| I4_donor_paired_loss_improvement_export | True |
| I5_inference_agreement_export | True |

The reference intervals and p-values were reproduced directly from original seeds, and primary point estimates were independently checked against the original Phase1 performance tables. Pseudovalue/SE/t algebra, all paired donor losses, separate BH families, input hashes and figure exports pass the required evidence gates. All 26 figure panels were visually inspected. Editable SVG text uses a minimum of 6 pt, 183 mm width and 300 dpi PNG previews. [verification.json](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1i/verification.json>) records source-validator exceptions: PDF is not generated under the user’s output policy; the prospective export is SVG/PNG, and actual SVG dimensions supersede the validator’s width-expression heuristic.

This audit cannot establish an irreducible cohort limitation in every future analysis. It can distinguish conditional interval constructions, refitted donor sensitivity and partition variation in the current19 donors. Extremely large jackknife pseudovalues can arise from nonlinear R² denominators and changes in training composition; the jackknife t construction is approximate. Disagreement between its bias-corrected point and conditional OOF effects does not prove a causal biological reversal. The paired-loss sensitivities provide a different view without replacing the primary estimand.

Current stage verdict is **HOLD_FINITE_DONOR_PRECISION**; finite-donor precision limitation is **HIGH**. The six effect-level consensus rows and target-specific classifications above determine this conclusion under the saved operational rules; no effect was rescued by choosing an interval or excluding a donor.

```text
DONORS_EXCLUDED = NO
TARGETS_CHANGED = NO
FOLDS_CHANGED = NO
ENCODER_CHANGED = NO
MODEL_FAMILY_CHANGED = NO
PRIMARY_POINT_ESTIMAND_CHANGED = NO
UNCERTAINTY_METHOD_SELECTED_BY_SIGNIFICANCE = NO
EXTERNAL_DATA_ACCESSED = NO
```

## Stop

**PHASE 1I VERDICT:** HOLD_FINITE_DONOR_PRECISION

**EXTERNAL VALIDATION AUTHORIZED = NO**

Phase1I is complete. No external cohort search, dataset download, external validation, encoder change, target expansion or manuscript expansion was performed. Await an independent human decision.
