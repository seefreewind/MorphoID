1. **PHASE 1 VERDICT:** `HOLD_MODEL_INSTABILITY`
2. **Common-universe SHA256:** 0d944e4cae613f5e2353e9ebd970d119dda7c39739c1a567405c7ff77e66340e
3. **Primary morphology encoder:** Phikon-v2
4. **Encoder revision:** 2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd
5. **Primary morphology scale:** 150 μm
6. **Independent donors:** 19
7. **Sections:** 26
8. **Common-universe regions:** 12144
9. **Epithelial injury eligible regions/donors:** 3219 / 19
10. **Fibroblast activation eligible regions/donors:** 4869 / 19
11. **Macrophage inflammatory eligible regions/donors:** 2661 / 19
12. **Median M0 R²:** 0.078
13. **Median M1 R²:** 0.175
14. **Median M2 R²:** 0.155
15. **Median M3 R²:** 0.270
16. **Median M4 R²:** 0.291
17. **Median DeltaR² composition:** 0.009
18. **Median DeltaR² neighborhood:** -0.020
19. **Median DeltaR² morphology:** 0.047
20. **Median cross-fitted residual R²:** -0.119
21. **Targets with positive morphology increment:** epithelial injury; fibroblast activation
22. **Targets with positive residual predictability:** fibroblast activation
23. **PF-only concordance:** 3/3 direction-concordant
24. **Macenko robustness:** median ΔR² = -0.060
25. **Grayscale robustness:** median ΔR² = -0.006
26. **Color-only contribution:** median R² = -0.056
27. **Coordinate-only contribution:** median R² = -0.076
28. **Strongest shortcut:** Phase 0D color→disease coupling (common-universe donor-held-out BA 0.939); among outcome controls: color-only
29. **Donor heterogeneity:** HIGH
30. **Primary FDR-significant morphology targets:** None
31. **Primary FDR-significant residual targets:** None
32. **Scientific interpretation:** 估计的 donor 间方向或不确定性不足以作稳定 GO/NO-GO 判断。
33. **Recommended next phase:** 维持 HOLD。建议下一轮先预注册 donor/fold 稳定性与 M2 高杠杆误差审计，再决定是否启动独立验证；本轮已停止。

## Scope and verdict basis

本报告完成 Phase 1A–1G minimum identifiability pilot，并按 stop rule 在 Phase 1G 后停止。分析只使用冻结 Vannan cohort、冻结 donor folds、冻结 150 μm grid、冻结 within-lineage targets 和一个冻结 encoder。未扩展到其他 cohort、encoder、单细胞/细胞核形态、全转录组或临床结局。

Phase 0D 的 58 项输入门控无不匹配；其中 29 项有历史 sidecar SHA256 并逐字节通过，另外 29 项没有历史 sidecar，记录为 Phase 0D baseline capture 并通过语义/跨输出校验。该差异不应表述成 58 项都与历史 SHA 对照。报告生成前核对到 58 条门控记录。

Verdict operationalization 在 Phase 1 performance outputs 生成前固定于报告脚本。支持性 morphology target 要求 disease-adjusted ΔR²>0、donor-bootstrap CI lower bound>0、三 target BH q<0.05、cross-fitted residual R²及其 CI lower bound均为正、PF-only ΔR²同向为正、Macenko/gray 至少一项增量为正（两项均完整报告），残差三 target BH q<0.05，target donor heterogeneity 非 HIGH，且 raw M4 R²高于 color-only 与 coordinate-only。GO_A 还要求至少两个 target 的 composition 增量及 donor CI 为正、全部 residual CI 跨零、median M4 R²为正、median composition+neighborhood 增量高于 morphology 增量、donor heterogeneity 非 HIGH 且没有 shortcut-dominant target。本 pilot 不自动触发 NO_GO；该判断还需要不可恢复 shortcut 或 estimand 无法解释的独立证据。未达 GO 条件的情形保留 HOLD。

Frozen primary signature coverage: epithelial__injury: 5/5 frozen genes present, PRIMARY; fibroblast__fibroblast_activation: 7/7 frozen genes present, PRIMARY; macrophage__inflammatory: 6/6 frozen genes present, PRIMARY. 每个 primary target 均保持 ≥20 target-lineage-cell 阈值；macrophage inflammatory score 按冻结 ontology 在 macrophage/monocyte lineage 中计算。未按表达结果或模型表现筛选基因或区域。

永久排除的 VUILD105MA1 与 VUILD48LA1 未重新纳入，原因仍为 `FIXED_COORDINATE_SUPPORT_FAILURE`。150 μm crop audit 记录了 77 个全局越界 region；partial crop 未补白或截取。target-specific、边界及低 lineage-count 排除均见 [`region_exclusion_log.tsv`](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/region_exclusion_log.tsv>)。

## Encoder and image preprocessing

使用单一 frozen encoder：Owkin Phikon-v2，DINOv2 ViT-L/16，1024 维 CLS embedding；不微调、不做 target supervision、不按 target 选择 patch。checkpoint hash、版本、许可、预处理、分辨率与 batch 参数记录于 [`morphology_encoder_manifest.tsv`](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/morphology_encoder_manifest.tsv>)。模型使用公开权重，官方模型卡声明 non-commercial license，并说明预训练使用 PANCAN-XL（含 TCGA/CPTAC/GTEx）；本 pilot 没有导入这些额外 cohort。详见 [Phikon-v2 官方模型卡](https://huggingface.co/owkin/phikon-v2)。

每个 crop 保持 150 μm 物理视野，处理器将正方形 crop 调整到 224 × 224，bicubic resampling、center crop、ImageNet mean/SD normalization。RAW 使用作者注册 H&E；GRAYSCALE 使用亮度灰度转换。每张 section 的 source stain summary 用覆盖全图的确定性等距块起点估计：每块连续 100 个 RGB 像素，默认 200 块（共 20,000 点）；若 OD 筛选后不足 500 个染色像素，则扩大到 1,000 块（共 100,000 点）。实际样本数和是否触发 fallback 均记录在 `macenko_section_source_parameters.tsv`。Macenko 在每个 outer fold 只用 training donors 的 section stain summaries 估计目标染色基底：先在 donor 内跨 section 取中位数，再跨 training donors 取中位数；held-out donor 不参与目标参考。目标参考矩阵与逐 section stain 参数均留档。PF-only Ridge 的训练和调参仅使用 PF donors；冻结的 Macenko outer-fold reference 仍由该 fold 全部 training donors 的图像汇总，可能含 control donors。PF-only Macenko 应解释为使用预先冻结、无 outer-test donor 泄漏的参考，而非 PF-only stain-reference 重估。

每个 region 的 tissue_fraction 是每 8 像素采样一次、以 OD sum > 0.15 且排除纯黑像素得到的 proxy，不是经验证的组织分割掩膜。

## Models and validation

M0 = Z；M1 = C + Z；M2 = C + N + Z；M3 = raw X + Z；M4 = C + N + Z + raw X。POOLED 的 frozen Z 使用 normalized x/y；DISEASE_ADJUSTED 在此基础上加入冻结 disease covariate；PF_ONLY 在自身 donors 内独立重新训练；可选 CONTROL_ONLY 未执行。三个 strata 均报告 OOF ΔR² = R²(M4) − R²(M2)，保留负值。颜色输入严格使用 D5 冻结的 40 个 RGB/HED/亮度/饱和度/直方图变量；coordinate-only 使用冻结 Phase 0D 的六项绝对/归一化坐标与图像尺寸项。

所有 outer folds 直接使用冻结 donor folds。每个 alpha 仅在 outer-training donors 内用 GroupKFold 调参，alpha grid 为 log-spaced 10⁻⁴–10⁵；每个 inner split 的 scaler 只在其 training donors 拟合。Ridge 的训练损失及 alpha 选择以 donor-equal MSE 汇总，以免 region 数较多的 donor 独占拟合权重；报告 R²仍按用户定义对全部 OOF regions pooled 计算。

Phase 1G 每个 outer fold 的 M2 training-side residual 由 outer-training donors 内再次 donor-grouped cross-fitting 得到；每个 training region 的 M2 base model 未见过其 donor。residual Ridge 仅用这些 OOF residual 训练，并在 outer-held-out donors 评估。训练侧 residual 明细见 [`phase1g_train_oof_residuals.tsv`](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/phase1g_train_oof_residuals.tsv>)。

不确定性使用 2,000 次 donor-cluster bootstrap；每次整 donor 连同其 sections/regions 一起重采样。ΔR² 的区间直接对同一次 donor resample 计算 R²(M4) − R²(M2)，没有用两个独立区间相减。Primary FDR family 仅是 DISEASE_ADJUSTED 的三个 primary ΔR² morphology；BH 校正基于 donor-bootstrap two-sided sign-tail approximation。未运行 permutation test；因此 q 值按 bootstrap-tail 解释，不是 donor-label permutation p 值。Residual R² 的三个 disease-adjusted 检验组成独立 BH family，采用同样的 bootstrap sign-tail approximation；其他 strata 保留未校正描述性 p 值，不作确认性显著性声明。

## Primary results by target

### epithelial injury

- Pooled M2 R² = 0.180; pooled M4 R² = 0.268; pooled ΔR² morphology = 0.088.
- Disease-adjusted M0–M4 R² = 0.247, 0.256, 0.243, 0.270, 0.291; ΔR² composition = 0.009, neighborhood = -0.012. M2 R² = 0.243; M4 R² = 0.291; ΔR² morphology = 0.047 (95% donor-bootstrap CI -0.109 to 0.204; BH q = 0.519).
- Cross-fitted residual R² = -0.163 (95% donor-bootstrap CI -0.376 to 0.004); residual Pearson r = 0.026; Spearman ρ = 0.066; bootstrap-tail p = 0.062, BH q = 0.093.
- PF-only ΔR² morphology = 0.124; Macenko ΔR² = -0.060; grayscale ΔR² = -0.006.
- Color-only R² = 0.034; coordinate-only R² = -0.038.
- Donor consistency: 13/19 eligible donors had lower MAE under M4; heterogeneity MODERATE (top two positive-error-gain donors account for 41.2% of positive donor error gain).
- Target verdict: `UNRESOLVED_OR_SHORTCUT_SENSITIVE`.

### fibroblast activation

- Pooled M2 R² = -0.233; pooled M4 R² = 0.438; pooled ΔR² morphology = 0.670.
- Disease-adjusted M0–M4 R² = 0.078, -0.155, -0.444, 0.371, 0.438; ΔR² composition = -0.233, neighborhood = -0.289. M2 R² = -0.444; M4 R² = 0.438; ΔR² morphology = 0.882 (95% donor-bootstrap CI 0.006 to 2.221; BH q = 0.067).
- Cross-fitted residual R² = 0.036 (95% donor-bootstrap CI -0.139 to 0.068); residual Pearson r = 0.192; Spearman ρ = 0.261; bootstrap-tail p = 0.470, BH q = 0.470.
- PF-only ΔR² morphology = 0.232; Macenko ΔR² = 0.813; grayscale ΔR² = 0.816.
- Color-only R² = -0.056; coordinate-only R² = -0.234.
- Donor consistency: 11/19 eligible donors had lower MAE under M4; heterogeneity HIGH (top two positive-error-gain donors account for 87.5% of positive donor error gain).
- Target verdict: `UNRESOLVED_OR_SHORTCUT_SENSITIVE`.

### macrophage inflammatory

- Pooled M2 R² = 0.154; pooled M4 R² = -0.081; pooled ΔR² morphology = -0.236.
- Disease-adjusted M0–M4 R² = -0.082, 0.175, 0.155, -0.174, -0.081; ΔR² composition = 0.257, neighborhood = -0.020. M2 R² = 0.155; M4 R² = -0.081; ΔR² morphology = -0.237 (95% donor-bootstrap CI -0.460 to -0.007; BH q = 0.067).
- Cross-fitted residual R² = -0.119 (95% donor-bootstrap CI -0.443 to -0.016); residual Pearson r = -0.089; Spearman ρ = -0.200; bootstrap-tail p = 0.023, BH q = 0.069.
- PF-only ΔR² morphology = -0.258; Macenko ΔR² = -0.121; grayscale ΔR² = -0.303.
- Color-only R² = -0.156; coordinate-only R² = -0.076.
- Donor consistency: 9/19 eligible donors had lower MAE under M4; heterogeneity HIGH (top two positive-error-gain donors account for 95.0% of positive donor error gain).
- Target verdict: `COMPOSITION_DOMINANT`.

## Cross-target tables

| Target | Stratum | Donors | Sections | Regions | Median lineage cells | M0 R² | M1 R² | M2 R² | M3 R² | M4 R² | ΔR² morphology | 95% CI | FDR q |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| epithelial injury | pooled | 19 | 26 | 3219 | 34.0 | -0.044 | 0.172 | 0.180 | 0.251 | 0.268 | 0.088 | -0.066 to 0.248 | NA |
| epithelial injury | PF-only | 13 | 19 | 2474 | 35.0 | -0.078 | -0.052 | -0.045 | 0.050 | 0.079 | 0.124 | -0.075 to 0.376 | NA |
| epithelial injury | disease-adjusted | 19 | 26 | 3219 | 34.0 | 0.247 | 0.256 | 0.243 | 0.270 | 0.291 | 0.047 | -0.109 to 0.204 | 0.519 |
| fibroblast activation | pooled | 19 | 26 | 4869 | 28.0 | -0.158 | 0.179 | -0.233 | 0.368 | 0.438 | 0.670 | -0.011 to 1.692 | NA |
| fibroblast activation | PF-only | 13 | 19 | 4384 | 29.0 | -0.009 | 0.193 | 0.232 | 0.405 | 0.464 | 0.232 | 0.065 to 0.404 | NA |
| fibroblast activation | disease-adjusted | 19 | 26 | 4869 | 28.0 | 0.078 | -0.155 | -0.444 | 0.371 | 0.438 | 0.882 | 0.006 to 2.221 | 0.067 |
| macrophage inflammatory | pooled | 19 | 26 | 2661 | 28.0 | -0.081 | 0.190 | 0.154 | -0.175 | -0.081 | -0.236 | -0.458 to -0.001 | NA |
| macrophage inflammatory | PF-only | 13 | 19 | 2117 | 28.0 | -0.142 | 0.109 | 0.118 | -0.239 | -0.140 | -0.258 | -0.653 to 0.075 | NA |
| macrophage inflammatory | disease-adjusted | 19 | 26 | 2661 | 28.0 | -0.082 | 0.175 | 0.155 | -0.174 | -0.081 | -0.237 | -0.460 to -0.007 | 0.067 |

Full tables: [Table 1 — analysis dataset](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/tables/table1_analysis_dataset.tsv>), [Table 2 — primary model performance](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/tables/table2_primary_model_performance.tsv>), [Table 3 — cross-fitted residual performance](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/tables/table3_crossfit_residual_performance.tsv>), [Table 4 — shortcut stress test](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/tables/table4_shortcut_stress_test.tsv>). All OOF records are in [`oof_predictions.tsv`](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/oof_predictions.tsv>); complete model scores and donor-bootstrap CIs are in [`donor_bootstrap_metrics.tsv`](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/donor_bootstrap_metrics.tsv>).

## Figures

P1, P2 and P4 show completed OOF estimates and 95% percentile intervals from 2,000 donor-cluster resamples. P3 shows all eligible disease-adjusted region pairs as points on symmetric-log axes, linear within ±0.25 and annotates residual R² and its donor-bootstrap interval; the diagonal is identity and both axes retain the full data extent. No tail regions are hidden. P5 shows each of the 19 held-out donors directly, without within-donor bootstrap error bars; marker area follows log(1 + region count). Pooled/disease-adjusted estimates use 19 donors; PF-only uses 13 independently retrained donors. Eligible region counts and all model arms are recorded in the dataset and OOF tables. Zero-reference lines and negative R² are retained. Figures are editable-text SVG plus 300-dpi PNG exports.

### Figure P1. Donor-held-out model ladder

![Figure 1: Model ladder](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/figures/figure1_model_ladder.png>)

### Figure P2. Morphology increment against the M2 oracle

![Figure 2: Morphology increment](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/figures/figure2_morphology_increment.png>)

### Figure P3. Cross-fitted residual prediction

![Figure 3: Residual prediction](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/figures/figure3_crossfit_residual.png>)

### Figure P4. Shortcut stress test

![Figure 4: Shortcut stress test](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/figures/figure4_shortcut_stress_test.png>)

### Figure P5. Donor heterogeneity

![Figure 5: Donor heterogeneity](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/figures/figure5_donor_heterogeneity.png>)

## Donor heterogeneity rule

Donor heterogeneity is classified from DISEASE_ADJUSTED donor-level ΔMAE and concentration of positive squared-error gains. HIGH indicates a near-even split of donors improving vs worsening (40–60%) or the top two donors contributing >50% of positive squared-error gain; MODERATE indicates a 20–40% directional minority or >30% top-two concentration; otherwise LOW. Current overall class: **HIGH**. Target-specific fractions and gains are in [`donor_heterogeneity.tsv`](</Volumes/硬盘2/H&E 图像上的"虚拟空间转录组"到底能看见什么——可识别性上限与捷径学习审计/morphoid/results/phase1/donor_heterogeneity.tsv>).

## Limitations and stop rule

This is a small frozen-cohort pilot with 19 independent donors, and primary target eligibility differs by lineage. Donor bootstrap intervals can remain wide with this donor count. Color and coordinate controls do not fully remove TMA/run effects; those remain a Phase 0D confounding qualification and constrain pooled interpretation. The encoder's non-commercial license also limits downstream use. No causal or clinical validation claim follows from these results.

Secondary programs use a separate family and never alter the primary verdict. Gene-level inference, pseudo-transcriptome controls, 100/200 μm sensitivity, and permutation testing were not used to alter the primary verdict. No gene-level pilot was run because no frozen gene panel was authorized. Phase 1 stops here; extensions listed in the request were not run.

## Reproducibility artifacts

- `scripts/phase1_input_hash_gate.py`, `scripts/phase1a_build_matrix.R`, `scripts/phase1a_crop_support_audit.py`, `scripts/phase1d_fetch_encoder.py`, `scripts/phase1d_extract_morphology.py`, `scripts/phase1e_g_fit_models.py`, `scripts/phase1_figures.py`, `scripts/phase1_report.py`.
- `results/phase1/input_hash_gate.tsv`, encoder manifest, section/fold Macenko parameters, embedding extraction manifest, alpha selections, complete OOF predictions, bootstrap metrics, donor heterogeneity, Phase 1G train-side residuals, exclusions, and Tables 1–4.
- Runtime logs are under `logs/phase1*`; the Phase 1A input gate passed before modeling authorization, and no post-gate frozen input was changed.

## Why the primary verdict remains HOLD

Fibroblast activation has a positive raw disease-adjusted increment, but its M2 baseline R² is negative and its paired interval is broad. The morphology increment fails the three-target BH threshold, and its residual interval crosses zero. Its independently trained PF-only result is stronger and remains supportive descriptive evidence. It does not override the disease-adjusted primary family.

The two donors with the largest positive squared-error gains account for 87.5% of fibroblast positive gains. Epithelial injury has a positive increment with an interval crossing zero; macrophage inflammatory morphology and residual R² are negative. Neither residual null nor a negative result alone is classified as NO_GO. These findings motivate stability review rather than an expanded atlas.

## Frozen secondary programs

These are the three pre-specified secondary-usable within-lineage programs. They were fitted only after all primary targets completed, with the same common universe, folds, models, independent PF-only retraining, and strict residual cross-fit. Disease-adjusted ΔR² and residual tests each use a separate three-program secondary BH family; they do not enter the primary verdict.

| Program | Stratum | Regions | Donors | M2 R² | M4 R² | ΔR² morphology | 95% CI | Secondary q | Residual R² | Residual q |
|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
| epithelial_transitional | POOLED | 3219 | 19 | 0.115 | 0.148 | 0.033 | -0.072 to 0.127 | NA | -0.159 | NA |
| epithelial_transitional | PF_ONLY | 2474 | 13 | -0.020 | 0.048 | 0.068 | -0.057 to 0.235 | NA | -0.018 | NA |
| epithelial_transitional | DISEASE_ADJUSTED | 3219 | 19 | 0.196 | 0.156 | -0.040 | -0.156 to 0.066 | 0.672 | -0.120 | 0.078 |
| fibroblast_myofibroblast | POOLED | 4869 | 19 | 0.085 | 0.207 | 0.122 | -0.148 to 0.405 | NA | -0.004 | NA |
| fibroblast_myofibroblast | PF_ONLY | 4384 | 13 | 0.102 | 0.182 | 0.081 | -0.106 to 0.269 | NA | 0.004 | NA |
| fibroblast_myofibroblast | DISEASE_ADJUSTED | 4869 | 19 | 0.017 | 0.205 | 0.188 | -0.158 to 0.636 | 0.672 | -0.023 | 0.640 |
| macrophage_profibrotic | POOLED | 2661 | 19 | 0.505 | 0.476 | -0.029 | -0.107 to 0.065 | NA | 0.021 | NA |
| macrophage_profibrotic | PF_ONLY | 2117 | 13 | 0.390 | 0.383 | -0.006 | -0.136 to 0.134 | NA | 0.006 | NA |
| macrophage_profibrotic | DISEASE_ADJUSTED | 2661 | 19 | 0.496 | 0.477 | -0.019 | -0.109 to 0.074 | 0.694 | 0.014 | 0.827 |

Full secondary OOF, sensitivity scores, alpha selections, bootstrap intervals, donor heterogeneity and training-side residual records are retained in results/phase1/secondary. The six frozen broad secondary programs were optional and were not fitted in this minimum within-lineage pilot. Their full-region grayscale/Macenko sensitivity coverage is not present in the frozen cache; no new embeddings or arm-specific region exclusions were introduced to extend this analysis.

## Low-coverage targets

The frozen low-coverage programs are described only; no confirmatory fits or claims are made.

| Lineage | Program | Panel/reference genes | Frozen coverage |
|---|---|---:|---|
| epithelial | proliferation | 12/200 | 0.060 |
| fibroblast | ECM remodeling | 34/321 | 0.106 |
| T_cell | IFN / cytotoxic state | 24/200 | 0.120 |
| Broad | proliferation | 12/200 | 0.060 |
| Broad | hypoxia | 15/200 | 0.075 |

## Frozen common universe and audit assertions

All model and sensitivity arms use the prospective R2 common universe: 12,144 of 12,260 regions. The exclusions are 77 unsupported coordinates and 39 zero-valid-pixel regions after the unchanged D5 black filter, with no overlap. Target-specific ≥20 lineage cells and finite frozen score eligibility is applied afterward; no arm-specific filtering or imputation is used. Original D5 and historical HOLD evidence remain intact. The post-modeling hash gate rechecks all 43 scientific input hashes.

```text
FULL_DATA_RESIDUALIZATION_USED = NO
MICROREGION_PRIMARY_BOOTSTRAP_USED = NO
TEST_DONOR_PREPROCESSING_LEAKAGE = NO
DONOR_CROSS_FOLD_LEAKAGE = NO
POSTHOC_TARGET_SELECTION = NO
POSTHOC_SECTION_EXCLUSION = NO
COMMON_UNIVERSE_CHANGED_AFTER_MODELING = NO
```

All opening medians and confirmatory target verdicts refer to DISEASE_ADJUSTED; pooled and independently retrained PF-only results remain fully tabulated. R² CIs condition on the completed OOF fits rather than re-fitting each bootstrap. Macenko reference estimation excludes each outer-test donor; its fixed outer-training reference is reused during inner tuning, so inner validation donors can contribute to this unsupervised stain reference. This does not expose outer-test donors but is a limitation of the frozen stain workflow.

The user-requested named tables and compressed primary OOF table are in results/phase1. M3 and M4 were both refitted for grayscale and Macenko; M3 sensitivity scores appear in primary_model_performance.tsv and shortcut_robustness.tsv. Figure audit: reports/PHASE1_FIGURE_QA.md.
