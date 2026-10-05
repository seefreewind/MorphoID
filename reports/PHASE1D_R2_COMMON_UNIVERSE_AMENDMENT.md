1. **VERDICT:** PHASE1D_R2_COMMON_UNIVERSE_PASS
2. **Protocol amendment rule:** ZERO_VALID_IMAGE_PIXELS_AFTER_FROZEN_BLACK_FILTER
3. **Outcome-blind rule:** YES
4. **Phase 0D total regions:** 12,260
5. **Coordinate-support exclusions:** 77
6. **Zero-valid-image exclusions:** 39
7. **Overlap between exclusion classes:** 0
8. **Unique technical exclusions:** 116
9. **Final common universe:** 12144
10. **Retained sections:** 26
11. **Retained donors:** 19
12. **Original D5 modified:** NO
13. **Imputation used:** NO
14. **Original 39 independently reproduced:** YES
15. **Epithelial injury eligible regions/donors:** 3219/19
16. **Fibroblast activation eligible regions/donors:** 4869/19
17. **Macrophage inflammatory eligible regions/donors:** 2661/19
18. **Color-only TMA macro-AUROC:** 0.7897550366300365
19. **Color-only disease BA/AUROC:** 0.9391025641025641/0.9551282051282052
20. **HIGH shortcut conclusion preserved:** YES
21. **Scientific parameters unchanged:** YES
22. **Common-universe SHA256:** 0d944e4cae613f5e2353e9ebd970d119dda7c39739c1a567405c7ff77e66340e
23. **Phase 1E–1G authorized:** YES
24. **Recommended next action:** STOP; wait for explicit human launch of Phase1E-G

# Phase 1D-R2 protocol amendment

这是输入完整性审计后、任何Phase1预测结果出现前加入的技术eligibility amendment；不声称该排除在原Phase0D方案中已存在。历史PHASE1D_R=TECHNICAL_HOLD_D5_REPAIR保留且报告字节未改。新规则不生成或填补颜色值，original D5的39行NaN保持原样。

独立重算使用全部12260个区域的原D5 whole-slide stride=8采样、0.2125μm/px、150μm网格、float32/dtype最大值归一化与max(RGB)<=0.02黑色过滤。没有把局部crop重置采样原点。只将计数=0作为图像输入无效，没有新增第二阈值。全量计数表与39个crop证据保留。

final_common_phase1_region_universe.tsv仅含common_universe=TRUE的保留区域；全12260行资格记录位于all_region_common_eligibility.tsv。共同集使用C、N、Z与图像输入完整性，不使用lineage门槛或target score，之后再按冻结20-cell和非缺失目标规则计算target资格。所有可直接比较的model arms必须使用同一membership。Primary raw morphology embedding决定共同集；grayscale/Macenko缓存按原Phase1D目标合格union提取，因此非目标合格区域可有空的sensitivity向量，不能据此在target-blind共同集中预删区域。三个目标合格行已逐一验证所有sensitivity向量可用。当前模型脚本尚未接入该新manifest，未来人工启动时须先校验hash并应用membership，不能直接按旧脚本运行。

| target                  |   eligible_regions |   eligible_donors |   eligible_sections |   pre_amendment_regions |   pre_amendment_donors |   difference_vs_pre_amendment |   lineage_min_cells |
|:------------------------|-------------------:|------------------:|--------------------:|------------------------:|-----------------------:|------------------------------:|--------------------:|
| epithelial injury       |               3219 |                19 |                  26 |                    3221 |                     19 |                            -2 |                  20 |
| fibroblast activation   |               4869 |                19 |                  26 |                    4879 |                     19 |                           -10 |                  20 |
| macrophage inflammatory |               2661 |                19 |                  26 |                    2664 |                     19 |                            -3 |                  20 |

## 排除结构

### disease

| disease            |   excluded |   retained |   total |   excluded_fraction |   share_of_all_excluded |
|:-------------------|-----------:|-----------:|--------:|--------------------:|------------------------:|
| control            |         25 |       2367 |    2392 |          0.0104515  |                0.215517 |
| pulmonary_fibrosis |         91 |       9777 |    9868 |          0.00922173 |                0.784483 |

### TMA

|   TMA |   excluded |   retained |   total |   excluded_fraction |   share_of_all_excluded |
|------:|-----------:|-----------:|--------:|--------------------:|------------------------:|
|     1 |         15 |       3039 |    3054 |          0.00491159 |               0.12931   |
|     2 |         52 |       2223 |    2275 |          0.0228571  |               0.448276  |
|     3 |         11 |       3856 |    3867 |          0.00284458 |               0.0948276 |
|     4 |         38 |       3026 |    3064 |          0.0124021  |               0.327586  |

### donor_id

| donor_id   |   excluded |   retained |   total |   excluded_fraction |   share_of_all_excluded |
|:-----------|-----------:|-----------:|--------:|--------------------:|------------------------:|
| THD0008    |          0 |        650 |     650 |          0          |              0          |
| THD0011    |          2 |        393 |     395 |          0.00506329 |              0.0172414  |
| TILD117    |          4 |        912 |     916 |          0.00436681 |              0.0344828  |
| TILD175    |          2 |        384 |     386 |          0.00518135 |              0.0172414  |
| VUHD069    |          6 |        303 |     309 |          0.0194175  |              0.0517241  |
| VUHD095    |          5 |        136 |     141 |          0.035461   |              0.0431034  |
| VUHD113    |          0 |        361 |     361 |          0          |              0          |
| VUHD116    |         12 |        524 |     536 |          0.0223881  |              0.103448   |
| VUILD102   |          0 |        956 |     956 |          0          |              0          |
| VUILD104   |         14 |        773 |     787 |          0.0177891  |              0.12069    |
| VUILD105   |          1 |        305 |     306 |          0.00326797 |              0.00862069 |
| VUILD106   |          1 |       1185 |    1186 |          0.00084317 |              0.00862069 |
| VUILD107   |          0 |        531 |     531 |          0          |              0          |
| VUILD110   |          5 |       1078 |    1083 |          0.00461681 |              0.0431034  |
| VUILD115   |          5 |        943 |     948 |          0.00527426 |              0.0431034  |
| VUILD48    |         26 |        345 |     371 |          0.0700809  |              0.224138   |
| VUILD78    |         16 |        692 |     708 |          0.0225989  |              0.137931   |
| VUILD91    |         14 |        645 |     659 |          0.0212443  |              0.12069    |
| VUILD96    |          3 |       1028 |    1031 |          0.0029098  |              0.0258621  |

### sample_id

| sample_id   |   excluded |   retained |   total |   excluded_fraction |   share_of_all_excluded |
|:------------|-----------:|-----------:|--------:|--------------------:|------------------------:|
| THD0008     |          0 |        650 |     650 |          0          |              0          |
| THD0011     |          2 |        393 |     395 |          0.00506329 |              0.0172414  |
| TILD117LA   |          4 |        425 |     429 |          0.00932401 |              0.0344828  |
| TILD117MA1  |          0 |        487 |     487 |          0          |              0          |
| TILD175MA   |          2 |        384 |     386 |          0.00518135 |              0.0172414  |
| VUHD069     |          6 |        303 |     309 |          0.0194175  |              0.0517241  |
| VUHD095     |          5 |        136 |     141 |          0.035461   |              0.0431034  |
| VUHD113     |          0 |        361 |     361 |          0          |              0          |
| VUHD116A    |          0 |        225 |     225 |          0          |              0          |
| VUHD116B    |         12 |        299 |     311 |          0.0385852  |              0.103448   |
| VUILD102LA  |          0 |        428 |     428 |          0          |              0          |
| VUILD102MA  |          0 |        528 |     528 |          0          |              0          |
| VUILD104MA1 |         14 |        335 |     349 |          0.0401146  |              0.12069    |
| VUILD104MA2 |          0 |        438 |     438 |          0          |              0          |
| VUILD105MA2 |          1 |        305 |     306 |          0.00326797 |              0.00862069 |
| VUILD106MA  |          1 |       1185 |    1186 |          0.00084317 |              0.00862069 |
| VUILD107MA  |          0 |        531 |     531 |          0          |              0          |
| VUILD110LA  |          5 |       1078 |    1083 |          0.00461681 |              0.0431034  |
| VUILD115MA  |          5 |        943 |     948 |          0.00527426 |              0.0431034  |
| VUILD48LA2  |         26 |        345 |     371 |          0.0700809  |              0.224138   |
| VUILD78LA   |          2 |        354 |     356 |          0.00561798 |              0.0172414  |
| VUILD78MA   |         14 |        338 |     352 |          0.0397727  |              0.12069    |
| VUILD91LA   |         14 |        296 |     310 |          0.0451613  |              0.12069    |
| VUILD91MA   |          0 |        349 |     349 |          0          |              0          |
| VUILD96LA   |          3 |        475 |     478 |          0.00627615 |              0.0258621  |
| VUILD96MA   |          0 |        553 |     553 |          0          |              0          |

**EXCLUSION_STRUCTURE_WARNING:** 39个zero-image区域集中于3个section；这是注册H&E图像的局部黑色区域分布。记录结构性集中，不据此改变规则。整体技术排除的donor/disease/TMA计数与比例见上表及TSV；未按任何标签设定eligibility。

## Phase0D regression check

| variant   | target   |   n_regions |   n_sections |   n_donors |   balanced_accuracy |   macro_auroc | risk   |
|:----------|:---------|------------:|-------------:|-----------:|--------------------:|--------------:|:-------|
| original  | TMA      |       12260 |           26 |         19 |            0.564583 |      0.801076 | HIGH   |
| original  | disease  |       12260 |           26 |         19 |            0.939103 |      0.955128 | HIGH   |
| common    | TMA      |       12144 |           26 |         19 |            0.564583 |      0.789755 | HIGH   |
| common    | disease  |       12144 |           26 |         19 |            0.939103 |      0.955128 | HIGH   |

诊断完全调用原section_level_predictions：section平均、donor-held-out frozen folds、训练侧StandardScaler、C=1/max_iter=500/class_weight=balanced logistic、donor-equal权重、1500 donor bootstrap。只重算TMA和disease；未重调参、未运行run预测或新permutation分析。原point metrics通过与历史表核对。此处是授权的Phase0D捷径诊断，不是Phase1 M0–M4或分子预测。

冻结文件、模型脚本、缓存与原HOLD报告before/after SHA相同，见scientific_freeze_assertions.tsv；新增唯一科学eligibility政策如上。报告与新configs manifest保留原HOLD历史。没有计算DeltaR²、residual、PF-only生物模型或图，半小时监测保持取消。

**STOP: PHASE1E_G_AUTHORIZED = YES; WAIT_FOR_HUMAN_LAUNCH.**

## Zero-image排除的分布补充

### donor_id

| donor_id   |   zero_image_excluded |   share_of_zero_image_excluded |
|:-----------|----------------------:|-------------------------------:|
| VUHD095    |                     4 |                       0.102564 |
| VUILD104   |                    14 |                       0.358974 |
| VUILD48    |                    21 |                       0.538462 |

### sample_id

| sample_id   |   zero_image_excluded |   share_of_zero_image_excluded |
|:------------|----------------------:|-------------------------------:|
| VUHD095     |                     4 |                       0.102564 |
| VUILD104MA1 |                    14 |                       0.358974 |
| VUILD48LA2  |                    21 |                       0.538462 |

### disease

| disease            |   zero_image_excluded |   share_of_zero_image_excluded |
|:-------------------|----------------------:|-------------------------------:|
| control            |                     4 |                       0.102564 |
| pulmonary_fibrosis |                    35 |                       0.897436 |

### TMA

|   TMA |   zero_image_excluded |   share_of_zero_image_excluded |
|------:|----------------------:|-------------------------------:|
|     2 |                    39 |                              1 |

39个zero-image排除中，Control为4个，Disease为35个；eligibility判定未读取这些标签。集中性仅作为结构警告记录，不据此更改规则。

最终文件核对：共同集12144行，原39集合完全一致，77/39 overlap=0，19 donors与26 sections保留。没有活动的Phase1模型或embedding重提取任务；未恢复定时监测。
