1. **VERDICT:** TECHNICAL_HOLD_D5_REPAIR
2. **Original D5 SHA256:** a17ef3770a411ac4916f87171d2da0441b7f878ba1f595911c6105eabdf7e2cb
3. **Original D5 rows:** 12260
4. **All-NaN D5 regions:** 39
5. **D5 root cause:** zero valid globally sampled pixels after frozen black filter; explicit counts=0 NaN outputs
6. **Positive-control regions recomputed:** 104
7. **Positive-control fidelity:** PASS
8. **D5 regions reconstructed:** 0
9. **D5 regions technically excluded:** 0
10. **D5 unresolved:** 39
11. **Imputation used:** NO
12. **Reconstructed D5 SHA256:** NOT GENERATED
13. **Phase 0D color shortcut regression check:** FAIL — gate not reached; no reconstructed input; no changed conclusion claimed
14. **Phase 0D regions:** 12,260
15. **Existing Phase 1D embeddings:** 12,183
16. **77-region difference fully reconciled:** YES
17. **Missing embeddings legitimately excluded:** 77
18. **Embeddings newly recovered:** 0
19. **Final Phase 1 region universe:** NOT FROZEN — unresolved D5 input
20. **Epithelial injury eligible regions/donors:** 3221/19 (pre-color requirement; final universe not frozen)
21. **Fibroblast activation eligible regions/donors:** 4879/19 (pre-color requirement; final universe not frozen)
22. **Macrophage inflammatory eligible regions/donors:** 2664/19 (pre-color requirement; final universe not frozen)
23. **Scientific parameters unchanged:** YES
24. **Phase 1E–1G authorized:** NO
25. **Recommended next action:** STOP; human review of frozen valid-pixel requirement and legitimate exclusion policy; no automatic imputation or threshold changes

# Phase 1D-R audit

原D5、矩阵、cohort、folds、缓存和提取manifest均保持原始字节；未运行Phase1模型、残差或生物学图。审计用原D5 `extract_one` 函数重算整节后选择预先确定对照，保持whole-slide stride采样原点。D5使用全球150μm网格内非黑色采样像素统计，并非直接对Phikon crop做局部采样。

`EMPTY_CROP` 在根因表仅表示过滤后统计像素集合为空；不是几何crop为空。几何有效不能使空像素集合的颜色统计变成已定义值。因此即使原函数重算一致，也不能以零、均值、去掉黑色过滤或更改采样填补39个NaN，亦未将这些区域擅自标成合法技术排除。

R8、R9、R13和R14未通过前置gate，未生成重建D5、最终universe或repaired-input manifest。R13上方目标数为冻结20-cell门槛与crop边界后的参考数，不是完整最终交集。R11若全部77个缺失属于out-of-bounds，R12不适用，无需重提取embedding。

代码来源：原D5脚本未被Git追踪；报告HEAD不是该文件生成时commit的证据。精确源文件SHA与代码定义保存在d5_feature_specification.yaml。脚本没有重写颜色定义。

所有审计表位于 `results/phase1/phase1d_r/`；日志 `logs/phase1d_r_color_repair.log`。原D5列数 50，region key唯一。26图像逐个SHA、可读、shape和dtype核对。

**STOP: PHASE1E_G_AUTHORIZED = NO**

## 定量证据与边界

- 104个deterministic对照覆盖 26 个section、19 位donor、PF及control、4 个TMA。逐特征最大绝对误差为 1.1102230246251563e-16，最大相对误差为 7.276153692151615e-13；绝对误差为机器精度量级。40项特征逐项误差见 fidelity summary。
- 39个几何crop均为 uint8 RGB，所有像素都为0；不是坐标错误、dtype范围错误、HED非有限数、解码异常或直方图bin改变。RGB和HED输入/输出均有限，所有采样像素被既定黑色过滤剔除。第一个未定义步骤是有效像素数为0后的统计量计算；原代码显式输出NaN。
- 根因已确定，但39个区域的颜色特征无法按冻结定义得到finite值，且纯黑crop不是当前已证实的冻结排除理由。`D5 unresolved = 39` 表示修复/合法排除决定尚未解决，不表示像素级根因未知。
- 77个未提取区域全部为 `FROZEN_COORDINATE_SUPPORT_FAILURE`，crop out of bounds；`should_have_embedding=YES` 的缺失区域数为 0。未补提取embedding，R12不适用。
- 当前完整颜色特征交集的目标区域数见 target_universe_audit.tsv，仅用于描述缺失影响，不用于静默缩减分析集。最终universe未冻结，未生成用户规定的final_phase1_region_universe.tsv。
- Phase0D shortcut regression check在开头记为FAIL仅表示该gate未完成；没有新shortcut结果，也没有证据称原HIGH结论已改变。没有重建文件时不运行该诊断。
- 历史input gate包含reports/CURRENT_PROJECT_STATUS.md，该状态报告已按用户明确要求更新，因此其字节与2026-09-29摘要不同。它是可变执行报告，不是科学输入；标记AUTHORIZED_REPORT_UPDATE。所有真正冻结输入按历史gate或审计起始hash核对，未发现不匹配。
- R14要求所有audit PASS才生成configs/phase1_repaired_input_manifest.yaml；本次HOLD未创建该文件，也没有覆盖任何原始冻结输入。

## 停止记录

Phase1D-R执行完成，审计脚本退出码0（表示审计执行完成，不表示修复PASS）。最终判定 `TECHNICAL_HOLD_D5_REPAIR`，`SCIENTIFIC_PARAMETERS_UNCHANGED = YES`，`PHASE1E_G_AUTHORIZED = NO`。不恢复半小时监测，不运行Phase1E–1G；等待人工决定。
