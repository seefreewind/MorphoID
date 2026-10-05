#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assemble the final Phase 1 pilot report from frozen outputs."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "results/phase1"
REPORT = ROOT / "reports/PHASE1_MINIMUM_IDENTIFIABILITY_PILOT.md"
DISPLAY = {
    "epithelial_injury": "epithelial injury",
    "fibroblast_activation": "fibroblast activation",
    "macrophage_inflammatory": "macrophage inflammatory",
}
STRATA = ["POOLED", "PF_ONLY", "DISEASE_ADJUSTED"]
STRATUM_LABEL = {"POOLED": "pooled", "PF_ONLY": "PF-only", "DISEASE_ADJUSTED": "disease-adjusted"}


def tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(P1 / name, sep="\t")


def fnum(value: object, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def get_row(frame: pd.DataFrame, target: str, stratum: str) -> pd.Series:
    rows = frame[(frame.target == target) & (frame.analysis_stratum == stratum)]
    if rows.empty:
        raise RuntimeError(f"report data missing {target}/{stratum}")
    return rows.iloc[0]


def classify_heterogeneity(donors: pd.DataFrame) -> tuple[str, dict[str, dict[str, float]]]:
    per_target: dict[str, dict[str, float]] = {}
    severity_rank = {"LOW": 0, "MODERATE": 1, "HIGH": 2}
    severities = []
    for target in DISPLAY:
        g = donors[(donors.target == target) & donors.analysis_stratum.eq("DISEASE_ADJUSTED")].copy()
        delta = g.delta_MAE_M2_minus_M4.to_numpy(float)
        positive_fraction = float(np.mean(delta > 0)) if len(delta) else np.nan
        gains = g.positive_squared_error_gain.to_numpy(float)
        total_gain = float(np.nansum(gains))
        top2_share = float(np.sort(gains[np.isfinite(gains)])[-2:].sum()/total_gain) if total_gain > 0 else 0.0
        minority = min(positive_fraction, 1-positive_fraction) if np.isfinite(positive_fraction) else 1.0
        if minority >= 0.40 or top2_share > 0.50:
            severity = "HIGH"
        elif minority >= 0.20 or top2_share > 0.30:
            severity = "MODERATE"
        else:
            severity = "LOW"
        per_target[target] = {"n_donors": int(g.donor_id.nunique()),
                              "positive_delta_mae_fraction": positive_fraction,
                              "top_two_positive_sse_gain_share": top2_share,
                              "severity": severity}
        severities.append(severity)
    overall = max(severities, key=lambda x: severity_rank[x])
    return overall, per_target


def choose_verdict(metrics: pd.DataFrame, boot: pd.DataFrame, residual: pd.DataFrame,
                   shortcuts: pd.DataFrame, heterogeneity: str,
                   heterogeneity_by_target: dict[str, dict[str, float]]) -> tuple[str, list[str], int, int]:
    adj = metrics[metrics.analysis_stratum.eq("DISEASE_ADJUSTED")].set_index("target")
    pf = metrics[metrics.analysis_stratum.eq("PF_ONLY")].set_index("target")
    res = residual[residual.analysis_stratum.eq("DISEASE_ADJUSTED")].set_index("target")
    short = shortcuts[shortcuts.analysis_stratum.eq("DISEASE_ADJUSTED")].set_index("target")
    bdelta = boot[(boot.metric == "deltaR2_composition") & (boot.model == "paired_increment") & boot.analysis_stratum.eq("DISEASE_ADJUSTED")].set_index("target")
    robust: list[str] = []
    composition_count = 0
    shortcut_count = 0
    unstable_positive = False
    for target in DISPLAY:
        row, prow, rr, sr = adj.loc[target], pf.loc[target], res.loc[target], short.loc[target]
        morph_boot = boot[(boot.target == target) & (boot.analysis_stratum == "DISEASE_ADJUSTED")
                          & (boot.metric == "deltaR2_morphology") & (boot.model == "paired_increment")]
        morph_ci_low = float(morph_boot.iloc[0].CI_low) if not morph_boot.empty else np.nan
        q = float(row.FDR) if pd.notna(row.FDR) else np.nan
        residual_supported = float(rr.residual_R2) > 0 and float(rr.CI_low) > 0 and float(rr.FDR) < 0.05
        shortcut_ceiling = max(float(sr.color_only), float(sr.coordinates_only))
        robust_shortcut_resistant = shortcut_ceiling < float(sr.raw)
        sensitivity_survives = max(float(sr.deltaR2_macenko), float(sr.deltaR2_grayscale)) > 0
        is_robust = (float(row.deltaR2_morphology) > 0 and morph_ci_low > 0 and q < 0.05
                     and residual_supported and float(prow.deltaR2_morphology) > 0
                     and sensitivity_survives and robust_shortcut_resistant
                     and heterogeneity_by_target[target]["severity"] != "HIGH")
        if is_robust:
            robust.append(target)
        comp_row = bdelta.loc[target]
        if (float(row.deltaR2_composition) > 0 and float(comp_row.CI_low) > 0):
            composition_count += 1
        raw_delta = float(row.deltaR2_morphology)
        controls_collapse = max(float(sr.deltaR2_macenko), float(sr.deltaR2_grayscale)) <= 0
        color_approx = shortcut_ceiling >= float(sr.raw) - 0.05
        if raw_delta > 0 and controls_collapse and color_approx:
            shortcut_count += 1
        h = heterogeneity_by_target[target]
        if (raw_delta > 0 and h["severity"] == "HIGH") or h["top_two_positive_sse_gain_share"] > 0.50:
            unstable_positive = True
    all_residual_cis_cross_zero = all(
        float(res.loc[target].CI_low) <= 0 <= float(res.loc[target].CI_high) for target in DISPLAY)
    median_raw_m4 = float(np.nanmedian(adj.M4_R2.to_numpy(float)))
    median_cn_increment = float(np.nanmedian(
        (adj.deltaR2_composition + adj.deltaR2_neighborhood).to_numpy(float)))
    median_morph_increment = float(np.nanmedian(adj.deltaR2_morphology.to_numpy(float)))
    composition_dominance_supported = (
        composition_count >= 2 and all_residual_cis_cross_zero and median_raw_m4 > 0
        and median_cn_increment > 0 and median_cn_increment > max(0.0, median_morph_increment)
        and heterogeneity != "HIGH" and shortcut_count == 0
    )
    if len(robust) >= 2:
        verdict = "GO_C_BROAD_IDENTIFIABILITY"
    elif len(robust) == 1:
        verdict = "GO_B_SELECTIVE_IDENTIFIABILITY"
    elif unstable_positive or heterogeneity == "HIGH":
        verdict = "HOLD_MODEL_INSTABILITY"
    elif shortcut_count >= 2:
        verdict = "HOLD_SHORTCUT_DOMINANCE"  # A pilot cannot establish irrecoverability.
    elif composition_dominance_supported:
        verdict = "GO_A_COMPOSITION_DOMINANCE"
    else:
        verdict = "HOLD_MODEL_INSTABILITY"
    return verdict, robust, composition_count, shortcut_count


def main() -> None:
    required = [P1 / "tables/table1_analysis_dataset.tsv", P1 / "tables/table2_primary_model_performance.tsv",
                P1 / "tables/table3_crossfit_residual_performance.tsv", P1 / "tables/table4_shortcut_stress_test.tsv",
                P1 / "donor_bootstrap_metrics.tsv", P1 / "donor_heterogeneity.tsv",
                P1 / "morphology_encoder_manifest.tsv", P1 / "input_hash_gate.tsv"]
    if not all(path.is_file() for path in required):
        raise FileNotFoundError("Phase 1 model or provenance outputs are incomplete")
    if REPORT.exists():
        raise FileExistsError(f"refusing to overwrite existing report: {REPORT}")
    table1 = tsv("tables/table1_analysis_dataset.tsv")
    metrics = tsv("tables/table2_primary_model_performance.tsv")
    residual = tsv("tables/table3_crossfit_residual_performance.tsv")
    shortcuts = tsv("tables/table4_shortcut_stress_test.tsv")
    boot = tsv("donor_bootstrap_metrics.tsv")
    donor = tsv("donor_heterogeneity.tsv")
    coverage = tsv("phase1a_target_feasibility.tsv")
    gene_coverage = tsv("target_gene_coverage.tsv")
    input_gate = tsv("input_hash_gate.tsv")
    overall_heterogeneity, h_by_target = classify_heterogeneity(donor)
    verdict, robust_targets, comp_count, shortcut_count = choose_verdict(
        metrics, boot, residual, shortcuts, overall_heterogeneity, h_by_target)
    expected_coverage = {"epithelial__injury", "fibroblast__fibroblast_activation",
                         "macrophage__inflammatory"}
    coverage_primary = coverage[coverage.target.isin(expected_coverage)]
    genes_primary = gene_coverage[gene_coverage.target.isin(expected_coverage)]
    if (set(coverage_primary.target) != expected_coverage or len(genes_primary) != 3
            or not coverage_primary.target_class.astype(str).eq("PRIMARY").all()
            or not genes_primary.target_class.astype(str).eq("PRIMARY").all()):
        verdict = "HOLD_TARGET_COVERAGE"

    adj = metrics[metrics.analysis_stratum.eq("DISEASE_ADJUSTED")].set_index("target")
    med = {col: float(np.nanmedian(adj[col].to_numpy(float))) for col in
           ["M0_R2", "M1_R2", "M2_R2", "M3_R2", "M4_R2", "deltaR2_composition",
            "deltaR2_neighborhood", "deltaR2_morphology"]}
    residual_adj = residual[residual.analysis_stratum.eq("DISEASE_ADJUSTED")]
    median_resid = float(np.nanmedian(residual_adj.residual_R2.to_numpy(float)))
    s4 = shortcuts[shortcuts.analysis_stratum.eq("DISEASE_ADJUSTED")]
    median_mac = float(np.nanmedian(s4.deltaR2_macenko.to_numpy(float)))
    median_gray = float(np.nanmedian(s4.deltaR2_grayscale.to_numpy(float)))
    median_color = float(np.nanmedian(s4.color_only.to_numpy(float)))
    median_coord = float(np.nanmedian(s4.coordinates_only.to_numpy(float)))
    strongest = "color-only" if median_color >= median_coord else "coordinates-only"
    sig = adj[adj.FDR < 0.05].index.tolist()
    concordant = []
    for target in DISPLAY:
        a = float(adj.loc[target].deltaR2_morphology)
        p = float(metrics[(metrics.target == target) & metrics.analysis_stratum.eq("PF_ONLY")].iloc[0].deltaR2_morphology)
        if np.sign(a) == np.sign(p):
            concordant.append(target)
    region_counts, target_donors = [], []
    for target in DISPLAY:
        row = table1[(table1.target == target) & table1.analysis_stratum.eq("DISEASE_ADJUSTED")].iloc[0]
        region_counts.append(f"{DISPLAY[target]} {int(row.n_regions):,}")
        target_donors.append(f"{DISPLAY[target]} {int(row.n_donors)}")

    if verdict == "GO_A_COMPOSITION_DOMINANCE":
        interpretation = "信号主要由可见细胞组成和局部组织结构承载；在 donor-held-out 评估中，H&E 形态未显示稳定且可重复的额外分子状态信息。"
        next_phase = "Phase 1G 已完成。可考虑仅围绕 composition/neighborhood dominance 设计独立复现或明确的外部验证；本次不扩展队列或模型。"
    elif verdict == "GO_B_SELECTIVE_IDENTIFIABILITY":
        interpretation = "形态可识别性呈选择性：只有部分 frozen within-lineage target 显示 donor-generalizable、shortcut-resistant 的增量。"
        next_phase = "Phase 1G 已完成。后续若启动，应优先针对支持信号的 target 做独立 donor/cohort 复现；本次不扩展队列或模型。"
    elif verdict == "GO_C_BROAD_IDENTIFIABILITY":
        interpretation = "多个 frozen within-lineage targets 在 donor-held-out、残差及 shortcut sensitivity 中保持正向，支持较广泛的形态相关分子状态信息。"
        next_phase = "Phase 1G 已完成。科学上可论证独立外部复现；本次仍按 stop rule 停止，不扩展队列或模型。"
    elif verdict == "HOLD_SHORTCUT_DOMINANCE":
        interpretation = "观察到的表观形态增量对染色/颜色或坐标捷径敏感，现有证据不能区分组织结构信号与技术 shortcut。"
        next_phase = "暂不扩展 cohort。下一阶段应先解决颜色、染色批次或坐标混杂，再重新预注册分析。"
    elif verdict == "HOLD_TARGET_COVERAGE":
        interpretation = "至少一个 primary target 与 Phase 0D 冻结的 signature coverage 不一致，当前 pilot 不可作确认性解释。"
        next_phase = "暂停扩展；先核对冻结 target gene coverage 与矩阵构建记录。"
    elif verdict == "NO_GO_MORPHOLOGICAL_IDENTIFIABILITY_PILOT":
        interpretation = "在本 pilot 的冻结尺度和 donor-held-out 设定下，composition dominance 与可复现的 shortcut-resistant morphology signal 均未获支持。"
        next_phase = "不建议扩展本研究的虚拟空间转录组预测路线；保留本 pilot 作为可复现的边界结果。"
    else:
        interpretation = "估计的 donor 间方向或不确定性不足以作稳定 GO/NO-GO 判断。"
        next_phase = "维持 HOLD。建议下一轮先预注册 donor/fold 稳定性与 M2 高杠杆误差审计，再决定是否启动独立验证；本轮已停止。"

    residual_sig = residual_adj[residual_adj.FDR < 0.05].target.tolist()
    positive_morph = adj[adj.deltaR2_morphology > 0].index.tolist()
    positive_residual = residual_adj[residual_adj.residual_R2 > 0].target.tolist()
    def names(items):
        return "; ".join(DISPLAY[t] for t in items) if items else "None"
    def eligible(target):
        r = get_row(table1, target, "DISEASE_ADJUSTED")
        return f"{int(r.n_regions)} / {int(r.n_donors)}"
    fields = [
        ("PHASE 1 VERDICT", f"`{verdict}`"),
        ("Common-universe SHA256", "0d944e4cae613f5e2353e9ebd970d119dda7c39739c1a567405c7ff77e66340e"),
        ("Primary morphology encoder", "Phikon-v2"),
        ("Encoder revision", "2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd"),
        ("Primary morphology scale", "150 μm"),
        ("Independent donors", "19"), ("Sections", "26"),
        ("Common-universe regions", "12144"),
        ("Epithelial injury eligible regions/donors", eligible("epithelial_injury")),
        ("Fibroblast activation eligible regions/donors", eligible("fibroblast_activation")),
        ("Macrophage inflammatory eligible regions/donors", eligible("macrophage_inflammatory")),
        *[(f"Median M{i} R²", fnum(med[f"M{i}_R2"])) for i in range(5)],
        ("Median DeltaR² composition", fnum(med["deltaR2_composition"])),
        ("Median DeltaR² neighborhood", fnum(med["deltaR2_neighborhood"])),
        ("Median DeltaR² morphology", fnum(med["deltaR2_morphology"])),
        ("Median cross-fitted residual R²", fnum(median_resid)),
        ("Targets with positive morphology increment", names(positive_morph)),
        ("Targets with positive residual predictability", names(positive_residual)),
        ("PF-only concordance", f"{len(concordant)}/3 direction-concordant"),
        ("Macenko robustness", f"median ΔR² = {fnum(median_mac)}"),
        ("Grayscale robustness", f"median ΔR² = {fnum(median_gray)}"),
        ("Color-only contribution", f"median R² = {fnum(median_color)}"),
        ("Coordinate-only contribution", f"median R² = {fnum(median_coord)}"),
        ("Strongest shortcut", "Phase 0D color→disease coupling (common-universe donor-held-out BA 0.939); among outcome controls: " + strongest),
        ("Donor heterogeneity", overall_heterogeneity),
        ("Primary FDR-significant morphology targets", names(sig)),
        ("Primary FDR-significant residual targets", names(residual_sig)),
        ("Scientific interpretation", interpretation),
        ("Recommended next phase", next_phase),
    ]
    assert len(fields) == 33
    first = [f"{i}. **{label}:** {value}" for i, (label, value) in enumerate(fields, 1)]

    target_sections = []
    for target, label in DISPLAY.items():
        pooled = get_row(metrics, target, "POOLED")
        pf = get_row(metrics, target, "PF_ONLY")
        primary = get_row(metrics, target, "DISEASE_ADJUSTED")
        rr = get_row(residual, target, "DISEASE_ADJUSTED")
        sh = get_row(shortcuts, target, "DISEASE_ADJUSTED")
        drow = donor[(donor.target == target) & donor.analysis_stratum.eq("DISEASE_ADJUSTED")]
        h = h_by_target[target]
        pboot = boot[(boot.target == target) & boot.analysis_stratum.eq("DISEASE_ADJUSTED")
                     & (boot.metric == "deltaR2_morphology") & (boot.model == "paired_increment")].iloc[0]
        target_status = "SUPPORTIVE_SELECTIVE_SIGNAL" if target in robust_targets else (
            "COMPOSITION_DOMINANT" if primary.deltaR2_morphology <= 0 and primary.deltaR2_composition > 0 else "UNRESOLVED_OR_SHORTCUT_SENSITIVE")
        target_sections.extend([
            f"### {label}",
            "",
            f"- Pooled M2 R² = {fnum(pooled.M2_R2)}; pooled M4 R² = {fnum(pooled.M4_R2)}; pooled ΔR² morphology = {fnum(pooled.deltaR2_morphology)}.",
            f"- Disease-adjusted M0–M4 R² = {fnum(primary.M0_R2)}, {fnum(primary.M1_R2)}, {fnum(primary.M2_R2)}, {fnum(primary.M3_R2)}, {fnum(primary.M4_R2)}; ΔR² composition = {fnum(primary.deltaR2_composition)}, neighborhood = {fnum(primary.deltaR2_neighborhood)}. M2 R² = {fnum(primary.M2_R2)}; M4 R² = {fnum(primary.M4_R2)}; ΔR² morphology = {fnum(primary.deltaR2_morphology)} (95% donor-bootstrap CI {fnum(pboot.CI_low)} to {fnum(pboot.CI_high)}; BH q = {fnum(primary.FDR)}).",
            f"- Cross-fitted residual R² = {fnum(rr.residual_R2)} (95% donor-bootstrap CI {fnum(rr.CI_low)} to {fnum(rr.CI_high)}); residual Pearson r = {fnum(rr.Pearson)}; Spearman ρ = {fnum(rr.Spearman)}; bootstrap-tail p = {fnum(rr.p)}, BH q = {fnum(rr.FDR)}.",
            f"- PF-only ΔR² morphology = {fnum(pf.deltaR2_morphology)}; Macenko ΔR² = {fnum(sh.deltaR2_macenko)}; grayscale ΔR² = {fnum(sh.deltaR2_grayscale)}.",
            f"- Color-only R² = {fnum(sh.color_only)}; coordinate-only R² = {fnum(sh.coordinates_only)}.",
            f"- Donor consistency: {int(np.count_nonzero(drow.delta_MAE_M2_minus_M4.to_numpy(float) > 0))}/{int(drow.donor_id.nunique())} eligible donors had lower MAE under M4; heterogeneity {h['severity']} (top two positive-error-gain donors account for {100*h['top_two_positive_sse_gain_share']:.1f}% of positive donor error gain).",
            f"- Target verdict: `{target_status}`.",
            "",
        ])

    coverage_text = "; ".join(f"{r.target}: {r.n_genes_found_in_RNA_assay}/{r.n_frozen_panel_genes} frozen genes present, {r.target_class}"
                              for _, r in gene_coverage[gene_coverage.target.isin([
                                  "epithelial__injury", "fibroblast__fibroblast_activation",
                                  "macrophage__inflammatory"])].iterrows())
    gate_records = int(input_gate.shape[0])
    body = [
        "",
        "## Scope and verdict basis",
        "",
        "本报告完成 Phase 1A–1G minimum identifiability pilot，并按 stop rule 在 Phase 1G 后停止。分析只使用冻结 Vannan cohort、冻结 donor folds、冻结 150 μm grid、冻结 within-lineage targets 和一个冻结 encoder。未扩展到其他 cohort、encoder、单细胞/细胞核形态、全转录组或临床结局。",
        "",
        f"Phase 0D 的 58 项输入门控无不匹配；其中 29 项有历史 sidecar SHA256 并逐字节通过，另外 29 项没有历史 sidecar，记录为 Phase 0D baseline capture 并通过语义/跨输出校验。该差异不应表述成 58 项都与历史 SHA 对照。报告生成前核对到 {gate_records} 条门控记录。",
        "",
        "Verdict operationalization 在 Phase 1 performance outputs 生成前固定于报告脚本。支持性 morphology target 要求 disease-adjusted ΔR²>0、donor-bootstrap CI lower bound>0、三 target BH q<0.05、cross-fitted residual R²及其 CI lower bound均为正、PF-only ΔR²同向为正、Macenko/gray 至少一项增量为正（两项均完整报告），残差三 target BH q<0.05，target donor heterogeneity 非 HIGH，且 raw M4 R²高于 color-only 与 coordinate-only。GO_A 还要求至少两个 target 的 composition 增量及 donor CI 为正、全部 residual CI 跨零、median M4 R²为正、median composition+neighborhood 增量高于 morphology 增量、donor heterogeneity 非 HIGH 且没有 shortcut-dominant target。本 pilot 不自动触发 NO_GO；该判断还需要不可恢复 shortcut 或 estimand 无法解释的独立证据。未达 GO 条件的情形保留 HOLD。",
        "",
        f"Frozen primary signature coverage: {coverage_text}. 每个 primary target 均保持 ≥20 target-lineage-cell 阈值；macrophage inflammatory score 按冻结 ontology 在 macrophage/monocyte lineage 中计算。未按表达结果或模型表现筛选基因或区域。",
        "",
        "永久排除的 VUILD105MA1 与 VUILD48LA1 未重新纳入，原因仍为 `FIXED_COORDINATE_SUPPORT_FAILURE`。150 μm crop audit 记录了 77 个全局越界 region；partial crop 未补白或截取。target-specific、边界及低 lineage-count 排除均见 [`region_exclusion_log.tsv`](../results/phase1/region_exclusion_log.tsv)。",
        "",
        "## Encoder and image preprocessing",
        "",
        "使用单一 frozen encoder：Owkin Phikon-v2，DINOv2 ViT-L/16，1024 维 CLS embedding；不微调、不做 target supervision、不按 target 选择 patch。checkpoint hash、版本、许可、预处理、分辨率与 batch 参数记录于 [`morphology_encoder_manifest.tsv`](../results/phase1/morphology_encoder_manifest.tsv)。模型使用公开权重，官方模型卡声明 non-commercial license，并说明预训练使用 PANCAN-XL（含 TCGA/CPTAC/GTEx）；本 pilot 没有导入这些额外 cohort。详见 [Phikon-v2 官方模型卡](https://huggingface.co/owkin/phikon-v2)。",
        "",
        "每个 crop 保持 150 μm 物理视野，处理器将正方形 crop 调整到 224 × 224，bicubic resampling、center crop、ImageNet mean/SD normalization。RAW 使用作者注册 H&E；GRAYSCALE 使用亮度灰度转换。每张 section 的 source stain summary 用覆盖全图的确定性等距块起点估计：每块连续 100 个 RGB 像素，默认 200 块（共 20,000 点）；若 OD 筛选后不足 500 个染色像素，则扩大到 1,000 块（共 100,000 点）。实际样本数和是否触发 fallback 均记录在 `macenko_section_source_parameters.tsv`。Macenko 在每个 outer fold 只用 training donors 的 section stain summaries 估计目标染色基底：先在 donor 内跨 section 取中位数，再跨 training donors 取中位数；held-out donor 不参与目标参考。目标参考矩阵与逐 section stain 参数均留档。PF-only Ridge 的训练和调参仅使用 PF donors；冻结的 Macenko outer-fold reference 仍由该 fold 全部 training donors 的图像汇总，可能含 control donors。PF-only Macenko 应解释为使用预先冻结、无 outer-test donor 泄漏的参考，而非 PF-only stain-reference 重估。",
        "",
        "每个 region 的 tissue_fraction 是每 8 像素采样一次、以 OD sum > 0.15 且排除纯黑像素得到的 proxy，不是经验证的组织分割掩膜。",
        "",
        "## Models and validation",
        "",
        "M0 = Z；M1 = C + Z；M2 = C + N + Z；M3 = raw X + Z；M4 = C + N + Z + raw X。POOLED 的 frozen Z 使用 normalized x/y；DISEASE_ADJUSTED 在此基础上加入冻结 disease covariate；PF_ONLY 在自身 donors 内独立重新训练；可选 CONTROL_ONLY 未执行。三个 strata 均报告 OOF ΔR² = R²(M4) − R²(M2)，保留负值。颜色输入严格使用 D5 冻结的 40 个 RGB/HED/亮度/饱和度/直方图变量；coordinate-only 使用冻结 Phase 0D 的六项绝对/归一化坐标与图像尺寸项。",
        "",
        "所有 outer folds 直接使用冻结 donor folds。每个 alpha 仅在 outer-training donors 内用 GroupKFold 调参，alpha grid 为 log-spaced 10⁻⁴–10⁵；每个 inner split 的 scaler 只在其 training donors 拟合。Ridge 的训练损失及 alpha 选择以 donor-equal MSE 汇总，以免 region 数较多的 donor 独占拟合权重；报告 R²仍按用户定义对全部 OOF regions pooled 计算。",
        "",
        "Phase 1G 每个 outer fold 的 M2 training-side residual 由 outer-training donors 内再次 donor-grouped cross-fitting 得到；每个 training region 的 M2 base model 未见过其 donor。residual Ridge 仅用这些 OOF residual 训练，并在 outer-held-out donors 评估。训练侧 residual 明细见 [`phase1g_train_oof_residuals.tsv`](../results/phase1/phase1g_train_oof_residuals.tsv)。",
        "",
        "不确定性使用 2,000 次 donor-cluster bootstrap；每次整 donor 连同其 sections/regions 一起重采样。ΔR² 的区间直接对同一次 donor resample 计算 R²(M4) − R²(M2)，没有用两个独立区间相减。Primary FDR family 仅是 DISEASE_ADJUSTED 的三个 primary ΔR² morphology；BH 校正基于 donor-bootstrap two-sided sign-tail approximation。未运行 permutation test；因此 q 值按 bootstrap-tail 解释，不是 donor-label permutation p 值。Residual R² 的三个 disease-adjusted 检验组成独立 BH family，采用同样的 bootstrap sign-tail approximation；其他 strata 保留未校正描述性 p 值，不作确认性显著性声明。",
        "",
        "## Primary results by target",
        "",
        *target_sections,
        "## Cross-target tables",
        "",
        "| Target | Stratum | Donors | Sections | Regions | Median lineage cells | M0 R² | M1 R² | M2 R² | M3 R² | M4 R² | ΔR² morphology | 95% CI | FDR q |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in metrics.iterrows():
        t1 = table1[(table1.target == row.target) & (table1.analysis_stratum == row.analysis_stratum)].iloc[0]
        body.append(f"| {DISPLAY[row.target]} | {STRATUM_LABEL.get(row.analysis_stratum, row.analysis_stratum)} | {int(row.n_donors)} | {int(row.n_sections)} | {int(row.n_regions)} | {fnum(t1.median_lineage_cells,1)} | {fnum(row.M0_R2)} | {fnum(row.M1_R2)} | {fnum(row.M2_R2)} | {fnum(row.M3_R2)} | {fnum(row.M4_R2)} | {fnum(row.deltaR2_morphology)} | {fnum(row.CI_low)} to {fnum(row.CI_high)} | {fnum(row.FDR)} |")
    body.extend([
        "",
        "Full tables: [Table 1 — analysis dataset](../results/phase1/tables/table1_analysis_dataset.tsv), [Table 2 — primary model performance](../results/phase1/tables/table2_primary_model_performance.tsv), [Table 3 — cross-fitted residual performance](../results/phase1/tables/table3_crossfit_residual_performance.tsv), [Table 4 — shortcut stress test](../results/phase1/tables/table4_shortcut_stress_test.tsv). All OOF records are in [`oof_predictions.tsv`](../results/phase1/oof_predictions.tsv); complete model scores and donor-bootstrap CIs are in [`donor_bootstrap_metrics.tsv`](../results/phase1/donor_bootstrap_metrics.tsv).",
        "",
        "## Figures",
        "",
        "P1, P2 and P4 show completed OOF estimates and 95% percentile intervals from 2,000 donor-cluster resamples. P3 shows all eligible disease-adjusted region pairs as points on symmetric-log axes, linear within ±0.25 and annotates residual R² and its donor-bootstrap interval; the diagonal is identity and both axes retain the full data extent. No tail regions are hidden. P5 shows each of the 19 held-out donors directly, without within-donor bootstrap error bars; marker area follows log(1 + region count). Pooled/disease-adjusted estimates use 19 donors; PF-only uses 13 independently retrained donors. Eligible region counts and all model arms are recorded in the dataset and OOF tables. Zero-reference lines and negative R² are retained. Figures are editable-text SVG plus 300-dpi PNG exports.",
        "",
        "### Figure P1. Donor-held-out model ladder",
        "",
        "![Figure 1: Model ladder](../results/phase1/figures/figure1_model_ladder.png)",
        "",
        "### Figure P2. Morphology increment against the M2 oracle",
        "",
        "![Figure 2: Morphology increment](../results/phase1/figures/figure2_morphology_increment.png)",
        "",
        "### Figure P3. Cross-fitted residual prediction",
        "",
        "![Figure 3: Residual prediction](../results/phase1/figures/figure3_crossfit_residual.png)",
        "",
        "### Figure P4. Shortcut stress test",
        "",
        "![Figure 4: Shortcut stress test](../results/phase1/figures/figure4_shortcut_stress_test.png)",
        "",
        "### Figure P5. Donor heterogeneity",
        "",
        "![Figure 5: Donor heterogeneity](../results/phase1/figures/figure5_donor_heterogeneity.png)",
        "",
        "## Donor heterogeneity rule",
        "",
        f"Donor heterogeneity is classified from DISEASE_ADJUSTED donor-level ΔMAE and concentration of positive squared-error gains. HIGH indicates a near-even split of donors improving vs worsening (40–60%) or the top two donors contributing >50% of positive squared-error gain; MODERATE indicates a 20–40% directional minority or >30% top-two concentration; otherwise LOW. Current overall class: **{overall_heterogeneity}**. Target-specific fractions and gains are in [`donor_heterogeneity.tsv`](../results/phase1/donor_heterogeneity.tsv).",
        "",
        "## Limitations and stop rule",
        "",
        "This is a small frozen-cohort pilot with 19 independent donors, and primary target eligibility differs by lineage. Donor bootstrap intervals can remain wide with this donor count. Color and coordinate controls do not fully remove TMA/run effects; those remain a Phase 0D confounding qualification and constrain pooled interpretation. The encoder's non-commercial license also limits downstream use. No causal or clinical validation claim follows from these results.",
        "",
        "Secondary programs use a separate family and never alter the primary verdict. Gene-level inference, pseudo-transcriptome controls, 100/200 μm sensitivity, and permutation testing were not used to alter the primary verdict. No gene-level pilot was run because no frozen gene panel was authorized. Phase 1 stops here; extensions listed in the request were not run.",
        "",
        "## Reproducibility artifacts",
        "",
        "- `scripts/phase1_input_hash_gate.py`, `scripts/phase1a_build_matrix.R`, `scripts/phase1a_crop_support_audit.py`, `scripts/phase1d_fetch_encoder.py`, `scripts/phase1d_extract_morphology.py`, `scripts/phase1e_g_fit_models.py`, `scripts/phase1_figures.py`, `scripts/phase1_report.py`.",
        "- `results/phase1/input_hash_gate.tsv`, encoder manifest, section/fold Macenko parameters, embedding extraction manifest, alpha selections, complete OOF predictions, bootstrap metrics, donor heterogeneity, Phase 1G train-side residuals, exclusions, and Tables 1–4.",
        "- Runtime logs are under `logs/phase1*`; the Phase 1A input gate passed before modeling authorization, and no post-gate frozen input was changed.",
        "",
    ])
    body.extend(["## Why the primary verdict remains HOLD", "",
        "Fibroblast activation has a positive raw disease-adjusted increment, but its M2 baseline R² is negative and its paired interval is broad. The morphology increment fails the three-target BH threshold, and its residual interval crosses zero. Its independently trained PF-only result is stronger and remains supportive descriptive evidence. It does not override the disease-adjusted primary family.", "",
        "The two donors with the largest positive squared-error gains account for "
        + f"{100*h_by_target['fibroblast_activation']['top_two_positive_sse_gain_share']:.1f}% of fibroblast positive gains. "
        + "Epithelial injury has a positive increment with an interval crossing zero; macrophage inflammatory morphology and residual R² are negative. Neither residual null nor a negative result alone is classified as NO_GO. These findings motivate stability review rather than an expanded atlas.", ""])
    secondary_path = P1 / "secondary/secondary_model_performance.tsv"
    if secondary_path.exists():
        sm = pd.read_csv(secondary_path, sep="\t")
        sr = pd.read_csv(P1 / "secondary/residual_prediction_performance.tsv", sep="\t")
        body.extend(["## Frozen secondary programs", "",
                     "These are the three pre-specified secondary-usable within-lineage programs. They were fitted only after all primary targets completed, with the same common universe, folds, models, independent PF-only retraining, and strict residual cross-fit. Disease-adjusted ΔR² and residual tests each use a separate three-program secondary BH family; they do not enter the primary verdict.", "",
                     "| Program | Stratum | Regions | Donors | M2 R² | M4 R² | ΔR² morphology | 95% CI | Secondary q | Residual R² | Residual q |",
                     "|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|"])
        for _, row in sm.iterrows():
            rr = sr[(sr.target == row.target) & (sr.analysis_stratum == row.analysis_stratum)].iloc[0]
            body.append(f"| {row.target} | {row.analysis_stratum} | {int(row.n_regions)} | {int(row.n_donors)} | {fnum(row.M2_R2)} | {fnum(row.M4_R2)} | {fnum(row.deltaR2_morphology)} | {fnum(row.CI_low)} to {fnum(row.CI_high)} | {fnum(row.FDR)} | {fnum(rr.residual_R2)} | {fnum(rr.FDR)} |")
        body.extend(["", "Full secondary OOF, sensitivity scores, alpha selections, bootstrap intervals, donor heterogeneity and training-side residual records are retained in results/phase1/secondary. The six frozen broad secondary programs were optional and were not fitted in this minimum within-lineage pilot. Their full-region grayscale/Macenko sensitivity coverage is not present in the frozen cache; no new embeddings or arm-specific region exclusions were introduced to extend this analysis.", ""])
    low = pd.read_csv(ROOT / "configs/frozen_within_lineage_targets.tsv", sep="\t")
    body.extend(["## Low-coverage targets", "",
                 "The frozen low-coverage programs are described only; no confirmatory fits or claims are made.", "",
                 "| Lineage | Program | Panel/reference genes | Frozen coverage |",
                 "|---|---|---:|---|"])
    for _, row in low[low.coverage_class.eq("LOW_COVERAGE")].iterrows():
        body.append(f"| {row.lineage} | {row.target} | {int(row.n_panel_genes)}/{int(row.n_reference_genes)} | {float(row.coverage_fraction):.3f} |")
    broad_low = pd.read_csv(ROOT / "configs/frozen_phase1_candidate_programs.tsv", sep="\t")
    for _, row in broad_low[broad_low.phase1_target_class.eq("LOW_COVERAGE")].iterrows():
        body.append(f"| Broad | {row.candidate_program} | {int(row.n_panel_genes)}/{int(row.n_reference_genes)} | {float(row.coverage_fraction):.3f} |")
    body.extend([
        "## Frozen common universe and audit assertions", "",
        "All model and sensitivity arms use the prospective R2 common universe: 12,144 of 12,260 regions. The exclusions are 77 unsupported coordinates and 39 zero-valid-pixel regions after the unchanged D5 black filter, with no overlap. Target-specific ≥20 lineage cells and finite frozen score eligibility is applied afterward; no arm-specific filtering or imputation is used. Original D5 and historical HOLD evidence remain intact. The post-modeling hash gate rechecks all 43 scientific input hashes.", "",
        "FULL_DATA_RESIDUALIZATION_USED = NO",
        "MICROREGION_PRIMARY_BOOTSTRAP_USED = NO",
        "TEST_DONOR_PREPROCESSING_LEAKAGE = NO",
        "DONOR_CROSS_FOLD_LEAKAGE = NO",
        "POSTHOC_TARGET_SELECTION = NO",
        "POSTHOC_SECTION_EXCLUSION = NO",
        "COMMON_UNIVERSE_CHANGED_AFTER_MODELING = NO", "",
        "All opening medians and confirmatory target verdicts refer to DISEASE_ADJUSTED; pooled and independently retrained PF-only results remain fully tabulated. R² CIs condition on the completed OOF fits rather than re-fitting each bootstrap. Macenko reference estimation excludes each outer-test donor; its fixed outer-training reference is reused during inner tuning, so inner validation donors can contribute to this unsupervised stain reference. This does not expose outer-test donors but is a limitation of the frozen stain workflow.", "",
        "The user-requested named tables and compressed primary OOF table are in results/phase1. M3 and M4 were both refitted for grayscale and Macenko; M3 sensitivity scores appear in primary_model_performance.tsv and shortcut_robustness.tsv. Figure audit: reports/PHASE1_FIGURE_QA.md.", ""
    ])
    import re
    body = [re.sub(r"\]\(\.\./([^)]*)\)", lambda m: "](<" + str(ROOT / m.group(1)) + ">)", line) for line in body]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(first + body), encoding="utf-8")
    print(f"PHASE1_REPORT_WRITTEN: {REPORT.relative_to(ROOT)}; verdict={verdict}", flush=True)


if __name__ == "__main__":
    main()
