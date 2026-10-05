#!/usr/bin/env python3
"""Five claim-first, donor-aware figures for the Phase 1 minimum pilot.

Figure contract
---------------
Core conclusion: at frozen 150-um scale, only an out-of-donor M4-vs-M2 gain
that persists in residual, PF-only, stain-normalized/grayscale and donor views
supports morphology-specific molecular-state information.
Archetype: quantitative grid.
Panel map: 1 model ladder; 2 morphology increment vs M2 baseline; 3 held-out
cross-fitted residual prediction; 4 shortcut stress tests; 5 donor heterogeneity.
Evidence hierarchy: donor-held-out primary estimates first, shortcut and
donor-resolved diagnostics second. Intervals are 95% donor-cluster bootstrap CIs.
Export: 183-mm-wide SVG (editable text) plus 300-dpi PNG preview; no PDF.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})  # No PDF export requested.
plt.rcParams["font.size"] = 7
plt.rcParams["axes.labelsize"] = 7
plt.rcParams["axes.titlesize"] = 8
plt.rcParams["xtick.labelsize"] = 6
plt.rcParams["ytick.labelsize"] = 6
plt.rcParams["legend.fontsize"] = 6
plt.rcParams["axes.linewidth"] = 0.65
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "results/phase1"
FIGDIR = P1 / "figures"
TARGET_LABELS = {
    "epithelial_injury": "Epithelial injury",
    "fibroblast_activation": "Fibroblast activation",
    "macrophage_inflammatory": "Macrophage inflammatory state",
}
STRATA = ["POOLED", "PF_ONLY", "DISEASE_ADJUSTED"]
STRATUM_LABELS = {"POOLED": "Pooled", "PF_ONLY": "PF-only", "DISEASE_ADJUSTED": "Disease-adjusted"}
STRATUM_COLORS = {"POOLED": "#767676", "PF_ONLY": "#C77728", "DISEASE_ADJUSTED": "#0F4D92"}
MODEL_ORDER = ["M0", "M1", "M2", "M3", "M4"]
MODEL_LABELS = {"M0": "Z", "M1": "C + Z", "M2": "C + N + Z", "M3": "X + Z", "M4": "C + N + Z + X"}
SHORTCUTS = ["raw", "macenko", "grayscale", "color_only", "coordinates_only"]
SHORTCUT_LABELS = {"raw": "Raw H&E + C/N/Z", "macenko": "Macenko + C/N/Z",
                   "grayscale": "Grayscale + C/N/Z", "color_only": "Color only",
                   "coordinates_only": "Coordinates only"}
SHORTCUT_COLORS = {"raw": "#0F4D92", "macenko": "#3775BA", "grayscale": "#42949E",
                   "color_only": "#C77728", "coordinates_only": "#9A4D8E"}


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(P1 / name, sep="\t")


def save(fig: plt.Figure, stem: str) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def ci_lookup(boot: pd.DataFrame, model: str, target: str, stratum: str,
              metric: str = "R2") -> tuple[float, float]:
    row = boot[(boot.metric == metric) & (boot.model == model)
               & (boot.target == target) & (boot.analysis_stratum == stratum)]
    if row.empty:
        return np.nan, np.nan
    return float(row.iloc[0].CI_low), float(row.iloc[0].CI_high)


def asymmetric_errors(estimate: float, low: float, high: float) -> np.ndarray:
    return np.asarray([[max(0.0, estimate-low)], [max(0.0, high-estimate)]])


def figure1_ladder() -> None:
    metrics, boot = read("tables/table2_primary_model_performance.tsv"), read("donor_bootstrap_metrics.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.35), sharey=True)
    x = np.arange(len(MODEL_ORDER))
    finite = []
    for target in TARGET_LABELS:
        sub = metrics[(metrics.target == target) & metrics.analysis_stratum.isin(STRATA)]
        finite.extend(sub[[f"{m}_R2" for m in MODEL_ORDER]].to_numpy(float).ravel())
        bsub = boot[(boot.target == target) & boot.analysis_stratum.isin(STRATA)
                    & boot.model.isin(MODEL_ORDER)]
        finite.extend(bsub[["CI_low", "CI_high"]].to_numpy(float).ravel())
    finite = np.asarray(finite)
    finite = finite[np.isfinite(finite)]
    ymin = min(-0.25, float(finite.min()) - 0.08) if len(finite) else -0.25
    ymax = max(0.8, float(finite.max()) + 0.10) if len(finite) else 0.8
    for ax, target in zip(axes, TARGET_LABELS):
        for stratum in STRATA:
            row = metrics[(metrics.target == target) & (metrics.analysis_stratum == stratum)]
            if row.empty:
                continue
            row = row.iloc[0]
            ys, lows, highs = [], [], []
            for model in MODEL_ORDER:
                value = float(row[f"{model}_R2"])
                low, high = ci_lookup(boot, model, target, stratum)
                ys.append(value); lows.append(low); highs.append(high)
            yerr = np.asarray([[max(0, v-lo) for v, lo in zip(ys, lows)],
                               [max(0, hi-v) for v, hi in zip(ys, highs)]])
            ax.errorbar(x, ys, yerr=yerr, marker="o", ms=3.0, lw=1.1, capsize=1.8,
                        color=STRATUM_COLORS[stratum], label=STRATUM_LABELS[stratum],
                        alpha=0.95, zorder=3)
        ax.axhline(0, color="#333333", lw=0.75, ls="--", zorder=1)
        ax.set_title(f"{chr(97 + list(TARGET_LABELS).index(target))}  {TARGET_LABELS[target]}\n19 donors; PF-only 13 donors", pad=7)
        ax.set_xticks(x, MODEL_ORDER)
        ax.set_ylim(ymin, ymax)
        ax.grid(axis="y", color="#E5E5E5", lw=0.45)
        ax.set_xlabel("Model")
    axes[0].set_ylabel("Donor-held-out out-of-fold R²")
    axes[-1].legend(frameon=False, loc="best", handlelength=1.5)
    fig.suptitle("Model ladder under donor-held-out validation", y=1.02, fontsize=9)
    fig.tight_layout(w_pad=0.8)
    save(fig, "figure1_model_ladder")


def figure2_increment() -> None:
    metrics, boot = read("tables/table2_primary_model_performance.tsv"), read("donor_bootstrap_metrics.tsv")
    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    marker = {"epithelial_injury": "o", "fibroblast_activation": "s", "macrophage_inflammatory": "^"}
    for _, row in metrics[metrics.analysis_stratum.isin(STRATA)].iterrows():
        target, stratum = row.target, row.analysis_stratum
        x, y = float(row.M2_R2), float(row.deltaR2_morphology)
        xlo, xhi = ci_lookup(boot, "M2", target, stratum)
        ylo, yhi = ci_lookup(boot, "paired_increment", target, stratum, "deltaR2_morphology")
        ax.errorbar(x, y, xerr=asymmetric_errors(x, xlo, xhi),
                    yerr=asymmetric_errors(y, ylo, yhi), marker=marker[target], ms=5,
                    color=STRATUM_COLORS[stratum], mec="white", mew=0.45,
                    ecolor=STRATUM_COLORS[stratum], capsize=2, lw=0.9, linestyle="none", zorder=3)
    ax.axhline(0, color="#303030", lw=0.8, ls="--")
    ax.set_xlabel("M2 out-of-fold R² (C + N + Z)")
    ax.set_ylabel("ΔR² morphology (M4 − M2)")
    ax.set_title("Morphology increment against the composition/neighborhood oracle")
    ax.grid(color="#E5E5E5", lw=0.45)
    handles = [plt.Line2D([0], [0], marker=marker[t], color="white", markerfacecolor="#555555",
                          markersize=5, label=TARGET_LABELS[t]) for t in TARGET_LABELS]
    handles += [plt.Line2D([0], [0], marker="o", color=STRATUM_COLORS[s], lw=0,
                           markersize=5, label=STRATUM_LABELS[s]) for s in STRATA]
    ax.legend(handles=handles, frameon=False, ncol=2, loc="best")
    fig.tight_layout()
    save(fig, "figure2_morphology_increment")


def figure3_residual() -> None:
    oof = read("oof_predictions.tsv")
    table = read("tables/table3_crossfit_residual_performance.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8))
    for ax, target in zip(axes, TARGET_LABELS):
        g = oof[(oof.target == target) & (oof.analysis_stratum == "DISEASE_ADJUSTED")]
        y = g.heldout_base_residual.to_numpy(float)
        p = g.morphology_predicted_residual.to_numpy(float)
        finite = np.isfinite(y) & np.isfinite(p)
        y, p = y[finite], p[finite]
        lim = float(np.max(np.abs(np.r_[y, p]))) * 1.03 if len(y) else 1.0
        lim = max(lim, 1e-3)
        ax.scatter(y, p, s=2, color="#0F4D92", alpha=0.22, rasterized=True,
                   linewidths=0)
        ax.set_xscale("symlog", linthresh=0.25)
        ax.set_yscale("symlog", linthresh=0.25)
        ax.axhline(0, color="#777777", lw=0.5)
        ax.axvline(0, color="#777777", lw=0.5)
        ax.plot([-lim, lim], [-lim, lim], color="#303030", lw=0.8, ls="--")
        row = table[(table.target == target) & table.analysis_stratum.eq("DISEASE_ADJUSTED")].iloc[0]
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        from matplotlib.ticker import FixedLocator, FuncFormatter
        ticks = [v for v in [-10, -1, 0, 1, 10] if abs(v) <= lim]
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_locator(FixedLocator(ticks))
            axis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
            axis.set_minor_locator(FixedLocator([]))
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{chr(97 + list(TARGET_LABELS).index(target))}  {TARGET_LABELS[target]}", pad=6)
        ax.text(0.04, 0.96, f"R² = {row.residual_R2:.2f}\n95% CI [{row.CI_low:.3f}, {row.CI_high:.3f}]\nr = {row.Pearson:.2f}\ndonors = {g.donor_id.nunique()}",
                transform=ax.transAxes, ha="left", va="top", fontsize=6,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5})
        ax.grid(color="#E5E5E5", lw=0.4)
    axes[0].set_ylabel("Predicted held-out residual")
    for ax in axes:
        ax.set_xlabel("Observed M2 residual (symlog; linear ±0.25)")
    fig.suptitle("Morphology prediction of expression variation left by cross-fitted M2", y=1.02, fontsize=9)
    fig.tight_layout(w_pad=0.8)
    save(fig, "figure3_crossfit_residual")


def figure4_shortcuts() -> None:
    metrics, boot = read("tables/table4_shortcut_stress_test.tsv"), read("donor_bootstrap_metrics.tsv")
    fig, axes = plt.subplots(3, 3, figsize=(7.2, 5.0), sharex=True)
    y = np.arange(len(SHORTCUTS))
    for i, target in enumerate(TARGET_LABELS):
        for j, stratum in enumerate(STRATA):
            ax = axes[i, j]
            row = metrics[(metrics.target == target) & (metrics.analysis_stratum == stratum)]
            if row.empty:
                ax.set_axis_off(); continue
            row = row.iloc[0]
            for k, model in enumerate(SHORTCUTS):
                estimate = float(row[model])
                # The Macenko and grayscale arms are M2 + the respective morphology embedding.
                model_key = {"raw": "M4", "macenko": "macenko", "grayscale": "grayscale",
                             "color_only": "color", "coordinates_only": "coordinate"}[model]
                low, high = ci_lookup(boot, model_key, target, stratum)
                err = asymmetric_errors(estimate, low, high)
                ax.errorbar(estimate, y[k], xerr=err, marker="o", ms=3.2, capsize=1.8,
                            color=SHORTCUT_COLORS[model], ecolor=SHORTCUT_COLORS[model], lw=0.9)
            ax.axvline(0, color="#303030", lw=0.65, ls="--")
            ax.grid(axis="x", color="#E5E5E5", lw=0.4)
            ax.set_title(f"{chr(97 + i * 3 + j)}  {TARGET_LABELS[target]}\n{STRATUM_LABELS[stratum]} ({13 if stratum == 'PF_ONLY' else 19} donors)", fontsize=6.5, pad=4)
            ax.set_yticks(y)
            ax.invert_yaxis()
            if j == 0:
                ax.set_yticklabels([SHORTCUT_LABELS[m] for m in SHORTCUTS], fontsize=5.5)
            else:
                ax.set_yticklabels([])
            if i == 2:
                ax.set_xlabel("Out-of-fold R²", fontsize=6)
    fig.suptitle("Shortcut stress test and morphology preprocessing sensitivity", y=1.01, fontsize=9)
    fig.tight_layout(h_pad=0.75, w_pad=0.45)
    save(fig, "figure4_shortcut_stress_test")


def figure5_heterogeneity() -> None:
    donors = read("donor_heterogeneity.tsv")
    donors = donors[donors.analysis_stratum.eq("DISEASE_ADJUSTED")].copy()
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 5.6))
    for ax, target in zip(axes, TARGET_LABELS):
        g = donors[donors.target == target].sort_values("delta_MAE_M2_minus_M4", ascending=True).reset_index(drop=True)
        assert g.n_regions.gt(0).all()  # log1p marker-size domain
        yy = np.arange(len(g))
        delta = g.delta_MAE_M2_minus_M4.to_numpy(float)
        ax.scatter(delta, yy, s=12 + 14*np.log1p(g.n_regions.to_numpy(float)),
                   color="#0F4D92", alpha=0.8, edgecolor="white", linewidth=0.35, zorder=3)
        ax.axvline(0, color="#303030", lw=0.8, ls="--")
        ax.set_yticks(yy)
        ax.set_yticklabels(g.donor_id.tolist(), fontsize=5.5)
        ax.set_title(f"{chr(97 + list(TARGET_LABELS).index(target))}  {TARGET_LABELS[target]}", pad=6)
        ax.set_xlabel("ΔMAE (M2 − M4)")
        ax.grid(axis="x", color="#E5E5E5", lw=0.45)
        positive = int(np.count_nonzero(delta > 0))
        ax.set_title(f"{chr(97 + list(TARGET_LABELS).index(target))}  {TARGET_LABELS[target]}\n{positive}/{len(g)} donors improve", pad=6)
    axes[0].set_ylabel("Held-out donor (size increases with log region count)")
    fig.suptitle("Donor-level heterogeneity in morphology error reduction", y=1.01, fontsize=9)
    fig.tight_layout(w_pad=0.8)
    save(fig, "figure5_donor_heterogeneity")


def main() -> None:
    required = [P1 / "oof_predictions.tsv", P1 / "tables/table2_primary_model_performance.tsv",
                P1 / "tables/table3_crossfit_residual_performance.tsv", P1 / "tables/table4_shortcut_stress_test.tsv",
                P1 / "donor_bootstrap_metrics.tsv", P1 / "donor_heterogeneity.tsv"]
    missing = [str(p.name) for p in required if not p.is_file()]
    if missing:
        raise FileNotFoundError("model outputs are not complete: " + ", ".join(missing))
    figure1_ladder(); figure2_increment(); figure3_residual(); figure4_shortcuts(); figure5_heterogeneity()
    print("PHASE1_FIGURES_PASS: 5 PNG + 5 SVG exports", flush=True)


if __name__ == "__main__":
    main()
