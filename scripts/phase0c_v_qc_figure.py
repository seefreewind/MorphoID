#!/usr/bin/env python3
"""Cohort and no-refit coordinate-bounds overview for Phase 0C-V."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/phase0/phase0c_v"
FIG = ROOT / "figures/phase0c_v"
FIG.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.7,
    "legend.frameon": False,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
})

cross = pd.read_csv(BASE / "vannan_sample_crosswalk.tsv", sep="\t")
cross = cross[cross.published_primary == "YES"].copy()
bounds = pd.read_csv(BASE / "spatial_bounds_diagnostic.tsv", sep="\t")
data = cross.merge(bounds, on=["sample_id", "donor_id", "GSM", "TMA", "disease_status"],
                   how="left", validate="one_to_one")
data = data.sort_values(["disease_status", "donor_id", "sample_id"], kind="stable").reset_index(drop=True)
data["row"] = np.arange(len(data))

blue, amber, fail, gray = "#4F7898", "#C08A4B", "#B84A4A", "#8C9298"
status_color = {"control": blue, "pulmonary_fibrosis": amber}

fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(8.6, 8.1), sharey=True,
    gridspec_kw={"height_ratios": [1.1, 1.0], "hspace": 0.20},
)
fig.patch.set_facecolor("white")

for r in data.itertuples(index=False):
    yy = r.row
    color = status_color.get(r.disease_status, gray)
    ax1.barh(yy, r.n_cells, height=0.64, color=color, edgecolor="none", zorder=2)
    ax1.text(r.n_cells + 1500, yy, f"{int(r.n_cells):,}", va="center", ha="left", fontsize=5.7,
             color="#444444", clip_on=False)
    ax1.text(1.055, yy, r.donor_id, transform=ax1.get_yaxis_transform(), va="center", ha="left",
             fontsize=5.6, color="#454545", clip_on=False)
    ax1.text(1.19, yy, f"TMA{r.TMA}", transform=ax1.get_yaxis_transform(), va="center", ha="left",
             fontsize=5.5, color="#686868", clip_on=False)

ax1.set_xlim(0, data.n_cells.max() * 1.13)
ax1.set_ylim(len(data) - 0.5, -0.5)
ax1.set_yticks(data.row)
ax1.set_yticklabels(data.sample_id, fontsize=5.7)
ax1.set_xlabel("Corrected author-RDS cells in section", labelpad=3)
ax1.set_title("A  Donor structure and section-level cell counts", loc="left", fontweight="bold", pad=6)
ax1.grid(axis="x", color="#E5E7E9", linewidth=0.5, zorder=0)
ax1.tick_params(axis="y", length=0, pad=3)
ax1.tick_params(axis="x", labelsize=6, length=2)
ax1.text(1.055, 1.012, "Donor", transform=ax1.transAxes, ha="left", va="bottom", fontsize=5.8,
         color="#454545", fontweight="bold")
ax1.text(1.19, 1.012, "TMA", transform=ax1.transAxes, ha="left", va="bottom", fontsize=5.8,
         color="#454545", fontweight="bold")
ax1.legend(handles=[
    Line2D([0], [0], color=blue, lw=5, label=f"Control ({cross[cross.disease_status == 'control'].donor_id.nunique()} donors)"),
    Line2D([0], [0], color=amber, lw=5, label=f"Pulmonary fibrosis ({cross[cross.disease_status == 'pulmonary_fibrosis'].donor_id.nunique()} donors)"),
], loc="lower right", fontsize=5.7, ncol=2, handlelength=1.3, columnspacing=1.2, borderaxespad=0.4)

for r in data.itertuples(index=False):
    yy = r.row
    if pd.notna(r.centroids_inside_image_fraction_R2):
        x = float(r.centroids_inside_image_fraction_R2)
        color = fail if r.R2_status == "FAIL" else blue
        ax2.scatter(x, yy, s=19, marker="o", color=color, edgecolor="white", linewidth=0.35, zorder=3)
    else:
        # An open triangle is a proxy-scale diagnostic for the two images whose
        # own TIFF metadata does not provide a physical unit. It is not R2.
        x = float(r.centroids_inside_image_fraction_consensus_diagnostic)
        ax2.scatter(x, yy, s=26, marker="^", facecolor="white", edgecolor=gray, linewidth=0.9, zorder=3)

ax2.axvline(0.95, color="#30363B", linewidth=0.9, linestyle=(0, (3, 2)), zorder=1)
ax2.axvspan(0, .95, color="#B84A4A", alpha=.045, zorder=0)
ax2.set_xlim(-0.025, 1.025)
ax2.set_xlabel("Fraction of cell centroids inside registered-H&E bounds", labelpad=3)
ax2.set_title("B  Direct identity-map bounds diagnostic (R2 threshold = 0.95)",
              loc="left", fontweight="bold", pad=6)
ax2.set_yticks(data.row)
ax2.set_yticklabels(data.sample_id, fontsize=5.7)
ax2.tick_params(axis="y", length=0, pad=3)
ax2.tick_params(axis="x", labelsize=6, length=2)
ax2.grid(axis="x", color="#E5E7E9", linewidth=0.5, zorder=0)
ax2.text(.95, -1.2, "0.95 gate", ha="right", va="bottom", fontsize=5.5, color="#30363B")

fig.suptitle("GSE250346 TMA1–4: donor roster is broad; image registration remains unresolved",
             x=0.31, y=0.994, ha="left", fontsize=9.5, fontweight="bold")
fig.text(0.31, 0.006,
         "28 sections / 19 donors; repeated cores remain donor-grouped. R2 uses the published common frame with no fitted offset, rotation or warp. "
         "Blue dots = R2 pass, red = fail, open triangles = missing-scale 0.2125 proxy only. R3–R6 remain incomplete; no marker denotes registration PASS.",
         ha="left", va="bottom", fontsize=5.6, color="#51565A", wrap=True)
fig.subplots_adjust(left=0.265, right=0.79, top=0.92, bottom=0.07)

out = FIG / "vannan_cohort_registration_bounds.png"
fig.savefig(out, dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(FIG / "vannan_cohort_registration_bounds.svg", bbox_inches="tight", facecolor="white")
fig.savefig(FIG / "vannan_cohort_registration_bounds.tiff", dpi=600, bbox_inches="tight",
            facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print(f"Wrote {out}")
