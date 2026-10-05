# Phase 1 Figure QA

**PILOT FIGURE QA: PASS_WITH_EXPORT_POLICY_EXCEPTION**

## Contract and source-data mapping

Backend: Python/matplotlib exclusively. Archetype: quantitative grid. The evidence chain is model ladder → paired increment → strict residual prediction → shortcut robustness → all-donor heterogeneity. Figures were generated only from the completed primary OOF predictions, performance tables, donor bootstrap intervals and donor diagnostic table. No simulated data, representative H&E selection or aesthetics-driven observation exclusions were used.

Pooled and disease-adjusted strata: 19 donors; independently refitted PF-only: 13 donors. Common universe: 12,144 regions, with target eligibility applied afterward (epithelial 3,219; fibroblast 4,869; macrophage 2,661). Every comparable model and control uses identical target/stratum rows. All 45 outer folds completed. Negative R² and paired ΔR² remain unmodified.

## Automated preflight and export policy

Nature-figure source validation completed. The remaining default FAIL is EXPORT-VECTOR because its generic submission bundle requires PDF. The user's output preference prohibits automatic PDF creation; this pilot therefore delivers editable SVG and 300-dpi PNG only. This export-policy exception is explicit and does not affect numerical or visual validation. No PDF was generated and no PDF glyph audit was needed. TIFF/600-dpi warnings refer to production submission formats, which were not requested. PNG files are pilot previews, not certified final journal production files.

Editable SVG text settings passed. Exported SVGs were independently parsed: all five contain selectable text, actual minimum text sizes are 5.5–6 pt, widths are approximately 183 mm, and PNG metadata reports approximately 300 dpi. Symmetric-log residual ticks were replaced with plain numeric labels to avoid undersized math superscripts. The source validator and export metadata are retained under results/phase1.

## Panel-by-panel visual inspection

Each final panel and its assembled figure was viewed after rendering. Titles, target names, donor counts, zero references, CI endpoints, negative values, label clearance, legend encoding and axes were checked against the source tables. No duplicated scientific panels were found.

| Panel | Unique claim | Center/summary | Spread/interval | Replicate unit | Labels/legend | Collision check | Verdict |
|---|---|---|---|---|---|---|---|
| P1a | M0–M4 held-out model ladder | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P1b | M0–M4 held-out model ladder | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P1c | M0–M4 held-out model ladder | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P2a | paired morphology increment vs M2 oracle | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P3a | strict cross-fitted residual prediction | completed OOF estimate | annotated donor-bootstrap residual R² CI; all region pairs displayed | donor | checked | checked | PASS |
| P3b | strict cross-fitted residual prediction | completed OOF estimate | annotated donor-bootstrap residual R² CI; all region pairs displayed | donor | checked | checked | PASS |
| P3c | strict cross-fitted residual prediction | completed OOF estimate | annotated donor-bootstrap residual R² CI; all region pairs displayed | donor | checked | checked | PASS |
| P4a | raw/stain/grayscale/color/coordinate controls | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P4b | raw/stain/grayscale/color/coordinate controls | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P4c | raw/stain/grayscale/color/coordinate controls | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P4d | raw/stain/grayscale/color/coordinate controls | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P4e | raw/stain/grayscale/color/coordinate controls | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P4f | raw/stain/grayscale/color/coordinate controls | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P4g | raw/stain/grayscale/color/coordinate controls | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P4h | raw/stain/grayscale/color/coordinate controls | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P4i | raw/stain/grayscale/color/coordinate controls | completed OOF estimate | 95% donor bootstrap CI (2,000 resamples) | donor | checked | checked | PASS |
| P5a | all held-out donors and direction heterogeneity | completed OOF estimate | none: individual donor diagnostics, no microregion uncertainty | individual held-out donor | checked | checked | PASS |
| P5b | all held-out donors and direction heterogeneity | completed OOF estimate | none: individual donor diagnostics, no microregion uncertainty | individual held-out donor | checked | checked | PASS |
| P5c | all held-out donors and direction heterogeneity | completed OOF estimate | none: individual donor diagnostics, no microregion uncertainty | individual held-out donor | checked | checked | PASS |

## Revisions made during QA

- P1 axis limits include every model CI, including the very negative fibroblast M2 lower bound.
- P2 crowded direct point annotations were removed; marker shape identifies target and color identifies stratum through a shared legend. Every point and both marginal donor-bootstrap intervals remain visible.
- P3 initially used density bins. The full fibroblast residual extent compressed its central distribution. The final panel shows every eligible region pair with symmetric-log axes (linear within ±0.25) and full extrema retained. This is a display transformation only; every metric was calculated on the untransformed residuals. CI endpoints use three decimals so a small positive bound is not rounded to zero.
- P4 preserves a common x scale, displays all five controls in the requested order, and includes raw/Macenko/grayscale M4 arms. Both M3 sensitivities are retained in the result tables.
- P5 shows every donor, with deterministic ordering by ΔMAE; no donor was selected or omitted. Improvement counts moved into titles, and a missing proportional-sign glyph was replaced by plain language. Marker area is 12 + 14 log(1 + n regions).
- No plotted point estimate falls outside its percentile interval; all plotted CI endpoints agree with the donor bootstrap source data.

## Uncertainty and interpretation

Percentile intervals use 2,000 complete-donor resamples, keeping all donor sections/regions together. Paired increment intervals recompute M4–M2 inside each resample. Intervals condition on the completed OOF fits. P3 intervals summarize donor uncertainty in residual R²; region points are descriptive. P5 presents each donor directly and intentionally has no region-level bootstrap interval. No figure treats microregions as independent inferential replicates.

## Export audit

| Figure | SVG size | Editable text elements | Minimum text | PNG pixels |
|---|---|---:|---:|---|
| figure1_model_ladder | 517.566094pt × 253.952906pt | 34 | 6.0 pt | [2157, 1059] |
| figure2_morphology_increment | 517.733438pt × 295.76625pt | 20 | 6.0 pt | [2157, 1233] |
| figure3_crossfit_residual | 517.626094pt × 221.90174pt | 42 | 6.0 pt | [2158, 925] |
| figure4_shortcut_stress_test | 517.524453pt × 372.774375pt | 61 | 5.5 pt | [2156, 1553] |
| figure5_donor_heterogeneity | 517.853828pt × 415.32825pt | 83 | 5.5 pt | [2158, 1731] |

No unit tests were run. The input gate, evidence assertions, syntax checks and required figure QA were executed. The original frozen inputs and 26 CPU embedding caches were preserved.
