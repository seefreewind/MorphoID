#!/usr/bin/env python3
"""Run only the three frozen secondary-usable within-lineage programs after primary completion."""
from pathlib import Path
import phase1e_g_fit_models as pilot
from threadpoolctl import threadpool_limits

def main():
    if not (pilot.P1 / "primary_model_performance.tsv").is_file():
        raise RuntimeError("all primary analyses must complete before secondary programs")
    threadpool_limits(limits=4)
    pilot.input_gate()
    pilot.TARGETS = {
        "epithelial_transitional": {"y": "Y_epithelial__transitional_epithelial",
            "n": "n_target_cells__epithelial__transitional_epithelial", "lineage": "epithelial"},
        "fibroblast_myofibroblast": {"y": "Y_fibroblast__myofibroblast",
            "n": "n_target_cells__fibroblast__myofibroblast", "lineage": "fibroblast"},
        "macrophage_profibrotic": {"y": "Y_macrophage__profibrotic",
            "n": "n_target_cells__macrophage__profibrotic", "lineage": "macrophage/monocyte"}}
    pilot.OUT = pilot.P1 / "secondary"
    pilot.OUT.mkdir(exist_ok=True)
    outputs = pilot.primary_analysis()
    names = ["oof_predictions_secondary.tsv.gz", "ridge_alpha_selection.tsv", "analysis_dataset.tsv",
             "secondary_model_performance.tsv", "donor_bootstrap_metrics.tsv",
             "residual_prediction_performance.tsv", "donor_heterogeneity.tsv",
             "phase1g_train_oof_residuals.tsv.gz", "residual_fold_metrics.tsv"]
    for name, frame in zip(names, outputs):
        frame.to_csv(pilot.OUT / name, sep="\t", index=False,
                     compression="gzip" if name.endswith(".gz") else None)
    pilot.input_gate()
    print("PHASE1_SECONDARY_PASS: three frozen programs; separate secondary BH families", flush=True)

if __name__ == "__main__":
    main()
