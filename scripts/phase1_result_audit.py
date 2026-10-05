#!/usr/bin/env python3
"""Required Phase 1 evidence QA; no model fitting or input mutation."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from phase1e_g_fit_models import input_gate, donor_cluster_bootstrap, r2_score, SEED

ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "results/phase1"
def read(name):
    return pd.read_csv(P1 / name, sep="\t")

def main():
    input_gate()
    oof = read("oof_predictions_primary.tsv.gz")
    common = read("phase1d_r2/final_common_phase1_region_universe.tsv")
    folds = pd.read_csv(ROOT / "configs/phase1_donor_folds.tsv", sep="\t")
    fmap = folds.set_index("donor_id").fold.to_dict()
    evidence = {}
    def claim(name, condition):
        evidence[name] = bool(condition)
        if not condition:
            raise RuntimeError("EVIDENCE_QA_FAILURE: " + name)
    claim("common_membership", set(oof.region_id).issubset(set(common.region_id)))
    claim("donor_frozen_folds", oof.fold.eq(oof.donor_id.map(fmap)).all())
    claim("no_duplicate_oof", not oof.duplicated(["region_id", "target", "stratum"]).any())
    predcols = [x for x in oof if x.endswith("_pred")] + ["heldout_base_residual", "morphology_predicted_residual"]
    claim("all_arms_same_rows_finite", np.isfinite(oof[predcols].to_numpy(float)).all())
    claim("exact_three_primary_three_strata", oof.target.nunique() == 3 and
          set(oof.stratum) == {"POOLED", "PF_ONLY", "DISEASE_ADJUSTED"})
    residual = read("phase1g_train_oof_residuals.tsv")
    claim("train_residual_outer_test_excluded", residual.outer_fold.ne(residual.fold).all())
    claim("train_residual_donor_crossfit", residual.train_side_oof_fold.eq(residual.fold).all())
    claim("residual_identity", np.allclose(oof.heldout_base_residual, oof.observed-oof.M2_pred))
    claim("training_residual_identity", np.allclose(residual.M2_train_oof_residual,
                                                   residual.observed-residual.M2_train_oof_pred))
    contract = json.loads((P1 / "phase1_execution_contract.json").read_text())
    import hashlib
    for path, expected_sha in contract["ancillary_input_sha256"].items():
        claim("unchanged_" + path, hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected_sha)
    reference = read("macenko_fold_target_parameters.tsv")
    for row in reference.itertuples():
        training = set(str(row.training_donors).split(";"))
        heldout = set(str(row.heldout_donors).split(";"))
        claim("macenko_" + row.outer_fold + "_outer_test_excluded", not training.intersection(heldout)
              and all(fmap[d] == row.outer_fold for d in heldout))
    alphas = read("ridge_alpha_selection.tsv")
    claim("tuning_in_frozen_grid", np.isclose(alphas.alpha.to_numpy(float)[:,None], np.logspace(-4,5,19)[None,:], rtol=1e-12, atol=0).any(axis=1).all())
    claim("complete_45_outer_folds", len([p for p in (P1/"model_fold_checkpoints").glob("*__*.tsv") if not p.name.startswith("._")]) == 45)
    expected = {"epithelial_injury":3219, "fibroblast_activation":4869,
                "macrophage_inflammatory":2661}
    m = read("primary_model_performance.tsv")
    for target, n in expected.items():
        for stratum in ["POOLED", "DISEASE_ADJUSTED"]:
            g = oof[(oof.target==target)&(oof.stratum==stratum)]
            claim(target+"_"+stratum+"_coverage", len(g)==n)
            row = m[(m.target==target)&(m.stratum==stratum)].iloc[0]
            claim(target+"_"+stratum+"_paired_delta", np.isclose(
                row.deltaR2_morph, r2_score(g.observed, g.M4_pred)-r2_score(g.observed,g.M2_pred)))
    # Paired sensitivity CIs are computed on the same completed OOF predictions.
    shortcut = read("shortcut_robustness.tsv")
    for j, ((target,stratum), g) in enumerate(oof.groupby(["target","stratum"], sort=False)):
        for arm, col in [("macenko","macenko_pred"), ("grayscale","grayscale_pred")]:
            b = donor_cluster_bootstrap(g.observed.to_numpy(float),
                {"M4":g[col].to_numpy(float),"M2":g.M2_pred.to_numpy(float)},
                g.donor_id.to_numpy(str), SEED+6000+j)["deltaR2_morphology"]
            sel=(shortcut.target==target)&(shortcut.stratum==stratum)
            shortcut.loc[sel,arm+"_deltaR2_CI_low"]=b["CI_low"]
            shortcut.loc[sel,arm+"_deltaR2_CI_high"]=b["CI_high"]
    shortcut.to_csv(P1/"shortcut_robustness.tsv",sep="\t",index=False)
    (P1/"phase1_result_audit.json").write_text(json.dumps(evidence,indent=2)+"\n")
    print("PHASE1_RESULT_AUDIT_PASS: "+str(len(evidence))+" evidence assertions",flush=True)

if __name__ == "__main__":
    main()
