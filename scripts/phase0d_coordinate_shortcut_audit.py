#!/usr/bin/env python3
"""Donor-held-out, coordinate-only shortcut audit (no expression input)."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / "results/phase0/phase0c_v"
VR = ROOT / "results/phase0/phase0c_vr"
OUT = ROOT / "results/phase0/phase0d"
CONFIG = ROOT / "configs"
SEED = 20260927
SCALE = 0.2125
GRID = 150


def metrics(y, pred, prob, classes, weights):
    ba = balanced_accuracy_score(y, pred, sample_weight=weights)
    try:
        if len(classes) == 2:
            auc = roc_auc_score(y, prob[:, 1], labels=classes, sample_weight=weights)
        else:
            auc = roc_auc_score(y, prob, labels=classes, multi_class="ovr", average="macro", sample_weight=weights)
    except Exception:
        auc = np.nan
    return float(ba), float(auc)


def donor_block_bootstrap_ba(y, pred, donors, n_boot=1500, seed=SEED + 20):
    """Fast donor-block bootstrap of weighted balanced accuracy.

    Each donor contributes total weight one, matching the region-level
    inverse-donor-count weights used for the point estimate. Precomputed donor
    confusion matrices make resampling equivalent without rescoring all
    regions or recomputing AUROC in every replicate.
    """
    y = np.asarray(y).astype(str)
    pred = np.asarray(pred).astype(str)
    donors = np.asarray(donors).astype(str)
    classes = np.array(sorted(np.unique(y)))
    donor_ids = np.array(sorted(np.unique(donors)))
    donor_cms = []
    for donor in donor_ids:
        mask = donors == donor
        weights = np.full(mask.sum(), 1.0 / mask.sum())
        donor_cms.append(confusion_matrix(y[mask], pred[mask], labels=classes,
                                          sample_weight=weights))
    donor_cms = np.asarray(donor_cms, dtype=float)
    rng = np.random.default_rng(seed)
    choices = rng.integers(0, len(donor_ids), size=(n_boot, len(donor_ids)))
    cms = donor_cms[choices].sum(axis=1)
    denom = cms.sum(axis=2)
    recalls = np.divide(np.diagonal(cms, axis1=1, axis2=2), denom,
                        out=np.full((n_boot, len(classes)), np.nan), where=denom > 0)
    with np.errstate(invalid="ignore"):
        ba = np.nanmean(recalls, axis=1)
    return ba[np.isfinite(ba)].tolist()


def grouped_oof(frame, target, features, folds):
    classes = np.array(sorted(frame[target].dropna().astype(str).unique()))
    work = frame.dropna(subset=[target, *features]).copy().reset_index(drop=True)
    donor_fold = folds.set_index("donor_id").fold
    work["fold"] = work.donor_id.map(donor_fold)
    X = work[features].to_numpy(float)
    y = work[target].astype(str).to_numpy()
    pred = np.full(len(work), "", dtype=object)
    prob = np.zeros((len(work), len(classes)), dtype=float)
    unsupported = set()
    donor_region_n = work.groupby("donor_id").size()
    sample_weight = work.donor_id.map(lambda d: 1.0/donor_region_n[d]).to_numpy(float)
    for fold in sorted(work.fold.unique()):
        te = work.fold.eq(fold).to_numpy(); tr = ~te
        train_classes = set(y[tr]); test_classes = set(y[te])
        unsupported.update(test_classes - train_classes)
        if len(train_classes) < 2:
            pred[te] = next(iter(train_classes))
            prob[np.ix_(te, np.flatnonzero(classes == pred[te][0]))] = 1.0
            continue
        scaler = StandardScaler().fit(X[tr])
        model = LogisticRegression(C=1.0, max_iter=500, class_weight="balanced", random_state=SEED)
        model.fit(scaler.transform(X[tr]), y[tr], sample_weight=sample_weight[tr])
        pp = model.predict_proba(scaler.transform(X[te]))
        pred[te] = model.classes_[np.argmax(pp, axis=1)]
        for j, cls in enumerate(model.classes_):
            prob[np.ix_(te, np.flatnonzero(classes == cls))] = pp[:, [j]]
    counts = work.groupby("donor_id").size()
    weights = work.donor_id.map(lambda d: 1.0/counts[d]).to_numpy(float)
    ba, auc = metrics(y, pred, prob, classes, weights)
    if unsupported:
        auc = np.nan
    boot = donor_block_bootstrap_ba(y, pred, work.donor_id.to_numpy())
    res = {"target": target, "n_regions": len(work), "n_donors": work.donor_id.nunique(),
           "n_classes": len(classes), "balanced_accuracy": ba, "macro_auroc": auc,
           "bootstrap_95ci_low": float(np.quantile(boot,.025)) if boot else np.nan,
           "bootstrap_95ci_high": float(np.quantile(boot,.975)) if boot else np.nan,
           "unseen_test_classes": ";".join(sorted(unsupported)),
           "folds": 5, "group": "donor_id", "features": ";".join(features)}
    predictions = work[["sample_id", "donor_id", "region_x", "region_y", target, "fold"]].copy()
    predictions["oof_prediction"] = pred
    for j, cls in enumerate(classes):
        predictions[f"p_{cls}"] = prob[:,j]
    return res, predictions


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cohort = pd.read_csv(CONFIG / "frozen_vannan_primary_cohort.tsv", sep="\t")
    cohort = cohort[cohort.primary_inclusion.eq("YES")]
    folds = pd.read_csv(CONFIG / "phase1_donor_folds.tsv", sep="\t")
    grid_meta = pd.read_csv(VR / "registered_he_150um_grid_audit.tsv", sep="\t")
    grid_meta = grid_meta.set_index("sample_id")
    usecols = ["sample_id", "donor_id", "region_x", "region_y", "grid_um", "total_cells", "disease_status", "affected_status"]
    regions = pd.read_csv(V / "microregion_feasibility.tsv", sep="\t", usecols=usecols)
    regions = regions[regions.sample_id.isin(cohort.sample_id) & regions.grid_um.eq(GRID)].copy()
    sec = cohort.set_index("sample_id")
    regions["TMA"] = regions.sample_id.map(sec.TMA).astype(str).map(lambda x: x if x.startswith("TMA") else "TMA"+x)
    regions["disease"] = regions.sample_id.map(sec.disease_status).astype(str).str.lower().map(
        lambda x: "control" if x in {"control", "unaffected"} else "pulmonary_fibrosis")
    regions["affected"] = regions.sample_id.map(sec.affected_status).astype(str).str.lower()
    regions["run"] = regions.sample_id.map(sec.run).astype(str).map(lambda x: x if x.startswith("Run") else "Run"+x)
    regions["x_um"] = (regions.region_x + .5) * GRID
    regions["y_um"] = (regions.region_y + .5) * GRID
    dims = grid_meta.loc[regions.sample_id]
    regions["width_um"] = dims.width_px.to_numpy(float) * SCALE
    regions["height_um"] = dims.height_px.to_numpy(float) * SCALE
    regions["x_norm"] = regions.x_um / regions.width_um
    regions["y_norm"] = regions.y_um / regions.height_um
    features = ["x_um", "y_um", "x_norm", "y_norm", "width_um", "height_um"]
    results, preds = [], []
    for target in ["TMA", "run", "disease", "affected"]:
        record, prediction = grouped_oof(regions, target, features, folds)
        results.append(record); preds.append(prediction)

    # Assign pathology region category from author cell-region assignments,
    # joined by exact cell ID to corrected-RDS Xenium coordinates. Use one
    # dominant category per 150-um region with >=10 author-assigned cells.
    cross = pd.read_csv(VR / "author_sample_run_crosswalk.tsv", sep="\t")
    alias_map = cross[cross.sample_id.isin(cohort.sample_id)][["sample_id", "author_id_list_alias"]]
    assignments_path = VR / "author_annotations/GSE250346_HE_annotations/cells_partitioned_by_annotation.csv"
    assignments = pd.read_csv(assignments_path, usecols=["sample", "cell_id", "annotation_type_instance"])
    assignments = assignments.merge(alias_map, left_on="sample", right_on="author_id_list_alias", how="inner", validate="many_to_one")
    cell_meta = pd.read_csv(V / "vannan_primary_cell_metadata.tsv.gz", sep="\t",
                            usecols=["sample_id", "cell_id_local", "x_um", "y_um"])
    cell_meta = cell_meta[cell_meta.sample_id.isin(cohort.sample_id)]
    joined = assignments.merge(cell_meta, left_on=["sample_id", "cell_id"],
                              right_on=["sample_id", "cell_id_local"], how="left", validate="many_to_one")
    joined = joined.dropna(subset=["x_um", "y_um"])
    joined["region_x"] = np.floor(joined.x_um / GRID).astype(int)
    joined["region_y"] = np.floor(joined.y_um / GRID).astype(int)
    joined["category"] = joined.annotation_type_instance.astype(str).str.replace(r"_\d+$", "", regex=True)
    counts = joined.groupby(["sample_id", "region_x", "region_y", "category"], as_index=False).size()
    totals = counts.groupby(["sample_id", "region_x", "region_y"]).size().rename("n_categories")
    counts = counts.merge(totals, on=["sample_id", "region_x", "region_y"])
    dominant = counts.sort_values(["sample_id", "region_x", "region_y", "size", "category"], ascending=[True,True,True,False,True])
    dominant = dominant.groupby(["sample_id", "region_x", "region_y"], as_index=False).first()
    # Exclude categories that cannot be represented in at least two distinct
    # outer folds; unseen classes are not scored as if they were predictable.
    cat_donors = dominant.groupby("category").donor_id.nunique() if "donor_id" in dominant else None
    # Attach donor/section metadata and count assigned cells per winning class.
    donor_by_sample = cohort.set_index("sample_id").donor_id
    dominant["donor_id"] = dominant.sample_id.map(donor_by_sample)
    dominant["assigned_cells_in_region"] = dominant["size"]
    category_support = dominant.groupby("category").agg(n_donors=("donor_id", "nunique"), n_regions=("sample_id", "size"))
    donor_fold = folds.set_index("donor_id").fold
    category_fold_count = dominant.assign(fold=dominant.donor_id.map(donor_fold)).groupby("category").fold.nunique()
    keep = category_support.index[(category_support.n_donors >= 3) & category_fold_count.reindex(category_support.index).fillna(0).ge(2)]
    dominant = dominant[dominant.category.isin(keep) & dominant.assigned_cells_in_region.ge(10)].copy()
    if not dominant.empty:
        dom_dims = grid_meta.loc[dominant.sample_id]
        dominant["x_um"] = (dominant.region_x + .5) * GRID
        dominant["y_um"] = (dominant.region_y + .5) * GRID
        dominant["width_um"] = dom_dims.width_px.to_numpy(float) * SCALE
        dominant["height_um"] = dom_dims.height_px.to_numpy(float) * SCALE
        dominant["x_norm"] = dominant.x_um/dominant.width_um
        dominant["y_norm"] = dominant.y_um/dominant.height_um
        # Add section-level labels.
        for col, source_col in [("TMA","TMA"),("disease","disease_status"),("affected","affected_status"),("run","run")]:
            if col == "TMA":
                dominant[col] = dominant.sample_id.map(sec.TMA).astype(str).map(lambda x: x if x.startswith("TMA") else "TMA"+x)
            elif col == "run":
                dominant[col] = dominant.sample_id.map(sec.run).astype(str).map(lambda x: x if x.startswith("Run") else "Run"+x)
            elif col == "disease":
                dominant[col] = dominant.sample_id.map(sec.disease_status).astype(str).str.lower().map(lambda x: "control" if x in {"control","unaffected"} else "pulmonary_fibrosis")
            else:
                dominant[col] = dominant.sample_id.map(sec.affected_status).astype(str).str.lower()
        record, prediction = grouped_oof(dominant, "category", features, folds)
        record["pathology_categories_scored"] = ";".join(sorted(keep))
        record["excluded_rare_category_count"] = int((~category_support.index.isin(keep)).sum())
        results.append(record); preds.append(prediction)
        category_support.reset_index().to_csv(OUT / "author_pathology_region_category_support.tsv", sep="\t", index=False)
        dominant.to_csv(OUT / "coordinate_pathology_region_labels.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(results).to_csv(OUT / "coordinate_shortcut_metrics.tsv", sep="\t", index=False)
    if preds:
        pd.concat(preds, ignore_index=True, sort=False).to_csv(OUT / "coordinate_shortcut_oof_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
    print(f"Coordinate-only audit complete: {len(regions)} regions; path classes={len(keep) if 'keep' in locals() else 0}", flush=True)


if __name__ == "__main__":
    main()
