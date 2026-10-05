#!/usr/bin/env python3
"""Extract permitted low-level H&E color summaries and audit shortcuts.

No texture, nuclei, deep embeddings, pathology embeddings, or expression data
are used. All outer evaluation folds hold out complete donors.
"""
from __future__ import annotations

import argparse
import gzip
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from skimage.color import rgb2hed
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / "results/phase0/phase0c_v"
VR = ROOT / "results/phase0/phase0c_vr"
OUT = ROOT / "results/phase0/phase0d"
RAW = ROOT / "data/raw/phase0c_v"
CONFIG = ROOT / "configs"
SCALE = 0.2125
GRID = 150.0
STRIDE = 8
SEED = 20260927
COLOR_FEATURES = [
    *[f"rgb_{c}_{s}" for c in "rgb" for s in ("mean", "sd")],
    *[f"hed_{c}_{s}" for c in ("hematoxylin", "eosin", "dab") for s in ("mean", "sd")],
    "brightness_mean", "brightness_sd", "saturation_mean", "saturation_sd",
    *[f"rgb_{c}_hist_{i}" for c in "rgb" for i in range(8)],
]


def load_image(path: Path) -> np.ndarray:
    with gzip.open(path, "rb") as fh:
        with tifffile.TiffFile(fh) as tif:
            arr = tif.pages[0].asarray()
    if arr.ndim == 2:
        arr = np.repeat(arr[..., None], 3, axis=2)
    if arr.ndim == 3 and arr.shape[0] in (3,4) and arr.shape[-1] not in (3,4):
        arr = np.moveaxis(arr, 0, -1)
    return arr[..., :3]


def aggregate_feature(values: np.ndarray, flat: np.ndarray, valid: np.ndarray, ngrid: int,
                      counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sums = np.bincount(flat[valid], weights=values[valid].astype(float), minlength=ngrid)
    sums2 = np.bincount(flat[valid], weights=np.square(values[valid].astype(float)), minlength=ngrid)
    mean = np.divide(sums, counts, out=np.full(ngrid, np.nan), where=counts > 0)
    var = np.divide(sums2, counts, out=np.full(ngrid, np.nan), where=counts > 0) - np.square(mean)
    return mean, np.sqrt(np.maximum(var, 0))


def extract_one(sample: str, image_file: str, micro: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    path = RAW / image_file
    if not path.exists():
        raise FileNotFoundError(path)
    arr = load_image(path)
    h, w = arr.shape[:2]
    sampled = arr[::STRIDE, ::STRIDE, :3]
    dtype_max = float(np.iinfo(arr.dtype).max) if np.issubdtype(arr.dtype, np.integer) else max(float(np.nanmax(sampled)), 1.0)
    rgb = np.clip(sampled.astype(np.float32) / dtype_max, 0, 1)
    nr, nc = rgb.shape[:2]
    yy, xx = np.indices((nr, nc))
    x_um = xx.ravel().astype(float) * STRIDE * SCALE
    y_um = yy.ravel().astype(float) * STRIDE * SCALE
    nx = int(math.ceil(w * SCALE / GRID)); ny = int(math.ceil(h * SCALE / GRID))
    gx = np.floor(x_um / GRID).astype(int); gy = np.floor(y_um / GRID).astype(int)
    inside = (gx < nx) & (gy < ny)
    flat = (gy[inside] * nx + gx[inside]).astype(int)
    flat_sample = sampled.reshape(-1, 3)[inside]
    rgb_sample = rgb.reshape(-1, 3)[inside]
    maxc = rgb_sample.max(axis=1)
    minc = rgb_sample.min(axis=1)
    black = maxc <= 0.02
    valid = ~black
    ngrid = nx * ny
    counts = np.bincount(flat[valid], minlength=ngrid)
    stats: dict[str, np.ndarray] = {}
    for j, channel in enumerate("rgb"):
        mean, sd = aggregate_feature(rgb_sample[:,j], flat, valid, ngrid, counts)
        stats[f"rgb_{channel}_mean"] = mean
        stats[f"rgb_{channel}_sd"] = sd
        bins = np.minimum((rgb_sample[:,j] * 8).astype(int), 7)
        hist = np.bincount(flat[valid] * 8 + bins[valid], minlength=ngrid*8).reshape(ngrid,8)
        hist = np.divide(hist, counts[:,None], out=np.full((ngrid,8),np.nan), where=counts[:,None]>0)
        for b in range(8):
            stats[f"rgb_{channel}_hist_{b}"] = hist[:,b]
    hed = rgb2hed(rgb.reshape(-1,3).reshape(nr,nc,3)).reshape(-1,3)[inside]
    for j, channel in enumerate(("hematoxylin", "eosin", "dab")):
        mean, sd = aggregate_feature(hed[:,j], flat, valid, ngrid, counts)
        stats[f"hed_{channel}_mean"] = mean
        stats[f"hed_{channel}_sd"] = sd
    brightness = maxc
    saturation = np.divide(maxc-minc, maxc, out=np.zeros_like(maxc), where=maxc>0)
    for name, values in (("brightness", brightness), ("saturation", saturation)):
        mean, sd = aggregate_feature(values, flat, valid, ngrid, counts)
        stats[f"{name}_mean"] = mean
        stats[f"{name}_sd"] = sd
    records = micro[micro.sample_id.eq(sample)].copy()
    grid_id = records.region_y.astype(int).to_numpy() * nx + records.region_x.astype(int).to_numpy()
    out = records[["sample_id", "region_x", "region_y"]].copy()
    for name in COLOR_FEATURES:
        out[name] = stats[name][grid_id]
    out["n_color_pixels_sampled_qc"] = counts[grid_id]
    out["black_fraction_sampled_qc"] = np.divide(
        np.bincount(flat[black], minlength=ngrid)[grid_id],
        np.maximum(np.bincount(flat, minlength=ngrid)[grid_id], 1))
    meta = {"sample_id": sample, "image_file": image_file, "width_px": w, "height_px": h,
            "sample_stride_px": STRIDE, "fixed_scale_um_per_px": SCALE,
            "color_feature_count": len(COLOR_FEATURES), "n_microregions": len(out),
            "regions_with_color_pixels": int(out.n_color_pixels_sampled_qc.gt(0).sum()),
            "image_sha_source": "registered-H&E SHA256 in VR author_sample_run_crosswalk.tsv"}
    return out, meta


def section_level_predictions(frame: pd.DataFrame, target: str, features: list[str], foldmap: pd.Series,
                              sample_weights: pd.Series | None = None, bootstrap_replicates: int = 1500) -> tuple[pd.DataFrame, dict]:
    # The labels (TMA/run/disease) are section-level. First average the allowed
    # 150-um region summaries within section, then fit one row per section so
    # adjacent regions cannot act as pseudo-replicates in the classifier.
    work = frame.dropna(subset=[target, *features]).copy()
    region_count = int(work.n_regions.sum()) if "n_regions" in work else len(work)
    work[target] = work[target].astype(str)
    work["fold"] = work.donor_id.map(foldmap)
    group_cols = ["sample_id", "donor_id", target, "fold"]
    sec = work.groupby(group_cols, as_index=False)[features].mean()
    if "n_regions" in work:
        nreg = work.groupby("sample_id").n_regions.sum()
        sec["n_regions"] = sec.sample_id.map(nreg).astype(int)
    sec["fold"] = sec.donor_id.map(foldmap)
    classes = sorted(sec[target].unique())
    X = sec[features].to_numpy(float); y = sec[target].to_numpy()
    oof = np.zeros((len(sec), len(classes)), dtype=float)
    pred = np.full(len(sec), "", dtype=object)
    unsupported = set()
    for fold in sorted(sec.fold.unique()):
        te = sec.fold.eq(fold).to_numpy(); tr = ~te
        train_classes = set(y[tr]); test_classes = set(y[te])
        unsupported.update(test_classes-train_classes)
        if len(train_classes) < 2:
            only = sorted(train_classes)[0]
            pred[te] = only; oof[np.ix_(te, np.flatnonzero(np.array(classes)==only))] = 1
            continue
        scaler = StandardScaler().fit(X[tr])
        model = LogisticRegression(C=1.0, max_iter=500, class_weight="balanced", random_state=SEED)
        nsec = sec.loc[tr].groupby("donor_id").sample_id.nunique()
        train_weights = sec.loc[tr].donor_id.map(lambda d: 1.0/nsec[d]).to_numpy(float)
        model.fit(scaler.transform(X[tr]), y[tr], sample_weight=train_weights)
        p = model.predict_proba(scaler.transform(X[te]))
        pred[te] = model.classes_[np.argmax(p,axis=1)]
        for j, cls in enumerate(model.classes_):
            oof[np.ix_(te, np.flatnonzero(np.array(classes)==cls))] = p[:,[j]]
    sec["oof_prediction"] = pred
    for j, cls in enumerate(classes):
        sec[f"p_{cls}"] = oof[:,j]
    # Each donor contributes equal score weight; donor blocks are bootstrap units.
    pcols = [f"p_{c}" for c in classes]
    sec["oof_prediction"] = sec[pcols].to_numpy().argmax(axis=1)
    sec["oof_prediction"] = [classes[i] for i in sec.oof_prediction]
    donor_sections = sec.groupby("donor_id").sample_id.nunique()
    weights = sec.donor_id.map(lambda d: 1.0/donor_sections[d]).to_numpy(float)
    ba = balanced_accuracy_score(sec[target], sec.oof_prediction, sample_weight=weights)
    try:
        if unsupported:
            raise ValueError("held-out class missing from at least one training fold")
        if len(classes) == 2:
            # sklearn's binary ROC-AUC API expects the score for one class;
            # passing the full two-column probability matrix is only valid for
            # a true multiclass call and otherwise leaves this metric blank.
            positive = classes[1]
            auc = roc_auc_score(
                sec[target].eq(positive).astype(int), sec[f"p_{positive}"],
                sample_weight=weights)
        else:
            auc = roc_auc_score(sec[target], sec[pcols].to_numpy(), labels=classes,
                                multi_class="ovr", average="macro", sample_weight=weights)
    except Exception:
        auc = np.nan
    rng = np.random.default_rng(SEED + 21)
    donor_ids = sec.donor_id.unique(); boot=[]
    for _ in range(bootstrap_replicates):
        chosen = rng.choice(donor_ids, len(donor_ids), replace=True)
        blocks = [np.flatnonzero(sec.donor_id.to_numpy()==d) for d in chosen]
        ix = np.concatenate(blocks)
        if sec[target].iloc[ix].nunique() < 2: continue
        ww = np.concatenate([np.full(len(block), 1/max(len(block),1)) for block in blocks])
        boot.append(balanced_accuracy_score(sec[target].iloc[ix],sec.oof_prediction.iloc[ix],sample_weight=ww))
    row = {"target": target, "n_regions": region_count, "n_sections": len(sec), "n_donors": sec.donor_id.nunique(),
           "n_classes": len(classes), "balanced_accuracy": float(ba), "macro_auroc": float(auc),
           "bootstrap_95ci_low": float(np.quantile(boot,.025)) if boot else np.nan,
           "bootstrap_95ci_high": float(np.quantile(boot,.975)) if boot else np.nan,
           "unseen_test_classes_across_folds": ";".join(sorted(unsupported)),
           "n_unseen_test_classes": len(unsupported), "outer_group": "donor_id",
           "model": "section mean of allowed 150-um region color features; standardized multinomial logistic, C=1, class_weight=balanced"}
    return sec, row


def block_permute_section_labels(sections: pd.DataFrame, target: str, rng: np.random.Generator) -> dict[str,str]:
    out: dict[str,str] = {}
    for nsec, block in sections.groupby("n_sections", sort=True):
        donor_ids = sorted(block.donor_id.astype(str).unique())
        patterns = {}
        for donor in donor_ids:
            d = sections[sections.donor_id.astype(str).eq(donor)].sort_values("sample_id")
            patterns[donor] = d[target].astype(str).tolist()
        source_donors = np.array(donor_ids, dtype=object)
        rng.shuffle(source_donors)
        for recipient, source in zip(donor_ids, source_donors):
            rec_samples = sections[sections.donor_id.astype(str).eq(recipient)].sort_values("sample_id").sample_id.astype(str).tolist()
            pattern = patterns[str(source)]
            for sample, label in zip(rec_samples, pattern):
                out[sample] = label
    return out


def evaluate_color_targets(regions: pd.DataFrame, features: list[str], folds: pd.DataFrame,
                           n_perm: int) -> tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    foldmap = folds.set_index("donor_id").fold
    regs = regions.copy()
    regs["TMA_label"] = regs.TMA.astype(str).map(lambda x: x if x.startswith("TMA") else "TMA"+x)
    regs["run_label"] = regs.run.astype(str).map(lambda x: x if x.startswith("Run") else "Run"+x)
    regs["disease_label"] = regs.disease_status.astype(str).str.lower().map(lambda x: "control" if x in {"control","unaffected"} else "pulmonary_fibrosis")
    sections = regs[["sample_id","donor_id","TMA_label","run_label","disease_label"]].drop_duplicates().copy()
    sections["n_sections"] = sections.groupby("donor_id").sample_id.transform("nunique")
    section_features = regs.groupby(
        ["sample_id", "donor_id", "TMA_label", "run_label", "disease_label"], as_index=False)[features].mean()
    section_features["n_regions"] = section_features.sample_id.map(regs.groupby("sample_id").size()).astype(int)
    summaries=[]; oof_all=[]; perm_rows=[]
    rng=np.random.default_rng(SEED+22)
    for target in ["TMA_label","run_label","disease_label"]:
        section_oof, observed = section_level_predictions(section_features,target,features,foldmap,None,1500)
        observed["target"] = target.removesuffix("_label")
        observed["class_counts_sections"] = ";".join(f"{k}:{v}" for k,v in sections[target].value_counts().sort_index().items())
        observed["color_only_risk"] = "HIGH" if observed["balanced_accuracy"] >= .70 or (np.isfinite(observed["macro_auroc"]) and observed["macro_auroc"] >= .75) else "NOT_HIGH_BY_PRESET_THRESHOLD"
        summaries.append(observed); oof_all.append(section_oof)
        # Donor-block label permutation within donor section-count strata.
        null=[]
        for b in range(n_perm):
            mapping=block_permute_section_labels(sections,target,rng)
            pframe=section_features.copy(); pframe[target]=pframe.sample_id.map(mapping)
            poof,_=section_level_predictions(pframe,target,features,foldmap,None,0)
            dct=poof.groupby("donor_id").size(); w=poof.donor_id.map(lambda d: 1.0/dct[d])
            null.append(float(balanced_accuracy_score(poof[target],poof.oof_prediction,sample_weight=w)))
        summaries[-1]["permutation_replicates"] = n_perm
        summaries[-1]["permutation_baseline_mean"] = float(np.mean(null))
        summaries[-1]["permutation_baseline_95pct"] = float(np.quantile(null,.95))
        summaries[-1]["permutation_p"] = (1+sum(v>=observed["balanced_accuracy"] for v in null))/(len(null)+1)
        perm_rows.append({"target": target.removesuffix("_label"), "replicates": n_perm,
                          "observed_balanced_accuracy": observed["balanced_accuracy"],
                          "null_mean": float(np.mean(null)), "null_95pct": float(np.quantile(null,.95)),
                          "permutation_p": summaries[-1]["permutation_p"], "permutation_unit": "whole donor block, stratified by number of sections"})
    # Donor class prediction has no evaluable held-out classes under donor CV.
    summaries.append({"target":"donor_id","n_regions":len(regs),"n_sections":regs.sample_id.nunique(),
                      "n_donors":regs.donor_id.nunique(),"n_classes":regs.donor_id.nunique(),
                      "balanced_accuracy":np.nan,"macro_auroc":np.nan,"bootstrap_95ci_low":np.nan,
                      "bootstrap_95ci_high":np.nan,"unseen_test_classes_across_folds":"all held-out donors",
                      "n_unseen_test_classes":regs.donor_id.nunique(),"outer_group":"donor_id",
                      "model":"NOT_RUN: donor-held-out classes are entirely unseen",
                      "color_only_risk":"NOT_ESTIMABLE_WITH_DONOR_HELD_OUT_SPLIT"})
    return pd.DataFrame(summaries),pd.concat(oof_all,ignore_index=True),pd.DataFrame(perm_rows)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--permutations",type=int,default=100)
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    cohort=pd.read_csv(CONFIG/"frozen_vannan_primary_cohort.tsv",sep="\t")
    cohort=cohort[cohort.primary_inclusion.eq("YES")]
    cross=pd.read_csv(VR/"author_sample_run_crosswalk.tsv",sep="\t")
    metadata=pd.read_csv(V/"vannan_sample_crosswalk.tsv",sep="\t")
    metadata=metadata[metadata.sample_id.isin(cohort.sample_id)].copy()
    micro=pd.read_csv(V/"microregion_feasibility.tsv",sep="\t",usecols=["sample_id","donor_id","region_x","region_y","grid_um","disease_status"])
    micro=micro[micro.sample_id.isin(cohort.sample_id) & micro.grid_um.eq(GRID)].copy()
    meta_cols=["sample_id","TMA","run","disease_status","affected_status"]
    regions=micro.merge(metadata[meta_cols],on="sample_id",suffixes=("_cells",""),validate="many_to_one")
    regions=regions.merge(cohort[["sample_id","donor_id"]],on="sample_id",suffixes=("","_cohort"),validate="many_to_one")
    # Keep the corrected donor label as the grouping unit.
    regions["donor_id"]=regions["donor_id_cohort"]
    del regions["donor_id_cohort"]
    fold=pd.read_csv(CONFIG/"phase1_donor_folds.tsv",sep="\t")
    feature_parts=[]; manifests=[]
    cache_dir = OUT / "d5_color_region_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    for i,row in enumerate(cohort.sort_values("sample_id").itertuples(index=False),1):
        match=cross[cross.sample_id.eq(row.sample_id)].iloc[0]
        cache_path = cache_dir / f"{row.sample_id}.tsv.gz"
        if cache_path.exists():
            feat = pd.read_csv(cache_path, sep="\t")
            manifest = {"sample_id": row.sample_id, "image_file": match.registered_HE_file,
                        "registered_HE_sha256": match.registered_HE_sha256,
                        "fixed_scale_um_per_px": SCALE, "sample_stride_px": STRIDE,
                        "color_feature_count": len(COLOR_FEATURES), "n_microregions": len(feat),
                        "cache_reused": "YES"}
        else:
            feat,manifest=extract_one(row.sample_id,match.registered_HE_file,regions)
            manifest["registered_HE_sha256"] = match.registered_HE_sha256
            manifest["cache_reused"] = "NO"
            feat.to_csv(cache_path,sep="\t",index=False,compression="gzip")
        feature_parts.append(feat); manifests.append(manifest)
        action = "loaded cached color features" if manifest["cache_reused"] == "YES" else "color features extracted"
        print(f"[{i}/{len(cohort)}] {action}: {row.sample_id} ({len(feat)} regions)",flush=True)
    color=pd.concat(feature_parts,ignore_index=True)
    color=color.merge(regions[["sample_id","region_x","region_y","TMA","run","disease_status","affected_status","donor_id"]].drop_duplicates(),
                      on=["sample_id","region_x","region_y"],how="left",validate="one_to_one")
    feat_path=OUT/"registered_he_region_color_features.tsv.gz"
    color.to_csv(feat_path,sep="\t",index=False,compression="gzip")
    pd.DataFrame(manifests).to_csv(OUT/"d5_color_feature_extraction_manifest.tsv",sep="\t",index=False)
    pd.DataFrame([{"feature_set":"D5_color_only","included_features":";".join(COLOR_FEATURES),
                   "feature_count":len(COLOR_FEATURES),"excluded_features":"texture;nuclei;deep embeddings;pathology embeddings;expression",
            "sample_stride_px":STRIDE,"grid_um":GRID,"scale_um_per_px":SCALE,
            "classifier_input": "mean of 150-um region color summaries within section; one row per section",
                   "excluded_pure_black_pixels_from_color_statistics":"YES","donor_grouped_outer_folds":5}]).to_csv(
        OUT/"d5_color_feature_specification.tsv",sep="\t",index=False)
    metrics_df,oof,perm=evaluate_color_targets(color,COLOR_FEATURES,fold,args.permutations)
    metrics_df.to_csv(OUT/"color_shortcut_metrics.tsv",sep="\t",index=False)
    oof.to_csv(OUT/"color_shortcut_oof_predictions.tsv",sep="\t",index=False)
    perm.to_csv(OUT/"color_shortcut_permutation_baseline.tsv",sep="\t",index=False)
    print(f"Color-only audit complete: {len(color)} regions; {args.permutations} donor-block permutations per target.",flush=True)


if __name__=="__main__":
    main()
