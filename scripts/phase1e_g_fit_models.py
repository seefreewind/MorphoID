#!/usr/bin/env python3
"""Fit the frozen donor-held-out Phase 1 minimum identifiability pilot."""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "results/phase1"
MATRIX = P1 / "region_analysis_matrix.tsv"
FOLDS = ROOT / "configs/phase1_donor_folds.tsv"
COHORT = ROOT / "configs/frozen_vannan_primary_cohort.tsv"
COLOR_FEATURES = ROOT / "results/phase0/phase0d/registered_he_region_color_features.tsv.gz"
COLOR_SPEC = ROOT / "results/phase0/phase0d/d5_color_feature_specification.tsv"
EMBEDDING_DIR = P1 / "embedding_cache"
OUT = P1
N_BOOT = 2000
SEED = 20260927
ALPHAS = np.logspace(-4, 5, 19)
INNER_MAX_SPLITS = 5
TARGETS = {
    "epithelial_injury": {
        "y": "Y_epithelial__injury",
        "n": "n_target_cells__epithelial__injury",
        "lineage": "epithelial",
    },
    "fibroblast_activation": {
        "y": "Y_fibroblast__fibroblast_activation",
        "n": "n_target_cells__fibroblast__fibroblast_activation",
        "lineage": "fibroblast",
    },
    "macrophage_inflammatory": {
        "y": "Y_macrophage__inflammatory",
        "n": "n_target_cells__macrophage__inflammatory",
        "lineage": "macrophage/monocyte",
    },
}
MAIN_STRATA = ["POOLED", "PF_ONLY", "DISEASE_ADJUSTED"]
ALL_STRATA = MAIN_STRATA  # Minimum mandatory pilot; optional control-only not executed.


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def input_gate() -> None:
    """Verify scientific inputs against the pre-modeling R2 snapshot."""
    import yaml
    cfg = yaml.safe_load((ROOT / "configs/phase1_common_universe.yaml").read_text())
    baseline = read_tsv(P1 / "phase1d_r2/scientific_freeze_assertions.tsv")
    checks = []
    def check(path, expected):
        path = Path(path)
        if not path.is_absolute():
            path = ROOT / path
        h = hashlib.sha256()
        if path.is_file():
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                    h.update(block)
            observed = h.hexdigest()
        else:
            observed = "MISSING"
        checks.append({"path": str(path), "expected_sha256": expected,
                       "observed_sha256": observed, "pass": observed == expected})
    for row in baseline.itertuples():
        # Scripts and narrative status are not frozen scientific inputs.
        if str(row.path).startswith(("configs/", "results/")) and not str(row.path).endswith(".py"):
            check(row.path, row.after_sha256)
    check(cfg["common_universe_file"], cfg["common_universe_sha256"])
    check(cfg["original_D5_path"], cfg["original_D5_sha256"])
    check("results/phase1/model_cache/owkin_phikon-v2/model.safetensors",
          "261ae680fa699b3b951597fd57aa19c02ef735805acb104b93af69b36d928569")
    gate = pd.DataFrame(checks)
    gate.to_csv(P1 / "phase1e_input_hash_gate.tsv", sep="\t", index=False)
    if not gate["pass"].all():
        raise RuntimeError("TECHNICAL_HOLD_PHASE1E_INPUT_MISMATCH")
    print(f"PHASE1E_INPUT_GATE_PASS: {len(gate)} checks", flush=True)


def donor_weights(donors: np.ndarray) -> np.ndarray:
    donors = np.asarray(donors).astype(str)
    counts = pd.Series(donors).value_counts().to_dict()
    return np.asarray([1.0 / counts[d] for d in donors], dtype=np.float64)


class WeightedRidge:
    """Standardized donor-equal Ridge solved by its weighted normal equations."""

    def __init__(self, scaler: StandardScaler, coefficient: np.ndarray, intercept: float):
        self.scaler = scaler
        self.coefficient = coefficient
        self.intercept = intercept

    def predict(self, x: np.ndarray) -> np.ndarray:
        return self.intercept + self.scaler.transform(np.asarray(x, dtype=np.float64)) @ self.coefficient


def fit_ridge(x: np.ndarray, y: np.ndarray, donors: np.ndarray, alpha: float) -> WeightedRidge:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    w = donor_weights(donors)
    scaler = StandardScaler().fit(x, sample_weight=w)
    xs = scaler.transform(x)
    ybar = float(np.average(y, weights=w))
    yw = y - ybar
    gram = xs.T @ (xs * w[:, None])
    rhs = xs.T @ (yw * w)
    coefficient = np.linalg.solve(gram + float(alpha) * np.eye(gram.shape[0]), rhs)
    return WeightedRidge(scaler, coefficient, ybar)


def donor_equal_mse(y: np.ndarray, pred: np.ndarray, donors: np.ndarray) -> float:
    _, inverse = np.unique(np.asarray(donors).astype(str), return_inverse=True)
    sqerr = (np.asarray(y, dtype=np.float64) - np.asarray(pred, dtype=np.float64)) ** 2
    sums = np.bincount(inverse, weights=sqerr)
    counts = np.bincount(inverse)
    return float(np.mean(sums / counts))


def choose_alpha(x: np.ndarray, y: np.ndarray, donors: np.ndarray) -> float:
    """Choose alpha by nested donor-grouped CV and donor-equal validation MSE."""
    donors = np.asarray(donors).astype(str)
    n_donors = len(np.unique(donors))
    if n_donors < 3:
        raise RuntimeError(f"Ridge tuning requires >=3 training donors; found {n_donors}")
    splitter = GroupKFold(n_splits=min(INNER_MAX_SPLITS, n_donors))
    splits = list(splitter.split(x, y, groups=donors))
    losses = np.zeros(len(ALPHAS), dtype=np.float64)
    n_validated_donors = np.zeros(len(ALPHAS), dtype=np.int64)
    for tr, va in splits:
        xv = np.asarray(x[va], dtype=np.float64)
        yv = np.asarray(y[va], dtype=np.float64)
        dtr = donors[tr]
        w = donor_weights(dtr)
        scaler = StandardScaler().fit(np.asarray(x[tr], dtype=np.float64), sample_weight=w)
        xt = scaler.transform(np.asarray(x[tr], dtype=np.float64))
        xval = scaler.transform(xv)
        ybar = float(np.average(y[tr], weights=w))
        yc = np.asarray(y[tr], dtype=np.float64) - ybar
        gram = xt.T @ (xt * w[:, None])
        rhs = xt.T @ (yc * w)
        eigval, eigvec = np.linalg.eigh(gram)
        projected = eigvec.T @ rhs
        for ai, alpha in enumerate(ALPHAS):
            beta = eigvec @ (projected / (eigval + alpha))
            pred = ybar + xval @ beta
            # Each donor contributes one mean squared error regardless of region count.
            n_va_donors = len(np.unique(donors[va]))
            val_loss = donor_equal_mse(yv, pred, donors[va])
            losses[ai] += val_loss * n_va_donors
            n_validated_donors[ai] += n_va_donors
    mean_loss = losses / np.maximum(n_validated_donors, 1)
    return float(ALPHAS[int(np.argmin(mean_loss))])


def r2_score(y: np.ndarray, pred: np.ndarray) -> float:
    y = np.asarray(y, dtype=np.float64)
    pred = np.asarray(pred, dtype=np.float64)
    denom = float(np.sum((y - y.mean()) ** 2))
    return float(1.0 - np.sum((y - pred) ** 2) / denom) if denom > 0 else np.nan


def regression_metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    y = np.asarray(y, dtype=np.float64)
    pred = np.asarray(pred, dtype=np.float64)
    from scipy.stats import pearsonr, spearmanr
    if len(y) < 2 or np.std(y) == 0 or np.std(pred) == 0:
        pearson, spearman = np.nan, np.nan
    else:
        pearson = float(pearsonr(y, pred).statistic)
        spearman = float(spearmanr(y, pred).statistic)
    err = y - pred
    return {"R2": r2_score(y, pred), "Pearson_r": pearson, "Spearman_rho": spearman,
            "MAE": float(np.mean(np.abs(err))), "RMSE": float(np.sqrt(np.mean(err ** 2)))}


def donor_cluster_bootstrap(y: np.ndarray, predictions: dict[str, np.ndarray], donors: np.ndarray,
                            seed: int, n_boot: int = N_BOOT) -> dict[str, dict[str, float]]:
    """Resample complete donors; compute model R2 and paired, direct DeltaR2."""
    donors = np.asarray(donors).astype(str)
    donor_ids = np.asarray(sorted(np.unique(donors)))
    rng = np.random.default_rng(seed)
    n_by_donor = np.asarray([np.count_nonzero(donors == d) for d in donor_ids], dtype=np.float64)
    sum_y_by_donor = np.asarray([np.sum(y[donors == d]) for d in donor_ids], dtype=np.float64)
    sum_y2_by_donor = np.asarray([np.sum(np.asarray(y[donors == d], dtype=np.float64) ** 2)
                                  for d in donor_ids], dtype=np.float64)
    sse_by_model = {
        name: np.asarray([np.sum((y[donors == d] - pred[donors == d]) ** 2)
                          for d in donor_ids], dtype=np.float64)
        for name, pred in predictions.items()
    }
    chosen = rng.integers(0, len(donor_ids), size=(n_boot, len(donor_ids)))
    total_n = n_by_donor[chosen].sum(axis=1)
    total_y = sum_y_by_donor[chosen].sum(axis=1)
    total_y2 = sum_y2_by_donor[chosen].sum(axis=1)
    sst = total_y2 - total_y ** 2 / total_n
    estimates: dict[str, np.ndarray] = {}
    for name, donor_sse in sse_by_model.items():
        sampled_sse = donor_sse[chosen].sum(axis=1)
        estimates[name] = np.divide(sst - sampled_sse, sst,
                                    out=np.full(n_boot, np.nan, dtype=np.float64), where=sst > 0)
    delta_pairs = {
        "deltaR2_composition": ("M1", "M0"),
        "deltaR2_neighborhood": ("M2", "M1"),
        "deltaR2_morphology": ("M4", "M2"),
    }
    deltas: dict[str, np.ndarray] = {
        name: estimates[upper] - estimates[lower]
        for name, (upper, lower) in delta_pairs.items()
        if upper in estimates and lower in estimates
    }
    output: dict[str, dict[str, float]] = {}
    for model, values in estimates.items():
        values = values[np.isfinite(values)]
        q = np.quantile(values, [0.025, 0.975])
        tail = min((np.count_nonzero(values <= 0) + 1) / (len(values) + 1),
                   (np.count_nonzero(values >= 0) + 1) / (len(values) + 1))
        output[model] = {"CI_low": float(q[0]), "CI_high": float(q[1]),
                         "bootstrap_two_sided_sign_p": float(min(1, 2 * tail))}
    for name, samples in deltas.items():
        samples = samples[np.isfinite(samples)]
        q = np.quantile(samples, [0.025, 0.975])
        arr = samples
        tail = min((np.count_nonzero(arr <= 0) + 1) / (len(arr) + 1),
                   (np.count_nonzero(arr >= 0) + 1) / (len(arr) + 1))
        output[name] = {
            "CI_low": float(q[0]), "CI_high": float(q[1]),
            "bootstrap_two_sided_sign_p": float(min(1.0, 2.0 * tail)),
        }
    return output


def bh_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan)
    good = np.isfinite(p)
    if not np.any(good):
        return q.tolist()
    vals = p[good]
    order = np.argsort(vals)
    ranked = vals[order] * len(vals) / (np.arange(len(vals)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    q_good = np.empty_like(vals)
    q_good[order] = ranked
    q[good] = q_good
    return q.tolist()


def load_embeddings() -> tuple[dict[str, np.ndarray], dict[str, np.ndarray],
                                dict[str, dict[str, np.ndarray]]]:
    if not EMBEDDING_DIR.is_dir():
        raise FileNotFoundError("Phase 1D embedding cache is missing")
    raw: dict[str, np.ndarray] = {}
    gray: dict[str, np.ndarray] = {}
    mac: dict[str, dict[str, np.ndarray]] = defaultdict(dict)
    files = sorted(p for p in EMBEDDING_DIR.glob("*.npz") if not p.name.startswith("._"))
    if len(files) != 26:
        raise RuntimeError(f"expected 26 section caches, found {len(files)}")
    for path in files:
        with np.load(path, allow_pickle=False) as z:
            if str(z["model_revision"].item()) != "2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd":
                raise RuntimeError(f"checkpoint revision mismatch in {path}")
            if str(z["inference_device"].item()) != "cpu":
                raise RuntimeError(f"backend mismatch in {path}")
            ids = z["region_ids"].astype(str)
            for i, rid in enumerate(ids):
                raw[rid] = z["raw"][i].astype(np.float32, copy=True)
                gray[rid] = z["grayscale"][i].astype(np.float32, copy=True)
                for key in z.files:
                    if key.startswith("macenko_fold_"):
                        fold = key.removeprefix("macenko_")
                        mac[fold][rid] = z[key][i].astype(np.float32, copy=True)
    return raw, gray, mac


def load_inputs() -> tuple[pd.DataFrame, list[str], dict[str, np.ndarray],
                            dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]:
    matrix = read_tsv(MATRIX)
    common = read_tsv(P1 / "phase1d_r2/final_common_phase1_region_universe.tsv")
    if common.region_id.duplicated().any() or len(common) != 12144:
        raise RuntimeError("invalid frozen common universe")
    if not set(common.region_id).issubset(set(matrix.region_id)):
        raise RuntimeError("common universe missing from matrix")
    matrix = matrix[matrix.region_id.isin(common.region_id)].copy()
    if len(matrix) != len(common):
        raise RuntimeError("common universe matrix cardinality mismatch")
    print(f"COMMON_UNIVERSE_LOADED: regions={len(matrix)}", flush=True)
    crop = read_tsv(P1 / "image_region_crop_audit.tsv")
    crop = crop[["region_id", "full_crop_in_bounds"]]
    matrix = matrix.merge(crop, on="region_id", validate="one_to_one",
                          suffixes=("", "_audit"))
    if "full_crop_in_bounds_audit" in matrix:
        if not matrix.full_crop_in_bounds.eq(matrix.full_crop_in_bounds_audit).all():
            raise RuntimeError("matrix and crop audit disagree on crop support")
        matrix = matrix.drop(columns="full_crop_in_bounds_audit")
    if not matrix.full_crop_in_bounds.eq("YES").any():
        raise RuntimeError("no supported image crops found")
    spec = read_tsv(COLOR_SPEC)
    color_cols = str(spec.iloc[0]["included_features"]).split(";")
    color = pd.read_csv(COLOR_FEATURES, sep="\t", compression="gzip",
                        usecols=["sample_id", "region_x", "region_y", *color_cols])
    color_keys = ["sample_id", "region_x", "region_y"]
    if color.duplicated(color_keys).any():
        raise RuntimeError("D5 color features contain duplicate sample/grid coordinates")
    matrix = matrix.merge(color, on=color_keys, how="left", validate="one_to_one")
    if matrix[color_cols].isna().any().any():
        raise RuntimeError("D5 frozen color features are missing for analysis regions")
    cohort = read_tsv(COHORT)
    cohort = cohort[cohort.primary_inclusion.eq("YES")]
    if matrix.sample_id.isin(cohort.sample_id).sum() != len(matrix):
        raise RuntimeError("matrix contains samples outside frozen admitted cohort")
    folds = read_tsv(FOLDS)
    if folds.donor_id.duplicated().any() or set(matrix.donor_id) - set(folds.donor_id):
        raise RuntimeError("donor folds fail frozen cohort correspondence")
    matrix["fold"] = matrix.donor_id.map(folds.set_index("donor_id").fold)
    matrix["disease_binary"] = matrix.disease.astype(str).eq("pulmonary_fibrosis").astype(float)
    raw, gray, mac = load_embeddings()
    return matrix, color_cols, raw, gray, mac


def build_model_features(frame: pd.DataFrame, color_cols: list[str], raw: dict[str, np.ndarray],
                         gray: dict[str, np.ndarray], mac: dict[str, dict[str, np.ndarray]],
                         outer_fold: str) -> dict[str, np.ndarray]:
    c_cols = [c for c in frame.columns if c.startswith("p_")]
    n_cols = [c for c in frame.columns if c.startswith("N_p_")]
    coords = ["x_center_um", "y_center_um", "normalized_x", "normalized_y", "width_um", "height_um"]
    z_cols = ["normalized_x", "normalized_y"]
    if set(c_cols) != set(f"p_{x}" for x in [
        "epithelial", "endothelial", "fibroblast", "smooth_muscle_pericyte", "macrophage",
        "monocyte", "neutrophil", "CD4_T", "CD8_T", "B", "plasma", "mast",
        "other_immune", "other_stromal", "other"]):
        raise RuntimeError("coarse composition feature set differs from frozen ontology")
    if set(n_cols) != {f"N_{x}" for x in c_cols}:
        raise RuntimeError("neighborhood composition features differ from composition ontology")
    if not set(coords + z_cols + color_cols).issubset(frame.columns):
        raise RuntimeError("frozen covariate, coordinate, or color feature columns missing")
    if frame.analysis_stratum.iloc[0] == "DISEASE_ADJUSTED":
        z_cols = z_cols + ["disease_binary"]
    ids = frame.region_id.astype(str).tolist()
    def embed(source: dict[str, np.ndarray], label: str) -> np.ndarray:
        missing = [rid for rid in ids if rid not in source or not np.isfinite(source[rid]).all()]
        if missing:
            raise RuntimeError(f"{label} embedding missing/invalid for {len(missing)} target-eligible regions")
        return np.stack([source[rid] for rid in ids]).astype(np.float64, copy=False)
    x_raw = embed(raw, "raw")
    x_gray = embed(gray, "grayscale")
    x_mac = embed(mac.get(outer_fold, {}), f"Macenko {outer_fold}")
    z = frame[z_cols].to_numpy(dtype=np.float64)
    c = frame[c_cols].to_numpy(dtype=np.float64)
    n = frame[n_cols].to_numpy(dtype=np.float64)
    color_x = frame[color_cols].to_numpy(dtype=np.float64)
    coord_x = frame[coords].to_numpy(dtype=np.float64)
    cnz = np.column_stack((c, n, z))
    return {
        "M0": z,
        "M1": np.column_stack((c, z)),
        "M2": cnz,
        "M3": np.column_stack((x_raw, z)),
        "M4": np.column_stack((cnz, x_raw)),
        "color": color_x,
        "coordinate": coord_x,
        "grayscale": np.column_stack((cnz, x_gray)),
        "macenko": np.column_stack((cnz, x_mac)),
        "M3_grayscale": np.column_stack((x_gray, z)),
        "M3_macenko": np.column_stack((x_mac, z)),
        "raw_embedding": x_raw,
    }


def fit_predict_one(x_train: np.ndarray, y_train: np.ndarray, d_train: np.ndarray,
                    x_test: np.ndarray) -> tuple[np.ndarray, float]:
    alpha = choose_alpha(x_train, y_train, d_train)
    model = fit_ridge(x_train, y_train, d_train, alpha)
    return model.predict(x_test), alpha


def crossfit_training_residuals(x: np.ndarray, y: np.ndarray, donors: np.ndarray,
                                frozen_folds: np.ndarray, target: str, stratum: str,
                                outer_fold: str, rows: pd.DataFrame) -> tuple[np.ndarray, list[dict]]:
    donors = np.asarray(donors).astype(str)
    folds = np.asarray(frozen_folds).astype(str)
    pred = np.full(len(y), np.nan, dtype=np.float64)
    records: list[dict] = []
    for inner_fold in sorted(np.unique(folds)):
        va = folds == inner_fold
        if not np.any(va):
            continue
        tr = ~va
        if len(np.unique(donors[tr])) < 3:
            raise RuntimeError(f"training-side M2 crossfit has <3 donors: {target}/{stratum}/{outer_fold}/{inner_fold}")
        yhat, alpha = fit_predict_one(x[tr], y[tr], donors[tr], x[va])
        pred[va] = yhat
        block = rows.loc[va, ["region_id", "sample_id", "donor_id", "fold"]].copy()
        block["target"] = target
        block["analysis_stratum"] = stratum
        block["outer_fold"] = outer_fold
        block["train_side_oof_fold"] = inner_fold
        block["observed"] = y[va]
        block["M2_train_oof_pred"] = yhat
        block["M2_train_oof_residual"] = y[va] - yhat
        block["M2_alpha"] = alpha
        records.extend(block.to_dict("records"))
    if not np.isfinite(pred).all():
        raise RuntimeError("incomplete training-side OOF residuals")
    return pred, records


def primary_analysis() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame,
                                pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix, color_cols, raw, gray, mac = load_inputs()
    folds_table = read_tsv(FOLDS)
    fold_by_donor = dict(zip(folds_table.donor_id.astype(str), folds_table.fold.astype(str)))
    prediction_rows: list[pd.DataFrame] = []
    alpha_rows: list[dict] = []
    training_residual_rows: list[dict] = []
    table1_rows: list[dict] = []
    residual_rows: list[dict] = []

    for target_i, (target, spec) in enumerate(TARGETS.items()):
        target_base = matrix[
            matrix["full_crop_in_bounds"].eq("YES")
            & matrix[spec["n"]].ge(20)
            & matrix[spec["y"]].notna()
        ].copy()
        if target_base.empty:
            raise RuntimeError(f"no regions available for primary target {target}")
        for stratum_i, stratum in enumerate(ALL_STRATA):
            frame = target_base.copy()
            if stratum == "PF_ONLY":
                frame = frame[frame.disease.eq("pulmonary_fibrosis")].copy()
            elif stratum == "CONTROL_ONLY":
                frame = frame[frame.disease.eq("control")].copy()
            frame = frame.sort_values(["donor_id", "sample_id", "region_id"]).reset_index(drop=True)
            if frame.empty:
                continue
            n_donors = frame.donor_id.nunique()
            if n_donors < 4:
                if stratum == "CONTROL_ONLY":
                    continue
                raise RuntimeError(f"insufficient donor count for {target}/{stratum}: {n_donors}")
            frame["analysis_stratum"] = stratum
            y = frame[spec["y"]].to_numpy(dtype=np.float64)
            donors = frame.donor_id.astype(str).to_numpy()
            folds = frame.fold.astype(str).to_numpy()
            if set(np.unique(folds)) != set(fold_by_donor[d] for d in np.unique(donors)):
                raise RuntimeError(f"incomplete frozen-fold representation for {target}/{stratum}")
            min_cells = int(frame[spec["n"]].min())
            table1_rows.append({"target": target, "analysis_stratum": stratum,
                                "n_donors": int(n_donors), "n_sections": int(frame.sample_id.nunique()),
                                "n_regions": int(len(frame)), "median_lineage_cells": float(frame[spec["n"]].median()),
                                "min_lineage_cells": min_cells,
                                "primary_minimum_lineage_cells": 20})
            oof = pd.DataFrame({"region_id": frame.region_id, "sample_id": frame.sample_id,
                                "donor_id": frame.donor_id, "fold": frame.fold,
                                "target": target, "analysis_stratum": stratum, "observed": y})
            for name in ["M0_pred", "M1_pred", "M2_pred", "M3_pred", "M4_pred",
                         "color_pred", "coordinate_pred", "grayscale_pred", "macenko_pred", "M3_grayscale_pred", "M3_macenko_pred",
                         "heldout_base_residual", "morphology_predicted_residual"]:
                oof[name] = np.nan
            for outer_fold in sorted(np.unique(folds)):
                print(f"OUTER_START {target}/{stratum}/{outer_fold}", flush=True)
                te = folds == outer_fold
                tr = ~te
                if len(np.unique(donors[tr])) < 3:
                    raise RuntimeError(f"outer training fold has <3 donors: {target}/{stratum}/{outer_fold}")
                outer_rows = frame.loc[tr].reset_index(drop=True)
                test_rows = frame.loc[te].reset_index(drop=True)
                features = build_model_features(frame, color_cols, raw, gray, mac, outer_fold)
                for model_name in ["M0", "M1", "M2", "M3", "M4", "color", "coordinate", "grayscale", "macenko", "M3_grayscale", "M3_macenko"]:
                    print(f"FIT {target}/{stratum}/{outer_fold}/{model_name}", flush=True)
                    x = features[model_name]
                    pred, alpha = fit_predict_one(x[tr], y[tr], donors[tr], x[te])
                    column = {"color": "color_pred", "coordinate": "coordinate_pred",
                              "grayscale": "grayscale_pred", "macenko": "macenko_pred",
                        "M3_grayscale": "M3_grayscale_pred", "M3_macenko": "M3_macenko_pred"}.get(model_name,
                                                                                              f"{model_name}_pred")
                    oof.loc[te, column] = pred
                    alpha_rows.append({"target": target, "analysis_stratum": stratum,
                                       "outer_fold": outer_fold, "model": model_name,
                                       "alpha": alpha, "n_train_regions": int(np.count_nonzero(tr)),
                                       "n_train_donors": int(np.unique(donors[tr]).size),
                                       "alpha_selection": "donor-grouped inner CV; donor-equal MSE; standardized weighted Ridge"})

                # Strict Phase 1G: outer-test M2 residuals and training-only cross-fit residuals.
                x_m2 = features["M2"]
                x_raw = features["raw_embedding"]
                base_test = oof.loc[te, "M2_pred"].to_numpy(dtype=np.float64)
                test_resid = y[te] - base_test
                train_base_oof, records = crossfit_training_residuals(
                    x_m2[tr], y[tr], donors[tr], folds[tr], target, stratum, outer_fold, outer_rows)
                training_residual_rows.extend(records)
                training_residual = y[tr] - train_base_oof
                residual_alpha = choose_alpha(x_raw[tr], training_residual, donors[tr])
                residual_model = fit_ridge(x_raw[tr], training_residual, donors[tr], residual_alpha)
                residual_pred = residual_model.predict(x_raw[te])
                oof.loc[te, "heldout_base_residual"] = test_resid
                oof.loc[te, "morphology_predicted_residual"] = residual_pred
                resid_alpha_row = {"target": target, "analysis_stratum": stratum,
                                   "outer_fold": outer_fold, "model": "PHASE1G_residual_RAW_X",
                                   "alpha": residual_alpha, "n_train_regions": int(np.count_nonzero(tr)),
                                   "n_train_donors": int(np.unique(donors[tr]).size),
                                   "alpha_selection": "donor-grouped inner CV on training-side cross-fitted M2 residuals"}
                alpha_rows.append(resid_alpha_row)
                residual_rows.append({"target": target, "analysis_stratum": stratum,
                                      "outer_fold": outer_fold,
                                      "n_train_donors": int(np.unique(donors[tr]).size),
                                      "n_test_donors": int(np.unique(donors[te]).size),
                                      "n_test_regions": int(np.count_nonzero(te)),
                                      "residual_alpha": residual_alpha,
                                      "residual_R2_fold": r2_score(test_resid, residual_pred)})
                checkpoint_dir = OUT / "model_fold_checkpoints"
                checkpoint_dir.mkdir(exist_ok=True)
                oof.loc[te].to_csv(checkpoint_dir / f"{target}__{stratum}__{outer_fold}.tsv",
                                  sep="\t", index=False)
                pd.DataFrame(alpha_rows).to_csv(checkpoint_dir / "alpha_progress.tsv",
                                               sep="\t", index=False)
                pd.DataFrame(training_residual_rows).to_csv(
                    checkpoint_dir / "training_residual_progress.tsv", sep="\t", index=False)
                print(f"OUTER_COMPLETE {target}/{stratum}/{outer_fold}", flush=True)
            if oof[["M0_pred", "M1_pred", "M2_pred", "M3_pred", "M4_pred", "color_pred",
                    "coordinate_pred", "grayscale_pred", "macenko_pred", "M3_grayscale_pred", "M3_macenko_pred",
                    "heldout_base_residual", "morphology_predicted_residual"]].isna().any().any():
                raise RuntimeError(f"incomplete OOF predictions for {target}/{stratum}")
            prediction_rows.append(oof)
            print(f"MODELLED {target}/{stratum}: donors={n_donors} regions={len(frame)}", flush=True)

    oof_all = pd.concat(prediction_rows, ignore_index=True)
    train_resid = pd.DataFrame(training_residual_rows)
    alpha_frame = pd.DataFrame(alpha_rows)
    table1 = pd.DataFrame(table1_rows)
    resid_fold = pd.DataFrame(residual_rows)
    summary_rows: list[dict] = []
    bootstrap_rows: list[dict] = []
    residual_summary_rows: list[dict] = []
    donor_rows: list[dict] = []
    metric_to_column = {"M0": "M0_pred", "M1": "M1_pred", "M2": "M2_pred", "M3": "M3_pred",
                        "M4": "M4_pred", "color": "color_pred", "coordinate": "coordinate_pred",
                        "grayscale": "grayscale_pred", "macenko": "macenko_pred",
                        "M3_grayscale": "M3_grayscale_pred", "M3_macenko": "M3_macenko_pred"}
    for (target, stratum), group in oof_all.groupby(["target", "analysis_stratum"], sort=False):
        yv = group.observed.to_numpy(float)
        dv = group.donor_id.astype(str).to_numpy()
        preds = {name: group[col].to_numpy(float) for name, col in metric_to_column.items()}
        m = {name: regression_metrics(yv, p) for name, p in preds.items()}
        boot = donor_cluster_bootstrap(yv, preds,
                                       dv, SEED + 100 * list(TARGETS).index(target) + ALL_STRATA.index(stratum))
        delta_comp = m["M1"]["R2"] - m["M0"]["R2"]
        delta_neigh = m["M2"]["R2"] - m["M1"]["R2"]
        delta_morph = m["M4"]["R2"] - m["M2"]["R2"]
        delta_ci = boot["deltaR2_morphology"]
        summary_rows.append({"target": target, "analysis_stratum": stratum,
                             "n_donors": int(group.donor_id.nunique()), "n_sections": int(group.sample_id.nunique()),
                             "n_regions": int(len(group)),
                             **{f"{name}_R2": m[name]["R2"] for name in metric_to_column},
                             **{f"{name}_{metric}": m[name][metric]
                                for name in metric_to_column
                                for metric in ["Pearson_r", "Spearman_rho", "MAE", "RMSE"]},
                             "deltaR2_composition": delta_comp,
                             "deltaR2_neighborhood": delta_neigh,
                             "deltaR2_morphology": delta_morph,
                             "deltaR2_morphology_CI_low": delta_ci["CI_low"],
                             "deltaR2_morphology_CI_high": delta_ci["CI_high"],
                             "CI_low": delta_ci["CI_low"], "CI_high": delta_ci["CI_high"],
                             "bootstrap_two_sided_sign_p": delta_ci["bootstrap_two_sided_sign_p"]})
        for model, metric in m.items():
            b = boot.get(model, {})
            bootstrap_rows.append({"target": target, "analysis_stratum": stratum, "metric": "R2",
                                   "model": model, "estimate": metric["R2"],
                                   "CI_low": b.get("CI_low", np.nan), "CI_high": b.get("CI_high", np.nan),
                                   "bootstrap_replicates": N_BOOT, "bootstrap_unit": "donor"})
        for metric_name, estimate in [("deltaR2_composition", delta_comp), ("deltaR2_neighborhood", delta_neigh),
                                      ("deltaR2_morphology", delta_morph)]:
            ci = boot[metric_name]
            bootstrap_rows.append({"target": target, "analysis_stratum": stratum, "metric": metric_name,
                                   "model": "paired_increment", "estimate": estimate,
                                   "CI_low": ci.get("CI_low", np.nan), "CI_high": ci.get("CI_high", np.nan),
                                   "bootstrap_two_sided_sign_p": ci.get("bootstrap_two_sided_sign_p", np.nan),
                                   "bootstrap_replicates": N_BOOT, "bootstrap_unit": "donor"})
        yr = group.heldout_base_residual.to_numpy(float)
        pr = group.morphology_predicted_residual.to_numpy(float)
        rm = regression_metrics(yr, pr)
        resid_boot = donor_cluster_bootstrap(yr, {"residual": pr}, dv,
                                             SEED + 3000 + 100 * list(TARGETS).index(target) + ALL_STRATA.index(stratum))
        residual_summary_rows.append({"target": target, "analysis_stratum": stratum,
                                      "n_donors": int(group.donor_id.nunique()), "n_regions": int(len(group)),
                                      "residual_R2": rm["R2"], "Pearson": rm["Pearson_r"],
                                      "Spearman": rm["Spearman_rho"], "MAE": rm["MAE"], "RMSE": rm["RMSE"],
                                      "CI_low": resid_boot["residual"]["CI_low"],
                                      "CI_high": resid_boot["residual"]["CI_high"],
                                      "bootstrap_replicates": N_BOOT,
                                      "p": resid_boot["residual"]["bootstrap_two_sided_sign_p"]})
        for donor, dgroup in group.groupby("donor_id", sort=True):
            yy = dgroup.observed.to_numpy(float)
            p2 = dgroup.M2_pred.to_numpy(float)
            p4 = dgroup.M4_pred.to_numpy(float)
            rb = dgroup.heldout_base_residual.to_numpy(float)
            rp = dgroup.morphology_predicted_residual.to_numpy(float)
            sst = float(np.sum((yy - yy.mean()) ** 2))
            donor_rows.append({"target": target, "analysis_stratum": stratum, "donor_id": donor,
                               "n_regions": int(len(dgroup)), "n_sections": int(dgroup.sample_id.nunique()),
                               "observed_variance": float(np.var(yy, ddof=1)) if len(yy) > 1 else np.nan,
                               "M2_R2_within_donor": 1 - np.sum((yy-p2)**2)/sst if sst > 0 else np.nan,
                               "M4_R2_within_donor": 1 - np.sum((yy-p4)**2)/sst if sst > 0 else np.nan,
                               "M2_MAE": float(np.mean(np.abs(yy-p2))), "M4_MAE": float(np.mean(np.abs(yy-p4))),
                               "delta_MAE_M2_minus_M4": float(np.mean(np.abs(yy-p2))-np.mean(np.abs(yy-p4))),
                               "residual_R2_within_donor": 1 - np.sum((rb-rp)**2)/np.sum((rb-rb.mean())**2)
                               if np.sum((rb-rb.mean())**2) > 0 else np.nan,
                               "residual_Pearson_within_donor": float(np.corrcoef(rb, rp)[0, 1])
                               if len(rb) > 1 and np.std(rb) > 0 and np.std(rp) > 0 else np.nan,
                               "residual_observed_variance": float(np.var(rb, ddof=1)) if len(rb) > 1 else np.nan,
                               "M2_SSE": float(np.sum((yy-p2)**2)), "M4_SSE": float(np.sum((yy-p4)**2)),
                               "positive_squared_error_gain": float(max(0, np.sum((yy-p2)**2)-np.sum((yy-p4)**2)))})
    summary = pd.DataFrame(summary_rows)
    primary_mask = summary.analysis_stratum.eq("DISEASE_ADJUSTED")
    q_values = bh_adjust(summary.loc[primary_mask, "bootstrap_two_sided_sign_p"].tolist())
    summary["FDR"] = np.nan
    summary.loc[primary_mask, "FDR"] = q_values
    boot_frame = pd.DataFrame(bootstrap_rows)
    residual_summary = pd.DataFrame(residual_summary_rows)
    residual_summary["FDR"] = np.nan
    residual_mask = residual_summary.analysis_stratum.eq("DISEASE_ADJUSTED")
    residual_summary.loc[residual_mask, "FDR"] = bh_adjust(
        residual_summary.loc[residual_mask, "p"].tolist())
    donors_out = pd.DataFrame(donor_rows)
    donor_resid_out = train_resid
    return (oof_all, alpha_frame, table1, summary, boot_frame,
            residual_summary, donors_out, donor_resid_out, resid_fold)


def main() -> None:
    from threadpoolctl import threadpool_limits
    threadpool_limits(limits=4)
    input_gate()
    needed = [MATRIX, FOLDS, COHORT, COLOR_FEATURES, COLOR_SPEC]
    missing = [str(p) for p in needed if not p.is_file()]
    if missing:
        raise FileNotFoundError("missing frozen Phase 1 inputs: " + "; ".join(missing))
    (oof, alphas, table1, metrics, bootstrap, residual, donor, train_resid,
     residual_folds) = primary_analysis()
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    oof.to_csv(OUT / "oof_predictions.tsv", sep="\t", index=False, na_rep="NA")
    alphas.to_csv(OUT / "ridge_alpha_selection.tsv", sep="\t", index=False)
    table1.to_csv(OUT / "tables/table1_analysis_dataset.tsv", sep="\t", index=False)
    metrics.to_csv(OUT / "tables/table2_primary_model_performance.tsv", sep="\t", index=False)
    residual.to_csv(OUT / "tables/table3_crossfit_residual_performance.tsv", sep="\t", index=False)
    bootstrap.to_csv(OUT / "donor_bootstrap_metrics.tsv", sep="\t", index=False)
    donor.to_csv(OUT / "donor_heterogeneity.tsv", sep="\t", index=False)
    train_resid.to_csv(OUT / "phase1g_train_oof_residuals.tsv", sep="\t", index=False)
    residual_folds.to_csv(OUT / "phase1g_residual_fold_metrics.tsv", sep="\t", index=False)
    # Table 4 is directly reproducible from the complete OOF prediction table.
    metrics["deltaR2_macenko"] = metrics["macenko_R2"] - metrics["M2_R2"]
    metrics["deltaR2_grayscale"] = metrics["grayscale_R2"] - metrics["M2_R2"]
    table4 = metrics[["target", "analysis_stratum", "M4_R2", "macenko_R2", "grayscale_R2",
                      "deltaR2_macenko", "deltaR2_grayscale", "color_R2", "coordinate_R2"]].rename(columns={
                          "M4_R2": "raw", "macenko_R2": "macenko", "grayscale_R2": "grayscale",
                          "color_R2": "color_only", "coordinate_R2": "coordinates_only"})
    table4.to_csv(OUT / "tables/table4_shortcut_stress_test.tsv", sep="\t", index=False)
    oof.rename(columns={"analysis_stratum": "stratum"}).to_csv(
        OUT / "oof_predictions_primary.tsv.gz", sep="\t", index=False, compression="gzip")
    metrics.rename(columns={"analysis_stratum": "stratum", "deltaR2_morphology": "deltaR2_morph",
        "deltaR2_morphology_CI_low": "deltaR2_morph_CI_low",
        "deltaR2_morphology_CI_high": "deltaR2_morph_CI_high",
        "bootstrap_two_sided_sign_p": "deltaR2_morph_p", "FDR": "deltaR2_morph_FDR"}).to_csv(
        OUT / "primary_model_performance.tsv", sep="\t", index=False)
    residual.rename(columns={"analysis_stratum": "stratum", "Pearson": "residual_Pearson",
        "Spearman": "residual_Spearman"}).to_csv(
        OUT / "residual_prediction_performance.tsv", sep="\t", index=False)
    donor.rename(columns={"analysis_stratum": "stratum"}).to_csv(
        OUT / "donor_heterogeneity_primary.tsv", sep="\t", index=False)
    shortcut = metrics[["target", "analysis_stratum", "deltaR2_morphology",
        "deltaR2_macenko", "deltaR2_grayscale", "color_R2", "coordinate_R2",
        "M3_R2", "M3_macenko_R2", "M3_grayscale_R2"]].rename(columns={
        "analysis_stratum": "stratum", "deltaR2_morphology": "raw_deltaR2",
        "deltaR2_macenko": "macenko_deltaR2", "deltaR2_grayscale": "grayscale_deltaR2",
        "color_R2": "color_only_R2", "coordinate_R2": "coordinate_only_R2"})
    shortcut["shortcut_interpretation"] = "MIXED_OR_NO_POSITIVE_INCREMENT"
    color_dom = ((shortcut.raw_deltaR2 > 0) & (shortcut.macenko_deltaR2 <= 0)
                 & (shortcut.grayscale_deltaR2 <= 0))
    shortcut.loc[color_dom, "shortcut_interpretation"] = "COLOR_SENSITIVE_INCREMENT"
    structure = ((shortcut.macenko_deltaR2 > 0) & (shortcut.grayscale_deltaR2 > 0)
                 & (shortcut.M3_R2 > shortcut.color_only_R2)
                 & (shortcut.M3_R2 > shortcut.coordinate_only_R2))
    shortcut.loc[structure, "shortcut_interpretation"] = "STRUCTURE_COMPATIBLE_SIGNAL"
    shortcut.to_csv(OUT / "shortcut_robustness.tsv", sep="\t", index=False)
    input_gate()
    print(f"PHASE1E_G_MODELS_PASS: predictions={len(oof)}; strata={oof.analysis_stratum.nunique()}; "
          f"targets={oof.target.nunique()}; donor_bootstrap={N_BOOT}", flush=True)


if __name__ == "__main__":
    main()
