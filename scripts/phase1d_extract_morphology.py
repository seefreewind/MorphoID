#!/usr/bin/env python3
"""Extract frozen Phikon-v2 features from fixed 150-um registered H&E crops."""
from __future__ import annotations

import csv
import gc
import gzip
import hashlib
import json
import os
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import torch
from PIL import Image
from transformers import AutoImageProcessor, AutoModel

ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "results/phase1"
MATRIX = P1 / "region_analysis_matrix.tsv"
CROP_AUDIT = P1 / "image_region_crop_audit.tsv"
COHORT = ROOT / "configs/frozen_vannan_primary_cohort.tsv"
FOLDS = ROOT / "configs/phase1_donor_folds.tsv"
MODEL_DIR = P1 / "model_cache/owkin_phikon-v2"
MODEL_REVISION = "2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd"
P1_TARGETS = ["epithelial__injury", "fibroblast__fibroblast_activation",
              "macrophage__inflammatory"]
BATCH = 4
PIXEL_SAMPLE_N = 20_000
PIXEL_SAMPLE_FALLBACK_N = 100_000
PIXEL_SAMPLE_BLOCK_SIZE = 100
BETA_OD = 0.15
ALPHA_PERCENTILE = 1.0
DECOMPRESS_CHUNK_BYTES = 16 * 1024 * 1024
DECOMPRESS_LOG_INTERVAL = 256 * 1024 * 1024
CROP_LOG_INTERVAL = 50


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def pixels_to_od(rgb: np.ndarray) -> np.ndarray:
    return -np.log((rgb.astype(np.float32) + 1.0) / 256.0)


def order_stains(stains: np.ndarray) -> np.ndarray:
    # Hematoxylin has the stronger blue-channel optical density; keep stain
    # columns consistently ordered H then E across sections and folds.
    stains = np.asarray(stains, dtype=np.float64)
    for j in range(stains.shape[1]):
        if stains[:, j].sum() < 0:
            stains[:, j] *= -1
        stains[:, j] /= max(np.linalg.norm(stains[:, j]), 1e-12)
    if stains[2, 0] < stains[2, 1]:
        stains = stains[:, ::-1]
    return stains


def estimate_stain_parameters(rgb: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """Macenko source stain vectors and 99th-percentile concentrations."""
    od = pixels_to_od(rgb.reshape(-1, 3))
    mask = (od.max(axis=1) > BETA_OD) & (od.max(axis=1) < 4.0)
    od = od[mask]
    if len(od) < 500:
        raise ValueError(f"too few stained pixels for Macenko fit: n={len(od)}")
    cov = np.cov(od, rowvar=False)
    eigval, eigvec = np.linalg.eigh(cov)
    plane = eigvec[:, np.argsort(eigval)[-2:]]
    projected = od @ plane
    angle = np.arctan2(projected[:, 1], projected[:, 0])
    lo, hi = np.percentile(angle, [ALPHA_PERCENTILE, 100.0 - ALPHA_PERCENTILE])
    raw = np.column_stack((plane @ [np.cos(lo), np.sin(lo)],
                           plane @ [np.cos(hi), np.sin(hi)]))
    stains = order_stains(raw)
    concentration = np.maximum(od @ np.linalg.pinv(stains.T), 0.0)
    p99 = np.percentile(concentration, 99.0, axis=0)
    return stains, np.maximum(p99, 1e-6), int(len(od))


def sample_slide_pixels(image: np.ndarray, n: int = PIXEL_SAMPLE_N) -> np.ndarray:
    h, w = image.shape[:2]
    total = h * w
    n_pixels = min(total, n)
    n_blocks = int(np.ceil(n_pixels / PIXEL_SAMPLE_BLOCK_SIZE))
    last_start = max(total - PIXEL_SAMPLE_BLOCK_SIZE, 0)
    starts = np.linspace(0, last_start, n_blocks, dtype=np.int64)
    flat = image.reshape(-1, 3)
    sampled = np.empty((n_pixels, 3), dtype=np.uint8)
    cursor = 0
    for start in starts:
        take = min(PIXEL_SAMPLE_BLOCK_SIZE, n_pixels - cursor)
        sampled[cursor:cursor + take] = flat[start:start + take]
        cursor += take
        if cursor >= n_pixels:
            break
    if cursor != n_pixels:
        raise RuntimeError(f"incomplete deterministic pixel sample: {cursor}/{n_pixels}")
    return sampled


def make_target_reference(source: dict[str, dict], donor_by_sample: dict[str, str],
                          training_donors: set[str]) -> tuple[np.ndarray, np.ndarray]:
    by_donor: dict[str, list[dict]] = defaultdict(list)
    for sample, values in source.items():
        donor = donor_by_sample[sample]
        if donor in training_donors:
            by_donor[donor].append(values)
    if set(by_donor) != training_donors:
        missing = sorted(training_donors - set(by_donor))
        raise RuntimeError(f"Macenko training reference missing donors: {missing}")
    donor_stains, donor_p99 = [], []
    for donor in sorted(by_donor):
        donor_stains.append(np.median(np.stack([x["stains"] for x in by_donor[donor]]), axis=0))
        donor_p99.append(np.median(np.stack([x["p99"] for x in by_donor[donor]]), axis=0))
    target_stains = order_stains(np.median(np.stack(donor_stains), axis=0))
    target_p99 = np.maximum(np.median(np.stack(donor_p99), axis=0), 1e-6)
    return target_stains, target_p99


def macenko_transform(rgb: np.ndarray, source_stains: np.ndarray, source_p99: np.ndarray,
                      target_stains: np.ndarray, target_p99: np.ndarray) -> np.ndarray:
    od = pixels_to_od(rgb.reshape(-1, 3))
    concentrations = np.maximum(od @ np.linalg.pinv(source_stains.T), 0.0)
    normalized = concentrations * (target_p99 / np.maximum(source_p99, 1e-6))
    reconstructed = normalized @ target_stains.T
    out = np.clip(256.0 * np.exp(-reconstructed) - 1.0, 0, 255)
    return out.reshape(rgb.shape).astype(np.uint8)


def grayscale_rgb(rgb: np.ndarray) -> np.ndarray:
    gray = np.rint(rgb[..., 0] * 0.2126 + rgb[..., 1] * 0.7152 + rgb[..., 2] * 0.0722)
    gray = np.clip(gray, 0, 255).astype(np.uint8)
    return np.repeat(gray[..., None], 3, axis=2)


def tissue_fraction_proxy(rgb: np.ndarray) -> float:
    # Same 8-pixel sampling interval as the Phase 0D color audit. This is an
    # OD-based stained-pixel proxy, not a validated histologic tissue mask.
    sub = rgb[::8, ::8, :]
    od = pixels_to_od(sub.reshape(-1, 3))
    nonblack = sub.reshape(-1, 3).mean(axis=1) > 0.05 * 255
    return float(np.mean((od.sum(axis=1) > 0.15) & nonblack))


def batch_embed(images: list[np.ndarray], processor, model, device: torch.device,
                stage: str) -> np.ndarray:
    print(f"  FEATURE_PROCESSOR_START {stage} n={len(images)}", flush=True)
    pil = [Image.fromarray(im) for im in images]
    batch = processor(images=pil, return_tensors="pt")["pixel_values"].to(device)
    print(f"  FEATURE_MODEL_FORWARD_START {stage} device={device}", flush=True)
    with torch.inference_mode():
        hidden = model(pixel_values=batch).last_hidden_state[:, 0, :]
    result = hidden.float().cpu().numpy().astype(np.float32, copy=False)
    print(f"  FEATURE_MODEL_FORWARD_COMPLETE {stage} shape={result.shape}", flush=True)
    return result


def open_registered_tiff_gz(path: Path) -> tuple[np.ndarray, Path]:
    temp_dir = P1 / "tiff_work"
    temp_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix="phase1_", suffix=".tif", dir=temp_dir,
                                     delete=False) as temp:
        temp_path = Path(temp.name)
    try:
        with gzip.open(path, "rb") as src, temp_path.open("wb") as dst:
            copied = 0
            next_log = DECOMPRESS_LOG_INTERVAL
            while chunk := src.read(DECOMPRESS_CHUNK_BYTES):
                dst.write(chunk)
                copied += len(chunk)
                if copied >= next_log:
                    print(f"  DECOMPRESS {path.name}: {copied / (1024 ** 2):.0f} MiB", flush=True)
                    next_log = ((copied // DECOMPRESS_LOG_INTERVAL) + 1) * DECOMPRESS_LOG_INTERVAL
        print(f"  DECOMPRESS_COMPLETE {path.name}: {copied / (1024 ** 3):.2f} GiB", flush=True)
        mapped = tifffile.memmap(temp_path, mode="r")
        if mapped.ndim != 3 or mapped.shape[2] != 3 or mapped.dtype != np.uint8:
            raise ValueError(f"unexpected registered TIFF array layout: shape={mapped.shape}, dtype={mapped.dtype}")
        print(f"  TIFF_RAM_COPY_START shape={mapped.shape} bytes={mapped.nbytes}", flush=True)
        image = np.array(mapped, copy=True, order="C")
        del mapped
        print(f"  TIFF_RAM_COPY_COMPLETE shape={image.shape} bytes={image.nbytes}", flush=True)
        return image, temp_path
    except BaseException:
        temp_path.unlink(missing_ok=True)
        raise


def crop_from_slide(image: np.ndarray, row: dict) -> np.ndarray:
    crop = np.asarray(image[int(row["y0_px"]):int(row["y1_px_exclusive"]),
                            int(row["x0_px"]):int(row["x1_px_exclusive"]), :])
    expected = (int(row["crop_height_px"]), int(row["crop_width_px"]), 3)
    if crop.shape != expected:
        raise ValueError(f"crop shape mismatch for {row['region_id']}: {crop.shape} != {expected}")
    return crop


def extract_batches(items: list[tuple[str, np.ndarray]], processor, model,
                    device: torch.device, stage: str) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    n_batches = (len(items) + BATCH - 1) // BATCH
    for start in range(0, len(items), BATCH):
        batch_items = items[start:start + BATCH]
        batch_id = start // BATCH + 1
        label = f"{stage} batch={batch_id}/{n_batches}"
        print(f"EMBED_BATCH_START {label} items={len(batch_items)}", flush=True)
        features = batch_embed([x[1] for x in batch_items], processor, model, device, label)
        for (region_id, _), feature in zip(batch_items, features):
            out[region_id] = feature
        print(f"EMBED_BATCH_COMPLETE {label}", flush=True)
    return out


def main() -> None:
    if not (MODEL_DIR / "model.safetensors").is_file():
        raise FileNotFoundError("run scripts/phase1d_fetch_encoder.py first")
    matrix = read_tsv(MATRIX)
    crop_audit = read_tsv(CROP_AUDIT)
    cohort = read_tsv(ROOT / "configs/frozen_vannan_primary_cohort.tsv")
    folds_df = read_tsv(ROOT / "configs/phase1_donor_folds.tsv")
    donor_by_sample = dict(zip(cohort.loc[cohort.primary_inclusion == "YES", "sample_id"],
                               cohort.loc[cohort.primary_inclusion == "YES", "donor_id"]))
    fold_by_donor = dict(zip(folds_df["donor_id"], folds_df["fold"]))
    matrix = matrix.merge(pd.DataFrame(crop_audit)[["region_id", "full_crop_in_bounds"]],
                          on="region_id", how="left", validate="one_to_one")
    if matrix["full_crop_in_bounds"].isna().any():
        raise RuntimeError("crop audit does not match frozen region matrix")
    valid = matrix[matrix.full_crop_in_bounds.eq("YES")].copy()
    primary_union = np.zeros(len(matrix), dtype=bool)
    for target in P1_TARGETS:
        primary_union |= matrix[f"n_target_cells__{target}"].to_numpy(int) >= 20
    primary_by_id = set(matrix.loc[primary_union & matrix.full_crop_in_bounds.eq("YES"), "region_id"])
    region_ids = matrix["region_id"].astype(str).tolist()
    all_index = {rid: i for i, rid in enumerate(region_ids)}

    source_path_by_sample = {}
    d5 = read_tsv(ROOT / "results/phase0/phase0d/d5_color_feature_extraction_manifest.tsv")
    d5_records = d5.to_dict("records")
    image_sha_by_sample = {str(r["sample_id"]): str(r["registered_HE_sha256"])
                           for r in d5_records}
    for r in d5_records:
        source_path_by_sample[r["sample_id"]] = ROOT / "data/raw/phase0c_v" / r["image_file"]
    if set(source_path_by_sample) != set(donor_by_sample):
        raise RuntimeError("registered-image manifest differs from frozen admitted cohort")

    # Reuse the deterministic source-stain table only when every row matches
    # the current admitted image SHA256 and the frozen sampling/fit rule.
    source_out = P1 / "macenko_section_source_parameters.tsv"
    sample_rule = (
        "deterministic equally spaced whole-slide block starts; 100 contiguous RGB pixels per block; "
        "primary 200 blocks/20000 pixels; fallback 1000 blocks/100000 pixels only if fewer than 500 pass OD filter"
    )
    source_params: dict[str, dict] = {}
    stain_rows = []
    stain_cache_reused = False
    if source_out.is_file():
        try:
            cached = read_tsv(source_out)
            cache_ok = (
                len(cached) == len(source_path_by_sample)
                and cached["sample_id"].astype(str).is_unique
                and set(cached["sample_id"].astype(str)) == set(source_path_by_sample)
                and cached["sample_rule"].astype(str).eq(sample_rule).all()
            )
            if cache_ok:
                for row in cached.to_dict("records"):
                    sample = str(row["sample_id"])
                    stains = np.asarray(json.loads(row["source_stain_matrix_H_then_E"]), dtype=np.float64)
                    p99 = np.asarray(json.loads(row["source_concentration_p99_H_E"]), dtype=np.float64)
                    cache_ok = (
                        str(row["donor_id"]) == str(donor_by_sample[sample])
                        and str(row["source_image_sha256"]) == str(image_sha_by_sample[sample])
                        and stains.shape == (3, 2) and p99.shape == (2,)
                        and np.isfinite(stains).all() and np.isfinite(p99).all()
                    )
                    if not cache_ok:
                        break
                    source_params[sample] = {"stains": stains, "p99": p99}
                if cache_ok and len(source_params) == len(source_path_by_sample):
                    stain_rows = cached.to_dict("records")
                    stain_cache_reused = True
                    print("STAIN_SOURCE_CACHE_REUSED rows=26; all image SHA256, donors, and sampling rules match",
                          flush=True)
        except Exception as exc:
            print(f"STAIN_SOURCE_CACHE_REJECTED {type(exc).__name__}: {exc}", flush=True)
    if not stain_cache_reused:
        source_params.clear()
        stain_rows.clear()
        print("STAIN_SOURCE_CACHE_REJECTED; recomputing all source stain parameters", flush=True)
        for i, sample in enumerate(sorted(source_path_by_sample), 1):
            print(f"STAIN_SOURCE_START [{i}/26] {sample}", flush=True)
            image, temp_path = open_registered_tiff_gz(source_path_by_sample[sample])
            try:
                print(f"  STAIN_PIXEL_SAMPLE_START n={PIXEL_SAMPLE_N}", flush=True)
                rgb = sample_slide_pixels(image)
                fallback_used = False
                try:
                    stains, p99, n_pixels = estimate_stain_parameters(rgb)
                except ValueError as exc:
                    if "too few stained pixels" not in str(exc):
                        raise
                    fallback_used = True
                    print(f"  STAIN_PIXEL_SAMPLE_FALLBACK n={PIXEL_SAMPLE_FALLBACK_N}", flush=True)
                    rgb = sample_slide_pixels(image, PIXEL_SAMPLE_FALLBACK_N)
                    stains, p99, n_pixels = estimate_stain_parameters(rgb)
                source_params[sample] = {"stains": stains, "p99": p99}
                stain_rows.append({"sample_id": sample, "donor_id": donor_by_sample[sample],
                                   "sampled_pixels": len(rgb), "stained_pixels_after_OD_filter": n_pixels,
                                   "sample_size_fallback_used": "YES" if fallback_used else "NO",
                                   "sample_rule": sample_rule,
                                   "source_stain_matrix_H_then_E": json.dumps(stains.tolist()),
                                   "source_concentration_p99_H_E": json.dumps(p99.tolist()),
                                   "source_image_sha256": image_sha_by_sample[sample],
                                   "fit_basis": "within-section source stain basis for transform only; no held-out data used in training-donor target reference"})
            finally:
                del image
                temp_path.unlink(missing_ok=True)
                gc.collect()
            print(f"STAIN_SOURCE [{i}/26] {sample}", flush=True)
        pd.DataFrame(stain_rows).to_csv(source_out, sep="\t", index=False)
    fold_rows = folds_df.to_dict("records")
    target_refs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    target_ref_rows = []
    for fold in sorted(folds_df.fold.unique()):
        held = set(folds_df.loc[folds_df.fold.eq(fold), "donor_id"])
        training = set(folds_df.donor_id) - held
        matrix_ref, p99_ref = make_target_reference(source_params, donor_by_sample, training)
        target_refs[fold] = (matrix_ref, p99_ref)
        target_ref_rows.append({"outer_fold": fold, "training_donors": ";".join(sorted(training)),
                                "heldout_donors": ";".join(sorted(held)),
                                "target_stain_matrix_H_then_E": json.dumps(matrix_ref.tolist()),
                                "target_concentration_p99_H_E": json.dumps(p99_ref.tolist()),
                                "reference_rule": "median within donor across sections, then median across training donors; held-out donors excluded"})
    pd.DataFrame(target_ref_rows).to_csv(P1 / "macenko_fold_target_parameters.tsv", sep="\t", index=False)
    source_json = {s: {"stains": p["stains"].tolist(), "p99": p["p99"].tolist()}
                   for s, p in source_params.items()}

    device_override = os.environ.get("MORPHOID_PHASE1D_DEVICE", "").strip().lower()
    if device_override and device_override not in {"cpu", "mps"}:
        raise ValueError("MORPHOID_PHASE1D_DEVICE must be 'cpu' or 'mps'")
    if device_override == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("MORPHOID_PHASE1D_DEVICE=mps requested but MPS is unavailable")
    device = torch.device(device_override or ("mps" if torch.backends.mps.is_available() else "cpu"))
    backend_choice = f"override={device_override}" if device_override else "override=AUTO"
    print(f"PHASE1D_DEVICE {device}; torch={torch.__version__}; batch={BATCH}; {backend_choice}", flush=True)
    processor = AutoImageProcessor.from_pretrained(MODEL_DIR, local_files_only=True)
    model = AutoModel.from_pretrained(MODEL_DIR, local_files_only=True, use_safetensors=True)
    model.eval().to(device)
    torch.set_num_threads(4)
    cache_dir = P1 / "embedding_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    file_manifest = []
    tissue_by_id: dict[str, float] = {}

    crop_rows_by_sample: dict[str, list[dict]] = defaultdict(list)
    crop_records = crop_audit.to_dict("records")
    for row in crop_records:
        if row["full_crop_in_bounds"] == "YES":
            crop_rows_by_sample[row["sample_id"]].append(row)

    for i, sample in enumerate(sorted(source_path_by_sample), 1):
        output = cache_dir / f"{sample}.npz"
        image_input = source_path_by_sample[sample]
        if output.exists():
            with np.load(output, allow_pickle=False) as z:
                if str(z["model_revision"].item()) != MODEL_REVISION:
                    raise RuntimeError(f"embedding cache revision mismatch: {output}")
                cache_device = (str(z["inference_device"].item())
                                if "inference_device" in z.files else "UNRECORDED")
                if cache_device not in {"UNRECORDED", str(device)}:
                    raise RuntimeError(f"embedding cache backend mismatch: {output}; "
                                       f"cache={cache_device}, requested={device}")
                for rid, fraction in zip(z["region_ids"].astype(str), z["tissue_fraction"]):
                    tissue_by_id[rid] = float(fraction)
            file_manifest.append({"sample_id": sample, "donor_id": donor_by_sample[sample],
                                  "image_sha256": image_sha_by_sample[sample],
                                  "embedding_file": str(output.relative_to(ROOT)),
                                  "embedding_file_sha256": sha256(output), "cache_reused": "YES",
                                  "inference_device": cache_device})
            print(f"EMBEDDING_CACHE_REUSED [{i}/26] {sample}", flush=True)
            continue

        rows = crop_rows_by_sample[sample]
        print(f"EMBEDDING_START [{i}/26] {sample}; crops={len(rows)}", flush=True)
        image, temp_path = open_registered_tiff_gz(image_input)
        try:
            raw_items, gray_items = [], []
            for crop_i, row in enumerate(rows, 1):
                rid = row["region_id"]
                crop = crop_from_slide(image, row)
                tissue_by_id[rid] = tissue_fraction_proxy(crop)
                raw_items.append((rid, crop))
                if rid in primary_by_id:
                    gray_items.append((rid, grayscale_rgb(crop)))
                if crop_i == 1 or crop_i % CROP_LOG_INTERVAL == 0 or crop_i == len(rows):
                    print(f"CROP_COLLECTION_PROGRESS [{i}/26] {sample}: {crop_i}/{len(rows)}",
                          flush=True)
            print(f"CROP_COLLECTION_COMPLETE [{i}/26] {sample}; raw={len(raw_items)} grayscale={len(gray_items)}",
                  flush=True)
            raw = extract_batches(raw_items, processor, model, device, f"{sample} raw")
            gray = extract_batches(gray_items, processor, model, device, f"{sample} grayscale")
            n_raw, n_gray = len(raw_items), len(gray_items)
            del raw_items, gray_items
            gc.collect()
            region_sample_ids = [r["region_id"] for r in rows]
            n = len(region_sample_ids)
            raw_arr = np.full((n, 1024), np.nan, dtype=np.float32)
            gray_arr = np.full((n, 1024), np.nan, dtype=np.float32)
            position = {rid: j for j, rid in enumerate(region_sample_ids)}
            for rid, feat in raw.items():
                raw_arr[position[rid]] = feat
            for rid, feat in gray.items():
                gray_arr[position[rid]] = feat
            mac_arrays = {fold: np.full((n, 1024), np.nan, dtype=np.float32)
                          for fold in target_refs}
            mac_rows = [r for r in rows if r["region_id"] in primary_by_id]
            for fold, (target_stains, target_p99) in target_refs.items():
                source = source_params[sample]
                items = []
                for row in mac_rows:
                    crop = crop_from_slide(image, row)
                    normalized = macenko_transform(crop, source["stains"], source["p99"],
                                                   target_stains, target_p99)
                    items.append((row["region_id"], normalized))
                features = extract_batches(items, processor, model, device, f"{sample} macenko {fold}")
                arr = mac_arrays[fold]
                for rid, feat in features.items():
                    arr[position[rid]] = feat
                print(f"  MACENKO {fold}: regions={len(features)}", flush=True)
                del items, features
                gc.collect()

            save_values = {"region_ids": np.asarray(region_sample_ids, dtype="U64"),
                           "raw": raw_arr, "grayscale": gray_arr,
                           "model_revision": np.asarray(MODEL_REVISION),
                           "inference_device": np.asarray(str(device)),
                           "tissue_fraction": np.asarray([tissue_by_id[rid] for rid in region_sample_ids], dtype=np.float32)}
            for fold, arr in mac_arrays.items():
                save_values[f"macenko_{fold}"] = arr
            np.savez(output, **save_values)
            file_manifest.append({"sample_id": sample, "donor_id": donor_by_sample[sample],
                                  "image_sha256": image_sha_by_sample[sample],
                                  "n_regions_in_bounds": n, "n_primary_union_regions": len(mac_rows),
                                  "embedding_file": str(output.relative_to(ROOT)),
                                  "embedding_file_sha256": sha256(output), "cache_reused": "NO",
                                  "inference_device": str(device)})
            print(f"EMBEDDED [{i}/26] {sample}: raw={n_raw} grayscale={n_gray} mac_fold_regions={len(mac_rows)}",
                  flush=True)
        finally:
            del image
            temp_path.unlink(missing_ok=True)
            gc.collect()

    pd.DataFrame(file_manifest).to_csv(P1 / "morphology_embedding_extraction_manifest.tsv",
                                       sep="\t", index=False)
    # Fill the Phase 1A matrix's image-derived tissue-fraction proxy and crop support.
    matrix["tissue_fraction"] = matrix["region_id"].map(tissue_by_id)
    crop_support = dict(zip(crop_audit["region_id"], crop_audit["full_crop_in_bounds"]))
    matrix["image_crop_supported"] = matrix["region_id"].map(crop_support)
    matrix["tissue_fraction_definition"] = "OD sum > 0.15 among every-8th-pixel samples; pure black excluded; proxy only"
    matrix.to_csv(MATRIX, sep="\t", index=False, na_rep="NA")
    print(f"PHASE1D_EMBEDDINGS_PASS: device={device}; samples={len(file_manifest)}; regions={len(valid)}; batch={BATCH}",
          flush=True)


if __name__ == "__main__":
    main()
