#!/usr/bin/env python3
"""Frozen no-refit H&E bounds and tissue-mask audit for GSE250346.

Reads each author-registered TIFF once, creates only a 3000-pixel thumbnail,
and evaluates the identity map implied by the published common coordinate
space plus per-image TIFF scale. No translation, orientation, scale, or warp
is searched or fitted. R4/R5/R6 remain separately gated.
"""
from __future__ import annotations

import csv
import gzip
import io
import os
import shutil
import sys
import tempfile
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import tifffile

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/phase0/phase0c_v"
MANIFEST = BASE / "phase0c_v_input_manifest.tsv"
CELL_FILE = BASE / "vannan_primary_cell_metadata.tsv.gz"
META_FILE = BASE / "registered_he_tiff_metadata.tsv"
OUT_FILE = BASE / "registration_qc.tsv"
SUMMARY_FILE = BASE / "registration_image_qc_summary.tsv"
CONTACT = ROOT / "figures/phase0c_v"
THUMB_MAX = 3000
R2_THRESHOLD = 0.95
R3_THRESHOLD = 0.90


def emit(message: str) -> None:
    print(message, flush=True)


def make_thumbnail(tmp_path: Path, width: int, height: int) -> np.ndarray:
    # These source TIFFs are single-strip, uncompressed RGB. Memmap the pixel
    # strip so peak RAM is the 3000-pixel thumbnail plus its HSV masks.
    thumb_w, thumb_h = (width, height)
    scale = min(1.0, THUMB_MAX / max(width, height))
    thumb_w = max(1, int(round(width * scale)))
    thumb_h = max(1, int(round(height * scale)))
    with tifffile.TiffFile(tmp_path, is_imagej=False) as tf:
        page = tf.pages[0]
        if page.compression != 1:
            raise ValueError(f"inner TIFF compression is {page.compression}, expected none")
        if page.samplesperpixel != 3 or page.planarconfig != 1:
            raise ValueError(f"expected chunky RGB; samples={page.samplesperpixel}; planar={page.planarconfig}")
        if len(page.dataoffsets) != 1 or len(page.databytecounts) != 1:
            raise ValueError(f"expected one pixel strip; found offsets={len(page.dataoffsets)}")
        offset = int(page.dataoffsets[0])
        rows_per_strip = int(page.rowsperstrip)
        if rows_per_strip < height:
            raise ValueError(f"multiple strip rows are not supported: {rows_per_strip} < {height}")
    pixels = np.memmap(tmp_path, mode="r", dtype=np.uint8, offset=offset,
                       shape=(height, width, 3), order="C")
    thumb = cv2.resize(pixels, (thumb_w, thumb_h), interpolation=cv2.INTER_AREA)
    return thumb


def primary_mask(rgb: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    raw10 = ((hsv[:, :, 1] > 10) & (hsv[:, :, 2] < 245)).astype(np.uint8)
    raw18 = ((hsv[:, :, 1] > 18) & (hsv[:, :, 2] < 245)).astype(np.uint8)
    raw10 = cv2.morphologyEx(raw10, cv2.MORPH_CLOSE, np.ones((3, 3), np.uint8))
    raw18 = cv2.morphologyEx(raw18, cv2.MORPH_CLOSE, np.ones((3, 3), np.uint8))
    nlabels, labels, stats, _ = cv2.connectedComponentsWithStats(raw18, 8)
    if nlabels <= 1:
        return np.zeros_like(raw10, dtype=bool), {"components": 0, "selected_components": 0}
    max_area = int(stats[1:, cv2.CC_STAT_AREA].max())
    selected = [k for k in range(1, nlabels) if stats[k, cv2.CC_STAT_AREA] >= 0.10 * max_area]
    support = np.zeros_like(raw18, dtype=np.uint8)
    kernel = np.ones((31, 31), dtype=np.uint8)  # frozen 15-pixel dilation radius
    for label in selected:
        comp = (labels == label).astype(np.uint8)
        support |= cv2.dilate(comp, kernel)
    tissue = (raw10 > 0) & (support > 0)
    return tissue, {"components": nlabels - 1, "selected_components": len(selected)}


def main() -> None:
    cv2.setNumThreads(1)
    assets = [r for r in csv.DictReader(MANIFEST.open(encoding="utf-8"), delimiter="\t")
              if r["role"] == "author-registered same-section H&E"]
    cells = pd.read_csv(CELL_FILE, sep="\t", usecols=["sample_id", "donor_id", "x_um", "y_um"])
    metadata = pd.read_csv(META_FILE, sep="\t").set_index("sample_id")
    existing = pd.read_csv(OUT_FILE, sep="\t").set_index("sample_id") if OUT_FILE.exists() else None
    out = []
    qc_dir = CONTACT / "registration_overlays"
    qc_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = BASE / ".registration_tiff_tmp"
    temp_dir.mkdir(exist_ok=True)
    try:
        for i, asset in enumerate(assets, start=1):
            p = ROOT / asset["file_path"]
            gsm = Path(p.name).name.split("_")[0]
            row = metadata[metadata.GSM == gsm]
            if len(row) != 1:
                raise ValueError(f"GSM {gsm} does not resolve to one TIFF metadata row")
            info = row.iloc[0]
            sample = str(info.sample_id)
            sub = cells[cells.sample_id == sample]
            if len(sub) == 0:
                raise ValueError(f"No corrected primary cells for {sample}")
            px = float(info.physical_pixel_size_x_um) if pd.notna(info.physical_pixel_size_x_um) else 0.2125
            py = float(info.physical_pixel_size_y_um) if pd.notna(info.physical_pixel_size_y_um) else 0.2125
            scale_known = info.ImageJ_spatial_unit in {"micron", "microns", "um"} and pd.notna(info.physical_pixel_size_x_um)
            xthumb = sub.x_um.to_numpy(float) / px * min(1.0, THUMB_MAX / max(int(info.width_px), int(info.height_px)))
            ythumb = sub.y_um.to_numpy(float) / py * min(1.0, THUMB_MAX / max(int(info.width_px), int(info.height_px)))
            thw = max(1, int(round(int(info.width_px) * min(1.0, THUMB_MAX / max(int(info.width_px), int(info.height_px))))))
            thh = max(1, int(round(int(info.height_px) * min(1.0, THUMB_MAX / max(int(info.width_px), int(info.height_px))))))
            inb = (xthumb >= 0) & (xthumb < thw) & (ythumb >= 0) & (ythumb < thh)
            r2 = float(inb.mean()) if scale_known else np.nan
            if existing is not None and sample in existing.index:
                prior_r2 = existing.loc[sample, "centroids_inside_image"]
                if pd.notna(prior_r2) and scale_known:
                    r2 = float(prior_r2)
            emit(f"[{i}/{len(assets)}] {sample}: decompressing registered H&E ({int(info.width_px)}x{int(info.height_px)})")
            tmp = temp_dir / f"{gsm}.tif"
            try:
                with gzip.open(p, "rb") as src, tmp.open("wb") as dst:
                    shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)
                thumb = make_thumbnail(tmp, int(info.width_px), int(info.height_px))
                tissue, components = primary_mask(thumb)
                ix = np.clip(np.rint(xthumb).astype(int), 0, thw - 1)
                iy = np.clip(np.rint(ythumb).astype(int), 0, thh - 1)
                r3 = float((inb & tissue[iy, ix]).mean())
                # Outcome-blind centroid overlay; no expression is loaded.
                overlay = thumb.copy()
                valid = np.flatnonzero(inb)
                for idx in valid:
                    cv2.circle(overlay, (int(ix[idx]), int(iy[idx])), 1, (255, 40, 40), -1,
                               lineType=cv2.LINE_AA)
                cv2.imwrite(str(qc_dir / f"{sample}_identity_map_centroids.png"),
                            cv2.cvtColor(overlay, cv2.COLOR_RGB2BGR))
                out.append({
                    "sample_id": sample, "donor_id": str(sub.donor_id.iloc[0]), "GSM": gsm,
                    "registered_HE": str(asset["file_path"]), "registered_HE_sha256": asset["sha256"],
                    "n_cells": len(sub), "R1_physical_scale_status": "PASS_SCALE_TAG" if scale_known else "UNRESOLVED_MISSING_PHYSICAL_UNIT",
                    "pixel_size_x_um": px if scale_known else np.nan,
                    "pixel_size_y_um": py if scale_known else np.nan,
                    "mapping_rule": "author-registered image; identity orientation and zero translation; no search or refit",
                    "centroids_inside_image": r2, "R2_threshold": R2_THRESHOLD,
                    "R2_status": ("PASS" if r2 >= R2_THRESHOLD else "FAIL") if scale_known else "UNRESOLVED",
                    "tissue_coverage": r3, "R3_threshold": R3_THRESHOLD,
                    "R3_status": "PASS" if r3 >= R3_THRESHOLD else "FAIL",
                    "n_large_saturation_components": components["selected_components"],
                    "all_saturation_components": components["components"],
                    "dapi_rho": np.nan, "R5_status": "UNAVAILABLE_DAPI_NOT_IN_BOUNDED_INPUTS",
                    "R4_status": "N_A_NO_WITHHELD_LANDMARKS_SUPPLIED",
                    "visual_pass": "NOT_REVIEWED", "n_visual_patches": 0,
                    "R6_status": "UNRESOLVED_NO_PATCH_REVIEW",
                    "registration_verdict": "FAIL" if (scale_known and (r2 < R2_THRESHOLD or r3 < R3_THRESHOLD)) else "UNRESOLVED",
                    "reason": "R4 landmarks absent and R5 DAPI concordance unavailable; frozen policy requires R6 as well. Direct identity-map R2/R3 are diagnostics only.",
                })
            finally:
                if tmp.exists():
                    tmp.unlink()
    finally:
        try:
            temp_dir.rmdir()
        except OSError:
            pass
    q = pd.DataFrame(out).sort_values("sample_id")
    q.to_csv(OUT_FILE, sep="\t", index=False, na_rep="NA")
    summary = q.groupby("registration_verdict", dropna=False).size().rename("sections").reset_index()
    summary.to_csv(SUMMARY_FILE, sep="\t", index=False)
    emit(f"DONE: {len(q)} sections; wrote {OUT_FILE}; overlays {qc_dir}")


if __name__ == "__main__":
    main()
