#!/usr/bin/env python3
"""Audit fixed 150-um crop support against frozen image bounds."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "results/phase1"
MATRIX = P1 / "region_analysis_matrix.tsv"
IMAGE_META = ROOT / "results/phase0/phase0c_v/registered_he_tiff_metadata.tsv"
FOLDS = ROOT / "configs/phase1_donor_folds.tsv"
LOG = P1 / "region_exclusion_log.tsv"
OUT = P1 / "image_region_crop_audit.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    images = {r["sample_id"]: r for r in read_tsv(IMAGE_META)}
    donor_fold = {r["donor_id"]: r["fold"] for r in read_tsv(FOLDS)}
    regions = read_tsv(MATRIX)
    audit: list[dict[str, object]] = []
    failing: list[dict[str, object]] = []
    for r in regions:
        sample = r["sample_id"]
        m = images[sample]
        mpp = float(r["fixed_scale_um_per_px"])
        gx, gy = int(r["region_x"]), int(r["region_y"])
        # Match the zero-origin 150-um grid and round only at pixel boundaries.
        x0 = int((gx * 150 / mpp) + 0.5)
        y0 = int((gy * 150 / mpp) + 0.5)
        x1 = int(((gx + 1) * 150 / mpp) + 0.5)
        y1 = int(((gy + 1) * 150 / mpp) + 0.5)
        width, height = int(m["width_px"]), int(m["height_px"])
        ok = x0 >= 0 and y0 >= 0 and x1 <= width and y1 <= height and x1 > x0 and y1 > y0
        row = {
            "region_id": r["region_id"], "sample_id": sample, "donor_id": r["donor_id"],
            "fold": donor_fold[r["donor_id"]], "region_x": gx, "region_y": gy, "image_file": r["image_file"],
            "x0_px": x0, "y0_px": y0, "x1_px_exclusive": x1, "y1_px_exclusive": y1,
            "crop_width_px": x1 - x0, "crop_height_px": y1 - y0,
            "target_fov_um": 150, "fixed_scale_um_per_px": mpp,
            "image_width_px": width, "image_height_px": height,
            "full_crop_in_bounds": "YES" if ok else "NO",
            "exclusion_reason": "" if ok else "FROZEN_COORDINATE_SUPPORT_FAILURE",
        }
        audit.append(row)
        if not ok:
            failing.append({
                "region_id": r["region_id"], "sample_id": sample, "donor_id": r["donor_id"],
                "target": "ALL_TARGETS", "observed_target_lineage_cells": "NA",
                "exclusion_stage": "1A", "exclusion_reason": "FROZEN_COORDINATE_SUPPORT_FAILURE",
                "sensitivity_only_10_cell_threshold": "NO",
            })
    write_tsv(OUT, audit, list(audit[0].keys()))
    previous = [r for r in read_tsv(LOG) if not (
        r.get("exclusion_stage") == "1A" and
        r.get("exclusion_reason") == "FROZEN_COORDINATE_SUPPORT_FAILURE")]
    fields = list(previous[0].keys()) if previous else [
        "region_id", "sample_id", "donor_id", "target", "observed_target_lineage_cells",
        "exclusion_stage", "exclusion_reason", "sensitivity_only_10_cell_threshold"]
    write_tsv(LOG, previous + failing, fields)
    primaries = defaultdict(lambda: [0, 0])
    audit_by_id = {r["region_id"]: r for r in audit}
    for r in regions:
        a = audit_by_id[r["region_id"]]
        for target in ("epithelial__injury", "fibroblast__fibroblast_activation",
                       "macrophage__inflammatory"):
            if int(r[f"n_target_cells__{target}"]) >= 20:
                primaries[target][0] += 1
                primaries[target][1] += int(a["full_crop_in_bounds"] == "YES")
    print("CROP_SUPPORT", {k: {"eligible_before_crop": v[0], "eligible_with_full_150um_crop": v[1]}
                          for k, v in primaries.items()})
    print("out_of_bounds_regions", sum(r["full_crop_in_bounds"] == "NO" for r in audit))


if __name__ == "__main__":
    main()
