#!/usr/bin/env python3
"""Verify Phase 1A frozen inputs against Phase 0D/Phase 0C provenance.

This gate is read-only with respect to source inputs. It records checksums in
results/phase1/input_hash_gate.tsv and exits non-zero on any mismatch. It does
not regenerate donor folds or construct modeling features.
"""
from __future__ import annotations

import csv
import hashlib
import sys
import time
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase1/input_hash_gate.tsv"
COHORT = ROOT / "configs/frozen_vannan_primary_cohort.tsv"
COHORT_SHA = ROOT / "configs/frozen_vannan_primary_cohort.sha256"
FOLDS = ROOT / "configs/phase1_donor_folds.tsv"
PH0C_INPUTS = ROOT / "results/phase0/phase0c_v/phase0c_v_input_manifest.tsv"
PH0D_IMAGE_MANIFEST = ROOT / "results/phase0/phase0d/d5_color_feature_extraction_manifest.tsv"

FROZEN_ARTIFACTS = [
    "configs/frozen_vannan_primary_cohort.tsv",
    "configs/frozen_vannan_primary_cohort.sha256",
    "configs/phase1_donor_folds.tsv",
    "configs/frozen_covariates_Z.tsv",
    "configs/frozen_image_preprocessing.yaml",
    "configs/frozen_within_lineage_targets.tsv",
    "configs/frozen_phase1_candidate_programs.tsv",
    "configs/frozen_phase1_estimand.yaml",
    "configs/frozen_morphology_scale.yaml",
    "configs/frozen_section_selection.yaml",
    "configs/random_seeds.yaml",
    "configs/current_project_status.yaml",
    "reports/CURRENT_PROJECT_STATUS.md",
    "reports/PHASE0D_METADATA_SHORTCUT_AUDIT.md",
    "results/phase0/phase0d/d5_color_feature_extraction_manifest.tsv",
    "results/phase0/phase0c_v/phase0c_v_input_manifest.tsv",
    "results/phase0/phase0c_vr/source_provenance_manifest.tsv",
    "results/phase0/phase0d/donor_fold_balance.tsv",
    "results/phase0/phase0d/within_disease_feasibility.tsv",
    "results/phase0/phase0d/contingency_alias_audit.tsv",
    "results/phase0/phase0d/color_shortcut_metrics.tsv",
    "results/phase0/phase0d/coordinate_shortcut_metrics.tsv",
    "results/phase0/phase0d/composition_disease_coupling.tsv",
    "results/phase0/phase0d/dependency_edges.tsv",
    "results/phase0/phase0d/spatial_dependence_by_section.tsv",
    "results/phase0/phase0d/microregion_independence_summary.tsv",
    "results/phase0/phase0d/author_pathology_region_category_support.tsv",
    "results/phase0/phase0d/frozen_cohort_summary.tsv",
    "results/phase0/phase0d/composition_metadata_shortcut.tsv",
]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(rows: list[dict[str, object]]) -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["input_id", "path", "expected_sha256", "observed_sha256",
            "comparison_basis", "status", "size_bytes"]
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def row(input_id: str, path: Path, expected: str, observed: str,
        basis: str, status: str) -> dict[str, object]:
    try:
        size = path.stat().st_size
    except OSError:
        size = "NA"
    return {"input_id": input_id, "path": str(path.relative_to(ROOT)),
            "expected_sha256": expected, "observed_sha256": observed,
            "comparison_basis": basis, "status": status, "size_bytes": size}


def main() -> int:
    started = time.time()
    records: list[dict[str, object]] = []
    failures: list[str] = []

    required = [ROOT / p for p in FROZEN_ARTIFACTS] + [COHORT_SHA, FOLDS,
                PH0C_INPUTS, PH0D_IMAGE_MANIFEST]
    for path in required:
        if not path.is_file():
            failures.append(f"missing frozen file: {path.relative_to(ROOT)}")

    status_text = (ROOT / "reports/CURRENT_PROJECT_STATUS.md").read_text(encoding="utf-8") if (ROOT / "reports/CURRENT_PROJECT_STATUS.md").exists() else ""
    project_status = yaml.safe_load((ROOT / "configs/current_project_status.yaml").read_text(encoding="utf-8"))
    if "PASS_PHASE0D_WITH_CONFOUNDING_QUALIFICATION" not in status_text or \
       project_status.get("phase0d_status") != "PASS_PHASE0D_WITH_CONFOUNDING_QUALIFICATION" or \
       project_status.get("phase1_authorized") not in ("YES", True):
        failures.append("current project status does not match the authorized Phase 0D state")

    # Phase 0D did not preserve sidecar digests for every frozen config/result.
    # Capture their current hashes and validate semantic values against the
    # Phase 0D report, frozen result tables, and cross-artifact fold mapping.
    for rel in FROZEN_ARTIFACTS:
        path = ROOT / rel
        if not path.is_file():
            continue
        observed = digest(path)
        records.append(row(f"frozen_artifact:{rel}", path, "NO_PHASE0D_SIDECAR_DIGEST",
                           observed,
                           "current SHA256 captured; semantic/cross-output validation below",
                           "CAPTURED_PHASE0D_BASELINE_NOT_SIDECARRED"))

    # Validate the exact frozen design fields before any large file reads.
    cfg = ROOT / "configs"
    zrows = read_tsv(cfg / "frozen_covariates_Z.tsv")
    primary_z = {r["variable"] for r in zrows if r.get("included_primary") == "YES"}
    if primary_z != {"normalized_x", "normalized_y"}:
        failures.append("primary Z differs from Phase 0D freeze")
    target_rows = read_tsv(cfg / "frozen_within_lineage_targets.tsv")
    primary_targets = {(r["lineage"], r["target"]) for r in target_rows
                       if r.get("coverage_class") == "PRIMARY_USABLE"}
    if primary_targets != {("epithelial", "injury"),
                           ("fibroblast", "fibroblast activation"),
                           ("macrophage", "inflammatory")}:
        failures.append("within-lineage primary target set differs from Phase 0D freeze")
    program_rows = read_tsv(cfg / "frozen_phase1_candidate_programs.tsv")
    if len(program_rows) != 8:
        failures.append("frozen broad-program set is not the Phase 0D eight-program set")
    estimand = yaml.safe_load((cfg / "frozen_phase1_estimand.yaml").read_text(encoding="utf-8"))
    if estimand.get("morphology_unit_um") != 150 or \
       estimand.get("outer_validation") != "donor_held_out" or \
       estimand.get("same_donor_crosses_folds") is not False:
        failures.append("frozen estimand scale or donor-grouped validation differs from Phase 0D")
    prep = yaml.safe_load((cfg / "frozen_image_preprocessing.yaml").read_text(encoding="utf-8"))
    if prep.get("primary", {}).get("normalization") != "none" or \
       prep.get("sensitivity", [{}])[0].get("name") != "Macenko":
        failures.append("frozen image preprocessing differs from Phase 0D")

    cohort_expected = COHORT_SHA.read_text(encoding="utf-8").split()[0]
    cohort_observed = digest(COHORT)
    cohort_status = "PASS" if cohort_expected == cohort_observed else "MISMATCH"
    records.append(row("frozen_cohort", COHORT, cohort_expected, cohort_observed,
                       "Phase 0D frozen .sha256 sidecar", cohort_status))
    if cohort_status != "PASS":
        failures.append("frozen cohort SHA256 mismatch")

    cohort_rows = read_tsv(COHORT)
    admitted = [r for r in cohort_rows if r.get("primary_inclusion") == "YES"]
    admitted_donors = {r.get("donor_id", "") for r in admitted}
    if len(cohort_rows) != 28 or len(admitted) != 26 or len(admitted_donors) != 19:
        failures.append("frozen cohort counts differ from Phase 0D: expected 28/26/19")
    if {r.get("sample_id") for r in cohort_rows if r.get("primary_inclusion") == "NO"} != {
            "VUILD105MA1", "VUILD48LA1"}:
        failures.append("permanent exclusion set differs from Phase 0D")

    # Validate the donor split itself before any outcome or embedding work.
    folds = read_tsv(FOLDS)
    if len(folds) != 19 or len({r.get("donor_id") for r in folds}) != 19 or \
       len({r.get("fold") for r in folds}) != 5:
        failures.append("donor-fold table does not contain 19 unique donors in five folds")
    cohort_disease = {r["donor_id"]: r["disease_status"] for r in admitted}
    for fold_name in sorted({r.get("fold", "") for r in folds}):
        fold_donors = [r["donor_id"] for r in folds if r.get("fold") == fold_name]
        if not fold_donors or any(d not in cohort_disease for d in fold_donors):
            failures.append(f"fold contains unknown donor: {fold_name}")
        elif len({cohort_disease[d] for d in fold_donors}) < 2:
            failures.append(f"frozen fold lacks a disease class: {fold_name}")

    # Verify folds were not regenerated: Phase 0D OOF outputs carry the fold
    # assignment actually used by both color and composition diagnostics.
    donor_fold = {r["donor_id"]: r["fold"] for r in folds}
    phase0d_oof = ROOT / "results/phase0/phase0d/color_shortcut_oof_predictions.tsv"
    if phase0d_oof.is_file():
        oof = read_tsv(phase0d_oof)
        if len({r["sample_id"] for r in oof}) != 26 or \
           any(donor_fold.get(r["donor_id"]) != r["fold"] for r in oof):
            failures.append("frozen donor folds differ from Phase 0D color OOF folds")
    else:
        failures.append("Phase 0D color OOF fold record is missing")
    comp_oof_path = ROOT / "results/phase0/phase0d/composition_disease_oof_by_donor.tsv"
    if comp_oof_path.is_file():
        comp_oof = read_tsv(comp_oof_path)
        if len(comp_oof) != 19 or any(donor_fold.get(r["donor_id"]) != r["fold"] for r in comp_oof):
            failures.append("frozen donor folds differ from Phase 0D composition OOF folds")
    else:
        failures.append("Phase 0D composition OOF fold record is missing")

    source_rows = read_tsv(PH0C_INPUTS)
    source_by_name = {Path(r["file_path"]).name: r for r in source_rows}
    d5_rows = read_tsv(PH0D_IMAGE_MANIFEST)
    expected_samples = {r["sample_id"] for r in admitted}
    d5_samples = {r.get("sample_id", "") for r in d5_rows}
    if len(d5_rows) != 26 or d5_samples != expected_samples:
        failures.append("Phase 0D color image manifest does not match the frozen 26 admitted sections")

    data_inputs: list[tuple[str, Path, str, str]] = []
    for name, kind in [
        ("GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds", "expression_object"),
        ("GSE250346_HE_annotations.tar.gz", "author_annotation_archive"),
    ]:
        source = source_by_name.get(name)
        if source is None:
            failures.append(f"Phase 0C source manifest lacks {name}")
        else:
            data_inputs.append((kind, ROOT / source["file_path"], source["sha256"],
                                "Phase 0C frozen input manifest"))

    # All 26 admitted registered images must agree across Phase 0C and Phase
    # 0D manifests before hashing their bytes.
    for image in d5_rows:
        name = image["image_file"]
        source = source_by_name.get(name)
        if source is None:
            failures.append(f"registered image absent from Phase 0C manifest: {name}")
            continue
        expected = source["sha256"]
        d5_expected = image["registered_HE_sha256"]
        if expected != d5_expected:
            failures.append(f"image SHA differs between Phase 0C/0D manifests: {name}")
        data_inputs.append((f"registered_HE:{image['sample_id']}",
                            ROOT / source["file_path"], expected,
                            "Phase 0C input manifest cross-checked with Phase 0D image manifest"))

    # Do not perform multi-gigabyte reads if any cheap/frozen comparison failed.
    if failures:
        for message in failures:
            records.append({"input_id": "GATE", "path": "NA", "expected_sha256": "NA",
                            "observed_sha256": "NA", "comparison_basis": message,
                            "status": "TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH", "size_bytes": "NA"})
        write_rows(records)
        print("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH")
        print("; ".join(failures))
        return 2

    write_rows(records)

    for input_id, path, expected, basis in data_inputs:
        if not path.is_file():
            observed = "MISSING"
            status = "MISMATCH"
        else:
            observed = digest(path)
            status = "PASS" if expected == observed else "MISMATCH"
        records.append(row(input_id, path, expected, observed, basis, status))
        write_rows(records)
        if status != "PASS":
            print("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH")
            print(f"{input_id}: expected={expected}; observed={observed}")
            return 2
        print(f"PASS {input_id} ({path.stat().st_size} bytes)", flush=True)

    print(f"PASS_PHASE1A_INPUT_HASH_GATE; checks={len(records)}; elapsed_seconds={time.time()-started:.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
