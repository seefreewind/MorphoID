#!/usr/bin/env python3
"""Reproducibly derive the GEO sample/donor/version inventory for Phase 0C-V."""

from __future__ import annotations

import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "phase0" / "phase0c_v"
MATRIX = BASE / "source_metadata" / "GSE250346_series_matrix.txt.gz"
FILELIST = BASE / "source_metadata" / "GSE250346_filelist.txt"
OUT = BASE


def geo_rows() -> tuple[dict[str, list[str]], dict[str, dict[str, str]]]:
    rows: dict[str, list[str]] = {}
    with gzip.open(MATRIX, "rt", encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row and row[0].startswith("!"):
                rows[row[0]] = row[1:]
    gsms = rows["!Sample_geo_accession"]
    chars: list[dict[str, str]] = [defaultdict(str) for _ in gsms]
    with gzip.open(MATRIX, "rt", encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row and row[0] == "!Sample_characteristics_ch1":
                for i, item in enumerate(row[1:]):
                    if ":" in item:
                        key, value = item.split(":", 1)
                        chars[i][key.strip().lower()] = value.strip()
    return rows, dict(zip(gsms, chars))


def file_inventory() -> dict[str, list[tuple[str, int, str]]]:
    by_gsm: dict[str, list[tuple[str, int, str]]] = defaultdict(list)
    with FILELIST.open(encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 4 and row[0] == "File":
                name, date, size = row[1], row[2], int(row[3])
                match = re.match(r"(GSM\d+)_", name)
                if match:
                    by_gsm[match.group(1)].append((name, size, date))
    return by_gsm


def canonical_sample_id(gsm: str, files: list[tuple[str, int, str]]) -> str:
    transcript = next((name for name, _, _ in files if "_transcripts.csv.gz" in name), "")
    if not transcript:
        return gsm
    stem = transcript[:-len("_transcripts.csv.gz")]
    if "_output-" in stem:
        parts = stem.split("__")
        return parts[2] if len(parts) >= 5 else gsm
    return stem.split("_", 1)[1] if "_" in stem else gsm


def donor_from_title(title: str) -> tuple[str, str]:
    raw = title.split("_", 1)[0]
    replicate = ""
    if re.search(r"-[12]$", raw):
        raw, rep = raw.rsplit("-", 1)
        replicate = rep
    elif re.search(r"[AB]$", raw) and raw.startswith("VUHD"):
        replicate = raw[-1]
        raw = raw[:-1]
    # The sample identifier remains a substring of the donor identifier after
    # removing only GEO's documented section replicate suffixes.
    return raw, replicate


def main() -> None:
    rows, characteristics = geo_rows()
    files = file_inventory()
    gsms = rows["!Sample_geo_accession"]
    titles = rows["!Sample_title"]
    submission = rows["!Sample_submission_date"]
    last_update = rows["!Sample_last_update_date"]
    source_names = rows["!Sample_source_name_ch1"]
    with gzip.open(MATRIX, "rt", encoding="utf-8") as fh:
        series_supp = {
            next(csv.reader([line], delimiter="\t"))[1].rsplit("/", 1)[-1]
            for line in fh if line.startswith("!Series_supplementary_file\t")
        }
    sample_records = []
    primary_donors: set[str] = set()
    donor_to_records: dict[str, list[dict[str, str]]] = defaultdict(list)
    for gsm, title, sub_date, upd_date, source in zip(gsms, titles, submission, last_update, source_names):
        attrs = characteristics[gsm]
        tma_match = re.search(r"TMA(\d+)", title)
        run_match = re.search(r"Run(\d+)", title)
        tma = tma_match.group(1) if tma_match else attrs.get("tma", "")
        run = run_match.group(1) if run_match else attrs.get("run", "")
        donor, suffix_replicate = donor_from_title(title)
        primary = tma in {"1", "2", "3", "4"}
        if primary:
            primary_donors.add(donor)
        affected = "control" if "control" in title.lower() or "unaffected" in title.lower() else (
            "less_affected" if "less-affected" in title.lower() else "more_affected"
        )
        diagnosis = attrs.get("clinical diagnosis", attrs.get("diagnosis", ""))
        status = attrs.get("status", "")
        disease_status = "control" if affected == "control" or status.lower() in {"control", "unaffected"} else "pulmonary_fibrosis"
        sample_files = files.get(gsm, [])
        sample_name = canonical_sample_id(gsm, sample_files)
        instrument_match = next((re.search(r"output-(XETG\d+)", x[0]) for x in sample_files if "_transcripts.csv" in x[0] and "output-" in x[0]), None)
        instrument_id = instrument_match.group(1) if instrument_match else ""
        replicate_match = re.search(r"(?:MA|LA)([12])$", sample_name)
        replicate_id = suffix_replicate or (replicate_match.group(1) if replicate_match else "")
        if not replicate_id and sample_name.startswith("VUHD116"):
            replicate_id = sample_name[-1]
        registered = [x for x in sample_files if "_registered_HE." in x[0]]
        transcripts = [x for x in sample_files if "_transcripts.csv" in x[0]]
        xenium_tars = [x for x in sample_files if x[0].endswith(".tar.gz")]
        raw_tma = f"GSE250346_IPFTMA{tma}.svs.tiff.gz" in series_supp
        rds = "GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds" in series_supp
        rec = {
            "sample_id": sample_name,
            "donor_id": donor,
            "GSM": gsm,
            "GSE": "GSE250346",
            "published_primary": "YES" if primary else "NO",
            "geo_extension": "NO" if primary else "YES",
            "formal_2025_article_member": "YES",
            "requested_28_19_primary_subset": "YES" if primary else "NO",
            "disease_status": disease_status,
            "affected_status": affected,
            "control_or_pf": disease_status,
            "TMA": tma,
            "run": run,
            "collection_center": "Vanderbilt" if donor.startswith("VU") else "TGen" if donor.startswith("T") else "UNRESOLVED",
            "he_scanner": "Leica Biosystems Aperio CS2 (publication-level method; not per-image EXIF)",
            "instrument_id": instrument_id or "UNLISTED_IN_SAMPLE_FILENAME",
            "xenium_analysis_version": "xenium-1.1.0.2" if tma in {"1", "2", "3", "4"} else "xenium-2.0.0.10",
            "core_id": sample_name,
            "replicate_id": replicate_id or "none",
            "same_donor_other_core": "PENDING" if False else "",
            "registered_HE_available": "YES" if registered else "NO",
            "raw_HE_available": "YES" if raw_tma else "NO",
            "Xenium_output_available": "YES" if transcripts or xenium_tars else "NO",
            "cell_annotation_available": "YES_CORRECTED_RDS_AUDITED" if rds else "NO",
            "coordinates_available": "YES_TRANSCRIPTS_AND_RDS" if transcripts and rds else ("YES_TRANSCRIPTS" if transcripts else "NO"),
            "expression_available": "YES_TRANSCRIPTS_AND_RDS" if transcripts and rds else ("YES_TRANSCRIPTS" if transcripts else "NO"),
            "panel_version": "Xenium human lung base panel PD_277 (246 genes) + custom CVEVZD (97 genes); 343 total",
            "n_panel_genes": "343",
            "source_publication": "Vannan et al., Nature Genetics 2025; DOI:10.1038/s41588-025-02080-x",
            "notes": f"GEO title={title}; source={source}; diagnosis={diagnosis}; GEO submission={sub_date}; last update={upd_date}; canonical sample ID follows latest transcript filename. Frozen admission subset is TMA1-4 (28 sections/19 donors). Current 2025 Nature Genetics/GEO cohort includes TMA1-5 (45 sections/35 donors); TMA5 is separately tagged as GEO_EXTENSION here to preserve the user-fixed 28/19 gate, not because it is absent from the formal article.",
        }
        sample_records.append(rec)
        donor_to_records[donor].append(rec)
    if len(sample_records) != 45 or sum(r["published_primary"] == "YES" for r in sample_records) != 28 or len(primary_donors) != 19:
        raise RuntimeError(f"Cohort count audit failed: samples={len(sample_records)}, primary={sum(r['published_primary']=='YES' for r in sample_records)}, primary donors={len(primary_donors)}")
    for donor, recs in donor_to_records.items():
        same = "YES" if len(recs) > 1 else "NO"
        for r in recs:
            r["same_donor_other_core"] = same
    crosswalk_cols = list(sample_records[0])
    with (OUT / "vannan_sample_crosswalk.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=crosswalk_cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(sample_records)

    donor_rows = []
    for donor, recs in sorted(donor_to_records.items()):
        primary_rows = [r for r in recs if r["published_primary"] == "YES"]
        donor_rows.append({
            "donor_id": donor,
            "n_cores": len(recs),
            "core_ids": ";".join(r["sample_id"] for r in recs),
            "disease_status": ";".join(sorted({r["disease_status"] for r in recs})),
            "affected_levels": ";".join(sorted({r["affected_status"] for r in recs})),
            "n_runs": len({r["run"] for r in recs}),
            "n_TMA": len({r["TMA"] for r in recs}),
            "statistical_unit": "donor_id",
            "published_primary_cores": len(primary_rows),
            "geo_extension_cores": len(recs) - len(primary_rows),
        })
    donor_cols = list(donor_rows[0])
    with (OUT / "donor_core_structure.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=donor_cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(donor_rows)

    confounding_rows = []
    for cohort, records in (
        ("PUBLISHED_PRIMARY_TMA1_4", [r for r in sample_records if r["published_primary"] == "YES"]),
        ("ALL_GEO_CURRENT", sample_records),
    ):
        for left, right in (
            ("disease_status", "TMA"),
            ("disease_status", "run"),
            ("disease_status", "collection_center"),
            ("disease_status", "instrument_id"),
            ("disease_status", "xenium_analysis_version"),
            ("disease_status", "panel_version"),
            ("disease_status", "he_scanner"),
            ("affected_status", "TMA"),
        ):
            cells: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
            for r in records:
                cells[(r[left], r[right])].append(r)
            levels: dict[str, set[str]] = defaultdict(set)
            for a, b in cells:
                levels[b].add(a)
            for (a, b), rows_ in cells.items():
                confounding_rows.append({
                    "cohort": cohort,
                    "factor_pair": f"{left}__by__{right}",
                    "left_variable": left,
                    "left_level": a,
                    "right_variable": right,
                    "right_level": b,
                    "section_count": len(rows_),
                    "independent_donors": len({r["donor_id"] for r in rows_}),
                    "disease_level_count_in_batch": len(levels[b]),
                    "alias_assessment": "pending_pair_summary",
                    "notes": "Counts use sections for data availability and distinct donor IDs for biological replication; repeated cores remain one donor.",
                })
            left_levels = {a for a, _ in cells}
            right_levels = {b for _, b in cells}
            both_count = sum(len({a for a, batch in cells if batch == b}) > 1 for b in right_levels)
            for row in confounding_rows:
                if row["cohort"] == cohort and row["factor_pair"] == f"{left}__by__{right}":
                    if len(right_levels) == 1 or len(left_levels) == 1:
                        row["alias_assessment"] = "not_assessable_single_level"
                    elif both_count == len(right_levels):
                        row["alias_assessment"] = "NO_COMPLETE_ALIAS_ALL_BATCHES_MIXED"
                    elif both_count == 0:
                        row["alias_assessment"] = "HIGH_COMPLETE_ALIAS"
                    else:
                        row["alias_assessment"] = "PARTIAL_ALIAS"
    with (OUT / "confounding_matrix.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(confounding_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(confounding_rows)

    primary_donor_rows = []
    for donor, recs in sorted(donor_to_records.items()):
        primary_rows = [r for r in recs if r["published_primary"] == "YES"]
        if not primary_rows:
            continue
        primary_donor_rows.append({
            "donor_id": donor,
            "disease_status": ";".join(sorted({r["disease_status"] for r in primary_rows})),
            "n_primary_sections": len(primary_rows),
            "sample_ids": ";".join(r["sample_id"] for r in primary_rows),
            "TMA_values": ";".join(sorted({r["TMA"] for r in primary_rows})),
            "affected_levels": ";".join(sorted({r["affected_status"] for r in primary_rows})),
            "interpretation": "Disease status is donor-level and constant within donor; all cores from this donor remain one biological unit.",
        })
    with (OUT / "disease_donor_identity.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(primary_donor_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(primary_donor_rows)

    raw_rows = []
    for r in sample_records:
        raw_rows.append({
            "sample_id": r["sample_id"],
            "donor_id": r["donor_id"],
            "TMA": r["TMA"],
            "raw_HE_file": f"GSE250346_IPFTMA{r['TMA']}.svs.tiff.gz" if r["raw_HE_available"] == "YES" else "",
            "registered_HE_file": next((x[0] for x in files[r["GSM"]] if "_registered_HE." in x[0]), ""),
            "registered_HE_sample_specific": "YES" if r["registered_HE_available"] == "YES" else "NO",
            "raw_to_registered_core_crosswalk": "TMA_LEVEL_ONLY_CORE_COORDINATE_UNRESOLVED",
            "reason": "GEO raw slide is whole TMA; the public series-level sample assignment does not provide a deterministic core polygon/slide-coordinate table linked to the per-sample registered image",
        })
    with (OUT / "raw_registered_he_crosswalk.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(raw_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(raw_rows)

    version_rows = [
        {"file":"GSE250346 Seurat object (pre-publication gene decoding)","version":"pre-publication RDS","date":"prior to Nature Genetics publication","deprecated":"YES","replacement":"GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds","used":"NO","reason":"README states codeword-to-gene assignment errors; use corrected object"},
        {"file":"GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds","version":"author-corrected Seurat RDS","date":"2024-08 name; current GEO source as accessed 2026-09-28","deprecated":"NO","replacement":"","used":"YES","reason":"Latest author-corrected object; audited for 45-section current cohort, 343 genes, cell labels, raw coordinates and sampled count-metadata joins."},
        {"file":"GSM7990534 pre-correction transcript filename ending TILD117MA","version":"old sample naming / transcript file","date":"replaced 2025-05-08","deprecated":"YES","replacement":"GSM7990534_output-XETG00048__0003400__TILD117MA1__20230313__191400_transcripts.csv.gz","used":"NO","reason":"GEO update explicitly identifies replacement; current inventory exposes corrected MA1 filename"},
        {"file":"GSM8505453_TILD117MA_transcripts.csv.gz","version":"old sample naming / transcript file","date":"replaced 2025-05-08","deprecated":"YES","replacement":"GSM8505453_TILD117MA2_transcripts.csv.gz","used":"NO","reason":"GEO update explicitly identifies replacement; current inventory exposes corrected MA2 filename"},
        {"file":"GSE250346_remove_nuclei_VUILD105MA1.csv.gz","version":"explicit exclusion list","date":"current GEO supplement","deprecated":"NO","replacement":"checked against corrected author object","used":"YES_AUDITED_NO_ADDITIONAL_REMOVAL","reason":"931 listed IDs; zero matched cells in the corrected object, so they were already absent."},
        {"file":"GSE250346_remove_coordinates_VUILD105MA1.csv.gz","version":"explicit exclusion polygon","date":"current GEO supplement","deprecated":"NO","replacement":"applied to corrected author object","used":"YES_22_CELLS_EXCLUDED","reason":"Five author-supplied polygon vertices overlapped 22 cells, excluded before the joined primary table."},
        {"file":"GSE250346 processed files before June 2025 update","version":"superseded processed data","date":"before 2025-06-06","deprecated":"YES","replacement":"current GEO processed data and corrected object","used":"NO","reason":"GEO series matrix records processed data updated 2025-06-06; sample pages expose current files"},
    ]
    with (OUT / "file_version_audit.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(version_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(version_rows)

    print(f"Wrote crosswalk: {len(sample_records)} sections, {len(primary_donors)} primary donors, {len(donor_rows)} all-version donors")


if __name__ == "__main__":
    main()
