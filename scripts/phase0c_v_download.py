#!/usr/bin/env python3
"""Download the bounded, corrected Phase 0C-V source set from GEO.

This intentionally excludes the 163.7-GB GEO RAW archive and all HESCAPE data.
Only TMA1-4 registered H&E files are downloaded for primary-cohort auditing;
TMA5 raw whole-slide scans are included only for raw/registered image mapping.
"""

from __future__ import annotations

import concurrent.futures
import csv
import datetime as dt
import gzip
import hashlib
import os
import re
import time
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECT = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250346/suppl/"
RAW_DIR = ROOT / "data" / "raw" / "phase0c_v"
OUT_DIR = ROOT / "results" / "phase0" / "phase0c_v"
FILELIST = OUT_DIR / "source_metadata" / "GSE250346_filelist.txt"
MATRIX = OUT_DIR / "source_metadata" / "GSE250346_series_matrix.txt.gz"
MANIFEST = OUT_DIR / "phase0c_v_input_manifest.tsv"
LOG = OUT_DIR / "download.log"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_primary_gsms() -> set[str]:
    with gzip.open(MATRIX, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession\t"):
                gsm = next(csv.reader([line], delimiter="\t"))[1:]
                break
        else:
            raise RuntimeError("GEO series matrix has no Sample_geo_accession row")

        fh.seek(0)
        for line in fh:
            if line.startswith("!Sample_title\t"):
                titles = next(csv.reader([line], delimiter="\t"))[1:]
                break
        else:
            raise RuntimeError("GEO series matrix has no Sample_title row")
    if len(gsm) != len(titles):
        raise RuntimeError("GEO sample and title row lengths differ")
    primary = set()
    for accession, title in zip(gsm, titles):
        match = re.search(r"TMA(\d+)", title)
        if match and match.group(1) in {"1", "2", "3", "4"}:
            primary.add(accession)
    if len(primary) != 28:
        raise RuntimeError(f"Expected 28 TMA1-4 primary sections, found {len(primary)}")
    return primary


def inventory() -> list[tuple[str, str, int | None, str, str]]:
    primary = read_primary_gsms()
    items: list[tuple[str, str, int | None, str, str]] = []
    for line in FILELIST.read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split("\t")
        if len(fields) < 4 or fields[0] != "File":
            continue
        name, size = fields[1], int(fields[3])
        match = re.match(r"(GSM\d+)_", name)
        if match and match.group(1) in primary and "_registered_HE." in name:
            gsm = match.group(1)
            bucket = re.match(r"(GSM\d{4})\d{3}$", gsm).group(1) + "nnn"
            url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{bucket}/{gsm}/suppl/{name}"
            items.append((name, url, size, gsm, "author-registered same-section H&E"))
    series_names: set[str] = set()
    with gzip.open(MATRIX, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("!Series_supplementary_file\t"):
                row = next(csv.reader([line], delimiter="\t"))
                series_names.add(row[1].rsplit("/", 1)[-1])
    wanted_series = {
        "GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds",
        "GSE250346_HE_annotations.tar.gz",
        "GSE250346_remove_coordinates_VUILD105MA1.csv.gz",
        "GSE250346_remove_nuclei_VUILD105MA1.csv.gz",
        *(f"GSE250346_IPFTMA{i}.svs.tiff.gz" for i in range(1, 6)),
    }
    for name in sorted(series_names & wanted_series):
        role = (
            "corrected author cell object" if name.endswith(".rds")
            else "author pathology annotations" if "annotations" in name
            else "explicitly corrected coordinate exclusion list" if "remove_coordinates" in name
            else "explicitly corrected nucleus exclusion list" if "remove_nuclei" in name
            else "raw TMA whole-slide H&E"
        )
        items.append((name, PROJECT + name, None, "GSE250346", role))
    got = sum("author-registered" in role for *_, role in items)
    if got != 28:
        raise RuntimeError(f"Expected registered H&E for all 28 primary sections, found {got}")
    expected_names = {
        "GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds",
        "GSE250346_HE_annotations.tar.gz",
        "GSE250346_remove_coordinates_VUILD105MA1.csv.gz",
        "GSE250346_remove_nuclei_VUILD105MA1.csv.gz",
        *(f"GSE250346_IPFTMA{i}.svs.tiff.gz" for i in range(1, 6)),
    }
    found_names = {name for name, *_ in items}
    if not expected_names <= found_names:
        raise RuntimeError(f"Missing expected bounded inputs: {sorted(expected_names - found_names)}")
    return items


def head_length(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "MorphoID-Phase0C-V/1.0"})
    with urllib.request.urlopen(req, timeout=120) as response:
        value = response.headers.get("Content-Length")
    if value is None:
        raise IOError(f"no Content-Length in HEAD response for {url}")
    return int(value)


def download_one(item: tuple[str, str, int | None, str, str]) -> tuple[str, int, str, str]:
    name, url, expected, accession, role = item
    final = RAW_DIR / name
    if expected is None:
        expected = head_length(url)
    if final.exists() and final.stat().st_size == expected:
        return name, final.stat().st_size, sha256(final), "reused"
    part = final.with_suffix(final.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(1, 4):
        offset = part.stat().st_size if part.exists() else 0
        h = hashlib.sha256()
        if offset:
            with part.open("rb") as old:
                for block in iter(lambda: old.read(8 * 1024 * 1024), b""):
                    h.update(block)
        headers = {"User-Agent": "MorphoID-Phase0C-V/1.0"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        req = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=300) as response:
                if offset and response.status == 206:
                    content_range = response.headers.get("Content-Range", "")
                    if not content_range.startswith(f"bytes {offset}-"):
                        raise IOError(f"invalid resume range for {name}: {content_range}")
                    mode = "ab"
                else:
                    offset = 0
                    h = hashlib.sha256()
                    mode = "wb"
                with part.open(mode) as out:
                    while True:
                        block = response.read(8 * 1024 * 1024)
                        if not block:
                            break
                        out.write(block)
                        h.update(block)
            total = part.stat().st_size
            if total != expected:
                raise IOError(f"size mismatch for {name}: expected {expected}, got {total}; partial retained for resume")
            os.replace(part, final)
            return name, total, h.hexdigest(), "downloaded"
        except Exception as exc:
            last_error = exc
            if attempt < 3:
                time.sleep(2 * attempt)
    raise IOError(f"{last_error}; partial_bytes={part.stat().st_size if part.exists() else 0}")


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    items = inventory()
    downloaded: dict[str, tuple[str, int, str, str]] = {}
    errors: list[str] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        future_to_item = {executor.submit(download_one, item): item for item in items}
        for future in concurrent.futures.as_completed(future_to_item):
            item = future_to_item[future]
            try:
                downloaded[item[0]] = future.result()
            except Exception as exc:
                errors.append(f"{item[0]}\t{exc}")
    with LOG.open("w", encoding="utf-8") as log:
        log.write(f"UTC completed: {dt.datetime.now(dt.timezone.utc).isoformat()}\n")
        log.write(f"Selected inputs: {len(items)}; successful: {len(downloaded)}; failed: {len(errors)}\n")
        for error in errors:
            log.write(f"ERROR\t{error}\n")
    if errors:
        raise RuntimeError(f"{len(errors)} GEO source downloads failed; see {LOG}")

    with MANIFEST.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["file_path", "source_url", "accession", "version", "download_date_utc", "sha256", "size_bytes", "role", "status"])
        for name, url, expected, accession, role in items:
            result_name, size, digest, status = downloaded[name]
            writer.writerow([
                str((Path("data/raw/phase0c_v") / result_name).as_posix()),
                url,
                accession,
                "GEO current as accessed 2026-09-28; file-list size verified",
                dt.datetime.now(dt.timezone.utc).date().isoformat(),
                digest,
                size,
                role,
                status,
            ])
        for path in [FILELIST, MATRIX, OUT_DIR / "source_metadata" / "GSE250346_RDS_README.txt"]:
            writer.writerow([
                str(path.relative_to(ROOT).as_posix()),
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250346/",
                "GSE250346",
                "GEO current as accessed 2026-09-28",
                dt.datetime.now(dt.timezone.utc).date().isoformat(),
                sha256(path),
                path.stat().st_size,
                "GEO source inventory/metadata",
                "downloaded",
            ])
    total = sum(result[1] for result in downloaded.values())
    print(f"Downloaded/reused {len(downloaded)} bounded source files ({total / 1e9:.2f} GB); manifest={MANIFEST}; log={LOG}")


if __name__ == "__main__":
    main()
