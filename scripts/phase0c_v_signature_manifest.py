#!/usr/bin/env python3
"""Hash the versioned, coverage-only MSigDB sources used in Phase 0C-V."""
import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results/phase0/phase0c_v/source_metadata"
OUT = ROOT / "results/phase0/phase0c_v/phase0c_v_signature_source_manifest.tsv"
BASE = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/"
FILES = [
    ("h.all.v2025.1.Hs.symbols.gmt", "MSigDB Hallmark collection"),
    ("c7.all.v2025.1.Hs.symbols.gmt", "MSigDB immunologic signatures collection"),
    ("c2.cp.reactome.v2025.1.Hs.symbols.gmt", "MSigDB Reactome collection"),
    ("c5.go.bp.v2025.1.Hs.symbols.gmt", "MSigDB GO Biological Process collection"),
]

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

rows = []
date = datetime.now(timezone.utc).date().isoformat()
for filename, role in FILES:
    path = SRC / filename
    if not path.is_file() or path.name.startswith("._"):
        raise FileNotFoundError(path)
    rows.append({
        "file_path": path.relative_to(ROOT).as_posix(),
        "source_url": BASE + filename,
        "collection": role,
        "release": "MSigDB 2025.1.Hs",
        "retrieval_date_utc": date,
        "sha256": sha256(path),
        "size_bytes": path.stat().st_size,
        "use": "predefined gene-set coverage only; no scoring or prediction",
    })
with OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(f"Wrote {len(rows)} versioned signature source records: {OUT}")
