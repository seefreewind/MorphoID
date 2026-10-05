#!/usr/bin/env python3
"""Read TIFF headers from author-registered H&E files without refitting transforms."""
import csv
import gzip
import io
import re
from fractions import Fraction
from pathlib import Path

import tifffile

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/phase0/phase0c_v"
MANIFEST = BASE / "phase0c_v_input_manifest.tsv"
OUT = BASE / "registered_he_tiff_metadata.tsv"

with MANIFEST.open(encoding="utf-8") as handle:
    assets = [r for r in csv.DictReader(handle, delimiter="\t")
              if r["role"] == "author-registered same-section H&E"]

rows = []
for asset in assets:
    p = ROOT / asset["file_path"]
    record = {"sample_id": "", "GSM": asset["accession"], "file": asset["file_path"],
              "size_bytes": p.stat().st_size, "sha256": asset["sha256"], "status": "OK"}
    name_match = re.match(r"(GSM\d+)_", p.name)
    record["GSM"] = name_match.group(1) if name_match else asset["accession"]
    try:
        # Tifffile seeks to EOF while opening a stream. Seeking to EOF through
        # gzip decompresses the entire multi-hundred-MB image just to inspect
        # its first IFD. The TIFF directory and scalar metadata live at the
        # start of these files, so parse a bounded 2 MiB prefix instead.
        with gzip.open(p, "rb") as stream:
            prefix = io.BytesIO(stream.read(2 * 1024 * 1024))
        # ImageJ's optional trailing metadata block is outside this bounded
        # prefix; TIFF scalar tags are complete, so skip parsing that block.
        with tifffile.TiffFile(prefix, is_imagej=False) as tf:
            page = tf.pages[0]
            tags = page.tags
            width, height = int(page.imagewidth), int(page.imagelength)
            desc = str(tags["ImageDescription"].value) if "ImageDescription" in tags else ""
            xres = tags["XResolution"].value if "XResolution" in tags else None
            yres = tags["YResolution"].value if "YResolution" in tags else None
            unit = str(tags["ResolutionUnit"].value) if "ResolutionUnit" in tags else ""
            orient = str(tags["Orientation"].value) if "Orientation" in tags else "1"
            if hasattr(xres, "numerator"):
                xres_v = float(Fraction(xres.numerator, xres.denominator))
            elif isinstance(xres, (tuple, list)) and len(xres) == 2:
                xres_v = float(xres[0]) / float(xres[1])
            else:
                xres_v = float(xres) if xres is not None else None
            if hasattr(yres, "numerator"):
                yres_v = float(Fraction(yres.numerator, yres.denominator))
            elif isinstance(yres, (tuple, list)) and len(yres) == 2:
                yres_v = float(yres[0]) / float(yres[1])
            else:
                yres_v = float(yres) if yres is not None else None
            unit_match = re.search(r"(?m)^unit=(.+)$", desc)
            spatial_unit = unit_match.group(1).strip() if unit_match else "UNSPECIFIED"
            px_x = 1.0 / xres_v if xres_v and spatial_unit.lower() in {"micron", "microns", "um"} else None
            px_y = 1.0 / yres_v if yres_v and spatial_unit.lower() in {"micron", "microns", "um"} else None
            record.update({
                "n_pages": len(tf.pages), "width_px": width, "height_px": height,
                "axes": str(tf.series[0].axes), "dtype": str(tf.series[0].dtype),
                "description": desc.replace("\n", "; ")[:500],
                "x_resolution_raw": str(xres), "y_resolution_raw": str(yres),
                "resolution_unit_tag": unit, "ImageJ_spatial_unit": spatial_unit,
                "physical_pixel_size_x_um": px_x, "physical_pixel_size_y_um": px_y,
                "physical_extent_x_um": width * px_x if px_x else None,
                "physical_extent_y_um": height * px_y if px_y else None,
                "orientation_tag": orient, "compression_tag": str(tags["Compression"].value) if "Compression" in tags else "",
                "photometric_tag": str(tags["PhotometricInterpretation"].value) if "PhotometricInterpretation" in tags else "",
            })
    except Exception as exc:
        record.update({"status": "ERROR", "error": repr(exc)})
    rows.append(record)

# Recover canonical section IDs only through the accession-based, versioned sample crosswalk.
with (BASE / "vannan_sample_crosswalk.tsv").open(encoding="utf-8") as handle:
    xwalk = {r["GSM"]: r["sample_id"] for r in csv.DictReader(handle, delimiter="\t")}
for row in rows:
    row["sample_id"] = xwalk.get(row["GSM"], "UNRESOLVED_GSM")

with OUT.open("w", newline="", encoding="utf-8") as handle:
    fields = list(dict.fromkeys(k for row in rows for k in row))
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
print(f"Parsed {len(rows)} registered H&E TIFF headers: {OUT}")
