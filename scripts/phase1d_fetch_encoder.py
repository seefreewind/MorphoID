#!/usr/bin/env python3
"""Fetch and verify the single frozen public pathology encoder fallback."""
from __future__ import annotations

import csv
import hashlib
import json
import platform
from pathlib import Path

import torch
import transformers
from huggingface_hub import snapshot_download
from PIL import __version__ as pillow_version

ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "results/phase1"
MODEL_ID = "owkin/phikon-v2"
REVISION = "2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd"
EXPECTED_SHA256 = "261ae680fa699b3b951597fd57aa19c02ef735805acb104b93af69b36d928569"
MODEL_DIR = P1 / "model_cache" / "owkin_phikon-v2"
CACHE_DIR = P1 / "huggingface_cache"
OUT = P1 / "morphology_encoder_manifest.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    P1.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.parent.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    snapshot_download(repo_id=MODEL_ID, revision=REVISION,
                      local_dir=str(MODEL_DIR), cache_dir=str(CACHE_DIR))
    weight = MODEL_DIR / "model.safetensors"
    if not weight.is_file():
        raise FileNotFoundError(f"expected checkpoint missing: {weight}")
    observed = sha256(weight)
    if observed != EXPECTED_SHA256:
        raise RuntimeError(f"CHECKPOINT_SHA256_MISMATCH: expected={EXPECTED_SHA256}; observed={observed}")
    config = json.loads((MODEL_DIR / "config.json").read_text(encoding="utf-8"))
    processor = json.loads((MODEL_DIR / "preprocessor_config.json").read_text(encoding="utf-8"))
    fields = [
        "encoder_name", "model_id", "revision_commit", "checkpoint_file", "checkpoint_sha256",
        "checkpoint_size_bytes", "source", "license", "fallback_reason", "architecture",
        "pretraining_provenance", "embedding_dimension", "input_resolution_px", "physical_fov_um",
        "patch_pixel_fov", "resize_interpolation", "processor_normalization", "input_color",
        "model_dtype", "extraction_device", "extraction_batch_size", "frozen", "fine_tuning",
        "target_supervision", "target_dependent_patch_selection", "python_version", "torch_version",
        "transformers_version", "Pillow_version", "checkpoint_local_path",
    ]
    row = {
        "encoder_name": "Phikon-v2 (Owkin)", "model_id": MODEL_ID,
        "revision_commit": REVISION, "checkpoint_file": "model.safetensors",
        "checkpoint_sha256": observed, "checkpoint_size_bytes": weight.stat().st_size,
        "source": f"https://huggingface.co/{MODEL_ID}/tree/{REVISION}",
        "license": "Owkin non-commercial license (research use)",
        "fallback_reason": "UNI2-h and UNI weights absent locally and gated; no credentials or access terms accepted on user's behalf; selected one publicly downloadable pathology encoder before Phase 1 predictions",
        "architecture": f"DINOv2 {config['model_type']} ViT-L/16; {config['num_hidden_layers']} layers; {config['num_attention_heads']} heads; {config['hidden_size']} hidden dimensions",
        "pretraining_provenance": "PANCAN-XL public histology collections including TCGA, CPTAC and GTEx; no additional cohort data imported for this pilot",
        "embedding_dimension": config["hidden_size"],
        "input_resolution_px": f"{processor['crop_size']['height']}x{processor['crop_size']['width']}",
        "physical_fov_um": 150,
        "patch_pixel_fov": "150-um grid cell; pixel bounds rounded from frozen 0.2125 um/px; resized to encoder input; no cell-centered/adaptive crop",
        "resize_interpolation": f"shortest edge {processor['size']['shortest_edge']} px; PIL resample code {processor['resample']} (bicubic); center crop {processor['crop_size']['height']}x{processor['crop_size']['width']}",
        "processor_normalization": json.dumps({"rescale_factor": processor["rescale_factor"], "mean": processor["image_mean"], "std": processor["image_std"]}, separators=(",", ":")),
        "input_color": "RGB, converted by frozen processor",
        "model_dtype": "float32",
        "extraction_device": "mps" if torch.backends.mps.is_available() else "cpu",
        "extraction_batch_size": 4,
        "frozen": "YES", "fine_tuning": "NO", "target_supervision": "NO",
        "target_dependent_patch_selection": "NO", "python_version": platform.python_version(),
        "torch_version": torch.__version__, "transformers_version": transformers.__version__,
        "Pillow_version": pillow_version, "checkpoint_local_path": str(weight.relative_to(ROOT)),
    }
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerow(row)
    print(f"ENCODER_CHECKPOINT_PASS model={MODEL_ID}@{REVISION} sha256={observed} bytes={weight.stat().st_size}")


if __name__ == "__main__":
    main()
