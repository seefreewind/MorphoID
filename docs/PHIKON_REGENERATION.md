# Frozen embedding regeneration

Official source owkin/phikon-v2; revision 2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd. model.safetensors SHA256: 261ae680fa699b3b951597fd57aa19c02ef735805acb104b93af69b36d928569.

The caller supplies the official model/processor files locally. Official config files are not copied into this candidate; their expected SHA256 values are recorded below. No model/processor reimplementation is introduced.

- config.json: `31b069a93e6b68741aba703e9b44b01c2e0d00a251023fde2af8e6fd80bc77fc`
- preprocessor_config.json: `dde0ea8d89652cd0c34dd25b5bf68ea9d948defc06edbeefed21aa1c103f2b51`

Frozen crops use the rounded bounds in image_region_crop_audit.tsv at 0.2125 µm/pixel, a 150 µm field, RGB uint8, official processor bicubic resize/center crop to 224×224, rescale 1/255, mean [0.485,0.456,0.406], std [0.229,0.224,0.225]. Inference uses CPU, float32 and batch 4, without fine-tuning. Features are last_hidden_state[:,0,:] CLS vectors of dimension 1024. The official processor config is authoritative.

NPZ schema: region_ids (Unicode ordered IDs), raw (N,1024 float32), grayscale (N,1024 float32), macenko_fold_1 through macenko_fold_5 (N,1024 float32; missing/non-selected entries preserved), model_revision and inference_device scalar strings, tissue_fraction float32. Region ordering follows full-support crop audit rows within each cohort section; no sorting or reindexing is added. The section manifest records section/donor/image hashes, counts and original NPZ hashes. Historic encoder metadata reports the initially selected MPS device; actual completed extraction provenance is CPU in the extraction manifest. Both receipts are retained unchanged.

Run from a separate writable study copy after reconstructing all historical prerequisites:

```sh
python scripts/regenerate_phikon_embeddings.py --model-dir /your/local/official/phikon-v2 --revision 2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd --accepted-provider-terms
```

The revision declaration does not independently prove provenance; provider download receipts should be retained. Hashes validate the locally supplied files. The wrapper does not download, rebundle or relax the scientific hash gates. It uses the original extraction main function and produces local files. Source images, prior hash baselines and historical prerequisites must be mapped/reconstructed; reproduction is not claimed to have been tested on another installation.

Historical analysis sequence: phase1a_build_matrix.R and crop support audit; original phase1d extraction via wrapper; phase1d_r2 common universe audit; phase1e_g_fit_models.py; phase1h_audit.py; phase1i_audit.py; phase1j precision/feasibility; original figure-generation script. These are reproducibility entry points, not authorization to expand the frozen study. Negative R², donor folds and training-side cross-fitting must be preserved.

The archived final region matrix already contains full_crop_in_bounds. The wrapper asserts that this column matches the frozen crop audit and supplies a view without that duplicate column to the original extractor's merge, avoiding suffix ambiguity. No file is rewritten by this input adapter, no rows/values are changed, and the original extraction code remains bytewise unchanged. This adapter was syntax checked, not executed.
