#!/usr/bin/env python3
"""Load the frozen extractor with a user-supplied official local model; no download."""
import argparse, hashlib, importlib.util, os
from pathlib import Path

def sha256(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--model-dir',type=Path,required=True)
    ap.add_argument('--revision',required=True)
    ap.add_argument('--accepted-provider-terms',action='store_true')
    a=ap.parse_args()
    if not a.accepted_provider_terms: ap.error('Independently obtain the model and accept applicable provider terms first.')
    if a.revision!='2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd': ap.error('Frozen revision mismatch')
    expected={'model.safetensors': '261ae680fa699b3b951597fd57aa19c02ef735805acb104b93af69b36d928569', 'config.json': '31b069a93e6b68741aba703e9b44b01c2e0d00a251023fde2af8e6fd80bc77fc', 'preprocessor_config.json': 'dde0ea8d89652cd0c34dd25b5bf68ea9d948defc06edbeefed21aa1c103f2b51'}
    for name,h in expected.items():
        p=a.model_dir/name
        if not p.is_file() or sha256(p)!=h: ap.error('Missing or mismatched official local file: '+name)
    os.environ['HF_HUB_OFFLINE']='1';os.environ['TRANSFORMERS_OFFLINE']='1'
    os.environ['MORPHOID_PHASE1D_DEVICE']='cpu'
    path=Path(__file__).resolve().with_name('phase1d_extract_morphology.py')
    spec=importlib.util.spec_from_file_location('frozen_morphoid_extractor',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.MODEL_DIR=a.model_dir.resolve()
    # The archived final matrix already contains the crop-audit column.
    # Supply an identical, audit-checked in-memory view to the original merge.
    original_read=module.read_tsv
    def read_frozen(path):
        df=original_read(path)
        if Path(path)==module.MATRIX and 'full_crop_in_bounds' in df.columns:
            audit=original_read(module.CROP_AUDIT).set_index('region_id')['full_crop_in_bounds']
            expected=df['region_id'].map(audit)
            if expected.isna().any() or not df['full_crop_in_bounds'].eq(expected).all():
                raise RuntimeError('Archived matrix/crop-audit disagreement')
            return df.drop(columns=['full_crop_in_bounds'])
        return df
    module.read_tsv=read_frozen
    module.main()

if __name__=='__main__': main()
