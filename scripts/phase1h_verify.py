#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Required evidence gates for completed Phase1H outputs, not unit tests."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from PIL import Image
import xml.etree.ElementTree as ET
import phase1e_g_fit_models as m
R=m.ROOT;O=R/'results/phase1h';CK=O/'checkpoints'
def read(n):return pd.read_csv(O/n,sep='\t')
def main():
 c=json.loads((O/'analysis_complete.json').read_text());checks={'analysis_complete':c['status']=='PHASE1H_ANALYSIS_COMPLETE','19_donors_completed':c['loo_donors']==19,'50_seeds_preserved':c['repeated_seeds']==50,'input_hashes_unchanged':read('input_hash_gate.tsv')['pass'].all()}
 original=pd.read_csv(m.P1/'oof_predictions_primary.tsv.gz',sep='\t');original=original.rename(columns={'analysis_stratum':'stratum'})
 folds=m.read_tsv(m.FOLDS).set_index('donor_id').fold.to_dict();loo=read('leave_one_donor_out_influence.tsv');rep=read('repeated_fold_metrics.tsv');part=read('repeated_partition_assignments.tsv')
 checks['loo_three_targets_three_strata_19_donors']=len(loo)==171 and loo.groupby(['target','stratum']).donor_id.nunique().eq(19).all()
 checks['50_partitions_or_explicit_invalid']=len(rep)+len(read('invalid_partitions.tsv'))==300
 checks['partition_donor_unique']=not part.duplicated(['seed','donor_id']).any() and part.groupby('seed').donor_id.nunique().eq(19).all()
 baseline=read('full_frozen_metrics.tsv');pm=m.read_tsv(m.P1/'primary_model_performance.tsv');rm=m.read_tsv(m.P1/'residual_prediction_performance.tsv')
 checks['baseline_matches_original_metrics']=all(np.isclose(b.DeltaR2,pm[pm.target.eq(b.target)&pm.stratum.eq(b.stratum)].deltaR2_morph.iloc[0],rtol=1e-10,atol=1e-10) and np.isclose(b.residual_R2,rm[rm.target.eq(b.target)&rm.stratum.eq(b.stratum)].residual_R2.iloc[0],rtol=1e-10,atol=1e-10) for b in baseline.itertuples())
 checks['original_bootstrap_intervals_reproduced']=all(np.isclose(b.delta_CI_low,pm[pm.target.eq(b.target)&pm.stratum.eq(b.stratum)].deltaR2_morph_CI_low.iloc[0],rtol=1e-8,atol=1e-8) and np.isclose(b.residual_CI_low,rm[rm.target.eq(b.target)&rm.stratum.eq(b.stratum)].CI_low.iloc[0],rtol=1e-8,atol=1e-8) for b in baseline.itertuples())
 unchanged=[]
 for row in loo.itertuples():
  if row.donor_not_in_stratum:continue
  p=CK/(row.run+'__'+folds[row.donor_id]+'.tsv.gz');g=pd.read_csv(p,sep='\t');b=original[original.target.eq(row.target)&original.stratum.eq(row.stratum)].set_index('region_id').loc[g.region_id]
  cols=['M2_pred','M4_pred','heldout_base_residual','morphology_predicted_residual']
  unchanged.append(np.allclose(g[cols].to_numpy(),b[cols].to_numpy(),rtol=1e-7,atol=1e-7))
 checks['loo_unchanged_training_fold_predictions_reproduced']=all(unchanged)
 af=[];crossfitok=True;finite=True;expectedcounts=True
 for p in CK.glob('*_alphas.tsv'):
  if p.name.startswith('._'):continue
  a=pd.read_csv(p,sep='\t');af.append(a)
  tag=p.name[:-len('_alphas.tsv')];tf=p.with_name(tag+'_train_residual.tsv.gz');fp=p.with_name(tag+'.tsv.gz');tr=pd.read_csv(tf,sep='\t');te=pd.read_csv(fp,sep='\t');crossfitok &= set(tr.donor_id).isdisjoint(te.donor_id) and not tr.region_id.duplicated().any() and np.isfinite(tr[['observed','M2_train_oof_pred','residual']].to_numpy(float)).all();crossfitok &= np.allclose(tr.observed-tr.M2_train_oof_pred,tr.residual)
  finite &= np.isfinite(te[[k for k in te if k.endswith('_pred') or k=='heldout_base_residual']].to_numpy(float)).all()
  run=a.run.iloc[0];t=a.target.iloc[0];s=a.stratum.iloc[0];base=original[original.target.eq(t)&original.stratum.eq(s)];drop=run[len('loo_'):].split('__')[0] if run.startswith('loo_') else None
  universe=base[base.donor_id.ne(drop)] if drop else base
  expectedcounts &= len(tr)+len(te)==len(universe) and set(tr.region_id).union(te.region_id)==set(universe.region_id)
  mapping=folds if run.startswith('loo_') else part[part.seed.eq(int(run.split('__')[0].replace('repeat_','')))].set_index('donor_id').fold.to_dict()
  for inner,ig in tr.groupby('crossfit_fold'):crossfitok &= all(mapping[d]==inner for d in ig.donor_id.unique()) and set(ig.donor_id).isdisjoint(set(tr[tr.crossfit_fold.ne(inner)].donor_id))
 a=pd.concat(af);checks['all_alphas_frozen_grid']=np.isclose(a.alpha.to_numpy()[:,None],m.ALPHAS[None,:],rtol=1e-12).any(1).all();checks['all_train_test_donor_overlap_zero']=a.donor_overlap.eq(0).all();checks['strict_training_side_crossfit']=bool(crossfitok);checks['all_oof_predictions_finite']=bool(finite);checks['no_region_losses_inside_runs']=bool(expectedcounts)
 write=a.to_csv(O/'all_alpha_and_split_audit.tsv',sep='\t',index=False)
 # SVG text remains editable and every literal rendered font is >=5 pt.
 svgs=sorted((O/'figures').glob('H*.svg'));pngs=sorted((O/'figures').glob('H*.png'));checks['five_svg_and_png_figures']=len(svgs)==5 and len(pngs)==5
 import re
 sizes=[];export=[]
 for p in svgs:
  root=ET.parse(p).getroot();text=root.findall('.//{http://www.w3.org/2000/svg}text');fonts=[]
  for element in text:
   style=element.attrib.get('style','');matches=re.findall(r'(?:font-size:|font:)\s*([\d.]+)px',style);fonts.extend(float(z) for z in matches)
  sizes.extend(fonts);width=float(root.attrib['width'].replace('pt',''))*25.4/72
  export.append(dict(file=p.name,width_mm=width,editable_text_elements=len(text),minimum_font_pt=min(fonts) if fonts else None))
  checks[p.stem+'_editable_text']=len(text)>0;checks[p.stem+'_width']=abs(width-183)<.1
 checks['svg_glyph_floor']=len(sizes)>0 and min(sizes)>=5
 for p in pngs:
  with Image.open(p) as im:checks[p.stem+'_raster_300dpi']=im.width>=2160 and abs(im.info.get('dpi',(0,0))[0]-300)<1
 (O/'figure_export_audit.json').write_text(json.dumps(export,indent=2))
 out={'checks':{k:bool(v) for k,v in checks.items()},'all_required_checks_pass':bool(all(checks.values())),'visual_review_complete':False,'visual_review_panels':[],'source_validator_export_exception':'PDF intentionally not generated under explicit user output policy; editable SVG font audit used.'}
 (O/'verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert out['all_required_checks_pass']
if __name__=='__main__':main()
