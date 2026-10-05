#!/usr/bin/env python3
"""Phase1D-R audit only; original derived inputs remain untouched."""
import sys, json, hashlib, importlib.util, subprocess, gc, traceback
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import yaml
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'results/phase1/phase1d_r'
sys.path.insert(0,str(ROOT/'scripts'))
import phase0d_color_shortcut_audit as d5

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
 return h.hexdigest()
def tsv(x,name): pd.DataFrame(x).to_csv(OUT/name,sep='\t',index=False,na_rep='NA')
def read(p): return pd.read_csv(ROOT/p,sep='\t')
original=ROOT/'results/phase0/phase0d/registered_he_region_color_features.tsv.gz'
mat=read('results/phase1/region_analysis_matrix.tsv'); crop=read('results/phase1/image_region_crop_audit.tsv')
color=pd.read_csv(original,sep='\t'); feats=d5.COLOR_FEATURES; keys=['sample_id','region_x','region_y']
assert len(feats)==40 and not color.duplicated(keys).any() and mat.region_id.is_unique and crop.region_id.is_unique
joined=mat.merge(color[keys+feats+['n_color_pixels_sampled_qc','black_fraction_sampled_qc']],on=keys,validate='one_to_one')
missing=joined[ joined[feats].isna().all(axis=1)].copy(); assert len(missing)==39
cohort=read('configs/frozen_vannan_primary_cohort.tsv'); cohort=cohort[cohort.primary_inclusion.eq('YES')]
assert len(cohort)==26 and cohort.donor_id.nunique()==19
assert not set(cohort.sample_id)&{'VUILD105MA1','VUILD48LA1'}
manifest=read('results/phase1/morphology_embedding_extraction_manifest.tsv')
meta=read('results/phase0/phase0c_v/registered_he_tiff_metadata.tsv').set_index('sample_id')
d5man=read('results/phase0/phase0d/d5_color_feature_extraction_manifest.tsv').set_index('sample_id')
targets={'epithelial':'epithelial__injury','fibroblast':'fibroblast__fibroblast_activation','macrophage':'macrophage__inflammatory'}
for label,t in targets.items(): missing['target_eligibility_'+label]=missing['n_target_cells__'+t].ge(20)&missing['Y_'+t].notna()
missing=missing.merge(crop[['region_id','x0_px','y0_px','x1_px_exclusive','y1_px_exclusive']],on='region_id',validate='one_to_one')
missing=missing.rename(columns={'x0_px':'crop_x0','y0_px':'crop_y0','x1_px_exclusive':'crop_x1','y1_px_exclusive':'crop_y1','total_cells':'n_cells'})
missing['all_40_features_nan']=True; missing['reason_initial']='UNDETERMINED_BEFORE_SOURCE_AUDIT'
cols=['region_id','sample_id','donor_id',*[ 'target_eligibility_'+k for k in targets],'x_center_um','y_center_um','crop_x0','crop_y0','crop_x1','crop_y1','full_crop_in_bounds','n_cells','tissue_fraction','all_40_features_nan','reason_initial','n_color_pixels_sampled_qc','black_fraction_sampled_qc']
tsv(missing[cols],'d5_missing_regions.tsv')
protected=list((ROOT/'configs').glob('frozen*'))+[ROOT/'configs/phase1_donor_folds.tsv',original,ROOT/'results/phase1/region_analysis_matrix.tsv',ROOT/'results/phase1/image_region_crop_audit.tsv',ROOT/'results/phase1/morphology_embedding_extraction_manifest.tsv']+sorted((ROOT/'results/phase1/embedding_cache').glob('[!._]*.npz'))
protected=[p for p in protected if p.is_file()]
before={str(p.relative_to(ROOT)):sha(p) for p in protected}
audit=[{'path':p,'sha256':s,'expected_sha256':s,'basis':'repair-start immutable snapshot','status':'CAPTURED'} for p,s in before.items()]
gate=read('results/phase1/input_hash_gate.tsv')
for r in gate.itertuples():
 p=ROOT/r.path
 if p.is_file():
  obs=sha(p)
  mutable_report = r.path == 'reports/CURRENT_PROJECT_STATUS.md'
  audit.append({'path':r.path,'sha256':obs,'expected_sha256':r.observed_sha256,'basis':'user-authorized mutable execution report; not scientific input' if mutable_report else 'Phase1A historical gate','status':'AUTHORIZED_REPORT_UPDATE' if mutable_report else ('PASS' if obs==r.observed_sha256 else 'FAIL')})
tsv(audit,'input_hash_audit.tsv')
if any(r['status']=='FAIL' for r in audit): raise RuntimeError('TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH')
commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
source=ROOT/'scripts/phase0d_color_shortcut_audit.py'
tracked=subprocess.run(['git','ls-files','--error-unmatch',str(source.relative_to(ROOT))],cwd=ROOT,capture_output=True).returncode==0
spec={'source_file':str(source.relative_to(ROOT)),'source_sha256':sha(source),'repository_HEAD':commit,'source_tracked_at_HEAD':tracked,'exact_source_commit':'UNTRACKED; HEAD does not attest source content' if not tracked else commit,'functions':['load_image','aggregate_feature','extract_one'],'feature_names':feats,'grid_um':d5.GRID,'scale_um_per_px':d5.SCALE,'stride_px':d5.STRIDE,'sampling_origin':'whole-slide pixel (0,0); NOT crop-local stride','histogram':'min(int(rgb_float32*8),7), 8 bins; counts normalized by nonblack count','rgb':'float32 / dtype maximum, clipped [0,1]','HED':'skimage.color.rgb2hed, exact imported function','brightness':'max(R,G,B)','saturation':'(max-min)/max; zero max gives 0','valid_pixel':'max RGB > 0.02','NaN_policy':'no valid pixels: means/SD/hist bins explicitly NaN via divide where counts>0','crop_definition':'D5 global 150um grid bin based on floor(pixel*0.2125/150); Phikon rounding bounds separately audited'}
(OUT/'d5_feature_specification.yaml').write_text(yaml.safe_dump(spec,sort_keys=False))
# Region reconciliation reads only existing caches, including arrays to ensure readability.
ids=set()
for p in sorted((ROOT/'results/phase1/embedding_cache').glob('[!._]*.npz')):
 with np.load(p,allow_pickle=False) as z:
  assert str(z['inference_device'].item())=='cpu'
  assert str(z['model_revision'].item())=='2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd'
  for k in z.files: a=z[k]
  ids.update(z['region_ids'].astype(str))
rec=mat[['region_id','sample_id','donor_id']].merge(crop[['region_id','full_crop_in_bounds','exclusion_reason']],on='region_id',validate='one_to_one')
rec['embedding_present']=rec.region_id.isin(ids)
rec['reason']=np.where(rec.embedding_present,'EMBEDDING_PRESENT',np.where(rec.full_crop_in_bounds.eq('NO'),'CROP_OUT_OF_BOUNDS','UNEXPLAINED'))
rec['should_have_embedding']=np.where(rec.full_crop_in_bounds.eq('YES'),'YES','NO')
rec['source_file']='results/phase1/image_region_crop_audit.tsv';rec['source_log']='logs/phase1d_embedding_extraction.log';rec['timestamp_if_available']='NOT_RECORDED_PER_REGION'
for label,t in targets.items(): rec['eligible_'+label]=mat['n_target_cells__'+t].ge(20).to_numpy()&mat['Y_'+t].notna().to_numpy()
tsv(rec,'phase0d_to_phase1d_region_reconciliation.tsv');tsv(rec[~rec.embedding_present],'missing_embedding_regions_audit.tsv')
print('INITIAL_GATE_PASS',len(mat),len(ids),'missing_embeddings',int((~rec.embedding_present).sum()),flush=True)
# Deterministic controls covering every admitted section: 4/section (104), no target/performance selection.
controls=pd.concat([g.sample(n=min(4,len(g)),random_state=20260927) for _,g in joined[joined[feats].notna().all(axis=1)].groupby('sample_id',sort=True)])
assert len(controls)>=100
root=[];diag=[];comparisons=[];recomputed=[];original_load=d5.load_image
for i,row in enumerate(manifest.itertuples(),1):
 sm=row.sample_id; image=ROOT/'data/raw/phase0c_v'/d5man.loc[sm,'image_file']
 obs=sha(image); expected=row.image_sha256
 assert expected==d5man.loc[sm,'registered_HE_sha256']
 audit.append({'path':str(image.relative_to(ROOT)),'sha256':obs,'expected_sha256':expected,'basis':'Phase1D + D5 image manifests','status':'PASS' if obs==expected else 'FAIL'})
 tsv(audit,'input_hash_audit.tsv')
 if obs!=expected: raise RuntimeError('TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH')
 print('IMAGE_HASH_PASS',i,26,sm,flush=True)
 arr=original_load(image)
 assert arr.shape[:2]==(int(meta.loc[sm,'height_px']),int(meta.loc[sm,'width_px'])) and str(arr.dtype)==str(meta.loc[sm,'dtype'])
 audit.append({'path':str(image.relative_to(ROOT)),'sha256':obs,'expected_sha256':expected,'basis':f'readable shape={arr.shape} dtype={arr.dtype}; historical metadata match','status':'PASS'})
 d5.load_image=lambda path:arr
 out,_=d5.extract_one(sm,image.name,mat)
 d5.load_image=original_load
 merged=controls[controls.sample_id.eq(sm)][['region_id']+keys+feats].merge(out[keys+feats],on=keys,suffixes=('_original','_recomputed'),validate='one_to_one')
 for r in merged.to_dict('records'):
  for f in feats:
   a=float(r[f+'_original']);b=float(r[f+'_recomputed']);e=abs(a-b)
   comparisons.append({'region_id':r['region_id'],'sample_id':sm,'feature':f,'original':a,'recomputed':b,'absolute_error':e,'relative_error':e/max(abs(a),1e-15)})
 for r in missing[missing.sample_id.eq(sm)].to_dict('records'):
  # D5 sampled pixels use global stride and physical-grid membership.
  ys=np.arange(0,arr.shape[0],d5.STRIDE);xs=np.arange(0,arr.shape[1],d5.STRIDE)
  ys=ys[np.floor(ys*d5.SCALE/d5.GRID).astype(int)==int(r['region_y'])] if 'region_y' in r else ys[np.floor(ys*d5.SCALE/d5.GRID).astype(int)==int(mat.set_index('region_id').loc[r['region_id'],'region_y'])]
  xs=xs[np.floor(xs*d5.SCALE/d5.GRID).astype(int)==int(mat.set_index('region_id').loc[r['region_id'],'region_x'])]
  pix=arr[np.ix_(ys,xs)].reshape(-1,3);rgb=np.clip(pix.astype(np.float32)/255,0,1);valid=rgb.max(axis=1)>0.02
  cp=arr[int(r['crop_y0']):int(r['crop_y1']),int(r['crop_x0']):int(r['crop_x1']),:]
  hed=d5.rgb2hed(rgb.reshape(1,-1,3)); rgbok=bool(np.isfinite(rgb).all()); hedok=bool(np.isfinite(hed).all())
  vals=out[(out.region_x==mat.set_index('region_id').loc[r['region_id'],'region_x'])&(out.region_y==mat.set_index('region_id').loc[r['region_id'],'region_y'])].iloc[0]
  finite=bool(np.isfinite(vals[feats].to_numpy(dtype=float)).all())
  root.append({'region_id':r['region_id'],'crop_valid':cp.size>0,'rgb_valid':rgbok,'hed_valid':hedok,'hist_valid':bool(valid.sum()>0),'first_failure_step':'aggregate_feature division: counts=0; explicit NaN output' if valid.sum()==0 else 'UNKNOWN','exception':'NONE','root_cause_class':'EMPTY_CROP' if valid.sum()==0 else 'UNKNOWN','root_cause_detail':'EMPTY_VALID_PIXEL_SET_AFTER_FROZEN_BLACK_FILTER; geometric crop nonempty','n_sampled_pixels':len(pix),'n_valid_pixels':int(valid.sum()),'repair_class':'RECONSTRUCTED_DERIVED_FEATURE' if finite else 'UNRESOLVED'})
  ds={'region_id':r['region_id'],'shape':str(cp.shape),'dtype':str(cp.dtype),'min':int(cp.min()),'max':int(cp.max()),'mean':float(cp.mean()),'SD':float(cp.std()),'fraction_zero':float((cp==0).mean()),'fraction_255':float((cp==255).mean()),'fraction_white':float((cp.min(axis=2)>=250).mean()),'fraction_black':float((cp.max(axis=2)<=5).mean()),'fraction_nonfinite':float((~np.isfinite(cp)).mean()),'channel_count':cp.shape[2],'empty_crop':cp.size==0,'alpha_channel':False,'read_error':'NONE','coordinate_rounding':'frozen Phikon bounds; D5 global stride bins retained','HED_finite_fraction':float(np.isfinite(hed).mean()),'zero_optical_density':bool((hed==0).all()),'conversion_range_min':float(rgb.min()),'conversion_range_max':float(rgb.max())}
  for j,ch in enumerate('rgb'):
   ds[ch+'_finite_fraction']=float(np.isfinite(rgb[:,j]).mean());ds[ch+'_variance']=float(np.var(rgb[:,j]));ds[ch+'_histogram_support']=int(np.unique(np.minimum((rgb[:,j]*8).astype(int),7)).size)
  diag.append(ds);recomputed.append({'region_id':r['region_id'],**{f:vals[f] for f in feats},'original_nan_flag':True,'recomputation_status':'FINITE' if finite else 'STILL_UNDEFINED_NO_VALID_COLOR_PIXELS'})
 tsv(root,'d5_nan_root_cause.tsv');tsv(diag,'d5_crop_diagnostics.tsv');tsv(comparisons,'d5_recompute_positive_controls.tsv')
 print('RECOMPUTED',sm,'controls',len(merged),'missing',int(missing.sample_id.eq(sm).sum()),flush=True)
 del arr,out,merged;gc.collect()
errors=pd.DataFrame(comparisons);summ=errors.groupby('feature').agg(max_absolute_error=('absolute_error','max'),median_absolute_error=('absolute_error','median'),max_relative_error=('relative_error','max')).reset_index();tsv(summ,'d5_positive_control_fidelity_summary.tsv')
fidelity=bool(np.isfinite(errors.recomputed).all() and np.allclose(errors.original,errors.recomputed,rtol=1e-12,atol=1e-12))
if fidelity: tsv(recomputed,'d5_missing_regions_recomputed.tsv')
unchanged=all(sha(ROOT/p)==s for p,s in before.items());tsv([{'path':p,'before_sha256':s,'after_sha256':sha(ROOT/p),'unchanged':sha(ROOT/p)==s} for p,s in before.items()],'scientific_freeze_assertions.tsv')
reconciled=not rec.reason.eq('UNEXPLAINED').any() and len(ids)==12183 and len(rec)==12260
unresolved=sum(r['repair_class']=='UNRESOLVED' for r in root)
verdict='TECHNICAL_HOLD_D5_REPAIR' if unresolved or not fidelity else 'TECHNICAL_HOLD_REGION_RECONCILIATION'
counts=[]
for label,t in targets.items():
 eligible=mat['n_target_cells__'+t].ge(20)&mat['Y_'+t].notna()&mat.full_crop_in_bounds.eq('YES')
 counts.append(f"{int(eligible.sum())}/{mat.loc[eligible,'donor_id'].nunique()} (pre-color requirement; final universe not frozen)")
fields=[('VERDICT',verdict),('Original D5 SHA256',before[str(original.relative_to(ROOT))]),('Original D5 rows',len(color)),('All-NaN D5 regions',39),('D5 root cause','zero valid globally sampled pixels after frozen black filter; explicit counts=0 NaN outputs'),('Positive-control regions recomputed',len(controls)),('Positive-control fidelity','PASS' if fidelity else 'FAIL'),('D5 regions reconstructed',39-unresolved if fidelity else 0),('D5 regions technically excluded',0),('D5 unresolved',unresolved),('Imputation used','NO'),('Reconstructed D5 SHA256','NOT GENERATED'),('Phase 0D color shortcut regression check','FAIL — gate not reached; no reconstructed input; no changed conclusion claimed'),('Phase 0D regions','12,260'),('Existing Phase 1D embeddings','12,183'),('77-region difference fully reconciled','YES' if reconciled else 'NO'),('Missing embeddings legitimately excluded',int((rec.reason=='CROP_OUT_OF_BOUNDS').sum())),('Embeddings newly recovered',0),('Final Phase 1 region universe','NOT FROZEN — unresolved D5 input'),('Epithelial injury eligible regions/donors',counts[0]),('Fibroblast activation eligible regions/donors',counts[1]),('Macrophage inflammatory eligible regions/donors',counts[2]),('Scientific parameters unchanged','YES' if unchanged else 'NO'),('Phase 1E–1G authorized','NO'),('Recommended next action','STOP; human review of frozen valid-pixel requirement and legitimate exclusion policy; no automatic imputation or threshold changes')]
report='\n'.join(f'{i}. **{k}:** {v}' for i,(k,v) in enumerate(fields,1))
report+='\n\n# Phase 1D-R audit\n\n原D5、矩阵、cohort、folds、缓存和提取manifest均保持原始字节；未运行Phase1模型、残差或生物学图。审计用原D5 `extract_one` 函数重算整节后选择预先确定对照，保持whole-slide stride采样原点。D5使用全球150μm网格内非黑色采样像素统计，并非直接对Phikon crop做局部采样。\n\n`EMPTY_CROP` 在根因表仅表示过滤后统计像素集合为空；不是几何crop为空。几何有效不能使空像素集合的颜色统计变成已定义值。因此即使原函数重算一致，也不能以零、均值、去掉黑色过滤或更改采样填补39个NaN，亦未将这些区域擅自标成合法技术排除。\n\nR8、R9、R13和R14未通过前置gate，未生成重建D5、最终universe或repaired-input manifest。R13上方目标数为冻结20-cell门槛与crop边界后的参考数，不是完整最终交集。R11若全部77个缺失属于out-of-bounds，R12不适用，无需重提取embedding。\n\n代码来源：原D5脚本未被Git追踪；报告HEAD不是该文件生成时commit的证据。精确源文件SHA与代码定义保存在d5_feature_specification.yaml。脚本没有重写颜色定义。\n\n所有审计表位于 `results/phase1/phase1d_r/`；日志 `logs/phase1d_r_color_repair.log`。原D5列数 '+str(len(color.columns))+'，region key唯一。26图像逐个SHA、可读、shape和dtype核对。\n\n**STOP: PHASE1E_G_AUTHORIZED = NO**\n'
(ROOT/'reports/PHASE1D_R_D5_COLOR_REPAIR_AUDIT.md').write_text(report)
tsv(audit,'input_hash_audit.tsv')
print('FINAL',verdict,'fidelity',fidelity,'unresolved',unresolved,'reconciled',reconciled,'unchanged',unchanged,flush=True)
