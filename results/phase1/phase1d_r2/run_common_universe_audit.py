#!/usr/bin/env python3
"""Prospective image eligibility amendment; no Phase1 predictive models."""
from pathlib import Path
import sys,hashlib,json,gc
from datetime import datetime
import numpy as np,pandas as pd,yaml
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'results/phase1/phase1d_r2';sys.path.insert(0,str(ROOT/'scripts'))
import phase0d_color_shortcut_audit as d5
RULE='ZERO_VALID_IMAGE_PIXELS_AFTER_FROZEN_BLACK_FILTER'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def read(p):return pd.read_csv(ROOT/p,sep='\t')
def save(x,p):pd.DataFrame(x).to_csv(OUT/p,sep='\t',index=False,na_rep='NA')
protected=list((ROOT/'configs').glob('frozen*'))+[ROOT/'configs/phase1_donor_folds.tsv',ROOT/'results/phase0/phase0d/registered_he_region_color_features.tsv.gz',ROOT/'reports/PHASE1D_R_D5_COLOR_REPAIR_AUDIT.md',ROOT/'results/phase1/region_analysis_matrix.tsv',ROOT/'results/phase1/image_region_crop_audit.tsv',ROOT/'results/phase1/morphology_embedding_extraction_manifest.tsv',ROOT/'scripts/phase1e_g_fit_models.py',ROOT/'scripts/phase0d_color_shortcut_audit.py']+list((ROOT/'results/phase1/embedding_cache').glob('[!._]*.npz'))
protected=[p for p in protected if p.is_file()];before={str(p.relative_to(ROOT)):sha(p) for p in protected}
m=read('results/phase1/region_analysis_matrix.tsv');c=pd.read_csv(ROOT/'results/phase0/phase0d/registered_he_region_color_features.tsv.gz',sep='\t');crop=read('results/phase1/image_region_crop_audit.tsv');co=read('configs/frozen_vannan_primary_cohort.tsv');co=co[co.primary_inclusion.eq('YES')];folds=read('configs/phase1_donor_folds.tsv');man=read('results/phase1/morphology_embedding_extraction_manifest.tsv')
assert len(m)==12260 and m.region_id.is_unique and len(co)==26 and co.donor_id.nunique()==19
assert not set(co.sample_id)&{'VUILD105MA1','VUILD48LA1'}
keys=['sample_id','region_x','region_y'];F=d5.COLOR_FEATURES
joined=m.merge(c[keys+F+['TMA','run']],on=keys,validate='one_to_one');old=set(joined.loc[joined[F].isna().all(axis=1),'region_id'])
variables={'image pixels':True,'crop coordinates':True,'frozen black threshold':True,'donor_id':False,'disease':False,'target score':False,'cell composition':False,'neighborhood':False,'morphology embedding':False,'model prediction':False,'model error':False,'residual':False}
save([{'variable':k,'used_in_rule':'YES' if v else 'NO','scientific_role':'technical image validity' if v else 'not used in exclusion decision'} for k,v in variables.items()],'exclusion_rule_audit.tsv')
# Selection independent of original NaN list and molecular/metadata outcomes.
def valid_counts(image,regions):
 sampled=image[::d5.STRIDE,::d5.STRIDE,:3]
 dtype_max=float(np.iinfo(image.dtype).max) if np.issubdtype(image.dtype,np.integer) else max(float(np.nanmax(sampled)),1.)
 rgb=np.clip(sampled.astype(np.float32)/dtype_max,0,1)
 yy,xx=np.indices(rgb.shape[:2]);gx=np.floor(xx.ravel().astype(float)*d5.STRIDE*d5.SCALE/d5.GRID).astype(int);gy=np.floor(yy.ravel().astype(float)*d5.STRIDE*d5.SCALE/d5.GRID).astype(int)
 nx=int(np.ceil(image.shape[1]*d5.SCALE/d5.GRID));ny=int(np.ceil(image.shape[0]*d5.SCALE/d5.GRID));inside=(gx<nx)&(gy<ny);flat=gy[inside]*nx+gx[inside]
 black=rgb.reshape(-1,3)[inside].max(axis=1)<=.02
 total=np.bincount(flat,minlength=nx*ny);count=np.bincount(flat[~black],minlength=nx*ny);gi=regions.region_y.astype(int).to_numpy()*nx+regions.region_x.astype(int).to_numpy()
 return total[gi],count[gi]
rows=[];zero=[];imageaudit=[];d5man=read('results/phase0/phase0d/d5_color_feature_extraction_manifest.tsv').set_index('sample_id');met=read('results/phase0/phase0c_v/registered_he_tiff_metadata.tsv').set_index('sample_id');cb=crop.set_index('region_id')
if (OUT/'all_region_image_validity.tsv').exists() and len(read('results/phase1/phase1d_r2/all_region_image_validity.tsv'))==len(m):
 rows=read('results/phase1/phase1d_r2/all_region_image_validity.tsv').to_dict('records');zero=read('results/phase1/phase1d_r2/zero_valid_pixel_regions.tsv').to_dict('records');imageaudit=read('results/phase1/phase1d_r2/registered_image_audit.tsv').to_dict('records')
 assert len(imageaudit)==26 and all(q['PASS'] for q in imageaudit)
 assert all(q['image_sha256']==man.set_index('sample_id').loc[q['sample_id'],'image_sha256'] for q in imageaudit)
 print('RESUME_FULL_IMAGE_RULE_AUDIT_COMPLETE: all 12260 pixel counts already independently recomputed',flush=True)
 remaining=[]
else:remaining=list(man.itertuples())
for i,r in enumerate(remaining,1):
 sm=r.sample_id;p=ROOT/'data/raw/phase0c_v'/d5man.loc[sm,'image_file'];hashval=sha(p)
 assert hashval==r.image_sha256==d5man.loc[sm,'registered_HE_sha256']
 arr=d5.load_image(p);assert arr.dtype==np.uint8 and arr.shape[:2]==(int(met.loc[sm,'height_px']),int(met.loc[sm,'width_px']))
 sub=m[m.sample_id.eq(sm)].copy();total,counts=valid_counts(arr,sub)
 for a,nt,nv in zip(sub.itertuples(),total,counts):
  rows.append({'region_id':a.region_id,'n_total_sampled_pixels':int(nt),'n_black_sampled_pixels':int(nt-nv),'n_valid_pixels':int(nv),'image_input_valid':bool(nv>0)})
  if nv==0:
   b=cb.loc[a.region_id];cp=arr[int(b.y0_px):int(b.y1_px_exclusive),int(b.x0_px):int(b.x1_px_exclusive)]
   zero.append({'region_id':a.region_id,'sample_id':sm,'donor_id':a.donor_id,'crop_bounds':f'{b.x0_px},{b.y0_px},{b.x1_px_exclusive},{b.y1_px_exclusive}','crop_shape':str(cp.shape),'dtype':str(cp.dtype),'min_pixel':int(cp.min()) if cp.size else np.nan,'max_pixel':int(cp.max()) if cp.size else np.nan,'mean_pixel':float(cp.mean()) if cp.size else np.nan,'n_total_pixels':int(nt),'n_black_pixels':int(nt-nv),'n_valid_pixels':int(nv),'n_crop_pixels':int(cp.shape[0]*cp.shape[1]),'original_D5_all_nan':a.region_id in old,'exclusion_reason':RULE})
 imageaudit.append({'sample_id':sm,'image_sha256':hashval,'expected_sha256':r.image_sha256,'shape':str(arr.shape),'dtype':str(arr.dtype),'PASS':True})
 save(rows,'all_region_image_validity.tsv');save(zero,'zero_valid_pixel_regions.tsv');save(imageaudit,'registered_image_audit.tsv')
 print('RULE_REAPPLIED',i,26,sm,'regions',len(sub),'zero',int((counts==0).sum()),flush=True);del arr;gc.collect()
a=pd.DataFrame(rows);new=set(a.loc[~a.image_input_valid,'region_id']);match=new==old
save([{'n_zero_valid_pixels':len(new),'overlap_with_original_39':len(new&old),'unexpected_new_regions':len(new-old),'missing_expected_regions':len(old-new),'exact_set_match':match}],'zero_set_comparison.tsv')
if not match:raise RuntimeError('HOLD_COMMON_UNIVERSE_MISMATCH')
# No embedding inference: read the frozen caches solely to audit presence and finite values.
emb={};section_by_id={};sensitivity_valid={}
for p in sorted((ROOT/'results/phase1/embedding_cache').glob('[!._]*.npz')):
 with np.load(p,allow_pickle=False) as z:
  assert str(z['inference_device'].item())=='cpu' and str(z['model_revision'].item())=='2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd'
  ids=z['region_ids'].astype(str);valid=np.isfinite(z['raw']).all(axis=1);svalid=np.isfinite(z['grayscale']).all(axis=1)
  for k in z.files:
   if k.startswith('macenko_fold_'):svalid &= np.isfinite(z[k]).all(axis=1)
  for rid,v,sv in zip(ids,valid,svalid):
   assert rid not in emb;emb[rid]=bool(v);sensitivity_valid[rid]=bool(sv);section_by_id[rid]=p.stem
x=joined.merge(a,on='region_id',validate='one_to_one');x['coordinate_support_failure']=~x.full_crop_in_bounds.eq('YES');x['zero_valid_image_pixels']=x.n_valid_pixels.eq(0);x['final_excluded']=x.coordinate_support_failure|x.zero_valid_image_pixels
x['exclusion_reason']=x.apply(lambda r:';'.join(([ 'FROZEN_COORDINATE_SUPPORT_FAILURE'] if r.coordinate_support_failure else [])+([RULE] if r.zero_valid_image_pixels else [])),axis=1);x['source_phase']=x.apply(lambda r:';'.join((['PHASE1A/PHASE1D'] if r.coordinate_support_failure else [])+(['PHASE1D_R2'] if r.zero_valid_image_pixels else [])),axis=1)
save(x.loc[x.final_excluded,['region_id','sample_id','donor_id','coordinate_support_failure','zero_valid_image_pixels','final_excluded','exclusion_reason','source_phase']],'final_technical_exclusion_registry.tsv')
C=[k for k in m if k.startswith('p_')];N=[k for k in m if k.startswith('N_p_')];assert len(C)==len(N)==15
coords=['x_center_um','y_center_um','normalized_x','normalized_y','width_um','height_um'];Z=['normalized_x','normalized_y']
x['coordinate_valid']=x.full_crop_in_bounds.eq('YES')&np.isfinite(x[coords]).all(axis=1)
x['image_valid']=x.image_input_valid;x['embedding_valid']=x.region_id.map(emb).fillna(False).astype(bool);x['color_valid']=np.isfinite(x[F]).all(axis=1);x['composition_valid']=np.isfinite(x[C]).all(axis=1);x['neighborhood_valid']=np.isfinite(x[N]).all(axis=1);x['Z_valid']=np.isfinite(x[Z]).all(axis=1)&x.disease.isin(['control','pulmonary_fibrosis'])
x['section_id']=x.sample_id
flags=['coordinate_valid','image_valid','embedding_valid','color_valid','composition_valid','neighborhood_valid','Z_valid']
x['common_universe']=x[flags].all(axis=1)&x.sample_id.isin(co.sample_id)&~x.final_excluded
assert (x.common_universe==~x.final_excluded).all(),'unexpected nontechnical input incompleteness'
fields=['region_id','sample_id','donor_id','section_id','disease','TMA','run',*flags,'common_universe','n_valid_pixels']
save(x[fields],'all_region_common_eligibility.tsv');common=x[x.common_universe].copy();save(common[fields],'final_common_phase1_region_universe.tsv')
assert common.donor_id.nunique()==19,'HOLD_COMMON_UNIVERSE_DONOR_LOSS';assert common.sample_id.nunique()==26
counts={'N_total':len(x),'N_coordinate_excluded':int(x.coordinate_support_failure.sum()),'N_zero_image_excluded':int(x.zero_valid_image_pixels.sum()),'N_overlap':int((x.coordinate_support_failure&x.zero_valid_image_pixels).sum()),'N_unique_excluded':int(x.final_excluded.sum()),'N_common_universe':len(common)};save([counts],'universe_size_check.tsv')
for field in ['donor_id','sample_id','disease','TMA','run']:
 q=x.groupby([field,'common_universe']).size().unstack(fill_value=0).rename(columns={False:'excluded',True:'retained'});q['total']=q.excluded+q.retained;q['excluded_fraction']=q.excluded/q.total;q['share_of_all_excluded']=q.excluded/x.final_excluded.sum();save(q.reset_index(),field+'_retention_distribution.tsv')
# Target eligibility is calculated only after serialization of target-blind common membership.
tr=[]
for label,t in {'epithelial injury':'epithelial__injury','fibroblast activation':'fibroblast__fibroblast_activation','macrophage inflammatory':'macrophage__inflammatory'}.items():
 eligible=x['n_target_cells__'+t].ge(20)&x['Y_'+t].notna()&x.donor_id.isin(folds.donor_id);previous=x[eligible&x.full_crop_in_bounds.eq('YES')];now=x[eligible&x.common_universe]
 assert now.region_id.map(sensitivity_valid).all(),'target-eligible sensitivity embedding missing'
 tr.append({'target':label,'eligible_regions':len(now),'eligible_donors':now.donor_id.nunique(),'eligible_sections':now.sample_id.nunique(),'pre_amendment_regions':len(previous),'pre_amendment_donors':previous.donor_id.nunique(),'difference_vs_pre_amendment':len(now)-len(previous),'lineage_min_cells':20})
assert all(r['eligible_donors']==19 for r in tr);save(tr,'primary_target_eligibility.tsv')
# Only the two explicitly authorized Phase0D diagnostics; preserve original section means, folds, C and weights.
base=read('results/phase0/phase0d/color_shortcut_metrics.tsv').set_index('target');diagnostics=[];oofs=[];foldmap=folds.set_index('donor_id').fold
for variant,frame in [('original',c),('common',c.merge(common[keys],on=keys,how='inner',validate='one_to_one'))]:
 regs=frame.copy();regs['TMA_label']=regs.TMA.astype(str).map(lambda a:a if a.startswith('TMA') else 'TMA'+a);regs['run_label']=regs.run.astype(str).map(lambda a:a if a.startswith('Run') else 'Run'+a);regs['disease_label']=regs.disease_status.astype(str).str.lower().map(lambda a:'control' if a in {'control','unaffected'} else 'pulmonary_fibrosis')
 sec=regs.groupby(['sample_id','donor_id','TMA_label','run_label','disease_label'],as_index=False)[F].mean();sec['n_regions']=sec.sample_id.map(regs.groupby('sample_id').size()).astype(int)
 for target in ['TMA','disease']:
  oof,result=d5.section_level_predictions(sec,target+'_label',F,foldmap,None,1500);result['variant']=variant;result['target']=target;result['risk']='HIGH' if result['balanced_accuracy']>=.70 or result['macro_auroc']>=.75 else 'NOT_HIGH_BY_PRESET_THRESHOLD';diagnostics.append(result);oof['variant']=variant;oof['diagnostic_target']=target;oofs.append(oof);print('PHASE0D_DIAGNOSTIC',variant,target,result['balanced_accuracy'],result['macro_auroc'],flush=True)
dg=pd.DataFrame(diagnostics);save(dg,'color_shortcut_regression_check.tsv');save(pd.concat(oofs,ignore_index=True),'color_shortcut_oof_predictions.tsv')
assert all(np.isclose(dg[(dg.variant=='original')&(dg.target==t)].iloc[0].macro_auroc,base.loc[t,'macro_auroc'],atol=1e-12) for t in ['TMA','disease']),'original diagnostic reproduction mismatch'
high=dg[dg.variant.eq('common')].risk.eq('HIGH').all()
unchanged=all(sha(ROOT/p)==h for p,h in before.items());save([{'path':p,'before_sha256':h,'after_sha256':sha(ROOT/p),'unchanged':sha(ROOT/p)==h} for p,h in before.items()],'scientific_freeze_assertions.tsv');assert unchanged
universefile=OUT/'final_common_phase1_region_universe.tsv';uh=sha(universefile);verdict='PHASE1D_R2_COMMON_UNIVERSE_PASS' if high else 'HOLD_COMMON_UNIVERSE_ALTERS_SHORTCUT_CONCLUSION';auth=bool(high)
settings={'protocol_amendment_id':'PHASE1D_R2','rule':RULE,'introduced_before_phase1_modeling':True,'deltaR2_seen_before_rule':False,'residual_results_seen_before_rule':False,'model_performance_seen_before_rule':False,'original_D5_modified':False,'imputation_used':False,'common_universe_file':str(universefile.relative_to(ROOT)),'common_universe_sha256':uh,'n_regions_total':len(x),'n_regions_excluded_coordinate':counts['N_coordinate_excluded'],'n_regions_excluded_zero_image':counts['N_zero_image_excluded'],'n_regions_exclusion_overlap':counts['N_overlap'],'n_regions_excluded_unique':counts['N_unique_excluded'],'n_regions_common':len(common),'n_donors_common':19,'n_sections_common':26,'original_D5_path':'results/phase0/phase0d/registered_he_region_color_features.tsv.gz','original_D5_sha256':before['results/phase0/phase0d/registered_he_region_color_features.tsv.gz'],'membership_applies_to':'all M0-M4 and directly compared shortcut/sensitivity/residual arms; targets gated afterward','black_threshold':0.02,'stride_px':8,'physical_grid_um':150,'fixed_scale_um_per_px':.2125,'sampling_origin':'original D5 whole-slide zero-origin sampling; not a new crop-local sampling origin','prior_hold':'PHASE1D_R=TECHNICAL_HOLD_D5_REPAIR','PHASE1E_G_AUTHORIZED':auth,'date':'2026-10-02 Asia/Shanghai'}
if auth:(ROOT/'configs/phase1_common_universe.yaml').write_text(yaml.safe_dump(settings,sort_keys=False))
else:(OUT/'unfrozen_common_universe_candidate.yaml').write_text(yaml.safe_dump(settings,sort_keys=False))
commonmetrics=dg[dg.variant.eq('common')].set_index('target');dis=commonmetrics.loc['disease'];tma=commonmetrics.loc['TMA']
values=[('VERDICT',verdict),('Protocol amendment rule',RULE),('Outcome-blind rule','YES'),('Phase 0D total regions','12,260'),('Coordinate-support exclusions',counts['N_coordinate_excluded']),('Zero-valid-image exclusions',len(new)),('Overlap between exclusion classes',counts['N_overlap']),('Unique technical exclusions',counts['N_unique_excluded']),('Final common universe',len(common)),('Retained sections',26),('Retained donors',19),('Original D5 modified','NO'),('Imputation used','NO'),('Original 39 independently reproduced','YES'),('Epithelial injury eligible regions/donors',f"{tr[0]['eligible_regions']}/{tr[0]['eligible_donors']}"),('Fibroblast activation eligible regions/donors',f"{tr[1]['eligible_regions']}/{tr[1]['eligible_donors']}"),('Macrophage inflammatory eligible regions/donors',f"{tr[2]['eligible_regions']}/{tr[2]['eligible_donors']}"),('Color-only TMA macro-AUROC',tma.macro_auroc),('Color-only disease BA/AUROC',f'{dis.balanced_accuracy}/{dis.macro_auroc}'),('HIGH shortcut conclusion preserved','YES' if high else 'NO'),('Scientific parameters unchanged','YES'),('Common-universe SHA256',uh),('Phase 1E–1G authorized','YES' if auth else 'NO'),('Recommended next action','STOP; wait for explicit human launch of Phase1E-G' if auth else 'STOP; human review')]
report='\n'.join(f'{i}. **{k}:** {v}' for i,(k,v) in enumerate(values,1))+'\n\n# Phase 1D-R2 protocol amendment\n\n'
report+='这是输入完整性审计后、任何Phase1预测结果出现前加入的技术eligibility amendment；不声称该排除在原Phase0D方案中已存在。历史PHASE1D_R=TECHNICAL_HOLD_D5_REPAIR保留且报告字节未改。新规则不生成或填补颜色值，original D5的39行NaN保持原样。\n\n'
report+='独立重算使用全部12260个区域的原D5 whole-slide stride=8采样、0.2125μm/px、150μm网格、float32/dtype最大值归一化与max(RGB)<=0.02黑色过滤。没有把局部crop重置采样原点。只将计数=0作为图像输入无效，没有新增第二阈值。全量计数表与39个crop证据保留。\n\n'
report+='final_common_phase1_region_universe.tsv仅含common_universe=TRUE的保留区域；全12260行资格记录位于all_region_common_eligibility.tsv。共同集使用C、N、Z与图像输入完整性，不使用lineage门槛或target score，之后再按冻结20-cell和非缺失目标规则计算target资格。所有可直接比较的model arms必须使用同一membership。Primary raw morphology embedding决定共同集；grayscale/Macenko缓存按原Phase1D目标合格union提取，因此非目标合格区域可有空的sensitivity向量，不能据此在target-blind共同集中预删区域。三个目标合格行已逐一验证所有sensitivity向量可用。当前模型脚本尚未接入该新manifest，未来人工启动时须先校验hash并应用membership，不能直接按旧脚本运行。\n\n'
report+=pd.DataFrame(tr).to_markdown(index=False)+'\n\n## 排除结构\n\n'
for field in ['disease','TMA','donor_id','sample_id']:
 q=pd.read_csv(OUT/(field+'_retention_distribution.tsv'),sep='\t');report+=f'### {field}\n\n'+q.to_markdown(index=False)+'\n\n'
report+='**EXCLUSION_STRUCTURE_WARNING:** 39个zero-image区域集中于3个section；这是注册H&E图像的局部黑色区域分布。记录结构性集中，不据此改变规则。整体技术排除的donor/disease/TMA计数与比例见上表及TSV；未按任何标签设定eligibility。\n\n## Phase0D regression check\n\n'
report+=dg[['variant','target','n_regions','n_sections','n_donors','balanced_accuracy','macro_auroc','risk']].to_markdown(index=False)+'\n\n'
report+='诊断完全调用原section_level_predictions：section平均、donor-held-out frozen folds、训练侧StandardScaler、C=1/max_iter=500/class_weight=balanced logistic、donor-equal权重、1500 donor bootstrap。只重算TMA和disease；未重调参、未运行run预测或新permutation分析。原point metrics通过与历史表核对。此处是授权的Phase0D捷径诊断，不是Phase1 M0–M4或分子预测。\n\n'
report+='冻结文件、模型脚本、缓存与原HOLD报告before/after SHA相同，见scientific_freeze_assertions.tsv；新增唯一科学eligibility政策如上。报告与新configs manifest保留原HOLD历史。没有计算DeltaR²、residual、PF-only生物模型或图，半小时监测保持取消。\n\n**STOP: PHASE1E_G_AUTHORIZED = '+('YES' if auth else 'NO')+'; WAIT_FOR_HUMAN_LAUNCH.**\n'
(ROOT/'reports/PHASE1D_R2_COMMON_UNIVERSE_AMENDMENT.md').write_text(report)
methods='''# Phase 1D-R2 methods amendment

Before any Phase 1 predictive modeling, regions with zero valid H&E pixels after application of the pre-specified black-pixel filter were classified as technical image-input failures and excluded from the common analysis universe. This rule depended only on frozen image pixels and crop coordinates and did not use molecular outcomes, disease labels, donor identity, model predictions, or residuals. No imputation or modification of the frozen color-feature definition was performed.

This eligibility rule was introduced as a protocol amendment after the input-completeness audit and before any Phase 1 predictive results were obtained; it was not part of the original Phase 0D exclusion policy. Valid sampled pixels were counted using the original whole-slide sampling origin and stride, physical grid assignment, RGB normalization and black-pixel threshold. The common universe was applied identically to all directly compared model arms; target eligibility was evaluated afterward using the unchanged lineage-cell threshold and target-availability criteria.

'''+f'Of 12,260 regions, {counts["N_unique_excluded"]} unique regions were technically excluded ({counts["N_coordinate_excluded"]} coordinate-support failures and {len(new)} zero-valid-image regions; overlap {counts["N_overlap"]}). The common universe contained {len(common)} regions from 26 sections and 19 donors. Original color rows and morphology caches were preserved.\n'
(ROOT/'reports/PHASE1D_R2_METHODS_AMENDMENT.md').write_text(methods)
status=ROOT/'reports/CURRENT_PROJECT_STATUS.md';s=status.read_text();s+='\n\n## 2026-10-02 Phase 1D-R2 execution update\n\n- Historical `PHASE1D_R = TECHNICAL_HOLD_D5_REPAIR` remains unchanged.\n- `PHASE1D_R2 = '+verdict+'`。\n- `PHASE1E_G_AUTHORIZED = '+('YES' if auth else 'NO')+'`，仅表示新共同集audit资格；未人工启动，当前仍停止。\n- common universe = '+str(len(common))+' regions / 26 sections / 19 donors。\n- 新技术规则 `'+RULE+'`；原D5和39行NaN不变。\n- 详见PHASE1D_R2_COMMON_UNIVERSE_AMENDMENT.md；半小时监测保持取消。\n';status.write_text(s)
print('FINAL',verdict,counts,'SHA256',uh,'PHASE1E_G_AUTHORIZED',auth,flush=True)
