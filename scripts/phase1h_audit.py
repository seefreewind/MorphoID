#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prospectively frozen donor/fold audit; no changes to primary outputs."""
from pathlib import Path
import json, hashlib, time
import numpy as np
import pandas as pd
import yaml
from threadpoolctl import threadpool_limits
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.covariance import LedoitWolf
from sklearn.neighbors import NearestNeighbors
import phase1e_g_fit_models as m
ROOT=m.ROOT; OUT=ROOT/'results/phase1h'; CK=OUT/'checkpoints'
CK.mkdir(exist_ok=True)
CFG=yaml.safe_load((ROOT/'configs/phase1h_repeated_fold_seeds.yaml').read_text())
TARGETS=list(m.TARGETS)

def write(df,name):df.to_csv(OUT/name,sep='\t',index=False)
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def gate():
 expected=json.loads((OUT/'pre_run_contract.json').read_text())['sha256']
 baseline=m.read_tsv(m.P1/'phase1d_r2/scientific_freeze_assertions.tsv')
 for row in baseline.itertuples():
  if str(row.path).startswith(('configs/','results/')) and not str(row.path).endswith('.py'):expected[row.path]=row.after_sha256
 common=yaml.safe_load((ROOT/'configs/phase1_common_universe.yaml').read_text())
 expected[common['common_universe_file']]=common['common_universe_sha256']
 expected[common['original_D5_path']]=common['original_D5_sha256']
 expected['results/phase1/model_cache/owkin_phikon-v2/model.safetensors']='261ae680fa699b3b951597fd57aa19c02ef735805acb104b93af69b36d928569'
 ancillary=json.loads((m.P1/'phase1_execution_contract.json').read_text())['ancillary_input_sha256']
 expected.update(ancillary)
 rows=[]
 for p,e in expected.items():
  path=Path(p); path=path if path.is_absolute() else ROOT/path
  a=sha(path);rows.append({'path':p,'expected_sha256':e,'observed_sha256':a,'pass':a==e})
 df=pd.DataFrame(rows);write(df,'input_hash_gate.tsv');assert df['pass'].all(),'FROZEN_INPUT_MISMATCH'
 print('PHASE1H_HASH_GATE_PASS',len(df),flush=True)

def efficient_embeddings():
 raw={};gray={};mac={}
 files=sorted(p for p in m.EMBEDDING_DIR.glob('*.npz') if not p.name.startswith('._'))
 assert len(files)==26
 for p in files:
  with np.load(p,allow_pickle=False) as z:
   assert z['inference_device'].item()=='cpu'
   assert z['model_revision'].item()=='2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd'
   ids=z['region_ids'].astype(str)
   raw.update(zip(ids,z['raw'].astype(np.float64)));gray.update(zip(ids,z['grayscale'].astype(np.float64)))
   for k in z.files:
    if k.startswith('macenko_fold_'):mac.setdefault(k[len('macenko_'):],{}).update(zip(ids,z[k].astype(np.float64)))
 return raw,gray,mac

def metrics(oof):
 y=oof.observed.to_numpy(float);d=oof.donor_id.to_numpy(str)
 a={k:m.regression_metrics(y,oof[k+'_pred'].to_numpy(float)) for k in ['M2','M4']}
 r=m.regression_metrics(oof.heldout_base_residual.to_numpy(float),oof.morphology_predicted_residual.to_numpy(float))
 b=m.donor_cluster_bootstrap(y,{k:oof[k+'_pred'].to_numpy(float) for k in ['M2','M4']},d,m.SEED+100*TARGETS.index(oof.target.iloc[0])+m.ALL_STRATA.index(oof.stratum.iloc[0]))['deltaR2_morphology']
 rb=m.donor_cluster_bootstrap(oof.heldout_base_residual.to_numpy(float),{'residual':oof.morphology_predicted_residual.to_numpy(float)},d,m.SEED+100*TARGETS.index(oof.target.iloc[0])+m.ALL_STRATA.index(oof.stratum.iloc[0])+3000)['residual']
 return dict(M2_R2=a['M2']['R2'],M4_R2=a['M4']['R2'],DeltaR2=a['M4']['R2']-a['M2']['R2'],residual_R2=r['R2'],residual_Pearson=r['Pearson_r'],residual_Spearman=r['Spearman_rho'],delta_CI_low=b['CI_low'],delta_CI_high=b['CI_high'],delta_p=b['bootstrap_two_sided_sign_p'],residual_CI_low=rb['CI_low'],residual_CI_high=rb['CI_high'],residual_p=rb['bootstrap_two_sided_sign_p'],n_donors=len(np.unique(d)),n_regions=len(y))

# Cache exact fit/predict calls for reused training subsets, bounded to avoid memory growth.
# No model algebra, alpha grid, scaling or target transformation is altered.
def refit(frame,features,folds,tag,complete):
 path=CK/(tag+'.tsv.gz');summary=CK/(tag+'.json')
 if summary.exists():return json.loads(summary.read_text())
 y=frame[m.TARGETS[frame.target.iloc[0]]['y']].to_numpy(float);d=frame.donor_id.to_numpy(str)
 oof=frame[['region_id','sample_id','donor_id','target','stratum']].copy();oof['fold']=folds;oof['observed']=y
 cols=list(features) if complete else ['M2','M4']
 cols=[k for k in cols if k!='raw_embedding' and k!='macenko_by_fold']
 for k in cols:oof[k+'_pred']=np.nan
 oof['heldout_base_residual']=np.nan;oof['morphology_predicted_residual']=np.nan
 audit=[]; training=[]
 for f in sorted(np.unique(folds)):
  fp=CK/(tag+'__'+f+'.tsv.gz');ap=CK/(tag+'__'+f+'_alphas.tsv');tp=CK/(tag+'__'+f+'_train_residual.tsv.gz')
  te=folds==f;tr=~te
  assert set(d[te]).isdisjoint(d[tr]) and len(np.unique(d[tr]))>=3
  if fp.exists() and ap.exists() and tp.exists():
   block=pd.read_csv(fp,sep='\t'); assert block.region_id.tolist()==oof.loc[te].region_id.tolist()
   for c in block.columns:
    if c.endswith('_pred') or c in ['heldout_base_residual','morphology_predicted_residual']:oof.loc[te,c]=block[c].to_numpy()
   audit.extend(pd.read_csv(ap,sep='\t').to_dict('records'));continue
  start=len(audit)
  for k in cols:
   x=features[k] if k not in ('macenko','M3_macenko') else features['macenko_by_fold'][f][k]
   pred,alpha=m.fit_predict_one(x[tr],y[tr],d[tr],x[te]);oof.loc[te,k+'_pred']=pred
   audit.append(dict(run=tag,target=frame.target.iloc[0],stratum=frame.stratum.iloc[0],fold=f,model=k,alpha=alpha,n_train_donors=len(np.unique(d[tr])),n_test_donors=len(np.unique(d[te])),donor_overlap=0))
  # Training-side donor cross-fit, using this run's outer partitions only.
  base=np.full(tr.sum(),np.nan);dt=d[tr];ft=folds[tr];xt=features['M2'][tr];yt=y[tr]
  rows=frame.loc[tr,['region_id','donor_id']].reset_index(drop=True)
  for inner in sorted(np.unique(ft)):
   va=ft==inner;itr=~va
   assert set(dt[va]).isdisjoint(dt[itr]) and len(np.unique(dt[itr]))>=3
   pred,alpha=m.fit_predict_one(xt[itr],yt[itr],dt[itr],xt[va]);base[va]=pred
   b=rows.loc[va].copy();b['outer_fold']=f;b['crossfit_fold']=inner;b['observed']=yt[va];b['M2_train_oof_pred']=pred;b['residual']=yt[va]-pred;b['alpha']=alpha;training.append(b)
  assert np.isfinite(base).all()
  test_resid=y[te]-oof.loc[te,'M2_pred'].to_numpy(float)
  pred,alpha=m.fit_predict_one(features['raw_embedding'][tr],yt-base,dt,features['raw_embedding'][te])
  oof.loc[te,'heldout_base_residual']=test_resid;oof.loc[te,'morphology_predicted_residual']=pred
  audit.append(dict(run=tag,target=frame.target.iloc[0],stratum=frame.stratum.iloc[0],fold=f,model='crossfit_residual',alpha=alpha,n_train_donors=len(np.unique(d[tr])),n_test_donors=len(np.unique(d[te])),donor_overlap=0))
  oof.loc[te].to_csv(fp,sep='\t',index=False,compression='gzip');pd.DataFrame(audit[start:]).to_csv(ap,sep='\t',index=False)
  pd.concat(training[-len(np.unique(ft)):]).to_csv(tp,sep='\t',index=False,compression='gzip')
 assert np.isfinite(oof[[k+'_pred' for k in cols]+['heldout_base_residual','morphology_predicted_residual']].to_numpy(float)).all()
 oof.to_csv(path,sep='\t',index=False,compression='gzip');res=metrics(oof);res.update(target=frame.target.iloc[0],stratum=frame.stratum.iloc[0],run=tag)
 summary.write_text(json.dumps(res,indent=2));print('RUN_COMPLETE',tag,flush=True);return res

def ci(a,b):return 'positive' if a>0 else 'negative' if b<0 else 'spans_zero'
def add_fdr(df):
 df=df.copy();df['delta_q']=np.nan;df['residual_q']=np.nan
 for s,g in df.groupby('stratum'):
  df.loc[g.index,'delta_q']=m.bh_adjust(g.delta_p.tolist());df.loc[g.index,'residual_q']=m.bh_adjust(g.residual_p.tolist())
 return df

def support(frame,x,original):
 # Distance definition and unsupported thresholds fixed before this function runs.
 d=frame.donor_id.to_numpy(str);fold=frame.fold.to_numpy(str);c=[q for q in frame if q.startswith('p_')];n=['N_'+q for q in c]
 allrows=[];regionrows=[]
 for f in sorted(np.unique(fold)):
  tr=fold!=f;te=~tr;w=m.donor_weights(d[tr]);sc=StandardScaler().fit(x[tr],sample_weight=w);a=sc.transform(x[tr]);b=sc.transform(x[te])
  lw=LedoitWolf().fit(a);mdtrain=np.sqrt(np.maximum(lw.mahalanobis(a),0));mdtest=np.sqrt(np.maximum(lw.mahalanobis(b),0))
  # Exact other-donor NN for both reference and held-out; training reference avoids self/within-donor neighbors.
  nntrain=np.full(len(a),np.inf);nnte=np.full(len(b),np.inf);td=d[tr]
  for donor in np.unique(td):
   nbr=NearestNeighbors(n_neighbors=1,algorithm='brute',n_jobs=1).fit(a[td==donor])
   eligible=td!=donor;nntrain[eligible]=np.minimum(nntrain[eligible],nbr.kneighbors(a[eligible],return_distance=True)[0][:,0]);nnte=np.minimum(nnte,nbr.kneighbors(b,return_distance=True)[0][:,0])
  ti=np.flatnonzero(te)
  for donor in np.unique(d[te]):
   mask=d[te]==donor;ind=ti[mask];base=dict(target=frame.target.iloc[0],stratum=frame.stratum.iloc[0],donor_id=donor,fold=f,n_regions=len(ind))
   for label,train_dist,test_dist in [('mahalanobis',mdtrain,mdtest),('nearest_neighbor',nntrain,nnte)]:
    v=test_dist[mask];q95,q99=np.quantile(train_dist,[.95,.99]);base.update({label+'_median':np.median(v),label+'_p90':np.quantile(v,.9),label+'_p95':np.quantile(v,.95),label+'_fraction_outside95':np.mean(v>q95),label+'_fraction_outside99':np.mean(v>q99),label+'_train95':q95,label+'_train99':q99})
   for label,columns in [('C',c),('N',n)]:
    v=frame.iloc[ind][columns].to_numpy(float);a0=frame.loc[tr,columns].to_numpy(float);outside=(v<a0.min(0))|(v>a0.max(0))
    base.update({label+'_fraction_features_outside':outside.mean(),label+'_fraction_regions_ge1':np.mean(outside.sum(1)>=1),label+'_fraction_regions_ge2':np.mean(outside.sum(1)>=2)})
   allrows.append(base)
  rr=frame.loc[te,['region_id','donor_id','target','stratum']].copy();rr['mahalanobis']=mdtest;rr['nearest_neighbor']=nnte;regionrows.append(rr)
 return allrows,regionrows

def main():
 gate();m.load_embeddings=efficient_embeddings
 matrix,color,raw,gray,mac=m.load_inputs();ft=m.read_tsv(m.FOLDS);co=m.read_tsv(m.COHORT);donors=sorted(ft.donor_id.astype(str))
 original=pd.read_csv(m.P1/'oof_predictions_primary.tsv.gz',sep='\t');original=original.rename(columns={'analysis_stratum':'stratum'})
 assert original.stratum.nunique()==3
 datasets={};baseline=[];errors=[];supports=[];regsupport=[]
 for t in TARGETS:
  spec=m.TARGETS[t];base=matrix[matrix.full_crop_in_bounds.eq('YES')&matrix[spec['n']].ge(20)&matrix[spec['y']].notna()].copy()
  for s in m.ALL_STRATA:
   frame=base[base.disease.eq('pulmonary_fibrosis')].copy() if s=='PF_ONLY' else base.copy()
   frame=frame.sort_values(['donor_id','sample_id','region_id']).reset_index(drop=True);frame['target']=t;frame['stratum']=s;frame['analysis_stratum']=s
   # Fixed embeddings retained; original Macenko references are conditional frozen scientific inputs.
   features=m.build_model_features(frame,color,raw,gray,mac,'fold_1')
   features['macenko_by_fold']={f:{k:v for k,v in m.build_model_features(frame,color,raw,gray,mac,f).items() if k in ['macenko','M3_macenko']} for f in sorted(ft.fold.unique())}
   datasets[t,s]=(frame,features)
   oof=original[original.target.eq(t)&original.stratum.eq(s)].copy();assert set(oof.region_id)==set(frame.region_id)
   v=metrics(oof);v.update(target=t,stratum=s);baseline.append(v)
   for donor,g in oof.groupby('donor_id'):
    y=g.observed.to_numpy(float);a=y-g.M2_pred.to_numpy(float);b=y-g.M4_pred.to_numpy(float);r=g.heldout_base_residual-g.morphology_predicted_residual
    errors.append(dict(target=t,stratum=s,donor_id=donor,n_regions=len(g),target_mean=y.mean(),target_variance=y.var(),M2_R2=m.r2_score(y,g.M2_pred),M4_R2=m.r2_score(y,g.M4_pred),M2_MAE=np.abs(a).mean(),M4_MAE=np.abs(b).mean(),M2_RMSE=np.sqrt(np.mean(a*a)),M4_RMSE=np.sqrt(np.mean(b*b)),mean_absolute_M2_residual=np.abs(a).mean(),mean_absolute_M4_residual=np.abs(b).mean(),improvement=np.abs(a).mean()-np.abs(b).mean(),M4_minus_M2_error=np.abs(b).mean()-np.abs(a).mean(),paired_SSE_gain=np.sum(a*a-b*b),residual_model_MAE=np.abs(r).mean(),residual_model_RMSE=np.sqrt(np.mean(r*r)),residual_R2=m.r2_score(g.heldout_base_residual,g.morphology_predicted_residual)))
   sr,rr=support(frame,features['M2'],oof);supports.extend(sr);regsupport.extend(rr)
 baseline=add_fdr(pd.DataFrame(baseline));write(baseline,'full_frozen_metrics.tsv')
 err=pd.DataFrame(errors)
 for _,g in err.groupby(['target','stratum']):err.loc[g.index,'LOW_TARGET_VARIANCE']=g.target_variance.le(g.target_variance.quantile(.1))
 write(err,'donor_error_profile.tsv');sp=pd.DataFrame(supports);write(sp,'m2_support_distance.tsv');rs=pd.concat(regsupport);write(rs,'region_support_distance.tsv.gz')
 paired=err.merge(sp,on=['target','stratum','donor_id','n_regions'],validate='one_to_one');paired['unsupported_donor']=paired.mahalanobis_fraction_outside99.ge(.1)
 for _,g in paired.groupby(['target','stratum']):
  bad2=g.M2_MAE.gt(g.M2_MAE.quantile(.75));bad4=g.M4_MAE.gt(g.M4_MAE.quantile(.75));better=g.M4_MAE.lt(g.M2_MAE);much=g.M4_MAE.le(.75*g.M2_MAE)
  case=np.select([bad2&bad4&~much,bad2&much,~bad2&better,~bad2&~better],['CASE1_HARD_DONOR','CASE2_POTENTIAL_RESCUE','CASE3_INCREMENTAL_SIGNAL','CASE4_MORPHOLOGY_HARM'],default='MIXED_BAD_M2_WITHOUT_LARGE_RESCUE');paired.loc[g.index,'paired_case']=case
 write(paired,'m2_m4_paired_leverage.tsv')
 associations=[]
 for (t,s),g in paired.groupby(['target','stratum']):
  rho,p=spearmanr(g.mahalanobis_median,g.M2_MAE);r=rs[rs.target.eq(t)&rs.stratum.eq(s)].merge(original[original.target.eq(t)&original.stratum.eq(s)][['region_id','observed','M2_pred']],on='region_id',validate='one_to_one');rr,rp=spearmanr(r.mahalanobis,np.abs(r.observed-r.M2_pred))
  associations.append(dict(target=t,stratum=s,donor_Spearman=rho,donor_p=p,n_donors=len(g),region_Spearman_exploratory=rr,region_p_exploratory=rp))
 write(pd.DataFrame(associations),'support_error_association.tsv')
 print('DIAGNOSTIC_SUPPORT_COMPLETE',flush=True)
 # Full original metrics are retained; every LOO refits all model/shortcut arms and nested residual pipeline.
 loo=[]
 for di,donor in enumerate(donors):
  rows=[]
  for (t,s),(frame,features) in datasets.items():
   keep=frame.donor_id.ne(donor).to_numpy();tag=f'loo_{donor}__{t}__{s}'
   if keep.all():v=baseline[baseline.target.eq(t)&baseline.stratum.eq(s)].iloc[0].to_dict();v['run']=tag;v['donor_not_in_stratum']=True
   else:
    ff={k:({f:{kk:vv[keep] for kk,vv in fs.items()} for f,fs in v.items()} if k=='macenko_by_fold' else v[keep]) for k,v in features.items()}
    sub=frame.loc[keep].reset_index(drop=True);v=refit(sub,ff,sub.fold.to_numpy(str),tag,True);v['donor_not_in_stratum']=False
   v['donor_id']=donor;rows.append(v)
  rows=add_fdr(pd.DataFrame(rows));loo.append(rows);write(pd.concat(loo),'loo_refit_metrics_progress.tsv');print('LOO_COMPLETE',di+1,donor,flush=True)
 loo=pd.concat(loo,ignore_index=True);infl=loo.merge(baseline,on=['target','stratum'],suffixes=('','_full'),validate='many_to_one')
 for k in ['M2_R2','M4_R2','DeltaR2','residual_R2','residual_Pearson','residual_Spearman']:infl['full_'+k]=infl[k+'_full'];infl['loo_'+k]=infl[k]
 infl['delta_influence']=infl.DeltaR2-infl.DeltaR2_full;infl['residual_influence']=infl.residual_R2-infl.residual_R2_full
 infl['delta_sign_flip']=np.sign(infl.DeltaR2)!=np.sign(infl.DeltaR2_full);infl['residual_sign_flip']=np.sign(infl.residual_R2)!=np.sign(infl.residual_R2_full)
 for k in ['delta','residual']:
  infl[k+'_CI_classification_change']=[ci(a,b)!=ci(c,d) for a,b,c,d in zip(infl[k+'_CI_low'],infl[k+'_CI_high'],infl[k+'_CI_low_full'],infl[k+'_CI_high_full'])]
  infl[k+'_FDR_classification_change']=infl[k+'_q'].lt(.05)!=infl[k+'_q_full'].lt(.05)
 infl['classification_change']=infl[['delta_CI_classification_change','residual_CI_classification_change','delta_FDR_classification_change','residual_FDR_classification_change']].any(axis=1)
 infl['verdict_change']=infl.classification_change|infl.delta_sign_flip|infl.residual_sign_flip
 for _,g in infl.groupby(['target','stratum']):
  extreme=g.verdict_change;mag=g.delta_influence.abs();high=mag.ge(mag.quantile(.9))|(mag/np.maximum(g.DeltaR2_full.abs(),1e-12)).ge(.5)
  infl.loc[g.index,'leverage_category']=np.select([extreme,high,mag.ge(mag.median())],['EXTREME_LEVERAGE','HIGH_LEVERAGE','MODERATE_LEVERAGE'],default='LOW_LEVERAGE')
 write(infl,'leave_one_donor_out_influence.tsv')
 # Clinical/technical annotations use frozen metadata and frozen predictor distributions only.
 annotations=[]
 for donor in sorted(infl[infl.leverage_category.isin(['HIGH_LEVERAGE','EXTREME_LEVERAGE'])].donor_id.unique()):
  g=matrix[matrix.donor_id.eq(donor)];c=co[co.donor_id.eq(donor)&co.primary_inclusion.eq('YES')]
  row=dict(donor_id=donor,disease=';'.join(sorted(g.disease.unique())),TMA=';'.join(map(str,sorted(c.TMA.unique()))),run=';'.join(map(str,sorted(c.run.unique()))),section_count=g.sample_id.nunique(),total_common_regions=len(g),LF_MF=';'.join(sorted(c.affected_status.astype(str).unique())))
  for col in [k for k in g if k.startswith(('p_','N_p_'))]:row[col+'_median']=g[col].median();row[col+'_p10']=g[col].quantile(.1);row[col+'_p90']=g[col].quantile(.9)
  annotations.append(row)
 write(pd.DataFrame(annotations),'high_leverage_frozen_metadata.tsv')
 # Deterministic disease-balanced partitions, written in full before fitting repetitions.
 partitions=[];disease=matrix.groupby('donor_id').disease.first()
 for seed in CFG['seeds']:
  rng=np.random.default_rng(seed);mapping={}
  for label in sorted(disease.unique()):
   ids=np.array(sorted(disease[disease.eq(label)].index));rng.shuffle(ids);offset=int(rng.integers(0,5))
   mapping.update({donor:'fold_'+str(1+(i+offset)%5) for i,donor in enumerate(ids)})
  partitions.extend(dict(seed=seed,donor_id=d,fold=mapping[d],disease=disease[d]) for d in donors)
 part=pd.DataFrame(partitions);write(part,'repeated_partition_assignments.tsv')
 repeats=[];invalid=[]
 for seed in CFG['seeds']:
  mapping=part[part.seed.eq(seed)].set_index('donor_id').fold.to_dict()
  for (t,s),(frame,features) in datasets.items():
   if s not in CFG['repeated_strata']:continue
   folds=frame.donor_id.map(mapping).to_numpy(str)
   valid=all(len(np.unique(frame.donor_id.to_numpy(str)[(folds!=f)&(folds!=inner)]))>=3 for f in np.unique(folds) for inner in np.unique(folds[folds!=f]))
   if not valid:invalid.append(dict(seed=seed,target=t,stratum=s,status='INVALID_MINIMUM_TRAINING_DONORS'));continue
   v=refit(frame,{k:features[k] for k in ['M2','M4','raw_embedding']},folds,f'repeat_{seed}__{t}__{s}',False);v['seed']=seed;repeats.append(v)
  write(pd.DataFrame(repeats),'repeated_fold_metrics.tsv');write(pd.DataFrame(invalid,columns=['seed','target','stratum','status']),'invalid_partitions.tsv');print('SEED_COMPLETE',seed,flush=True)
 gate();(OUT/'analysis_complete.json').write_text(json.dumps({'status':'PHASE1H_ANALYSIS_COMPLETE','loo_donors':len(donors),'repeated_seeds':len(CFG['seeds']),'n_repeated_target_strata':len(repeats),'n_invalid':len(invalid),'scientific_inputs_unchanged':True},indent=2))
 print('PHASE1H_ANALYSIS_PASS',flush=True)

if __name__=='__main__':
 with threadpool_limits(limits=4):main()
