#!/usr/bin/env python3
"""Prospective augmentation precision projections from frozen donor blocks only."""
from pathlib import Path
import json,hashlib,datetime
import numpy as np
import pandas as pd
import yaml
from scipy.stats import t
R=Path(__file__).resolve().parents[1]; O=R/'results/phase1j'
C=yaml.safe_load((R/'configs/phase1j_precision_plan.yaml').read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def gate():
 p=O/'pre_run_contract.json'
 if not p.exists():
  old=json.loads((R/'results/phase1i/pre_run_contract.json').read_text())['sha256']
  for root in ['configs','results/phase1i']:
   for f in sorted((R/root).rglob('*')):
    if f.is_file() and not f.name.startswith('._') and '.mplconfig' not in str(f):old[str(f.relative_to(R))]=sha(f)
  for s in ['reports/PHASE1I_FINITE_DONOR_INFERENCE_ROBUSTNESS.md','reports/PHASE1H_DONOR_FOLD_STABILITY_AUDIT.md','scripts/phase1j_precision.py']:old[s]=sha(R/s)
  p.write_text(json.dumps({'created_before_projection_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':old,'planning_contract':C},indent=2))
 rows=[]
 for s,h in json.loads(p.read_text())['sha256'].items():
  f=Path(s);f=f if f.is_absolute() else R/f
  v=sha(f);rows.append(dict(path=s,expected_sha256=h,observed_sha256=v,pass_gate=v==h))
 d=pd.DataFrame(rows);d.to_csv(O/'input_hash_gate.tsv',sep='\t',index=False)
 assert d.pass_gate.all(),'Frozen inputs changed'
 print('PHASE1J_HASH_GATE_PASS',len(d),flush=True)
def wilson(k,n):
 z=1.959963984540054;p=k/n;den=1+z*z/n
 return (p+z*z/(2*n)-z*np.sqrt(p*(1-p)/n+z*z/(4*n*n)))/den
def effect(counts,b,e):
 # Whole donor sums reconstruct the original region-weighted R² algebra.
 N=counts@b[:,0]
 if e=='DeltaR2_morph':sy=counts@b[:,1];sy2=counts@b[:,2];num=counts@(b[:,3]-b[:,4])
 else:sy=counts@b[:,5];sy2=counts@b[:,6];num=None
 sst=sy2-sy*sy/N
 if e=='Residual_R2':num=sst-counts@b[:,7]
 return np.divide(num,sst,out=np.full(sst.shape,np.nan),where=sst>0)
def main():
 O.mkdir(parents=True,exist_ok=True);gate()
 rng=np.random.default_rng(C['seed']);T=C['targets'];E=C['effects'];K=C['outer_projections'];B=C['inner_oof_bootstraps']
 df=pd.read_csv(R/'results/phase1i/donor_prediction_blocks.tsv',sep='\t');pv=pd.read_csv(R/'results/phase1i/jackknife_pseudovalues.tsv',sep='\t');agree=pd.read_csv(R/'results/phase1i/inference_agreement_matrix.tsv',sep='\t');jk=pd.read_csv(R/'results/phase1i/jackknife_summary.tsv',sep='\t')
 donors=sorted(df.donor_id.unique());assert len(donors)==19
 cols=['n_regions','sum_y','sum_y2','SSE_M2','SSE_M4','sum_r','sum_r2','SSE_residual']
 blocks={a:df[df.target==a].set_index('donor_id').reindex(donors)[cols].to_numpy(float) for a in T}
 pseudos={};baselines={};current=[]
 basecounts=rng.multinomial(19,np.full(19,1/19),size=C['baseline_inner_oof_bootstraps'])
 for a in T:
  for e in E:
   pp=pv[(pv.target==a)&(pv.effect==e)].set_index('donor_id').reindex(donors).pseudovalue.to_numpy(float);pseudos[a,e]=pp
   row=agree[(agree.target==a)&(agree.effect==e)&(agree.method=='REFERENCE_BOOTSTRAP')].iloc[0];j=jk[(jk.target==a)&(jk.effect==e)].iloc[0]
   bca=agree[(agree.target==a)&(agree.effect==e)&(agree.method=='BCA_INFERENCE_SENSITIVITY')].iloc[0]
   vals=effect(basecounts,blocks[a],e);raw=np.diff(np.nanquantile(vals,[.025,.975]))[0]
   baselines[a,e]=dict(reference=row.CI_width,jackknife=j.CI_high-j.CI_low,calibration=row.CI_width/raw,raw=raw)
   assert np.isclose(2*t.ppf(.975,18)*pp.std(ddof=1)/np.sqrt(19),baselines[a,e]['jackknife'])
   current.append(dict(target=a,effect=e,original_estimate=row.original_estimate,reference_CI_width=row.CI_width,jackknife_CI_width=j.CI_high-j.CI_low,BCa_CI_width=bca.CI_width if bca.valid else np.nan,pseudovalue_SD=pp.std(ddof=1),jackknife_SE=j.jackknife_SE,top2_influence_fraction=j.top2_fraction_abs_influence,baseline_simulated_OOF_CI_width=raw,calibration=baselines[a,e]['calibration']))
 pd.DataFrame(current).to_csv(O/'current_precision.tsv',sep='\t',index=False)
 # Pair all targets/effects and nested donor-count additions through common draws.
 additions=rng.integers(0,19,size=(K,40));outer=np.ones((K,19),dtype=np.int16);rows=[];distributions=[]
 for n in range(19,60):
  if n>19:np.add.at(outer,(np.arange(K),additions[:,n-20]),1)
  assert np.all(outer.sum(axis=1)==n) and np.all(outer>=1)
  widths={(a,e):[] for a in T for e in E}
  invalid={(a,e):[] for a in T for e in E}
  for k in range(K if n>19 else 1):
   c=rng.multinomial(n,outer[k].astype(float)/n,size=B)
   for a in T:
    for e in E:
     v=effect(c,blocks[a],e);bad=1-np.isfinite(v).mean();w=np.diff(np.nanquantile(v,[.025,.975]))[0] if bad<=.05 else np.nan
     widths[a,e].append(w);invalid[a,e].append(bad)
  for a in T:
   for e in E:
    base=baselines[a,e];p=pseudos[a,e]
    mean=outer@p/n;variance=(outer@(p*p)-n*mean*mean)/(n-1)
    jw=2*t.ppf(.975,n-1)*np.sqrt(np.maximum(variance,0)/n)
    ow=np.asarray(widths[a,e])*base['calibration']
    if n==19:ow=np.full(K,base['reference']) # frozen current, not a stochastic redesign
    for method,w,key in [('OOF_BLOCK_EMPIRICAL',ow,'reference'),('REFIT_PSEUDOVALUE_PROXY',jw,'jackknife')]:
     relative=w/base[key];d=dict(target=a,effect=e,method=method,total_donors=n,added_donors=n-19,projected_CI_width=np.nanmedian(w),width_p10=np.nanquantile(w,.1),width_p90=np.nanquantile(w,.9),relative_width=np.nanmedian(relative),current_CI_width=base[key],n_projections=K,unreliable_projection_fraction=float(np.mean(~np.isfinite(w))),asymptotic_relative_width=np.sqrt(19/n),max_inner_invalid_fraction=max(invalid[a,e]))
     for ratio,label in [(.8,'80'),(.67,'67'),(.5,'50')]:
      hits=np.count_nonzero(np.isfinite(relative)&(relative<=ratio));d['probability_width_le_'+label+'pct_current']=hits/K;d['probability_'+label+'pct_Wilson95_lower']=wilson(hits,K)
     rows.append(d)
     distributions.extend(dict(target=a,effect=e,method=method,total_donors=n,projection=i+1,CI_width=float(v)) for i,v in enumerate(w))
  if n in C['total_donor_grid']:print('PRECISION_PROJECTION_COMPLETE n='+str(n),flush=True)
 out=pd.DataFrame(rows);out.to_csv(O/'precision_projection.tsv',sep='\t',index=False);pd.DataFrame(distributions).to_csv(O/'planning_width_distributions.tsv.gz',sep='\t',index=False,compression='gzip')
 minima=[]
 for a in T:
  for label,effect_filter,method_filter in [('CONSERVATIVE_BOTH_EFFECTS_BOTH_METHODS',E,['OOF_BLOCK_EMPIRICAL','REFIT_PSEUDOVALUE_PROXY'])]+[(m+'_'+e,[e],[m]) for e in E for m in ['OOF_BLOCK_EMPIRICAL','REFIT_PSEUDOVALUE_PROXY']]:
   d=out[(out.target==a)&out.effect.isin(effect_filter)&out.method.isin(method_filter)]
   good=d.groupby('total_donors').probability_67pct_Wilson95_lower.min()>=.8
   ns=[n for n in good.index if n>19 and good.loc[n:].all()]
   minima.append(dict(target=a,criterion=label,minimum_total_donors=min(ns) if ns else np.nan,MINIMUM_ADDITIONAL_DONORS=min(ns)-19 if ns else np.nan,status='FEASIBLE_UNDER_EMPIRICAL_PROJECTION' if ns else 'PRECISION_RESCUE_NOT_FEASIBLE_AT_PRACTICAL_N'))
 pd.DataFrame(minima).to_csv(O/'minimum_information_threshold.tsv',sep='\t',index=False)
 gate();(O/'precision_complete.json').write_text(json.dumps({'status':'PHASE1J_PRECISION_PLANNING_COMPLETE','no_model_fits':True,'seed':C['seed'],'outer':K,'inner':B,'n_grid':list(range(19,60))},indent=2));print('PHASE1J_PRECISION_COMPLETE',flush=True)
if __name__=='__main__':main()
