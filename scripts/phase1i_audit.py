#!/usr/bin/env python3
"""Frozen prediction-block inference audit; no fitting and no external data."""
from pathlib import Path
import hashlib,json,datetime
import numpy as np
import pandas as pd
import yaml
from scipy.stats import t,norm,binomtest,rankdata
import phase1e_g_fit_models as m
R=m.ROOT;O=R/'results/phase1i';H=R/'results/phase1h';CFG=R/'configs/phase1i_inference_seeds.yaml'
C=yaml.safe_load(CFG.read_text());TS=C['targets'];EFFECTS=['DeltaR2_morph','Residual_R2'];S='DISEASE_ADJUSTED';N=19

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def read(p):return pd.read_csv(p,sep='\t')
def write(d,n):d.to_csv(O/n,sep='\t',index=False)
def freeze_gate():
 cp=O/'pre_run_contract.json'
 if not cp.exists():
  prior=read(H/'input_hash_gate.tsv');assert prior['pass'].all()
  expected=dict(zip(prior.path,prior.expected_sha256))
  inputs=['configs/phase1i_inference_seeds.yaml','results/phase1h/leave_one_donor_out_influence.tsv','results/phase1h/full_frozen_metrics.tsv','results/phase1h/repeated_fold_metrics.tsv','results/phase1h/repeated_partition_assignments.tsv','results/phase1h/analysis_complete.json','results/phase1h/verification.json','results/phase1h/pre_run_contract.json','results/phase1h/human_verdict_authorization.json','results/phase1h/verdict.json','results/phase1h/phase1h_execution_status.json','reports/PHASE1H_DONOR_FOLD_STABILITY_AUDIT.md','scripts/phase1i_audit.py']
  inputs += [str(p.relative_to(R)) for p in sorted(m.EMBEDDING_DIR.glob('*.npz')) if not p.name.startswith('._')]
  loo=read(H/'leave_one_donor_out_influence.tsv')
  inputs += ['results/phase1h/checkpoints/'+run+'.json' for run in loo.run if (H/'checkpoints'/(run+'.json')).is_file()]
  for p in inputs:
   assert (R/p).is_file(),p
   expected[p]=sha(R/p)
  cp.write_text(json.dumps({'created_before_computation_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':expected,'request':C,'original_phase1h_verdict':'HOLD_INFERENCE_INSTABILITY'},indent=2))
 expected=json.loads(cp.read_text())['sha256'];rows=[]
 for p,v in expected.items():
  path=Path(p);path=path if path.is_absolute() else R/path
  observed=sha(path);rows.append({'path':p,'expected_sha256':v,'observed_sha256':observed,'pass':observed==v})
 d=pd.DataFrame(rows);write(d,'input_hash_gate.tsv');assert d['pass'].all(),'Frozen input mismatch; do not refit automatically'
 print('PHASE1I_HASH_GATE_PASS',len(d),flush=True)

def sign(x):return 'POSITIVE' if x>0 else 'NEGATIVE' if x<0 else 'ZERO'
def classification(est,lo,hi,p=np.nan):
 if np.isfinite(lo) and np.isfinite(hi):return 'POSITIVE_SUPPORT' if lo>0 else 'NEGATIVE_SUPPORT' if hi<0 else 'UNCERTAIN'
 return sign(est)+'_SUPPORT' if np.isfinite(p) and p<.05 and est!=0 else 'UNCERTAIN'
def accel(v):
 u=np.mean(v)-v;s=np.sum(u*u)
 return float(np.sum(u**3)/(6*s**1.5)) if s>0 else np.nan

def block_stats(g):
 rows=[]
 for d,x in g.groupby('donor_id',sort=True):
  y=x.observed.to_numpy(float);a=x.M2_pred.to_numpy(float);b=x.M4_pred.to_numpy(float);r=x.heldout_base_residual.to_numpy(float);pr=x.morphology_predicted_residual.to_numpy(float)
  rows.append(dict(donor_id=d,n_regions=len(x),sum_y=y.sum(),sum_y2=(y*y).sum(),SSE_M2=np.sum((y-a)**2),SSE_M4=np.sum((y-b)**2),MAE_M2=np.mean(abs(y-a)),MAE_M4=np.mean(abs(y-b)),sum_r=r.sum(),sum_r2=np.sum(r*r),SSE_residual=np.sum((r-pr)**2),MAE_zero=np.mean(abs(r)),MAE_residual=np.mean(abs(r-pr))))
 return pd.DataFrame(rows)
def block_effect(counts,b,e):
 n=counts@b.n_regions.to_numpy(float)
 if e=='DeltaR2_morph':sy=counts@b.sum_y.to_numpy(float);sy2=counts@b.sum_y2.to_numpy(float);num=counts@(b.SSE_M2-b.SSE_M4).to_numpy(float)
 else:sy=counts@b.sum_r.to_numpy(float);sy2=counts@b.sum_r2.to_numpy(float);num=sy2-counts@b.SSE_residual.to_numpy(float)
 sst=sy2-sy*sy/n
 if e=='Residual_R2':num=sst-counts@b.SSE_residual.to_numpy(float)
 return np.divide(num,sst,out=np.full(np.shape(sst),np.nan),where=sst>0)
def count_draws(chosen):
 return np.array([np.bincount(v,minlength=N) for v in chosen],dtype=np.int16)
def tail_p(v):
 v=v[np.isfinite(v)];return min(1,2*min((np.count_nonzero(v<=0)+1)/(len(v)+1),(np.count_nonzero(v>=0)+1)/(len(v)+1)))
def w_exact(v):
 v=np.asarray(v);v=v[v!=0]
 if not len(v):return 1.,0.
 ranks=np.rint(rankdata(abs(v),method='average')*2).astype(int);obs=int(ranks[v>0].sum());dp=np.zeros(ranks.sum()+1,dtype=np.int64);dp[0]=1
 for rank in ranks:dp[rank:]+=dp[:-rank].copy()
 p=min(1.,2*min(dp[:obs+1].sum(),dp[obs:].sum())/(2**len(v)))
 return p,float(min(obs,ranks.sum()-obs)/2)
def bca(v,theta,jack):
 finite=v[np.isfinite(v)];a=accel(jack);prop=(np.count_nonzero(finite<theta)+.5*np.count_nonzero(finite==theta))/len(finite);z0=norm.ppf(prop);reason=[]
 invalid=1-len(finite)/len(v)
 if invalid>.05:reason.append('more_than_5pct_invalid')
 if len(np.unique(finite))<2:reason.append('degenerate_bootstrap_distribution')
 if not np.isfinite(a) or not np.isfinite(z0):reason.append('undefined_acceleration_or_bias')
 def endpoints(values):
  prop=(np.count_nonzero(values<theta)+.5*np.count_nonzero(values==theta))/len(values);z=norm.ppf(prop)
  den=1-a*(z+norm.ppf([.025,.975]));probs=norm.cdf(z+(z+norm.ppf([.025,.975]))/den)
  return np.quantile(values,probs),probs,den
 lo=hi=np.nan;probs=[np.nan,np.nan];spread=np.nan
 if not reason:
  endpoints_full,probs,den=endpoints(finite);lo,hi=endpoints_full
  if np.any(den<=0) or np.any(~np.isfinite(probs)):reason.append('unstable_adjustment_denominator')
  if min(probs[0],1-probs[1])*len(finite)<10:reason.append('adjusted_tail_under_10_replicates')
  chunkends=[]
  for chunk in np.array_split(v,4):
   chunk=chunk[np.isfinite(chunk)]
   try:ends,_,_=endpoints(chunk);chunkends.append(ends)
   except (ValueError,FloatingPointError):reason.append('chunk_endpoint_undefined')
  if len(chunkends)==4:
   spread=np.max(abs(np.asarray(chunkends)-endpoints_full))
   if hi<=lo or spread>.25*(hi-lo):reason.append('endpoint_chunk_deviation_over_25pct_width')
 return dict(BCa_valid=not bool(reason),BCa_status='BCA_VALID' if not reason else 'BCA_UNRELIABLE',BCa_reason=';'.join(reason),CI_low=float(lo),CI_high=float(hi),acceleration=float(a),bias_correction_z0=float(z0),adjusted_quantile_low=float(probs[0]),adjusted_quantile_high=float(probs[1]),chunk_endpoint_max_deviation=float(spread),invalid_fraction=invalid)

def main():
 freeze_gate()
 pred=read(m.P1/'oof_predictions_primary.tsv.gz');loo=read(H/'leave_one_donor_out_influence.tsv');rep=read(H/'repeated_fold_metrics.tsv');pm=read(m.P1/'primary_model_performance.tsv');rm=read(m.P1/'residual_prediction_performance.tsv')
 common=read(m.P1/'phase1d_r2/final_common_phase1_region_universe.tsv');assert len(common)==12144
 admitted=read(m.COHORT);assert admitted.primary_inclusion.eq('YES').sum()==26
 fold=read(m.FOLDS);assert len(fold)==19 and fold.fold.nunique()==5
 pseudo=[];influence=[];js=[];blocks=[];lossrows=[];signrows=[];equal=[];agree=[];bootrows=[];bcadiag=[];directions=[]
 for ti,target in enumerate(TS):
  g=pred[pred.target.eq(target)&pred.stratum.eq(S)].copy();assert g.donor_id.nunique()==19 and not g.region_id.duplicated().any() and set(g.region_id)<=set(common.region_id)
  assert g.fold.eq(g.donor_id.map(fold.set_index('donor_id').fold)).all()
  b=block_stats(g);b['target']=target;blocks.append(b);l=loo[loo.target.eq(target)&loo.stratum.eq(S)].sort_values('donor_id');assert len(l)==19 and l.n_donors.eq(18).all() and l.donor_id.tolist()==b.donor_id.tolist()
  for rr in l.itertuples():
   cached=json.loads((H/'checkpoints'/f'{rr.run}.json').read_text())
   assert np.isclose(rr.DeltaR2,cached['DeltaR2']) and np.isclose(rr.residual_R2,cached['residual_R2'])
  rng=np.random.default_rng(C['balanced_seeds'][target]);bag=np.repeat(np.arange(N),C['balanced_replicates']);rng.shuffle(bag);draws=bag.reshape(C['balanced_replicates'],N);assert np.bincount(bag).tolist()==[10000]*19;balanced_counts=count_draws(draws)
  delete_counts=np.ones((N,N));np.fill_diagonal(delete_counts,0)
  for ei,e in enumerate(EFFECTS):
   theta=float(block_effect(np.ones((1,N)),b,e)[0]);r0=(pm[pm.target.eq(target)&pm.stratum.eq(S)] if ei==0 else rm[rm.target.eq(target)&rm.stratum.eq(S)]).iloc[0]
   reference_est=float(r0.deltaR2_morph if ei==0 else r0.residual_R2);assert np.isclose(theta,reference_est,atol=1e-12)
   reflo=float(r0.deltaR2_morph_CI_low if ei==0 else r0.CI_low);refhi=float(r0.deltaR2_morph_CI_high if ei==0 else r0.CI_high);refp=float(r0.deltaR2_morph_p if ei==0 else r0.p);refq=float(r0.deltaR2_morph_FDR if ei==0 else r0.FDR)
   refseed=m.SEED+100*ti+m.ALL_STRATA.index(S)+(3000 if ei else 0);refdraws=np.random.default_rng(refseed).integers(0,N,size=(2000,N));refvalues=block_effect(count_draws(refdraws),b,e);refinterval=np.quantile(refvalues[np.isfinite(refvalues)],[.025,.975]);assert np.allclose(refinterval,[reflo,refhi],rtol=1e-8,atol=1e-8) and np.isclose(tail_p(refvalues),refp)
   lm=l.DeltaR2.to_numpy() if ei==0 else l.residual_R2.to_numpy();pv=N*theta-(N-1)*lm;jmean=pv.mean();bias=(N-1)*(lm.mean()-theta);se=pv.std(ddof=1)/np.sqrt(N);jklo,jkhi=jmean+np.array([-1,1])*t.ppf(.975,18)*se;jp=2*t.sf(abs(jmean/se),18) if se>0 else np.nan
   assert np.isclose(jmean,theta-bias)
   for idx,d in enumerate(b.donor_id):
    others=np.delete(pv,idx);sd=others.std(ddof=1);z=(pv[idx]-jmean)/pv.std(ddof=1) if pv.std(ddof=1)>0 else np.nan;student=(pv[idx]-others.mean())/(sd*np.sqrt(N/(N-1))) if sd>0 else np.nan
    pseudo.append(dict(target=target,effect=e,donor_id=d,n_donors=N,theta_full=theta,theta_minus_i=lm[idx],pseudovalue=pv[idx],jackknife_mean=jmean,jackknife_bias=bias,jackknife_SE=se,df=18))
    influence.append(dict(target=target,effect=e,donor_id=d,delete1_influence=lm[idx]-theta,centered_pseudovalue=pv[idx]-jmean,standardized_influence=z,studentized_influence=student))
   absin=abs(lm-theta);order=np.argsort(-absin);total=absin.sum()
   js.append(dict(target=target,effect=e,original_estimate=theta,JACKKNIFE_BIAS_CORRECTED_ESTIMATE=jmean,jackknife_bias=bias,jackknife_SE=se,CI_low=jklo,CI_high=jkhi,p_value=jp,df=18,max_absolute_standardized_influence=np.max(abs((pv-jmean)/pv.std(ddof=1))),max_absolute_studentized_influence=max(abs(x['studentized_influence']) for x in influence if x['target']==target and x['effect']==e),median_absolute_influence=np.median(absin),p90_absolute_influence=np.quantile(absin,.9),top_donor=b.donor_id.iloc[order[0]],top2_donors=';'.join(b.donor_id.iloc[order[:2]]),top2_combined_abs_influence=absin[order[:2]].sum(),top2_fraction_abs_influence=absin[order[:2]].sum()/total))
   base=dict(target=target,effect=e,loss='PRIMARY_R2',original_estimate=theta,uncertainty='OOF conditional donor-block sampling; models fixed',n_donors=19)
   agree.append(dict(base,method='REFERENCE_BOOTSTRAP',estimate=theta,CI_low=reflo,CI_high=refhi,p_value=refp,FDR=refq,p_value_status='Original bootstrap sign-tail approximation retained verbatim',valid=True))
   agree.append(dict(base,method='JACKKNIFE_NORMAL_CI',estimate=jmean,CI_low=jklo,CI_high=jkhi,p_value=jp,FDR=np.nan,p_value_status='Approximate t_18 pivot for pseudovalue mean; refitted delete1 sensitivity',valid=True,uncertainty='Refitted delete1 donor sensitivity with t_18 approximation'))
   values=block_effect(balanced_counts,b,e);finite=values[np.isfinite(values)];blo,bhi=np.quantile(finite,[.025,.975]);agree.append(dict(base,method='BALANCED_DONOR_BOOTSTRAP',estimate=theta,CI_low=blo,CI_high=bhi,p_value=tail_p(values),FDR=np.nan,p_value_status='Bootstrap sign-tail approximation; not null-centered test',valid=len(finite)/len(values)>=.95))
   blockjack=block_effect(delete_counts,b,e);bd=bca(values,theta,blockjack);bd.update(target=target,effect=e,refit_jackknife_acceleration=accel(lm),acceleration_source='Frozen OOF-block delete1, matched to conditional bootstrap');bcadiag.append(bd)
   agree.append(dict(base,method='BCA_INFERENCE_SENSITIVITY',estimate=theta,CI_low=bd['CI_low'] if bd['BCa_valid'] else np.nan,CI_high=bd['CI_high'] if bd['BCa_valid'] else np.nan,p_value=np.nan,FDR=np.nan,p_value_status='No calibrated BCa p-value constructed',valid=bd['BCa_valid']))
   for idx,v in enumerate(values):bootrows.append(dict(target=target,effect=e,method='BALANCED_DONOR_BOOTSTRAP',replicate=idx+1,estimate=v))
   # Paired donor loss improvements are sensitivity estimands, never primary R2 replacements.
   gains={'SSE_GAIN':(b.SSE_M2-b.SSE_M4).to_numpy() if ei==0 else (b.sum_r2-b.SSE_residual).to_numpy(),'MAE_GAIN':(b.MAE_M2-b.MAE_M4).to_numpy() if ei==0 else (b.MAE_zero-b.MAE_residual).to_numpy()}
   for loss,v in gains.items():
    positives=int(np.count_nonzero(v>0));negatives=int(np.count_nonzero(v<0));ties=int(np.count_nonzero(v==0));sp=binomtest(positives,positives+negatives,.5,alternative='two-sided').pvalue if positives+negatives else 1.;wp,stat=w_exact(v);median=float(np.median(v));mean=float(np.mean(v))
    for d,nr,value in zip(b.donor_id,b.n_regions,v):lossrows.append(dict(target=target,effect=e,loss=loss,donor_id=d,n_regions=nr,paired_gain=value,improved=value>0,independent_weight=1))
    signrows.append(dict(target=target,effect=e,loss=loss,analysis_role='SENSITIVITY',n_donors=19,donors_improved=positives,donors_worsened=negatives,ties=ties,exact_sign_p=sp,exact_Wilcoxon_p=wp,Wilcoxon_statistic=stat,mean=mean,median=median))
    equal.append(dict(target=target,effect=e,loss=loss,n_donors=19,independent_weight_per_donor=1,mean=mean,median=median,unit='sum squared error units per donor' if loss=='SSE_GAIN' else 'mean absolute error units',summary_selected_by_significance=False))
    for method,pval in [('EXACT_SIGN_TEST',sp),('WILCOXON_SENSITIVITY',wp)]:agree.append(dict(target=target,effect=e,loss=loss,original_estimate=theta,uncertainty='Donor-equal paired-loss sensitivity; different quantity from R2',n_donors=19,method=method,estimate=median,CI_low=np.nan,CI_high=np.nan,p_value=pval,FDR=np.nan,p_value_status='Exact two-sided sign test' if method=='EXACT_SIGN_TEST' else 'Exact conditional signed-rank randomization; symmetry assumption',valid=True))
   part=rep[rep.target.eq(target)&rep.stratum.eq('POOLED')];assert len(part)==50
   partvals=part.DeltaR2.to_numpy() if ei==0 else part.residual_R2.to_numpy();pmid=np.median(partvals);sgn=np.sign(theta);sse=gains['SSE_GAIN'];mae=gains['MAE_GAIN'];major=int(np.sign(theta));loofrac=np.mean(np.sign(lm)==major);partfrac=np.mean(np.sign(partvals)==major)
   summarydirs=[sign(theta),sign(jmean),sign(np.median(lm)),sign(pmid),sign(np.mean(sse)),sign(np.median(sse))];consistent=len(set(summarydirs))==1 and summarydirs[0]!='ZERO' and loofrac>=.8 and partfrac>=.8
   directions.append(dict(target=target,effect=e,full_estimate=theta,full_direction=sign(theta),jackknife_mean=jmean,jackknife_direction=sign(jmean),LOO_median=np.median(lm),LOO_positive=int(np.sum(lm>0)),LOO_negative=int(np.sum(lm<0)),LOO_zero=int(np.sum(lm==0)),LOO_same_full_fraction=loofrac,partition_median=pmid,partition_median_direction=sign(pmid),partition_same_full_fraction=partfrac,partition_p5=np.quantile(partvals,.05),partition_p95=np.quantile(partvals,.95),partition_distribution_role='PARTITION_VARIABILITY_NOT_CI',donor_equal_SSE_mean=np.mean(sse),donor_equal_SSE_median=np.median(sse),donor_equal_MAE_mean=np.mean(mae),donor_equal_MAE_median=np.median(mae),donor_equal_SSE_mean_direction=sign(np.mean(sse)),donor_equal_SSE_median_direction=sign(np.median(sse)),donor_equal_MAE_mean_direction=sign(np.mean(mae)),donor_equal_MAE_median_direction=sign(np.median(mae)),direction_consensus='DIRECTION_CONSISTENT' if consistent else 'DIRECTION_MIXED'))
  print('INFERENCE_COMPLETE',target,flush=True)
 agreement=pd.DataFrame(agree)
 for (method,e,loss),g in agreement.groupby(['method','effect','loss']):
  if method=='REFERENCE_BOOTSTRAP':continue
  ix=g[g.valid&g.p_value.notna()].index
  if len(ix):agreement.loc[ix,'FDR']=m.bh_adjust(agreement.loc[ix,'p_value'].tolist())
 agreement['direction']=agreement.estimate.map(sign);agreement['zero_excluded']=agreement.valid&((agreement.CI_low>0)|(agreement.CI_high<0))
 agreement['classification']=[classification(x.estimate,x.CI_low,x.CI_high,x.p_value) if x.valid else 'UNCERTAIN' for x in agreement.itertuples()]
 agreement['FDR_classification']=[(sign(x.estimate)+'_SUPPORT') if x.valid and np.isfinite(x.FDR) and x.FDR<.05 else 'UNCERTAIN' if np.isfinite(x.FDR) else 'NOT_DEFINED' for x in agreement.itertuples()]
 agreement['CI_width']=agreement.CI_high-agreement.CI_low;agreement['relative_precision_original']=np.where(agreement.original_estimate.abs()>=C['precision_near_zero'],agreement.CI_width/agreement.original_estimate.abs(),np.nan)
 for n,df in [('jackknife_pseudovalues.tsv',pd.DataFrame(pseudo)),('studentized_donor_influence.tsv',pd.DataFrame(influence)),('jackknife_summary.tsv',pd.DataFrame(js)),('donor_prediction_blocks.tsv',pd.concat(blocks)),('donor_paired_losses.tsv',pd.DataFrame(lossrows)),('donor_sign_wilcoxon_results.tsv',pd.DataFrame(signrows)),('donor_equal_effect_summary.tsv',pd.DataFrame(equal)),('bca_feasibility.tsv',pd.DataFrame(bcadiag)),('inference_agreement_matrix.tsv',agreement),('effect_direction_consensus.tsv',pd.DataFrame(directions))]:write(df,n)
 write(pd.DataFrame(bootrows),'balanced_bootstrap_distribution.tsv.gz')
 consensus=[]
 for d in pd.DataFrame(directions).itertuples():
  a=agreement[agreement.target.eq(d.target)&agreement.effect.eq(d.effect)&agreement.valid];primary=a[a.loss.eq('PRIMARY_R2')];signcheck=a[a.method.eq('EXACT_SIGN_TEST')&a.loss.eq('SSE_GAIN')]
  switches=len(set(pd.concat([primary,signcheck]).classification))>1
  fdrset=set(primary[primary.FDR.notna()].FDR_classification);fdrswitch=len(fdrset)>1
  pointdirs=set(primary.direction)|{d.donor_equal_SSE_mean_direction,d.donor_equal_SSE_median_direction};pointdirs.discard('ZERO');directionconflict=len(pointdirs)>1
  supported=set(primary.classification)-{'UNCERTAIN'};opposed={'POSITIVE_SUPPORT','NEGATIVE_SUPPORT'}<=supported
  ic='INFERENCE_DIRECTION_CONFLICT' if directionconflict else 'INFERENCE_METHOD_SENSITIVE' if switches or fdrswitch or d.direction_consensus=='DIRECTION_MIXED' else 'INFERENCE_ROBUST'
  consensus.append(dict(target=d.target,effect=d.effect,direction_consensus=d.direction_consensus,inference_consensus=ic,support_uncertain_switching=switches,FDR_switching=fdrswitch,opposing_supported_CI_methods=opposed,method_point_direction_conflict=len(pointdirs)>1,primary_supported_CI_methods=int(primary.zero_excluded.sum()),valid_BCa=bool(primary.method.eq('BCA_INFERENCE_SENSITIVITY').any())))
 con=pd.DataFrame(consensus);write(con,'effect_inference_consensus.tsv')
 targetrows=[]
 for target in TS:
  dx=pd.DataFrame(directions).query('target == @target').set_index('effect');cx=con.query('target == @target').set_index('effect');aa=agreement[agreement.target.eq(target)&agreement.loss.eq('PRIMARY_R2')&agreement.valid];delta=dx.loc[EFFECTS[0]];resid=dx.loc[EFFECTS[1]]
  null=(aa.CI_low.ge(-.05)&aa.CI_high.le(.05)).all() and aa.estimate.abs().le(.01).all() and not cx.method_point_direction_conflict.any()
  if null:tc='ROBUST_NULL'
  elif cx.direction_consensus.eq('DIRECTION_MIXED').any():tc='DIRECTION_UNSTABLE'
  else:
   direction=delta.full_direction;robust=cx.loc[EFFECTS[0],'inference_consensus']=='INFERENCE_ROBUST' and cx.loc[EFFECTS[0],'primary_supported_CI_methods']>=2
   compatible=resid.full_direction=='POSITIVE' and cx.loc[EFFECTS[1],'inference_consensus']=='INFERENCE_ROBUST' and cx.loc[EFFECTS[1],'primary_supported_CI_methods']>=2
   tc=('ROBUST_POSITIVE' if robust and compatible else 'POSITIVE_BUT_INFERENCE_SENSITIVE') if direction=='POSITIVE' else ('ROBUST_NEGATIVE' if robust else 'NEGATIVE_BUT_INFERENCE_SENSITIVE')
  ll=loo[loo.target.eq(target)&loo.stratum.eq(S)];depend=int(ll.qualitative_interpretation_change.sum());eligible=tc=='ROBUST_POSITIVE' and not (0<depend<=2)
  targetrows.append(dict(target=target,target_classification=tc,delta_direction_consensus=delta.direction_consensus,residual_direction_consensus=resid.direction_consensus,delta_inference=cx.loc[EFFECTS[0],'inference_consensus'],residual_inference=cx.loc[EFFECTS[1],'inference_consensus'],interpretation_flip_donors=depend,external_positive_eligible=eligible))
 tc=pd.DataFrame(targetrows);write(tc,'target_inference_classification.tsv')
 rel=[]
 for target in TS:
  for e in EFFECTS:
   pp=agreement[agreement.target.eq(target)&agreement.effect.eq(e)&agreement.method.isin(['REFERENCE_BOOTSTRAP','JACKKNIFE_NORMAL_CI'])]
   rel.append(bool(pp.relative_precision_original.ge(2).any() or ((pp.original_estimate.abs()<.01)&pp.CI_width.gt(.1)).any()))
 precision='HIGH' if sum(rel)>=2 else 'MODERATE' if sum(rel)==1 else 'LOW';positive=tc.external_positive_eligible.sum()
 if positive:
  verdict='INFERENCE_ROBUST_GO_EXTERNAL_VALIDATION' if tc.target_classification.isin(['ROBUST_POSITIVE','ROBUST_NEGATIVE','ROBUST_NULL']).all() else 'SELECTIVE_INFERENCE_ROBUSTNESS'
 elif con.opposing_supported_CI_methods.any() and con.method_point_direction_conflict.any() and con.direction_consensus.eq('DIRECTION_MIXED').any():verdict='NO_GO_INFERENCE_IDENTIFIABILITY'
 elif precision=='HIGH':verdict='HOLD_FINITE_DONOR_PRECISION'
 else:verdict='HOLD_INFERENCE_METHOD_SENSITIVITY'
 result={'phase1i_verdict':verdict,'external_validation_authorized':'YES' if positive else 'NO','finite_donor_precision_limitation':precision,'primary_point_estimates_changed':False,'original_models_refit':False,'optional_precision_curve':'NOT_RUN','BCa_all_valid':bool(pd.DataFrame(bcadiag).BCa_valid.all()),'n_donors':19,'n_targets':3,'methods_with_direction_conflict':con[con.method_point_direction_conflict].loc[:,['target','effect']].to_dict('records'),'effects_with_support_uncertain_switching':con[con.support_uncertain_switching|con.FDR_switching].loc[:,['target','effect']].to_dict('records'),'recommended_next_phase':'Separate pre-registered independent-donor precision/transportability design, subject to human authorization. No external analysis performed.'}
 (O/'verdict.json').write_text(json.dumps(result,indent=2))
 freeze_gate();(O/'analysis_complete.json').write_text(json.dumps({'status':'PHASE1I_ANALYSIS_COMPLETE','completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_donors':19,'balanced_replicates':10000,'primary_stratum':S,'hash_gate_pass':True,'original_models_refit':False},indent=2));print('PHASE1I_ANALYSIS_PASS',json.dumps(result),flush=True)
if __name__=='__main__':main()
