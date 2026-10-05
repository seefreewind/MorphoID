#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Evaluate original frozen Phase1 decision logic on complete LOO refits."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import phase1e_g_fit_models as m
import phase1_report as policy
O=m.ROOT/'results/phase1h';CK=O/'checkpoints'

def prepare(blocks):
 metrics=[];boots=[];residual=[];shortcuts=[];donors=[]
 for (t,s),g in blocks.groupby(['target','stratum'],sort=False):
  y=g.observed.to_numpy(float);d=g.donor_id.to_numpy(str);pred={k:g[k+'_pred'].to_numpy(float) for k in ['M0','M1','M2','M3','M4','color','coordinate','grayscale','macenko','M3_grayscale','M3_macenko']}
  scores={k:m.r2_score(y,v) for k,v in pred.items()};seed=m.SEED+100*list(m.TARGETS).index(t)+m.ALL_STRATA.index(s)
  b=m.donor_cluster_bootstrap(y,pred,d,seed)
  metrics.append(dict(target=t,analysis_stratum=s,M4_R2=scores['M4'],deltaR2_morphology=scores['M4']-scores['M2'],deltaR2_composition=scores['M1']-scores['M0'],deltaR2_neighborhood=scores['M2']-scores['M1'],p=b['deltaR2_morphology']['bootstrap_two_sided_sign_p']))
  for key in ['deltaR2_composition','deltaR2_neighborhood','deltaR2_morphology']:boots.append(dict(target=t,analysis_stratum=s,metric=key,model='paired_increment',**b[key]))
  rb=m.donor_cluster_bootstrap(g.heldout_base_residual.to_numpy(float),{'residual':g.morphology_predicted_residual.to_numpy(float)},d,seed+3000)['residual']
  residual.append(dict(target=t,analysis_stratum=s,residual_R2=m.r2_score(g.heldout_base_residual,g.morphology_predicted_residual),CI_low=rb['CI_low'],CI_high=rb['CI_high'],p=rb['bootstrap_two_sided_sign_p']))
  shortcuts.append(dict(target=t,analysis_stratum=s,color_only=scores['color'],coordinates_only=scores['coordinate'],raw=scores['M4'],deltaR2_macenko=scores['macenko']-scores['M2'],deltaR2_grayscale=scores['grayscale']-scores['M2']))
  for donor,v in g.groupby('donor_id'):
   y0=v.observed.to_numpy();a=y0-v.M2_pred.to_numpy();b0=y0-v.M4_pred.to_numpy()
   donors.append(dict(target=t,analysis_stratum=s,donor_id=donor,delta_MAE_M2_minus_M4=np.mean(abs(a))-np.mean(abs(b0)),positive_squared_error_gain=max(0,np.sum(a*a-b0*b0))))
 df=pd.DataFrame(metrics);rf=pd.DataFrame(residual)
 for s in m.ALL_STRATA:
  ix=df.analysis_stratum.eq(s);df.loc[ix,'FDR']=m.bh_adjust(df.loc[ix,'p'].tolist());ir=rf.analysis_stratum.eq(s);rf.loc[ir,'FDR']=m.bh_adjust(rf.loc[ir,'p'].tolist())
 hetero,ht=policy.classify_heterogeneity(pd.DataFrame(donors))
 verdict,robust,comp,shortcut=policy.choose_verdict(df,pd.DataFrame(boots),rf,pd.DataFrame(shortcuts),hetero,ht)
 return dict(verdict=verdict,robust_targets=robust,composition_count=comp,shortcut_count=shortcut,heterogeneity=hetero)

def main():
 assert (O/'analysis_complete.json').exists()
 old=pd.read_csv(m.P1/'oof_predictions_primary.tsv.gz',sep='\t').rename(columns={'analysis_stratum':'stratum'})
 full=prepare(old);assert full['verdict']=='HOLD_MODEL_INSTABILITY'
 influence=pd.read_csv(O/'leave_one_donor_out_influence.tsv',sep='\t');records=[]
 for donor in sorted(influence.donor_id.unique()):
  p=CK/('phase1_verdict_without_'+donor+'.json')
  if p.exists():result=json.loads(p.read_text())
  else:
   blocks=[]
   for t in m.TARGETS:
    for s in m.ALL_STRATA:
     path=CK/(f'loo_{donor}__{t}__{s}.tsv.gz')
     blocks.append(pd.read_csv(path,sep='\t') if path.exists() else old[old.target.eq(t)&old.stratum.eq(s)].copy())
   result=prepare(pd.concat(blocks,ignore_index=True));p.write_text(json.dumps(result,indent=2))
  records.append(dict(donor_id=donor,full_phase1_verdict=full['verdict'],loo_phase1_verdict=result['verdict'],verdict_change=result['verdict']!=full['verdict'],full_robust_targets=';'.join(full['robust_targets']),loo_robust_targets=';'.join(result['robust_targets']),loo_heterogeneity=result['heterogeneity']))
  print('PHASE1_VERDICT_DIAGNOSTIC',donor,result['verdict'],flush=True)
 # Separate actual stage-verdict change from the prospectively defined target CI/FDR interpretation change.
 influence['qualitative_interpretation_change']=influence.classification_change|influence.delta_sign_flip|influence.residual_sign_flip
 influence=influence.drop(columns='verdict_change').merge(pd.DataFrame(records),on='donor_id',validate='many_to_one')
 influence.to_csv(O/'leave_one_donor_out_influence.tsv',sep='\t',index=False);pd.DataFrame(records).to_csv(O/'loo_original_phase1_verdict.tsv',sep='\t',index=False)
if __name__=='__main__':main()
