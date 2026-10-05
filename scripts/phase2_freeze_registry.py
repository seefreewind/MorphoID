"""Copy frozen estimates and source hashes; never fit models or alter evidence."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
R=Path(__file__).resolve().parents[1]; O=R/'results/final';O.mkdir(exist_ok=True)
def read(f):return pd.read_csv(R/f,sep='\t')
def sha(f):return hashlib.sha256((R/f).read_bytes()).hexdigest()
P='results/phase1/tables/table2_primary_model_performance.tsv';Q='results/phase1/tables/table3_crossfit_residual_performance.tsv';C='results/phase1/tables/table4_shortcut_stress_test.tsv';H='results/phase1h/repeated_fold_stability.tsv';L='results/phase1h/leave_one_donor_out_influence.tsv';I='results/phase1i/inference_agreement_matrix.tsv'
p,q,c,h,l,i=map(read,[P,Q,C,H,L,I]);levels={'epithelial_injury':'C','fibroblast_activation':'B','macrophage_inflammatory':'D'}
rows=[]
def add(phase,analysis,f,r,metric,value,lo=np.nan,hi=np.nan,pv=np.nan,qv=np.nan,stratum=None):
 t=r.get('target','ALL');rows.append(dict(phase=phase,analysis=analysis,target=t,stratum=stratum or r.get('analysis_stratum',r.get('stratum','DIAGNOSTIC')),metric=metric,estimate=value,CI_low=lo,CI_high=hi,p=pv,FDR=qv,n_donors=r.get('n_donors',np.nan),n_sections=r.get('n_sections',np.nan),n_regions=r.get('n_regions',np.nan),evidence_level=levels.get(t,'DIAGNOSTIC'),claim_status='DESCRIPTIVE_ONLY' if analysis not in ['PRIMARY_INCREMENT','RESIDUAL'] else 'SUPPORTED_WITH_QUALIFICATION',source_file=f,source_hash=sha(f)))
for _,r in p.iterrows():
 for m in [x for x in p.columns if x.endswith(('_R2','_Pearson_r','_Spearman_rho','_MAE','_RMSE'))]+['deltaR2_composition','deltaR2_neighborhood']:
  add('1','MODEL_LADDER',P,r,m,r[m])
 add('1','PRIMARY_INCREMENT',P,r,'deltaR2_morphology',r.deltaR2_morphology,r.CI_low,r.CI_high,r.bootstrap_two_sided_sign_p,r.FDR)
for _,r in q.iterrows():
 add('1','RESIDUAL',Q,r,'residual_R2',r.residual_R2,r.CI_low,r.CI_high,r.p,r.FDR)
 for m in ['Pearson','Spearman','MAE','RMSE']:add('1','RESIDUAL_DIAGNOSTIC',Q,r,m,r[m])
for _,r in c.iterrows():
 for m in c.columns[2:]:add('1','SHORTCUT_CONTROL',C,r,m,r[m])
for _,r in h.iterrows():
 for m in h.columns[3:]:
  if isinstance(r[m],(int,float,np.number)):add('1H','REPEATED_PARTITION',H,r,m,r[m])
for _,r in read('results/phase1h/repeated_fold_metrics.tsv').iterrows():
 for m in ['M2_R2','M4_R2','DeltaR2','residual_R2']:add('1H','PARTITION_'+str(r['run']),'results/phase1h/repeated_fold_metrics.tsv',r,m,r[m])
for _,r in l.iterrows():
 if str(r.donor_not_in_stratum).lower()=='true':continue
 for m,low,high,pv,qv in [('loo_DeltaR2','delta_CI_low','delta_CI_high','delta_p','delta_q'),('loo_residual_R2','residual_CI_low','residual_CI_high','residual_p','residual_q')]:add('1H','REFIT_LOO_'+str(r.donor_id),L,r,m,r[m],r[low],r[high],r[pv],r[qv])
for _,r in i.iterrows():add('1I',str(r.method),I,r,str(r.effect)+'_'+str(r.loss),r.estimate,r.CI_low,r.CI_high,r.p_value,r.FDR,'DISEASE_ADJUSTED')
for f in ['results/phase0/phase0d/color_shortcut_metrics.tsv','results/phase0/phase0d/coordinate_shortcut_metrics.tsv','results/phase0/phase0d/composition_disease_coupling.tsv']:
 for _,r in read(f).iterrows():
  for m in ['balanced_accuracy','macro_auroc']:
   add('0D','SHORTCUT_DIAGNOSTIC',f,r,m,r.get(m,np.nan),r.get('bootstrap_95ci_low',np.nan) if m=='balanced_accuracy' else np.nan,r.get('bootstrap_95ci_high',np.nan) if m=='balanced_accuracy' else np.nan,r.get('permutation_p',np.nan) if m=='balanced_accuracy' else np.nan)
# Preserve all numerical summaries from the completed precision/feasibility phase.
for f in sorted((R/'results/phase1j').glob('*.tsv')):
 if f.name.startswith('._') or f.name=='input_hash_gate.tsv':continue
 d=pd.read_csv(f,sep='\t');rel=str(f.relative_to(R))
 for idx,r in d.iterrows():
  for m in d.select_dtypes(include='number').columns:
   add('1J','FROZEN_'+f.stem+'_row'+str(idx),rel,r,m,r[m],stratum=str(r.get('stratum','PRECISION_OR_FEASIBILITY_DIAGNOSTIC')))
reg=pd.DataFrame(rows);reg.to_csv(O/'master_result_registry.tsv',sep='\t',index=False,na_rep='NA')
classrows=[]
s=read('results/phase1h/target_stability_classification.tsv')
for t,level in levels.items():
 rr=p[(p.target==t)&(p.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0];res=q[(q.target==t)&(q.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0]
 v=dict(target=t,FINAL_EVIDENCE_LEVEL=level,FINAL_CLAIM_STATUS='SUPPORTED_WITH_QUALIFICATION',scope='Incremental direction only; not selective identifiability' if level=='B' else 'Limited positive tendency; residual unsupported' if level=='C' else 'Bounded negative result; not universal impossibility',primary_stratum='DISEASE_ADJUSTED',n_donors=rr.n_donors,n_sections=rr.n_sections,n_regions=rr.n_regions)
 for m in ['M2_R2','M4_R2','deltaR2_morphology','CI_low','CI_high','bootstrap_two_sided_sign_p','FDR','color_R2','coordinate_R2','grayscale_R2','macenko_R2']:v['disease_adjusted_'+m]=rr[m]
 for m in ['residual_R2','CI_low','CI_high','p','FDR']:v['disease_adjusted_residual_'+m]=res[m]
 for st in ['POOLED','PF_ONLY']:
  hr=h[(h.target==t)&(h.stratum==st)].iloc[0];pr=p[(p.target==t)&(p.analysis_stratum==st)].iloc[0];qr=q[(q.target==t)&(q.analysis_stratum==st)].iloc[0]
  for m in ['median_DeltaR2','DeltaR2_p5','DeltaR2_p95','delta_fraction_positive','delta_fraction_negative','median_residual_R2','residual_p5','residual_p95','residual_fraction_positive','residual_fraction_negative']:v[st+'_repeat_'+m]=hr[m]
  for m in ['M2_R2','M4_R2','deltaR2_morphology','CI_low','CI_high','FDR','macenko_R2','grayscale_R2','color_R2','coordinate_R2']:v[st+'_'+m]=pr[m]
  for m in ['residual_R2','CI_low','CI_high']:v[st+'_residual_'+m]=qr[m]
 sr=s[s.target==t].iloc[0]
 for m in s.columns[1:]:v['LOO_or_H_class_'+m]=sr[m]
 v['phase1i_all_methods_json']=i[i.target==t].to_json(orient='records',double_precision=15)
 v['source_files_json']=json.dumps({f:sha(f) for f in [P,Q,C,H,L,I,'results/phase1h/target_stability_classification.tsv']})
 classrows.append(v)
pd.DataFrame(classrows).to_csv(O/'target_evidence_classification.tsv',sep='\t',index=False,na_rep='NA')
print('Registry rows',len(reg),'targets',len(classrows))
