#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarize the prospective audit and render required source-linked figures."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy.stats import percentileofscore
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];O=ROOT/'results/phase1h';FIG=O/'figures';FIG.mkdir(exist_ok=True)
TS=['epithelial_injury','fibroblast_activation','macrophage_inflammatory'];NAMES=['Epithelial injury','Fibroblast activation','Macrophage inflammatory']
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'legend.frameon':False})
def read(name):return pd.read_csv(O/name,sep='\t')
def write(df,name):df.to_csv(O/name,sep='\t',index=False)
def stability(v):
 pos=np.mean(v>0);neg=np.mean(v<0);s=max(pos,neg);return s,'HIGH' if s>=.8 else 'MODERATE' if s>=.65 else 'LOW',pos,neg

def summarize():
 assert json.loads((O/'analysis_complete.json').read_text())['status']=='PHASE1H_ANALYSIS_COMPLETE'
 repeat=read('repeated_fold_metrics.tsv');full=read('full_frozen_metrics.tsv');loo=read('leave_one_donor_out_influence.tsv');paired=read('m2_m4_paired_leverage.tsv');assoc=read('support_error_association.tsv')
 summaries=[]
 for (t,s),g in repeat.groupby(['target','stratum']):
  b=full[full.target.eq(t)&full.stratum.eq(s)].iloc[0];row=dict(target=t,stratum=s,n_partitions=len(g),median_DeltaR2=g.DeltaR2.median(),DeltaR2_IQR=g.DeltaR2.quantile(.75)-g.DeltaR2.quantile(.25),DeltaR2_p5=g.DeltaR2.quantile(.05),DeltaR2_p95=g.DeltaR2.quantile(.95),median_residual_R2=g.residual_R2.median(),residual_p5=g.residual_R2.quantile(.05),residual_p95=g.residual_R2.quantile(.95),median_M2_R2=g.M2_R2.median(),median_M4_R2=g.M4_R2.median(),median_residual_Pearson=g.residual_Pearson.median(),median_residual_Spearman=g.residual_Spearman.median(),full_DeltaR2=b.DeltaR2)
  for metric,label in [('DeltaR2','delta'),('residual_R2','residual')]:
   st,cat,pos,neg=stability(g[metric]);row.update({label+'_sign_stability':st,label+'_stability_class':cat,label+'_fraction_positive':pos,label+'_fraction_negative':neg})
   # CV omitted when mean is close to zero to avoid meaningless inflation.
   row[label+'_CV']=g[metric].std()/abs(g[metric].mean()) if abs(g[metric].mean())>.02 else np.nan
  row['original_fold_percentile']=percentileofscore(g.DeltaR2,b.DeltaR2,kind='mean');row['original_fold_extreme']=row['original_fold_percentile']<5 or row['original_fold_percentile']>95;summaries.append(row)
 summary=pd.DataFrame(summaries);write(summary,'repeated_fold_stability.tsv')
 classes=[]
 for t in TS:
  g=loo[loo.target.eq(t)&loo.stratum.eq('DISEASE_ADJUSTED')];b=full[full.target.eq(t)&full.stratum.eq('DISEASE_ADJUSTED')].iloc[0];s=summary[summary.target.eq(t)&summary.stratum.eq('POOLED')].iloc[0];p=paired[paired.target.eq(t)&paired.stratum.eq('DISEASE_ADJUSTED')];a=assoc[assoc.target.eq(t)&assoc.stratum.eq('DISEASE_ADJUSTED')].iloc[0]
  flips=g[g.qualitative_interpretation_change];unsupported=set(p[p.unsupported_donor].donor_id);positive=p.paired_SSE_gain.clip(lower=0);concentration=positive[p.donor_id.isin(unsupported)].sum()/positive.sum() if positive.sum()>0 else 0
  reduction=g[g.donor_id.isin(unsupported)].delta_influence.min() if unsupported else np.nan
  extrap=bool(a.donor_Spearman>=.5 and a.donor_p<.05 and concentration>=.5 and np.isfinite(reduction) and reduction<=-.5*abs(b.DeltaR2))
  null=bool(s.DeltaR2_p5>=-.02 and s.DeltaR2_p95<=.02 and s.residual_p5>=-.02 and s.residual_p95<=.02 and abs(s.median_DeltaR2)<=.01 and abs(s.median_residual_R2)<=.01 and b.delta_CI_low>=-.05 and b.delta_CI_high<=.05 and b.residual_CI_low>=-.05 and b.residual_CI_high<=.05 and g.delta_influence.abs().max()<=.02 and g.residual_influence.abs().max()<=.02)
  stable=bool(s.delta_sign_stability>=.8 and s.residual_sign_stability>=.8 and len(flips)==0 and not s.original_fold_extreme and not extrap)
  donor=bool(0<len(flips)<=2 and min(s.delta_sign_stability,s.residual_sign_stability)<.8 and not extrap)
  # A persistent sign alone cannot establish positive predictive evidence.
  if extrap:classification='M2_EXTRAPOLATION_DRIVEN'
  elif null:classification='STABLE_NULL'
  elif stable:classification='STABLE_SIGNAL'
  elif donor:classification='DONOR_LEVERAGE_DRIVEN'
  elif min(s.delta_sign_stability,s.residual_sign_stability)<.8:classification='FOLD_INSTABILITY'
  else:classification='STABILITY_NOT_ESTABLISHED' # None of the specified positive classifications is established; do not falsely attribute to <=2 donors.
  high=g.loc[g.delta_influence.abs().idxmax()]
  classes.append(dict(target=t,classification=classification,full_DeltaR2=b.DeltaR2,loo_DeltaR2_min=g.DeltaR2.min(),loo_DeltaR2_max=g.DeltaR2.max(),delta_sign_flips=int(g.delta_sign_flip.sum()),residual_sign_flips=int(g.residual_sign_flip.sum()),interpretation_flips=len(flips),highest_leverage_donor=high.donor_id,maximum_abs_delta_influence=abs(high.delta_influence),median_DeltaR2=s.median_DeltaR2,DeltaR2_p5=s.DeltaR2_p5,DeltaR2_p95=s.DeltaR2_p95,delta_sign_stability=s.delta_sign_stability,residual_sign_stability=s.residual_sign_stability,median_residual_R2=s.median_residual_R2,donor_support_Spearman=a.donor_Spearman,donor_support_p=a.donor_p,unsupported_positive_gain_fraction=concentration,unsupported_LOO_max_reduction=reduction,original_fold_percentile=s.original_fold_percentile,original_fold_extreme=s.original_fold_extreme,stable_direction=stable,positive_signal=stable and s.median_DeltaR2>0 and s.median_residual_R2>0,extrapolation_criteria=extrap,stable_null_criteria=null,low_repeat_signs=min(s.delta_sign_stability,s.residual_sign_stability)<.65))
 cl=pd.DataFrame(classes)
 positive=cl.positive_signal.sum()
 if positive==3:verdict='STABILITY_PASS_GO_EXTERNAL_VALIDATION'
 elif positive>0:verdict='SELECTIVE_SIGNAL_WITH_HETEROGENEITY'
 elif cl.extrapolation_criteria.any():verdict='HOLD_M2_EXTRAPOLATION'
 elif cl.classification.eq('DONOR_LEVERAGE_DRIVEN').any():verdict='HOLD_DONOR_LEVERAGE'
 elif cl.stable_null_criteria.sum()>=2:verdict='STABLE_NULL_GO_COMPOSITION_DOMINANCE'
 elif cl.classification.eq('FOLD_INSTABILITY').any():verdict='HOLD_FOLD_INSTABILITY'
 else:verdict=None # No supplied stage category is established; require a human decision instead of forcing fold instability.
 write(cl,'target_stability_classification.tsv')
 for name,df in [('tableH1_donor_influence.tsv',loo),('tableH2_m2_support.tsv',read('m2_support_distance.tsv')),('tableH3_repeated_stability.tsv',summary),('tableH4_target_classification.tsv',cl)]:write(df,name)
 # Low-variance co-occurrence is descriptive and never used for exclusion.
 error=read('donor_error_profile.tsv');joint=error.merge(loo[['target','stratum','donor_id','delta_influence','residual_influence','leverage_category']],on=['target','stratum','donor_id'],validate='one_to_one')
 for _,g in joint.groupby(['target','stratum']):
  joint.loc[g.index,'extreme_M2_R2']=g.M2_R2.le(g.M2_R2.quantile(.1))|g.M2_R2.ge(g.M2_R2.quantile(.9));joint.loc[g.index,'extreme_delta_influence']=g.delta_influence.abs().ge(g.delta_influence.abs().quantile(.9));joint.loc[g.index,'extreme_residual_influence']=g.residual_influence.abs().ge(g.residual_influence.abs().quantile(.9))
 write(joint,'low_variance_influence_cooccurrence.tsv')
 result={'verdict':verdict,'external_validation_authorized':'YES' if positive else 'NO','positive_stable_targets':cl[cl.positive_signal].target.tolist(),'n_extreme_donors':int(loo[loo.stratum.eq('DISEASE_ADJUSTED')&loo.leverage_category.eq('EXTREME_LEVERAGE')].donor_id.nunique()),'highest_leverage_donor':cl.loc[cl.maximum_abs_delta_influence.idxmax(),'highest_leverage_donor'],'next_phase':{'HOLD_M2_EXTRAPOLATION':'Pre-register support-aware estimand / restricted-support sensitivity in a separate phase.','HOLD_DONOR_LEVERAGE':'Pre-register donor-leverage and transportability design; retain all current donors.','HOLD_FOLD_INSTABILITY':'Pre-register a fold-robust estimand and uncertainty design before validation.','STABLE_NULL_GO_COMPOSITION_DOMINANCE':'Separate composition-dominance evidence phase.','SELECTIVE_SIGNAL_WITH_HETEROGENEITY':'Consider independently authorized target-specific external validation.','STABILITY_PASS_GO_EXTERNAL_VALIDATION':'Consider independently authorized external validation.'}.get(verdict,'Human decision required: HIGH repeat sign consistency coexists with donor-sensitive CI/FDR interpretation; supplied categories do not cover this combination.'),'decision_status':'DECIDED' if verdict else 'HUMAN_CLASSIFICATION_REQUIRED'}
 (O/'verdict.json').write_text(json.dumps(result,indent=2));return full,loo,paired,summary,cl,result

def save(fig,name):
 fig.savefig(FIG/(name+'.svg'));fig.savefig(FIG/(name+'.png'),dpi=300);plt.close(fig)
def label(ax,i,title):
 if title in NAMES:title=title.replace('Fibroblast activation','Fibroblast\nactivation').replace('Macrophage inflammatory','Macrophage\ninflammatory')
 ax.set_title(chr(97+i)+'  '+title,loc='left',fontweight='bold',pad=8)
def figures(full,loo,paired,summary):
 # H1 raw diagnostic influences for all donors; zero reference, no confidence intervals on deterministic refits.
 fig,axes=plt.subplots(3,1,figsize=(183/25.4,6.5));fig.subplots_adjust(left=.11,right=.98,top=.95,bottom=.09,hspace=.8)
 for i,(t,ax) in enumerate(zip(TS,axes)):
  g=loo[loo.target.eq(t)&loo.stratum.eq('DISEASE_ADJUSTED')].sort_values('donor_id');x=np.arange(len(g));ax.axhline(0,color='.6',lw=.7);ax.vlines(x,0,g.delta_influence,color='#537C99',lw=1);ax.scatter(x,g.delta_influence,c=['#AD655A' if k=='EXTREME_LEVERAGE' else '#537C99' for k in g.leverage_category],s=17,zorder=3)
  ax.set_xticks(x);ax.set_xticklabels(g.donor_id,rotation=65,ha='right',rotation_mode='anchor',fontsize=6);ax.set_ylabel('LOO minus full ΔR²');label(ax,i,NAMES[i]+' | 19 diagnostic refits')
 save(fig,'H1_leave_one_donor_influence')
 for name,mode in [('H2_m2_support_donor_error','support'),('H5_paired_m2_m4_error','paired')]:
  fig,axes=plt.subplots(1,3,figsize=(183/25.4,2.75));fig.subplots_adjust(left=.09,right=.98,top=.8,bottom=.23,wspace=.5)
  for i,(t,ax) in enumerate(zip(TS,axes)):
   g=paired[paired.target.eq(t)&paired.stratum.eq('DISEASE_ADJUSTED')];x=g.mahalanobis_median if mode=='support' else g.M2_MAE;y=g.M2_MAE if mode=='support' else g.M4_MAE
   ax.scatter(x,y,s=25,c='#537C99',edgecolor='white',linewidth=.4);ax.set_xlabel('Median M2 support distance' if mode=='support' else 'M2 MAE');ax.set_ylabel('M2 MAE' if mode=='support' else 'M4 MAE');label(ax,i,NAMES[i]);
   if mode=='paired':
    lim=max(x.max(),y.max())*1.12;ax.plot([0,lim],[0,lim],color='.55',lw=.7,ls='--');ax.set_xlim(0,lim);ax.set_ylim(0,lim)
   # Label the two largest errors using complete donor IDs; all donor points retained.
   for rank,j in enumerate(sorted(np.argsort(y.to_numpy())[-2:],key=lambda j:x.iloc[j])):
    inward=x.iloc[j]>=x.max()-.15*(x.max()-x.min());left=True if mode=='paired' else (rank==0 or inward)
    ax.annotate(g.donor_id.iloc[j],(x.iloc[j],y.iloc[j]),xytext=(-4 if left else 4,(9 if rank==0 else -5) if mode=='paired' else (9 if rank==0 else 15)),ha='right' if left else 'left',textcoords='offset points',fontsize=5.5,arrowprops={'arrowstyle':'-','lw':.4,'color':'.5'})
   ax.margins(x=.2,y=.25)
  save(fig,name)
 for name,metric,ylabel in [('H3_repeated_fold_delta','DeltaR2','Morphology increment ΔR²'),('H4_repeated_fold_residual','residual_R2','Cross-fit residual R²')]:
  rep=read('repeated_fold_metrics.tsv');fig,axes=plt.subplots(1,3,figsize=(183/25.4,3.0));fig.subplots_adjust(left=.1,right=.98,top=.8,bottom=.2,wspace=.5)
  for i,(t,ax) in enumerate(zip(TS,axes)):
   ax.axhline(0,color='.65',lw=.7)
   for j,s in enumerate(['POOLED','PF_ONLY']):
    g=rep[rep.target.eq(t)&rep.stratum.eq(s)].sort_values('seed');v=g[metric].to_numpy();x=j+np.linspace(-.15,.15,len(v));ax.scatter(x,v,s=9,color='#819FAF' if j==0 else '#B0A09B',alpha=.7)
    q=np.quantile(v,[.05,.25,.5,.75,.95]);ax.errorbar(j,q[2],yerr=[[q[2]-q[0]],[q[4]-q[2]]],fmt='o',color='#355F79',ms=4,capsize=3,lw=1.2);ax.plot([j,j],[q[1],q[3]],lw=4,color='#355F79')
    b=full[full.target.eq(t)&full.stratum.eq(s)].iloc[0];ax.scatter(j+.24,b[metric],marker='D',s=22,facecolor='white',edgecolor='#AD655A',linewidth=1.2,zorder=5)
   ax.set_xticks([0,1]);ax.set_xticklabels(['Pooled','PF only']);ax.set_ylabel(ylabel);label(ax,i,NAMES[i]);ax.margins(x=.23,y=.1)
  fig.text(.5,.96,'50 partitions: median, IQR and 5–95% range; open diamond = original frozen fold',ha='center',fontsize=7)
  save(fig,name)

if __name__=='__main__':
 full,loo,paired,summary,cl,res=summarize();figures(full,loo,paired,summary)
 print(json.dumps(res,indent=2))
