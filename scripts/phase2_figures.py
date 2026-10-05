"""Render manuscript figures from frozen tables only. No model fitting or resampling."""
from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
R=Path(__file__).resolve().parents[1];O=R/'manuscript/figures';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6,'svg.fonttype':'none','pdf.fonttype':42,'savefig.dpi':300,'axes.spines.top':False,'axes.spines.right':False})
TS=['epithelial_injury','fibroblast_activation','macrophage_inflammatory']; names=['Epithelial injury','Fibroblast activation','Macrophage inflammatory']; colors=['#3B6FB6','#B06624','#6E4A8B']; W=183/25.4
sources={};qa=[]
def rd(f):sources[f]=hashlib.sha256((R/f).read_bytes()).hexdigest();return pd.read_csv(R/f,sep='\t')
p=rd('results/phase1/tables/table2_primary_model_performance.tsv');q=rd('results/phase1/tables/table3_crossfit_residual_performance.tsv');c=rd('results/phase1/tables/table4_shortcut_stress_test.tsv');h=rd('results/phase1h/repeated_fold_metrics.tsv');l=rd('results/phase1h/leave_one_donor_out_influence.tsv');i=rd('results/phase1i/inference_agreement_matrix.tsv')
def label(ax,letter,title):ax.set_title(title,loc='left');ax.text(-.08,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=8)
def save(fig,name,claims):
 fig.savefig(O/(name+'.svg'));fig.savefig(O/(name+'.png'),dpi=300);plt.close(fig)
 for a,claim in enumerate(claims):qa.append({'figure':name,'panel':chr(97+a),'claim':claim,'font_floor_pt':6,'replicate_unit':'Donor; partitions are repeated dependent splits','source_preflight':'PENDING','visual_QA':'PENDING'})
# Contract exists before plotting; schematic is a diagram, not biological mechanistic evidence.
(R/'results/final/figure_contract.json').write_text(json.dumps({'backend':'Python matplotlib','width_mm':183,'font_floor_pt':6,'exports':['editable SVG','300 dpi PNG'],'PDF':'not requested','new_effects_or_resampling':False,'main_figures':{'F1':'Framework and frozen estimand','F2':'Cohort and distinct technical diagnostics','F3':'Target-specific ladder including negative R²','F4':'Increment and cross-fitted residual intervals','F5':'Absolute control performance without attributing causality','F6':'Repeated partitions, refit LOO and uncertainty separated by stratum'},'supplement':'Reuse frozen microscopy/dependency assets without edits; replot all LOO; link full source tables','display':'No fabricated data; bootstrap intervals only where frozen; no inferential error bars for single frozen point ladder/control values'},indent=2))
fig,ax=plt.subplots(figsize=(W,4.3));ax.axis('off')
def box(x,y,w,hh,txt,color='#F2F4F7'):
 ax.add_patch(FancyBboxPatch((x,y),w,hh,boxstyle='round,pad=0.008',facecolor=color,edgecolor='#64748B',lw=.7));ax.text(x+w/2,y+hh/2,txt,ha='center',va='center',fontsize=7)
box(.02,.75,.27,.19,'Same-section registered H&E\n150 μm microregion\nX: frozen Phikon-v2')
box(.36,.75,.27,.19,'Measured Xenium cell labels\nC: 15-lineage composition\nN: cardinal-neighbor composition')
box(.70,.75,.27,.19,'Y: within-lineage programs\nZ: spatial coordinates\n+ disease in adjusted analysis')
box(.05,.39,.39,.27,'Frozen donor-held-out ladder\nM0 = Z       M1 = C + Z\nM2 = C + N + Z       M3 = X + Z\nM4 = C + N + Z + X\nΔR² = R²(M4) − R²(M2)')
box(.54,.39,.39,.27,'Training-side donor cross-fit\nM2 out-of-fold residuals → X-only ridge\nTest residual = Y − M2(test prediction)\nHeld-out residual R²\nNo full-data residualization')
box(.05,.06,.39,.23,'Controls and strata\nRaw / training-donor Macenko / grayscale\nColor-only / coordinates-only\nPooled / disease-adjusted / PF-only')
box(.54,.06,.39,.23,'Donor-level uncertainty\n19 donors, 26 sections\n50 partitions: pooled and PF-only\nRefit LOO / jackknife / balanced / BCa\nNo external validation')
for x in [.20,.50,.82]:ax.annotate('',xy=(x,.67),xytext=(x,.74),arrowprops={'arrowstyle':'->','color':'#64748B'})
ax.annotate('',xy=(.49,.53),xytext=(.44,.53),arrowprops={'arrowstyle':'->'})
ax.annotate('',xy=(.25,.30),xytext=(.25,.38),arrowprops={'arrowstyle':'->'})
ax.annotate('',xy=(.74,.30),xytext=(.74,.38),arrowprops={'arrowstyle':'->'})
fig.subplots_adjust(left=.02,right=.98,top=.98,bottom=.02);save(fig,'F1_framework',['Audit inputs, model ladder, residual estimand and limits'])
fig,axs=plt.subplots(1,3,figsize=(W,3.05),layout='constrained')
ax=axs[0];label(ax,'a','Independent units');ax.bar([0,1],[6,13],color=['#777777','#B06624']);ax.bar([0,1],[7,19],width=.35,facecolor='none',edgecolor='black',hatch='//');ax.set_xticks([0,1],['Control','PF']);ax.set_ylabel('Count');ax.set_ylim(0,23);ax.text(.5,21,'Solid: donors; hatch: sections',ha='center',fontsize=6)
ax=axs[1];label(ax,'b','Registration and common universe');ax.axis('off');ax.text(0,.95,'28 audited sections\n↓ fixed coordinate support\n26 admitted sections / 19 donors\n\n12,260 initial regions\n−77 full-crop coordinate failures\n−39 zero valid image-pixel regions\n=12,144 common regions\n\nTPS-registered author H&E\n150 μm qualification; no single-cell claim',va='top',fontsize=7)
ax=axs[2];label(ax,'c','Distinct shortcut diagnostics')
cf=rd('results/phase0/phase0d/color_shortcut_metrics.tsv');co=rd('results/phase0/phase0d/coordinate_shortcut_metrics.tsv');cp=rd('results/phase0/phase0d/composition_disease_coupling.tsv')
ss=[cf[cf.target=='disease'].iloc[0],cf[cf.target=='TMA'].iloc[0],co[co.target=='TMA'].iloc[0],cp.iloc[0]]
for jj,r in enumerate(ss):ax.errorbar(r.balanced_accuracy,jj,xerr=[[r.balanced_accuracy-r.bootstrap_95ci_low],[r.bootstrap_95ci_high-r.balanced_accuracy]],fmt='o',color='#3B6FB6',capsize=2)
ax.set_yticks(range(4),['Color → disease','Color → TMA','Coordinates → TMA','Composition → disease']);ax.invert_yaxis();ax.set_xlim(0,1.03);ax.set_xlabel('Balanced accuracy (95% donor CI)');save(fig,'F2_cohort_shortcuts',['Donors and sections differ','Transparent pre-model universe amendment','Strong nuisance/biological association; assays and metrics distinct'])
fig,axs=plt.subplots(1,3,figsize=(W,2.8),layout='constrained')
for a,t,col,n in zip(axs,TS,colors,names):
 r=p[(p.target==t)&(p.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0];a.plot(range(5),[r['M'+str(j)+'_R2'] for j in range(5)],'o-',color=col,lw=1.2);a.axhline(0,color='#888888',lw=.7);a.set_xticks(range(5),['M0','M1','M2','M3','M4']);a.set_ylim(-.53,.53);a.set_ylabel('Donor-held-out R²');label(a,chr(97+TS.index(t)),n)
save(fig,'F3_model_ladder',['Disease-adjusted five-model ladder; fixed point estimates']*3)
fig,axs=plt.subplots(1,2,figsize=(W,3),layout='constrained')
for ax,f,m,title in [(axs[0],p,'deltaR2_morphology','Morphology increment ΔR²'),(axs[1],q,'residual_R2','Cross-fitted residual R²')]:
 for j,(t,col) in enumerate(zip(TS,colors)):
  r=f[(f.target==t)&(f.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0];ax.errorbar(r[m],j,xerr=[[r[m]-r.CI_low],[r.CI_high-r[m]]],fmt='o',color=col,capsize=3)
 ax.axvline(0,color='#888888',lw=.7);ax.set_yticks(range(3),names);ax.invert_yaxis();ax.set_xlabel('Frozen estimate and 95% donor-bootstrap CI');label(ax,'a' if ax is axs[0] else 'b',title)
save(fig,'F4_increment_residual',['Primary adjusted increments; no target q < .05','Residual support remains qualified; no target q < .05'])
fig,axs=plt.subplots(1,3,figsize=(W,2.8),layout='constrained')
for a,t,col,n in zip(axs,TS,colors,names):
 r=c[(c.target==t)&(c.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0];vals=[r.raw,r.macenko,r.grayscale,r.color_only,r.coordinates_only];a.plot(range(5),vals,'o',color=col);a.axhline(0,color='#888888',lw=.7);a.set_xticks(range(5),['Raw','Macenko','Gray','Color','Coord'],rotation=35,ha='right',rotation_mode='anchor');a.set_ylim(-.28,.53);a.set_ylabel('Donor-held-out R²');label(a,chr(97+TS.index(t)),n)
save(fig,'F5_shortcut_controls',['Raw/Macenko/gray are augmented M4; color/coord are standalone controls']*3)
fig,axs=plt.subplots(3,3,figsize=(W,7.4),layout='constrained')
methods=['REFERENCE_BOOTSTRAP','JACKKNIFE_NORMAL_CI','BALANCED_DONOR_BOOTSTRAP','BCA_INFERENCE_SENSITIVITY'];mt=['Reference','Jackknife','Balanced','BCa']
for row,(t,n,col) in enumerate(zip(TS,names,colors)):
 ax=axs[row,0]
 for yy,st,mk in [(0,'POOLED','o'),(1,'PF_ONLY','s')]:
  z=h[(h.target==t)&(h.stratum==st)];ax.scatter(z.DeltaR2,np.full(len(z),yy),s=10,marker=mk,color=col,alpha=.5)
 ax.set_yticks([0,1],['Pooled','PF-only']);ax.set_ylim(-.5,1.5);ax.axvline(0,color='#888888',lw=.7);ax.set_xlabel('ΔR²; all 50 partitions');label(ax,chr(97+row*3),n+' / partitions')
 ax=axs[row,1];z=l[(l.target==t)&(l.stratum=='DISEASE_ADJUSTED')];ax.scatter(z.loo_DeltaR2,range(len(z)),s=12,color=col);ax.axvline(0,color='#888888',lw=.7);ax.set_yticks([]);ax.set_ylabel('19 omitted donors');ax.set_xlabel('Adjusted refit LOO ΔR²');label(ax,chr(98+row*3),'Donor deletion')
 ax=axs[row,2]
 for yy,method in enumerate(methods):
  z=i[(i.target==t)&(i.effect=='DeltaR2_morph')&(i.method==method)]
  if len(z)!=1:raise ValueError((t,method,len(z)))
  r=z.iloc[0]
  if pd.notna(r.CI_low) and pd.notna(r.CI_high):ax.hlines(yy,r.CI_low,r.CI_high,color=col);ax.plot(r.estimate,yy,'o',color=col,ms=3)
  else:ax.text(0,yy,'Invalid',fontsize=6)
 ax.set_yticks(range(4),mt);ax.invert_yaxis();ax.axvline(0,color='#888888',lw=.7);ax.set_xlabel('Adjusted ΔR² method-specific CI');label(ax,chr(99+row*3),'Inference sensitivity')
save(fig,'F6_donor_inference_stability',['Dependent repeated partitions, not new donors','All 19 adjusted refit omissions','Method estimates/intervals, not one common center']*3)
# Supplementary assets preserved byte-for-byte, not edited or filtered.
reuse={'S1_registration_example':'figures/phase0c_vr/THD0008_fixed_patches.png','S2_donor_dependency':'figures/phase0d/vannan_dependency_graph.png','S4_all_partition_increments':'results/phase1h/figures/H3_repeated_fold_delta.png','S5_all_partition_residuals':'results/phase1h/figures/H4_repeated_fold_residual.png','S8_m2_support':'results/phase1h/figures/H2_m2_support_donor_error.png','S9_jackknife_pseudovalues':'results/phase1i/figures/I3_jackknife_donor_pseudovalues.png','S10_residual_inference':'results/phase1i/figures/I2_residual_inference_comparison.png','S11_inference_agreement':'results/phase1i/figures/I5_inference_agreement.png','S12_precision_projection':'figures/phase1j/J1_precision_projection.png'}
for name,src in reuse.items():
 f=R/src;shutil.copyfile(f,O/(name+'.png'));sources[src]=hashlib.sha256(f.read_bytes()).hexdigest()
 if f.with_suffix('.svg').exists():shutil.copyfile(f.with_suffix('.svg'),O/(name+'.svg'));sources[str(f.with_suffix('.svg').relative_to(R))]=hashlib.sha256(f.with_suffix('.svg').read_bytes()).hexdigest()
 qa.append({'figure':name,'panel':'ALL inherited panels','claim':'Frozen audit asset; see original report and supplemental legend','font_floor_pt':'Inherited','replicate_unit':'Original frozen assay','source_preflight':'Inherited source available','visual_QA':'PENDING'})
fig,ax=plt.subplots(figsize=(W,2.5),layout='constrained');ax.barh(['Common universe','Coordinate support excluded','Black-filter excluded'],[12144,77,39],color=['#777777','#3B6FB6','#B06624']);ax.set_xlabel('Regions; 116 excluded, no overlap');ax.set_title('Pre-model amendment: same universe for all model comparisons');save(fig,'S3_common_universe_QC',['Fixed technical exclusions before predictive results'])
for metric,name in [('loo_DeltaR2','S6_all_donor_LOO_increment'),('loo_residual_R2','S7_all_donor_LOO_residual')]:
 fig,axs=plt.subplots(3,3,figsize=(W,7.7),layout='constrained')
 for row,t in enumerate(TS):
  for cc,st in enumerate(['POOLED','DISEASE_ADJUSTED','PF_ONLY']):
   ax=axs[row,cc];z=l[(l.target==t)&(l.stratum==st)&(~l.donor_not_in_stratum.astype(str).str.lower().eq('true'))].sort_values('donor_id');ax.scatter(z[metric],range(len(z)),s=10,color=colors[row]);ax.set_yticks(range(len(z)),z.donor_id,fontsize=6);ax.axvline(0,color='#888888',lw=.7);ax.set_xlabel('ΔR²' if metric=='loo_DeltaR2' else 'Residual R²');label(ax,chr(97+row*3+cc),names[row]+'\n'+st.replace('_',' ').title())
 save(fig,name,['All actual refit donor omissions; PF-only n=13, pooled/adjusted n=19']*9)
(R/'results/final/figure_source_hashes.json').write_text(json.dumps(sources,indent=2));pd.DataFrame(qa).to_csv(R/'results/final/figure_panel_QA.tsv',sep='\t',index=False)
print('Rendered 6 main and 12 supplement figures; no models run.')
