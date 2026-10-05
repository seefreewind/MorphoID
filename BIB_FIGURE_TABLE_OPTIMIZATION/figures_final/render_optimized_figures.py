"""Presentation-only redraw of all frozen rows for Figures 1–5; no fitting or inference."""
from pathlib import Path
import hashlib, json, shutil
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
R=Path(__file__).resolve().parents[2]; O=Path(__file__).resolve().parent
W=183/25.4
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
 'font.size':9.5,'axes.titlesize':10,'axes.labelsize':9.5,'xtick.labelsize':9.5,'ytick.labelsize':9.5,
 'legend.fontsize':9.5,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,
 'axes.spines.right':False,'axes.linewidth':.8,'savefig.dpi':601})
TS=['epithelial_injury','fibroblast_activation','macrophage_inflammatory']
NAMES=['Epithelial injury','Fibroblast activation','Macrophage inflammatory']
COLORS=['#3B6FB6','#B06624','#6E4A8B']; sources={}; exports=[]; plotted=[]
def read(f):
 p=R/f; sources[f]=hashlib.sha256(p.read_bytes()).hexdigest()
 shutil.copyfile(p,O.parent/'source_data'/p.name)
 return pd.read_csv(p,sep='\t')
p=read('results/phase1/tables/table2_primary_model_performance.tsv')
q=read('results/phase1/tables/table3_crossfit_residual_performance.tsv')
c=read('results/phase1/tables/table4_shortcut_stress_test.tsv')
cf=read('results/phase0/phase0d/color_shortcut_metrics.tsv')
co=read('results/phase0/phase0d/coordinate_shortcut_metrics.tsv')
cp=read('results/phase0/phase0d/composition_disease_coupling.tsv')
def label(ax,letter,title):
 ax.set_title(title,loc='left',pad=10); ax.text(-.1,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=10)
def save(fig,name,contract):
 fig.savefig(O/(name+'.svg')); fig.savefig(O/(name+'.pdf')); fig.savefig(O/(name+'.png'),dpi=601)
 text=[t.get_fontsize() for ax in fig.axes for t in ax.texts+[ax.title,ax.xaxis.label,ax.yaxis.label]+ax.get_xticklabels()+ax.get_yticklabels() if t.get_text()]
 exports.append({'figure':name,'width_mm':183,'height_mm':fig.get_size_inches()[1]*25.4,
                 'raster_dpi':601,'font_floor_pt':min(text),'editable':'SVG text elements',
                 'contract':contract,'new_estimates':False,'resampling':False})
 plt.close(fig)
fig,ax=plt.subplots(figsize=(W,5.1));ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1)
def box(x,y,w,h,text,fc='#F2F4F7'):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008',facecolor=fc,edgecolor='#64748B',lw=.8))
 ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=9.5,linespacing=1.5)
ax.text(.5,.98,'What does morphology add beyond C + N + Z?',ha='center',fontsize=11,fontweight='bold')
box(.02,.77,.27,.16,'Histology\nX: fixed morphology\nrepresentation')
box(.36,.77,.27,.16,'Measured context\nC: composition\nN: neighborhood/context')
box(.70,.77,.27,.16,'Target and covariates\nY: predefined target\nZ: prespecified covariates')
box(.04,.41,.41,.28,'Conditional model comparison\nM0 = Z     M1 = C + Z\nM2 = C + N + Z     M3 = X + Z\nM4 = C + N + Z + X\nΔR² = R²(M4) − R²(M2)')
box(.55,.41,.41,.28,'Cross-fitted residual test\nTraining-side grouped cross-fitting\nM2 out-of-fold residuals → X-only model\nTest residual = Y − M2(test prediction)\nHeld-out residual R²\nNo full-data residualization')
box(.04,.05,.41,.26,'Technical shortcut audit\nStain normalization\nGrayscale\nColor / spatial controls')
box(.55,.05,.41,.26,'Independent-unit uncertainty\nGrouped resampling\nRepeated grouped partitions\nRefit unit omission\nInference-sensitivity audit')
for x in [.16,.50,.84]:ax.annotate('',xy=(x,.70),xytext=(x,.76),arrowprops={'arrowstyle':'->','color':'#64748B','lw':.9})
ax.annotate('',xy=(.54,.55),xytext=(.46,.55),arrowprops={'arrowstyle':'->','lw':.9})
ax.text(.5,.35,'Interpretation safeguards',ha='center',fontsize=9.5,fontweight='bold')
for x in [.25,.75]:ax.annotate('',xy=(x,.32),xytext=(x,.40),arrowprops={'arrowstyle':'->','lw':.9})
fig.subplots_adjust(left=.015,right=.985,top=.985,bottom=.015)
save(fig,'Figure1_MorphoID_Framework','General conditional protocol; arrows are workflow, not causal edges. No cohort-specific counts.')
fig=plt.figure(figsize=(W,5.6));gs=fig.add_gridspec(2,2,width_ratios=[1,1.2],height_ratios=[1,1],hspace=.58,wspace=.30)
a=fig.add_subplot(gs[0,0]);a.axis('off');label(a,'a','Empirical cohort')
for y,name,donors,sections in [(.68,'Control',6,7),(.23,'PF',13,19)]:
 a.text(0,y,name,fontweight='bold',fontsize=10)
 a.text(0,y-.14,f'{donors} donors / {sections} sections',fontsize=9.5)
a.text(0,.01,'Donors and sections are distinct units',fontsize=9.5)
b=fig.add_subplot(gs[:,1]);b.axis('off');label(b,'b','Section and region inclusion')
flow=[(.92,'SECTION FLOW'),(.85,'28 audited sections'),(.75,'−2 fixed-coordinate-support failures'),(.65,'26 admitted sections / 19 donors'),(.49,'REGION FLOW'),(.42,'12,260 initial nonempty regions'),(.33,'−77 crop-support failures'),(.25,'−39 zero-valid-image regions'),(.15,'12,144 common regions')]
for yy,txt in flow:b.text(.5,yy,txt,ha='center',va='center',fontsize=9.5,fontweight='bold' if yy in [.92,.65,.49,.15] else 'normal')
for high,low in [(.82,.78),(.72,.68),(.39,.36),(.30,.28),(.22,.18)]:b.annotate('',xy=(.5,low),xytext=(.5,high),arrowprops={'arrowstyle':'->','lw':.8,'color':'#64748B'})
b.text(.5,.015,'Technical exclusions; no overlap\n150 μm regional analysis',ha='center',fontsize=9.5)
a=fig.add_subplot(gs[1,0]);label(a,'c','Shortcut diagnostics')
rows=[cf[cf.target=='disease'].iloc[0],cp.iloc[0],cf[cf.target=='TMA'].iloc[0],co[co.target=='TMA'].iloc[0]]
ys=[3,2,0,-1]
for j,(r,y) in enumerate(zip(rows,ys)):
 a.errorbar(r.balanced_accuracy,y,xerr=[[r.balanced_accuracy-r.bootstrap_95ci_low],[r.bootstrap_95ci_high-r.balanced_accuracy]],fmt='o',color='#4C5866',capsize=3,lw=1.2,ms=4)
 plotted.append({'figure':'F2','panel':'c','row':j,'estimate':r.balanced_accuracy,'low':r.bootstrap_95ci_low,'high':r.bootstrap_95ci_high})
a.set_yticks(ys,['Color → disease','Composition → disease','Color → TMA','Coordinates → TMA']);a.set_xlim(0,1.03);a.set_ylim(-1.6,4.2)
a.text(0,3.7,'Disease-associated',fontsize=9.5,fontweight='bold');a.text(0,.7,'Technical identity',fontsize=9.5,fontweight='bold')
a.set_xlabel('Balanced accuracy\n95% donor-bootstrap CI')
fig.subplots_adjust(left=.29,right=.985,top=.93,bottom=.13)
save(fig,'Figure2_Cohort_Shortcuts','Distinct classification tasks grouped by disease vs technical identity; no chance values inferred.')
fig,axs=plt.subplots(1,3,figsize=(W,3.6),layout='constrained')
for j,(a,t,col,name) in enumerate(zip(axs,TS,COLORS,NAMES)):
 r=p[(p.target==t)&(p.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0]
 vals=[r[f'M{k}_R2'] for k in range(5)]
 a.plot(range(5),vals,'o',color=col,ms=5,linestyle='none');a.axhline(0,color='#606060',lw=1)
 a.axvspan(-.45,2.45,color='#F1F2F3',zorder=0);a.text(.02,-.34,'Measured context: M0–M2\nMorphology-only: M3\nCombined: M4',transform=a.transAxes,fontsize=9.5)
 a.set_xticks(range(5),['M0','M1','M2','M3','M4']);a.set_ylim(-.53,.53);a.set_ylabel('Donor-held-out R²');label(a,chr(97+j),name)
 plotted.append({'figure':'F3','target':t,'values':vals})
save(fig,'Figure3_Model_Comparison','Five independent fixed estimates per target, no M2–M3 connectivity or invented confidence intervals.')
fig,axs=plt.subplots(1,2,figsize=(W,3.9),layout='constrained')
for ax,f,m,title,lim in [(axs[0],p,'deltaR2_morphology','Increment test: M4 − M2',(-.6,2.4)),(axs[1],q,'residual_R2','Residual test',(-.5,.12))]:
 for j,(t,col) in enumerate(zip(TS,COLORS)):
  r=f[(f.target==t)&(f.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0]
  ax.errorbar(r[m],j,xerr=[[r[m]-r.CI_low],[r.CI_high-r[m]]],fmt='o',color=col,capsize=3,lw=1.3,elinewidth=1.3,ms=5)
  plotted.append({'figure':'F4','metric':m,'target':t,'estimate':r[m],'low':r.CI_low,'high':r.CI_high})
 ax.axvline(0,color='#404040',lw=1.3);ax.set_yticks(range(3),NAMES,fontsize=9.5);ax.invert_yaxis();ax.set_xlim(*lim)
 ax.set_xlabel('Estimate and 95%\ndonor-bootstrap CI');label(ax,'a' if ax is axs[0] else 'b',title)
save(fig,'Figure4_Increment_Residual','Existing adjusted point estimates and fixed-OOF donor-block intervals, separate BH families, n=19.')
fig,axs=plt.subplots(1,3,figsize=(W,3.8),layout='constrained')
for j,(a,t,col,name) in enumerate(zip(axs,TS,COLORS,NAMES)):
 r=c[(c.target==t)&(c.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0]
 vals=[r.raw,r.macenko,r.grayscale,r.color_only,r.coordinates_only]
 a.plot(range(3),vals[:3],'o',color=col,ms=5);a.plot([3,4],vals[3:],'s',color=col,ms=5,fillstyle='none')
 a.axvline(2.5,color='#888888',lw=.8,linestyle='--');a.axhline(0,color='#606060',lw=1)
 a.set_xticks(range(5),['Raw','Macenko','Gray','Color','Coord'],rotation=45,ha='right',rotation_mode='anchor')
 baseline=p[(p.target==t)&(p.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0].M2_R2
 assert f'{baseline:.3f}'==['0.243','-0.444','0.155'][j]
 a.axhline(baseline,color='#444444',lw=1,linestyle=(0,(5,3)))
 a.text(.98,baseline,'M2 baseline',transform=a.get_yaxis_transform(),ha='right',va='bottom',fontsize=9.5)
 a.set_ylim(-.53,.53);a.set_xlim(-.4,4.4);a.set_ylabel('Donor-held-out R²');label(a,chr(97+j),name)
 a.text(1,-.39,'Augmented\nM4',transform=a.get_xaxis_transform(),ha='center',fontsize=9.5)
 a.text(3.5,-.39,'Standalone\ncontrols',transform=a.get_xaxis_transform(),ha='center',fontsize=9.5)
 plotted.append({'figure':'F5','target':t,'values':vals,'M2_baseline':baseline})
save(fig,'Figure5_Shortcut_Controls','Augmented raw/Macenko/gray M4 separated from standalone color/coordinate controls; no new uncertainty.')
(O/'figure_export_manifest.json').write_text(json.dumps(exports,indent=2))
(O/'figure_source_hashes.json').write_text(json.dumps(sources,indent=2))
(O/'plotted_frozen_values.json').write_text(json.dumps(plotted,indent=2))
print('Five figures exported as editable SVG and genuine 601 dpi PNG; no model or resampling executed.')
