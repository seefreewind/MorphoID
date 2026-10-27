"""Presentation-only redraw of all frozen rows for Figures 1–5; no fitting or inference."""
from pathlib import Path
import hashlib, json, shutil
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=Path('reproduced'))
args=parser.parse_args()
R=Path(__file__).resolve().parents[1]; O=args.output.resolve(); O.mkdir(parents=True,exist_ok=True)
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
 expected=json.loads((R/'manifests/frozen_display_sources.json').read_text())[f]
 assert sources[f]==expected, 'Frozen source checksum mismatch: '+f
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
 fig.savefig(O/(name+'.svg')); fig.savefig(O/(name+'.png'),dpi=601)
 text=[t.get_fontsize() for ax in fig.axes for t in ax.texts+[ax.title,ax.xaxis.label,ax.yaxis.label]+ax.get_xticklabels()+ax.get_yticklabels() if t.get_text()]
 exports.append({'figure':name,'width_mm':183,'height_mm':fig.get_size_inches()[1]*25.4,
                 'raster_dpi':601,'font_floor_pt':min(text),'editable':'SVG text elements',
                 'contract':contract,'new_estimates':False,'resampling':False})
 plt.close(fig)
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

(O/'plotted_frozen_values.json').write_text(json.dumps(plotted,indent=2))
(O/'figure_source_hashes.json').write_text(json.dumps(sources,indent=2))
(O/'figure_export_manifest.json').write_text(json.dumps(exports,indent=2))
rows=[]
for t,name in zip(TS,NAMES):
 a=p[(p.target==t)&(p.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0]
 b=q[(q.target==t)&(q.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0]
 pf=p[(p.target==t)&(p.analysis_stratum=='PF_ONLY')].iloc[0]
 d=c[(c.target==t)&(c.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0]
 rows.append(dict(target=t,n_regions=a['n_regions'],M2_R2=a.M2_R2,M4_R2=a.M4_R2,delta_R2=a.deltaR2_morphology,delta_CI_low=a.CI_low,delta_CI_high=a.CI_high,delta_q=a.FDR,residual_R2=b.residual_R2,residual_CI_low=b.CI_low,residual_CI_high=b.CI_high,PF_delta_R2=pf.deltaR2_morphology,Macenko_M4_R2=d.macenko,Grayscale_M4_R2=d.grayscale))
pd.DataFrame(rows).to_csv(O/'Table1_frozen_numeric.tsv',sep='\t',index=False,float_format='%.17g')
print('PASS: Figure 4, Figure 5 and Table 1 reproduced from checksum-verified frozen tables; no estimation or resampling.')
