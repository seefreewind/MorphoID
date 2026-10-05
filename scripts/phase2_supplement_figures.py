"""Presentation-only revisions of frozen registration and dependency evidence."""
from pathlib import Path
import json,hashlib
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
R=Path(__file__).resolve().parents[1];O=R/'manuscript/figures'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':7,'svg.fonttype':'none'})
# Original asset retained in place. Each original patch represents 150 micrometres;
# retain all panels and state calibration in the legend rather than infer image margins.
c=pd.read_csv(R/'configs/frozen_vannan_primary_cohort.tsv',sep='\t');c=c[c.primary_inclusion=='YES']
f,ax=plt.subplots(figsize=(183/25.4,7.2),layout='constrained');ax.axis('off')
ax.set_title('Frozen donor–section dependency structure',loc='left')
rows=[]
for donor,g in c.groupby('donor_id',sort=True):
 rows.append([donor,('PF' if g.disease_status.iloc[0]=='pulmonary_fibrosis' else 'Control'),', '.join(g.sample_id),', '.join(map(str,sorted(g.TMA.unique()))),', '.join(map(str,sorted(g.run.unique())))])
t=ax.table(cellText=rows,colLabels=['Donor','Disease','Admitted sections','TMA','Run'],loc='center',cellLoc='left',colWidths=[.18,.12,.46,.12,.12]);t.auto_set_font_size(False);t.set_fontsize(6);t.scale(1,1.42)
for (r,k),cell in t.get_celld().items():
 cell.set_edgecolor('#D1D5DB');cell.set_linewidth(.4)
 if r==0:cell.set_facecolor('#E5E7EB')
ax.text(0,.02,'19 independent donors → 26 sections → 12,144 common regions.\nSections/regions from the same donor remain in one outer fold. TMA/run can span donors.',transform=ax.transAxes,fontsize=7)
f.savefig(O/'S2_donor_dependency.svg');f.savefig(O/'S2_donor_dependency.png',dpi=300);plt.close(f)
s=json.loads((R/'results/final/figure_source_hashes.json').read_text());p='configs/frozen_vannan_primary_cohort.tsv';s[p]=hashlib.sha256((R/p).read_bytes()).hexdigest();(R/'results/final/figure_source_hashes.json').write_text(json.dumps(s,indent=2))
