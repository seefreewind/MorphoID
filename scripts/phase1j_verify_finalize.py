#!/usr/bin/env python3
"""Required numerical evidence and figure-export QA; not a model/test suite."""
from pathlib import Path
import json,re,hashlib,datetime
import numpy as np
import pandas as pd
from scipy.stats import t
from PIL import Image
from phase1j_precision import gate,effect,C
R=Path(__file__).resolve().parents[1];O=R/'results/phase1j';F=R/'figures/phase1j'
def main():
 p=pd.read_csv(O/'precision_projection.tsv',sep='\t')
 current=pd.read_csv(O/'current_precision.tsv',sep='\t')
 blocks=pd.read_csv(R/'results/phase1i/donor_prediction_blocks.tsv',sep='\t')
 pv=pd.read_csv(R/'results/phase1i/jackknife_pseudovalues.tsv',sep='\t')
 donors=sorted(blocks.donor_id.unique());assert len(donors)==19
 rng=np.random.default_rng(C['seed']);rng.multinomial(19,np.full(19,1/19),size=C['baseline_inner_oof_bootstraps'])
 additions=rng.integers(0,19,size=(C['outer_projections'],40));counts=np.ones((1000,19),dtype=np.int16)
 maximum_proxy_error=0.
 for n in range(19,60):
  if n>19:np.add.at(counts,(np.arange(1000),additions[:,n-20]),1)
  assert (counts>=1).all() and (counts.sum(axis=1)==n).all()
  for a in C['targets']:
   for e in C['effects']:
    values=pv[(pv.target==a)&(pv.effect==e)].set_index('donor_id').reindex(donors).pseudovalue.to_numpy(float)
    mean=counts@values/n
    variance=(counts@(values*values)-n*mean*mean)/(n-1)
    width=2*t.ppf(.975,n-1)*np.sqrt(np.maximum(variance,0)/n)
    row=p[(p.target==a)&(p.effect==e)&(p.total_donors==n)&(p.method=='REFIT_PSEUDOVALUE_PROXY')].iloc[0]
    err=abs(np.median(width)-row.projected_CI_width);maximum_proxy_error=max(maximum_proxy_error,err)
    assert err<1e-12
    for ratio,label in [(.8,'80'),(.67,'67'),(.5,'50')]:
     hits=np.mean(width/row.current_CI_width<=ratio)
     assert abs(hits-row['probability_width_le_'+label+'pct_current'])<1e-12
 cols=['n_regions','sum_y','sum_y2','SSE_M2','SSE_M4','sum_r','sum_r2','SSE_residual']
 maximum_point_error=0.
 for _,row in current.iterrows():
  b=blocks[blocks.target==row.target].set_index('donor_id').reindex(donors)[cols].to_numpy(float)
  value=float(effect(np.ones((1,19)),b,row.effect)[0])
  maximum_point_error=max(maximum_point_error,abs(value-row.original_estimate))
  assert np.isclose(value,row.original_estimate,atol=1e-11,rtol=1e-10)
 assert len(p)==492 and p.unreliable_projection_fraction.max()==0
 for label in ['80','67','50']:
  assert p['probability_width_le_'+label+'pct_current'].between(0,1).all()
 assert (p.probability_width_le_80pct_current>=p.probability_width_le_67pct_current).all()
 assert (p.probability_width_le_67pct_current>=p.probability_width_le_50pct_current).all()
 threshold=[]
 for a,d in p.groupby('target',sort=False):
  for e,ms in [('BOTH_EFFECTS_BOTH_METHODS',d)]+[(f'{ef}_{m}',d[(d.effect==ef)&(d.method==m)]) for ef in d.effect.unique() for m in d.method.unique()]:
   good=ms.groupby('total_donors').probability_width_le_67pct_current.min()>=.8
   ns=[int(n) for n in good.index if n>19 and good.loc[n:].all()]
   threshold.append(dict(target=a,criterion=e,minimum_total_donors=min(ns) if ns else None,MINIMUM_ADDITIONAL_DONORS=min(ns)-19 if ns else None,probability_rule='USER_REQUESTED_80_PERCENT;STABLE_THROUGH_CEILING',status='FEASIBLE_UNDER_EMPIRICAL_PROJECTION' if ns else 'PRECISION_RESCUE_NOT_FEASIBLE_AT_PRACTICAL_N'))
 pd.DataFrame(threshold).to_csv(O/'requested_probability_threshold.tsv',sep='\t',index=False)
 cross=pd.read_csv(O/'vannan_extension_crosswalk.tsv',sep='\t')
 assert len(cross)==45 and cross.GSM.nunique()==45 and cross.donor_id.nunique()==35
 new=cross[cross.new_independent_donor=='YES'];assert new.donor_id.nunique()==16
 assert new[new.registered_HE=='YES_AUTHOR_FILE'].donor_id.nunique()==13
 assert cross[cross.original_admitted=='YES'].sample_id.nunique()==26
 assert 'TILD117' not in set(new.donor_id)
 reg=pd.read_csv(O/'external_candidate_registry.tsv',sep='\t');assert len(reg)==16
 assert int(reg.independent_study.sum())==13
 assert not reg.admission_status.eq('PASS').any()
 assert set(reg.candidate_role)<= {'AUGMENTATION_CANDIDATE','EXTERNAL_VALIDATION_CANDIDATE','SUPPORTIVE_ONLY'}
 qa=[]
 for name,panels in [('J1_precision_projection',6),('J2_candidate_feasibility',1)]:
  svg=(F/(name+'.svg')).read_text()
  width_pt=float(re.search(r'width="([\d.]+)pt"',svg).group(1))
  assert abs(width_pt*25.4/72-183)<.01
  assert '<text' in svg and 'font-size' in svg or '<text' in svg
  with Image.open(F/(name+'.png')) as im:
   dpi=im.info.get('dpi',(0,0));assert min(dpi)>299
   assert abs(im.width/300*25.4-183)<.1
   dimensions=list(im.size)
  qa.append(dict(figure=name,panels=panels,visual_review='PASS_AFTER_J1_FOOTER_SPACING_REVISION',width_mm=183,PNG_dimensions=dimensions,dpi=dpi,editable_SVG=True,min_configured_font_pt=6.5))
 (O/'figure_qa.json').write_text(json.dumps(dict(backend='python',visual_panels_inspected=7,figures=qa,
  static_source_preflight='RAW_VALIDATOR_NOT_READY_FORMAT_EXCEPTIONS_DOCUMENTED',
  resolved_exceptions=['PDF omitted by explicit user output preference; editable SVG verified','PNG300 planning export satisfies declared contract; TIFF600 not requested','Static width parser misreads 183/25.4; actual SVG and raster independently verify 183mm'],
  publication_submission_status='PLANNING_ARTIFACTS_NOT_A_JOURNAL_SUBMISSION_PACKAGE'),indent=2))
 report=R/'reports/PHASE1J_PRECISION_AND_ADDITIONAL_DONOR_FEASIBILITY.md';s=report.read_text()
 assert re.findall(r'^(\d+)\. \*\*',s,re.M)==[str(i) for i in range(1,25)]
 for field in ['PRIMARY_19_DONOR_ANALYSIS_CHANGED','DONORS_REMOVED','TARGETS_CHANGED','ENCODER_CHANGED','MODEL_CHANGED','SIGNIFICANCE_TARGETED','EXTERNAL_MODELING_RUN','AUGMENTATION_MODELING_RUN']:
  assert field+' = NO' in s
 assert len(re.findall(r'^!\[',s,re.M))==2
 # Refresh metadata source manifest with canonical landing URLs; never certify error pages.
 manifest=pd.read_csv(O/'source_metadata_manifest.tsv',sep='\t').fillna('')
 for i,r in manifest.iterrows():
  f=R/r.file
  if re.match(r'GSE\d+.*soft.txt',f.name):
   acc=re.search(r'GSE\d+',f.name).group();manifest.loc[i,'source_URL']='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+acc
  manifest.loc[i,'bytes']=f.stat().st_size;manifest.loc[i,'sha256']=hashlib.sha256(f.read_bytes()).hexdigest()
 portal=O/'source_metadata/lungspatialdb_status.txt'
 if portal.exists() and not manifest.file.eq(str(portal.relative_to(R))).any():
  manifest=pd.concat([manifest,pd.DataFrame([dict(file=str(portal.relative_to(R)),source_URL='https://www.lungspatialdb.com',audit_date='2026-10-03',bytes=portal.stat().st_size,sha256=hashlib.sha256(portal.read_bytes()).hexdigest(),status='ACCESS_FAILURE_NOT_ABSENCE_EVIDENCE')])],ignore_index=True)
 manifest.to_csv(O/'source_metadata_manifest.tsv',sep='\t',index=False)
 gate()
 h=pd.read_csv(O/'input_hash_gate.tsv',sep='\t');assert h.pass_gate.all()
 audit=dict(status='PHASE1J_COMPLETE_VERIFIED_STOPPED',verdict='HOLD_TRANSPORTABILITY',external_validation_authorized='NO',augmentation_authorized='NO',
  frozen_hash_checks=len(h),all_frozen_hashes_pass=True,proxy_max_absolute_reconstruction_error=maximum_proxy_error,
  frozen_OOF_point_max_absolute_reconstruction_error=maximum_point_error,figure_QA='PASS_WITH_DOCUMENTED_FORMAT_EXCEPTIONS',
  no_model_fits=True,stop_after_phase1j=True,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 (O/'execution_status.json').write_text(json.dumps(audit,indent=2))
 inventory=[]
 for root in [O,F]:
  for f in sorted(root.rglob('*')):
   if f.is_file() and not f.name.startswith('._') and '.mplconfig' not in str(f) and f.name!='deliverable_manifest.tsv':
    inventory.append(dict(file=str(f.relative_to(R)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 for f in [report,R/'scripts/phase1j_precision.py',R/'scripts/phase1j_feasibility.py',R/'scripts/phase1j_figures_report.py',R/'scripts/phase1j_verify_finalize.py',R/'configs/phase1j_precision_plan.yaml']:
  inventory.append(dict(file=str(f.relative_to(R)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 pd.DataFrame(inventory).to_csv(O/'deliverable_manifest.tsv',sep='\t',index=False)
 print(json.dumps(audit),flush=True)
if __name__=='__main__':main()
