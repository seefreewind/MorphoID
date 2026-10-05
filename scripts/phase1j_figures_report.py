#!/usr/bin/env python3
"""Phase J planning figures and single evidence report; no biological model fit.
Figure contract: J1, six equal panels, precision planning curves, no effect signs.
J2, gate matrix, metadata support distinct from actual admission. 183 mm width.
Source: frozen donor contributions and public metadata. n=19 original donors,
1000 hypothetical augmentation draws, 500 nested donor bootstraps. Ranges are
10th–90th planning percentiles, not sampling confidence intervals. No tests,
P values or significance thresholds. Python-only. SVG + PNG300; PDF not requested.
"""
from pathlib import Path
import os,json,hashlib,sys
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1];O=R/'results/phase1j';F=R/'figures/phase1j'
F.mkdir(parents=True,exist_ok=True)
os.environ['MPLCONFIGDIR']=str(O/'.mplconfig')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':7,'axes.labelsize':7,'axes.titlesize':8,'xtick.labelsize':6.5,'ytick.labelsize':6.5,
 'legend.fontsize':6.5,'axes.linewidth':.6,'svg.fonttype':'none','pdf.fonttype':42,
 'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
T=['epithelial_injury','fibroblast_activation','macrophage_inflammatory']
NAMES=['Epithelial injury','Fibroblast activation','Macrophage inflammatory']
E=['DeltaR2_morph','Residual_R2'];METHOD=['OOF_BLOCK_EMPIRICAL','REFIT_PSEUDOVALUE_PROXY']
COLORS=['#0072B2','#D55E00']
def export(fig,name):
 fig.savefig(F/(name+'.svg'))
 fig.savefig(F/(name+'.png'),dpi=300)
 plt.close(fig)
def table(d):
 def fmt(v):
  if pd.isna(v):return 'NA'
  if isinstance(v,(float,np.floating)):return f'{v:.4g}'
  return str(v).replace('|',' / ').replace('\n',' ')
 cols=list(d.columns);s='| '+' | '.join(cols)+' |\n| '+' | '.join(['---']*len(cols))+' |\n'
 return s+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in d.itertuples(index=False,name=None))+'\n'
def main():
 p=pd.read_csv(O/'precision_projection.tsv',sep='\t')
 current=pd.read_csv(O/'current_precision.tsv',sep='\t')
 minima=pd.read_csv(O/'minimum_information_threshold.tsv',sep='\t')
 cross=pd.read_csv(O/'vannan_extension_crosswalk.tsv',sep='\t').fillna('')
 reg=pd.read_csv(O/'external_candidate_registry.tsv',sep='\t').fillna('UNKNOWN')
 trans=pd.read_csv(O/'transportability_matrix.tsv',sep='\t').fillna('UNKNOWN')
 summary=json.loads((O/'feasibility_complete.json').read_text())
 assert len(p)==492 and set(p.total_donors)==set(range(19,60))
 assert p.unreliable_projection_fraction.max()==0
 fig,axes=plt.subplots(2,3,figsize=(183/25.4,142/25.4))
 fig.subplots_adjust(left=.085,right=.98,bottom=.22,top=.87,wspace=.38,hspace=.5)
 for row,e in enumerate(E):
  for col,a in enumerate(T):
   ax=axes[row,col]
   for m,color in zip(METHOD,COLORS):
    d=p[(p.target==a)&(p.effect==e)&(p.method==m)].sort_values('total_donors')
    x=d.total_donors.to_numpy();y=d.projected_CI_width.to_numpy();lo=d.width_p10.to_numpy();hi=d.width_p90.to_numpy()
    assert np.all(np.isfinite(y)) and np.all(lo<=y) and np.all(y<=hi)
    # No smoothing or imposed monotonicity; points join observed integer-grid projections.
    ax.fill_between(x,lo,hi,color=color,alpha=.13,linewidth=0)
    ax.plot(x,y,color=color,lw=1.4)
    ax.plot(x,d.current_CI_width.to_numpy()*np.sqrt(19/x),color=color,lw=.7,ls=':')
   ax.set_xlim(19,59);ax.set_xticks([19,29,39,49,59]);ax.set_ylim(bottom=0)
   ax.set_title(NAMES[col],pad=9)
   ax.text(-.23,1.08,chr(97+row*3+col),transform=ax.transAxes,fontweight='bold',fontsize=8)
   ax.set_ylabel(('Delta R²' if row==0 else 'Residual R²')+' CI width')
   if row==1:ax.set_xlabel('Total independent donors (hypothetical)')
   ax.tick_params(width=.6,length=3)
 fig.suptitle('Precision projections retain current donor heterogeneity',fontsize=10,y=.975)
 handles=[Line2D([0],[0],color=COLORS[0],lw=1.4,label='Frozen OOF block projection'),
  Line2D([0],[0],color=COLORS[1],lw=1.4,label='Refit pseudovalue dispersion proxy'),
  Line2D([0],[0],color='#666666',ls=':',lw=.8,label='1 / sqrt(n) reference')]
 fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.025),ncol=2,frameon=False)
 fig.text(.5,.09,'Shading: 10th–90th percentiles across hypothetical cohorts; no new model fits',ha='center',fontsize=6.5)
 export(fig,'J1_precision_projection')
 # Metadata gates: + documented support, ? required evidence pending, x incompatible.
 # A metadata + is not a numerical registration PASS.
 tokens=[['+','?','+','?','+','?','?'],['+','+','?','?','x','x','x'],
  ['+','?','x','?','x','?','x'],['+','+','?','?','x','?','?'],
  ['+','x','?','?','+','x','x'],['+','+','?','x','+','?','x'],
  ['+','?','?','?','x','?','?'],['+','?','?','?','x','?','?'],
  ['+','+','?','?','x','?','?'],['+','+','?','?','x','?','?'],
  ['?','x','?','?','x','x','x'],['?','?','?','?','x','?','?'],
  ['+','?','?','?','x','?','?'],['+','?','?','?','x','?','x'],
  ['x','+','+','+','+','+','x'],['x','?','?','?','?','?','x']]
 assert len(tokens)==len(reg)
 labels=['Vannan TMA5: 16 new','CHUV NSCLC: 10','COPD Xenium: 38','CosMx COVID: 22',
  'Mayr IPF: 7 (study)','PF Visium: 8','Adult MERFISH: 6','NSCLC MERFISH: 46',
  'NeuMap HD: 8','EGFR HD: 4','Asthma: 8 confirmed*','AML Xenium: unknown',
  'CosMx COVID: 6','Vendor CosMx: 3','Vannan-derived: 0 new','Single-donor demo: 1']
 code={'+':0,'?':1,'x':2};matrix=np.array([[code[t] for t in row] for row in tokens])
 fig,ax=plt.subplots(figsize=(183/25.4,156/25.4))
 fig.subplots_adjust(left=.32,right=.98,top=.86,bottom=.15)
 ax.imshow(matrix,aspect='auto',cmap=ListedColormap(['#D5EAF2','#F1E8CE','#ECE0DE']),vmin=-.5,vmax=2.5)
 ax.set_xticks(range(7));ax.set_xticklabels(['Independent\ndonors','Same H&E\nsection','Frozen\ngenes','C / N\nstructure','Frozen\nZ structure','150 µm\nadmission','Overall\nadmission'])
 ax.xaxis.tick_top();ax.tick_params(axis='both',length=0,pad=6)
 ax.set_yticks(range(len(labels)));ax.set_yticklabels(labels,fontsize=7)
 for i,row in enumerate(tokens):
  for j,t in enumerate(row):ax.text(j,i,t,ha='center',va='center',fontsize=8,fontweight='bold')
 for s in ax.spines.values():s.set_visible(False)
 fig.suptitle('No candidate clears the complete frozen admission gate',fontsize=10,y=.975)
 fig.text(.33,.095,'+ Metadata supported     ? Pending / unverified     x Incompatible / excluded',fontsize=7)
 fig.text(.33,.06,'* Eight confirmed in cohort 1; 16 specimen IDs across two cohorts.\nVannan: 13 author-registered new donors; numerical admission remains pending.',fontsize=6.5)
 export(fig,'J2_candidate_feasibility')
 save=pd.DataFrame(tokens,columns=['Independent donors','Same HE section','Frozen genes','C/N structure','Z structure','150um admission','Overall admission'])
 save.insert(0,'candidate_id',reg.candidate_id);save.to_csv(O/'J2_figure_source.tsv',sep='\t',index=False)
 if '--figures-only' in sys.argv:
  print('PHASE1J_FIGURES_REVISED',flush=True)
  return
 # Table J1 preserves every requested grid and all three design probabilities.
 rows=[]
 for a in T:
  for e in E:
   for n in [19,24,29,34,39,49,59]:
    d=p[(p.target==a)&(p.effect==e)&(p.total_donors==n)].set_index('method')
    def pair(k):return ' / '.join(f'{d.loc[m,k]:.3f}' for m in METHOD)
    rows.append({'Target':a,'Effect':e,'Total / added':f'{n} / {n-19}',
     'Width OOF / proxy':pair('projected_CI_width'),'Relative width':pair('relative_width'),
     'Pr ≤80%':pair('probability_width_le_80pct_current'),'Pr ≤67%':pair('probability_width_le_67pct_current'),
     'Pr ≤50%':pair('probability_width_le_50pct_current')})
 J1=pd.DataFrame(rows);J1.to_csv(O/'Table_J1_display.tsv',sep='\t',index=False)
 def cw(a):
  d=current[current.target==a].set_index('effect')
  return '; '.join(e+f": reference {d.loc[e,'reference_CI_width']:.6f}, jackknife {d.loc[e,'jackknife_CI_width']:.6f}, BCa "+(f"{d.loc[e,'BCa_CI_width']:.6f}" if pd.notna(d.loc[e,'BCa_CI_width']) else 'UNRELIABLE / NOT USED') for e in E)
 now='2026-10-03'
 report=R/'reports/PHASE1J_PRECISION_AND_ADDITIONAL_DONOR_FEASIBILITY.md'
 assert not report.exists(), 'Inspect existing final report before any revision'
 lines=[
 '1. **PHASE 1J VERDICT:** HOLD_TRANSPORTABILITY',
 '2. **Current independent donors:** 19',
 '3. **Primary cohort changed:** NO',
 '4. **Primary point estimates changed:** NO',
 '5. **Significance used as planning target:** NO',
 '6. **Current epithelial CI width:** '+cw(T[0]),
 '7. **Current fibroblast CI width:** '+cw(T[1]),
 '8. **Current macrophage CI width:** '+cw(T[2]),
 '9. **Added donors for 33% precision gain — epithelial:** 29 at the requested ≥80% probability threshold (total 48); 32 with the conservative Monte Carlo safeguard (total 51)',
 '10. **Added donors for 33% precision gain — fibroblast:** PRECISION_RESCUE_NOT_FEASIBLE_AT_PRACTICAL_N (evaluated through 40 added / 59 total)',
 '11. **Added donors for 33% precision gain — macrophage:** 29 at the requested ≥80% probability threshold (total 48); 30 with the conservative Monte Carlo safeguard (total 49)',
 '12. **Vannan extension new independent donors:** 16',
 '13. **Vannan extension technically admissible donors:** 0 verified under frozen admission gates; 13 have author-registered per-core H&E metadata',
 '14. **Independent external cohorts identified:** 13 study sources screened; 0 admitted independent validation cohorts',
 '15. **Best external candidate:** NONE ADMITTED. Same-section image-access follow-up priority: Columbia CosMx COVID22; suitable only for supportive lung replication unless unchanged estimand transport is established',
 '16. **Same-section registration feasible:** NO — complete frozen admission not established for a new cohort; same-slide acquisition is documented for some candidates',
 '17. **Frozen targets transportable:** NO — no external candidate has cleared all unchanged gene, measurement-unit and lineage gates',
 '18. **Composition oracle transportable:** NO — no candidate has cleared the complete frozen oracle/unit/ontology gate',
 '19. **150 μm morphology transportable:** NO — new-donor numerical coordinate, tissue and full-crop admission not established',
 '20. **Augmentation feasible:** NO — strict admission remains pending; a 13-donor metadata-supported source exists',
 '21. **Independent validation feasible:** NO',
 '22. **External validation authorized:** NO',
 '23. **Augmentation authorized:** NO',
 '24. **Recommended next phase:** Await explicit authorization for a bounded Vannan TMA5 admission audit of the 13 author-registered new donors; original discovery remains frozen. Do not start augmentation or validation models.',
 '', '# Phase 1J precision and additional-donor feasibility', '',
 f'Audit date: {now}. Scope: planning and public metadata only. The frozen Phase 1I verdict remains `HOLD_FINITE_DONOR_PRECISION`. This audit neither changes that inference nor creates an external effect estimate.',
 '', '## Decision', '',
 'There are plausible additional donor sources, so project closure for lack of any source is premature. Their availability does not establish admission. The closest augmentation source is Vannan TMA5: 16 donors outside the original 19, including 13 with author-level per-core registered H&E. Coordinate validity, tissue coverage, 150 μm crop support and the nuclear transcript unit have not been cleared for those new donors. No independent cohort currently passes the complete frozen estimand gate.',
 '', 'The empirical projections indicate that 13–16 potential extension donors would not meet the conservative 33% width-reduction target for any of the three targets: the minima are 32 and 30 for epithelial and macrophage, and unattained for fibroblast through 40 added donors. This does not mean the extension would add no information; it means the prespecified practical precision target is not established by that source alone.',
 '', '## Frozen design and planning contract', '',
 'The primary analysis retains 19 donors, 26 admitted sections, the common universe, ≥20 lineage cells per eligible region, the three signatures, 150 μm morphology, C/N/Z, Phikon-v2, Ridge model family, original donor folds, OOF predictions and all Phase 1H/1I results. New observations would be labeled AUGMENTATION or INDEPENDENT_VALIDATION and could never replace the original discovery results.',
 '', 'The planning contract was written to `configs/phase1j_precision_plan.yaml` before simulation and captured with SHA256 in `results/phase1j/pre_run_contract.json`. Both start and end checks compare the frozen scientific inputs, Phase 1I outputs and reports; the final gate is recorded in `input_hash_gate.tsv`. Model training and candidate image/count processing were not run.',
 '', '### Current precision and donor dispersion', '',
 table(current[['target','effect','original_estimate','reference_CI_width','jackknife_CI_width','BCa_CI_width','pseudovalue_SD','top2_influence_fraction']]),
 '', 'Fibroblast has the greatest relative precision need by the preregistered maximum jackknife/reference-width ratio (approximately 2.82 for residual R²); its Delta R² also has the largest absolute width. The original target family is unchanged. The unreliable fibroblast residual BCa interval remains unavailable. The top-two influence fraction and pseudovalue dispersion are descriptive summaries; no donor was removed.',
 '', '### Heterogeneity-preserving projections', '',
 'Each hypothetical cohort retains all 19 original donors exactly once and appends 0–40 draws from their empirical joint donor distribution. The same donor identities are used across all targets and both estimands. Repeated historical blocks are a planning approximation to future independent observations; they are not additional real donors and do not increase the actual observed donor count.',
 '', 'For the OOF projection, 1,000 hypothetical augmentation compositions are evaluated at every integer total n=19–59. Each composition receives 500 whole-donor bootstrap draws; donor sums reconstruct the original region-weighted Delta R² and residual R² algebra. These intervals are conditional on frozen predictions and do not include future model-training uncertainty. A 10,000-draw n=19 baseline supplies a multiplicative calibration to the existing frozen reference width. Calibration factors (0.955–1.007) and simulated baseline widths are retained in `current_precision.tsv`; no observed effect is recentered.',
 '', 'The second method appends donor-paired Phase 1I refit jackknife pseudovalues and projects width as 2 × t(0.975,n−1) × sample SD / √n. This exactly reproduces the original jackknife width at n=19. It is a dispersion proxy for refitted-estimator uncertainty, not a new jackknife refit or a prediction of future training stability. The dotted 1/√n scaling is a reference only.',
 '', 'P1, P2 and P3 require width ≤80%, ≤67% and ≤50% of each method’s own frozen current width. Probabilities use all 1,000 outer projections; invalid projections cannot count as successes. Nonpositive target variance invalidates a bootstrap draw and >5% invalid inner draws invalidates its projection. No unreliable projections were observed. The two methods answer different precision questions and their widths are not pooled.',
 '', 'The requested threshold of ≥80% probability for P2 is reported separately by effect and method. The conservative target threshold also requires the Wilson 95% Monte Carlo lower bound ≥0.80 for both effects and both methods, with the criterion retained at every larger evaluated n. This additional safeguard was specified before simulation. It measures numerical planning reliability, not biological certainty. No criterion uses p values, FDR, effect direction or whether an interval crosses zero.',
 '', table(minima),
 '', 'At 59 total donors, fibroblast Delta R² meets P2 in only 63.4% of OOF projections (Wilson lower bound 60.4%), despite a median relative width of 0.639. The simpler square-root rule would conceal this tail risk. `PRECISION_LIMIT_NOT_PRACTICALLY_RESCUABLE` applies to its conservative family-level rescue under this empirical planning contract, not as proof that all future data could never improve precision. No n>59 was simulated.',
 '', '## Table J1. Precision projections', '',
 'Every paired cell is OOF / pseudovalue proxy. Widths are medians across hypothetical cohorts; probabilities are the proportions meeting each reduction target. The complete continuous-n table includes planning percentiles, Wilson bounds and invalid-draw diagnostics in `results/phase1j/precision_projection.tsv`.',
 '', table(J1), '', '## Figure J1. CI width versus hypothetical total donor count', '',
 '![J1 precision planning](../figures/phase1j/J1_precision_projection.png)', '',
 '**Legend.** Columns represent the three retained targets; upper and lower rows show Delta R² and residual R² interval widths. Blue curves are calibrated frozen OOF block projections; orange curves are refit-pseudovalue dispersion proxies. Shading is the 10th–90th planning percentile across 1,000 hypothetical augmentation compositions, not a sampling confidence band. Dotted curves show the corresponding 1/√n reference. The original 19 donors remain fixed; appended donors are sampled from their joint empirical distribution. No biological model was fitted and no statistical significance test was used. Source: `precision_projection.tsv`.',
 '', '## Vannan extension provenance', '',
 'The current GSE250346 record contains 45 sample records and 35 donor IDs, matching the historical complete crosswalk. The frozen discovery subset was drawn from the original 28 cores/19 donors, of which 26 sections passed admission. The 17 TMA5 records contain 16 new donors plus another core of the original donor TILD117. TMA5 is already part of the final 2025 publication; it is an extension of the frozen preprint-era subset, not a newly discovered post-final-publication cohort. No later additional GSM was found in this audit. [Current GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346), [final author article](https://www.nature.com/articles/s41588-025-02080-x).',
 '', 'The shared author panel is 343 genes and contains all 18 frozen signature genes. TMA5 uses Xenium 2.0.0.10 and multi-tissue staining, versus 1.1.0.2 in the original TMAs. Author annotations for all 17 TMA5 sample IDs already exist locally. The author analysis restricts cell-aware expression to nuclear transcripts; extension admission must preserve that unit instead of silently substituting updated full-cell counts. Public registered per-core H&E is absent for TILD028LA, TILD167LA and VUHD049. Whole-TMA or alternate-modality images are not counted as an admitted per-core replacement.',
 '', 'Historical corrections remain provenance events: the August 2024 codeword correction affects ACTA2/MRC1; the corrected annotation/expression source is retained. The May 2025 TILD117 rename separates original MA1 (GSM7990534) from TMA5 MA2 (GSM8505453); it does not create a donor. The June 2025 processed-data update and current February 2026 GEO update are tracked without changing the primary inputs. Same-section author provenance is supported; numerical registration PASS has not been inferred from filenames.',
 '', '### Table J2. TMA5 donor audit', '',
 table(cross[cross.TMA==5][['sample_id','donor_id','GSM','new_independent_donor','registered_HE','coordinate_validity','morphology_150um']]),
 '', '`VANNAN_EXTENSION_INSUFFICIENT` under the complete technical-admission rule: 16 new donor IDs, 13 metadata-supported registered partners, 0 newly technically admitted. This verdict records unresolved admission evidence, not a finding that 16 donors are numerically below five. The 45-row crosswalk retains all requested fields, correction notes, source URLs and original admission membership.',
 '', '## Public-cohort search and independence', '',
 'A bounded public search through 2026-10-03 screened Xenium first, followed by CosMx, MERFISH, Visium HD and conventional Visium. Thirteen independent study sources were examined, alongside Vannan extension, derived repositories and a single-donor platform contingency. Unknown specimen duplication, donor crosswalk, gene coverage or registration is kept unknown. Repeated cores, modalities, visits and panel releases are never counted as independent donors. The registry is a screening record, not a declaration that all listed studies are compatible. No exhaustive-absence claim is made.',
 '', 'The most direct independent same-slide acquisition lead is the 22-donor Columbia CosMx study. Its protocol explicitly stains the imaged slides with H&E after CosMx. Public raw spatial files, registration transforms, exact frozen gene coverage and FOV support remain unverified. Its GEO superseries points to sequence experiments; GSE287111 is scRNA-seq and must not be mislabeled as the CosMx donor dataset. COVID/control is not the frozen PF/control adjustment. [Author spatial study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12463091/).',
 '', 'The COPD custom panel was checked directly from all four small public panel JSON files: only 10 of 18 frozen genes are present. Missing genes are '+', '.join(summary['COPD_missing_frozen_genes'])+'. All three target definitions lose genes. They are `TARGET_NOT_TRANSPORTABLE`; no subset signature or proxy gene can rescue admission. Gene-level evidence is in `frozen_gene_coverage.tsv`. Linked participant metadata are controlled through BioLINCC; fluorescence morphology files are not H&E. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313006), [author article](https://www.nature.com/articles/s41588-025-02480-z).',
 '', '### Table J3. Candidate registry', '',
 table(reg[['candidate_id','accession','n_donors','same_section','candidate_role','admission_status','reason','source_URL']]),
 '', 'Roles are prospective suitability classes. Candidates missing a required gate have no permission for modeling. SUPPORTIVE_ONLY includes partial tissue/platform replication that cannot apply the complete frozen estimand. No independent study currently qualifies for EXTERNAL_VALIDATION_CANDIDATE admission. Existing CHUV registration failures are reused from the frozen Phase 0C-S evidence, without rerunning transform optimization. Asthma has eight confirmed cohort-1 donors and 16 specimen IDs across historical cohorts; independent cross-cohort duplication is not certified. Its H&E uses additional sections, which fails strict same-section admission.',
 '', '## Table J4. Frozen-estimand transportability', '',
 table(trans[['candidate_id','target_definitions','C_ontology','N_structure','Z_structure','morphology_150um','encoder_preprocessing','admission_status']]),
 '', '“NO” in the opening transport fields means no complete new-donor PASS has been established. Pending evidence is distinguished from biological incompatibility in Table J4. Commercial availability of high-plex genes is not exact gene coverage. Cell labels must map to the frozen coarse C ontology and lineage units; deconvolved spot mixtures or image-inferred labels are not automatically true composition oracles. N requires the original physical neighborhood, Z retains normalized x/y and PF/control status, and H&E support must permit full 150 μm crops with unchanged encoder preprocessing. DAPI registration or RNA/protein panel coregistration cannot substitute for H&E admission.',
 '', '## Figure J2. Candidate feasibility gate map', '',
 '![J2 candidate feasibility](../figures/phase1j/J2_candidate_feasibility.png)', '',
 '**Legend.** Rows represent screened study sources, the Vannan augmentation source and excluded duplicate/single-donor contingencies. Columns summarize independent donor evidence, same H&E section, frozen genes, C/N, frozen Z, 150 μm admission and overall admission. “+” denotes metadata support, “?” unresolved evidence, and “x” incompatibility or exclusion. Metadata support does not establish a numerical registration PASS. No biological effects or predicted outcomes enter the map. Asthma’s eight confirmed cohort-1 donors and 16 historical specimen IDs are distinguished. Source: `J2_figure_source.tsv` and the full candidate registry.',
 '', '## Next-stage options — defined, not executed', '',
 '**Option A: independent validation.** Admission would require ≥5 independently identified donors from an external cohort, unchanged signatures and measurement units, supported C/N/Z, same-section H&E registration PASS and valid 150 μm morphology. Keep the same encoder, scale, Ridge family and donor-held-out logic. Define an external cohort’s own frozen donor split without rewriting the original split. Fit/evaluate within that separate cohort only after explicit authorization; report positive, null and negative effects equally. No candidate is selected for its anticipated fibroblast result.',
 '', '**Option B: augmentation.** If the 13 author-registered TMA5 donors clear admission, retain the original 19 as the primary discovery cohort and identify new donors as a prospective augmentation set. Prespecify: original discovery fixed, extension-only estimation, then combined sensitivity. Any combined output remains a sensitivity result and cannot replace the frozen discovery result. Annotation/unit/geometry admission comes before model execution. The 16-donor count cannot be treated as an assured solution to the 30–32 donor precision needs.',
 '', '**Option C: no further compatible data.** If the bounded admission audit cannot establish new compatible donors, retain FINITE_DONOR_PRECISION_LIMITATION and move to manuscript interpretation only with authorization. Closure is not triggered here because a plausible augmentation source remains. No larger simulation, new target or encoder comparison is proposed.',
 '', '### Future cohort-level synthesis', '',
 'Only after independently authorized compatible cohorts exist: estimate Delta R² and residual R² separately per cohort, then perform fixed-effect descriptive synthesis with explicit estimator and variance compatibility. Add random-effects sensitivity only with ≥3 cohorts; with two, do not give heterogeneity variance a strong interpretation. Cohort directions must all be reported. No meta-analysis, concatenation, batch correction or embedding harmonization was performed in Phase 1J.',
 '', '## Evidence, QA and limitations', '',
 'Source metadata are retained under `results/phase1j/source_metadata/` with URL, bytes, SHA256 and response status in `source_metadata_manifest.tsv`. These are article XML, GEO descriptive/platform records, small panel files and one participant/specimen metadata spreadsheet, not full expression matrices or H&E datasets. Large GEO platform SOFT records are metadata; their size does not imply a downloaded candidate count matrix. Unsupported GEO `targ=samples` responses are explicitly flagged and not used as evidence. Public portal access failure is an access limitation, not proof that no files exist.',
 '', 'The projections assume future donor contributions resemble the present 19-donor empirical distribution. They cannot identify unseen disease/site heterogeneity, future refit instability, measurement shifts or the information gain of genuinely new donors outside that support. OOF and pseudovalue methods preserve the observed heterogeneity in different ways and must not be read as guaranteed future CIs. The strict target minimum combines two effects and two methods and is more conservative than a single-method threshold.',
 '', 'Figure QA uses the nature-figure Python workflow: source preflight, editable SVG, 300-dpi PNG and complete six-panel/matrix visual inspection at the declared 183-mm size. No PDF was requested. The raw static validator remains NOT READY under its default submission bundle: it requires PDF/TIFF and its 600-dpi default, and misreads the 183/25.4 width expression. These are documented planning-export exceptions in figure_qa.json; actual editable SVG width (183 mm), raster size and 300-dpi resolution are independently verified. Visual and numerical checks are separate. The final completion record is written only after visual inspection and the final frozen-input hash gate.',
 '', '## Required assertions', '',
 '```text',
 'PRIMARY_19_DONOR_ANALYSIS_CHANGED = NO',
 'DONORS_REMOVED = NO',
 'TARGETS_CHANGED = NO',
 'ENCODER_CHANGED = NO',
 'MODEL_CHANGED = NO',
 'SIGNIFICANCE_TARGETED = NO',
 'EXTERNAL_MODELING_RUN = NO',
 'AUGMENTATION_MODELING_RUN = NO',
 'EXTERNAL_VALIDATION_AUTHORIZED = NO',
 'AUGMENTATION_AUTHORIZED = NO',
 '```', '',
 'Phase 1J stops here. No candidate full dataset download, new embedding extraction, predictive refit, external Delta R² calculation or biological signal validation is authorized by this report.'
 ]
 report.write_text('\n'.join(lines)+'\n')
 print('PHASE1J_FIGURES_REPORT_CREATED',report,flush=True)
if __name__=='__main__':main()
