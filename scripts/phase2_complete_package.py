"""Integrate frozen tables into manuscript documents. No model execution."""
from pathlib import Path
import pandas as pd,json,hashlib,re,sys,platform,subprocess,importlib.metadata as md
R=Path(__file__).resolve().parents[1];M=R/'manuscript';F=R/'results/final'
def read(p):return pd.read_csv(R/p,sep='\t')
def table(d):
 def val(x):
  if pd.isna(x):return 'NA'
  if isinstance(x,float):return f'{x:.5g}'
  return str(x).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(map(str,d.columns))+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(val(x) for x in row)+' |' for row in d.itertuples(index=False,name=None))+'\n'
def write(p,s):(R/p).write_text(s,encoding='utf-8')
p=read('results/phase1/tables/table2_primary_model_performance.tsv');q=read('results/phase1/tables/table3_crossfit_residual_performance.tsv');e=read('results/final/target_evidence_classification.tsv')
rows=[]
for _,r in p[p.analysis_stratum=='DISEASE_ADJUSTED'].iterrows():
 t=r.target;rr=q[(q.target==t)&(q.analysis_stratum=='DISEASE_ADJUSTED')].iloc[0];pf=p[(p.target==t)&(p.analysis_stratum=='PF_ONLY')].iloc[0]
 rows.append([t,int(r.n_regions),r.M2_R2,r.M4_R2,r.deltaR2_morphology,f'[{r.CI_low:.4f}, {r.CI_high:.4f}]',r.FDR,rr.residual_R2,f'[{rr.CI_low:.4f}, {rr.CI_high:.4f}]',pf.deltaR2_morphology,r.macenko_R2,r.grayscale_R2,e[e.target==t].FINAL_EVIDENCE_LEVEL.iloc[0]])
t2=table(pd.DataFrame(rows,columns=['Target','Eligible regions','M2 R²','M4 R²','ΔR²','ΔR² 95% CI','ΔR² q','Residual R²','Residual 95% CI','PF-only ΔR²','Macenko M4 R²','Gray M4 R²','Class']))
t1=table(pd.DataFrame([['Dataset','GSE250346; author-registered paired H&E–Xenium'],['Independent biological units','19 donors: 6 control, 13 PF'],['Sections','26: 7 control, 19 PF'],['Morphology unit','150 × 150 μm'],['Common analysis universe','12,144 regions'],['Pre-model exclusions','77 coordinate-support; 39 zero valid image pixels; no overlap'],['Targets','3 within-lineage programs; target-specific eligibility'],['Primary comparison','Disease-adjusted M4 versus M2; frozen donor outer folds'],['Sensitivity strata','Pooled; PF-only independently retrained on 13 donors'],['Frozen encoder','Phikon-v2, 1,024-dimensional CLS'],['Validation boundary','Donor-held-out within one cohort; no independent cohort']],columns=['Item','Frozen structure']))
t3=table(pd.DataFrame([
 ['Epithelial injury','Pooled C increment 0.21624','Adjusted Δ positive, CI spans zero','Adjusted negative; partition support mixed','Macenko/gray adjusted increments not positive','98% pooled; 96% PF Δ positive','Unstable; donor deletion can reverse direction','C: limited positive tendency; residual unsupported'],
 ['Fibroblast activation','Pooled C increment 0.33708; adjusted M2 negative','Positive across partitions and adjusted omissions','Small adjusted value; 4 LOO residual reversals; PF positive','Positive augmented increments in Macenko/gray; no causal attribution','100% pooled/PF Δ positive','Magnitude, CI/FDR sensitive; residual BCa invalid','B for incremental direction; selective identifiability not established'],
 ['Macrophage inflammatory','Pooled C increment 0.27078','Adjusted Δ negative; all adjusted omissions negative','Adjusted and PF residual negative','Selected preprocessing controls do not restore increment','92% pooled; 90% PF Δ negative','Negative but inference-sensitive','D under this cohort/panel/encoder/150 μm design']
],columns=['Target','Composition contribution','Morphology increment','Residual support','Shortcut sensitivity','Partition direction','Inference robustness','Conclusion']))
f=M/'MorphoID_Main_Manuscript.md';s=f.read_text();refs=json.loads((F/'literature_metadata.json').read_text());bib='\n'.join(f"{r['id']}. {r.get('authors','Primary resource')}. {r['title']} *{r.get('journal','Primary resource')}* ({r.get('year','n.d.')}). [{r.get('doi') or 'Primary source'}]({r['url']})." for r in refs)
s=s.replace('{{TABLE1}}','**Table 1. Cohort and analysis structure.**\n\n'+t1).replace('{{TABLE2}}','**Table 2. Frozen primary results.**\n\n'+t2+'\nPrimary intervals are 2,000 fixed-OOF donor-block percentile intervals. q values belong to the three-target adjusted morphology-increment family. PF-only is a separate sensitivity fit. Macenko/gray columns show augmented M4 performance, not isolated morphology performance. Full precision, residual q values and all strata are in Supplementary Tables.').replace('{{TABLE3}}','**Table 3. Target-specific interpretation matrix.**\n\n'+t3).replace('{{REFERENCES}}',bib)
s=s.replace('expression benchmarks and tissue-level inference frameworks','expression benchmarks, communication modeling and tissue-level inference frameworks')
f.write_text(s)
# Chronology and methods are source-grounded; complete earlier reports remain intact.
sm='''# MorphoID Supplementary Methods

## Scope and evidence freeze

This package integrates completed evidence. No model, donor, section, region, target, fold or uncertainty method was added in Phase 2. The scientific configurations and prior results retain their byte hashes. The current administrative project status is STOP_FINITE_DONOR_LIMITATION_ACCEPTED. Primary discovery analysis is frozen; finite-donor precision limitations are accepted; external validation and augmentation are unavailable and unauthorized; no computational rescue is authorized.

## Registration reconstruction and raw–registered provenance

The discovery dataset is GSE250346. The audited H&E images were author-supplied post-Xenium images of the same section, with registration reconstructed from deposited images, coordinate support and author methodology. The author procedure used Aperio CS2 20× imaging, ImageJ/BigWarp and approximately 200 landmarks with thin-plate-spline transformation. Original manual landmarks and a complete raw-core crosswalk were unavailable. The author-registered TIFF is not a new native scanner image generated in this project. Qualification applies to 150 μm microregions; single-cell or nuclear correspondence was not validated. The provenance and registration reconstruction reports preserve unresolved measurement limits.

The original L8 candidate remained on HOLD_AUTHOR_ALIGNMENT_REQUIRED. The official pancreas positive control was a registration control, not a lung discovery donor. The replacement lung cohort was admitted after the documented provenance and coordinate-support audit. Of 28 audited sections, VUILD105MA1 and VUILD48LA1 were permanently excluded for fixed coordinate-support failures. The 26 admitted sections represent 19 donors. Section aliases in extraction logs must be resolved with the frozen cohort table rather than interpreted as extra samples.

## D5 missing-color audit and common-universe amendment

The frozen D5 table contained 39 regions with all 40 color features missing despite matching merge keys and nominal full-crop coordinate support. The repair audit established zero valid image pixels after the fixed black-pixel rule for these regions. The D5 source was preserved; no imputation or outcome-dependent exclusion was performed. A separately authorized pre-model common-universe amendment removed these 39 regions and 77 full-crop coordinate failures from all model comparisons, with no overlap. Thus 12,260 original regions yielded 12,144 common supported regions. The image-pixel audit used black threshold 0.02, stride 8, grid offset 0 and image scale 0.2125. Target eligibility then selected 3,219 epithelial, 4,869 fibroblast and 2,661 macrophage regions in the pooled/adjusted analyses. These target sets overlap conceptually and are not summed as an independent sample count.

## Molecular programs and measured covariates

The frozen within-lineage target and candidate-program tables specify gene membership, cell-lineage assignment and eligibility rules. Their complete contents and panel coverage are included in Supplementary Tables. Composition comprises 15 measured cell lineages. Neighborhood features summarize measured composition in the four cardinal neighboring grid regions; missing neighbors contribute zero according to the frozen construction. Covariates contain normalized coordinates, and the adjusted stratum also contains the PF indicator. Outcome definitions, panel coverage and eligibility thresholds were not revised during manuscript integration.

## Embeddings and execution provenance

Phikon-v2 is fixed to official Owkin revision 2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd and checkpoint SHA256 261ae680fa699b3b951597fd57aa19c02ef735805acb104b93af69b36d928569. The representation is the 1,024-dimensional CLS vector after the frozen processor. The extraction used CPU backend and batch size 4. All 26 section caches were preserved. Extraction covered 12,183 regions before application of the amended common universe; this is not the final modeling population. Whole-slide stain sampling uses 200 equally spaced block starts × 100 contiguous RGB pixels; the fallback is 1,000 × 100 when fewer than 500 pixels remain after optical-density filtering. TIFF decompression followed by sequential RAM copying changed I/O only, preserving image layout, uint8 values and crop coordinates.

## Model ladder, tuning and conditional residual estimand

M0=Z; M1=C+Z; M2=C+N+Z; M3=X+Z; M4=C+N+Z+X. ΔR² is M4 minus M2 held-out R². Ridge tuning uses 19 logarithmically spaced penalties from 10⁻⁴ to 10⁵, with training-donor GroupKFold inner splits up to five folds. Transformations and tuning are training-side. Inverse donor region counts provide equal total donor weight during fitting/tuning. Reported primary R² is the region-weighted pooled out-of-fold score; training weights do not convert it into a donor-average estimand.

For each outer test fold, M2 training-side predictions are generated by donor cross-fitting within the outer training set, retaining the frozen fold grouping and nested tuning. The X-only ridge model fits these training residuals. Test residuals are measured outcomes minus the outer-training M2 test predictions. Test outcomes are never used to fit a residualizer. This defines a model-conditioned predictive residual; it is neither a causal effect nor an information-theoretic proof of biological independence.

## Strata and shortcut controls

The primary inferential stratum is disease-adjusted. Pooled and PF-only estimates are separate sensitivity analyses. PF-only retrains on 13 PF donors and 19 sections rather than subsetting pooled predictions. Raw, training-donor-only Macenko and grayscale embeddings enter augmented M4; standalone color and coordinate controls are distinct models. Grayscale uses rounded uint8 0.2126R+0.7152G+0.0722B replicated to three channels. Metadata shortcuts and molecular outcome controls answer different questions. Phase 0D diagnostics preceded the common-universe amendment and must not be described as newly refitted on 12,144 regions.

## Uncertainty and partition sensitivity

The primary 2,000-replicate donor-block bootstrap resamples fixed out-of-fold predictions; it does not refit the model. Percentile intervals and two-sided bootstrap sign-tail p values are reported, with BH correction separately across three adjusted ΔR² tests and three adjusted residual tests. The sign-tail calculation is not a fully null-calibrated permutation test.

Phase 1H evaluates 50 frozen donor partitions for pooled and PF-only strata only. These reuse the same donors and are dependent split sensitivities, not 50 biological replications. No 50-partition disease-adjusted series was executed. Actual refit donor omissions total 19 pooled, 19 adjusted and 13 PF-only per target; the six control-donor PF-only no-op rows per target are retained in the original log but excluded from the displayed refit series. The support-distance/error diagnostics are descriptive.

Phase 1I retains the reference bootstrap, delete-one-donor refit jackknife normal/t intervals, balanced donor bootstrap with 10,000 draws, and BCa sensitivity where valid. Jackknife pseudo-values have their own mean and t₁₈ interval, which need not equal the original estimate. BCa acceleration uses matched delete-one fixed-OOF blocks; the fibroblast residual BCa result was unreliable and is not replaced by a favorable interval. Exact sign and Wilcoxon sensitivities concern donor-level loss gains; their p values are not R² effect tests. All 48 method/effect rows, including invalid entries, remain available.

## Historical HOLD transitions and closure

Phase 0 registration HOLD/PASS histories and Phase 0D confounding qualification remain preserved. The primary pilot concluded HOLD_MODEL_INSTABILITY; Phase 1H concluded HOLD_INFERENCE_INSTABILITY; Phase 1I concluded HOLD_FINITE_DONOR_PRECISION. Phase 1J assessed additional-donor feasibility and frozen precision projections without admitting an external or augmentation cohort. Its transportability HOLD is not overwritten. The administrative closure accepts these limits, rather than converting historical HOLD verdicts into scientific PASS results.

The hypothetical precision curves are already-completed Phase 1J projections conditional on current donor effects, not forecasts of successfully recruited donors or external replication. No new data acquisition, modeling, projection or resampling was performed for this package.

## Reproduction and figure interpretation

All frozen table hashes, configuration hashes, result registry and figure sources are indexed in FINAL_REPRODUCIBILITY_MANIFEST.md. Reproduction order is historical QC → approved common-universe amendment → extraction → primary models → Phase 1H → Phase 1I → Phase 1J → Phase 2 presentation. This describes provenance, not an instruction to rerun the closed project. The checked-out commit alone does not capture uncommitted scientific files; the hash contract is the controlling snapshot. Current observed software is recorded separately from the declared environment. Original images, tables, caches and negative estimates remain intact.
'''
write('manuscript/MorphoID_Supplementary_Methods.md',sm)
# Full evidence tables, not selection by favorable inference.
items=[('S1 Frozen cohort','configs/frozen_vannan_primary_cohort.tsv'),('S2 Targets','configs/frozen_within_lineage_targets.tsv'),('S3 Candidate programs','configs/frozen_phase1_candidate_programs.tsv'),('S4 Panel coverage','results/phase1/target_gene_coverage.tsv'),('S5 All primary results','results/phase1/tables/table2_primary_model_performance.tsv'),('S6 All residual results','results/phase1/tables/table3_crossfit_residual_performance.tsv'),('S7 Shortcut controls','results/phase1/tables/table4_shortcut_stress_test.tsv'),('S8 All 50 partitions','results/phase1h/repeated_fold_metrics.tsv'),('S9 Partition summary','results/phase1h/repeated_fold_stability.tsv'),('S10 All donor omissions including labeled PF no-ops','results/phase1h/leave_one_donor_out_influence.tsv'),('S11 All uncertainty methods','results/phase1i/inference_agreement_matrix.tsv')]
parts=['# MorphoID Supplementary Tables\n\nFull frozen rows are reproduced below. NA retains undefined/not-applicable values; no negative estimates are censored. Source TSVs provide machine-readable full precision. Registry contains 10,981 evidence records, with source hashes. PF-only control-donor no-op rows are present in S10 for transparency and are not counted as actual refits.\n']
for title,path in items:
 if not (R/path).exists():raise FileNotFoundError(path)
 parts.append('## '+title+'\n\nSource: ['+path+']('+str(R/path)+'). SHA256 `'+hashlib.sha256((R/path).read_bytes()).hexdigest()+'`.\n\n'+table(read(path)))
parts.append('## S12 Exclusion, amendment and hash evidence\n\n')
for path in ['configs/phase1_common_universe.yaml','reports/PHASE1D_R_D5_COLOR_REPAIR_AUDIT.md','reports/PHASE1D_R2_COMMON_UNIVERSE_AMENDMENT.md','reports/PHASE1D_R2_METHODS_AMENDMENT.md','reports/PHASE0C_VR_REGISTRATION_EVIDENCE_RECONSTRUCTION.md','reports/PHASE0C_V_IMAGE_PROVENANCE.md','reports/PHASE1J_PRECISION_AND_ADDITIONAL_DONOR_FEASIBILITY.md','results/final/phase2_preintegration_contract.json','results/final/master_result_registry.tsv','results/final/target_evidence_classification.tsv']:
 parts.append('- ['+path+']('+str(R/path)+'), SHA256 `'+hashlib.sha256((R/path).read_bytes()).hexdigest()+'`.\n')
write('manuscript/MorphoID_Supplementary_Tables.md','\n'.join(parts))
legends={
'F1_framework':'MorphoID framework. X is frozen Phikon-v2 morphology, C measured 15-lineage composition, N cardinal-neighbor composition, Z coordinates plus disease in the adjusted stratum, and Y a within-lineage program. Arrows depict analysis flow, not biological causation. M4−M2 and training-side cross-fitted residual prediction answer distinct conditional questions. Independent units are donors.',
'F2_cohort_shortcuts':'Cohort and shortcut landscape. (a) Solid bars show donors (6 control, 13 PF); hatched bars sections (7,19). (b) Registration qualification and technical amendment counts. (c) Four different Phase 0D metadata/coupling diagnostics: balanced accuracy with frozen donor-bootstrap 95% intervals. These diagnostic models predate the common-universe amendment. They do not measure attribution of molecular-target prediction to color or position.',
'F3_model_ladder':'Disease-adjusted donor-held-out ladder. Panels a–c show epithelial, fibroblast and macrophage programs. M0–M4 point estimates are frozen pooled out-of-fold R² from 19 donors, with target-specific n=3,219/4,869/2,661 eligible regions. No inferential error bar is attached to individual ladder estimates. Negative scores are retained. Colors consistently denote targets.',
'F4_increment_residual':'Conditional prediction estimates. (a) ΔR²=M4−M2; (b) cross-fitted residual R². Dots are adjusted frozen estimates, bars primary 95% percentile donor-block intervals (2,000 fixed-OOF resamples; 19 donors). The zero line is the no-gain reference. None passes its respective three-target BH q<0.05 family. A positive increment does not establish stable residual identifiability.',
'F5_shortcut_controls':'Disease-adjusted outcome controls. Panels a–c correspond to the three targets. Raw, Macenko and gray denote augmented M4 with different frozen embeddings; Color and Coord are standalone controls. Points are fixed donor-held-out R², with no newly computed uncertainty. Training-donor stain references avoid outer-test fitting. These controls narrow selected shortcut explanations without assigning causal contributions.',
'F6_donor_inference_stability':'Split, donor and inference sensitivity. Rows are epithelial/fibroblast/macrophage. Left column shows all 50 pooled and PF-only partition ΔR² estimates; these are dependent partitions of 19 and 13 donors, not independent experiments. Middle shows all 19 adjusted actual donor refits. Right shows method-specific adjusted ΔR² estimates/95% intervals: reference, jackknife, balanced and BCa. Method centers and estimands can differ; no common center is imposed. Zero lines and negative estimates are retained. Complete residual and loss-test inference is in S11 of Supplementary Tables.',
'S1_registration_example':'Registered H&E crop example. All twelve frozen THD0008 patches from the registration audit are retained, including low-tissue areas. Each original square patch represents 150 × 150 μm; this stated field width is the physical calibration, rather than a scale bar inferred from composite-image margins. Panel labels and coordinates are inherited. This visual QC example establishes regional usability, not measured single-cell registration accuracy. Original source is hash-linked.',
'S2_donor_dependency':'Donor dependency structure. Every admitted donor is shown with its sections, TMA and acquisition run. Sections and their regions remain grouped within donors. TMA/run groupings can cross donors. The 28-section historical audit graph is preserved at its original source; this readable table diagram displays only the frozen 26 admitted sections, without changing selection.',
'S3_common_universe_QC':'Technical common-universe amendment. Counts are 12,144 retained, 77 coordinate failures and 39 zero-valid-pixel exclusions, with no overlap. Exclusions were fixed before predictive results and applied to all comparisons; no imputation of frozen D5 values was performed. Bar sizes are counts, not biological effect sizes.',
'S4_all_partition_increments':'All frozen partition increments. Three target panels show pooled and PF-only distributions over 50 donor partitions. Original diamonds, interval bars and summary marks retain their source definitions (median, interquartile and 5th–95th percentile spread). Split spread is not a confidence interval from independent experiments. Source Phase 1H and full S8/S9 tables control interpretation.',
'S5_all_partition_residuals':'All frozen partition residual scores. Same strata and dependent split structure as S4. Negative scores are retained. Residual sign consistency varies by target and stratum; pooled/PF evidence is not an adjusted 50-partition result.',
'S6_all_donor_LOO_increment':'Actual refit donor omissions for ΔR². Rows are targets; columns pooled, disease-adjusted and PF-only. Each labeled point is one omitted donor, with 19/19/13 actual refits. Control-donor PF-only no-op rows are omitted from this graphic but preserved in the source table. There is no selection by effect direction.',
'S7_all_donor_LOO_residual':'Actual refit donor omissions for residual R². Layout, donor labels and denominators match S6. The zero line exposes residual sign changes; this refit sensitivity is distinct from resampling fixed predictions.',
'S8_m2_support':'Frozen donor-level M2 support diagnostics. Target panels relate training-support distances to held-out M2 errors across donors, retaining original annotated donors. These are descriptive associations; support distance does not establish the causal source of a morphology gain. The particularly poor adjusted fibroblast M2 baseline must be considered when interpreting its large ΔR².',
'S9_jackknife_pseudovalues':'Delete-one-donor refit pseudo-values. Six inherited panels show target-specific morphology increments and residual scores for the adjusted 19-donor analysis. Dashed lines denote original full-cohort effects. Dispersion and influential donors limit precision; pseudo-values are not additional independent participants.',
'S10_residual_inference':'Method-specific residual inference. All three adjusted targets and four uncertainty methods are retained. Intervals and points come from Phase 1I. Fibroblast residual BCa is explicitly unreliable and is not silently dropped, replaced or presented as confirmatory.',
'S11_inference_agreement':'Frozen direction/interval agreement. The two matrices show morphology increments and residual effects under each uncertainty method. Symbols encode the source interval support classes: positive, negative, ambiguous or unavailable. This is interval/direction classification, not a matrix of BH significance. Complete p, q and invalid-method details remain in S11 of Supplementary Tables.',
'S12_precision_projection':'Frozen Phase 1J precision projections. Six panels display target/effect-specific hypothetical donor totals, including fixed-OOF donor-block and pseudo-value dispersion proxies. Shading and reference curves retain their original definitions, including 10th–90th percentile bands and inverse-square-root reference. These conditional projections do not represent recruited donors, independent validation, a guaranteed power calculation, or a new Phase 2 analysis.'}
write('manuscript/MorphoID_Figure_Legends.md','# MorphoID Figure Legends\n\n'+ '\n\n'.join('## '+name+'\n\n'+txt+'\n\n![Figure '+name+']('+str(M/'figures'/f'{name}.png')+')' for name,txt in legends.items()))
write('manuscript/MorphoID_Cover_Letter_Draft.md','''# Cover letter draft — author review required

Dear Editors,

We submit for your consideration “Decomposing morphology-associated information in virtual spatial transcriptomics” [authors, journal and article type to be confirmed].

Histology-to-expression prediction is often summarized by performance scores. MorphoID asks a different methodological question: what predictive information does morphology add after measured cell composition and neighborhood structure are available? We audited paired H&E–Xenium data from 19 independent donors and 26 sections at a qualified 150 μm scale, using a frozen donor-held-out model ladder, training-side cross-fitted residual prediction, and stain/color/position controls.

The findings require distinct conclusions for three molecular programs. Fibroblast activation showed reproducible positive morphology-associated increments across repeated pooled/PF-only partitions, with finite-donor uncertainty and a small, unstable adjusted residual signal. Epithelial injury showed a positive tendency without consistent residual support. Macrophage inflammatory state showed no reliable incremental gain in the current cohort, panel, encoder and spatial unit. Strong technical associations and donor-deletion/inference sensitivity constrain interpretation. None of the primary adjusted effects met its respective three-target q<0.05 criterion.

The contribution is an auditable conditional evaluation framework that separates prediction from molecular identifiability while preserving negative results. Fifty partitions reuse the same donors and are not independent replication. The manuscript states the 19-donor limitation, targeted panel, registration qualification and absence of independent validation. Historical HOLD verdicts remain part of the evidence. No further computational rescue or cohort expansion was performed.

Before submission, the authors must confirm authorship, corresponding contact, declarations, data/code deposition and exclusive submission. This draft does not assert that those checks have been completed.

Sincerely,
[Corresponding author and affiliation]
''')
write('manuscript/PROJECT_MANUSCRIPT_RULES.md','''# Frozen manuscript rules

- Primary evidence: 19 donors, 26 sections, 150 μm regions, three programs; frozen results only.
- No external validation, augmentation, new fitting, folds or resampling.
- No Level A; fibroblast B applies to incremental direction, epithelial C, macrophage D within the tested design.
- Preserve negative R², unstable intervals, invalid BCa and historical HOLD verdicts.
- No citations in Abstract, Methods or Results. Max three references per sentence.
- No causal, universal impossibility, clinical-validation or selective-identifiability claim.
- Main text separates metadata shortcuts from molecular outcome controls and pooled/PF partitions from adjusted inference.
- Authors must complete declarations and deposition. No automatic submission or PDF generation.
''')
print('Integrated main tables, references, methods, full supplementary tables, legends and cover letter.')
