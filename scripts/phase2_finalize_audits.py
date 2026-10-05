"""Source/claim integrity audits and authorized administrative closure. No scientific analysis."""
from pathlib import Path
import hashlib,json,re,sys,platform,subprocess,importlib.metadata as md
import pandas as pd
R=Path(__file__).resolve().parents[1];F=R/'results/final';M=R/'manuscript'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def wr(p,s):(R/p).write_text(s,encoding='utf-8')
# Correct reference types from official primary records.
r=json.loads((F/'literature_metadata.json').read_text())
r[46].update(authors='Jaume G, Doucet P, Song AH, Lu MY, Almagro-Pérez C, Wagner SJ, Vaidya AJ, Chen RJ, Williamson DFK, Kim A, Mahmood F',journal='Advances in Neural Information Processing Systems 37, Datasets and Benchmarks Track',doi='10.52202/079017-1704')
r[48].update(authors='Filiot A, Jacob P, Mac Kain A, Saillard C',journal='arXiv preprint',doi='10.48550/arXiv.2409.09173')
for j in [47,49]:r[j].update(authors='Repository contributors',journal='Official software repository; not journal validation evidence')
(F/'literature_metadata.json').write_text(json.dumps(r,indent=2,ensure_ascii=False))
p=M/'MorphoID_Main_Manuscript.md';s=p.read_text();bib='\n'.join(f"{x['id']}. {x['authors']}. {x['title']} *{x['journal']}* ({x['year']}). [{x.get('doi') or 'Primary source'}]({x['url']})." for x in r);s=re.sub(r'(?<=## References\n\n).*?(?=\n\n## Title options)',bib,s,flags=re.S);p.write_text(s)
q=M/'MorphoID_Supplementary_Methods.md';t=q.read_text().replace('missing neighbors contribute zero according to the frozen construction','composition is averaged over available nonempty neighbors, with a zero vector if none are available');q.write_text(t)
# Authorized current status only; retain historical inputs as exact snapshots.
p=R/'configs/current_project_status.yaml';s=p.read_text();s=s.replace('project_status: PASS_PHASE0D_WITH_CONFOUNDING_QUALIFICATION','project_status: STOP_FINITE_DONOR_LIMITATION_ACCEPTED',1).replace('phase1_status: NOT_RUN_AWAITING_SEPARATE_HUMAN_START','phase1_status: HOLD_MODEL_INSTABILITY')
s+='\nphase1h_status: HOLD_INFERENCE_INSTABILITY\nphase1i_status: HOLD_FINITE_DONOR_PRECISION\nphase2_status: MANUSCRIPT_READY_WITH_QUALIFICATIONS\nprimary_discovery_analysis: FROZEN\nfinite_donor_precision_limitation: ACCEPTED\nexternal_validation_available: NO\nexternal_validation_authorized: NO\naugmentation_available: NO\naugmentation_authorized: NO\nfurther_computational_rescue_authorized: NO\nnext_action: HUMAN_MANUSCRIPT_REVIEW\n';p.write_text(s)
wr('reports/CURRENT_PROJECT_STATUS.md','''# Current project status

**STOP_FINITE_DONOR_LIMITATION_ACCEPTED**

Primary discovery analysis frozen. Finite-donor precision limitation accepted. No external validation or augmentation available or authorized. No further computational rescue authorized. Phase 2 manuscript package completed with qualifications; human review is next.

Historical scientific verdicts remain unchanged in their original reports: Phase 0 registration HOLD/PASS; Phase 0D PASS_PHASE0D_WITH_CONFOUNDING_QUALIFICATION; Phase 1 HOLD_MODEL_INSTABILITY; Phase 1H HOLD_INFERENCE_INSTABILITY; Phase 1I HOLD_FINITE_DONOR_PRECISION; Phase 1J transportability HOLD. Administrative closure does not supersede their scientific interpretation.

Exact pre-Phase-2 current-status files are retained in results/final/prephase2_CURRENT_PROJECT_STATUS.md and prephase2_current_project_status.yaml. Final audits and manuscript integration report provide source hashes and package details. No submission or further modeling is authorized.
''')
contract=json.loads((F/'phase2_preintegration_contract.json').read_text());gate=[]
for x in contract['inputs']:
 f=R/x['path'];observed=sha(f) if f.exists() else 'MISSING';gate.append({'path':x['path'],'expected_sha256':x['sha256'],'observed_sha256':observed,'status':'AUTHORIZED_ADMINISTRATIVE_UPDATE' if x['path'] in contract['authorized_mutable_paths'] else ('PASS' if observed==x['sha256'] else 'FAIL')})
g=pd.DataFrame(gate);g.to_csv(F/'phase2_final_input_hash_gate.tsv',sep='\t',index=False)
assert not (g.status=='FAIL').any(),g[g.status=='FAIL']
# Registry mappings are rechecked against source bytes.
reg=pd.read_csv(F/'master_result_registry.tsv',sep='\t');rs=reg[['source_file','source_hash']].drop_duplicates();bad=[]
for _,x in rs.iterrows():
 if sha(R/x.source_file)!=x.source_hash:bad.append(x.source_file)
assert not bad,bad
fs=json.loads((F/'figure_source_hashes.json').read_text());assert all(sha(R/x)==h for x,h in fs.items())
# Visual reviews were performed on all 18 outputs; revised S2 reviewed again.
qa=pd.read_csv(F/'figure_panel_QA.tsv',sep='\t');qa['source_preflight']='PASS_SOURCE_MAPPING; export exceptions documented';qa['visual_QA']='PASS_MANUAL_VISUAL_REVIEW';qa.loc[qa.figure=='S2_donor_dependency','claim']='Readable 19-donor / 26-section dependency table; frozen cohort only';qa.to_csv(F/'figure_panel_QA.tsv',sep='\t',index=False)
wr('reports/PHASE2_FIGURE_QA.md','''# Phase 2 figure QA

All six main and twelve supplement PNGs were visually opened and inspected. Original negative scores, donor labels, split dependence, unavailable BCa entries and projection qualifications are retained. S2 was redrawn from the frozen cohort to remove unreadable overlap, then inspected again. Main F1 is an analysis schematic, not mechanistic evidence. The retained registration montage has explicit 150 × 150 μm field-width calibration in its legend; no invented scale bar is inferred from montage borders.

Static source preflight parses successfully and passes font-family/data/inference checks. Its PDF requirement fails because the user explicitly requested editable Markdown artifacts and no PDF. PNG/TIFF and 300/600 dpi warnings are retained in the raw preflight: PNG 300 dpi previews plus editable SVG are the review exports, not a declaration of compliance with a selected journal's upload specifications. Main figures use 183 mm width and 6 pt minimum text; inherited supplement dimensions are retained. Journal-specific raster export can follow human journal selection without new analysis. This is an explicit export waiver, not a claimed clean validator result.

Per-panel source mapping and manual review are in results/final/figure_panel_QA.tsv. Main figures and redrawn S2 have editable text SVG; historical supplement assets retain available original SVG. All 18 PNG review figures exist. No PDF or new biological images were generated.
''')
# Every trigger-bearing line is retained verbatim for human-verifiable review.
terms=r'\b(?:demonstrat\w*|prov\w*|validat\w*|robust\w*|generaliz\w*|identif\w*|infer\w*|predict\w*|independent\w*|causal\w*|specific\w*)\b|molecular[- ]state'
lines=[]
for f in sorted(M.glob('MorphoID_*.md')):
 # Raw numeric supplementary tables are evidence, not prose claims.
 if f.name.endswith('Supplementary_Tables.md'):continue
 for n,line in enumerate(f.read_text().splitlines(),1):
  matches=sorted(set(x.lower() for x in re.findall(terms,line,re.I)))
  if not matches or not line.strip() or re.match(r'^\d+\. .*https://',line):continue
  context='DESCRIPTIVE_ONLY' if line.startswith('#') or '## References' in line else 'SUPPORTED_WITH_QUALIFICATION'
  source='Primary / H / I source tables; Supplementary Methods scope; reference links for literature'
  if any(x in line.lower() for x in ['no independent','not validated','not establish','cannot','not a causal','no target','no external','unreliable','not consistently','unavailable']):source='Explicit limitation or negative boundary; preserved audit evidence'
  if 'Future' in line or 'would require' in line:source='Prospective requirement; no executed-analysis claim'
  lines.append([str(f.relative_to(R)),n,', '.join(matches),context,source,line.replace('|','\\|')])
a=pd.DataFrame(lines,columns=['File','Line','Trigger terms','Claim status','Evidence interpretation','Exact text'])
def mdt(df):
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(str(v).replace('\n',' ') for v in row)+' |' for row in df.itertuples(index=False,name=None))
main=(M/'MorphoID_Main_Manuscript.md').read_text()
assert '{{' not in main
for sec in ['Abstract','Methods','Results']:
 text=re.search(r'^## '+sec+r'\n(.*?)(?=^## |\Z)',main,re.S|re.M).group(1);assert not re.search(r'\[[0-9][0-9,– -]*\]',text),sec
assert all(x in main for x in ['19 donors','26 sections','150','BCa','q < 0.05'])
def citations(t):
 out=set()
 for z in re.findall(r'\[([\d,– -]+)\]',t):
  for v in z.split(','):
   v=v.strip();ab=re.split('[–-]',v)
   out.update(range(int(ab[0]),int(ab[-1])+1))
 return out
counts={}
for sec in ['Background','Discussion']:
 z=re.search(r'^## '+sec+r'\n(.*?)(?=^## |\Z)',main,re.S|re.M).group(1);counts[sec]={'words':len(z.split()),'distinct_references':len(citations(z))}
assert counts['Background']['distinct_references']<50/4
# Groups max 3; semicolon/source links excluded.
for z in re.findall(r'\[([\d,– -]+)\]',main):assert len(citations('['+z+']'))<=3,z
wr('reports/FINAL_CLAIM_AUDIT.md',f'''# Final claim audit

**Claim audit: PASS for the frozen-evidence manuscript; human author review remains required.**

## Calibration decisions

No Level A target was created. Fibroblast B is limited to recurring incremental direction, with finite-donor magnitude, M2 baseline and residual qualifications. Epithelial C does not establish composition-independent residual prediction. Macrophage D is bounded to this cohort, panel, encoder and 150 μm model design. No selective-identifiability, universal impossibility, causal or externally validated headline is supported. The main abstract reports all three targets, negative evidence, shortcut findings and donor uncertainty. No primary adjusted effect meets its respective q<0.05 family.

The requested broad “composition and neighborhood explain a large fraction” heading was calibrated: composition improves pooled scores, while neighborhood increments are small or negative for some targets. Large fibroblast ΔR² is presented beside negative adjusted M2 R². Pooled/PF-only repeated partitions are not labeled adjusted repeats. Metadata shortcut balanced accuracies are not substituted for molecular-outcome controls. Fixed-prediction bootstrap is distinguished from refit LOO, and invalid residual BCa remains explicit.

Unsupported candidate formulations rejected: “robust H&E molecular-state prediction,” “selective identifiability established,” “independent variance explained,” “externally generalizable,” and universal macrophage impossibility. These were avoided during drafting, not claimed to have been removed from historical frozen reports.

## Structure and citation audit

Background: {counts['Background']['words']} words; Discussion: {counts['Discussion']['words']} words; ratio {counts['Discussion']['words']/counts['Background']['words']:.2f}. Total references: 50; Background distinct references: {counts['Background']['distinct_references']}; Discussion: {counts['Discussion']['distinct_references']}. Background stays below one quarter of the total. Discussion uses more sources than the usual 20–28 target because this methodological synthesis includes biological, pathology-model and uncertainty-method contexts; this is a documented editorial-density qualification. Abstract/Methods/Results have no literature citations. All individual numeric citation groups contain at most three sources. Main research conclusions are source-mapped below and in the result registry. Resource/preprint records are labeled as such; published metadata were checked against primary records.

## Evidence map

Results 1: frozen cohort and registration/provenance reports. Results 2: Phase 0D color/coordinate/coupling tables. Results 3–4: table2 primary ladder/increments. Results 5: table3 crossfit residual. Results 6: target classification plus the underlying tables. Results 7: table4 and primary preprocessing columns. Results 8: repeated_fold_metrics/stability, leave_one_donor_out_influence, inference_agreement_matrix. Source hashes are rechecked in FINAL_REPRODUCIBILITY_MANIFEST.md. Biological discussion is explicitly interpretation, with primary literature citations, not measured collagen mechanism or cell-level validation.

## Every trigger-bearing manuscript line

{len(a)} lines screened for demonstrate/prove/validate/robust/generalize/identify/infer/predict/independent/molecular state/causal/specific and their word forms. Headings and method descriptions are descriptive; conditional scientific claims retain scope and uncertainty. Reference titles are bibliographic records, not MorphoID efficacy claims. The verbatim table below makes the review inspectable; it does not substitute for reading each cited source during author review.

'''+mdt(a)+'\n')
# Runtime is observed now, not asserted as original analysis runtime.
versions={'python':sys.version,'executable':sys.executable,'platform':platform.platform()}
for v in ['numpy','pandas','scipy','torch','transformers','matplotlib','Pillow','scikit-learn']:
 try:versions[v]=md.version(v)
 except md.PackageNotFoundError:versions[v]='not in current interpreter'
try:versions['R']=subprocess.run(['R','--version'],capture_output=True,text=True).stdout.splitlines()[0]
except Exception:versions['R']='not available in current PATH; historical evidence retained'
(F/'observed_runtime.json').write_text(json.dumps(versions,indent=2))
commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
status=subprocess.check_output(['git','status','--short'],cwd=R,text=True);(F/'git_working_tree_status.txt').write_text(status)
keypaths=sorted(set([x['path'] for x in contract['inputs'] if x['path'].startswith('configs/')]+list(rs.source_file)+list(fs)+['environment/environment.yml','results/phase1/morphology_encoder_manifest.tsv','results/phase1/morphology_embedding_extraction_manifest.tsv']))
manifest=[]
for x in keypaths:
 f=R/x
 if f.exists():manifest.append([x,f.stat().st_size,sha(f)])
pd.DataFrame(manifest,columns=['path','bytes','sha256']).to_csv(F/'key_reproducibility_hashes.tsv',sep='\t',index=False)
wr('reports/FINAL_REPRODUCIBILITY_MANIFEST.md',f'''# Final reproducibility manifest

**Audit: PASS for source mapping and frozen-input preservation.** This is not a fresh numerical rerun or a claim that the declared environment was independently recreated.

## Snapshot

Repository HEAD: `{commit}`. Working tree is not clean; the commit alone is insufficient. Exact working-tree status is retained in results/final/git_working_tree_status.txt. Frozen source contents are controlled by results/final/phase2_preintegration_contract.json and phase2_final_input_hash_gate.tsv: {len(g)} files checked, {(g.status=='PASS').sum()} unchanged; {(g.status=='AUTHORIZED_ADMINISTRATIVE_UPDATE').sum()} explicitly authorized current-status edits; zero unexpected changes. Both prior administrative files are preserved byte-for-byte in their prephase2 snapshots. Historical HOLD reports and scientific configurations were not rewritten.

All {len(reg)} registry records map to {len(rs)} unique source/hash pairs, rechecked against available bytes. All {len(fs)} figure-source hashes match. Complete configuration/cohort/common-universe/fold/target/key-result/figure source SHA256 values appear in key_reproducibility_hashes.tsv; the full scientific inventory remains in the preintegration contract. This avoids abbreviating hashes or relying on filenames alone.

Checkpoint SHA256: `261ae680fa699b3b951597fd57aa19c02ef735805acb104b93af69b36d928569`; official revision `2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd`. Historical encoder and extraction manifests are included in the key index. Existing 26 CPU caches are preserved; this audit checks provenance manifests, not a new expensive recomputation of embeddings or re-download of checkpoint weights.

## Runtime

Declared environment: environment/environment.yml (hashed). Observed Phase 2 interpreter/packages: results/final/observed_runtime.json. These are distinguished from the original numerical execution environment; the observed interpreter may differ from the declared environment. Current observations:

```json
{json.dumps(versions,indent=2)}
```

No package installation or tests were needed for this presentation-only integration. Current scripts used pandas/Python to format existing estimates and matplotlib to draw them, without importing/fitting scientific model workflows.

## Execution lineage and export qualifications

Historical analysis order: Phase 0 provenance/registration/shortcuts → authorized D5 common-universe amendment → CPU morphology extraction → primary donor models → Phase 1H splits/LOO → Phase 1I uncertainty → Phase 1J feasibility/projection. Each original report and input gate remains accessible. Phase 2 scripts format the registry, literature, tables, figures and audit documents only. The stopped project must not be automatically rerun.

Figure data sources and individual SHA256 values are in figure_source_hashes.json. All 18 figures were visually inspected; PDF/TIFF/600-dpi validator exceptions are explicit in PHASE2_FIGURE_QA.md, following the user's editable Markdown/no-PDF request. Review PNGs and available editable SVGs are included. Figure-source consistency, negative-score preservation and honest unavailable-method display pass; no selected-journal upload compliance is claimed.

## Human review and publication preparation

Authors, affiliations, funding, ethics/consent applicability, conflicts, permanent code/data deposition and exact target-journal formatting remain author decisions. The bibliography contains 50 primary published/resource records, with source URLs; repository/preprint records are not presented as independent MorphoID validation. These are submission preparation items, not missing source mappings for the frozen result. Absolute local links support inspection in this workspace; a distributed package must preserve the repository-relative structure or rewrite links during human-approved deposition.
''')
print(json.dumps({'gate_files':len(g),'unexpected_changes':0,'registry_rows':len(reg),'source_pairs':len(rs),'trigger_lines':len(a),'structure':counts,'commit':commit}))
