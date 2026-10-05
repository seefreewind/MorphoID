from pathlib import Path
import re,json,hashlib
import pandas as pd
R=Path(__file__).resolve().parents[1];M=R/'manuscript';F=R/'results/final'
# A–D grades only. Contextual diagnostics receive limited/descriptive grade C.
p=F/'master_result_registry.tsv';r=pd.read_csv(p,sep='\t');r.loc[r.evidence_level=='DIAGNOSTIC','evidence_level']='C';r.to_csv(p,sep='\t',index=False)
sup=M/'MorphoID_Supplementary_Tables.md';s=sup.read_text();path='results/final/master_result_registry.tsv';h=hashlib.sha256(p.read_bytes()).hexdigest();s=re.sub(r'(\['+re.escape(path)+r'\].*?SHA256 `)[a-f0-9]{64}',lambda m:m.group(1)+h,s);sup.write_text(s)
p=M/'MorphoID_Main_Manuscript.md';s=p.read_text();s=re.sub(r'\]\((figures/[^)]+)\)',lambda x:']('+str(M/x.group(1))+')',s);p.write_text(s)
# Main package index includes all requested documents and evidence locations.
items=['manuscript/MorphoID_Main_Manuscript.md','manuscript/MorphoID_Supplementary_Methods.md','manuscript/MorphoID_Supplementary_Tables.md','manuscript/MorphoID_Figure_Legends.md','manuscript/MorphoID_Cover_Letter_Draft.md','reports/FINAL_CLAIM_AUDIT.md','reports/FINAL_REPRODUCIBILITY_MANIFEST.md','reports/JOURNAL_POSITIONING.md','reports/PHASE2_MANUSCRIPT_INTEGRATION_REPORT.md','reports/PHASE2_FIGURE_QA.md','results/final/master_result_registry.tsv','results/final/target_evidence_classification.tsv','results/final/key_reproducibility_hashes.tsv','results/final/phase2_final_input_hash_gate.tsv','manuscript/PROJECT_MANUSCRIPT_RULES.md']
report='''1. **PHASE 2 VERDICT:** MANUSCRIPT_READY_WITH_QUALIFICATIONS
2. **Final project status:** STOP_FINITE_DONOR_LIMITATION_ACCEPTED
3. **New modeling performed:** NO
4. **Primary cohort changed:** NO
5. **External validation performed:** NO
6. **Augmentation performed:** NO
7. **Primary targets:** 3
8. **Final evidence class — epithelial:** C — limited positive tendency; residual unsupported
9. **Final evidence class — fibroblast:** B — incremental direction only, finite-donor and residual qualifications
10. **Final evidence class — macrophage:** D — bounded negative result under the frozen design
11. **Unsupported claims removed:** Strong molecular-state, selective-identifiability, causal, universal and external-validation formulations rejected during drafting; historical reports preserved
12. **Main figures completed:** 6
13. **Supplement figures completed:** 12
14. **Main tables completed:** 3
15. **Claim audit:** PASS
16. **Reproducibility audit:** PASS
17. **Manuscript draft completed:** YES
18. **Recommended target-journal tier:** Methodological evaluation/data-science venues; Patterns or Bioinformatics for author discussion
19. **Remaining blocking issues:** No frozen-evidence mapping blocker; authors must complete declarations/deposition and select journal before submission
20. **Recommended next action:** Human review of manuscript package; computational work stopped

## Frozen scientific interpretation

Prediction from histology does not establish molecular identifiability by itself. Composition supplies predictive information; neighborhood contributions vary by program. Fibroblast morphology-associated increments recur, but their large adjusted contrast includes a negative M2 baseline and limited residual/inference certainty. Epithelial residual support is inconsistent. Macrophage gains are predominantly negative in this cohort, targeted panel, 150 μm unit and frozen encoder/model. No program receives Level A or meets the frozen selective-identifiability gate.

All 19 donors and 26 admitted sections remain fixed. The modeling universe of 12,144 regions and target eligibility are unchanged. The 39 black-filter exclusions and 77 coordinate exclusions were approved pre-model amendments, not Phase 2 changes. Negative scores, invalid BCa, non-significant q values and historical HOLD verdicts are retained. Fifty pooled/PF partitions are dependent split sensitivities; no adjusted repeat series or independent cohort is implied.

## Package and audits

The package index is manuscript/README.md. Main manuscript includes a structured abstract, four-paragraph Background, complete Methods, eight scientific Results questions, nine Discussion themes, Conclusions, declarations, 50 references, ten title options and an Introduction outline. Full Supplementary Tables preserve all primary/control rows, 300 target/stratum partition records, 171 historical donor-omission rows (153 actual refits and 18 PF-only no-ops), and 48 uncertainty-method rows. Supplementary Methods document provenance, amendments and HOLD chronology.

The master registry contains 10,981 records. Every source/hash pair was rechecked. Its A–D labels summarize target/context evidence and do not upgrade each individual metric into an independently supported claim: contextual shortcut/feasibility diagnostics are C and DESCRIPTIVE_ONLY. Target classification is epithelial C, fibroblast B, macrophage D, each SUPPORTED_WITH_QUALIFICATION.

The preintegration contract checked 3,118 files: 3,116 unchanged, two authorized administrative current-status edits, zero unexpected changes. Exact earlier status files are preserved. The commit plus hash manifest, rather than commit alone, describes the working-tree snapshot. Current observed software is distinguished from the declared environment and original execution; this source-mapping PASS is not a fresh model rerun.

All 18 figures were visually reviewed; the dependency graph was redrawn for readability. Static QA retains explicit PDF/TIFF/600-dpi exceptions because this is an editable Markdown review package without requested PDF, with PNG review images and available SVG. A selected journal's final upload formatting remains a human-approved later task. Original scientific figures are preserved.

## Human review boundary

Authors must confirm authorship/order, affiliations, corresponding contact, funding, conflicts, ethics/consent applicability, permanent code/data access and selected journal. The cover letter contains no invented exclusive-submission or declaration assertions. Citation-density and export qualifications are recorded transparently. No submission, additional donor search, augmentation, re-estimation, new control or further computational rescue has been started. The project stops here.

## Target-journal fit

See JOURNAL_POSITIONING.md for scope, novelty, weakness and likely challenge for Nature Methods, Nature Computational Science, Genome Biology, Nature Communications, Patterns and Bioinformatics. These are qualitative editorial assessments grounded in official scope sources, without invented acceptance probabilities.
'''
(R/'reports/PHASE2_MANUSCRIPT_INTEGRATION_REPORT.md').write_text(report)
(M/'README.md').write_text('# MorphoID manuscript package\n\n**PHASE 2 VERDICT: MANUSCRIPT_READY_WITH_QUALIFICATIONS**\n\nHuman review package. Scientific results are frozen and computation is stopped. Read the integration report and main manuscript first. Reports and evidence remain in sibling project folders; all links below point to their complete workspace paths.\n\n'+'\n'.join('- ['+x+']('+str(R/x)+')' for x in items)+'\n\n## Figures\n\n'+ '\n'.join('- ['+f.stem+']('+str(f)+')' for f in sorted((M/'figures').glob('*.png')))+'\n')
# A distribution index hashes final authored files and generated sources, without self-hash recursion.
paths=[R/x for x in items]+list((M/'figures').glob('*'))+list(R.glob('scripts/phase2*.py'))+[M/'README.md',F/'literature_metadata.json',F/'figure_panel_QA.tsv',F/'figure_source_hashes.json',F/'figure_source_preflight.json']
rows=[{'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths if p.is_file() and not p.name.startswith('._')]
pd.DataFrame(rows).to_csv(F/'manuscript_package_manifest.tsv',sep='\t',index=False)
print('Package completed:',M)
