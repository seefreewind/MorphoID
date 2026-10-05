#!/usr/bin/env python3
"""Metadata-only Phase J audit. No image/count download or model execution."""
from pathlib import Path
import re, json, gzip, hashlib, datetime
import numpy as np
import pandas as pd

R=Path(__file__).resolve().parents[1]
O=R/'results/phase1j'; S=O/'source_metadata'
GENES={'epithelial_injury':['ITGB6','KRT18','KRT8','MMP7','SOX4'],
       'fibroblast_activation':['ACTA2','COL1A1','COL1A2','COL3A1','CTHRC1','FN1','POSTN'],
       'macrophage_inflammatory':['CCL2','IL1B','S100A8','S100A9','SPP1','TNF']}
def save(d,name): d.to_csv(O/name,sep='\t',index=False)
def main():
 old=pd.read_csv(R/'results/phase0/phase0c_v/vannan_sample_crosswalk.tsv',sep='\t').fillna('UNKNOWN')
 frozen=pd.read_csv(R/'configs/frozen_vannan_primary_cohort.tsv',sep='\t')
 primary=set(frozen.donor_id); admitted=set(frozen.loc[frozen.primary_inclusion=='YES','sample_id'])
 assert len(primary)==19 and len(admitted)==26
 text=(S/'GSE250346_current.soft.txt').read_text()
 samples={b.splitlines()[0].strip():b for b in text.split('^SAMPLE = ')[1:]}
 assert set(samples)==set(old.GSM), 'Current GEO samples differ from historical full crosswalk'
 panel=set(pd.read_csv(R/'results/phase0/phase0c_v/gene_panel_registry.tsv',sep='\t').gene_symbol)
 assert set(sum(GENES.values(),[]))<=panel
 rows=[]
 for _,r in old.iterrows():
  b=samples[r.GSM]; files=re.findall(r'!Sample_supplementary_file(?:_\d+)? = (.+)',b)
  image=[u.strip() for u in files if 'registered' in u.lower() and any(w in u.lower() for w in ['.tif','.png','.jpeg'])]
  new=r.donor_id not in primary
  ann=R/'results/phase0/phase0c_vr/author_annotations/GSE250346_HE_annotations/IPFTMA5_annotations'/r.sample_id
  annfiles=[p.name for p in ann.iterdir() if p.is_file() and not p.name.startswith('._')] if ann.exists() else []
  rows.append(dict(sample_id=r.sample_id,donor_id=r.donor_id,GSM=r.GSM,publication_primary=r.published_primary,
   post_publication_extension='NO_FINAL_2025_ARTICLE_MEMBER',same_donor_as_primary='YES' if not new else 'NO',
   new_independent_donor='YES' if new else 'NO',TMA=r.TMA,run=r.run,registered_HE='YES_AUTHOR_FILE' if image else 'NO_PER_CORE_REGISTERED_FILE',
   expression=r.expression_available,coordinates=r.coordinates_available,
   annotation=r.cell_annotation_available,panel='PASS_343_GENES_18_OF_18_FROZEN_GENES',
   same_section_provenance='YES_AUTHOR_POST_XENIUM_HE_PROTOCOL',original_admitted='YES' if r.sample_id in admitted else 'NO',
   historical_preprint_extension=r.geo_extension,formal_2025_article_member=r.formal_2025_article_member,
   registered_HE_URL=';'.join(image),annotation_local_files=';'.join(annfiles),
   coordinate_validity='PASS_FROZEN_PRIMARY' if r.sample_id in admitted else 'PENDING_FROZEN_ADMISSION_QA',
   morphology_150um='PASS_FROZEN_PRIMARY' if r.sample_id in admitted else 'PENDING_TISSUE_CROP_QA',
   nuclear_transcript_unit='PASS_FROZEN_PRIMARY' if r.sample_id in admitted else 'AUTHOR_NUCLEAR_ANALYSIS_SUPPORTED_EXTENSION_UNIT_CHECK_PENDING',
   source_URL='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+r.GSM,
   notes=str(r.notes)+' | J audit: TMA5 is already included in the formal 2025 article. Independent of discovery donors, same study/technical context. Registered filenames alone are not numerical registration PASS.'))
 cross=pd.DataFrame(rows); save(cross,'vannan_extension_crosswalk.tsv')
 ext=cross[cross.new_independent_donor=='YES']; ready=ext[ext.registered_HE=='YES_AUTHOR_FILE']
 assert ext.donor_id.nunique()==16 and ready.donor_id.nunique()==13
 coverage=[]
 for i in range(1,5):
  j=json.loads(gzip.decompress((S/f'COPD_TMA{i}_gene_panel.json.gz').read_bytes()))
  g={x['type']['data'].get('name','') for x in j['payload']['targets'] if x['type'].get('descriptor')=='gene'}
  assert len(g)==480
  for t,genes in GENES.items():
   for gene in genes:coverage.append(dict(study='COPD_XENIUM',sample=f'TMA{i}',target=t,gene=gene,in_panel=gene in g,source=f'source_metadata/COPD_TMA{i}_gene_panel.json.gz'))
 for t,genes in GENES.items():
  for gene in genes:coverage.append(dict(study='VANNAN_EXTENSION',sample='ALL_TMA5_AUTHOR_SHARED_PANEL',target=t,gene=gene,in_panel=gene in panel,source='results/phase0/phase0c_v/gene_panel_registry.tsv'))
 cov=pd.DataFrame(coverage);save(cov,'frozen_gene_coverage.tsv')
 assert cov[cov.study=='COPD_XENIUM'].groupby('sample').in_panel.sum().eq(10).all()
 # Gate strings reflect observed metadata; UNKNOWN never means PASS.
 entries=[]
 def add(key,study,accession,platform,disease,n,same,he,expr,coords,ann,panel,composition,reg,access,role,status,reason,url,C,N,Z,scale,preprocess,independent=True):
  entries.append(dict(candidate_id=key,study=study,accession=accession,platform=platform,tissue='human lung',disease=disease,n_donors=n,
   same_section=same,**{'H&E':he},expression=expr,coordinates=coords,annotation=ann,panel_compatible=panel,
   composition_possible=composition,registration_evidence=reg,public_access=access,candidate_role=role,admission_status=status,
   reason=reason,source_URL=url,independent_study=independent,target_definitions=panel,C_ontology=C,N_structure=N,Z_structure=Z,
   morphology_150um=scale,encoder_preprocessing=preprocess,metadata_audit_date='2026-10-03'))
 add('VANNAN_TMA5','Vannan extension','GSE250346','Xenium','PF/control',16,'YES','13_NEW_DONORS_AUTHOR_REGISTERED;3_MISSING_PER_CORE',
  'YES','YES','AUTHOR_SHARED_ONTOLOGY','PASS_18_OF_18','PENDING_EXTENSION_NUCLEAR_UNIT_CHECK',
  'AUTHOR_MAPPING_PRESENT_13;NUMERICAL_QA_PENDING','PUBLIC','AUGMENTATION_CANDIDATE','PENDING',
  '16 independent discovery-external donors; same study. 13 have per-core registered H&E; 0 newly admitted under frozen numerical gates.',
  'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346','AUTHOR_MAPPING_AVAILABLE_PENDING_UNIT_CHECK','PENDING_COORDINATE_QA','PF_CONTROL_METADATA_AVAILABLE','PENDING_CROP_QA','RGB_HE_AVAILABLE_PENDING_QA',False)
 add('CHUV_XENIUM','Bilous CHUV lung panel cohort','GSE311609','Xenium','NSCLC',10,'YES_AUTHOR_PROTOCOL','AVAILABLE_8_DONORS;L9_L10_MISSING',
  'YES','YES','CELL_ANNOTATIONS_AVAILABLE','5K_EXACT_GENES_NOT_VERIFIED;OTHER_PANELS_INCOMPLETE','POSSIBLE_MAPPING_PENDING',
  'FROZEN_PHASE0_S_0_PASS;FAILED_OR_UNRESOLVED_FRAGMENT_MAPPING','PUBLIC','SUPPORTIVE_ONLY','NOT_ADMITTED',
  '22 panel/section records are 10 donors; 5K has six distinct donors. Existing registration audit admitted none; PF/control disease structure absent.',
  'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311609','PENDING','PENDING_REGISTERED_COORDINATES','FAIL_PF_CONTROL_ABSENT','FAIL_OR_UNRESOLVED_REGISTRATION','HE_AVAILABLE_GEOMETRY_UNRESOLVED')
 add('COPD_XENIUM','Sauler COPD spatial cohort','GSE313006','Xenium','COPD',38,'HE_PAIR_UNKNOWN','NO_REGISTERED_HE_IN_CURRENT_SUPPLEMENT_LIST',
  'YES','YES','CELL_ANNOTATIONS_PAPER;DONOR_LINKAGE_CONTROLLED','FAIL_10_OF_18','PENDING_CONTROLLED_DONOR_LINKAGE',
  'MORPHOLOGY_OME_IS_FLUORESCENCE;HE_TRANSFORM_NOT_FOUND','COUNTS_PUBLIC;LINKED_METADATA_CONTROLLED_BIOLINCC','SUPPORTIVE_ONLY','NOT_ADMITTED',
  'Four TMA records represent 38 participants, not four donors. All four panels lack eight frozen genes. No signature substitution allowed.',
  'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313006','PENDING','PENDING','FAIL_COPD_NOT_FROZEN_PF_STATUS','UNKNOWN_HE_SCALE','NO_PUBLIC_REGISTERED_RGB_HE')
 add('COSMX_COVID22','Columbia immuno-lung spatial cohort','LungSpatialDB;GSE287114_SEQUENCE_ONLY','CosMx','COVID explant/autopsy/control',22,
  'YES_POST_COSMX_HE_PROTOCOL','METHOD_CONFIRMED;DOWNLOAD_AND_TRANSFORM_UNKNOWN','PAPER_YES_PUBLIC_PORTAL_UNVERIFIED','PAPER_YES_PUBLIC_FILE_UNVERIFIED',
  'AUTHOR_CELLTYPES;PUBLIC_CROSSWALK_UNVERIFIED','UNKNOWN_EXACT_18_GENES','POSSIBLE_PENDING_ANNOTATION_AND_SAMPLING',
  'SAME_SLIDE_STAINING_CONFIRMED;NUMERICAL_TRANSFORM_UNKNOWN','PORTAL_NOT_ACCESSIBLE_IN_THIS_AUDIT','SUPPORTIVE_ONLY','PENDING',
  '22 donors, 400 FOVs. GEO GSE287111 is organoid scRNA-seq, not the 22-donor CosMx matrix. Sparse FOV sampling and PF/control absence block unchanged estimand.',
  'https://pmc.ncbi.nlm.nih.gov/articles/PMC12463091/','PENDING','PENDING_FOV_SAMPLING','FAIL_COVID_NOT_PF','PENDING_FOV_BOUNDARIES','POST_STAIN_HE_EXISTS_PUBLIC_MAPPING_UNKNOWN')
 add('MAYR_IPF','Mayr IPF spatial study','Zenodo10012934;10015169','Xenium/Visium','IPF/control',7,'XENIUM_ADJACENT_VISIUM_HE','VISIUM_HE_PRESENT',
  'YES','YES','AUTHOR_ANNOTATIONS','289_PANEL_EXACT_18_NOT_VERIFIED','XENIUM_POSSIBLE;VISIUM_DECONVOLUTION_NOT_ORACLE',
  'XENIUM_TO_HE_ADJACENT_SECTION','PUBLIC','SUPPORTIVE_ONLY','NOT_ADMITTED',
  'Seven donors across Visium; do not infer seven distinct Xenium donors. Adjacent-section pairing fails strict same-section Xenium requirement.',
  'https://pmc.ncbi.nlm.nih.gov/articles/PMC11313858/','PENDING','VISIUM_SPOT_NEIGHBORHOODS_DIFFER','PF_CONTROL_AVAILABLE','FAIL_XENIUM_ADJACENT_HE','VISIUM_HE_SUPPORTIVE_ONLY')
 add('PF_VISIUM8','Pulmonary fibrosis Visium cohort','S-BSST1410','Visium','IPF/control',8,'VISIUM_SAME_SECTION','YES',
  'YES','YES','CELL2LOCATION_INFERRED_COMPOSITION','WHOLE_TRANSCRIPTOME;EXACT_ASSAY_COVERAGE_UNVERIFIED','NO_OBSERVED_WITHIN_LINEAGE_ORACLE',
  'SPACE_RANGER_IMAGE_MAPPING;NO_FROZEN_ADMISSION_QA','PUBLIC','SUPPORTIVE_ONLY','NOT_ADMITTED',
  'Four IPF and four control donors. Spot mixtures cannot supply unchanged single-cell within-lineage targets and true composition oracle.',
  'https://www.nature.com/articles/s41588-024-01819-2','FAIL_INFERRED_COMPOSITION','FAIL_SPOT_UNIT','PF_CONTROL_AVAILABLE','PENDING_NUMERICAL_QA','HE_AVAILABLE')
 add('ADULT_MERFISH6','Adult lung multimodal atlas','HuBMAP/LungMAP publication a10041ad9ebae0b42d3c7f602ba37b82','MERFISH','healthy',6,
  'EXACT_HE_MERFISH_PAIR_UNVERIFIED','SERIAL_HISTOLOGY_DOCUMENTED','YES','YES','REFERENCE_TRANSFER_AUTHOR_LABELS','503_GENES_EXACT_18_UNKNOWN',
  'POSSIBLE_PENDING_MAPPING','SERIAL_HE_NOT_SAME_SECTION_TRANSFORM','PUBLIC_PORTAL','SUPPORTIVE_ONLY','PENDING',
  'Six spatial donors; eleven overall atlas donors are not six additional spatial donors. Same-section H&E not established; healthy-only cannot reproduce PF-adjusted comparison.',
  'https://pmc.ncbi.nlm.nih.gov/articles/PMC12140004/','PENDING','PENDING','FAIL_HEALTHY_ONLY','UNKNOWN_SAME_HE_SECTION','UNKNOWN_REGISTERED_RGB_HE')
 add('NSCLC_MERFISH46','Chen NSCLC immunity hubs','SCP2510;Zenodo11198494','MERFISH','NSCLC',46,'EXACT_HE_PAIR_UNKNOWN','UNKNOWN_PUBLIC_REGISTERED_HE',
  'YES','YES','AUTHOR_LABELS_PUBLIC','484_GENES_EXACT_18_UNKNOWN','POSSIBLE_PENDING_MAPPING',
  'PANEL_COREGISTRATION_IS_NOT_HE_REGISTRATION','PUBLIC','SUPPORTIVE_ONLY','PENDING',
  '46 processed patient datasets documented by author deposit. RNA/protein coregistration does not establish a registered H&E partner or PF/control transport.',
  'https://zenodo.org/records/11198494','PENDING','PENDING','FAIL_PF_CONTROL_ABSENT','UNKNOWN_HE_PAIR','UNKNOWN_REGISTERED_RGB_HE')
 add('NEUMAP_HD8','NeuMap human lung TMA','GSE266680;GSE298732;GSM9021740','Visium HD','NSCLC/adjacent lung',8,'YES_VISIUM_HD_HE_WORKFLOW','YES_AUTHOR_IMAGE_AND_REGISTRATION_METADATA',
  'YES','YES','BIN2CELL_INFERRED_CELLS_TYPES','WHOLE_TRANSCRIPTOME_PROBE_ASSAY;EXACT_18_UNVERIFIED','RECOVERABILITY_PENDING_NOT_TRUE_ORACLE_PASS',
  'REGISTRATION_JSON_SCALEFACTORS;NUMERICAL_QA_NOT_RUN','PUBLIC','SUPPORTIVE_ONLY','PENDING',
  'Eight patients, three cores per patient, not 24 donors. Core-to-patient mapping and cell-state recovery need verification; PF/control Z absent.',
  'https://www.nature.com/articles/s41586-025-09807-0','PENDING_BIN2CELL_ONTOLOGY','PENDING_RECOVERED_CELLS','FAIL_PF_CONTROL_ABSENT','AUTHOR_SCALE_SUPPORTED_QA_PENDING','HE_AVAILABLE')
 add('EGFR_HD4','EGFR lung Visium HD response study','GSE288758','Visium HD','NSCLC EGFR treatment',4,'YES_VISIUM_HD_WORKFLOW','YES_METADATA',
  'YES','YES','CELL_COMPOSITION_RECOVERY_UNVERIFIED','WHOLE_TRANSCRIPTOME_PROBE_ASSAY;EXACT_18_UNVERIFIED','PENDING',
  'AUTHOR_WORKFLOW_ONLY','PUBLIC','SUPPORTIVE_ONLY','PENDING',
  'Four named patients ER, SR, LR1, LR2; pre/post specimens are repeated donors. Below preferred five, no frozen PF/control comparison.',
  'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288758','PENDING','PENDING','FAIL_PF_CONTROL_ABSENT','PENDING','HE_METADATA_AVAILABLE')
 add('ASTHMA_XENIUM','Asthma endobronchial spatial cohort','GSE269354','Xenium','asthma/control','8_CONFIRMED_COHORT1;16_SPECIMEN_IDS_UPPER_BOUND','HE_ADDITIONAL_SECTION','HE_HISTOLOGY_ADJACENT',
  'YES','YES','AUTHOR_LABELS','339_GENES_EXACT_18_UNVERIFIED','POSSIBLE_PENDING','ADDITIONAL_SLIDES_FOR_HE_NO_SAME_SECTION_MAPPING',
  'PUBLIC','SUPPORTIVE_ONLY','NOT_ADMITTED',
  'Metadata lists eight distinct biopsy IDs in each of two historical cohorts; patient cross-cohort duplication cannot be independently ruled out. Sixteen specimen IDs are a donor upper bound, not a certified independent count.',
  'https://www.nature.com/articles/s41590-025-02161-3','PENDING','PENDING','FAIL_ASTHMA_NOT_PF','FAIL_ADJACENT_HE','FAIL_SAME_SECTION_PAIR')
 add('AML_XENIUM','AML lung infiltration spatial study','GSE319763','Xenium','AML infiltration','UNKNOWN','UNKNOWN_HE_PAIR','NO_HE_FOUND_IN_SCREENED_METADATA',
  'YES','YES','PENDING_PUBLIC_LABELS','UNKNOWN_EXACT_18','PENDING','NO_VERIFIED_HE_TRANSFORM',
  'PUBLIC_GEO','SUPPORTIVE_ONLY','PENDING',
  'Three GSM records do not certify three independent donors. Disease context differs; unknown donor crosswalk and registered H&E block admission.',
  'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319763','PENDING','PENDING','FAIL_AML_NOT_PF','UNKNOWN','UNKNOWN')
 add('COSMX_COVID6','Severe COVID macrophage spatial niches','GSE327879','CosMx','COVID/control',6,'UNKNOWN_HE_PAIR','UNKNOWN_PUBLIC_REGISTERED_HE',
  'YES','YES_METADATA','PENDING','1000_GENES_EXACT_18_UNKNOWN','PENDING','NO_VERIFIED_HE_TRANSFORM',
  'PUBLIC_GEO','SUPPORTIVE_ONLY','PENDING',
  'GEO design identifies six patients on two slides; no strict H&E pairing/registration or unchanged PF/control Z established.',
  'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327879','PENDING','PENDING','FAIL_COVID_NOT_PF','UNKNOWN','UNKNOWN')
 add('VENDOR_COSMX3','Bruker NSCLC FFPE example','Bruker CosMx NSCLC FFPE','CosMx','NSCLC',3,'UNKNOWN_HE_PAIR','FLUORESCENCE_MORPHOLOGY_NOT_HE_PROOF',
  'YES','YES','LABELS_AVAILABLE','UNKNOWN_EXACT_18','PENDING','NO_VERIFIED_REGISTERED_HE',
  'PUBLIC_VENDOR','SUPPORTIVE_ONLY','NOT_ADMITTED',
  'Three patients; repeated tissue samples are not independent donors. Below preferred five and no unchanged disease structure.',
  'https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/','PENDING','PENDING','FAIL_PF_CONTROL_ABSENT','UNKNOWN','UNKNOWN')
 add('VANNAN_DERIVED','CAMEO-Lung / HESCAPE / overlapping HEST entries','CAMEO-Lung','derived Xenium','PF/control',19,'DERIVED','DERIVED_REGISTERED_HE',
  'DERIVED','DERIVED','DERIVED','SAME_PANEL','DERIVED','INHERITS_DISCOVERY','PUBLIC','SUPPORTIVE_ONLY','EXCLUDED_DUPLICATE_DISCOVERY',
  'Derived current Vannan donors add zero independent information; do not add repository counts across releases.',
  'https://huggingface.co/datasets/theislab/CAMEO-Lung','DERIVED','DERIVED','DERIVED','DERIVED','DERIVED',False)
 add('SINGLE_DONOR_DEMO','10x single human lung demonstration','10x public lung examples','Xenium','demonstration',1,'VERSION_DEPENDENT','VERSION_DEPENDENT',
  'YES','YES','VERSION_DEPENDENT','VERSION_DEPENDENT','PENDING','DATASET_SPECIFIC','PUBLIC_VENDOR','SUPPORTIVE_ONLY','EXCLUDED_SINGLE_DONOR',
  'Single-donor demonstration cannot meet the independent donor gate; no count aggregation across panel demos without donor provenance.',
  'https://www.10xgenomics.com/datasets','PENDING','PENDING','UNKNOWN','PENDING','PENDING',False)
 registry=pd.DataFrame(entries);save(registry,'external_candidate_registry.tsv')
 transport=registry[['candidate_id','target_definitions','C_ontology','N_structure','Z_structure','morphology_150um','encoder_preprocessing','admission_status','reason','source_URL']].copy()
 transport['frozen_estimand_transport_PASS']='NO';save(transport,'transportability_matrix.tsv')
 first_panel=json.loads(gzip.decompress((S/'COPD_TMA1_gene_panel.json.gz').read_bytes()))
 present={x['type']['data'].get('name','') for x in first_panel['payload']['targets'] if x['type'].get('descriptor')=='gene'}
 missing=sorted(set(sum(GENES.values(),[]))-present)
 assert len(missing)==8
 summary=dict(audit_date='2026-10-03',verdict='HOLD_TRANSPORTABILITY',vannan_extension_new_donors=16,
  vannan_extension_author_registered_new_donors=13,vannan_extension_technically_admitted_new_donors=0,
  vannan_extension_verdict='VANNAN_EXTENSION_INSUFFICIENT',independent_external_study_sources_screened=int(registry.independent_study.sum()),
  admitted_external_cohorts=0,external_validation_authorized='NO',augmentation_authorized='NO',COPD_missing_frozen_genes=missing,
  no_models_run=True,no_full_candidate_datasets_downloaded=True)
 (O/'feasibility_complete.json').write_text(json.dumps(summary,indent=2))
 # HTTP errors or unsupported GEO targ=samples responses remain flagged artifacts.
 urls={
  'GSE250346_current.soft.txt':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346&targ=all&form=text&view=full',
  'GSE269354_Samples_list_xenium_metadata_v2.xlsx':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269354/suppl/GSE269354_Samples_list_xenium_metadata_v2.xlsx'}
 manifests=[]
 for f in sorted(S.iterdir()):
  if not f.is_file() or f.name.startswith('._'):continue
  url=urls.get(f.name,'')
  if f.name.startswith('PMC'):url='https://www.ebi.ac.uk/europepmc/webservices/rest/'+f.stem+'/fullTextXML'
  elif re.match(r'GSE\d+.*soft.txt',f.name) and not url:
   acc=re.search(r'GSE\d+',f.name).group();url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+acc+'&targ='+('samples' if '_samples' in f.name else 'self')+'&form=text&view=full'
  elif f.name.startswith('COPD_TMA'):
   i=re.search(r'TMA(\d)',f.name).group(1);gsm='GSM935887'+i
   url=f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9358nnn/{gsm}/suppl/{gsm}_COPD_TMA{i}_gene_panel.json.gz'
  elif f.name=='lungspatialdb_status.txt':url='https://www.lungspatialdb.com'
  status='UNSUPPORTED_GEO_QUERY_NOT_EVIDENCE' if '_samples.soft' in f.name else 'METADATA_ONLY'
  manifests.append(dict(file=str(f.relative_to(R)),source_URL=url,audit_date='2026-10-03',bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),status=status))
 save(pd.DataFrame(manifests),'source_metadata_manifest.tsv')
 print('PHASE1J_METADATA_AUDIT_COMPLETE',json.dumps(summary),flush=True)
if __name__=='__main__':main()
