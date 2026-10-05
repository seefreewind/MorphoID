"""Retrieve bibliographic metadata only; no analysis data or model execution."""
import json,urllib.request,urllib.parse,concurrent.futures,difflib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
titles=[
'Visualization and analysis of gene expression in tissue sections by spatial transcriptomics',
'High-definition spatial transcriptomics for in situ tissue profiling',
'Slide-seq: A scalable technology for measuring genome-wide expression at high spatial resolution',
'Slide-seqV2: sensitive spatial genome-wide expression profiling at cellular resolution',
'High-spatial-resolution multi-omics sequencing via deterministic barcoding in tissue',
'Integrating spatial gene expression and breast tumour morphology via deep learning',
'A deep learning model to predict RNA-Seq expression of tumours from whole slide images',
'Spatial gene expression at single-cell resolution from histology using deep learning with GHIST',
'Systematic inference of super-resolution cell spatial profiles from histology images',
'A practical solution to pseudoreplication bias in single-cell studies',
'Confronting false discoveries in single-cell differential expression',
'Spatial transcriptomics identifies molecular niche dysregulation associated with distal lung remodeling in pulmonary fibrosis',
'Robust decomposition of cell type mixtures in spatial transcriptomics',
'Cell2location maps fine-grained cell types in spatial transcriptomics',
'Spatial deconvolution of heterogeneous tissue with RCTD',
'A single-cell atlas of the human healthy airways',
'Single-cell RNA-seq reveals ectopic and aberrant lung-resident cell populations in idiopathic pulmonary fibrosis',
'Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis',
'Alveolar fibroblast lineage orchestrates lung inflammation and fibrosis',
'Mapping spatially resolved transcriptomes in human and mouse pulmonary fibrosis',
'Persistence of a regeneration-associated, transitional alveolar epithelial cell state in pulmonary fibrosis',
'Inflammatory signals induce AT2 cell-derived damage-associated transient progenitors that mediate alveolar regeneration',
'Single-cell transcriptomics reveals that lung epithelial cells can enter a keratin 8-positive transitional state during alveolar regeneration',
'Monocyte-derived alveolar macrophages drive lung fibrosis and persist in the lung over the life span',
'A molecular cell atlas of the human lung from single-cell RNA sequencing',
'Integrated single-cell atlas of endothelial cells of the human lung',
'Single-cell RNA sequencing identifies diverse roles of epithelial cells in idiopathic pulmonary fibrosis',
'Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure',
'Using and understanding cross-validation strategies. Perspectives on Saeb et al.',
'Data leakage inflates prediction performance in connectome-based machine learning models',
'Batch effects in digital histopathology are not trivial',
'Quantifying the effects of data augmentation and stain color normalization in convolutional neural networks for computational pathology',
'A foundation model for generalizable cancer diagnosis and survival prediction from histopathological images',
'Towards a general-purpose foundation model for computational pathology',
'A visual-language foundation model for computational pathology',
'A whole-slide foundation model for digital pathology from real-world data',
'Prediction of molecular pathway activities from histopathology images using deep learning',
'A systematic comparison of algorithms for spatial transcriptomics',
'Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution',
'A new method for non-parametric multivariate analysis of variance',
'Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing',
'Better Bootstrap Confidence Intervals',
'Double/debiased machine learning for treatment and structural parameters',
'Genomic and transcriptomic correlates of immunotherapy response within the tumor microenvironment',
'EcoTyper: A framework for discovering and characterizing cell states and ecosystems',
]
def retrieve(item):
 i,t=item
 q='TITLE:"'+t+'"'
 url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':q,'format':'json','resultType':'core','pageSize':3})
 try:
  data=json.load(urllib.request.urlopen(url,timeout=35));rs=data.get('resultList',{}).get('result',[])
  if not rs:
   url='https://api.crossref.org/works?'+urllib.parse.urlencode({'query.title':t,'rows':2})
   c=json.load(urllib.request.urlopen(url,timeout=35))['message']['items']
   r=max(c,key=lambda x:difflib.SequenceMatcher(None,t.lower(),x.get('title',[''])[0].lower()).ratio())
   return {'id':i,'query':t,'title':r.get('title',[''])[0],'authors':'; '.join(x.get('family','') for x in r.get('author',[])[:6]),'year':r.get('published',r.get('issued',{})).get('date-parts',[[None]])[0][0],'journal':r.get('container-title',[''])[0],'doi':r.get('DOI'),'url':'https://doi.org/'+r.get('DOI',''),'match':difflib.SequenceMatcher(None,t.lower(),r.get('title',[''])[0].lower()).ratio(),'metadata_source':url}
  r=max(rs,key=lambda x:difflib.SequenceMatcher(None,t.lower(),x.get('title','').lower()).ratio())
  return {'id':i,'query':t,'title':r.get('title'),'authors':r.get('authorString'),'year':r.get('pubYear'),'journal':r.get('journalInfo',{}).get('journal',{}).get('title'),'doi':r.get('doi'),'url':'https://doi.org/'+r['doi'] if r.get('doi') else 'https://europepmc.org/article/'+r['source']+'/'+r['id'],'match':difflib.SequenceMatcher(None,t.lower(),r.get('title','').lower()).ratio(),'metadata_source':url}
 except Exception as e:return {'id':i,'query':t,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(retrieve,enumerate(titles,1)))
(ROOT/'results/final/literature_metadata.json').write_text(json.dumps(rows,indent=2))
for r in rows:print(r['id'],r.get('match'),r.get('title',r.get('error')),r.get('doi'))
