from pathlib import Path
import json,urllib.request,urllib.parse,time
R=Path(__file__).resolve().parents[1];p=R/'results/final/literature_metadata.json';rows=json.loads(p.read_text())
fix={3:'10.1126/science.aaw1219',4:'10.1038/s41587-020-0739-1',8:'10.1038/s41592-025-02795-z',11:'10.1038/s41467-021-25960-2',13:'10.1038/s41587-021-00830-w',15:'10.1038/s41592-021-01264-7',16:'10.1164/rccm.201911-2199OC',23:'10.1038/s41467-020-17358-3',31:'10.1038/s41467-021-24698-1',33:'10.1038/s41586-024-07894-z',37:'10.1038/s43018-020-0085-8',38:'10.1038/s41467-025-56618-y',40:'10.1214/aos/1013699998',42:'10.1080/01621459.1987.10478410',43:'10.1111/ectj.12097',44:'10.1038/s41592-019-0667-5',45:'10.1038/s41587-019-0114-2',46:'10.1093/bib/bbac297'}
for ident,doi in fix.items():
 url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':'DOI:'+doi,'format':'json','resultType':'core','pageSize':1})
 try:
  rs=json.load(urllib.request.urlopen(url,timeout=25))['resultList']['result']
  if rs:
   x=rs[0];r={'id':ident,'title':x['title'],'authors':x.get('authorString'),'year':x.get('pubYear'),'journal':x.get('journalInfo',{}).get('journal',{}).get('title'),'doi':doi,'url':'https://doi.org/'+doi,'match':1,'metadata_source':url}
  else:
   url='https://api.crossref.org/works/'+urllib.parse.quote(doi,safe='');x=json.load(urllib.request.urlopen(url,timeout=25))['message'];r={'id':ident,'title':x['title'][0],'authors':'; '.join(a.get('family','') for a in x.get('author',[])[:6]),'year':x.get('published',x['issued'])['date-parts'][0][0],'journal':x.get('container-title',[''])[0],'doi':doi,'url':'https://doi.org/'+doi,'match':1,'metadata_source':url}
  rows=[r if x['id']==ident else x for x in rows]
  if ident>45:rows.append(r)
  print(ident,r['title'])
 except Exception as e:print('ERROR',ident,doi,e)
 time.sleep(.3)
for ident,title,url,yr in [(47,'HEST-1k: A Dataset For Spatial Transcriptomics and Histology Image Analysis','https://proceedings.neurips.cc/paper_files/paper/2024/hash/60a899cc31f763be0bde781a75e04458-Abstract-Datasets_and_Benchmarks_Track.html',2024),(48,'GHIST+ official implementation: tissue-wide reconstruction of single-cell molecular states','https://github.com/SydneyBioX/GHIST_plus',2026),(49,'Phikon-v2, A large and public feature extractor for biomarker prediction','https://arxiv.org/abs/2409.09173',2024),(50,'STP-BENCH: official benchmark release','https://github.com/NEXGEM/STP-Bench',2026)]:rows.append({'id':ident,'title':title,'authors':'See primary source','year':yr,'journal':'Proceedings / preprint / software record (not peer-reviewed journal evidence)','url':url,'match':1,'metadata_source':'Official primary record verified with web'})
p.write_text(json.dumps(sorted(rows,key=lambda x:x['id']),indent=2))
