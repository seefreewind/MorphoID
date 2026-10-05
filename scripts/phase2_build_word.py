"""Integrate frozen manuscript material into an editable Word review package."""
from pathlib import Path
import re,json,hashlib
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.section import WD_SECTION_START,WD_ORIENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
from PIL import Image
R=Path(__file__).resolve().parents[1];M=R/'manuscript';O=M/'MorphoID_Integrated_Manuscript.docx';D=Document()
for name in ['Normal','Title','Subtitle','Heading 1','Heading 2','Heading 3','Caption']:
 st=D.styles[name];st.font.name='Times New Roman';st.font.color.rgb=RGBColor(0,0,0);st.font.size=Pt(11 if name=='Normal' else {'Title':20,'Subtitle':11,'Heading 1':16,'Heading 2':13,'Heading 3':11,'Caption':9}[name]);st.paragraph_format.space_after=Pt(6)
D.styles['Normal'].paragraph_format.line_spacing=1.15
for name in ['Heading 1','Heading 2','Heading 3']:D.styles[name].paragraph_format.keep_with_next=True
D.styles['Caption'].font.bold=False
for st in D.styles:
 if st.type==1:
  pr=st.element.find(qn('w:pPr'))
  if pr is not None:
   for b in list(pr.findall(qn('w:pBdr'))):pr.remove(b)
  rp=st.element.find(qn('w:rPr'))
  if rp is not None:
   ff=rp.find(qn('w:rFonts'))
   if ff is not None:
    for attr in ['asciiTheme','hAnsiTheme','eastAsiaTheme','cstheme']:
     ff.attrib.pop(qn('w:'+attr),None)

sec=D.sections[0];sec.page_width=Inches(8.27);sec.page_height=Inches(11.69);sec.top_margin=sec.bottom_margin=Inches(.7);sec.left_margin=sec.right_margin=Inches(.8)
def footer(sec):
 p=sec.footer.paragraphs[0];p.alignment=WD_ALIGN_PARAGRAPH.CENTER;r=p.add_run('MorphoID  |  ');r.font.size=Pt(8)
 f=OxmlElement('w:fldSimple');f.set(qn('w:instr'),'PAGE');p._p.append(f)
footer(sec)
def section(land=False):
 s=D.add_section(WD_SECTION_START.NEW_PAGE);s.orientation=WD_ORIENT.LANDSCAPE if land else WD_ORIENT.PORTRAIT;s.page_width=Inches(11.69 if land else 8.27);s.page_height=Inches(8.27 if land else 11.69);s.left_margin=s.right_margin=Inches(.65 if land else .8);s.top_margin=s.bottom_margin=Inches(.65 if land else .7)
 return s
# Inline Markdown emphasis and clickable source links.
def inline(p,txt):
 pat=r'(\[[^\]]+\]\((?:<[^>]+>|[^)]+)\)|\*\*[^*]+\*\*|`[^`]+`|(?<!\*)\*[^*]+\*(?!\*))'
 for part in re.split(pat,txt):
  if not part:continue
  m=re.fullmatch(r'\[([^\]]+)\]\((?:<([^>]+)>|([^)]*))\)',part)
  if m:
   url=m.group(2) or m.group(3);h=OxmlElement('w:hyperlink');rid=p.part.relate_to(url,'http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink',is_external=True);h.set(qn('r:id'),rid);rr=OxmlElement('w:r');pr=OxmlElement('w:rPr');col=OxmlElement('w:color');col.set(qn('w:val'),'244A73');pr.append(col);rr.append(pr);tt=OxmlElement('w:t');tt.text=m.group(1);rr.append(tt);h.append(rr);p._p.append(h)
  else:
   rr=p.add_run(part.strip('*`') if part.startswith(('*','`')) else part)
   if part.startswith('**'):rr.bold=True
   elif part.startswith('*'):rr.italic=True
   elif part.startswith('`'):rr.font.name='Consolas';rr.font.size=Pt(8)
def para(txt,style=None):
 p=D.add_paragraph(style=style);inline(p,txt);return p
abbr={'epithelial_injury':'Epithelial','fibroblast_activation':'Fibroblast','macrophage_inflammatory':'Macrophage','DISEASE_ADJUSTED':'Adjusted','PF_ONLY':'PF only','POOLED':'Pooled','REFERENCE_BOOTSTRAP':'Reference','JACKKNIFE_NORMAL_CI':'Jackknife','BALANCED_DONOR_BOOTSTRAP':'Balanced','BCA_INFERENCE_SENSITIVITY':'BCa','EXACT_SIGN_TEST':'Sign test','WILCOXON_SENSITIVITY':'Wilcoxon'}
def clean(v):
 if v.startswith('repeat_'):
  v=v.split('_')[1]
 for a,b in abbr.items():v=v.replace(a,b)
 return v.replace('_',' ').replace('\\|','|')
def grid(headers,rows,small=False,transpose=False):
 if transpose:
  keys=[r[0] for r in rows];rows=[[h]+[r[j] for r in rows] for j,h in enumerate(headers[1:],1)];headers=[headers[0]]+keys
 # Split wide computational tables, preserving a stable row index in every column group.
 if len(headers)>11:
  keys=[j for j,h in enumerate(headers) if h in ['target','analysis_stratum','stratum','donor_id','run','method','effect','loss']][:3]
  fields=[j for j in range(len(headers)) if j not in keys]
  for k in range(0,len(fields),7):
   ix=keys+fields[k:k+7];cp=para('Column group '+str(k//7+1),'Caption');cp.paragraph_format.keep_with_next=True;grid(['Record']+[headers[j] for j in ix],[[str(n+1)]+[r[j] for j in ix] for n,r in enumerate(rows)],True)
  return
 t=D.add_table(rows=1,cols=len(headers));t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.autofit=False
 width=10.39 if D.sections[-1].orientation==WD_ORIENT.LANDSCAPE else 6.67
 widths=[width/len(headers)]*len(headers)
 if len(headers)==2:widths=[width*.32,width*.68]
 if len(headers)==4:widths=[width*.24]+[width*.76/3]*3
 for c,w in zip(t.columns,widths):c.width=Inches(w)
 for c,h in zip(t.rows[0].cells,headers):inline(c.paragraphs[0],clean(h))
 h=OxmlElement('w:tblHeader');t.rows[0]._tr.get_or_add_trPr().append(h)
 for vals in rows:
  cells=t.add_row().cells
  for c,v in zip(cells,vals):inline(c.paragraphs[0],clean(v))
 for ni,row in enumerate(t.rows):
  trpr=row._tr.get_or_add_trPr();n=OxmlElement('w:cantSplit');trpr.append(n)
  for ci,c in enumerate(row.cells):
   c.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER;pr=c._tc.get_or_add_tcPr();marg=OxmlElement('w:tcMar')
   for side in ['top','left','bottom','right']:
    x=OxmlElement('w:'+side);x.set(qn('w:w'),'65');x.set(qn('w:type'),'dxa');marg.append(x)
   pr.append(marg);borders=OxmlElement('w:tcBorders')
   for side in ['top','left','bottom','right']:
    x=OxmlElement('w:'+side);x.set(qn('w:val'),'single');x.set(qn('w:sz'),'4');x.set(qn('w:color'),'D9D9D9');borders.append(x)
   pr.append(borders)
   if ni==0 or ni%2==0:
    sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'E3E8ED' if ni==0 else 'F6F7F8');pr.append(sh)
   for p in c.paragraphs:
    p.paragraph_format.space_after=Pt(2);p.paragraph_format.space_before=Pt(2);p.paragraph_format.line_spacing=1.05;p.alignment=WD_ALIGN_PARAGRAPH.CENTER if small and ci>0 else WD_ALIGN_PARAGRAPH.LEFT
    for r in p.runs:r.font.size=Pt(8 if small else 9);r.bold=ni==0
 D.add_paragraph().paragraph_format.space_after=Pt(2)
legendtext=(M/'MorphoID_Figure_Legends.md').read_text();leg={}
for name,body in re.findall(r'^## (\S+)\n(.*?)(?=^## |\Z)',legendtext,re.M|re.S):leg[name]=body.split('![')[0].strip()
images=[]
def picture(path,caption=None,number=None):
 path=Path(str(path).strip('<>'));im=Image.open(path);w,h=im.size
 limitw=6.67;limith=6.5;scale=min(limitw/w,limith/h);p=D.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.paragraph_format.keep_with_next=True
 run=p.add_run();run.add_picture(str(path),width=Inches(w*scale));images.append(str(path))
 if caption:para(caption,'Caption')
 if number in leg:para(leg[number],'Caption')
# Main tables 2/3 transpose for readable numeric/categorical comparison.
def markdown(text,main=False,levelshift=0):
 lines=text.splitlines();n=0
 while n<len(lines):
  line=lines[n].strip();n+=1
  if not line:continue
  if line.startswith('```'):
   block=[]
   while n<len(lines) and not lines[n].startswith('```'):block.append(lines[n]);n+=1
   n+=1;para('\n'.join(block),'Normal');continue
  if line.startswith('|') and n<len(lines) and re.match(r'^\|[ :|\-]+$',lines[n].strip()):
   headers=[x.strip() for x in line.strip('|').split('|')];n+=1;rows=[]
   while n<len(lines) and lines[n].lstrip().startswith('|'):
    rows.append([x.strip() for x in re.split(r'(?<!\\)\|',lines[n].strip().strip('|'))]);n+=1
   grid(headers,rows,small=not main,transpose=main and len(headers)>7);continue
  match=re.fullmatch(r'!\[([^\]]*)\]\((?:<([^>]+)>|([^)]*))\)',line)
  if match:
   path=match.group(2) or match.group(3)
   if not path.startswith('/'):path=str(M/path)
   picture(path,match.group(1),Path(path).stem);continue
  hm=re.match(r'^(#{1,6})\s+(.*)',line)
  if hm:
   rank=len(hm.group(1));title=hm.group(2)
   if rank==1 and main:para(title,'Title')
   else:para(title,'Heading '+str(min(3,max(1,rank-1+levelshift))))
   continue
  if line.startswith('- '):para(line[2:],'List Bullet')
  elif re.match(r'^\d+\. ',line):para(line)
  else:para(line)
# Preserve the original main manuscript text and place six figures near their Results.
markdown((M/'MorphoID_Main_Manuscript.md').read_text(),main=True)
section();para('Supplementary methods','Title');markdown((M/'MorphoID_Supplementary_Methods.md').read_text().split('\n',1)[1])
section();para('Supplementary figures','Title')
for name in sorted(leg,key=lambda x: int(re.match(r'[FS](\d+)',x).group(1))):
 if not name.startswith('S'):continue
 if len(images)>6:D.add_page_break()
 para(name.replace('_',' '),'Heading 1');picture(M/'figures'/f'{name}.png',number=name)
# Supplementary numeric tables: retain every source row; display scientifically relevant
# nonduplicated fields, with direct full-source hyperlinks and hashes retained.
section(True);para('Supplementary evidence tables','Title')
para('All frozen rows are retained in the displayed evidence series. Computationally duplicated fields remain in the linked full-precision source TSVs. Target labels are abbreviated to Epithelial, Fibroblast and Macrophage; strata remain Pooled, Adjusted and PF only. NA denotes unavailable or inapplicable values. No model or estimate was changed.')
text=(M/'MorphoID_Supplementary_Tables.md').read_text();blocks=re.split(r'(?=^## )',text,flags=re.M)
selected={
'S5':['target','analysis_stratum','n_donors','n_sections','n_regions','M0_R2','M1_R2','M2_R2','M3_R2','M4_R2','deltaR2_morphology','CI_low','CI_high','FDR'],
'S8':['target','stratum','run','seed','M2_R2','M4_R2','DeltaR2','residual_R2','delta_CI_low','delta_CI_high','residual_CI_low','residual_CI_high','n_donors','n_regions'],
'S9':['target','stratum','n_partitions','median_DeltaR2','DeltaR2_p5','DeltaR2_p95','delta_fraction_positive','delta_fraction_negative','median_residual_R2','residual_p5','residual_p95','residual_fraction_positive'],
'S10':['target','stratum','donor_id','donor_not_in_stratum','full_DeltaR2','loo_DeltaR2','full_residual_R2','loo_residual_R2','delta_sign_flip','residual_sign_flip','classification_change','delta_FDR_classification_change','residual_FDR_classification_change','loo_phase1_verdict'],
'S11':['target','effect','loss','method','original_estimate','estimate','CI_low','CI_high','p_value','FDR','valid','classification','FDR_classification']}
for block in blocks:
 if not block.startswith('## '):continue
 key=block.split()[1]
 if key in selected:
  ls=block.splitlines();start=next(j for j,x in enumerate(ls) if x.startswith('|'));hdr=[x.strip() for x in ls[start].strip('|').split('|')];ix=[hdr.index(c) for c in selected[key] if c in hdr]
  new=['| '+' | '.join(hdr[j] for j in ix)+' |','| '+' | '.join(['---']*len(ix))+' |']
  for line in ls[start+2:]:
   if not line.startswith('|'):continue
   cells=[x.strip() for x in re.split(r'(?<!\\)\|',line.strip('|'))];new.append('| '+' | '.join(cells[j] for j in ix)+' |')
  block='\n'.join(ls[:start]+new)
 markdown(block)
section();para('Author review and submission materials','Title');markdown((M/'MorphoID_Cover_Letter_Draft.md').read_text().split('\n',1)[1])
for name,title,truncate in [('PHASE2_MANUSCRIPT_INTEGRATION_REPORT.md','Integration report',False),('FINAL_CLAIM_AUDIT.md','Claim calibration audit',True),('FINAL_REPRODUCIBILITY_MANIFEST.md','Reproducibility manifest',False),('JOURNAL_POSITIONING.md','Journal positioning',False),('PHASE2_FIGURE_QA.md','Figure quality audit',False)]:
 D.add_page_break();para(title,'Title');text=(R/'reports'/name).read_text()
 if truncate:text=text.split('## Every trigger-bearing manuscript line')[0]+f'\nFull line-by-line audit: [{name}](<{R/"reports"/name}>).\n'
 markdown(text)
D.core_properties.title='MorphoID integrated manuscript';D.core_properties.subject='Frozen results manuscript and supplementary evidence';D.core_properties.author='';D.core_properties.keywords='MorphoID; donor-held-out; spatial transcriptomics'
# Table and figure quantities are checked structurally; no scientific tests.
for par in D.paragraphs:
 pr=par._p.find(qn('w:pPr'))
 if pr is not None:
  for b in list(pr.findall(qn('w:pBdr'))):pr.remove(b)
D.save(O)
(M/'word_qa').mkdir(exist_ok=True)
(M/'word_qa/build_manifest.json').write_text(json.dumps({'output':str(O),'embedded_figures':len(images),'tables':len(D.tables),'sources':['MorphoID_Main_Manuscript.md','MorphoID_Supplementary_Methods.md','MorphoID_Supplementary_Tables.md','MorphoID_Figure_Legends.md','MorphoID_Cover_Letter_Draft.md'],'new_analysis':False,'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in M.glob('*.md') if not p.name.startswith('._')}},indent=2))
print(json.dumps({'output':str(O),'figures':len(images),'tables':len(D.tables)}))
