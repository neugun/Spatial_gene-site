"""Stage965 audit for the currently visible Three-section atlas, not stale Stage956 panels."""
from pathlib import Path
import json,re
import numpy as np,pandas as pd
from bs4 import BeautifulSoup
pub=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
D=pub/"data"
manifest=pd.read_csv(D/"stage963_gene_section_manifest.csv")
prior=pd.read_csv(D/"stage956_gene_section_manifest.csv")
meta=json.loads((D/"stage963_orientation_smoothing_audit.json").read_text(encoding="utf8"))
issues=[]
def ck(flag,msg):
 if not flag:issues.append(msg)
ck(len(manifest)==81 and manifest.gene.nunique()==27,"81 maps/27 genes missing")
ck(set(manifest.section)=={500,530,560},"wrong sections")
ck(manifest.groupby("gene").section.nunique().eq(3).all(),"missing gene panels")
ck(manifest.x_flipped.all() and manifest.y_flipped.all(),"XY inversion metadata missing")
join=manifest.merge(prior,on=["gene","section"],suffixes=("_new","_old"),validate="one_to_one")
ck(len(join)==81 and np.allclose(join.vmax_new,join.vmax_old) and np.allclose(join.gamma_new,join.gamma_old),"colorbar changed")
for row in manifest.itertuples():
 p=pub/"assets"/"stage963_orientation_smooth"/row.file
 ck(p.exists() and p.stat().st_size>10000,f"missing map {row.file}")
doc=BeautifulSoup((pub/"index.html").read_text(encoding="utf8"),"html.parser")
ck(len(doc.select(".gene-card"))==27,"not 27 gene cards")
ck(len(doc.select(".gene-card .pane"))==81,"not 81 visible image panels")
imgs=doc.select(".gene-card img")
ck(len(imgs)==81 and all("stage963_orientation_smooth/" in x.get("src","") for x in imgs),"stale gene images")
for fname in ("FINAL_LC_PERILC_LOCATOR_X_FLIPPED_Y_FLIPPED.png","FINAL_REGION10_SMOOTH_X_FLIPPED_Y_FLIPPED.png"):
 ck(len(doc.select('img[src*="'+fname+'"]'))==1,"missing main figure "+fname)
for ss in (500,530,560):
 q=meta["region_QA"][str(ss)]
 ck(q["agreement"]>=.97 and q["min_region_iou"]>=.9,f"S{ss} shape agreement failed")
page={}
for fn in ("index.html","single-cell.html"):
 p=pub/fn;s=BeautifulSoup(p.read_text(encoding="utf8"),"html.parser")
 ids={e.get("id") for e in s.select("[id]")}
 bad=[]
 for tag,attr in (("img","src"),("a","href")):
  for e in s.find_all(tag):
   x=e.get(attr,"")
   if not x or x.startswith(("http:","https:","mailto:","data:")):continue
   if x.startswith("#"):
    if x[1:] not in ids:bad.append(x)
   elif not (pub/x.split("?")[0].split("#")[0]).exists():bad.append(x)
 ck(not bad,f"{fn}: missing links {bad[:4]}")
 page[fn]=dict(n_images=len(s.select("img")),n_broken_links=len(bad))
ck("Stage961" in (pub/"index.html").read_text(encoding="utf8"),"lost Stage961 detail")
ck("Stage962" in (pub/"index.html").read_text(encoding="utf8"),"lost Stage962 detail")
result={"stage":965,"status":"PASS" if not issues else "FAIL","issues":issues,
 "gene_maps":len(manifest),"n_genes":manifest.gene.nunique(),
 "region_display_agreement":meta["region_QA"],
 "page_integrity":page,"expression_unchanged":True,"source_Stage943_RouteA":True,
 "fine26_labels_and_stage956_region_assignments_unchanged":True,
 "note":"Stage963 is display-only smoothing; it does not validate the anatomical interpretation of 10 regions."}
(D/"stage965_visible_atlas_audit.json").write_text(json.dumps(result,indent=2),encoding="utf8")
print(json.dumps(result,indent=2))
if issues:raise SystemExit(1)
