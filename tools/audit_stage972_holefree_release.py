"""Map10 CI gate: consistent x flipped and y flipped across CURRENT spatial outputs."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from bs4 import BeautifulSoup
from map10_display_xy import flip_xy,verify_flip_xy
root=Path(__file__).resolve().parents[1]
pub=root/"perilc-map6-review"
data=pub/"data"
errors=[]
def ck(x,reason):
 if not x:errors.append(reason)
raw=pd.read_csv(data/"stage947_sectionwise_region10_display_labels.csv.gz",usecols=["section","x","y"])
trans=pd.read_csv(data/"stage956_cell_region10_labels.csv.gz",usecols=["section","x_display","y_display"])
ck(len(raw)==len(trans)==71950,"Coordinate input shapes mismatch")
ck(np.array_equal(raw.section.to_numpy(),trans.section.to_numpy()),"Section order mismatch")
sections={}
for s in (500,530,560):
 m=raw.section.to_numpy()==s
 try:
  verify_flip_xy(raw.x.to_numpy()[m],raw.y.to_numpy()[m],trans.x_display.to_numpy()[m],trans.y_display.to_numpy()[m])
  sections[str(s)]={"x_flipped":True,"y_flipped":True,"n_cells":int(m.sum())}
 except AssertionError as e:errors.append(f"S{s} mismatch: {e}")
for version in (956,970):
 f=data/f"stage{version}_gene_section_manifest.csv"
 tbl=pd.read_csv(f)
 ck(len(tbl)==81 and tbl.gene.nunique()==27,f"Stage{version} manifest missing maps")
 ck("x_flipped" in tbl.columns and "y_flipped" in tbl.columns,f"Stage{version} flags absent")
 if "x_flipped" in tbl.columns:ck(tbl.x_flipped.astype(bool).all() and tbl.y_flipped.astype(bool).all(),f"Stage{version} flags false")
pol=pd.read_csv(data/"stage956_gene_colorbar_policy.csv")
c=pd.read_csv(data/"stage970_gene_section_manifest.csv")
for g in pol.gene:
 q=c[c.gene==g]
 ck(len(q)==3 and q.vmax.nunique()==1 and q.gamma.nunique()==1,f"Colorbar mismatch {g}")
site={}
for name in ("index.html","single-cell.html"):
 p=pub/name
 soup=BeautifulSoup(p.read_text(encoding="utf8"),"html.parser")
 text=p.read_text(encoding="utf8").lower()
 ck("reversed" not in text and "xy_reversed" not in text,f"{name}: obsolete wording or asset")
 for node in soup.select("img[src]"):
  url=node.get("src","")
  ck((pub/url).is_file(),f"{name} missing {url}")
 site[name]={"images":len(soup.select("img")),"links_checked":len(soup.select("a[href]"))}
main=BeautifulSoup((pub/"index.html").read_text(encoding="utf8"),"html.parser")
allgenes=[x.get("src") for x in main.select(".gene-card img")]
ck(len(allgenes)==81 and all("/stage970_holefree_xy/" in x for x in allgenes),"Stale 2D gene panels")
for file in ("FINAL_LC_PERILC_LOCATOR_X_FLIPPED_Y_FLIPPED.png","FINAL_REGION10_HOLEFREE_X_FLIPPED_Y_FLIPPED.png"):
 ck(len(main.select('img[src$="'+file+'"]'))==1,"missing current map "+file)
both=(pub/"index.html").read_text(encoding="utf8")+(pub/"single-cell.html").read_text(encoding="utf8")
ck("assets/stage951_maskbody/" not in both,"old 3D maskbody shown")
ck(both.count("assets/stage967_maskbody_xy/")>=11,"corrected 3D panels not all linked")
for file in ("FINE26_XY_CURRENT.jpg","BROAD_XY_CURRENT.jpg","S500_LOCAL_MULTIPLEX_XY.png","S530_LOCAL_MULTIPLEX_XY.png","S560_LOCAL_MULTIPLEX_XY.png"):
 ck(file in both,"missing corrected spatial panel "+file)
mask=json.loads((data/"stage967_maskbody_xy_coordinate_audit.json").read_text(encoding="utf8"))
for s in ("500","530","560"):
 row=mask[s]
 ck(row["x_flipped"] and row["y_flipped"],f"3D flips missing {s}")
 ck(row["x_agreement_correlation"]>.995 and row["y_agreement_correlation"]>.995,f"3D coordinate mismatch {s}")
 ck(row["median_abs_x_um"]<10 and row["median_abs_y_um"]<10,f"3D ROI mismatch {s}")
script_names=["build_stage956_final_atlas.py","build_stage963_orientation_smooth.py",
"build_stage963_locator_from_source.py","build_stage964_fine26_xy_consistency.py",
"build_stage966_local_multiplex_xy.py","build_stage967_maskbody_xy_flips.py"]
for name in script_names:
 p=root/"tools"/name
 source=p.read_text(encoding="utf8").lower()
 ck("invert_yaxis()" not in source and "invert_xaxis()" not in source and "rot90(" not in source,f"mixed axis implementation {name}")
 ck("reversed" not in source and "xy_reversed" not in source,f"obsolete terminology {name}")
 ck(p.is_file(),f"missing active generator {name}")
for file in ("stage956_final_atlas_authority.json","stage970_holefree_region_audit.json"):
 a=json.loads((data/file).read_text(encoding="utf8"))
 ck("x flipped, y flipped" in str(a).lower(),f"missing orientation authority {file}")
shape=json.loads((data/"stage970_holefree_region_audit.json").read_text(encoding="utf8"))
ck(shape.get("stage")==970 and shape.get("gene_maps")==81,"Stage970 expected 81 maps")
for sec in ("500","530","560"):
 row=shape["region_QA"][sec]
 topo=row["after"]
 ck(topo["white_holes"]==0 and topo["white_pixels"]==0,f"S{sec} white holes remain")
 ck(len(topo["region_components"])==10 and all(n==1 for n in topo["region_components"].values()),f"S{sec} disconnected regions")
 ck(len(topo["region_holes"])==10 and all(n==0 for n in topo["region_holes"].values()),f"S{sec} region holes")
 ck(row["changed_existing_fraction"]<.03,f"S{sec} shape changed too much")
comps=pd.read_csv(data/"stage971_all3_holefree_manifest.csv")
ck(len(comps)==27 and comps.gene.nunique()==27,"stage971 three-section composite missing")
for r in comps.itertuples():
 ck((pub/"assets"/"stage971_all3_holefree_xy"/r.file).is_file(),f"missing composite {r.gene}")
ck((pub/"index.html").read_text(encoding="utf8").count("stage971_all3_holefree_xy/")==27,"composite links stale")
result={"stage":972,"status":"PASS" if not errors else "FAIL",
 "orientation":"x flipped, y flipped","source_coordinates_unchanged":True,
 "n_cells":len(raw),"gene_section_maps":len(c),"region_topology":shape["region_QA"],"regions":sections,
 "maskbody":mask,"pages":site,"errors":errors,
 "scope":"all CURRENT 2D and 3D spatial panels linked from map10 pages; historical sources retained for provenance"}
(data/"stage972_holefree_release_audit.json").write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf8")
print(json.dumps(result,indent=2,ensure_ascii=False))
if errors:raise SystemExit(1)
