"""Stage976 independent check against original acquisition frame and Stage393/838."""
from pathlib import Path
import json,numpy as np,pandas as pd
from bs4 import BeautifulSoup
root=Path(r"G:\Spatial_gene_site_publish")
pub=root/"perilc-map6-review";D=pub/"data"
source=Path(r"G:\Map6_recover_all\N5_ARRAY_TO_STAGE631_AFFINE.csv")
aff=pd.read_csv(source)
raw=pd.read_csv(D/"stage947_sectionwise_region10_display_labels.csv.gz")
atlas=pd.read_csv(Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926\stage631_current_3d_atlas\CURRENT_3D_CELL_ATLAS.csv.gz"),usecols=["section","x_um","y_um"])
assert len(raw)==len(atlas)==71950
assert np.array_equal(raw.section,atlas.section)
assert np.allclose(raw.x,atlas.x_um,rtol=0,atol=.001)
assert np.allclose(raw.y,atlas.y_um,rtol=0,atol=.001)
lineage={}
for sec in (500,530,560):
 q=aff[(aff.section==sec)&aff.target.isin(["x_um","y_um"])]
 assert len(q)==2
 x=q[q.target=="x_um"].iloc[0];y=q[q.target=="y_um"].iloc[0]
 assert abs(x.coef_array_x+.46)<1e-6 and abs(y.coef_array_y+.46)<1e-6
 assert abs(x.coef_array_y)<1e-8 and abs(y.coef_array_x)<1e-8
 assert x.rmse<1e-5 and y.rmse<1e-5
 lineage[str(sec)]={"x_um_per_raw_array_x_px":float(x.coef_array_x),"y_um_per_raw_array_y_px":float(y.coef_array_y),
                   "x_intercept_um":float(x.intercept),"y_intercept_um":float(y.intercept),
                   "x_flipped_from_raw_array":True,"y_flipped_from_raw_array":True,"extra_xy_flip":False,
                   "cells":int((atlas.section==sec).sum())}
man=pd.read_csv(D/"stage973_gene_section_manifest.csv")
prior=pd.read_csv(D/"stage970_gene_section_manifest.csv")
assert len(man)==81 and man.gene.nunique()==27 and man.section.nunique()==3
assert man.x_flipped.all() and man.y_flipped.all() and not man.extra_flip.any()
compare=man.merge(prior,on=["gene","section"],validate="one_to_one",suffixes=("_new","_prior"))
assert np.allclose(compare.vmax_new,compare.vmax_prior)
assert np.allclose(compare.gamma_new,compare.gamma_prior)
shape=json.loads((D/"stage973_correct_orientation_region_audit.json").read_text(encoding="utf8"))
assert shape["stage"]==973 and shape["gene_maps"]==81
for sec in (500,530,560):
 t=shape["region_QA"][str(sec)]["after"]
 assert t["white_holes"]==0 and t["white_pixels"]==0
 assert all(v==1 for v in t["region_components"].values())
 assert all(v==0 for v in t["region_holes"].values())
ref=json.loads((D/"stage974_maskbody_coordinate_audit.json").read_text(encoding="utf8"))
for s in ["500","530","560"]:
 assert ref[s]["x_agreement_correlation"]>.999 and ref[s]["y_agreement_correlation"]>.999
html=(pub/"index.html").read_text(encoding="utf8")
sc=(pub/"single-cell.html").read_text(encoding="utf8")
assert html.count("assets/stage973_correct_orientation_roundedge/")>=164
assert html.count("assets/stage975_all3_correct_orientation/")==27
assert html.count("assets/stage974_preflipped_reference/")>=3
assert sc.count("assets/stage974_preflipped_reference/")>=15
for name,doc in [("index.html",html),("single-cell.html",sc)]:
 soup=BeautifulSoup(doc,"html.parser")
 for item in soup.select("img[src]"):
  assert (pub/item["src"]).is_file(),(name,item["src"])
 for item in soup.select("a[href]"):
  link=item["href"]
  if link.startswith("assets/") or link.startswith("data/"):
   assert (pub/link).is_file(),(name,link)
record={"stage":976,"status":"PASS","authority":"Stage393 / Stage838 current pre-flipped physical XY, verified by Stage631 affine from raw acquisition array",
 "input_coordinate_columns":["stage631.x_um","stage631.y_um"],"no_double_flip":True,
 "gene_maps":81,"composites":27,"cell_universe":71950,
 "raw_array_affine":lineage,"colorbar_unmodified":True,"region_topology_no_holes":True,
 "stage974_maskbody_roi_alignment":ref,
 "note":"All figure generation uses Stage631 cell coordinates and frozen Stage943 gene values; region visual mask is a remapping of frozen Stage956 categorical labels, not a raster-image flip."}
(D/"stage976_coordinate_origin_audit.json").write_text(json.dumps(record,indent=2),encoding="utf8")
print("STAGE976_COORDINATE_ORIGIN_PASS",json.dumps({"cells":71950,"section_affine":lineage,"maps":81,"composites":27}),flush=True)
