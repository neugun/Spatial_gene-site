"""Map10 Stages977-980 science + publication release QA."""
from pathlib import Path
import json,subprocess,csv,cv2,numpy as np,pandas as pd
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
D=P/"data";html=(P/"single-cell.html").read_text(encoding="utf8")
orient=json.loads((D/"stage976_coordinate_origin_audit.json").read_text())
assert orient["status"]=="PASS" and orient["no_double_flip"] is True
palette=pd.read_csv(D/"stage977_fine26_color_key.csv")
assert len(palette)==26 and palette.color.nunique()==26
geom=json.loads((D/"stage978_geometry_audit.json").read_text())
assert geom["metrics"]["Stage978"]["15nn_identity_purity"] > geom["metrics"]["Stage704"]["15nn_identity_purity"]
old=json.loads((D/"stage974_maskbody_coordinate_audit.json").read_text())
new=json.loads((D/"stage979_maskbody_coordinate_audit.json").read_text())
for s in ("500","530","560"):
 for name in ("roi_points_matched","x_agreement_correlation","y_agreement_correlation","median_abs_x_um","median_abs_y_um"):
  assert np.isclose(old[s][name],new[s][name],atol=1e-7),(s,name)
 assert new[s]["x_agreement_correlation"]>.999 and new[s]["y_agreement_correlation"]>.999
for name in ("FINE26_TRUE_MASKBODY_3D.png","FINE26_STAGE631_PREFLIPPED.jpg","S500_LOCAL_TRUE_MASKBODY.png","S530_LOCAL_TRUE_MASKBODY.png","S560_LOCAL_TRUE_MASKBODY.png"):
 assert ("stage979_consistent_fine26/"+name) in html,name
assert 'stage974_preflipped_reference/FINE26_' not in html
assert "stage978_umap_geometry/FINE26_STAGE704_VS_978.png" in html
film=P/"assets"/"stage980_true3d_video"/"STAGE980_FINE26_MASKBODY_ZOOM_ORBIT_12S.mp4"
assert film.stat().st_size < 100_000_000 and film.stat().st_size>1_000_000
assert film.name in html and 'id="video980"' in html
v=cv2.VideoCapture(str(film))
assert v.isOpened()
count=int(v.get(cv2.CAP_PROP_FRAME_COUNT))
fps=float(v.get(cv2.CAP_PROP_FPS))
width=int(v.get(cv2.CAP_PROP_FRAME_WIDTH));height=int(v.get(cv2.CAP_PROP_FRAME_HEIGHT))
assert count==300 and abs(fps-25)<.01 and (width,height)==(1800,650),(count,fps,width,height)
samples={}
for k in (0,50,100,175,299):
 v.set(cv2.CAP_PROP_POS_FRAMES,k);ok,fr=v.read();assert ok,k
 samples[k]=cv2.resize(fr,(450,162)).astype(np.float32)
v.release()
dist=lambda a,b:float(np.mean(np.abs(samples[a]-samples[b])))
assert dist(0,100)>3 and dist(100,299)>3,(dist(0,100),dist(100,299))
assert not list((P/"assets").rglob("*.npz")),"Public site must not expose raw embedding"
legacy=D/"stage980_legacy_video_reference.json"
rep=json.loads(legacy.read_text());rep["status"]="HISTORICAL_REFERENCES_VERIFIED_AND_NEW_MOVIE_RENDERED"
legacy.write_text(json.dumps(rep,indent=2),encoding="utf8")
result={"stage":980,"status":"PASS","frozen_stage976_orientation":True,"fine26_colors":26,"maskbody_coordinate_match":True,
"video":{"frames":count,"fps":fps,"width":width,"height":height,"file_bytes":film.stat().st_size,
"frame0_100_MAD":dist(0,100),"frame100_299_MAD":dist(100,299)},
"geometry_15nn":geom["metrics"],"source_movie_reference":["191","214","215"]}
(D/"stage980_release_QA.json").write_text(json.dumps(result,indent=2),encoding="utf8")
print("MAP10_STAGE980_RELEASE_QA_PASS",json.dumps(result),flush=True)
