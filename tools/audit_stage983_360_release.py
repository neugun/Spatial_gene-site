from pathlib import Path
import numpy as np,cv2,json,pandas as pd
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
D=P/"data";A=P/"assets"/"stage983_true3d_360_video"
v=A/"STAGE983_FINE26_MASKBODY_ZOOM_360_ZOOMOUT_27S.mp4"
cap=cv2.VideoCapture(str(v));assert cap.isOpened()
n=int(cap.get(cv2.CAP_PROP_FRAME_COUNT));fps=float(cap.get(cv2.CAP_PROP_FPS))
w=int(cap.get(cv2.CAP_PROP_FRAME_WIDTH));h=int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
assert (n,fps,w,h)==(675,25.,1650,600),(n,fps,w,h)
images={}
for k in [0,100,205,310,415,535,674]:
 cap.set(cv2.CAP_PROP_POS_FRAMES,k);ok,frame=cap.read();assert ok
 images[k]=cv2.resize(frame,(550,200))
cap.release()
mad=lambda x,y:float(np.mean(np.abs(images[x].astype(float)-images[y].astype(float))))
assert mad(100,310)>5 and mad(310,535)>5 and mad(0,205)>5
assert (A/"STAGE983_SEVEN_FRAME_QC.jpg").is_file()
ori=json.loads((D/"stage976_coordinate_origin_audit.json").read_text())
assert ori["status"]=="PASS" and ori["no_double_flip"]
s=(P/"single-cell.html").read_text(encoding="utf8")
for text in ["id=\"threshold981\"","STAGE983_FINE26_MASKBODY_ZOOM_360_ZOOMOUT_27S.mp4","stage984_h5ad_vs_routeA_source_audit.csv"]:
 assert text in s
X=pd.read_csv(D/"stage981_vgat_vglut2_four_class_by_section.csv")
assert len(X)==54
D.joinpath("stage983_360_video_QA.json").write_text(json.dumps({"stage":983,"status":"PASS","frames":n,"fps":fps,"duration_seconds":n/fps,"dimensions":[w,h],"bytes":v.stat().st_size,"motion_MAD":[mad(100,310),mad(310,535),mad(0,205)],"orientation_stage976":"PASS; no double flip","source":"Stage141 true maskcloud -> Stage631 matched ROI -> 26 Fine26 colors","camera":"Full 360-degree smooth orbit; zoom in then zoom out"},indent=2))
print("STAGE983_QA_PASS",n,fps,"bytes",v.stat().st_size,flush=True)


