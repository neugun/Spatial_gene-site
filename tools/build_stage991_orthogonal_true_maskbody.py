"""Stage991: actual Stage141 label-mask-body samples, orthogonal XY/XZ/YZ with XY dominant.
All sections included; no region gate. Physical x/y Stage631 already flipped (no double flip).
"""
from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
MC=Path(r"G:\PeriLC_current\D_mirror\11_PPTsummary\PeriLC_CHATGPT_RESULTS_20260914\141_DOWNSEG_MASKCLOUD")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
O=P/"assets"/"stage991_orthogonal_maskbody";O.mkdir(exist_ok=True);D=P/"data"
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["section","roi_id","fine26_stable_id","x_um","y_um","depth_um"])
K=pd.read_csv(D/"stage977_fine26_color_key.csv")
pal=dict(zip(K.fine26,K.color));bmap=dict(zip(K.fine26,K.broad))
BC={"E":"#ec862f","I":"#397cbd","NE":"#2f9b75","ChAT":"#b75fb3","Other":"#999999"}
rng=np.random.default_rng(991)
coords=[];summary=[]
for sec in (500,530,560):
 a=A[A.section==sec];q=np.load(MC/f"maskcloud_s{sec}_xy5_z7_cap45.npz");rid=q["roi_id"].astype(int)
 sz=max(int(rid.max()),int(a.roi_id.max()))+1
 lx=np.full(sz,np.nan);ly=np.full(sz,np.nan);colors=np.full(sz,-1,int);bl=np.full(sz,-1,int)
 lx[a.roi_id.to_numpy(int)]=a.x_um;ly[a.roi_id.to_numpy(int)]=a.y_um
 st=dict(zip(K.fine26,range(26)))
 colors[a.roi_id.to_numpy(int)]=[st.get(v,-1) for v in a.fine26_stable_id]
 B=["E","I","NE","ChAT","Other"]
 bl[a.roi_id.to_numpy(int)]=[B.index(bmap.get(v,"Other")) for v in a.fine26_stable_id]
 cx=q["x_um"].astype(float)/2;cy=q["y_flipped_um"].astype(float)/2;cz=q["z_um"].astype(float)/2
 ok=np.flatnonzero(np.isfinite(lx[rid])&np.isfinite(ly[rid])&(colors[rid]>=0))
 dx=np.median(lx[rid[ok]]+cx[ok]);dy=np.median(ly[rid[ok]]-cy[ok])
 xx=dx-cx;yy=dy+cy;zz=cz
 agree=(float(np.corrcoef(xx[ok],lx[rid[ok]])[0,1]),float(np.corrcoef(yy[ok],ly[rid[ok]])[0,1]))
 assert min(agree)>.995,agree
 # over 100k true mask voxels in every section, independent of centre identities.
 if len(ok)>200000:ok=np.sort(rng.choice(ok,200000,replace=False))
 coords.append((sec,xx[ok],yy[ok],zz[ok],colors[rid[ok]],bl[rid[ok]]))
 summary.append(dict(section=sec,n_cells=len(a),sampled_mask_voxels=len(ok),x_alignment=agree[0],y_alignment=agree[1],physical_x_extent=[float(xx.min()),float(xx.max())],physical_y_extent=[float(yy.min()),float(yy.max())],depth_extent=[float(zz.min()),float(zz.max())]))
 print("SECTION_MASK",sec,len(ok),agree,flush=True)
for mode in ("Fine26","Broad"):
 fig=plt.figure(figsize=(18,13),facecolor="white")
 gs=fig.add_gridspec(3,3,width_ratios=[2.2,1.0,1.0],wspace=.19,hspace=.33,left=.06,right=.985,bottom=.08,top=.92)
 for row,(sec,x,y,z,fine,broad) in enumerate(coords):
  cols=np.array([pal[K.fine26.iloc[int(t)]] for t in fine] if mode=="Fine26" else [list(BC.values())[int(t)] for t in broad])
  for col,(u,v,title,limx,limy) in enumerate([
    (x,y,"XY · primary coronal view",(-35,1065),(-35,1065)),
    (x,z,"XZ · depth profile",(-35,1065),(0,320)),
    (y,z,"YZ · depth profile",(-35,1065),(0,320))
   ]):
   ax=fig.add_subplot(gs[row,col])
   ax.scatter(u,v,c=cols,s=.22 if col==0 else .26,alpha=.73,linewidths=0,rasterized=True)
   ax.set_xlim(*limx);ax.set_ylim(*limy)
   ax.set_aspect("equal",adjustable="box")
   ax.set_xticks([]);ax.set_yticks([])
   ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
   ax.set_title(f"S{sec}  {title}",fontsize=11)
   if col==0:ax.set_ylabel("Y (flipped) · entire XY section",fontsize=9)
   elif col==1:ax.set_ylabel("Z depth",fontsize=9)
   else:ax.set_ylabel("Z depth",fontsize=9)
 fig.suptitle(f"Whole-section {mode} on authentic segmented-body voxel samples | XY, XZ, YZ orthogonal anatomy",fontsize=15)
 fig.savefig(O/f"STAGE991_{mode.upper()}_XY_XZ_YZ_TRUE_MASKBODY.png",dpi=205)
 fig.savefig(O/f"STAGE991_{mode.upper()}_XY_XZ_YZ_TRUE_MASKBODY.pdf")
 plt.close(fig);print("STAGE991_DONE",mode,flush=True)
(D/"stage991_orthogonal_authority.json").write_text(json.dumps({"stage":991,"status":"COMPLETE","section_rows":summary,"coordinate_frame":"Stage631 already flipped both XY; Stage141 source x sign corrected via ROI fit, same Stage979 transform","source":"Stage141 downseg.tif-derived true segmentation voxels, sampled, not reconstructed surfaces or cell-centroids","panels":"whole section XY primary, orthogonal XZ and YZ for depth"},indent=2))
