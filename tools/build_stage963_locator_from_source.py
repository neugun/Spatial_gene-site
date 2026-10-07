"""Stage963: plot LC polygons and periLC field in exactly the gene-atlas XY orientation.
Input annotation directory supplied by command line, never serialized to public outputs.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
raw=Path(sys.argv[1])
A=pd.read_csv(raw/"CURRENT_CELL_MEDIAL_PERILC_ANNOTATION.csv.gz")
P=pd.read_csv(raw/"LC_CORE_SAMPLE_POLYGONS.csv")
out=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review\assets\stage963_orientation_smooth")
out.mkdir(exist_ok=True)
ref=pd.read_csv(Path(r"G:\Spatial_gene_site_publish\perilc-map6-review\data\stage947_sectionwise_region10_display_labels.csv.gz"))
fig,axs=plt.subplots(1,3,figsize=(13.2,4.9))
for ax,s in zip(axs,(500,530,560)):
 q=A[A.section==s]
 t=ref[ref.section==s]
 assert len(q)==len(t)
 xr=q.x_um.to_numpy(float);yr=q.y_um.to_numpy(float)
 # Verify source annotation coordinate bounds correspond to gene-atlas bounds
 assert abs(xr.min()-t.x.min())<1 and abs(xr.max()-t.x.max())<1
 assert abs(yr.min()-t.y.min())<1 and abs(yr.max()-t.y.max())<1
 xx=xr.min()+xr.max()-xr
 yy=yr.min()+yr.max()-yr
 ax.scatter(xx,yy,s=.35,c="#d7d7d7",alpha=.65,linewidths=0,rasterized=True)
 flag=q.medial_wrap_450um.astype(bool).to_numpy()
 gx=np.linspace(xx.min(),xx.max(),220);gy=np.linspace(yy.min(),yy.max(),220)
 XX,YY=np.meshgrid(gx,gy)
 tree=cKDTree(np.c_[xx,yy])
 dist,idx=tree.query(np.c_[XX.ravel(),YY.ravel()],k=1)
 F=gaussian_filter(flag[idx].reshape(XX.shape).astype(float),2.2)
 Z=np.where(dist.reshape(XX.shape)<35,F,np.nan)
 ax.contour(gx,gy,Z,levels=[.5],colors=["#d97706"],linewidths=1.8)
 pp=P[P.section==s].sort_values("vertex_order")
 if len(pp):
  px=xr.min()+xr.max()-pp.x_um.to_numpy(float)
  py=yr.min()+yr.max()-pp.y_um.to_numpy(float)
  ax.plot(px,py,c="#2563eb",lw=1.9)
  ax.text(px.mean(),py.mean(),"LC",fontsize=9,ha="center",va="center",color="#1d4ed8",fontweight="bold")
 if flag.any():ax.text(xx[flag].mean(),yy[flag].mean(),"periLC",fontsize=9,ha="center",va="center",color="#b45309",fontweight="bold")
 ax.set_title(f"S{s}",fontsize=11)
 ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
 for sp in ax.spines.values():sp.set_visible(False)
fig.suptitle("LC / periLC locator · XY jointly reversed from source",fontsize=14)
fig.tight_layout(rect=[0,0,1,.95])
fig.savefig(out/"FINAL_LC_PERILC_LOCATOR_XY_REVERSED.png",dpi=190,bbox_inches="tight")
plt.close(fig)
print("LOCATOR_FROM_SOURCE_PASSED",len(A))
