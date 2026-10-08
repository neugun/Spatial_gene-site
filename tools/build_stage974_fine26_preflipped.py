"""Stage964: Fine26 and broad-class spatial panels in current atlas XY orientation."""
from pathlib import Path
import sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from map10_display_xy import flip_xy
root=Path(sys.argv[1])
A=pd.read_csv(root/"CURRENT_CELL_MEDIAL_PERILC_ANNOTATION.csv.gz")
P=pd.read_csv(root/"LC_CORE_SAMPLE_POLYGONS.csv")
pub=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
key=pd.read_csv(pub/"data"/"stage833_fine26_color_key.csv")
types=key.fine26.astype(str).tolist()
palette=dict(zip(key.fine26,key.color))
broad=dict(zip(key.fine26,key.broad))
broad_col={"E":"#f39c34","I":"#4aa3df","NE":"#48c774","ChAT":"#d36bc5","Other":"#b8b8b8"}
out=pub/"assets"/"stage974_preflipped_reference"
def setup(ax,s):
 q=A[A.section==s].copy()
 x=q.x_um.to_numpy(float);y=q.y_um.to_numpy(float)
 xx,yy=x.copy(),y.copy()
 q["x_show"]=xx;q["y_show"]=yy
 ax.set_facecolor("black")
 ax.scatter(xx,yy,s=.28,c="#222222",alpha=.38,linewidths=0,rasterized=True)
 pp=P[P.section==s].sort_values("vertex_order")
 if len(pp):
  xp,yp=pp.x_um.to_numpy(float),pp.y_um.to_numpy(float)
  ax.plot(xp,yp,c="#6cd4ff",lw=1.4)
 ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
 for sp in ax.spines.values():sp.set_visible(False)
 ax.set_title(f"S{s}",color="white",fontsize=12)
 return q
for mode in ("Fine26","Broad"):
 fig,axs=plt.subplots(1,3,figsize=(15.8,5.6),facecolor="black")
 for ax,s in zip(axs,(500,530,560)):
  q=setup(ax,s)
  q=q[q.primary_periLC_medial_450.astype(bool)].copy()
  q["fid"]=q.fine26_stable_id.astype(str)
  if mode=="Fine26":
   order=q.fid.value_counts().index.tolist()
   for t in order:
    w=q[q.fid==t]
    ax.scatter(w.x_show,w.y_show,s=1.8,c=palette.get(t,"#cccccc"),alpha=.90,linewidths=0,rasterized=True)
  else:
   q["bc"]=q.fid.map(broad).fillna("Other")
   for b in ("E","I","Other","ChAT","NE"):
    w=q[q.bc==b]
    ax.scatter(w.x_show,w.y_show,s=1.5,c=broad_col[b],alpha=.83,linewidths=0,rasterized=True)
 legend=types if mode=="Fine26" else ["E","I","Other","ChAT","NE"]
 colors=palette if mode=="Fine26" else broad_col
 handles=[plt.Line2D([0],[0],marker="o",linestyle="",markersize=5,markerfacecolor=colors[t],markeredgewidth=0,label=t) for t in legend]
 fig.legend(handles=handles,loc="lower center",ncol=13 if mode=="Fine26" else 5,frameon=False,labelcolor="white",fontsize=7,bbox_to_anchor=(.5,.005),columnspacing=.8,handletextpad=.25)
 fig.suptitle(("Fine26 populations" if mode=="Fine26" else "Broad molecular classes")+" in the medial periLC · already Stage631-flipped XY",color="white",fontsize=14,y=.98)
 fig.tight_layout(rect=[0,.08,1,.94])
 fname="FINE26_STAGE631_PREFLIPPED.jpg" if mode=="Fine26" else "BROAD_STAGE631_PREFLIPPED.jpg"
 fig.savefig(out/fname,dpi=160,facecolor="black",bbox_inches="tight")
 plt.close(fig)
 print("REBUILT",fname,flush=True)
