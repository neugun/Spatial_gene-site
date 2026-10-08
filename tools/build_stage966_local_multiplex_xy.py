"""Stage966 local Fine26 multiplex, corrected display XY only; selected windows frozen."""
from pathlib import Path
import sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from map10_display_xy import flip_xy
root=Path(sys.argv[1])
pub=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
D=pub/"data";out=pub/"assets"/"stage963_orientation_smooth"
A=pd.read_csv(root/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["section","x_um","y_um","fine26","fine26_stable_id"])
C=np.load(D/"stage943_routeA_all27.npy",mmap_mode="r")
genes=pd.read_csv(D/"stage956_gene_section_manifest.csv").gene.drop_duplicates().tolist()
wins=pd.read_csv(D/"stage950_local_singlecell_windows.csv")
assert len(A)==len(C)==71950 and C.shape[1]==len(genes)==27
sec=A.section.to_numpy(int)
cmap=plt.get_cmap("turbo")
for s in (500,530,560):
 ss=A.loc[sec==s].reset_index(drop=True);cs=np.asarray(C[sec==s],dtype=float)
 info=wins[wins.section==s].iloc[0]
 x0=float(info.x0);y0=float(info.y0);size=float(info.window_um)
 selected=(ss.x_um>=x0)&(ss.x_um<x0+size)&(ss.y_um>=y0)&(ss.y_um<y0+size)
 assert int(selected.sum())==int(info.n_cells)
 W=ss[selected].copy();CW=cs[selected.to_numpy()]
 burden=CW.sum(1);order=np.argsort(-burden)[:min(96,len(burden))]
 SW=W.iloc[order].copy();SC=CW[order]
 z=np.log1p(SC);z=(z-z.mean(0))/(z.std(0)+1e-9)
 order2=np.argsort(SW.fine26_stable_id.astype(str).to_numpy())
 z=z[order2]
 fig=plt.figure(figsize=(12,6.4),facecolor="black")
 gs=fig.add_gridspec(1,2,width_ratios=[.9,1.1],wspace=.24)
 ax=fig.add_subplot(gs[0,0]);ax.set_facecolor("black")
 xr=ss.x_um.to_numpy(float);yr=ss.y_um.to_numpy(float)
 originx=xr.min()+xr.max();originy=yr.min()+yr.max()
 xshow,yshow=flip_xy(W.x_um.to_numpy(float),W.y_um.to_numpy(float),x_bounds=(xr.min(),xr.max()),y_bounds=(yr.min(),yr.max()))
 for typ in sorted(W.fine26.unique()):
  sel=W.fine26.to_numpy()==typ
  ax.scatter(xshow[sel],yshow[sel],s=3.5,color=cmap(typ/25),alpha=.8,linewidths=0,rasterized=True)
 ax.set_xlim(originx-(x0+size),originx-x0)
 ax.set_ylim(originy-(y0+size),originy-y0)
 ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
 ax.set_title(f"S{s} local field | {len(W)} cells, {W.fine26.nunique()} Fine26",color="white",fontsize=10)
 ax=fig.add_subplot(gs[0,1]);ax.set_facecolor("black")
 ax.imshow(z.T,aspect="auto",cmap="coolwarm",vmin=-2.5,vmax=2.5)
 ax.set_yticks(range(len(genes)),genes,fontsize=6,color="white");ax.set_xticks([])
 ax.set_title("Selected cells x 27 genes",color="white",fontsize=10)
 fig.suptitle(f"Stage943 local Fine26 multiplex | x flipped, y flipped | S{s}",color="white",fontsize=13)
 fig.tight_layout()
 fig.savefig(out/f"S{s}_LOCAL_MULTIPLEX_XY.png",dpi=190,bbox_inches="tight",facecolor="black")
 plt.close(fig)
 print("STAGE966_LOCAL_OK",s,len(W),flush=True)
