"""Stage989: replace all principal single-cell embedding UMAP displays with frozen full-data t-SNE.
Unsupervised PCA15 multiscale 30/300, cannot use label leakage to force islands.
"""
from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
EXP_CMAP=LinearSegmentedColormap.from_list('pale_cyan_to_deep_red',['#d9f4f7','#acdce9','#f4e6d5','#ea9c73','#b72b37','#670013'],N=256)
from matplotlib.lines import Line2D
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data"
O=P/"assets"/"stage989_fine26_tsne";O.mkdir(exist_ok=True)
z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz")
xy=np.load(Path(r"G:\Map6_recover_all\stage995_tsne\STAGE995_CURRENT_LOGZ_PCA20_TSNE_30_100.npz"))["embedding"]
assert xy.shape==(71950,2)
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0)
G=list(z["genes"].astype(str));A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["fine26_stable_id","section"])
K=pd.read_csv(D/"stage977_fine26_color_key.csv");palette=dict(zip(K.fine26,K.color));bmap=dict(zip(K.fine26,K.broad))
BC={"E":"#ec862f","I":"#397cbd","NE":"#2f9b75","ChAT":"#b75fb3","Other":"#999999"}
typ=A.fine26_stable_id.astype(str).to_numpy()
bro=np.array([bmap.get(t,"Other") for t in typ]);x,y=xy.T
def axis(ax):
 ax.set_aspect("equal",adjustable="box")
 ax.set_xticks([]);ax.set_yticks([])
 for spine in ax.spines.values():spine.set_visible(False)
def save(fig,name):
 fig.savefig(O/(name+".png"),dpi=190,bbox_inches="tight")
 fig.savefig(O/(name+".pdf"),bbox_inches="tight")
 plt.close(fig)
fig,ax=plt.subplots(figsize=(11.0,8.3))
rng=np.random.default_rng(989)
for t in K.fine26:
 inds=np.flatnonzero(typ==t);rng.shuffle(inds)
 ax.scatter(x[inds],y[inds],s=1.1,c=palette[t],alpha=.82,rasterized=True,linewidths=0)
axis(ax)
handles=[Line2D([0],[0],marker="o",linestyle="",markersize=6,color=palette[t],label=f"{t}: {(typ==t).sum():,}") for t in K.fine26]
fig.legend(handles=handles,loc="center right",bbox_to_anchor=(.985,.51),ncol=2,fontsize=8,frameon=False)
fig.suptitle("Fine26 on full-cell t-SNE · 26-class palette",fontsize=16,x=.40)
fig.subplots_adjust(left=.03,right=.73,top=.92,bottom=.05)
save(fig,"FINE26_CATEGORICAL_TSNE")
fig,ax=plt.subplots(figsize=(9,8))
for b in ["Other","I","E","NE","ChAT"]:
 m=bro==b;ax.scatter(x[m],y[m],s=.7,alpha=.70,c=BC[b],label=f"{b} ({m.sum():,})",rasterized=True,linewidths=0)
axis(ax);ax.legend(loc="lower center",bbox_to_anchor=(.5,-.06),frameon=False,ncol=5,fontsize=9)
fig.suptitle("Frozen broad molecular classes · t-SNE",fontsize=15)
save(fig,"BROAD_CLASSES_TSNE")
genes=[("Slc32a1","VGAT"),("Slc17a6","VGLUT2"),("Slc6a2","NE transporter"),("Slc5a7","Cholinergic transporter")]
total=np.maximum(C.sum(axis=1),1)
for normalized in [False,True]:
 fig,axes=plt.subplots(2,2,figsize=(12.5,11),layout="constrained")
 for ax,(gene,name) in zip(axes.flat,genes):
  v=np.log1p(C[:,G.index(gene)]) if not normalized else np.log1p(C[:,G.index(gene)]*10000/total)
  mm=v>0
  q=np.quantile(v[mm],.995) if mm.any() else 1
  ax.scatter(x[~mm],y[~mm],s=.25,c="#eff3f4",alpha=.20,linewidths=0,rasterized=True)
  p=ax.scatter(x[mm],y[mm],s=.70,c=v[mm],vmin=0,vmax=q,cmap=EXP_CMAP,alpha=.78,linewidths=0,rasterized=True)
  fig.colorbar(p,ax=ax,fraction=.037,pad=.02,shrink=.78)
  axis(ax);ax.set_title(f"{name} ({gene}) · {(mm).sum():,} detected\n{100*mm.mean():.1f}% of 71,950 cells",fontsize=11)
 fig.suptitle("Gene-specific "+("cell-depth-normalized" if normalized else "raw-corrected")+" counts on t-SNE",fontsize=14)
 save(fig,"VGAT_VGLUT2_NE_CHAT_MARKER_TSNE_"+("NORM" if normalized else "RAW"))
fig,axes=plt.subplots(2,2,figsize=(12.5,11),layout="constrained")
for ax,(gene,name) in zip(axes.flat,genes):
 v=np.log1p(C[:,G.index(gene)]*10000/total)
 th=float(np.quantile(v,.95));hi=v>=th
 ax.scatter(x[~hi],y[~hi],s=.25,c="#eff3f4",alpha=.20,linewidths=0,rasterized=True)
 sc=ax.scatter(x[hi],y[hi],s=2.0,c=v[hi],vmin=th,vmax=float(np.quantile(v,.999)),cmap=EXP_CMAP,alpha=.82,linewidths=0,rasterized=True)
 axis(ax);fig.colorbar(sc,ax=ax,fraction=.035,pad=.02)
 ax.set_title(f"{name} · top 5% within all 71,950\nThreshold log1p norm {th:.2f}",fontsize=11)
fig.suptitle("Relative marker enrichment · all 4 transmitters on identical t-SNE",fontsize=14)
save(fig,"VGAT_VGLUT2_NE_CHAT_TOP5PCT_TSNE")
# complete coexpression overlays for proposed thresholds and depth-percentile
fig,axs=plt.subplots(1,3,figsize=(17,6))
VG=C[:,G.index("Slc32a1")];VL=C[:,G.index("Slc17a6")]
clscols=["#4387BE","#E37B34","#A362BC","#BABABA"]
for ax,vcut in zip(axs,[2,3,10]):
 a=VG>10;b=VL>vcut
 labels=np.full(len(a),3,int);labels[a&~b]=0;labels[~a&b]=1;labels[a&b]=2
 for k in [3,2,0,1]:
  i=labels==k;ax.scatter(x[i],y[i],s=.6,color=clscols[k],alpha=.17 if k==3 else .8,linewidths=0,rasterized=True)
 axis(ax);ax.set_title(f"VGAT >10  /  VGLUT2 >{vcut} corrected counts",fontsize=11)
fig.legend(handles=[Line2D([0],[0],marker="o",linestyle="",color=c,markersize=8,label=n) for c,n in zip(clscols,["VGAT-only","VGLUT2-only","Double-positive","Neither"])],loc="lower center",bbox_to_anchor=(.5,.01),ncol=4,frameon=False)
fig.suptitle("Threshold-sensitive four-class calls on FULL-data t-SNE",fontsize=14)
fig.subplots_adjust(bottom=.13,top=.87,wspace=.08)
save(fig,"FOURCLASS_VGAT10_VGLUT2_2_3_10_TSNE")
(D/"stage989_tsne_display_authority.json").write_text(json.dumps({"stage":989,"status":"BUILT","n_cells":71950,"input":"Stage995 log-count gene-wise StandardScaler PCA20, unsupervised multiscale 30/100 t-SNE","cluster_truth":"Stage631 Fine26 frozen; no label used to train t-SNE","outputs":["FINE26_CATEGORICAL_TSNE","BROAD_CLASSES_TSNE","VGAT_VGLUT2_NE_CHAT_MARKER_TSNE_RAW","VGAT_VGLUT2_NE_CHAT_MARKER_TSNE_NORM","VGAT_VGLUT2_NE_CHAT_TOP5PCT_TSNE","FOURCLASS_VGAT10_VGLUT2_2_3_10_TSNE"]},indent=2))
print("STAGE989_TSNE_DISPLAY_COMPLETE",flush=True)
