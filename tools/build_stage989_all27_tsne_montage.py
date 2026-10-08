"""Stage989b: complete 27-gene feature plots on one frozen t-SNE, not only 4 markers."""
from pathlib import Path
import numpy as np
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
O=P/"assets"/"stage989_fine26_tsne"
z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz")
xy=np.load(Path(r"G:\Map6_recover_all\stage982_tsne\multiscale_30_300.npz"))["embedding"]
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0).astype(float)
G=list(z["genes"].astype(str));n=len(C)
fig,axes=plt.subplots(6,5,figsize=(18,20),constrained_layout=True)
for j,ax in enumerate(axes.ravel()):
 if j>=len(G):ax.set_visible(False);continue
 raw=C[:,j];val=np.log1p(raw);pos=raw>0
 ax.scatter(xy[~pos,0],xy[~pos,1],s=.21,c="#dddddd",alpha=.25,linewidths=0,rasterized=True)
 vmax=float(np.quantile(val[pos],.99)) if pos.any() else 1.
 sc=ax.scatter(xy[pos,0],xy[pos,1],s=.55,c=val[pos],cmap="magma",vmin=0,vmax=max(.01,vmax),alpha=.84,linewidths=0,rasterized=True)
 ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
 for sp in ax.spines.values():sp.set_visible(False)
 ax.set_title(f"{G[j]} · {100*pos.mean():.0f}% detected",fontsize=10)
 cb=fig.colorbar(sc,ax=ax,fraction=.032,pad=.008,shrink=.63);cb.ax.tick_params(labelsize=7)
fig.suptitle("All 27 genes · same complete 71,950-cell t-SNE · gene-specific log(1+raw corrected count) colorbars",fontsize=17)
fig.savefig(O/"ALL27_GENE_TSNE_MONTAGE.png",dpi=165,bbox_inches="tight")
plt.close(fig)
print("STAGE989_ALL27_COMPLETE",len(G),n,flush=True)
