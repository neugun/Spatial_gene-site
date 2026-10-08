from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz")
full=np.load(Path(r"G:\Map6_recover_all\stage982_tsne\multiscale_30_300.npz"))["embedding"]
assert full.shape==(71950,2)
genes=list(z["genes"].astype(str))
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0)
va=C[:,genes.index("Slc32a1")]>10;ve=C[:,genes.index("Slc17a6")]>10
q=np.full(len(va),3,int);q[va&~ve]=0;q[~va&ve]=1;q[va&ve]=2
colors=["#4a7db6","#ec892b","#a45ac2","#bdbdbd"]
names=["VGAT-only","VGLUT2-only","Co-positive","Neither"]
fig,axs=plt.subplots(1,2,figsize=(14,6))
for ax,(lab,X) in zip(axs,[("Current UMAP",z["embedding"]),("Multiscale t-SNE",full)]):
 for k in (3,2,0,1):
  m=q==k
  ax.scatter(X[m,0],X[m,1],s=.6,c=colors[k],alpha=.12 if k==3 else .80,linewidths=0,rasterized=True)
 ax.set_title(lab+" · same >10 raw-count rule",fontsize=12);ax.set_xticks([]);ax.set_yticks([])
 ax.set_aspect("equal",adjustable="box")
 for sp in ax.spines.values():sp.set_visible(False)
handles=[plt.Line2D([0],[0],marker="o",linestyle="",color=c,markersize=7,label=f"{nm}: {(q==k).sum():,} ({(q==k).mean()*100:.1f}%)") for k,(c,nm) in enumerate(zip(colors,names))]
fig.legend(handles=handles,loc="lower center",ncol=4,frameon=False,fontsize=10)
fig.suptitle("All cells, paper-faithful >10 spots/cell VGAT/VGLUT2 four-class overlays\nNo artificial co-positive suppression",fontsize=14)
fig.tight_layout(rect=[0,.07,1,.94])
out=P/"assets"/"stage982_tsne_comparison"
fig.savefig(out/"STAGE982_PAPER_GT10_ALLCELL_FOURCLASS_UMAP_TSNE.png",dpi=195,bbox_inches="tight")
fig.savefig(out/"STAGE982_PAPER_GT10_ALLCELL_FOURCLASS_UMAP_TSNE.pdf",bbox_inches="tight")
plt.close(fig)
(P/"data"/"stage982_paper_gt10_fourclass_counts.csv").write_text(pd.DataFrame({"category":names,"n":[int((q==k).sum()) for k in range(4)],"percent":[float(100*(q==k).mean()) for k in range(4)]}).to_csv(index=False))
print("STAGE982_PAPER_GT10",dict(zip(names,np.bincount(q,minlength=4).tolist())),flush=True)
