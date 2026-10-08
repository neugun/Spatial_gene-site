"""Map10 Stage977b: relative marker enrichment versus frozen broad classes."""
from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
O=P/"assets"/"stage977_fine26_umap"
z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz")
E=z["embedding"];genes=list(z["genes"].astype(str))
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0)
assert C.shape==(len(E),27)
a=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["fine26_stable_id"])
key=pd.read_csv(P/"data"/"stage977_fine26_color_key.csv")
broad=a.fine26_stable_id.map(dict(zip(key.fine26,key.broad))).to_numpy()
assert pd.notna(broad).all()
geneset=[("Slc32a1","VGAT"),("Slc17a6","VGLUT2"),("Slc6a2","NE transporter"),("Slc5a7","Cholinergic transporter")]
total=np.maximum(C.sum(axis=1),1)
fig,axs=plt.subplots(2,2,figsize=(12.5,11),constrained_layout=True)
records=[]
for ax,(g,title) in zip(axs.ravel(),geneset):
 v=np.log1p(1e4*C[:,genes.index(g)]/total)
 th=float(np.quantile(v,.95))
 high=v>=th
 assert .04<high.mean()<.08,(g,high.mean())
 ax.scatter(E[~high,0],E[~high,1],s=.32,c="#d9d9d9",alpha=.20,linewidths=0,rasterized=True)
 sc=ax.scatter(E[high,0],E[high,1],c=v[high],s=2.5,cmap="magma",vmin=th,vmax=float(np.quantile(v,.999)),alpha=.9,linewidths=0,rasterized=True)
 ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect("equal")
 for sp in ax.spines.values():sp.set_visible(False)
 ax.set_title(f"{title} · {g}\nTop 5% of cell-depth-normalized counts (n={high.sum():,})",fontsize=12)
 cb=fig.colorbar(sc,ax=ax,fraction=.04,pad=.025,shrink=.77);cb.set_label("log(1 + 10k × count / all 27 counts)",fontsize=8)
 for b in ("E","I","NE","ChAT","Other"):
  m=(broad==b)
  p=float(np.mean(high[m])); baseline=float(np.mean(high))
  records.append({"gene":g,"broad_class":b,"n_class":int(m.sum()),"top5_n":int((high&m).sum()),"top5_fraction_in_class":p,"fold_enrichment_vs_global":p/baseline,"threshold_log_normalized":th})
fig.suptitle("Marker-enriched cells on shared Fine26 molecular geometry",fontsize=16,fontweight="bold")
fig.savefig(O/"VGAT_VGLUT2_NE_CHAT_TOP5PCT.png",dpi=185,bbox_inches="tight")
fig.savefig(O/"VGAT_VGLUT2_NE_CHAT_TOP5PCT.pdf",bbox_inches="tight")
plt.close(fig)
d=pd.DataFrame(records);d.to_csv(P/"data"/"stage977_marker_top5_broad_enrichment.csv",index=False)
pivot=d.pivot(index="gene",columns="broad_class",values="fold_enrichment_vs_global").reindex(index=[x[0] for x in geneset],columns=["E","I","NE","ChAT","Other"])
fig,ax=plt.subplots(figsize=(8.2,4.4))
im=ax.imshow(np.log2(np.maximum(pivot.to_numpy(float),1e-3)),cmap="coolwarm",vmin=-3,vmax=3,aspect="auto")
ax.set_yticks(range(4),pivot.index);ax.set_xticks(range(5),pivot.columns)
for i in range(4):
 for j in range(5):
  val=pivot.iloc[i,j]
  ax.text(j,i,f"{val:.1f}×",ha="center",va="center",fontsize=9,color="white" if abs(np.log2(max(val,1e-3)))>1.6 else "#181818")
for sp in ax.spines.values():sp.set_visible(False)
fig.colorbar(im,ax=ax,fraction=.034,pad=.03,label="log2 enrichment / all 71,950 cells")
ax.set_title("Top-5%-marker enrichment within frozen broad classes\nDescriptive cell counts, not biological replicates",fontsize=12)
fig.tight_layout()
fig.savefig(O/"TRANSMITTER_TOP5_BROAD_ENRICHMENT.png",dpi=200,bbox_inches="tight")
fig.savefig(O/"TRANSMITTER_TOP5_BROAD_ENRICHMENT.pdf",bbox_inches="tight")
plt.close(fig)
print("STAGE977B_COMPLETE",pivot.round(2).to_string(),flush=True)
