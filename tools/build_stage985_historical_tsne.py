"""Stage985: recover historical Yuhan-era tSNE on EXACT matched cells and score fairly."""
from pathlib import Path
import numpy as np,pandas as pd,json
from sklearn.neighbors import NearestNeighbors
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
PUB=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
O=PUB/"assets"/"stage985_historical_tsne";O.mkdir(parents=True,exist_ok=True)
old=pd.read_csv(r"G:\PeriLC_current\D_mirror\11_PPTsummary\PeriLC_CHATGPT_RESULTS_20260914\127_EITV2_TSNE\EITv2_tsne_coordinates.csv")
z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz")
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["roi_id","fine26_stable_id","section"])
idx=pd.Index(A.section.astype(str)+"_allroi_"+A.roi_id.astype(int).astype(str)).get_indexer(old.cell_id.astype(str))
assert len(old)==35249 and np.all(idx>=0) and len(set(idx.tolist()))==len(idx)
lab=A.fine26_stable_id.astype(str).to_numpy()[idx]
key=pd.read_csv(PUB/"data"/"stage977_fine26_color_key.csv")
palette=dict(zip(key.fine26,key.color));bmap=dict(zip(key.fine26,key.broad))
broad=np.asarray([bmap[t] for t in lab])
old2=old[["tsne1","tsne2"]].to_numpy(float)
new=np.asarray(z["embedding"])[idx]
full_new_tsne=np.load(Path(r"G:\Map6_recover_all\stage982_tsne")/"multiscale_30_300.npz")["embedding"][idx]
pca=np.asarray(z["pca15"])[idx]
rng=np.random.default_rng(985);sel=np.sort(rng.choice(len(idx),min(4000,len(idx)),replace=False))
ref=NearestNeighbors(n_neighbors=16).fit(pca[sel]).kneighbors(pca[sel],return_distance=False)[:,1:]
def score(E):
 nn=NearestNeighbors(n_neighbors=16).fit(E[sel]).kneighbors(E[sel],return_distance=False)[:,1:]
 return {"Fine26_15NN":float((lab[sel][nn]==lab[sel,None]).mean()),"Broad_15NN":float((broad[sel][nn]==broad[sel,None]).mean()),"highD_15NN_overlap":float(np.mean([len(set(a)&set(b))/15 for a,b in zip(ref,nn)]))}
scores={"Historical_EITv2_tSNE_matched":score(old2),"Current_UMAP_matched":score(new),"Current_multiscale_tSNE_matched":score(full_new_tsne)}
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0)[idx]
g=list(z["genes"].astype(str))
v=C[:,g.index("Slc32a1")]>10;e=C[:,g.index("Slc17a6")]>10
classes=np.full(len(v),3,int);classes[v&~e]=0;classes[~v&e]=1;classes[v&e]=2
col4=["#4a7db6","#ec892b","#a45ac2","#b8b8b8"]
fig,axs=plt.subplots(2,3,figsize=(19.5,12))
for j,(name,X) in enumerate([("Historical EITv2 tSNE",old2),("Current UMAP (matched cells)",new),("Current full-data multiscale tSNE",full_new_tsne)]):
 for st,c in palette.items():
  m=lab==st
  axs[0,j].scatter(X[m,0],X[m,1],c=c,s=.65,alpha=.78,linewidths=0,rasterized=True)
 for k in [3,2,0,1]:
  m=classes==k
  axs[1,j].scatter(X[m,0],X[m,1],c=col4[k],s=.7,alpha=.12 if k==3 else .8,linewidths=0,rasterized=True)
 for i in range(2):
  ax=axs[i,j];ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
  for sp in ax.spines.values():sp.set_visible(False)
 axs[0,j].set_title(f"{name} | frozen Fine26\nmatched n={len(idx):,}")
 axs[1,j].set_title(name+" | paper >10 spots/cell VGAT/VGLUT2")
fig.suptitle("Matched-cell historical tSNE versus current UMAP (no sample-universe confound)",fontsize=14)
fig.tight_layout()
fig.savefig(O/"HISTORICAL_TSNE_VS_CURRENT_EMBEDDINGS_MATCHED.png",dpi=190,bbox_inches="tight")
fig.savefig(O/"HISTORICAL_TSNE_VS_CURRENT_EMBEDDINGS_MATCHED.pdf",bbox_inches="tight")
plt.close(fig)
data={"stage":985,"cells_matched":len(idx),"source":"Stage127 EITv2 tSNE; 35,249 historical cells matched by cell_id to Stage631","scores":scores,"paper_threshold":">10 spots/cell","boundary":"Historical embedding computed with earlier EITv2 labels, not an independent new tSNE; matched current frozen labels used for evaluation"}
(PUB/"data"/"stage985_historical_tsne_quality.json").write_text(json.dumps(data,indent=2))
print("STAGE985_MATCHED_TSNE",json.dumps(scores),flush=True)
