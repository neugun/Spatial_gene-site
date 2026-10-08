"""Stage982 reference-aligned t-SNE and UMAP controls, frozen cell universe."""
from pathlib import Path
import numpy as np,pandas as pd,json,sys,time
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from openTSNE import affinity,initialization,TSNEEmbedding
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
PUB=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
D=PUB/"data";OUT=PUB/"assets"/"stage982_tsne_comparison";OUT.mkdir(parents=True,exist_ok=True)
LOCAL=Path(r"G:\Map6_recover_all\stage982_tsne");LOCAL.mkdir(parents=True,exist_ok=True)
z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz",allow_pickle=True)
X=np.ascontiguousarray(z["pca15"].astype("float32"))
assert X.shape==(71950,15)
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["fine26_stable_id","fine26"])
K=pd.read_csv(D/"stage977_fine26_color_key.csv")
labels=A.fine26_stable_id.astype(str).to_numpy()
broad=A.fine26_stable_id.map(dict(zip(K.fine26,K.broad))).to_numpy()
colors=dict(zip(K.fine26,K.color))
rng=np.random.default_rng(982)
sel=np.sort(rng.choice(len(X),4000,replace=False))
high=NearestNeighbors(n_neighbors=16).fit(X[sel]).kneighbors(X[sel],return_distance=False)[:,1:]
def metric(E):
 low=NearestNeighbors(n_neighbors=16).fit(E[sel]).kneighbors(E[sel],return_distance=False)[:,1:]
 return {"15nn_Fine26_purity":float((labels[sel][low]==labels[sel,None]).mean()),
 "15nn_Broad_purity":float((broad[sel][low]==broad[sel,None]).mean()),
 "15nn_highD_neighbor_overlap":float(np.mean([len(set(a)&set(b))/15 for a,b in zip(high,low)]))}
colors4=["#4a7db6","#ec892b","#a45ac2","#bbbbbb"]
gene=list(z["genes"].astype(str))
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0)
total=np.maximum(C.sum(1),1)
T=json.loads((D/"stage981_vgat_vglut2_threshold_authority.json").read_text())["selected"]
v=np.log1p(10000*C[:,gene.index("Slc32a1")]/total)
e=np.log1p(10000*C[:,gene.index("Slc17a6")]/total)
a=v>=T["VGAT_log_norm"];b=e>=T["VGLUT2_log_norm"]
q=np.full(len(a),3,int);q[a&~b]=0;q[~a&b]=1;q[a&b]=2
modes={"Stage704_UMAP":np.asarray(z["embedding"])}
if (LOCAL/"STAGE978_ALTERNATE_GEOMETRY.npz").exists(): modes["Stage978_UMAP"]=np.load(LOCAL/"STAGE978_ALTERNATE_GEOMETRY.npz")["embedding"]
cfg=sys.argv[1] if len(sys.argv)>1 else "multiscale_30_300"
perps={"multiscale_30_300":[30,300],"multiscale_30_100":[30,100],"single_50":[50]}[cfg]
path=LOCAL/(cfg+".npz")
if path.exists():E=np.load(path)["embedding"]
else:
 ts=time.time()
 print("TSNE_AFFINITIES",cfg,perps,"n",len(X),flush=True)
 aff=affinity.Multiscale(X,perplexities=perps,metric="euclidean",method="approx",n_jobs=8) if len(perps)>1 else affinity.PerplexityBasedNN(X,perplexity=perps[0],metric="euclidean",method="approx",n_jobs=8)
 init=initialization.pca(X,random_state=982)
 emb=TSNEEmbedding(init,aff,negative_gradient_method="fft",n_jobs=8,random_state=982)
 print("TSNE_OPT_EARLY",cfg,flush=True)
 emb=emb.optimize(n_iter=250,exaggeration=12,momentum=.5,learning_rate=len(X)/12,inplace=True)
 print("TSNE_OPT_LATE",cfg,"elapsed",round(time.time()-ts,1),flush=True)
 emb=emb.optimize(n_iter=450,exaggeration=1,momentum=.8,learning_rate=len(X)/12,inplace=True)
 E=np.asarray(emb,dtype=np.float32)
 np.savez_compressed(path,embedding=E)
 print("TSNE_FINISHED",cfg,"elapsed",round(time.time()-ts,1),flush=True)
modes["Stage982_"+cfg]=E
metrics={name:metric(arr) for name,arr in modes.items()}
print("METRICS",json.dumps(metrics),flush=True)
pd.DataFrame([{"embedding":name,**val} for name,val in metrics.items()]).to_csv(D/f"stage982_{cfg}_embedding_quality.csv",index=False)
fig,axs=plt.subplots(2,len(modes),figsize=(6.4*len(modes),11))
for j,(name,coords) in enumerate(modes.items()):
 ax=axs[0,j]
 for st,c in colors.items():
  m=labels==st;ax.scatter(coords[m,0],coords[m,1],s=.48,color=c,alpha=.75,linewidths=0,rasterized=True)
 ax.set_title(name+" | Fine26\n15NN purity "+f"{metrics[name]['15nn_Fine26_purity']:.3f}")
 ax=axs[1,j]
 for k in (3,2,0,1):
  m=q==k;ax.scatter(coords[m,0],coords[m,1],s=.60,color=colors4[k],alpha=.12 if k==3 else .72,linewidths=0,rasterized=True)
 ax.set_title(name+" | VGAT-only/VGLUT2-only/co/neither")
 for ax in (axs[0,j],axs[1,j]):
  ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect("equal",adjustable="box")
  for sp in ax.spines.values():sp.set_visible(False)
fig.suptitle("Full 71,950-cell Fine26 + transmitter comparison; identical markers/thresholds",fontsize=15)
fig.tight_layout()
fig.savefig(OUT/f"{cfg}_COMPARISON.png",dpi=175,bbox_inches="tight")
fig.savefig(OUT/f"{cfg}_COMPARISON.pdf",bbox_inches="tight")
plt.close(fig)
(D/f"stage982_{cfg}_audit.json").write_text(json.dumps({"stage":982,"config":cfg,"perplexities":perps,"n_cells":len(X),"metrics":metrics,"fine26_source":"frozen Stage631","input":"Stage704 PCA15 from current RouteA 27 genes","initialization":"PCA","lr":len(X)/12,"iterations":[250,450],"caveat":"Visualization-only, Broad cell classes potentially informed by the tested genes"},indent=2))
print("STAGE982_DONE",cfg,flush=True)
