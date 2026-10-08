"""Stage995 blind embedding improvement: log-count per-gene robust scaled PCA20, openTSNE multi30/100.
NO Fine26/broad labels used to train 2D projection. Later compare frozen label preservation
and high-dimensional neighbors separately.
"""
from pathlib import Path
import json,time
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from openTSNE import affinity,initialization,TSNEEmbedding
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data";O=P/"assets"/"stage995_tsne_preprocess";O.mkdir(exist_ok=True)
CACHE=Path(r"G:\Map6_recover_all\stage995_tsne");CACHE.mkdir(exist_ok=True)
Z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz")
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0).astype(float)
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["fine26_stable_id"])
K=pd.read_csv(D/"stage977_fine26_color_key.csv");lab=A.fine26_stable_id.astype(str).to_numpy();br=np.array([dict(zip(K.fine26,K.broad))[x] for x in lab])
# log adjusted gene count; robust within-gene quantile scaling prevents extreme hotspots.
L=np.log1p(C);clip=np.quantile(L,.997,axis=0);L=np.minimum(L,clip)
X=StandardScaler().fit_transform(L).astype(np.float32)
pca=PCA(n_components=20,random_state=995).fit_transform(X).astype(np.float32)
print("INPUT_PCA",pca.shape,"explained20",float(PCA(n_components=20,random_state=995).fit(X).explained_variance_ratio_.sum()),flush=True)
rng=np.random.default_rng(995);sel=np.sort(rng.choice(len(C),4000,False))
old=np.load(Path(r"G:\Map6_recover_all\stage982_tsne\multiscale_30_300.npz"))["embedding"]
def metrics(emb):
 m=NearestNeighbors(n_neighbors=16).fit(emb[sel]).kneighbors(emb[sel],return_distance=False)[:,1:]
 return float((lab[sel][m]==lab[sel,None]).mean()),float((br[sel][m]==br[sel,None]).mean())
print("REFERENCE",metrics(old),flush=True)
path=CACHE/"STAGE995_CURRENT_LOGZ_PCA20_TSNE_30_100.npz"
if path.is_file():Y=np.load(path)["embedding"]
else:
 start=time.time()
 aff=affinity.Multiscale(pca,perplexities=[30,100],metric="euclidean",method="approx",n_jobs=8)
 init=initialization.pca(pca,random_state=995)
 emb=TSNEEmbedding(init,aff,negative_gradient_method="fft",n_jobs=8,random_state=995)
 print("EARLY_TSNE",flush=True)
 emb=emb.optimize(n_iter=240,exaggeration=12,momentum=.5,learning_rate=len(C)/12,inplace=True)
 print("LATE_TSNE",round(time.time()-start,1),flush=True)
 emb=emb.optimize(n_iter=400,exaggeration=1,momentum=.8,learning_rate=len(C)/12,inplace=True)
 Y=np.asarray(emb,np.float32);np.savez_compressed(path,embedding=Y)
 print("STAGE995_TSNE_FINISHED",round(time.time()-start,1),flush=True)
metric=metrics(Y)
fig,axs=plt.subplots(1,2,figsize=(14,7))
pal=dict(zip(K.fine26,K.color))
for ax,(coords,title) in zip(axs,[(old,"Stage982 PCA15 30/300"),(Y,"Stage995 logZ PCA20 30/100")]):
 for ty,col in pal.items():
  m=lab==ty;ax.scatter(coords[m,0],coords[m,1],s=.8,c=col,linewidths=0,alpha=.72,rasterized=True)
 ax.set_title(title);ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect("equal")
 for sp in ax.spines.values():sp.set_visible(False)
fig.suptitle("Two UNSUPERVISED t-SNE embeddings · 71,950 cells · same frozen Fine26 labels")
fig.savefig(O/"STAGE995_ORIGINAL_VS_LOGZ_TSNE.png",dpi=190);plt.close(fig)
meta=dict(stage=995,cells=71950,input="log1p(RouteA gene counts) winsorized at each gene 99.7% then gene-wise StandardScaler then PCA20",embedding="tSNE multiscale 30/100 PCA init, lr n/12",original_label_purity=metrics(old),new_label_purity=metric,blind_training=True,caveat="local label concordance descriptive, must also evaluate expression-space neighborhood structure")
(D/"stage995_preprocessing_qa.json").write_text(json.dumps(meta,indent=2))
print("STAGE995_METRICS",meta,flush=True)
