from pathlib import Path
import json,numpy as np,pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data"
old=np.load(Path(r"G:\Map6_recover_all\stage982_tsne\multiscale_30_300.npz"))["embedding"]
new=np.load(Path(r"G:\Map6_recover_all\stage995_tsne\STAGE995_CURRENT_LOGZ_PCA20_TSNE_30_100.npz"))["embedding"]
H0=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz")["pca15"]
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0).astype(float)
L=np.log1p(C);L=np.minimum(L,np.quantile(L,.997,axis=0));X=StandardScaler().fit_transform(L).astype(np.float32)
H1=PCA(n_components=20,random_state=995).fit_transform(X).astype(np.float32)
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["fine26_stable_id"])
K=pd.read_csv(D/"stage977_fine26_color_key.csv")
lab=A.fine26_stable_id.astype(str).to_numpy();bm=dict(zip(K.fine26,K.broad));br=np.array([bm[x] for x in lab])
sel=np.sort(np.random.default_rng(995).choice(len(C),4000,False))
m={}
for name,emb in [("OldPCA15_multiscale_30_300",old),("NewLogZPCA20_multiscale_30_100",new)]:
 nn=NearestNeighbors(n_neighbors=16).fit(emb[sel]).kneighbors(emb[sel],return_distance=False)[:,1:]
 m[name]={"Fine26_purity":float((lab[sel][nn]==lab[sel,None]).mean()),"Broad_purity":float((br[sel][nn]==br[sel,None]).mean())}
 for hname,h in [("Old_PCA15",H0),("Current_logZ_PCA20",H1)]:
  true=NearestNeighbors(n_neighbors=16).fit(h[sel]).kneighbors(h[sel],return_distance=False)[:,1:]
  m[name]["knn_overlap_with_"+hname]=float(np.mean([len(set(a)&set(b))/15 for a,b in zip(nn,true)]))
print("STAGE996_METRICS",json.dumps(m,indent=2),flush=True)
(D/"stage996_tsne_preprocess_quality_comparison.json").write_text(json.dumps({"stage":996,"n_subsampled":len(sel),"reference_input_1":"Stage704 PCA15 embedding input","reference_input_2":"Direct log1p corrected counts winsorized each gene p99.7; StandardScaler genes; PCA20","metrics":m,"decision":"Prefer Stage995 if frozen-label neighborhood, lowD-to-highD fidelity, and broad preservation all improve; never use identities in embedding fit","scientific_limitations":["E/I/Fine26 labels may have gene-input circularity","No independent animal biological replicate","Silhouette separated islands alone are not a new molecular cell class"]},indent=2))
