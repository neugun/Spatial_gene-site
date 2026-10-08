"""Map10 Stage978 alternative Fine26 geometry with frozen cell identities."""
from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score
import umap
ROOT=Path(r"G:\Spatial_gene_site_publish");PUB=ROOT/"perilc-map6-review";BASE=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
out=PUB/"assets"/"stage978_umap_geometry";out.mkdir(parents=True,exist_ok=True)
z=np.load(BASE/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz")
P=z["pca15"].astype(np.float32);old=z["embedding"].astype(np.float32)
A=pd.read_csv(BASE/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["fine26","fine26_stable_id"])
assert P.shape==(71950,15) and len(A)==71950
key=pd.read_csv(PUB/"data"/"stage977_fine26_color_key.csv")
names=A.fine26_stable_id.astype(str).to_numpy();lab=A.fine26.to_numpy()
print("STAGE978_UMAP_BEGIN",flush=True)
model=umap.UMAP(n_neighbors=55,min_dist=.09,spread=1.3,metric="cosine",n_components=2,random_state=978,n_epochs=220,low_memory=True,n_jobs=1)
new=model.fit_transform(P).astype(np.float32)
private=Path(r"G:\Map6_recover_all\stage978_geometry");private.mkdir(parents=True,exist_ok=True)
np.savez_compressed(private/"STAGE978_ALTERNATE_GEOMETRY.npz",embedding=new,labels=lab)
r=np.random.default_rng(978);sel=np.sort(r.choice(len(lab),5000,replace=False))
def score(x):
 q=x[sel]
 kn=NearestNeighbors(n_neighbors=16).fit(q).kneighbors(q,return_distance=False)[:,1:]
 return {"15nn_identity_purity":float((lab[sel][kn]==lab[sel,None]).mean()),"silhouette_5k":float(silhouette_score(q,lab[sel],sample_size=3000,random_state=978))}
metrics={"Stage704":score(old),"Stage978":score(new)}
# Show the two geometries using the SAME class colors, without modifying underlying taxonomies.
fig,axes=plt.subplots(1,2,figsize=(17,8.5))
for ax,e,name in zip(axes,[old,new],["Stage704 / current geometry","Stage978 / alternative geometry"]):
 for _,row in key.iterrows():
  s=names==row.fine26
  ax.scatter(e[s,0],e[s,1],s=.8,color=row.color,alpha=.76,linewidths=0,rasterized=True)
 ax.set_aspect("equal",adjustable="box")
 ax.set_xticks([]);ax.set_yticks([])
 for spine in ax.spines.values():spine.set_visible(False)
 q=metrics[name.split(" / ")[0]]
 ax.set_title(name+f"\n15-NN identity purity {q['15nn_identity_purity']:.3f} | silhouette {q['silhouette_5k']:.3f}",fontsize=12)
fig.suptitle("Fine26 visualization geometry comparison · identical cells and frozen labels",fontweight="bold",fontsize=15)
fig.tight_layout(rect=[0,.03,1,.93])
fig.savefig(out/"FINE26_STAGE704_VS_978.png",dpi=190,bbox_inches="tight");fig.savefig(out/"FINE26_STAGE704_VS_978.pdf",bbox_inches="tight");plt.close(fig)
verdict={"stage":978,"status":"BUILT_AWAITING_VISUAL_REVIEW","sample_n":5000,"metrics":metrics,
"source":"Stage704 PCA15 of Stage934 Route A 27-gene counts",
"candidate_params":{"n_neighbors":55,"min_dist":.09,"spread":1.3,"metric":"cosine","n_epochs":220,"random_state":978},
"caution":"UMAP geometry is visualization and cannot serve as independent evidence for cluster identity or molecular function."}
(PUB/"data"/"stage978_geometry_audit.json").write_text(json.dumps(verdict,indent=2))
print("STAGE978_DONE",json.dumps(verdict),flush=True)
