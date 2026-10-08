from pathlib import Path
import numpy as np,pandas as pd,json
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from openTSNE import TSNE
from sklearn.neighbors import NearestNeighbors
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review');D=P/'data';O=P/'assets'/'stage1043_strict_transmitter_gate'
R=Path(r'Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926')
M=np.load(Path(r'G:\Map6_recover_all')/'stage1043_private_filter_masks.npz')
pick=M['selected'];classes=M['class3'][pick]
X=np.concatenate([np.load(Path(r'G:\Map6_recover_all\CURRENT_ROUTEA')/f'S{s}_cell_gene27_current.npy') for s in (500,530,560)]).astype(np.float32)
A=pd.read_csv(R/'stage631_current_3d_atlas'/'CURRENT_3D_CELL_ATLAS.csv.gz',usecols=['section'])
section=A.section.to_numpy(int)[pick]
Y=np.log1p(np.clip(X[pick],0,None))
Y=np.minimum(Y,np.quantile(Y,.997,axis=0)[None,:])
Z=PCA(n_components=20,random_state=1046).fit_transform(StandardScaler().fit_transform(Y))
print('PCA_DONE',Z.shape,flush=True)
E=np.asarray(TSNE(n_components=2,perplexity=100,initialization='pca',
 metric='euclidean',random_state=1046,n_jobs=6,verbose=False).fit(Z))
np.savez_compressed(Path(r'G:\Map6_recover_all')/'stage1046_strict_gate_tsne_private.npz',embedding=E,classes=classes,section=section)
print('TSNE_DONE',E.shape,flush=True)
colors=['#3f7fad','#dd8b52','#9559a6']
fig,axs=plt.subplots(1,2,figsize=(11.5,5.8))
for k in range(3):
 a=classes==k;axs[0].scatter(E[a,0],E[a,1],s=.7,alpha=.53,c=colors[k],linewidths=0,rasterized=True,label=['VGAT-only','VGLUT2-only','Dual'][k])
for i,(sec,c) in enumerate(zip((500,530,560),('#3f7fad','#dd8b52','#9559a6'))):
 a=section==sec;axs[1].scatter(E[a,0],E[a,1],s=.6,c=c,alpha=.5,linewidth=0,rasterized=True,label=f'S{sec}')
for ax in axs:
 ax.legend(frameon=False,markerscale=4);ax.set_xticks([]);ax.set_yticks([])
 ax.set_aspect('equal')
 for sp in ax.spines.values():sp.set_visible(False)
axs[0].set_title('Three marker classes');axs[1].set_title('Section dependence')
fig.suptitle('New label-free t-SNE fitted ONLY to 31,115 selected cells (27 gene panel) - illustrative, no independent subtype validation')
fig.tight_layout();fig.savefig(O/'STAGE1046_NEW_TSNE_USER_SELECTED_CELLS.png',dpi=200);plt.close(fig)
old=np.load(Path(r'G:\Map6_recover_all\stage995_tsne\STAGE995_CURRENT_LOGZ_PCA20_TSNE_30_100.npz'))['embedding'][pick]
rows=[]
for name,coords in [('previous_full_71950_then_filter',old),('refit_only_31115',E)]:
 idx=NearestNeighbors(n_neighbors=16).fit(coords).kneighbors(coords,return_distance=False)[:,1:]
 for feature,labels in [('transmitter_3class',classes),('section',section)]:
  rows.append(dict(embedding=name,feature=feature,knn15_label_agreement=float(np.mean(labels[:,None]==labels[idx]))))
pd.DataFrame(rows).to_csv(D/'stage1046_strict_gate_old_vs_new_embedding_knn.csv',index=False)
(D/'stage1046_strict_gate_tsne_authority.json').write_text(json.dumps({'stage':1046,'method':'Log1p corrected 27-gene counts of 31115 user-selected cells; winsorize per-gene 99.7%; each gene StandardScaler then PCA20; openTSNE Euclidean perplexity100 PCA initialization, random_state1046; no Fine26, broad, or transmitter labels are passed into fit. Compare source Stage995 all-cell embedding filtered post hoc.','caveat':'Mask was selected USING VGAT/VGLUT2 and Snap25; t-SNE class separability may be selection-induced and circular. Same 3 sections one biological specimen. kNN class agreement not cross-validated phenotype prediction.'},indent=2))
print('STAGE1046_COMPLETE',rows,flush=True)
