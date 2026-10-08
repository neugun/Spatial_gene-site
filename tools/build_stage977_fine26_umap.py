"""Map10 Stage977: auditable Fine26 UMAP and transmitter overlays.
No changes to frozen Fine26 labels, gene counts, or spatial XY.
"""
from pathlib import Path
import json, shutil, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import colorcet as cc
from sklearn.neighbors import NearestNeighbors
ROOT=Path(r"G:\Spatial_gene_site_publish")
PUB=ROOT/"perilc-map6-review"
BASE=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
DEST=PUB/"assets"/"stage977_fine26_umap";DEST.mkdir(parents=True,exist_ok=True)
DATA=PUB/"data"
source=BASE/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz"
z=np.load(source,allow_pickle=True)
E=z["embedding"].astype(np.float32)
A=pd.read_csv(BASE/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["fine26","fine26_stable_id","fine26_name","section"])
assert E.shape==(71950,2) and len(A)==71950
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in [500,530,560]],axis=0)
assert C.shape==(71950,27)
genes=[str(g) for g in z["genes"].tolist()]
assert all(g in genes for g in ["Slc32a1","Slc17a6","Slc6a2","Slc5a7"]), genes
palkey=pd.read_csv(DATA/"stage833_fine26_color_key.csv")
types=palkey.fine26.tolist()
assert len(types)==26 and A.fine26_stable_id.nunique()==26 and set(types)==set(A.fine26_stable_id)
# Glasbey was explicitly developed for many high-contrast categorical labels.
# Do not use nearby blue/orange gradients to distinguish the 26 subpopulations.
palette=dict(zip(types, cc.glasbey[:26]))
palkey["color_old"]=palkey["color"]
palkey["color"]=[palette[t] for t in palkey.fine26]
palkey.to_csv(DATA/"stage977_fine26_color_key.csv",index=False)
ids=A.fine26_stable_id.to_numpy(str)
classes=A.fine26.to_numpy(int)
plt.rcParams.update({"font.size":10,"font.family":"DejaVu Sans","savefig.facecolor":"white"})
xl,yl=E[:,0],E[:,1]
xmin,xmax=np.quantile(xl,[.001,.999]);ymin,ymax=np.quantile(yl,[.001,.999])
padx=(xmax-xmin)*.07;pady=(ymax-ymin)*.07
lims=(xmin-padx,xmax+padx,ymin-pady,ymax+pady)
def axis(ax):
    ax.set_xlim(*lims[:2]);ax.set_ylim(*lims[2:])
    ax.set_aspect("equal",adjustable="box")
    ax.set_xticks([]);ax.set_yticks([])
    for sp in ax.spines.values():sp.set_visible(False)
def write(fig, stem):
    fig.savefig(DEST/(stem+".png"),dpi=200,bbox_inches="tight")
    fig.savefig(DEST/(stem+".pdf"),bbox_inches="tight")
    plt.close(fig)
# Main: full Fine26 with separate legend and reduced occlusion.
fig,ax=plt.subplots(figsize=(10.8,9.2))
rng=np.random.default_rng(977)
for k,t in enumerate(types):
    loc=np.flatnonzero(ids==t)
    loc=rng.permutation(loc)
    ax.scatter(xl[loc],yl[loc],c=palette[t],s=1.15,alpha=.80,linewidths=0,rasterized=True)
axis(ax)
handles=[plt.Line2D([0],[0],linestyle="",marker="o",markersize=7,markerfacecolor=palette[t],markeredgewidth=0,label=f"{t}  (n={np.sum(ids==t):,})") for t in types]
fig.legend(handles=handles,loc="center right",bbox_to_anchor=(1.02,.51),ncol=2,fontsize=8,frameon=False,labelspacing=1.0)
fig.suptitle("Fine26 molecular populations · 26 distinguishable colors",fontweight="bold",fontsize=16,x=.40)
fig.text(.06,.035,"Stage704 27-gene embedding · Stage631 Fine26 labels frozen · n=71,950 cells",fontsize=9)
fig.subplots_adjust(right=.72,left=.035,top=.92,bottom=.07)
write(fig,"FINE26_CATEGORICAL_GLASBEY")
# Overlay direct gene readout, not broad class prediction. Same UMAP and global axes.
targets=[("Slc32a1","VGAT / inhibitory"),("Slc17a6","VGLUT2 / excitatory"),("Slc6a2","NE / noradrenergic"),("Slc5a7","ChAT / cholinergic")]
fig,axs=plt.subplots(2,2,figsize=(13,11),constrained_layout=True)
summ=[]
for ax,(g,title) in zip(axs.flat,targets):
    v=C[:,genes.index(g)].astype(float)
    pos=v>0
    vmax=float(np.quantile(np.log1p(v[pos]),.995)) if pos.any() else 1.
    ax.scatter(xl[~pos],yl[~pos],s=.27,color="#d6d6d6",alpha=.20,linewidths=0,rasterized=True)
    sc=ax.scatter(xl[pos],yl[pos],c=np.log1p(v[pos]),s=.8,alpha=.8,cmap="magma",vmin=0,vmax=max(vmax,.001),linewidths=0,rasterized=True)
    axis(ax)
    ax.set_title(f"{title} ({g})\n{int(pos.sum()):,} positive cells / {100*pos.mean():.1f}%",fontsize=12,fontweight="bold",pad=9)
    cb=fig.colorbar(sc,ax=ax,fraction=.035,pad=.02,shrink=.74)
    cb.set_label("log(1 + spot count)",fontsize=8)
    summ.append({"gene":g,"positive":int(pos.sum()),"positive_fraction":float(pos.mean()),"positive_log_q995":vmax})
fig.suptitle("Direct transmitter-marker distribution on the same Fine26 UMAP",fontsize=16,fontweight="bold")
write(fig,"VGAT_VGLUT2_NE_CHAT_MARKER_UMAP")
# Broad label check kept separate from observed single-gene overlays.
bmap=dict(zip(palkey.fine26,palkey.broad))
bcol={"E":"#df6d24","I":"#3b70b7","NE":"#1c9e67","ChAT":"#a64daf","Other":"#969696"}
fig,ax=plt.subplots(figsize=(9,8))
for b in ("Other","I","E","NE","ChAT"):
    ix=np.flatnonzero(np.array([bmap[t]==b for t in ids]))
    ax.scatter(xl[ix],yl[ix],s=.65,alpha=.7,color=bcol[b],linewidths=0,label=f"{b} (n={len(ix):,})",rasterized=True)
axis(ax)
ax.legend(loc="lower center",bbox_to_anchor=(.5,-.04),ncol=5,frameon=False,fontsize=10)
fig.suptitle("Frozen broad-class labels on the molecular embedding",fontsize=15,fontweight="bold")
write(fig,"BROAD_CLASSES_UMAP")
# Quantify marker sensitivity to class labels; avoid equating gene positivity and predicted identities.
tab=[]
for g,title in targets:
    v=C[:,genes.index(g)]>0
    for b in ("E","I","NE","ChAT","Other"):
        m=np.array([bmap[t]==b for t in ids])
        tab.append({"gene":g,"class":b,"n_class":int(m.sum()),"positive_n":int((v&m).sum()),"within_class_positive_fraction":float(v[m].mean())})
pd.DataFrame(tab).to_csv(DATA/"stage977_transmitter_per_broad.csv",index=False)
nn=NearestNeighbors(n_neighbors=16,algorithm="kd_tree").fit(E)
ix=nn.kneighbors(E,return_distance=False)[:,1:]
purity=float(np.mean(classes[ix]==classes[:,None]))
q={"stage":977,"n_cells":71950,"n_fine26":26,"gene_source":"Stage934 Route A counts, section-concatenated row order S500/S530/S560","embedding_source":"Stage704 current UMAP27 unmodified","label_source":"Stage631 frozen Fine26","palette":"colorcet.glasbey, 26 categorical colors","NN15_fine26_purity":purity,"targets":summ,"status":"BUILT_PENDING_VISUAL_QC"}
(DATA/"stage977_fine26_umap_authority.json").write_text(json.dumps(q,indent=2))
print("STAGE977_COMPLETED",json.dumps({"n_cells":71950,"nn15_purity":purity,"targets":summ}),flush=True)
