"""Stage1015: original user's own 4,701-neuron raw-count scRNA reference.
Count matrix already frozen, 25 clusters, NO re-clustering. Multiple neurotransmitter
gene checks and UMI-depth stratification, not paper-specific arbitrary threshold.
"""
from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data";O=P/"assets"/"stage1015_scrna_gaba_glu";O.mkdir(exist_ok=True)
source=Path(r"G:\PeriLC_current\D_mirror\11_PPTsummary\Slides_periLCFISH\periLCFISH\scRNA_spatial_mapping\scRNA_REFERENCE_PIPELINE_FROZEN_20260917\results\scrna_reference\scrna_reference_25cluster.h5ad")
a=ad.read_h5ad(source)
genes=["Snap25","Slc32a1","Gad1","Gad2","Slc17a6","Slc17a7","Slc17a8","Slc6a2","Dbh","Th","Slc18a2","Chat","Slc5a7","Slc18a3"]
exists=[g for g in genes if g in a.var_names]
X=np.asarray(a[:,exists].layers["counts"].toarray() if hasattr(a[:,exists].layers["counts"],"toarray") else a[:,exists].layers["counts"],float)
cs=dict(zip(exists,X.T));cluster=a.obs.cluster.astype(str).to_numpy()
umi=a.obs.nCount_RNA.to_numpy(float);sample=a.obs["orig.ident"].astype(str).to_numpy()
v=cs["Slc32a1"];e=cs["Slc17a6"];gat=(cs["Gad1"]+cs["Gad2"])>0
v2=(cs["Slc17a6"]+cs["Slc17a7"]+cs["Slc17a8"])>0
rows=[]
for ss in ["ALL"]+sorted(pd.unique(sample)):
 for c in ["ALL"]+sorted(pd.unique(cluster)):
  for gt in [0,1,2,3,5,10]:
   m=((sample==ss) if ss!="ALL" else np.ones(len(X),bool))&((cluster==c) if c!="ALL" else np.ones(len(X),bool))
   if not m.any():continue
   vg=v>gt;gl=e>gt
   n=int(m.sum());q=int(np.sum(m&vg&gl))
   rows.append(dict(source="Rong_frozen_25cluster_scRNA",specimen=ss,cluster=c,threshold="both genes >"+str(gt)+" raw UMI",
     n=n,vgat_pos=int(np.sum(m&vg)),vglut2_pos=int(np.sum(m&gl)),double=q,
     pct_double=100*q/n,pct_double_given_any=100*q/max(1,int(np.sum(m&(vg|gl)))),
     pct_gad_any=100*np.mean(gat[m]),pct_vglut_any=100*np.mean(v2[m]),
     pct_both_full_programs=100*np.mean((gat&v2)[m]),
     median_total_umi=float(np.median(umi[m]))))
T=pd.DataFrame(rows);T.to_csv(D/"stage1015_Rong_4701_cluster25_GABA_GLU_rawUMI_sensitivity.csv",index=False)
print("RONg_REFERENCE_4701",T.query('specimen=="ALL" and cluster=="ALL"').to_string(index=False),flush=True)
# per-cluster summary
summary=T[(T.specimen=="ALL")&(T.threshold=="both genes >0 raw UMI")&(T.cluster!="ALL")].copy()
summary=summary.sort_values(["pct_double","n"],ascending=[False,False])
fig,axs=plt.subplots(1,2,figsize=(14.5,7.2))
ax=axs[0];yy=np.arange(len(summary))
ax.barh(yy,summary.pct_double,color="#834aa5");ax.set_yticks(yy,summary.cluster,fontsize=8);ax.invert_yaxis()
ax.set(xlabel="Slc32a1 and Slc17a6 both >0 raw UMI (%)",title="Rong own scRNA: 25 frozen clusters, true raw UMI")
for i,r in enumerate(summary.itertuples(index=False)):
 ax.text(r.pct_double+.4,i,f"{r.n}",va="center",fontsize=7)
ax=axs[1]
x0=np.log1p(v);y0=np.log1p(e)
ax.scatter(x0,y0,s=4,alpha=.18,color="#ad9fb6",linewidths=0)
ax.set(xlabel="log1p Slc32a1 raw UMI",ylabel="log1p Slc17a6 raw UMI",title="Raw UMI co-detection, all 4,701 own neurons")
for ax in axs:
 for sp in ("top","right"):ax.spines[sp].set_visible(False)
fig.suptitle("Own RNA is the internal independent full-transcriptome reference; do not merge with 27-gene FISH thresholds",fontsize=12)
fig.subplots_adjust(top=.90,left=.10,right=.98,bottom=.11,wspace=.29)
fig.savefig(O/"STAGE1015_RONG_25CLUSTER_GABA_GLU_DOUBLE_RAW_UMI.png",dpi=190);plt.close(fig)
# Multi-marker status and depth deciles, empirical doublet-like dependency 
bins=pd.qcut(umi,10,labels=False,duplicates="drop")
tab=[]
for k in sorted(np.unique(bins)):
 m=bins==k
 tab.append(dict(umi_decile=int(k),n=int(m.sum()),umi_median=float(np.median(umi[m])),pct_Slc32a1_Slc17a6_double=100*np.mean((v[m]>0)&(e[m]>0)),
   pct_Gad1or2_VGlut_any=100*np.mean(gat[m]&v2[m])))
pd.DataFrame(tab).to_csv(D/"stage1015_Rong_scRNA_double_positive_vs_UMI_depth.csv",index=False)
fig,ax=plt.subplots(figsize=(7,5))
u=pd.DataFrame(tab)
ax.plot(u.umi_median,u.pct_Slc32a1_Slc17a6_double,"o-",color="#8d3baf",label="Slc32a1 & Slc17a6 raw>0")
ax.plot(u.umi_median,u.pct_Gad1or2_VGlut_any,"s-",color="#d97937",label="Gad1/2 & Slc17a6/7/8 any raw>0")
ax.set(xlabel="Median total UMI / cell",ylabel="Dual-detected scRNA cells (%)",title="Raw UMI depth dependence")
ax.legend(frameon=False)
for sp in ("top","right"):ax.spines[sp].set_visible(False)
fig.savefig(O/"STAGE1015_RONG_DOUBLE_UMI_DEPTH_QC.png",dpi=205);plt.close(fig)
meta={"stage":1015,"status":"COMPLETE","cells":len(X),"clusters":len(np.unique(cluster)),"specimens":sorted(pd.unique(sample).tolist()),"genes":exists,
"unit":"unmodified counts layer RAW UMI; 0/1/2/3/5/10 thresholds all descriptive; no molecule-per-cell threshold transfer to FISH",
"interpretation":"Dropout differs from smFISH; Slc32a1+Slc17a6 co-UMI detection not direct dual neurotransmitter release; doublet and ambient must be assessed prior to biological claim"}
(D/"stage1015_rong_reference_authority.json").write_text(json.dumps(meta,indent=2),encoding="utf8")
print("STAGE1015_RONG_DONE",len(X),len(np.unique(cluster)),flush=True)
