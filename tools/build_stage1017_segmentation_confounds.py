"""Stage1017: test segmentation area, per-neuron depth and local cell packing
as potential sources of apparent Slc32a1/Slc17a6 co-detection.
Metrics are descriptive; one mouse, Fine26 labels not independent of markers.
"""
from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad
from scipy.stats import spearmanr
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");O=P/"assets"/"stage1017_segmentation_confounds";O.mkdir(exist_ok=True);D=P/"data"
B=Path(r"G:\Map6_recover_all\stage1014_spatial_doublepositive_private")
src=Path(r"G:\PeriLC_current\D_mirror\cluster\Desktop\EASIFISH\PeriLC\_mapmycells_work\query\periLC_500_530_560_27gene_ALLCELL_FINE26_V2_BROADV4_XYFLIP_APDIRECT.h5ad")
a=ad.read_h5ad(src,backed="r")
m=a.obs[["section","roi_id","area","minor_axis_length","major_axis_length","aspect_ratio"]].copy()
a.file.close()
m["section"]=m["section"].astype(int)
all=[]
for ss in (500,530,560):
 q=pd.read_csv(B/f"raw_VGAT10_VGLUT2_3_S{ss}_cells.csv.gz")
 mm=m[m.section==ss]
 out=q.merge(mm,on=["section","roi_id"],validate="one_to_one",how="left")
 print("AREA_MISSING",ss,int(out.area.isna().sum()),len(out),flush=True)
 old=pd.read_csv(Path(r"G:\Map6_recover_all\stage1001_neuron_gate_cell_flags_PRIVATE.csv.gz"),usecols=["section","roi_id","reference_neurons"])
 out=out.merge(old[old.section==ss],on=["section","roi_id"],how="left",validate="one_to_one")
 all.append(out[out.reference_neurons].copy())
T=pd.concat(all,ignore_index=True)
T["double"]=(T.class4==2)
T["snap25_bin"]=pd.qcut(T.snap25,10,duplicates="drop")
T["area_bin"]=np.nan
T["local_dense_bin"]=np.nan
fields=["area","major_axis_length","minor_axis_length","z_um","depth27","snap25","nearest_both_marker_um","mixed_exclusive_neighborhood"]
summary=[]
for ss in (0,500,530,560):
 q=T if ss==0 else T[T.section==ss]
 for field in fields:
  good=q[field].notna()
  a=q.loc[good,field]
  if len(a)<40:continue
  b=q.loc[good,"double"]
  rankq=pd.qcut(a.rank(method="first"),10,labels=False,duplicates="drop")
  for k in range(10):
   m=rankq==k
   summary.append(dict(section=ss,feature=field,decile=k,n=int(m.sum()),feature_median=float(a[m].median()),
       positive=int(b[m].sum()),positive_pct=100*float(b[m].mean())))
S=pd.DataFrame(summary);S.to_csv(D/"stage1017_segment_area_depth_cellpacking_dual_fraction_deciles.csv",index=False)
cols=[("area","Segmentation cell area"),("depth27","Total 27-gene counts"),("nearest_both_marker_um","Distance to both exclusive classes"),("mixed_exclusive_neighborhood","Mixed VGAT/VGLUT neighborhood")]
fig,axs=plt.subplots(2,2,figsize=(12.6,10),layout="constrained")
colors={500:"#397ebb",530:"#df8a43",560:"#8d3baf"}
for ax,(f,title) in zip(axs.flat,cols):
 for ss in (500,530,560):
  q=S[(S.section==ss)&(S.feature==f)].sort_values("decile")
  ax.plot(q.decile+1,q.positive_pct,marker="o",ms=4,lw=1.4,label=f"S{ss}",color=colors[ss])
 ax.set(xlabel="Within-section feature decile (1 low → 10 high)",ylabel="Reference neurons double+ (%)",title=title)
 ax.legend(frameon=False)
 for sp in ["top","right"]:ax.spines[sp].set_visible(False)
fig.suptitle("Are dual-markers segmentation size / expression depth / adjacent-cell composition artifacts?",fontsize=14)
fig.savefig(O/"STAGE1017_DUAL_MARKER_SEGMENTATION_DEPTH_NEIGHBOR_QC.png",dpi=190);plt.close(fig)
f=T.groupby(["Fine26","section"]).agg(n=("double","size"),n_double=("double","sum"),pct_double=("double","mean"),median_area=("area","median"),median_depth27=("depth27","median")).reset_index()
f.pct_double*=100;f.to_csv(D/"stage1017_fine26_section_double_positive_confounds.csv",index=False)
pivot=f.pivot(index="Fine26",columns="section",values="pct_double")
fig,ax=plt.subplots(figsize=(10,9))
im=ax.imshow(pivot.to_numpy(),cmap="PuRd",vmin=0,vmax=75,aspect="auto")
ax.set_xticks(np.arange(3),["S500","S530","S560"])
ax.set_yticks(np.arange(len(pivot.index)),pivot.index)
ax.set(title="VGAT>10 and VGLUT2>3 co-positive (%) in each frozen Fine26 type")
fig.colorbar(im,ax=ax,fraction=.035,pad=.02,label="Double-positive among reference neurons (%)")
fig.tight_layout();fig.savefig(O/"STAGE1017_FINE26_DUAL_POSITIVE_PER_SECTION.png",dpi=180);plt.close(fig)
overall={}
for f in fields:
 q=T[[f,"double"]].dropna()
 rho,p=spearmanr(q[f],q.double.astype(int))
 overall[f]={"spearman_r":float(rho),"n":len(q)}
meta={"stage":1017,"n_reference_neurons":len(T),"double_definition":"VGAT>10 VGLUT2>3 corrected RouteA","reference_cell_shape":"H5AD obs.area from Stage631-matched original cell masks; physical area units source must be verified","summary":overall,
"interpretation":"Descriptive conditional features only; correlation not proof of segmentation contamination; all sections one biological specimen"}
(D/"stage1017_segmentation_confound_authority.json").write_text(json.dumps(meta,indent=2))
print("STAGE1017_DONE",len(T),json.dumps(overall),flush=True)



