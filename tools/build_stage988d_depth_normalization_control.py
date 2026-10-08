"""Stage988d: test whether 27-panel normalization creates apparent E/I VGLUT2 specificity.
Compare original H5AD, corrected RouteA, raw percentile and normalized percentile.
"""
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,json
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data";O=P/"assets"/"stage988_gene_percentile"
z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz");G=list(z["genes"].astype(str))
stage=A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["section","fine26_stable_id","roi_id"])
K=pd.read_csv(D/"stage977_fine26_color_key.csv");bro=A.fine26_stable_id.map(dict(zip(K.fine26,K.broad))).to_numpy()
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0)
h5=Path(r"G:\PeriLC_current\D_mirror\cluster\Desktop\EASIFISH\PeriLC\_mapmycells_work\query\periLC_500_530_560_27gene_ALLCELL_FINE26_V2_BROADV4_XYFLIP_APDIRECT.h5ad")
ah=ad.read_h5ad(h5,backed="r");assert np.array_equal(ah.obs_names.to_numpy(),(A.section.astype(str)+"_allroi_"+A.roi_id.astype(str)).to_numpy())
H=np.asarray(ah.X);ah.file.close()
rows=[]
for name,X in [("Current corrected Route A",C),("Historical H5AD",H)]:
 dep=X.sum(1);raw=X[:,G.index("Slc17a6")];N=np.divide(raw*10000,dep,out=np.zeros_like(raw,dtype=float),where=dep>0)
 for kind,tgt in [("raw",raw),("panel_normalized",N)]:
  q70=float(np.quantile(tgt,.70))
  for ss in [0,500,530,560]:
   for broad in ["E","I","NE","ChAT","Other"]:
    m=(bro==broad)&((A.section==ss) if ss else True)
    if not np.any(m):continue
    rows.append(dict(source=name,metric=kind,section=int(ss),broad=broad,n=int(m.sum()),global_p70_cutoff=q70,
       pct_positive=float(100*(tgt[m]>q70).mean()),depth_median=float(np.median(dep[m])),
       raw_gene_median=float(np.median(raw[m])),
       normalized_gene_median=float(np.median(N[m]))))
T=pd.DataFrame(rows);T.to_csv(D/"stage988_h5ad_routeA_depth_normalization_sensitivity.csv",index=False)
current=T[(T.source=="Current corrected Route A")&(T.section==0)&T.broad.isin(["E","I"])]
print("CURRENT_NORM_CONTROL",current.to_string(index=False,float_format=lambda x:f"{x:.2f}"),flush=True)
historic=T[(T.source=="Historical H5AD")&(T.section==0)&T.broad.isin(["E","I"])]
print("HISTORICAL_NORM_CONTROL",historic.to_string(index=False,float_format=lambda x:f"{x:.2f}"),flush=True)
fig,axs=plt.subplots(1,3,figsize=(14.5,5),layout="constrained")
for i,metric in enumerate(["raw","panel_normalized"]):
 ax=axs[i]
 d=T[(T.source=="Current corrected Route A")&(T.metric==metric)&(T.section==0)]
 for j,b in enumerate(["E","I"]):
  v=float(d[d.broad==b].pct_positive.iloc[0]);ax.bar(j,v,color=["#E37B34","#4387BE"][j],width=.65)
  ax.text(j,v+1,f"{v:.1f}%",ha="center")
 ax.set_ylim(0,80);ax.set_xticks([0,1],["E","I"])
 ax.set_title("VGLUT2 global P70: "+("raw counts" if i==0 else "27-gene-depth normalized"))
 ax.set_ylabel("Within-class positive (%)")
di=T[(T.source=="Current corrected Route A")&(T.metric=="panel_normalized")&(T.section==0)&T.broad.isin(["E","I"])]
vals=[float(di[di.broad==b].depth_median.iloc[0]) for b in ["E","I"]]
axs[2].bar([0,1],vals,color=["#E37B34","#4387BE"],width=.65);axs[2].set_xticks([0,1],["E","I"])
axs[2].set_title("Total 27 measured genes/cell · median");axs[2].set_ylabel("Sum of 27 corrected counts")
for j,v in enumerate(vals):axs[2].text(j,v+7,f"{v:.1f}",ha="center")
fig.suptitle("27-gene-panel denominator is highly class dependent; normalization is not an independent sequencing-depth control",fontsize=13)
for ax in axs:
 for sp in ("top","right"):ax.spines[sp].set_visible(False)
fig.savefig(O/"STAGE988_P70_DEPTH_DENOMINATOR_CONTROL.png",dpi=205);plt.close(fig)
jsonout={"stage":988,"control":"Raw vs panel-normalized P70 VGLUT2 in Route A and independent earlier H5AD",
 "E_median_27gene_total":vals[0],"I_median_27gene_total":vals[1],"class_depth_ratio_I_to_E":vals[1]/vals[0],
 "warning":"Total counts from a selected 27-gene biological marker panel are NOT unbiased total transcriptomic sequencing depth; class-dependent denominators may create dramatic normalized E-I specificity.","methodology":"Use this as a negative/sensitivity control before calling same-percentile VGLUT2 a selective transmitter marker."}
(D/"stage988_depth_denominator_control.json").write_text(json.dumps(jsonout,indent=2),encoding="utf8")
print("STAGE988D_DENOMINATOR_QA_COMPLETE",flush=True)
