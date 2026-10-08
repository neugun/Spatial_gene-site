"""Stage988c choose a transparent shared-percentile illustrative threshold closest
to user's raw-count target VGAT=15 VGLUT2=5, not optimized against E/I labels.
Scan p60..p85, report all alternatives and section-conditioned overlaps.
"""
from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data";O=P/"assets"/"stage988_gene_percentile"
z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz");g=list(z["genes"].astype(str));sec=z["section"].astype(int)
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0).astype(float)
u=np.maximum(C.sum(1),1);rv=C[:,g.index("Slc32a1")];re=C[:,g.index("Slc17a6")]
nv=10000*rv/u;ne=10000*re/u
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["fine26_stable_id"])
K=pd.read_csv(D/"stage977_fine26_color_key.csv");b=A.fine26_stable_id.map(dict(zip(K.fine26,K.broad))).to_numpy()
rows=[];matrix=[]
for p in np.arange(60,86,.5):
 p=float(p);vthr=float(np.percentile(nv,p));ethr=float(np.percentile(ne,p))
 rawv=float(np.percentile(rv,p));rawe=float(np.percentile(re,p))
 # Objective independent of class labels: absolute fractional log error of corresponding
 # corrected-count empirical quantiles to raw spot equivalent targets 15 and 5.
 relative=float((np.log(rawv/15)**2+np.log(rawe/5)**2)**.5)
 v=nv>vthr;e=ne>ethr
 row=dict(percentile=p,VGAT_norm_v=vthr,VGLUT_norm_v=ethr,
   VGAT_raw_quantile=rawv,VGLUT_raw_quantile=rawe,
   raw_target_relative_log_distance=relative,
   count_vgat_only=int((v&~e).sum()),count_vglut_only=int((~v&e).sum()),
   count_both=int((v&e).sum()),count_neither=int((~v&~e).sum()),
   positive_fraction_vgat=float(v.mean()),positive_fraction_vglut=float(e.mean()))
 rows.append(row)
best=min(rows,key=lambda r:r["raw_target_relative_log_distance"])
p=best["percentile"]
# "best continuous" is not automatically the simplest easily communicated 70th percentile;
# record both and assess p70 in advance for illustrative transparency.
print("BEST_LOG_TARGET_MATCH",json.dumps(best),flush=True)
tab=pd.DataFrame(rows);tab.to_csv(D/"stage988_continuous_percentile_near_spot_targets.csv",index=False)
pview=70.
vthr=np.percentile(nv,pview);ethr=np.percentile(ne,pview)
v=nv>vthr;e=ne>ethr
cls=np.full(len(v),3,int);cls[v&~e]=0;cls[~v&e]=1;cls[v&e]=2
parts=[]
for ss in [0,500,530,560]:
 for bb in ["ALL","E","I","NE","ChAT","Other"]:
  m=(sec==ss if ss else np.ones(len(v),bool))&((b==bb) if bb!="ALL" else True)
  if not np.any(m):continue
  n=int(m.sum());count=np.bincount(cls[m],minlength=4)
  parts.append(dict(section=ss,broad=bb,n=n,vgat_norm_threshold=float(vthr),
   vglut_norm_threshold=float(ethr),VGAT_only=int(count[0]),VGLUT_only=int(count[1]),double=int(count[2]),neither=int(count[3]),
   VGAT_only_pct=100*count[0]/n,VGLUT_only_pct=100*count[1]/n,double_pct=100*count[2]/n,neither_pct=100*count[3]/n))
q=pd.DataFrame(parts);q.to_csv(D/"stage988_shared_p70_classification_by_section.csv",index=False)
xy=np.load(Path(r"G:\Map6_recover_all\stage995_tsne\STAGE995_CURRENT_LOGZ_PCA20_TSNE_30_100.npz"))["embedding"]
fig,axs=plt.subplots(1,2,figsize=(12.6,5.3),layout="constrained")
col=["#4387BE","#E37B34","#A362BC","#B8B8B8"]
for k in (3,2,0,1):
 m=cls==k;axs[0].scatter(xy[m,0],xy[m,1],s=.65,c=col[k],alpha=.18 if k==3 else .72,rasterized=True,linewidths=0)
axs[0].set_title(f"Candidate shared P{pview:g} · VGAT norm>{vthr:.1f}, VGLUT2 norm>{ethr:.1f}",fontsize=11)
axs[0].set_aspect("equal");axs[0].set_xticks([]);axs[0].set_yticks([])
for k,nm in enumerate(["VGAT-only","VGLUT2-only","Both positive","Neither"]):
 axs[1].bar(k,(cls==k).mean()*100,color=col[k],width=.68)
 axs[1].text(k,(cls==k).mean()*100+1,f"{(cls==k).mean()*100:.1f}%",ha="center",fontsize=10)
axs[1].set_ylim(0,68);axs[1].set_xticks(range(4),["VGAT-only","VGLUT2-only","Both","Neither"],rotation=20)
axs[1].set_ylabel("Fraction of all 71,950 cells (%)")
axs[1].set_title("Same-ranked gene-specific thresholds: selection tradeoff")
for ax in axs:
 for side in ax.spines.values():side.set_visible(False)
fig.savefig(O/"STAGE988_P70_SHARED_PERCENTILE_TSNE.png",dpi=200);plt.close(fig)
q0=q[(q.section==0)&q.broad.isin(["ALL","E","I"])]
print("P70_RESULT",q0.to_string(index=False,float_format=lambda x:f"{x:.1f}"),flush=True)
auth={"stage":988,"p70_illustrative":True,"selected_percentile":pview,"closest_full_scan_percentile":best["percentile"],
 "selection_criterion":"Approximate raw count targets 15(VGAT)/5(VGLUT2) by normalized per-gene identical percentile; no E/I labels used for cutoff selection",
 "vgat_normalized_cutoff":float(vthr),"vglut_normalized_cutoff":float(ethr),
 "note":"A percentile cutoff is a marginal-rank visualization tool; independently validate FISH spot detection and neuronal class labels."}
(D/"stage988_p70_selection_authority.json").write_text(json.dumps(auth,indent=2),encoding="utf8")
