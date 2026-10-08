"""Stage988: gene-specific empirical percentiles of DEPTH-normalized and raw spots.
Count metric is 10,000*gene/total27 (no log), with log1p also reported.
Full-cell CDF and positive-only CDF are distinct, do not confuse them.
"""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata
ROOT=Path(r"G:\Spatial_gene_site_publish")
P=ROOT/"perilc-map6-review";D=P/"data";O=P/"assets"/"stage988_gene_percentile";O.mkdir(parents=True,exist_ok=True)
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz");g=list(z["genes"].astype(str));sect=z["section"].astype(int)
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0).astype(float)
assert C.shape==(71950,27)
library=C.sum(axis=1);dep=np.divide(C*10000,library[:,None],out=np.zeros_like(C),where=library[:,None]>0)
gen={"VGAT":"Slc32a1","VGLUT2":"Slc17a6"}
records=[];curves={}
for name,gg in gen.items():
 j=g.index(gg);r=C[:,j];norm=dep[:,j];log=np.log1p(norm)
 for population,m in [("all",np.ones(len(r),bool)),("positive",r>0)]+[(f"S{s}",sect==s) for s in (500,530,560)]:
  a=norm[m];raw=r[m]
  targets=[15.] if name=="VGAT" else [5.]
  for t in targets:
   records.append(dict(gene=name,population=population,kind="target_norm",normalised_value=t,
    log1p_normalised_value=float(np.log1p(t)),percentile_lt=float(100*(a<t).mean()),
    percentile_le=float(100*(a<=t).mean()),percentile_midrank=float(50*((a<t).mean()+(a<=t).mean())),
    raw_count_median_among_above=float(np.median(raw[a>=t])) if np.any(a>=t) else np.nan,
    raw_count_q10_among_above=float(np.quantile(raw[a>=t],.1)) if np.any(a>=t) else np.nan,
    raw_count_q90_among_above=float(np.quantile(raw[a>=t],.9)) if np.any(a>=t) else np.nan,
    raw_median_near_target=float(np.median(raw[np.abs(a-t)<=max(1,.1*t)])) if np.any(np.abs(a-t)<=max(1,.1*t)) else np.nan,
    n_population=len(a)))
  for pct in [50,60,65,70,75,80,85,90,92,95,97,98,99]:
   nt=float(np.percentile(a,pct))
   records.append(dict(gene=name,population=population,kind="common_percentile",normalised_value=nt,
    log1p_normalised_value=float(np.log1p(nt)),percentile_lt=float(100*(a<nt).mean()),
    percentile_le=float(100*(a<=nt).mean()),percentile_midrank=float(pct),
    raw_count_median_among_above=float(np.median(raw[a>=nt])),
    raw_count_q10_among_above=float(np.quantile(raw[a>=nt],.1)),
    raw_count_q90_among_above=float(np.quantile(raw[a>=nt],.9)),
    raw_median_near_target=float(np.median(raw[np.abs(a-nt)<=max(1,.1*nt)])) if np.any(np.abs(a-nt)<=max(1,.1*nt)) else np.nan,
    n_population=len(a)))
 curves[name]=(r,norm)
T=pd.DataFrame(records);T.to_csv(D/"stage988_gene_specific_percentile_alignment.csv",index=False)
print("TARGET_ROWS",T[(T.population=="all")&(T.kind=="target_norm")].to_string(index=False),flush=True)
pairs=T[(T.kind=="common_percentile")&(T.population=="all")].pivot(index="percentile_midrank",columns="gene",values="normalised_value")
print("COMMON_PERCENTILES",pairs.to_string(float_format=lambda f:f"{f:.2f}"),flush=True)
# Curves: x = empirical percentile including zeros; two genes use same X scale.
fig,axs=plt.subplots(1,2,figsize=(12.2,5.4),layout="constrained")
for name,col in [("VGAT","#3a82bd"),("VGLUT2","#dc8631")]:
 raw,norm=curves[name];p=np.linspace(0,100,1001)
 ax=axs[0];y=np.quantile(norm,p/100)
 ax.plot(p,y,color=col,lw=2,label=name)
 target=15 if name=="VGAT" else 5
 pc=float(100*(norm<=target).mean())
 ax.scatter(pc,target,color=col,s=55,zorder=5)
 ax.annotate(f"{name} {target:g}: ≤ {pc:.1f}th %ile",xy=(pc,target),xytext=(-10,18),textcoords="offset points",ha="right",fontsize=9,color=col)
 axs[1].plot(p,np.quantile(raw,p/100),color=col,lw=2,label=name)
 axs[1].scatter(pc,np.quantile(raw,pc/100),s=55,color=col,zorder=5)
for ax in axs:
 ax.set_xlim(0,99.5);ax.set_yscale("symlog",linthresh=5)
 ax.legend(frameon=False);ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False)
 ax.set_xlabel("Empirical percentile, ALL 71,950 cells (zeros included)")
axs[0].set_ylabel("10,000 × gene / sum(all 27 genes)");axs[0].set_title("Gene-specific normalized-count CDF")
axs[1].set_ylabel("Corrected gene spot-count equivalents");axs[1].set_title("Same percentile -> raw/corrected count CDF")
fig.savefig(O/"STAGE988_PERCENTILE_NORMALIZED_AND_RAW_CURVES.png",dpi=205);fig.savefig(O/"STAGE988_PERCENTILE_NORMALIZED_AND_RAW_CURVES.pdf");plt.close(fig)
# Common pct threshold: biology is NOT guaranteed equal accuracy by CDF alignment
def classify(a,b):
 return [float(np.mean(a&~b)),float(np.mean(~a&b)),float(np.mean(a&b)),float(np.mean(~a&~b))]
qrows=[]
v=curves["VGAT"][1];e=curves["VGLUT2"][1]
for pc in [65,70,75,80,85,90,92,95,97,98,99]:
 vt=np.quantile(v,pc/100);et=np.quantile(e,pc/100)
 q=classify(v>vt,e>et)
 qrows.append(dict(percentile=pc,vgat_norm_threshold=float(vt),vglut_norm_threshold=float(et),
   vgat_log_threshold=float(np.log1p(vt)),vglut_log_threshold=float(np.log1p(et)),
   vgat_only=q[0],vglut_only=q[1],co=q[2],neither=q[3],
   vgat_raw_at_percentile=float(np.quantile(curves["VGAT"][0],pc/100)),
   vglut_raw_at_percentile=float(np.quantile(curves["VGLUT2"][0],pc/100))))
Q=pd.DataFrame(qrows);Q.to_csv(D/"stage988_shared_percentile_coexpression_tradeoff.csv",index=False)
fig,axs=plt.subplots(1,2,figsize=(12.5,5),layout="constrained")
axs[0].plot(Q.percentile,Q.vgat_norm_threshold,label="VGAT cutoff",marker="o",color="#3a82bd")
axs[0].plot(Q.percentile,Q.vglut_norm_threshold,label="VGLUT2 cutoff",marker="o",color="#dc8631")
axs[0].set_yscale("symlog",linthresh=10);axs[0].set(ylabel="Gene-specific 10k-normalized cutoff",xlabel="Shared percentile across all cells",title="One percentile -> two gene-specific cutoffs");axs[0].legend(frameon=False)
for k,label,col in [("co","Double-positive","#a45ac2"),("neither","Neither","#999999"),("vgat_only","VGAT-only","#3a82bd"),("vglut_only","VGLUT2-only","#dc8631")]:
 axs[1].plot(Q.percentile,Q[k]*100,label=label,marker="o",markersize=4,color=col)
axs[1].set(xlabel="Shared percentile across all cells",ylabel="Percent of ALL cells",title="Coexpression and unassigned tradeoff");axs[1].legend(frameon=False)
for ax in axs:ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
fig.savefig(O/"STAGE988_COMMON_PERCENTILE_TRADEOFF.png",dpi=200);plt.close(fig)
jsonout={"stage":988,"status":"COMPLETE","cells":len(v),
 "definition":"cell depth corrected to library size: 10000*(corrected gene spot count / sum of all 27 corrected genes). No log transform for requested 15 and 5. log1p reported separately.",
 "all_cell_percentile_includes_zero_spot_cells":True,
 "target":T[(T.population=="all")&(T.kind=="target_norm")].to_dict("records"),
 "shared_percentile_range":[65,99],
 "scientific_limitations":["One equal percentile is a rank/marginal-frequency normalization, not measurement of ambient puncta FPR or neurotransmitter specificity","Gene target 15 and 5 refer to normalized values, not raw spot counts","Raw/corrected counts per threshold vary across cells with total gene depth","S500/S530/S560 are spatial sections of one specimen"]}
(D/"stage988_gene_percentile_authority.json").write_text(json.dumps(jsonout,indent=2),encoding="utf8")
print("STAGE988_COMPLETE",len(T),len(Q),flush=True)
