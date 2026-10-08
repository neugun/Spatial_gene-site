"""Stage988b: disambiguate raw 15/5 spot-equivalents from 10k-normalized 15/5.
Provide empirical percentiles, conditional depth, and fixed-library-depth conversions.
"""
from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
D=P/"data";O=P/"assets"/"stage988_gene_percentile"
z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz");g=list(z["genes"].astype(str))
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0).astype(float)
depth=C.sum(1)
rows=[];curves={}
for name,gene,rawtarget in [("VGAT","Slc32a1",15),("VGLUT2","Slc17a6",5)]:
 raw=C[:,g.index(gene)]
 norm=np.divide(raw*10000,depth,out=np.zeros_like(raw),where=depth>0)
 for pop,m in (("all",np.ones(len(raw),bool)),("positive",raw>0)):
  rr=raw[m];nn=norm[m];dd=depth[m]
  below=float((rr<rawtarget).mean()*100);at=float((rr<=rawtarget).mean()*100)
  near=(rr>=rawtarget-.5)&(rr<=rawtarget+.5)
  if near.sum()<30:near=(rr>=rawtarget-1)&(rr<=rawtarget+1)
  depth_ref=float(np.median(dd));norm_ref=10000*rawtarget/max(1,depth_ref)
  mapped=float(100*(nn<=norm_ref).mean())
  near_median=float(np.median(nn[near])) if near.any() else np.nan
  near_perc=float(100*(nn<=near_median).mean())
  rows.append(dict(gene=name,population=pop,raw_count_target=rawtarget,
    raw_count_percentile_lt=below,raw_count_percentile_le=at,
    depth_total27_median=depth_ref,
    normalized_equiv_at_population_median_depth=norm_ref,
    normalized_percentile_of_median_depth_equiv=mapped,
    near_raw_target_cell_n=int(near.sum()),
    normalized_median_among_near_raw_target=near_median,
    normalized_q10_among_near_raw_target=float(np.quantile(nn[near],.1)) if near.any() else np.nan,
    normalized_q90_among_near_raw_target=float(np.quantile(nn[near],.9)) if near.any() else np.nan,
    normalized_percentile_of_near_target_median=near_perc,
    raw_median_above_target=float(np.median(rr[rr>rawtarget]))))
 curves[name]=(raw,norm)
T=pd.DataFrame(rows);T.to_csv(D/"stage988_raw15_raw5_normalized_percentile_mapping.csv",index=False)
print("RAW_THRESHOLDS_PERCENTILES",T[T.population=="all"].to_string(index=False),flush=True)
# show threshold locations on cumulative percentile curves
fig,axs=plt.subplots(1,2,figsize=(12,5.5),layout="constrained")
for name,col in (("VGAT","#3675b4"),("VGLUT2","#db8135")):
 raw,norm=curves[name]
 p=np.linspace(0,.9995,1400)
 y=np.quantile(raw,p);x=p*100
 axs[0].plot(x,y,label=name,color=col,lw=2)
 t=15 if name=="VGAT" else 5
 pc=100*(raw<=t).mean()
 axs[0].scatter(pc,t,color=col,s=55,zorder=5)
 axs[0].annotate(f"{name} count>{t}: {pc:.1f}th pct",xy=(pc,t),xytext=(-22,18),textcoords="offset points",color=col,fontsize=9)
 axs[1].plot(x,np.quantile(norm,p),label=name,color=col,lw=2)
 row=T[(T.population=="all")&(T.gene==name)].iloc[0]
 n=row.normalized_equiv_at_population_median_depth
 pc2=row.normalized_percentile_of_median_depth_equiv
 axs[1].scatter(pc2,n,color=col,s=55)
 axs[1].annotate(f"{name} raw>{t} @median depth: {pc2:.1f}th pct",xy=(pc2,n),xytext=(-25,18),textcoords="offset points",color=col,fontsize=8)
for ax in axs:
 ax.set(xlabel="Percentile across ALL cells (zeros retained)")
 ax.set_yscale("symlog",linthresh=6)
 ax.legend(frameon=False);ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False)
axs[0].set(ylabel="Raw corrected spot-count equivalents",title="Requested RAW thresholds: VGAT 15 / VGLUT2 5")
axs[1].set(ylabel="10k-normalized count",title="Equivalent cutoff at MEDIAN library depth")
fig.savefig(O/"STAGE988_RAW15_RAW5_PERCENTILE_MAPPING.png",dpi=210);plt.close(fig)
