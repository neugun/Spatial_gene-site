"""Stage987: conditional specificity and ratio: never reinterpret VGAT-gated VGLUT calls as molecular exclusion."""
from pathlib import Path
import numpy as np,pandas as pd,json
from sklearn.metrics import roc_auc_score
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
Z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz");G=list(Z["genes"].astype(str))
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0)
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["section","fine26_stable_id"])
K=pd.read_csv(P/"data"/"stage977_fine26_color_key.csv")
b=A.fine26_stable_id.map(dict(zip(K.fine26,K.broad))).fillna("Other").to_numpy()
v=C[:,G.index("Slc32a1")];e=C[:,G.index("Slc17a6")];section=A.section.to_numpy()
ei=np.isin(b,["E","I"]);E=b=="E";I=b=="I"
ratio=np.log2((e+1)/(v+1))
auc_v=roc_auc_score(E[ei],-v[ei]);auc_e=roc_auc_score(E[ei],e[ei]);auc_ratio=roc_auc_score(E[ei],ratio[ei])
rows=[]
for gcut in [2,3,10]:
 pv=v>10;pe=e>gcut
 Vonly=pe&~pv;Iclass=pv # intentional GABA precedence, not mutually exclusive gene expression
 for sec in [0,500,530,560]:
  m=ei&((section==sec) if sec else True)
  ev=m&Vonly;iv=m&Iclass
  rows.append({"vglut_gt":gcut,"section":sec,"n_EI":int(m.sum()),
   "E_gated_VGLUT_only_n":int(ev.sum()),"E_gated_VGLUT_only_purity_pct":float(100*E[ev].mean()) if ev.sum() else float("nan"),
   "I_gated_VGAT_positive_n":int(iv.sum()),"I_gated_VGAT_positive_purity_pct":float(100*I[iv].mean()) if iv.sum() else float("nan"),
   "assigned_pct":float(100*(Vonly|pv)[m].mean()),"unassigned_pct":float(100*(~Vonly&~pv)[m].mean()),
   "true_dual_positive_pct":float(100*(pv&pe)[m].mean()),
   "auc_E_from_vglut_alone_global":float(auc_e),"auc_E_from_negative_vgat_alone_global":float(auc_v),
   "auc_E_from_log2_vglut1_over_vgat1_global":float(auc_ratio)})
df=pd.DataFrame(rows);df.to_csv(P/"data"/"stage987_conditional_gating_and_marker_ratio.csv",index=False)
fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout="constrained")
for bs,c in (("E","#E37B34"),("I","#4387BE")):
 z=ratio[b==bs];vals,edges=np.histogram(np.clip(z,-7,7),bins=80,density=True)
 axs[0].plot((edges[:-1]+edges[1:])/2,vals,color=c,label=f"Broad {bs} (n={len(z):,})",lw=1.5)
axs[0].set(xlabel="log2[(VGLUT2 + 1)/(VGAT + 1)]",ylabel="Probability density",title="Gene expression ratio; cannot independently confirm Broad labels")
axs[0].legend(frameon=False,fontsize=9)
for cat,c in (("E","#E37B34"),("I","#4387BE")):
 m=b==cat
 axs[1].plot(np.arange(11),[(e[m]>cut).mean()*100 for cut in range(11)],color=c,marker="o",ms=3,label=f"Broad {cat}")
axs[1].set(xlabel="VGLUT2 threshold (> corrected counts)",ylabel="Within-class positive (%)",title=f"VGLUT alone AUC(E vs I) = {auc_e:.3f}")
axs[1].legend(frameon=False)
fig.suptitle(f"Specificity audit · VGAT-alone AUC(E vs I): {auc_v:.3f}; combined log ratio AUC: {auc_ratio:.3f}",fontsize=11)
O=P/"assets"/"stage986_marker_threshold"
fig.savefig(O/"STAGE987_GENE_RATIO_AND_VGLUT_SPECIFICITY.png",dpi=190);plt.close(fig)
print("STAGE987_AUC",round(auc_e,4),round(auc_v,4),round(auc_ratio,4),flush=True)
print("STAGE987_GATED",df[df.section==0].to_string(index=False,float_format=lambda x:f"{x:.2f}"),flush=True)
