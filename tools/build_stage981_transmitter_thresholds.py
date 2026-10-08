"""Map10 Stage981: blinded-by-section threshold audit for inhibitory/excitatory markers."""
from pathlib import Path
import numpy as np,pandas as pd,json
from sklearn.metrics import f1_score,balanced_accuracy_score
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
PUB=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
D=PUB/"data";O=PUB/"assets"/"stage981_vgat_vglut2";O.mkdir(parents=True,exist_ok=True)
Z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz",allow_pickle=True)
G=list(Z["genes"].astype(str));S=np.asarray(Z["section"]).astype(int)
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0).astype(float)
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["fine26_stable_id","fine26","section"])
assert C.shape==(71950,27) and np.array_equal(S,A.section.to_numpy(int))
B=A.fine26_stable_id.map(dict(zip(pd.read_csv(D/"stage977_fine26_color_key.csv").fine26,pd.read_csv(D/"stage977_fine26_color_key.csv").broad))).to_numpy()
a=C[:,G.index("Slc32a1")];b=C[:,G.index("Slc17a6")]
total=np.maximum(C.sum(1),1)
nA=np.log1p(10000*a/total);nB=np.log1p(10000*b/total)
isI=B=="I";isE=B=="E";IE=isI|isE
# Historical positive threshold >10 from Wang et al. eLife 2023, retained as fixed-reference baseline
cuts=np.round(np.arange(.5,7.01,.25),2)
train=S==500;val=S==530;test=S==560
rows=[]
for ta in cuts:
 pa=nA>=ta
 for tb in cuts:
  pb=nB>=tb
  ca=pa&~pb;cb=pb&~pa;cc=pa&pb;cn=~pa&~pb
  # score on validation is descriptive because the Broad class may depend on these genes.
  ie=IE&val
  ci=float(np.mean(ca[ie&isI])) if (ie&isI).sum() else 0
  ce=float(np.mean(cb[ie&isE])) if (ie&isE).sum() else 0
  crossi=float(np.mean(cb[ie&isI])) if (ie&isI).sum() else 0
  crosse=float(np.mean(ca[ie&isE])) if (ie&isE).sum() else 0
  cop=float(np.mean(cc[ie]));nei=float(np.mean(cn[ie]))
  # Prespecified multiobjective: high target sensitivity, low opposite-class contamination;
  # mild penalties for co-detection and missing cells, no requirement to force co to zero.
  objective=.5*(ci+ce)-.6*.5*(crossi+crosse)-.22*cop-.12*nei
  rows.append((ta,tb,objective,ci,ce,crossi,crosse,cop,nei))
Q=pd.DataFrame(rows,columns=["VGAT_threshold","VGLUT2_threshold","validation_objective","I_only_recall","E_only_recall","I_to_E_wrong","E_to_I_wrong","IE_cofraction","IE_neither"])
Q.to_csv(D/"stage981_normalized_threshold_sweep.csv",index=False)
best=Q.sort_values(["validation_objective","IE_cofraction"],ascending=[False,True]).iloc[0]
ta,tb=float(best.VGAT_threshold),float(best.VGLUT2_threshold)
ref={"paper_raw_gt10":(a>10,b>10),"normalized_selected":(nA>=ta,nB>=tb),
     "raw_gt0":(a>0,b>0)}
summ=[]
colors=["#4a7db6","#ec892b","#a45ac2","#bbbbbb"]
names=["VGAT-only","VGLUT2-only","Co-positive","Neither"]
for label,(pa,pb) in ref.items():
 for sec in (500,530,560):
  m=S==sec
  for broad in ("E","I","NE","ChAT","Other","ALL"):
   mm=m&(B==broad) if broad!="ALL" else m
   pct=[float(np.mean(v[mm])) for v in (pa&~pb,pb&~pa,pa&pb,~pa&~pb)]
   summ.append(dict(method=label,section=sec,broad=broad,n=int(mm.sum()),VGAT_only=pct[0],VGLUT2_only=pct[1],co=pct[2],neither=pct[3]))
T=pd.DataFrame(summ);T.to_csv(D/"stage981_vgat_vglut2_four_class_by_section.csv",index=False)
fig,ax=plt.subplots(figsize=(7.5,6.2))
z=Q.pivot(index="VGLUT2_threshold",columns="VGAT_threshold",values="validation_objective").sort_index()
im=ax.imshow(z.values,origin="lower",aspect="auto",extent=[cuts.min(),cuts.max(),cuts.min(),cuts.max()],cmap="viridis")
ax.plot(ta,tb,marker="o",ms=8,color="white",mec="black",mew=1)
ax.set_xlabel("VGAT normalized log threshold");ax.set_ylabel("VGLUT2 normalized log threshold")
ax.set_title("Section-held-out threshold scan · validation S530",fontsize=12)
fig.colorbar(im,ax=ax,label="multiobjective validation score");fig.tight_layout();fig.savefig(O/"VGAT_VGLUT2_THRESHOLD_SCAN.png",dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(16,5.4))
sel=ref["normalized_selected"]
cls=np.full(len(S),3,int);cls[sel[0]&~sel[1]]=0;cls[sel[1]&~sel[0]]=1;cls[sel[0]&sel[1]]=2
for ax,sec in zip(axes,(500,530,560)):
 ix=S==sec
 count=np.bincount(cls[ix],minlength=4)
 ax.barh(np.arange(4),count/count.sum(),color=colors,height=.7)
 ax.set_yticks(np.arange(4),names);ax.invert_yaxis();ax.set_xlim(0,1)
 ax.set_title(f"S{sec} n={ix.sum():,}")
 for j,val in enumerate(count/count.sum()):ax.text(val+.01,j,f"{val:.1%}",va="center",fontsize=9)
 for sp in ("right","top"):ax.spines[sp].set_visible(False)
fig.suptitle(f"VGAT≥{ta:g} / VGLUT2≥{tb:g}: full-cell four-class composition",fontsize=14)
fig.tight_layout();fig.savefig(O/"VGAT_VGLUT2_FOUR_CLASS_SECTIONS.png",dpi=185);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(13,6))
for ax,(label,(pa,pb)) in zip(axes,[("Paper raw >10",ref["paper_raw_gt10"]),("Normalized threshold",ref["normalized_selected"])]):
 a0=np.full(len(S),3,int);a0[pa&~pb]=0;a0[pb&~pa]=1;a0[pa&pb]=2
 for k in (3,2,0,1):
  idx=a0==k
  ax.scatter(Z["embedding"][idx,0],Z["embedding"][idx,1],s=.5,color=colors[k],alpha=.55 if k!=3 else .13,linewidths=0,rasterized=True)
 ax.set_title(label);ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
 for sp in ax.spines.values():sp.set_visible(False)
fig.suptitle("Same molecular embedding; alternative transparent four-class marker rules")
fig.tight_layout();fig.savefig(O/"VGAT_VGLUT2_FOUR_CLASS_ON_UMAP.png",dpi=195);plt.close(fig)
metadata=dict(stage=981,status="COMPLETE",selected={"VGAT_log_norm":ta,"VGLUT2_log_norm":tb,"validation_S530_score":float(best.validation_objective)},
train_section=500,validation_section=530,untouched_test_section=560,
primary_check="Compare fixed paper raw-count >10 to normalized two-gene thresholds, all sections and broad types",
scientific_limitations=["Broad E/I labels may use these genes and are not independent ground truth","High positivity thresholds must not be tuned merely to eliminate co-positive cells","Three sections are one specimen and not biological replicates"],
genes=G,reference="Wang et al. eLife 2023 Supp4 used >10 spots/cell for co-expression")
(D/"stage981_vgat_vglut2_threshold_authority.json").write_text(json.dumps(metadata,indent=2))
print("STAGE981_DONE",json.dumps(metadata["selected"]),flush=True)
print("S560_HELDOUT",T[(T.method=="normalized_selected")&(T.section==560)&T.broad.isin(["I","E","ALL"])].to_string(index=False),flush=True)
