"""Stage1020: first-round and second-round Rong Cell Ranger filtered matrices,
separate provenance and non-overlap cautions. No mixing cell universes.
"""
from pathlib import Path
import scanpy as sc,numpy as np,pandas as pd,json
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
base=Path(r"Z:\sternsonlab\Past Group Members\Rong\RNAseq")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data";O=P/"assets"/"stage1020_rong_rounds";O.mkdir(exist_ok=True)
files=[
 ("round1","Female_Neuron","1st round/Sample_F_Neuron/outs/filtered_feature_bc_matrix.h5"),
 ("round1","Male_Neuron","1st round/Sample_M_Neuron/outs/filtered_feature_bc_matrix.h5"),
 ("round2","Female_Mix","2nd round/Sample_F_Mix/outs/filtered_feature_bc_matrix.h5"),
 ("round2","Male_Mix","2nd round/Sample_M_Mix/outs/filtered_feature_bc_matrix.h5")]
rows=[];meta=[]
for rnd,name,loc in files:
 path=base/loc;print("READ_RONG",rnd,name,flush=True)
 a=sc.read_10x_h5(path)
 # CellRanger feature matrix can have duplicate symbols, use sum of all same symbol columns
 g=np.array(a.var_names.astype(str));X=a.X
 def gene(sym):
  j=np.flatnonzero(g==sym)
  if not len(j):return np.zeros(a.shape[0])
  return np.asarray(X[:,j].sum(axis=1)).ravel().astype(float)
 v=gene("Slc32a1");e=gene("Slc17a6");gad=gene("Gad1")+gene("Gad2")
 glufam=gene("Slc17a6")+gene("Slc17a7")+gene("Slc17a8");snap=gene("Snap25")
 totals=np.asarray(X.sum(axis=1)).ravel()
 for gate,m in [("all",np.ones(len(v),bool)),("Snap25_rawUMI_gt0",snap>0)]:
  n=int(m.sum())
  for t in [0,1,2,3,5]:
   pos=(v>t)&(e>t)
   rows.append(dict(round=rnd,cohort=name,neuron_gate=gate,threshold_gt=t,n_cells=n,
    n_dual=int(np.sum(pos&m)),pct_dual=100*float(np.mean(pos[m])) if n else np.nan,
    pct_vgat=100*float(np.mean((v>t)[m])) if n else np.nan,
    pct_vglut2=100*float(np.mean((e>t)[m])) if n else np.nan,
    pct_both_programs=100*float(np.mean(((gad>0)&(glufam>0))[m])) if n else np.nan,
    median_total_umi=float(np.median(totals[m])) if n else np.nan))
 meta.append(dict(round=rnd,sample=name,filtered_cells=len(v),total_features=a.shape[1],raw_both_gt0=int(np.sum((v>0)&(e>0))),
   neuron_marker_positive=int(np.sum(snap>0))))
 del a,X
 pd.DataFrame(rows).to_csv(D/"stage1020_rong_round1_round2_cellranger_direct_GABA_GLU.csv",index=False)
 print("FINISHED",rnd,name,meta[-1],flush=True)
S=pd.DataFrame(rows)
fig,ax=plt.subplots(figsize=(9,5.4))
c={"Female_Neuron":"#3a7fb4","Male_Neuron":"#a653af","Female_Mix":"#e38747","Male_Mix":"#459d86"}
for name in S.cohort.unique():
 t=S[(S.cohort==name)&(S.neuron_gate=="all")].sort_values("threshold_gt")
 ax.plot(t.threshold_gt,t.pct_dual,"o-",label=name,color=c[name])
ax.set(xlabel="Both Slc32a1 & Slc17a6 exceed raw UMI cutoff",ylabel="Dual-positive in filtered CellRanger cells (%)",
 title="Rong first/second round: raw UMI double-detection by assay and sample (NOT same cell universe)")
ax.legend(frameon=False)
for sp in ("top","right"):ax.spines[sp].set_visible(False)
fig.savefig(O/"STAGE1020_RONG_ROUND1_ROUND2_GABA_GLU_DOUBLE.png",dpi=190);plt.close(fig)
(D/"stage1020_rong_rounds_authority.json").write_text(json.dumps({"stage":1020,"status":"COMPLETE","sources":meta,
"interpretation":"First-round neuron-enriched vs second-round Mix filtered Cell Ranger are different sample definitions. These may contain overlapping/different capture processing. Never add their n or treat four cohorts as four independent animals without source mouse IDs.",
"critical":"12k variants are reprocessings and intentionally excluded from the primary comparison; first-round 4701 curated 25-cluster differs from raw female+male 5161."},indent=2))
print("STAGE1020_ALL_FINISHED",flush=True)
