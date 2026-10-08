"""Stage1016: actual ALLEN Pons 143661x32285 raw 10X UMI, selected marker columns,
chunked H5AD backed sparse slice. No use of mapped 27-gene Fine26 to define G/I truth.
"""
from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");O=P/"assets"/"stage1016_allen_pons";O.mkdir(exist_ok=True);D=P/"data"
source=Path(r"Z:\sternsonlab\ABC_Atlas_Downloads\Pons\abc_pipeline_outputs\expression_matrices\WMB-10Xv3\20230630\WMB-10Xv3-P-raw.h5ad")
a=ad.read_h5ad(source,backed="r")
symbols=a.var.gene_symbol.astype(str).to_numpy()
targets=["Snap25","Slc32a1","Gad1","Gad2","Slc17a6","Slc17a7","Slc17a8","Slc6a2","Dbh","Th","Chat","Slc5a7","Slc18a2"]
idx=[int(np.flatnonzero(symbols==g)[0]) for g in targets if np.any(symbols==g)]
targets=[g for g in targets if np.any(symbols==g)]
assert len(targets)>=10
S=np.zeros((len(a),len(idx)),np.float32)
print("ALLEN_PONS_START",a.shape,"genemap",dict(zip(targets,idx)),flush=True)
for lo in range(0,len(a),5000):
 hi=min(len(a),lo+5000)
 z=a[lo:hi,idx].X
 S[lo:hi,:]=np.asarray(z.toarray() if hasattr(z,"toarray") else z,float)
 print("ALLEN_CHUNK",hi,"of",len(a),flush=True) if lo%25000==0 else None
meta=a.obs[["anatomical_division_label","library_label"]].reset_index(drop=True).copy()
a.file.close()
G=dict(zip(targets,S.T))
v=G["Slc32a1"];e=G["Slc17a6"];gad=(G["Gad1"]+G["Gad2"])>0
glu=((G["Slc17a6"]+G["Slc17a7"]+(G["Slc17a8"] if "Slc17a8" in G else 0))>0)
nr=G["Snap25"]>0
rows=[]
for label,m in [("ALL_Pons",np.ones(len(S),bool)),("Snap25_detected",nr)]:
 for cutoff in [0,1,2,3,5,10]:
  j=(v>cutoff)&(e>cutoff)&m
  n=int(m.sum())
  rows.append(dict(population=label,raw_umi_cutoff_gt=cutoff,n=n,
    pct_Slc32a1_detect=float(100*np.mean((v>cutoff)[m])),
    pct_Slc17a6_detect=float(100*np.mean((e>cutoff)[m])),
    n_both=int(j.sum()),pct_both=float(100*j.sum()/max(1,n)),
    pct_both_Gad_or_VGlu_family=float(100*np.mean((gad&glu)[m]))))
T=pd.DataFrame(rows);T.to_csv(D/"stage1016_allen_pons_GABA_GLU_rawUMI_threshold_sweep.csv",index=False)
print("ALLEN_PONS_COUNTS",T.to_string(index=False,float_format=lambda x:f"{x:.2f}"),flush=True)
# within Pons anatomical subdivisions and libraries
rows=[]
for ana,rr in meta.groupby("anatomical_division_label"):
 ix=rr.index.to_numpy();nn=len(ix)
 if nn<40:continue
 for cutoff in [0,1]:
  rows.append(dict(anatomical_division=ana,n=nn,umi_cutoff=cutoff,
   pct_VGAT=float(100*np.mean(v[ix]>cutoff)),
   pct_VGLUT=float(100*np.mean(e[ix]>cutoff)),
   pct_dual=float(100*np.mean((v[ix]>cutoff)&(e[ix]>cutoff))),
   pct_GadGlu=float(100*np.mean((gad&glu)[ix]))))
pd.DataFrame(rows).to_csv(D/"stage1016_allen_pons_anatomical_division_dual_umi.csv",index=False)
fig,axs=plt.subplots(1,2,figsize=(12,4.6),layout="constrained")
for label,c in [("ALL_Pons","#516ab0"),("Snap25_detected","#b45384")]:
 t=T[T.population==label]
 axs[0].plot(t.raw_umi_cutoff_gt,t.pct_both,marker="o",label=label,color=c)
 axs[1].plot(t.raw_umi_cutoff_gt,t.pct_Slc32a1_detect,marker="o",label=label+" VGAT",color=c)
 axs[1].plot(t.raw_umi_cutoff_gt,t.pct_Slc17a6_detect,marker="s",label=label+" VGLUT2",color=c,ls=":")
axs[0].set(xlabel="Raw UMI threshold strictly greater than",ylabel="VGAT and VGLUT2 double-detected (%)",title="Allen Pons 10Xv3 raw transcriptome")
axs[1].set(xlabel="Raw UMI threshold strictly greater than",ylabel="Gene-detected cells (%)",title="Individual marker detection")
for ax in axs:
 ax.legend(frameon=False,fontsize=8);ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
fig.savefig(O/"STAGE1016_ALLEN_PONS_SLC32A1_SLC17A6_RAW_UMI.png",dpi=190);plt.close(fig)
(D/"stage1016_allen_pons_authority.json").write_text(json.dumps({"stage":1016,"status":"COMPLETE","cells":len(S),"genes":targets,"data":"Allen Pons WMB-10Xv3 raw h5ad, original 143661 cells","comparable_scRNA":"Rong 4701 own raw UMI counts, dropout and doublet may differ","not_comparable_directly":"27-gene per-cell corrected FISH spot counts vs 10X single-cell raw UMI are different assay/detection modalities","specimen":"Allen Pons reference release, n independent mice not inferred"},indent=2))
