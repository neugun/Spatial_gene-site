import numpy as np,pandas as pd,anndata as ad,json
from pathlib import Path
root=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
p=Path(r"G:\PeriLC_current\D_mirror\cluster\Desktop\EASIFISH\PeriLC\_mapmycells_work\query\periLC_500_530_560_27gene_ALLCELL_FINE26_V2_BROADV4_XYFLIP_APDIRECT.h5ad")
A=ad.read_h5ad(p)
X=np.asarray(A.X,dtype=float)
gene=list(A.var_names)
R=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0)
assert X.shape==R.shape==(71950,27)
B=A.obs["fine26_broad_v4"].astype(str).to_numpy()
if not set(["E","I"]).issubset(set(B)):print("BROAD_VALUES",pd.Series(B).value_counts().to_dict(),flush=True)
out=[]
for g in ["Slc17a6","Slc32a1","Slc6a2","Slc5a7"]:
 j=gene.index(g);h=X[:,j];r=R[:,j]
 for q in ["E","I","NE","ChAT","Other"]:
  mm=B==q
  if not mm.any():continue
  out.append({"gene":g,"class":q,"n":int(mm.sum()),"historic_h5ad_gt10":float((h[mm]>10).mean()),"current_routeA_gt10":float((r[mm]>10).mean()),"historic_q50":float(np.median(h[mm])),"current_q50":float(np.median(r[mm])),"mae":float(np.mean(abs(h[mm]-r[mm]))),"pearson":float(np.corrcoef(h[mm],r[mm])[0,1])})
d=pd.DataFrame(out)
dest=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review\data\stage984_h5ad_vs_routeA_source_audit.csv");d.to_csv(dest,index=False)
print("STAGE984_H5AD_VS_ROUTEA",d[d.gene.isin(["Slc17a6","Slc32a1"])].to_string(index=False),flush=True)
