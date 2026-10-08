"""Stage1019: adjust apparent acquisition tile-overlap association for
cell size, detected 27-gene depth, Fine26 class, section.
Within one animal, purely descriptive conditional association.
"""
from pathlib import Path
import json,ast
import numpy as np,pandas as pd,anndata as ad
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import LogisticRegression
from scipy.special import expit
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data";O=P/"assets"/"stage1019_edge_adjusted";O.mkdir(exist_ok=True)
private=Path(r"G:\Map6_recover_all\stage1014_spatial_doublepositive_private")
geo=pd.read_csv(D/"stage1018_tile_geometry_and_Stage631_affines.csv")
a=ad.read_h5ad(Path(r"G:\PeriLC_current\D_mirror\cluster\Desktop\EASIFISH\PeriLC\_mapmycells_work\query\periLC_500_530_560_27gene_ALLCELL_FINE26_V2_BROADV4_XYFLIP_APDIRECT.h5ad"),backed="r")
shape=a.obs[["section","roi_id","area","major_axis_length","minor_axis_length"]].reset_index(drop=True)
shape.section=shape.section.astype(int);a.file.close()
ff=pd.read_csv(Path(r"G:\Map6_recover_all\stage1001_neuron_gate_cell_flags_PRIVATE.csv.gz"),usecols=["section","roi_id","reference_neurons"])
blocks=[]
for ss in (500,530,560):
 z=pd.read_csv(private/f"raw_VGAT10_VGLUT2_3_S{ss}_cells.csv.gz")
 z=z.merge(shape[shape.section==ss],on=["section","roi_id"],how="left",validate="one_to_one")
 z=z.merge(ff[ff.section==ss],on=["section","roi_id"],validate="one_to_one")
 assert z.area.notna().all()
 q=geo[geo.section==ss].iloc[0]
 xmid=np.array(ast.literal_eval(q.overlap_centers_x_um))
 ymid=np.array(ast.literal_eval(q.overlap_centers_y_um))
 xhalf=np.array(ast.literal_eval(q.overlap_widths_x_um))/2
 yhalf=np.array(ast.literal_eval(q.overlap_widths_y_um))/2
 dx=np.abs(z.x_um.to_numpy()[:,None]-xmid[None,:])
 dy=np.abs(z.y_um.to_numpy()[:,None]-ymid[None,:])
 near=((dx<=xhalf[None,:]).any(1)|(dy<=yhalf[None,:]).any(1))
 z["native_overlap"]=near.astype(int)
 z["near_15um"]=((dx<=15).any(1)|(dy<=15).any(1)).astype(int)
 z=z[z.reference_neurons].copy()
 z["dual"]=(z.class4==2).astype(int)
 z["logarea"]=np.log1p(z.area.clip(lower=0))
 z["logdepth"]=np.log1p(z.depth27.clip(lower=0))
 blocks.append(z)
T=pd.concat(blocks,ignore_index=True)
assert len(T)==68572
cols=["logarea","logdepth","native_overlap","near_15um"]
print("STAGE1019_COVARIATE_STDS",T[cols].std().to_dict(),flush=True)
# Cohort-matched g-computation: predict each cell twice, near=0/1.
from sklearn.preprocessing import OneHotEncoder
def robust_model(edge,use_fine):
 C=pd.DataFrame({"logarea":T.logarea,"logdepth":T.logdepth,"seam":T[edge]})
 C["S530"]=(T.section==530).astype(float);C["S560"]=(T.section==560).astype(float)
 if use_fine:
  dummy=pd.get_dummies(T.Fine26,prefix="fine",drop_first=True,dtype=float)
  C=pd.concat([C.reset_index(drop=True),dummy.reset_index(drop=True)],axis=1)
 # scaling continuous only; leave boundary indicator for direct OR in coefficient
 scale=StandardScaler();C[["logarea","logdepth"]]=scale.fit_transform(C[["logarea","logdepth"]])
 model=LogisticRegression(penalty="l2",C=1.0,solver="lbfgs",max_iter=500)
 model.fit(C,T.dual)
 x0=C.copy();x1=C.copy();x0["seam"]=0.;x1["seam"]=1.
 p0=model.predict_proba(x0)[:,1].mean();p1=model.predict_proba(x1)[:,1].mean()
 coef=float(model.coef_[0][C.columns.get_loc("seam")])
 return dict(n=len(T),seam_variable=edge,adjust_Fine26=use_fine,OR_adjusted=float(np.exp(coef)),
     standardized_risk_inside_pct=100*float(p1),standardized_risk_outside_pct=100*float(p0),
     standardized_diff_pct=100*float(p1-p0),crude_inside_pct=100*T.dual[T[edge]==1].mean(),
     crude_outside_pct=100*T.dual[T[edge]==0].mean(),logarea_coef=float(model.coef_[0][C.columns.get_loc("logarea")]),
     logdepth_coef=float(model.coef_[0][C.columns.get_loc("logdepth")]))
rows=[]
for edge in ["native_overlap","near_15um"]:
 for fine in [False,True]:
  rows.append(robust_model(edge,fine));print("MODEL_RESULT",rows[-1],flush=True)
Z=pd.DataFrame(rows);Z.to_csv(D/"stage1019_tile_seam_size_depth_adjusted_association.csv",index=False)
# Nonparametric area-depth stratified matched prevalence (20 quantiles global)
for t in ["area","depth27"]:
 T[t+"_decile"]=pd.qcut(T[t].rank(method="first"),10,labels=False)
for typ in ["native_overlap","near_15um"]:
 strat=[]
 for (ss,da,dc),gr in T.groupby(["section","area_decile","depth27_decile"],observed=True):
  if gr[typ].nunique()<2 or len(gr)<20:continue
  near=gr[gr[typ]==1];far=gr[gr[typ]==0]
  if len(near)<3 or len(far)<3:continue
  strat.append(dict(section=ss,area_decile=da,depth_decile=dc,
   n_near=len(near),n_far=len(far),near_rate=float(near.dual.mean()),far_rate=float(far.dual.mean())))
 S=pd.DataFrame(strat);S["test"]=typ
 S.to_csv(D/f"stage1019_{typ}_matched_area_depth_strata.csv",index=False)
 print("MATCHED",typ,"strata",len(S),"weighted_diff",float(np.average(S.near_rate-S.far_rate,weights=S.n_near+S.n_far)),flush=True)
# Presentation: crude vs adjusted risk difference, section-wise near/far
fig,axs=plt.subplots(1,2,figsize=(11.8,5),layout="constrained")
q=Z
for j,edge in enumerate(["native_overlap","near_15um"]):
 ax=axs[j];s=q[q.seam_variable==edge].reset_index(drop=True)
 crude=float(s.crude_inside_pct.iloc[0]-s.crude_outside_pct.iloc[0])
 for k,(lab,val,col) in enumerate([("Crude",crude,"#ae62bd"),("Area/depth/section",float(s[s.adjust_Fine26==False].standardized_diff_pct.iloc[0]),"#dd8947"),
  ("+ frozen Fine26",float(s[s.adjust_Fine26==True].standardized_diff_pct.iloc[0]),"#397ebb")]):
  ax.bar(k,val,width=.66,color=col)
  ax.text(k,val+(0.25 if val>=0 else -0.25),f"{val:+.2f} pp",ha="center",va="bottom" if val>=0 else "top")
 ax.axhline(0,color="#777777",lw=.8)
 ax.set_xticks([0,1,2],["Raw","Area/depth/section","+ Fine26"],rotation=19,ha="right")
 ax.set(title=edge.replace("_"," ").capitalize(),ylabel="Tile-edge minus away double+ prevalence (percentage points)")
 for sp in ("top","right"):ax.spines[sp].set_visible(False)
fig.suptitle("Acquisition tile overlap: apparent double positivity after correcting for cell size, detection depth and composition",fontsize=12)
fig.savefig(O/"STAGE1019_TILE_SEAM_ADJUSTED_FOR_CELL_SIZE_AND_DEPTH.png",dpi=190);plt.close(fig)
(D/"stage1019_size_depth_adjustment_authority.json").write_text(json.dumps({"stage":1019,"status":"COMPLETE","cells":len(T),
 "endpoint":"VGAT>10 VGLUT2>3 double positive, old reference neurons",
 "exposure":"Actual tiles.json overlap stripes transformed to Stage631, not final pixel blend seam",
 "adjustments":["log cell mask area","log total corrected27 counts","3 sections","optional frozen Fine26"],"estimators":"penalized logistic standardization (g-computation) and 10x10 area/depth strata sensitivity",
 "caveats":["One biological specimen: associations descriptive; no inferential animal-level p values","27 gene depth includes endpoint VGAT and VGLUT2 and may induce overadjustment (report both adjusted and crude)","Fine26 labels partially depend on endpoint markers; source of circularity; do not interpret as causal mediation","Measured cell area correlation can be expected even with true separate transcripts; cannot infer segmentation artifact without raw puncta count per mask"]},indent=2))
print("STAGE1019_DONE",flush=True)
