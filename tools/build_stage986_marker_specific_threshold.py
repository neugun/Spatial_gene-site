"""Map10 Stage986: per-gene raw-count thresholds, independent source and Fine26/tSNE QC.
Conservatively keep explicit > versus >=; no cluster forcing.
"""
from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
O=P/"assets"/"stage986_marker_threshold";O.mkdir(parents=True,exist_ok=True)
D=P/"data"
Z=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz")
G=list(Z["genes"].astype(str))
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)],axis=0)
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["section","roi_id","fine26_stable_id","x_um","y_um"])
K=pd.read_csv(D/"stage977_fine26_color_key.csv")
bmap=dict(zip(K.fine26,K.broad))
broad=A.fine26_stable_id.map(bmap).fillna("Other").to_numpy()
S=A.section.to_numpy(int)
assert len(A)==C.shape[0]==len(Z["embedding"])==71950 and np.array_equal(S,Z["section"].astype(int))
v=C[:,G.index("Slc32a1")]
e=C[:,G.index("Slc17a6")]
E=broad=="E";I=broad=="I";N=E|I
cuts=[0,1,2,3,4,5,7,10,15,20]
vgat_cuts=[8,10,12]
rows=[]
for vc in vgat_cuts:
 for ec in cuts:
  pv=v>vc;pe=e>ec
  cls=np.full(len(v),3,np.uint8)
  cls[pv&~pe]=0;cls[~pv&pe]=1;cls[pv&pe]=2
  for sec in [0,500,530,560]:
   msec=(S==sec) if sec else np.ones(len(v),bool)
   for cat,cm in (("ALL",np.ones(len(v),bool)),("neuronal_EI",N),("E",E),("I",I),("NE",broad=="NE"),("ChAT",broad=="ChAT"),("Other",broad=="Other")):
    m=msec&cm;n=int(m.sum())
    if not n:continue
    q=np.bincount(cls[m],minlength=4)
    rows.append(dict(vgat_thr=vc,vglut_thr=ec,section=sec,broad=cat,n=n,
      vgat_detect=float((pv[m]).mean()),vglut_detect=float((pe[m]).mean()),
      vgat_only=int(q[0]),vglut_only=int(q[1]),co_positive=int(q[2]),neither=int(q[3]),
      vgat_only_pct=100*q[0]/n,vglut_only_pct=100*q[1]/n,
      co_positive_pct=100*q[2]/n,neither_pct=100*q[3]/n,
      co_of_any_pct=100*q[2]/max(1,(q[0]+q[1]+q[2])),
      vglut_purity_E_pct=100*np.sum(pe[m]&E[m])/max(1,np.sum(pe[m]&N[m])) if cat=="ALL" else np.nan))
T=pd.DataFrame(rows)
T.to_csv(D/"stage986_vgat_vglut_threshold_all_sections_broad.csv",index=False)
# Source old raw expression counts from stored H5AD, row order matched by roi identity.
import anndata as ad
h5=Path(r"G:\PeriLC_current\D_mirror\cluster\Desktop\EASIFISH\PeriLC\_mapmycells_work\query\periLC_500_530_560_27gene_ALLCELL_FINE26_V2_BROADV4_XYFLIP_APDIRECT.h5ad")
ad0=ad.read_h5ad(h5,backed="r")
assert list(ad0.var_names.astype(str))==G and len(ad0)==len(C)
his=np.asarray(ad0[:,["Slc32a1","Slc17a6"]].X,dtype=float)
# Current row order validated by cell ids, required.
expect=(A.section.astype(str)+"_allroi_"+A.roi_id.astype(int).astype(str)).to_numpy()
assert np.array_equal(expect,ad0.obs_names.to_numpy()),"H5AD and RouteA order mismatch"
ad0.file.close()
oldrows=[]
for name,vv,ee in (("CURRENT_ROUTEA",v,e),("HISTORICAL_H5AD",his[:,0],his[:,1])):
 for vg_cut in (10,):
  for gl_cut in (0,1,2,3,4,5,7,10):
   pv=vv>vg_cut; pe=ee>gl_cut
   for sec in [500,530,560]:
    for cat,mcat in (("E",E),("I",I),("ALL",np.ones(len(v),bool))):
     m=(S==sec)&mcat
     oldrows.append(dict(data=name,section=sec,broad=cat,vgat_gt=vg_cut,vglut_gt=gl_cut,n=int(m.sum()),
       vglut_positive_pct=100*pe[m].mean(),vgat_positive_pct=100*pv[m].mean(),
       double_pct=100*(pv[m]&pe[m]).mean(),vglut_only_pct=100*(~pv[m]&pe[m]).mean(),
       neither_pct=100*(~pv[m]&~pe[m]).mean()))
Q=pd.DataFrame(oldrows);Q.to_csv(D/"stage986_vgat10_vglut_threshold_h5ad_routeA_reproducibility.csv",index=False)
# Summaries computed after all counts
total=T[(T.vgat_thr==10)&(T.section==0)&(T.vglut_thr.isin([2,3,10]))&T.broad.isin(["ALL","neuronal_EI","E","I"])]
print("MAIN_COUNTS",total[["vglut_thr","broad","n","vglut_detect","vgat_detect","vgat_only_pct","vglut_only_pct","co_positive_pct","neither_pct"]].to_string(index=False,float_format=lambda x:f"{x:.3f}"),flush=True)
# Fig 1: sensitivity with percentages, co-fraction and E vs I cutoffs
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
fig,axs=plt.subplots(2,2,figsize=(12,10),layout="constrained")
colorE="#E37B34";colorI="#4387BE"
for lab,clr in (("E",colorE),("I",colorI)):
 t=T[(T.vgat_thr==10)&(T.section==0)&(T.broad==lab)].sort_values("vglut_thr")
 axs[0,0].plot(t.vglut_thr,100*t.vglut_detect,marker="o",ms=4,label=lab,color=clr)
 axs[0,1].plot(t.vglut_thr,t.vglut_only_pct,marker="o",ms=4,label=lab,color=clr)
 axs[1,0].plot(t.vglut_thr,t.co_positive_pct,marker="o",ms=4,label=lab,color=clr)
 axs[1,1].plot(t.vglut_thr,t.neither_pct,marker="o",ms=4,label=lab,color=clr)
titles=["VGLUT2 positive (%)","VGLUT2-only (%)","VGAT/VGLUT2 double-positive (%)","Neither marker (%)"]
for ax,title in zip(axs.flat,titles):
 ax.set(xlabel="VGLUT2 count threshold (> N spots/cell)",ylabel="Cells within broad class (%)",title=title)
 ax.set_ylim(bottom=0);ax.set_xlim(0,20)
 for x in [2,3,10]:ax.axvline(x,color="#888888",lw=.7,alpha=.45,ls="--")
 ax.legend(frameon=False)
fig.suptitle("Stage986: VGAT >10 fixed; VGLUT2 >0–20 sweep (71,950 cells, 1 specimen)",size=14)
fig.savefig(O/"STAGE986_VGLUT2_THRESHOLD_SENSITIVITY.png",dpi=200);fig.savefig(O/"STAGE986_VGLUT2_THRESHOLD_SENSITIVITY.pdf");plt.close(fig)
# Fig 2: co-class counts for 10/2, 10/3, 10/10, broken down by section and all.
palette=["#4387BE","#E37B34","#A362BC","#BABABA"]
names=["VGAT-only","VGLUT2-only","Double+","Neither"]
fig,axs=plt.subplots(3,4,figsize=(15,10),layout="constrained")
for i,gcut in enumerate([2,3,10]):
 for j,section in enumerate([0,500,530,560]):
  ax=axs[i,j]
  t=T[(T.vgat_thr==10)&(T.vglut_thr==gcut)&(T.section==section)&(T.broad=="ALL")].iloc[0]
  portions=np.array([t.vgat_only_pct,t.vglut_only_pct,t.co_positive_pct,t.neither_pct])
  ax.bar(np.arange(4),portions,color=palette,width=.65)
  ax.set_ylim(0,100)
  ax.set_xticks(range(4),["VGAT-only","VGLUT2-only","Double","Neither"],rotation=48,ha="right")
  ax.set_title(f"VGAT >10 / VGLUT2 >{gcut} · "+("ALL" if section==0 else f"S{section}")+f" (n={t.n:,})",fontsize=10)
  for k,y in enumerate(portions):
   if y>1:ax.text(k,y+.8,f"{y:.1f}%",ha="center",fontsize=8)
fig.suptitle("Four classes across sections · exact percentage denominators displayed",fontsize=13)
fig.savefig(O/"STAGE986_FOUR_CLASS_10_2_3_10_PER_SECTION.png",dpi=185);plt.close(fig)
# Fig3: current frozen embedding comparison for both proposed thresholds
X1=np.asarray(Z["embedding"]);X2=np.load(Path(r"G:\Map6_recover_all\stage982_tsne\multiscale_30_300.npz"))["embedding"]
fig,axs=plt.subplots(2,3,figsize=(18,12),layout="constrained")
for j,glcut in enumerate([2,3,10]):
 pv=v>10;pe=e>glcut;cl=np.full(len(v),3,int);cl[pv&~pe]=0;cl[~pv&pe]=1;cl[pv&pe]=2
 for i,(label,coords) in enumerate([("UMAP",X1),("multiscale t-SNE",X2)]):
  ax=axs[i,j]
  for k in (3,2,0,1):
   m=cl==k
   ax.scatter(coords[m,0],coords[m,1],s=.55,color=palette[k],alpha=.13 if k==3 else .66,linewidths=0,rasterized=True)
  ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
  for side in ax.spines.values():side.set_visible(False)
  ax.set_title(f"{label} · VGAT >10, VGLUT2 >{glcut}",fontsize=11)
handles=[Line2D([0],[0],marker="o",linestyle="",markersize=8,color=c,label=n) for c,n in zip(palette,names)]
fig.legend(handles=handles,loc="lower center",ncol=4,frameon=False,bbox_to_anchor=(.5,-.01))
fig.suptitle("Same cell identities, same embedding · marker rule changes only",fontsize=14)
fig.savefig(O/"STAGE986_VGAT10_VGLUT2_2_3_10_UMAP_TSNE.png",dpi=170,bbox_inches="tight")
fig.savefig(O/"STAGE986_VGAT10_VGLUT2_2_3_10_UMAP_TSNE.pdf",bbox_inches="tight")
plt.close(fig)
# Section reproducibility of VGLUT2 positivity using all broad E/I and heldout.
fig,axs=plt.subplots(1,2,figsize=(12,5),layout="constrained")
for j,cat in enumerate(["E","I"]):
 for sec,c in [(500,"#3273C3"),(530,"#D46B30"),(560,"#54A879")]:
  t=T[(T.vgat_thr==10)&(T.section==sec)&(T.broad==cat)].sort_values("vglut_thr")
  axs[j].plot(t.vglut_thr,100*t.vglut_detect,label=f"S{sec}",color=c,marker="o",ms=3)
 axs[j].set(title=f"Broad {cat}: VGLUT2+ % by section",ylabel="Positive cells (%)",xlabel="VGLUT2 count threshold (> N)")
 axs[j].legend(frameon=False)
fig.savefig(O/"STAGE986_VGLUT_SECTION_ROBUSTNESS.png",dpi=185);plt.close(fig)
# Explicit conclusion in JSON, no tuned thresholds to artificially match expectations.
def summary(gcut):
 m=(T.vgat_thr==10)&(T.vglut_thr==gcut)&(T.section==0)
 y=T[m].set_index("broad")
 return {b:{k:float(y.loc[b,k]) for k in ("vglut_detect","vgat_detect","co_positive_pct","neither_pct","vglut_only_pct","vgat_only_pct")} for b in ("ALL","E","I","neuronal_EI")}
audit={"stage":986,"status":"COMPLETE","n_cells":int(len(v)),"n_sections":3,"biological_replicates":1,
"strict_rule":"VGAT count >10; VGLUT2 count >2 or >3; >10/10 historic comparator",
"model":"No reclassification of frozen Fine26 or broad labels. Merely marker-positive overlays.",
"candidate_10_2":summary(2),"candidate_10_3":summary(3),"historic_10_10":summary(10),
"raw_counts":"CURRENT_ROUTEA post-correction; HISTORICAL_H5AD sensitivity separately audited",
"limitations":["A low per-gene detection threshold requires single-molecule background/spot QC; do not label as proven transmitter class without specificity evidence.","Broad E/I identities may depend on same Slc32a1 and Slc17a6 probes (circular validation).","Coexpression does not imply dual transmitter release.","S500/S530/S560 are sections from one mouse, not independent biological repeats."]}
(D/"stage986_marker_threshold_authority.json").write_text(json.dumps(audit,indent=2),encoding="utf8")
print("STAGE986_FINISHED",json.dumps({k:audit[k] for k in ("candidate_10_2","candidate_10_3","historic_10_10")}),flush=True)
