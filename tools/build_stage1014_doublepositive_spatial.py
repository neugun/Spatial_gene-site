"""Stage1014: map VGAT / VGLUT2 double-positive across COMPLETE S500/S530/S560.
Test 5x5 estimated tile grid (with phase sensitivity), tissue perimeter, and
cross-transmitter cell-neighborhood border. Never claim inferred grid exact physical seam.
Raw corrected count rules >10/>3, >10/>10, >15/>5; normalized P70 shown as artifact control.
Same specimen three sections, not animal replication; within-section results descriptive.
"""
from pathlib import Path
import os,json
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.api as sm
ROOT=Path(r"G:\Spatial_gene_site_publish")
P=ROOT/"perilc-map6-review";O=P/"assets"/"stage1014_doublepositive_spatial";O.mkdir(exist_ok=True)
D=P/"data";R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
G=np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz");gene=list(G["genes"].astype(str))
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["section","roi_id","x_um","y_um","depth_um","fine26_stable_id"])
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)]).astype(float)
assert len(A)==len(C)==71950 and np.array_equal(A.section.to_numpy(int),G["section"].astype(int))
vg=C[:,gene.index("Slc32a1")];gl=C[:,gene.index("Slc17a6")]
snap=C[:,gene.index("Snap25")] if "Snap25" in gene else np.full(len(C),np.nan)
depth=C.sum(axis=1)
K=pd.read_csv(D/"stage977_fine26_color_key.csv");bmap=dict(zip(K.fine26,K.broad))
broad=A.fine26_stable_id.map(bmap).fillna("Other").to_numpy()
# old conservative gate from stage1001 outputs counts without exposing per-cell in public
prior=Path(r"G:\Map6_recover_all\stage1001_neuron_gate_cell_flags_PRIVATE.csv.gz")
flags=pd.read_csv(prior,usecols=["section","roi_id","reference_neurons","reference_neurons_and_Snap25_ge5"])
assert len(flags)==len(A) and np.array_equal(flags.roi_id.to_numpy(int),A.roi_id.to_numpy(int))
keepref=flags.reference_neurons.to_numpy(bool);keepstrict=flags.reference_neurons_and_Snap25_ge5.to_numpy(bool)
x=A.x_um.to_numpy(float);y=A.y_um.to_numpy(float);z=A.depth_um.to_numpy(float);section=A.section.to_numpy(int)
p70a=np.quantile(np.divide(vg*1e4,depth,out=np.zeros_like(vg),where=depth>0),.7)
p70b=np.quantile(np.divide(gl*1e4,depth,out=np.zeros_like(gl),where=depth>0),.7)
defs={
"raw_VGAT10_VGLUT2_3":(vg>10,gl>3),
"raw_VGAT10_VGLUT2_10":(vg>10,gl>10),
"raw_VGAT15_VGLUT2_5":(vg>15,gl>5),
"normalized_geneP70":(vg*1e4/np.maximum(depth,1)>p70a,gl*1e4/np.maximum(depth,1)>p70b)
}
col=["#397ebb","#df8a43","#8d3baf","#d3d6d8"];catnames=["VGAT-only","VGLUT2-only","Double-positive","Neither"]
rows=[];nearrows=[];audit={}
# 5x5 physical mosaic extents; exact tile seams unknown; include offset phase grid and declare proxy
pitch=1030/5  # downseg cell mosaic ~1029 um wide (2238 px * .46)
print("PIXEL_TO_UM_ESTIMATED_TILE_PITCH",pitch,flush=True)
for dn,(a,b) in defs.items():
 cls=np.full(len(A),3,np.uint8)
 cls[a&~b]=0;cls[~a&b]=1;cls[a&b]=2
 for gate_name,gate in [("all",np.ones(len(A),bool)),("old_reference_neurons",keepref),("reference_neurons_and_Snap255",keepstrict)]:
  for ss in (0,500,530,560):
   m=gate&((section==ss) if ss else np.ones(len(A),bool))
   counts=np.bincount(cls[m],minlength=4);rows.append({"definition":dn,"gate":gate_name,"section":ss,
     "n":int(m.sum()),"VGATonly":int(counts[0]),"VGLUTonly":int(counts[1]),"double":int(counts[2]),"neither":int(counts[3]),
     "double_pct":float(100*counts[2]/max(1,m.sum())),"dual_of_any_pct":float(100*counts[2]/max(1,counts[:3].sum()))})
 # Compute neighborhood border proxy independent of Fine26 label; only exclusive classes
 for ss in (500,530,560):
  m=np.flatnonzero(section==ss)
  p=np.column_stack((x[m],y[m]))
  tree=cKDTree(p)
  dd,nn=tree.query(p,k=11,workers=4)
  local_classes=cls[m]
  ex=((local_classes==0)|(local_classes==1))
  # near other transmitter exclusive neighbor(s): fraction of ten other cells that are a different exclusive class
  n0=(local_classes[nn[:,1:]]==0).mean(axis=1)
  n1=(local_classes[nn[:,1:]]==1).mean(axis=1)
  mixed=2*np.minimum(n0,n1) # zero if no mixture, up to 1
  opp_nearest=np.full(len(m),np.nan)
  ix0=np.flatnonzero(local_classes==0);ix1=np.flatnonzero(local_classes==1)
  if len(ix0)>0 and len(ix1)>0:
   d0=cKDTree(p[ix0]).query(p,k=1,workers=4)[0];d1=cKDTree(p[ix1]).query(p,k=1,workers=4)[0]
   opp_nearest=np.maximum(d0,d1) # distance needed to encounter BOTH exclusive transmitter classes
  # approximate distance to OUTER tissue rectangle; bounding pixels may be missing tissue, flagged.
  wall=np.minimum.reduce([x[m]-0,1030-x[m],y[m]-0,1030-y[m]])
  # 5x5 grid proxy phase from physical array corner, and grid phase sensitivity
  prox_dist=lambda xx,phase:np.abs(((xx-phase+pitch/2)%pitch)-pitch/2)
  seam=np.minimum(prox_dist(x[m],0),prox_dist(y[m],0))
  enrich=local_classes==2
  for basegate,subset in [("all",np.ones(len(m),bool)),("old_reference_neurons",keepref[m]),("reference_neurons_and_Snap255",keepstrict[m])]:
   for feature,vals,cut in [
    ("estimated_tile_seam_distance_um",seam,15),
    ("outer_rect_edge_distance_um",wall,40),
    ("opposing_marker_nearest_both_um",opp_nearest,25),
    ("cross_exclusive_neighbor_mix",mixed,.2)]:
    good=subset&np.isfinite(vals)
    for cat,t in [("near",vals<=cut),("far",vals>cut)]:
     mm=good&t
     nearrows.append(dict(definition=dn,section=ss,gate=basegate,feature=feature,threshold=cut,comparison=cat,n=int(mm.sum()),
        dual=int(np.sum(enrich&mm)),dual_fraction=float(np.mean(enrich[mm])) if mm.any() else np.nan))
  # per-cell internal QC saved PRIVATE for auditable seam hypotheses
  private=Path(r"G:\Map6_recover_all\stage1014_spatial_doublepositive_private");private.mkdir(parents=True,exist_ok=True)
  pd.DataFrame(dict(section=ss,roi_id=A.roi_id.to_numpy()[m],x_um=x[m],y_um=y[m],z_um=z[m],
    VGAT=vg[m],VGLUT2=gl[m],snap25=snap[m],depth27=depth[m],Fine26=A.fine26_stable_id.to_numpy()[m],
    broad=broad[m],class4=local_classes,seam_proxy_um=seam,outer_rect_edge_um=wall,
    mixed_exclusive_neighborhood=mixed,nearest_both_marker_um=opp_nearest)).to_csv(private/f"{dn}_S{ss}_cells.csv.gz",index=False,compression="gzip")
 # Three true FULL section panels, points above colored by class, purple drawn on top
 fig,axs=plt.subplots(1,3,figsize=(17,6.6),facecolor="white")
 for ax,ss in zip(axs,[500,530,560]):
  m=section==ss
  for k in (3,0,1,2):
   j=m&(cls==k)
   ax.scatter(x[j],y[j],s=1.0 if k!=2 else 3.4,c=col[k],linewidths=0,alpha=.42 if k==3 else .83,rasterized=True)
  # approx grid, explicitly not ground-truth stitched tile metadata
  for t in np.arange(pitch,1030,pitch):
   ax.axvline(t,c="#666666",ls="--",lw=.55,alpha=.65);ax.axhline(t,c="#666666",ls="--",lw=.55,alpha=.65)
  ax.set(xlim=(0,1030),ylim=(0,1030),xlabel="Stage631 X already flipped (µm)",ylabel="Stage631 Y already flipped (µm)")
  ax.set_aspect("equal");ax.set_title(f"S{ss} · n={m.sum():,}; dual={100*np.mean(cls[m]==2):.1f}%",fontsize=12)
  for sp in ["top","right"]:ax.spines[sp].set_visible(False)
 fig.suptitle(f"Complete XY sections · {dn.replace('_',' ')} · dashed = estimated 5×5 tile boundaries (NOT independently aligned)",fontsize=14)
 fig.legend(handles=[Line2D([0],[0],linestyle="",marker="o",markersize=8,color=c,label=n) for c,n in zip(col,catnames)],loc="lower center",ncol=4,frameon=False,bbox_to_anchor=(.5,.025))
 fig.subplots_adjust(left=.065,right=.99,top=.84,bottom=.20,wspace=.24)
 fig.savefig(O/f"STAGE1014_{dn}_FULL_SECTION_DUAL_SPATIAL.png",dpi=180,bbox_inches="tight")
 plt.close(fig)
 # per-section class4 overlay without grid for unbiased view
 if dn=="raw_VGAT10_VGLUT2_3":
  fig,axs=plt.subplots(1,3,figsize=(16.5,5.8),layout="constrained")
  for ax,ss in zip(axs,[500,530,560]):
   m=section==ss
   ax.scatter(x[m],y[m],s=.35,color="#e2e5e7",alpha=.40,linewidths=0,rasterized=True)
   dual=m&(cls==2)
   ax.scatter(x[dual],y[dual],s=2.1,color="#8d3baf",alpha=.86,linewidths=0,rasterized=True)
   ax.set(xlim=(0,1030),ylim=(0,1030),title=f"S{ss} · dual {dual.sum():,}");ax.set_aspect("equal")
   ax.set_xticks([]);ax.set_yticks([])
  fig.suptitle("Authentic double-positive cell locations (no inferred tile overlay; full measured sections)")
  fig.savefig(O/"STAGE1014_DUAL_ONLY_FULL_SECTION_NO_GRID.png",dpi=185);plt.close(fig)
pd.DataFrame(rows).to_csv(D/"stage1014_double_positive_counts_all_definitions_gates.csv",index=False)
N=pd.DataFrame(nearrows);N.to_csv(D/"stage1014_double_positive_tile_and_cell_border_proxies.csv",index=False)
# Build relative-risk plots of 3 edge definitions and all thresholds, no p-values implying animals
fig,axs=plt.subplots(1,3,figsize=(15,4.5),layout="constrained")
for ax,feature in zip(axs,["estimated_tile_seam_distance_um","outer_rect_edge_distance_um","opposing_marker_nearest_both_um"]):
 for definition,c in [("raw_VGAT10_VGLUT2_3","#8d3baf"),("raw_VGAT10_VGLUT2_10","#397ebb"),("raw_VGAT15_VGLUT2_5","#df8a43"),("normalized_geneP70","#279b83")]:
  t=N[(N.definition==definition)&(N.gate=="old_reference_neurons")&(N.feature==feature)]
  v=[]
  for ss in (500,530,560):
   q=t[t.section==ss];n=q[q.comparison=="near"].iloc[0];f=q[q.comparison=="far"].iloc[0]
   v.append(n.dual_fraction/max(f.dual_fraction,1e-8))
  ax.plot([500,530,560],v,marker="o",color=c,label=definition.replace("raw_","")[:27])
 ax.axhline(1,color="#777777",ls="--",lw=.8)
 ax.set(xlabel="Section",ylabel="Double-positive near/far fraction ratio",title=feature.replace("_"," "))
 ax.legend(fontsize=7,frameon=False)
fig.suptitle("Predefined edge proxies (estimated 5×5 grid phase; NOT verified actual stitch seams)",fontsize=12)
fig.savefig(O/"STAGE1014_EDGE_AND_NEIGHBOR_ENRICHMENT.png",dpi=190);plt.close(fig)
# tile proxy phase sensitivity: test x/y common offset on 4x regular internal grid, not cherry-picked
phaserows=[]
for dn,(a,b) in defs.items():
 cl=a&b
 for ss in (500,530,560):
  m=section==ss;u=x[m];v=y[m];good=keepref[m]
  for ph in np.arange(0,pitch,5.):
   # nearest vertical OR horizontal estimated seam
   wx=np.abs(((u-ph+pitch/2)%pitch)-pitch/2)
   wy=np.abs(((v-ph+pitch/2)%pitch)-pitch/2)
   near=(np.minimum(wx,wy)<=15)&good
   far=(np.minimum(wx,wy)>15)&good
   phaserows.append(dict(definition=dn,section=ss,offset_um=ph,n_near=int(near.sum()),n_far=int(far.sum()),
      dual_near=float(np.mean(cl[m][near])) if near.any() else np.nan,
      dual_far=float(np.mean(cl[m][far])) if far.any() else np.nan,
      enrichment=(float(np.mean(cl[m][near]))/max(1e-9,float(np.mean(cl[m][far])))) if near.any() and far.any() else np.nan))
phase=pd.DataFrame(phaserows);phase.to_csv(D/"stage1014_tile_phase_sensitivity_scan.csv",index=False)
fig,axs=plt.subplots(1,3,figsize=(15,4.6),layout="constrained")
for ax,ss in zip(axs,[500,530,560]):
 for dn,c in [("raw_VGAT10_VGLUT2_3","#8d3baf"),("raw_VGAT10_VGLUT2_10","#397ebb"),("raw_VGAT15_VGLUT2_5","#df8a43")]:
  t=phase[(phase.definition==dn)&(phase.section==ss)]
  ax.plot(t.offset_um,t.enrichment,label=dn.replace("raw_",""),color=c,lw=1.6)
 ax.axhline(1,color="#777777",ls="--",lw=.8)
 ax.set(title=f"S{ss}",xlabel="Unknown tile-grid phase offset (µm)",ylabel="Dual enrichment within 15µm")
 ax.legend(frameon=False,fontsize=8)
fig.suptitle("Crucial robustness: inferred grid phase could change apparent seam enrichment",fontsize=13)
fig.savefig(O/"STAGE1014_UNCERTAIN_TILE_PHASE_SENSITIVITY.png",dpi=190);plt.close(fig)
# concentration near opposite-type exclusive neighborhoods beyond chance, coarse matched strata
summ=N[(N.definition=="raw_VGAT10_VGLUT2_3")&(N.gate=="old_reference_neurons")]
for feature in summ.feature.unique():
 t=summ[summ.feature==feature]
 print("EDGE_PROXY",feature,[(int(ss),round(float(t[(t.section==ss)&(t.comparison=="near")].dual_fraction.iloc[0]),4),round(float(t[(t.section==ss)&(t.comparison=="far")].dual_fraction.iloc[0]),4)) for ss in (500,530,560)],flush=True)
print("DUAL_COUNTS",pd.DataFrame(rows).query('gate=="old_reference_neurons"').to_string(index=False,float_format=lambda x:f"{x:.3f}"),flush=True)
audit={"stage":1014,"status":"COMPLETE","cell_universe":71950,"per_section":{"500":31461,"530":20513,"560":19976},
"marker_pair":"Slc32a1/Slc17a6 corrected RouteA count equivalents",
"primary":"VGAT>10 and VGLUT2>3; no threshold claimed validated experimentally",
"checks":["VGAT>10/VGLUT2>10","VGAT>15/VGLUT2>5","gene-normalized shared P70 (artifact-sensitive)"],
"edge_geometry":"Estimated 5x5 physical acquisition GRID with 206um spacing and phase sensitivity ±one period; not tile metadata-validated physical stitching seams. True tile-transform/spot-level alignment still required.",
"additional_borders":["outer section bounding rectangle, not tissue actual perimeter","nearest simultaneous exclusive V and G classes","10NN opposing-type mixture"],
"statistical_unit":"three sections from same animal; edge contrasts are descriptive, not independent-animal p-values",
"private_cell_level":"G:/Map6_recover_all/stage1014_spatial_doublepositive_private (not published)",
"biological_limit":"two-gene co-detection does not establish bona fide co-releasing neuron; underlying segmentation/correction needs original punctum validation"}
(D/"stage1014_dual_spatial_authority.json").write_text(json.dumps(audit,indent=2),encoding="utf8")
print("STAGE1014_FINISHED",flush=True)

