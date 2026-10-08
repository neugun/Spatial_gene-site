"""Stage1018: actual per-section acquisition tiles.json geometry → registered
Stage631 full-section double-positive near real tile OVERLAP stripes, not guessed grid.
Geometry: 25 tile starts+width, 4 neighboring overlap strips per axis.
Affine x/y fitted from Stage631 cell coordinate against matched H5AD array_x/y_px.
Stitch fusion details may shift pixel seams within measured overlap, so test whole strips
and strip center ±10/15/25 and never call optimizer output an exact seam cut.
"""
from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");O=P/"assets"/"stage1018_acquisition_tile_overlaps";O.mkdir(exist_ok=True);D=P/"data"
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["section","roi_id","x_um","y_um","fine26_stable_id"])
base=Path(r"Z:\sternsonlab\Zhenggang\2acq\outputs")
acq={500:"M5L_500_4channels_5X5tile_1",530:"M5L_530_4channels_5X5tile_1",560:"M5L_560_4channels_5X5tile_3"}
hpath=Path(r"G:\PeriLC_current\D_mirror\cluster\Desktop\EASIFISH\PeriLC\_mapmycells_work\query\periLC_500_530_560_27gene_ALLCELL_FINE26_V2_BROADV4_XYFLIP_APDIRECT.h5ad")
h=ad.read_h5ad(hpath,backed="r");xy=h.obs[["section","roi_id","array_x_px","array_y_px"]].reset_index(drop=True);h.file.close()
xy["section"]=xy["section"].astype(int)
M=A.merge(xy,on=["section","roi_id"],how="left",validate="one_to_one")
print("VALID_XPIXELS",int(M.array_x_px.notna().sum()),len(M),flush=True);assert M.array_x_px.notna().sum()>20000
F=pd.read_csv(Path(r"G:\Map6_recover_all\stage1001_neuron_gate_cell_flags_PRIVATE.csv.gz"),usecols=["section","roi_id","reference_neurons"])
M=M.merge(F,on=["section","roi_id"],validate="one_to_one",how="left")
C=np.concatenate([np.load(Path(r"G:\Map6_recover_all\CURRENT_ROUTEA")/f"S{s}_cell_gene27_current.npy") for s in (500,530,560)]).astype(float)
g=list(np.load(R/"stage704_current_umap27"/"CURRENT_UMAP27_STAGE704.npz")["genes"].astype(str))
assert len(M)==C.shape[0]==71950
vg=C[:,g.index("Slc32a1")];glu=C[:,g.index("Slc17a6")]
definitions={"VGAT10_GLU3":(vg>10,glu>3),"VGAT10_GLU10":(vg>10,glu>10),"VGAT15_GLU5":(vg>15,glu>5)}
# true x/y overlap from tile origins, 1920 original mosaic pixels per tile
def strips_from_start(tiles,axis):
 starts=sorted(set(round(float(t["position"][axis]),6) for t in tiles))
 assert len(starts)==5,starts
 widths=[float(t["size"][axis]) for t in tiles]
 assert np.max(widths)-np.min(widths)<1e-6
 width=widths[0]
 overlap=[]
 for i in range(4):
  lower=starts[i+1]-starts[0]
  upper=starts[i]+width-starts[0]
  assert 0<upper-lower<300,(starts,lower,upper)
  overlap.append((lower,upper))
 return overlap
rows=[];geoms={};models=[]
for ss in (500,530,560):
 m=M[M.section==ss].copy()
 slopes={}
 for coord,pixel in [("x_um","array_x_px"),("y_um","array_y_px")]:
  good=m[pixel].notna()
  fit=np.polyfit(m.loc[good,pixel],m.loc[good,coord],1)
  pred=np.polyval(fit,m.loc[good,pixel]);err=np.abs(pred-m.loc[good,coord])
  slopes[coord]=fit
  print("STAGE1018_AFFINE",ss,coord,fit.tolist(),"p99abs",float(np.quantile(err,.99)),"max",float(err.max()),flush=True)
  assert float(np.quantile(err,.99))<1.5,"stage XY and raw array pixel coordinates not aligned"
 tiles=json.loads((base/acq[ss]/"stitching"/"tiles.json").read_text(encoding="utf8"))
 assert len(tiles)==25
 seam={}
 for axis,key,coord in [(0,"x","x_um"),(1,"y","y_um")]:
  raw=strips_from_start(tiles,axis)
  # acquisition output is downsampled 4-fold relative to tile JSON pixel coordinates.
  fit=slopes[coord]
  # Stitch origin may be translated and cropped. Pixel origin offset is determined by
  # installed Stage631~array affine, NOT by tile JSON alone. Report sensitivity below.
  # Using full mosaic origin in array px=0; tile strip raw pixel /4.
  converted=[tuple(sorted(np.polyval(fit,np.array([r0,r1])/4).tolist())) for r0,r1 in raw]
  seam[key]=converted
 geoms[ss]=seam
 models.append(dict(section=ss,n_cells=len(m),tiles=25,
  x_slope=float(slopes["x_um"][0]),x_intercept=float(slopes["x_um"][1]),
  y_slope=float(slopes["y_um"][0]),y_intercept=float(slopes["y_um"][1]),
  overlap_centers_x_um=[float(np.mean(v)) for v in seam["x"]],
  overlap_centers_y_um=[float(np.mean(v)) for v in seam["y"]],
  overlap_widths_x_um=[float(b-a) for a,b in seam["x"]],
  overlap_widths_y_um=[float(b-a) for a,b in seam["y"]]))
 print("STAGE1018_TILE_BOUNDARIES",ss,seam,flush=True)
 for dn,(v,g) in definitions.items():
  vv=v[M.section.to_numpy()==ss];gg=g[M.section.to_numpy()==ss]
  pos=vv&gg;neuron=m.reference_neurons.to_numpy(bool)
  xx=m.x_um.to_numpy(float);yy=m.y_um.to_numpy(float)
  union_stripe=np.zeros(len(m),bool)
  for axis,vals in [("x",xx),("y",yy)]:
   for lo,hi in seam[axis]:
    union_stripe|=(vals>=lo)&(vals<=hi)
  # centers; includes a 4-neighbor seam phase jitter ±10 because original tile optimizer
  dx=np.min(np.abs(xx[:,None]-np.array([np.mean(v) for v in seam["x"]])[None,:]),axis=1)
  dy=np.min(np.abs(yy[:,None]-np.array([np.mean(v) for v in seam["y"]])[None,:]),axis=1)
  center_distance=np.minimum(dx,dy)
  for mode,near in [("native_overlap_strip",union_stripe)]+[(f"center_pm_{v}um",center_distance<=v) for v in (8,15,25)]:
   for gname,gmask in [("all",np.ones(len(m),bool)),("reference_neuron",neuron)]:
    near=near&gmask;far=(~near)&gmask if False else (gmask&~(union_stripe if mode=="native_overlap_strip" else (center_distance<=float(mode.split("_")[2].replace("um","")))))
    if not near.any() or not far.any():continue
    p_near=float(pos[near].mean());p_far=float(pos[far].mean())
    rows.append(dict(section=ss,definition=dn,gate=gname,test=mode,
      n_near=int(near.sum()),n_far=int(far.sum()),both_near=int((pos&near).sum()),both_far=int((pos&far).sum()),
      pct_near=100*p_near,pct_far=100*p_far,risk_ratio=p_near/max(p_far,1e-9),
      diff_pct=(p_near-p_far)*100))
# full section visualization with empirically registered tile-strip overlays
for dn,(v,g) in definitions.items():
 cl=np.full(len(M),3,np.uint8);cl[v&~g]=0;cl[~v&g]=1;cl[v&g]=2
 fig,axs=plt.subplots(1,3,figsize=(16.4,6.2))
 for ax,ss in zip(axs,(500,530,560)):
  mk=M.section.to_numpy()==ss
  colors=["#397ebb","#df8a43","#8d3baf","#d3d6d8"]
  for z in (3,0,1,2):
   mm=mk&(cl==z);ax.scatter(M.x_um[mm],M.y_um[mm],s=1 if z!=2 else 2.3,color=colors[z],alpha=.32 if z==3 else .75,linewidths=0,rasterized=True)
  for lo,hi in geoms[ss]["x"]:
   ax.axvspan(lo,hi,color="#73ad77",alpha=.15,lw=0)
   ax.axvline((lo+hi)/2,color="#3f8b4d",lw=.65,ls="--",alpha=.7)
  for lo,hi in geoms[ss]["y"]:
   ax.axhspan(lo,hi,color="#73ad77",alpha=.15,lw=0)
   ax.axhline((lo+hi)/2,color="#3f8b4d",lw=.65,ls="--",alpha=.7)
  ax.set(xlim=(0,1030),ylim=(0,1030),title=f"S{ss} · {np.sum(cl[mk]==2):,} double+",
       xlabel="Stage631 x (already flipped)",ylabel="Stage631 y (already flipped)")
  ax.set_aspect("equal")
  for sp in ("top","right"):ax.spines[sp].set_visible(False)
 fig.legend(handles=[
    Line2D([0],[0],marker="o",linestyle="",color=c,markersize=7,label=n)
    for c,n in zip(colors,["VGAT-only","VGLUT2-only","Double+","Neither"])]+
    [Line2D([0],[0],linestyle="--",color="#3f8b4d",label="Native tile overlap midline")],
    loc="lower center",ncol=5,bbox_to_anchor=(.5,.04),frameon=False)
 fig.suptitle(f"Native 5×5 tiles.json-derived OVERLAP stripes · {dn} · all cells, full sections",y=.96)
 fig.subplots_adjust(left=.06,right=.985,top=.85,bottom=.22,wspace=.27)
 fig.savefig(O/f"STAGE1018_{dn}_NATIVE_TILE_OVERLAP_SPATIAL.png",dpi=195)
 plt.close(fig)
Z=pd.DataFrame(rows);Z.to_csv(D/"stage1018_native_tile_overlap_dual_rates.csv",index=False)
pd.DataFrame(models).to_csv(D/"stage1018_tile_geometry_and_Stage631_affines.csv",index=False)
fig,axs=plt.subplots(1,3,figsize=(14.5,5.5),layout="constrained")
for ax,ss in zip(axs,(500,530,560)):
 for dn,c in [("VGAT10_GLU3","#8d3baf"),("VGAT10_GLU10","#397ebb"),("VGAT15_GLU5","#df8a43")]:
  t=Z[(Z.section==ss)&(Z.definition==dn)&(Z.gate=="reference_neuron")].set_index("test")
  vals=[t.loc[k,"risk_ratio"] for k in ["native_overlap_strip","center_pm_8um","center_pm_15um","center_pm_25um"]]
  ax.plot(np.arange(4),vals,marker="o",color=c,label=dn)
 ax.axhline(1,color="#777777",ls="--")
 ax.set_xticks(range(4),["Overlap\nstrip","center ±8","center ±15","center ±25"])
 ax.set(title=f"S{ss}",ylabel="Double+ risk near / away from tile seams")
 ax.legend(frameon=False,fontsize=8)
 for sp in ("top","right"):ax.spines[sp].set_visible(False)
fig.suptitle("Acquisition tiles.json-derived overlap, Stage631 physical registration; neuron-reference gated",fontsize=13)
fig.savefig(O/"STAGE1018_NATIVE_TILE_ENRICHMENT_BY_SECTION_AND_THRESHOLD.png",dpi=185);plt.close(fig)
a={"stage":1018,"status":"GEOMETRY_DERIVED_FROM_ACTUAL_TILE_METADATA",
"tile_files":"25 tile starts per section from original stitching/tiles.json, optimizer-final available",
"coordinate_transform":"tile starts normalized to earliest physical position, pixel strip width1920 original; /4 to array px; fit Stage631 x/y affine to original array_x/y_px separately by section",
"important_uncertainty":"tiles.json positions are acquisition initial registration, not optimizer-final per-tile fusion seams; absolute origin translation between mosaic and registered Stage631 array requires full comparison to original stitched TIFF. 8/15/25µm variations are sensitivity controls.",
"scope":"All 71950 cells, original 3 sections one specimen","tests":["whole overlap strip","center±8µm","center±15µm","center±25µm"],
"no_causal_claim":"An edge association is observational; validated individual spot/ROI masks required"}
(D/"stage1018_native_tile_seam_authority.json").write_text(json.dumps(a,indent=2),encoding="utf8")
print("STAGE1018_OVERLAPS",Z.query('definition=="VGAT10_GLU3" and gate=="reference_neuron"').to_string(index=False,float_format=lambda x:f"{x:.3f}"),flush=True)


