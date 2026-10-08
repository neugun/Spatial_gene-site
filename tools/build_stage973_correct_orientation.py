"""Stage973: Stage631 pre-flipped coordinates; smooth inner and OUTER anatomy.
Display contours are smoothed only; Stage943 expression and Stage956 region assignments stay frozen.
"""
from pathlib import Path
import json
import numpy as np
from scipy.ndimage import gaussian_filter, label as cc_label, binary_fill_holes, distance_transform_edt, zoom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from PIL import Image
from map10_display_xy import flip_xy,verify_flip_xy
from map10_region_geometry import refine_region10,topology_counts
import pandas as pd
root=Path(r"G:\Spatial_gene_site_publish")
pub=root/"perilc-map6-review"; data=pub/"data"
out=pub/"assets"/"stage973_correct_orientation_roundedge"
out.mkdir(parents=True,exist_ok=True)
script=(root/"tools"/"build_stage956_final_atlas.py").read_text(encoding="utf8")
init=script.split('pd.DataFrame(sweep_rows).to_csv')[0]
context={"__name__":"stage963"}
exec(compile(init,"stage956_frozen_reconstruct","exec"),context)
infos=context["region_info"]
cells=context["cells"]
X=context["X"]; genes=context["genes"]; secs=(500,530,560)
xraw=context["x0"]; yraw=context["y0"]; section=context["section"]
xd=context["xd"]; yd=context["yd"]
# Important: Stage631 coordinates are already calibrated XY-flipped
# physical coordinates; prior Stage956 used a second flip.
# Retain region identity by reprojecting its numeric categorical grid.
# This is NOT a raster-image transform.
xd=xraw.copy();yd=yraw.copy()
display={}
diagnostics={}
for s in secs:
    original=infos[s]["local"][::-1,::-1].copy()
    support=infos[s]["mask"][::-1,::-1].copy()
    cleaned,mask,qa=refine_region10(original,support)
    display[s]=dict(region=cleaned,mask=mask,original=original,sigma=qa["smoothing_sigma"])
    diagnostics[str(s)]=qa
    assert qa["after"]["white_holes"]==0
    assert all(n==1 for n in qa["after"]["region_components"].values())
    assert all(n==0 for n in qa["after"]["region_holes"].values())
print("STAGE973_TOPOLOGY_PASSED",json.dumps(diagnostics),flush=True)
tab=plt.get_cmap("tab10")
def draw_bounds(ax,s,fill=False):
    # The outer support is smoothed by a Gaussian contour, not jagged imshow
    # edge pixels. Internal category fields are evaluated at 4x subgrid resolution.
    info=infos[s];G=display[s]["region"];xx=info["xc"];yy=info["yc"]
    n=4; dx=float(xx[1]-xx[0]);dy=float(yy[1]-yy[0])
    xe=(float(xx[0]-dx/2),float(xx[-1]+dx/2))
    ye=(float(yy[0]-dy/2),float(yy[-1]+dy/2))
    mask=display[s]["mask"]
    smooth_support=gaussian_filter(mask.astype(float),3.7)
    H=zoom(smooth_support,n,order=3)
    high_support=H>=.50
    high_support=binary_fill_holes(high_support)
    xs=np.linspace(xe[0]+dx/(2*n),xe[1]-dx/(2*n),H.shape[1])
    ys=np.linspace(ye[0]+dy/(2*n),ye[1]-dy/(2*n),H.shape[0])
    if fill:
        prob=np.stack([zoom(gaussian_filter((G==k).astype(float),1.05),n,order=3) for k in range(1,11)])
        categorical=prob.argmax(axis=0).astype(int)
        rgba=np.empty((*categorical.shape,4),dtype=float)
        palette=np.array([tab(k) for k in range(10)])
        rgba[:,:,:3]=palette[categorical,:3]
        # Soft anti-aliased cutout of external border. No pixel-cell imshow edges.
        alpha=np.clip(gaussian_filter(high_support.astype(float),.8)*2.2-.65,0,1)
        rgba[:,:,3]=alpha*.75
        ax.imshow(rgba,extent=(*xe,*ye),origin="lower",interpolation="bilinear",aspect="equal")
    # Single smooth external contour, internal labels retain source-cell geometry.
    ax.contour(xs,ys,H,levels=[.50],colors=["#444"],linewidths=.8 if fill else .46,alpha=.9 if fill else .42)
    for k in range(1,11):
        z=gaussian_filter((G==k).astype(float),1.35)
        z=np.where(smooth_support>.92,z,np.nan)
        if np.nanmax(z)>.50:
            ax.contour(xx,yy,z,levels=[.5],colors=["#444"],linewidths=.65 if fill else .4,alpha=.78 if fill else .26)
def clean_axes(ax):
 ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
 for sp in ax.spines.values():sp.set_visible(False)
fig,axs=plt.subplots(1,3,figsize=(13.2,4.95))
for ax,s in zip(axs,secs):
 draw_bounds(ax,s,True)
 G=display[s]["region"]
 for k in range(1,11):
  ys,xs=np.where(G==k)
  if len(xs):ax.text(float(infos[s]["xc"][xs].mean()),float(infos[s]["yc"][ys].mean()),str(k),ha="center",va="center",fontsize=9,bbox=dict(boxstyle="circle,pad=.2",fc="white",ec="#666",lw=.6))
 ax.set_title(f"S{s} · local regions 1–10",fontsize=11);clean_axes(ax)
fig.suptitle("Marker-guided region10 | hole-free domains | x flipped, y flipped",fontsize=14)
fig.tight_layout(rect=[0,0,1,.94])
fig.savefig(out/"FINAL_REGION10_STAGE631_PREFLIPPED_SMOOTH_OUTER.png",dpi=180,bbox_inches="tight")
plt.close(fig)
print("STAGE973_REGION_IMAGE_READY",flush=True)
oldman=pd.read_csv(data/"stage956_gene_section_manifest.csv")
policy=pd.read_csv(data/"stage956_gene_colorbar_policy.csv").set_index("gene")
rows=[]
for gi,g in enumerate(genes):
 vmax=float(policy.loc[g,"vmax"]); gamma=float(policy.loc[g,"gamma"])
 norm=PowerNorm(gamma=gamma,vmin=0,vmax=vmax)
 for s in secs:
  mask=section==s
  v=X[mask,gi];xp=xd[mask];yp=yd[mask]
  fig,ax=plt.subplots(figsize=(5.1,4.95))
  ax.scatter(xp,yp,s=.45,c="#e7e7e7",alpha=.55,linewidths=0,rasterized=True)
  active=v>0
  col=ax.scatter(xp[active],yp[active],s=1.55,c=np.clip(v[active],0,vmax),cmap="Reds",norm=norm,alpha=.90,linewidths=0,rasterized=True)
  draw_bounds(ax,s,False)
  clean_axes(ax);ax.set_title(f"{g} · S{s}",fontsize=11)
  cb=fig.colorbar(col,ax=ax,fraction=.048,pad=.026)
  cb.set_ticks([0,vmax/2,vmax])
  cb.set_ticklabels(["0",f"{vmax/2:.0f}" if vmax>=10 else f"{vmax/2:.1f}",f"{vmax:.0f}" if vmax>=10 else f"{vmax:.1f}"])
  cb.set_label("spot count",fontsize=9);cb.ax.tick_params(labelsize=8)
  fig.text(.5,.015,f"display 0–{vmax:.1f} spot count | x flipped, y flipped",ha="center",fontsize=8,color="#555")
  fig.tight_layout(rect=[0,.035,1,1])
  temp=out/"temp_working.png"
  fig.savefig(temp,dpi=180,bbox_inches="tight")
  plt.close(fig)
  fname=f"{g}_S{s}_expression_preview.jpg"
  im=Image.open(temp).convert("RGB");im.thumbnail((850,850));im.save(out/fname,quality=87)
  temp.unlink()
  rows.append(dict(gene=g,section=s,file=fname,x_flipped=True,y_flipped=True,extra_flip=False,vmax=vmax,gamma=gamma,smooth_sigma=display[s]["sigma"]))
 print("STAGE973_GENE",g,flush=True)
pd.DataFrame(rows).to_csv(data/"stage973_gene_section_manifest.csv",index=False)
auth=dict(stage=973,status="PREFLIPPED_COORDINATE_AND_SMOOTH_OUTER_QA_PASSED",
 expression_authority="Stage943 Route A",biological_region_authority="Stage956 unchanged",
 display="already pre-flipped Stage631 x_um/y_um; no second flip of coordinates",
 smooth_method="numeric reproject of frozen Stage956 label grid into Stage631-preflipped coordinates, hole-free segmentation, Gaussian outer contour sigma 3.7; per-cell identities unchanged",
 region_QA=diagnostics,gene_maps=len(rows))
(data/"stage973_correct_orientation_region_audit.json").write_text(json.dumps(auth,indent=2),encoding="utf8")
print("STAGE973_COMPLETE",json.dumps(auth),flush=True)
