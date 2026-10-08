"""Stage963: regenerate full 3-section atlas with a single dual-axis display transform.
Display contours are smoothed only; Stage943 expression and Stage956 region assignments stay frozen.
"""
from pathlib import Path
import json
import numpy as np
from scipy.ndimage import gaussian_filter, label as cc_label, binary_fill_holes, distance_transform_edt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from PIL import Image
from map10_display_xy import flip_xy,verify_flip_xy
import pandas as pd
root=Path(r"G:\Spatial_gene_site_publish")
pub=root/"perilc-map6-review"; data=pub/"data"
out=pub/"assets"/"stage963_orientation_smooth"
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
for s in secs:
 m=section==s
 verify_flip_xy(xraw[m],yraw[m],xd[m],yd[m])
display={}
diagnostics={}
for s in secs:
 original=infos[s]["local"]; support=infos[s]["mask"].copy()
 cc,n=cc_label(support)
 sizes=np.bincount(cc.ravel())
 main=np.argmax(sizes[1:])+1
 mainmask=(cc==main)
 removed=int((support & ~mainmask).sum())
 # tiny interior holes in display only; keeps any non-tiny holes unchanged
 holes=binary_fill_holes(mainmask)&(~mainmask)
 cc_h,nh=cc_label(holes)
 sizes_h=np.bincount(cc_h.ravel())
 small=np.isin(cc_h,np.where((sizes_h<=16)&(np.arange(len(sizes_h))>0))[0])
 mask=mainmask|small
 labels=original.copy()
 dist,ind=distance_transform_edt(labels==0,return_indices=True)
 labels[(labels==0)&mask]=labels[tuple(i[(labels==0)&mask] for i in ind)]
 best=None
 for sig in (0.8,1.2,1.6,2.0,2.4,3.0):
  prob=np.stack([gaussian_filter(((labels==k)&mask).astype(float),sig) for k in range(1,11)])
  new=prob.argmax(axis=0)+1
  agree=float((new[mainmask]==original[mainmask]).mean())
  ious=[]
  for k in range(1,11):
   old=(original==k)&mainmask
   newk=(new==k)&mainmask
   ious.append(float((old&newk).sum()/max((old|newk).sum(),1)))
  ok=(agree>=.97 and min(ious)>=.89)
  if ok:best=(sig,new,agree,min(ious))
 assert best is not None,f"S{s}: no shape-preserving smoothing"
 sig,new,agree,minio=best
 new[~mask]=0
 display[s]=dict(region=new,mask=mask,original=original,sigma=sig)
 diagnostics[str(s)]=dict(sigma=sig,agreement=agree,min_region_iou=minio,removed_detached_mask_pixels=removed,filled_small_hole_pixels=int(small.sum()))
print("SMOOTH_DIAGNOSTICS",json.dumps(diagnostics))
tab=plt.get_cmap("tab10")
def draw_bounds(ax,s,fill=False):
 info=infos[s];G=display[s]["region"];xx=info["xc"];yy=info["yc"]
 for k in range(1,11):
  zz=gaussian_filter((G==k).astype(float),.8)
  if fill:ax.contourf(xx,yy,zz,levels=[.5,1.5],colors=[tab(k-1)],alpha=.45)
  ax.contour(xx,yy,zz,levels=[.5],colors=["#444"],linewidths=.8 if fill else .42,alpha=.85 if fill else .27)
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
fig.suptitle("Marker-guided region10 | smooth contours | x flipped, y flipped",fontsize=14)
fig.tight_layout(rect=[0,0,1,.94])
fig.savefig(out/"FINAL_REGION10_SMOOTH_X_FLIPPED_Y_FLIPPED.png",dpi=180,bbox_inches="tight")
plt.close(fig)
print("STAGE963_REGION_IMAGE_READY",flush=True)
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
  rows.append(dict(gene=g,section=s,file=fname,x_flipped=True,y_flipped=True,vmax=vmax,gamma=gamma,smooth_sigma=display[s]["sigma"]))
 print("STAGE963_GENE",g,flush=True)
pd.DataFrame(rows).to_csv(data/"stage963_gene_section_manifest.csv",index=False)
auth=dict(stage=963,status="DISPLAY_READY_PENDING_VISUAL_QA",
 expression_authority="Stage943 Route A",biological_region_authority="Stage956 unchanged",
 display="x flipped, y flipped by section: x_display=xmin+xmax-x_raw; y_display=ymin+ymax-y_raw",
 smooth_method="display-only smoothing, unchanged Stage956 cell labels",
 region_QA=diagnostics,gene_maps=len(rows))
(data/"stage963_orientation_smoothing_audit.json").write_text(json.dumps(auth,indent=2),encoding="utf8")
print("STAGE963_COMPLETE",json.dumps(auth),flush=True)
