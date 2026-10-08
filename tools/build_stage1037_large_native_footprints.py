"""Stage1037 true N5 full-resolution LARGE 140um three-gene RS-FISH overlays.
Original N5 /c/s0 voxels; support footprint estimated locally from raw intensity,
NOT RS-FISH-detected radius (RS-FISH exports only centroid position).
Per-round source image coordinates, NOT cross-round registered cell co-localization.
"""
from pathlib import Path
import json,gzip,struct,functools
import numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
from scipy.ndimage import gaussian_filter,label as cc_label
from skimage.measure import find_contours
R=Path(r"Z:\sternsonlab\Zhenggang\2acq")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");O=P/"assets"/"stage1037_large_native_footprints";O.mkdir(parents=True,exist_ok=True);D=P/"data"
S=json.loads((D/"stage1033_nativeN5_overlay_provenance.json").read_text(encoding="utf8"))["samples"]
regions=[S[0],S[1],S[4]]
inputs={
"Snap25":(R/"outputs/M5L_500_4channels_5X5tile_1/stitching/export.n5/c1/s0",R/"map6_allsections_slurm_20260926/baseline/S500/R1_Snap25_c1/spots_c1.txt"),
"VGAT":(R/"outputs/M5L_500_4channels_5X5tile_R77/stitching/export.n5/c0/s0",R/"map6_allsections_slurm_20260926/baseline/S500/R7_Slc32a1_c0/spots_c0.txt"),
"VGLUT2":(R/"outputs/M5L_500_4channels_5X5tile_R99/stitching/export.n5/c3/s0",R/"map6_allsections_slurm_20260926/corrected_spots/S500/R9_Slc17a6_c3/spots_c3.txt")}
colors={"Snap25":"#24debd","VGAT":"#ffc052","VGLUT2":"#f16bac"}
@functools.lru_cache(maxsize=48)
def block(base,ix,iy,iz):
 f=Path(base)/str(ix)/str(iy)/str(iz)
 if not f.is_file():return None
 a=f.read_bytes();mode,n=struct.unpack(">HH",a[:4]);ds=struct.unpack(">"+str(n)+"I",a[4:4+4*n])
 q=np.frombuffer(gzip.decompress(a[4+4*n:]),dtype=">u2")
 return q.reshape((ds[2],ds[1],ds[0]))
def crop(base,cx,cy,cz,rad=300,zhalf=7):
 cx=int(round(cx));cy=int(round(cy));cz=int(round(cz))
 x0=cx-rad;y0=cy-rad;z0=cz-zhalf;x1=cx+rad;y1=cy+rad;z1=cz+zhalf+1
 vol=np.zeros((2*zhalf+1,2*rad,2*rad),dtype=np.uint16)
 for zk in range(z0//70,(z1-1)//70+1):
  for yk in range(y0//128,(y1-1)//128+1):
   for xk in range(x0//128,(x1-1)//128+1):
    v=block(str(base),xk,yk,zk)
    if v is None:continue
    xx0=max(x0,xk*128);xx1=min(x1,xk*128+v.shape[2])
    yy0=max(y0,yk*128);yy1=min(y1,yk*128+v.shape[1])
    zz0=max(z0,zk*70);zz1=min(z1,zk*70+v.shape[0])
    if min(xx1-xx0,yy1-yy0,zz1-zz0)<=0:continue
    vol[zz0-z0:zz1-z0,yy0-y0:yy1-y0,xx0-x0:xx1-x0]=v[zz0-zk*70:zz1-zk*70,yy0-yk*128:yy1-yk*128,xx0-xk*128:xx1-xk*128]
 return vol.astype("float32"),(x0,y0,z0)
# Record ALL source detections in big xy neighborhoods and broad z range, then select focal bands per gene.
spotsets={}
for gene,(n5,f) in inputs.items():
 chunks=[]
 for a in pd.read_csv(f,header=None,usecols=[0,1,2],dtype=np.float32,chunksize=275000):
  xx=a.iloc[:,0].to_numpy()/.23; yy=a.iloc[:,1].to_numpy()/.23; zz=a.iloc[:,2].to_numpy()/.42
  use=np.zeros(len(a),bool)
  for reg in regions:
   cx=reg["cx"]*8;cy=reg["cy"]*8;cz=reg["cz"]*4
   use|=(np.abs(xx-cx)<340)&(np.abs(yy-cy)<340)&(np.abs(zz-cz)<100)
  if use.any():chunks.append(np.column_stack([xx[use],yy[use],zz[use]]))
 spotsets[gene]=np.vstack(chunks) if len(chunks) else np.zeros((0,3))
 print("LOADED_GENE",gene,len(spotsets[gene]),flush=True)
def boundary(img,x,y,zcenter,z0):
 """Find local compact raw-intensity support near RS-FISH centroid, at spot's source z."""
 iy=int(round(y));ix=int(round(x));iz=int(round(zcenter-z0))
 if not(12<=ix<img.shape[2]-12 and 12<=iy<img.shape[1]-12 and 1<=iz<img.shape[0]-2):return None
 roi=img[iz-1:iz+2,iy-11:iy+12,ix-11:ix+12].max(axis=0)
 # local patch 23x23, central 5x5 peak; annulus for background
 yy,xx=np.indices(roi.shape)
 r=((yy-11)**2+(xx-11)**2)**.5
 ann=roi[(r>=8)&(r<=11)]
 bg=float(np.median(ann));mad=float(1.4826*np.median(abs(ann-bg))+1.)
 small=roi[9:14,9:14]; py,px=np.unravel_index(np.argmax(small),small.shape);py+=9;px+=9
 amp=float(roi[py,px]-bg)
 if amp<max(4*mad,10.):return None
 thr=bg+max(2.2*mad,.22*amp)
 b=(roi>=thr)&(r<=9.)
 seg,n=cc_label(b,np.ones((3,3)))
 lab=int(seg[py,px])
 if lab==0:return None
 component=seg==lab
 if component.sum()<2 or component.sum()>160:return None
 curves=find_contours(component.astype(float),.5)
 if not curves:return None
 # polygon in crop local native pixels, contours points (y,x)
 curve=max(curves,key=len)
 xv=ix-11+curve[:,1];yv=iy-11+curve[:,0]
 radius=(component.sum()/np.pi)**.5
 return np.column_stack([xv,yv]),radius,float(amp/max(mad,1e-6))
rows=[]
for ri,reg in enumerate(regions):
 for gene,(n5,f) in inputs.items():
  cx=reg["cx"]*8;cy=reg["cy"]*8;cz=reg["cz"]*4
  z=spotsets[gene];eligible=(np.abs(z[:,0]-cx)<290)&(np.abs(z[:,1]-cy)<290)&(np.abs(z[:,2]-cz)<75)
  nearby=z[eligible]
  # Choose z with most spots inside a ±7 native-z slab; require same source gene round
  if len(nearby)>3:
   candidate=np.arange(max(9,cz-60),cz+60,4)
   density=np.array([np.sum(np.abs(nearby[:,2]-zz)<=7) for zz in candidate])
   zi=float(candidate[np.argmax(density)])
  else:zi=float(cz)
  vol,(x0,y0,z0)=crop(n5,cx,cy,zi,rad=300,zhalf=7)
  raw=vol.max(axis=0)
  # RAW display has no background subtraction. High-pass for diagnostics ONLY.
  clean=np.maximum(0,raw-gaussian_filter(raw,6))
  hi=float(np.quantile(raw,.999));lo=float(np.quantile(raw,.08))
  rawnorm=np.clip((raw-lo)/max(1,hi-lo),0,1)
  hhi=float(np.quantile(clean,.996));enh=np.clip(clean/max(1,hhi),0,1)
  good=(nearby[:,0]>=x0+10)&(nearby[:,0]<x0+590)&(nearby[:,1]>=y0+10)&(nearby[:,1]<y0+590)&(np.abs(nearby[:,2]-zi)<=7.1)
  pts=nearby[good]
  masks=[];radii=[];snrs=[]
  for x,y,zsp in pts:
   obj=boundary(vol,x-x0,y-y0,zsp,z0)
   if obj is not None:
    xy,radius,snr=obj;masks.append(xy);radii.append(radius);snrs.append(snr)
  fig,axs=plt.subplots(1,3,figsize=(16.0,5.7))
  axs[0].imshow(rawnorm,cmap="gray",vmin=0,vmax=1,interpolation="nearest")
  axs[1].imshow(enh,cmap="gray",vmin=0,vmax=1,interpolation="nearest")
  axs[2].imshow(rawnorm,cmap="gray",vmin=0,vmax=1,interpolation="nearest")
  cc=colors[gene]
  for poly in masks:
   axs[2].add_patch(Polygon(poly,closed=True,facecolor=(*matplotlib.colors.to_rgb(cc),.22),edgecolor=cc,linewidth=1.4))
  # represent below-threshold / merged detections by thin ring only, never imply physical segmentation.
  # Show all detected source centroids as white one-pixel dots separately from footprints
  if len(pts):
   axs[2].scatter(pts[:,0]-x0,pts[:,1]-y0,s=3,c="white",alpha=.80,linewidth=0,rasterized=True)
  for ax in axs:
   ax.set(xlim=(0,600),ylim=(600,0));ax.set_aspect("equal")
   ax.set_xticks([]);ax.set_yticks([])
   for sp in ax.spines.values():sp.set_visible(False)
  axs[0].set_title("Unmodified N5 raw (display contrast only)")
  axs[1].set_title("Local-background high-pass (diagnostic)")
  axs[2].set_title(f"RS-FISH centers + image-derived full punctum support\n{len(pts)} detections; {len(masks)} supported footprints")
  fig.suptitle(f"S500 · {gene} · ~138 × 138 µm raw FOV · source z={zi:.0f}±7 voxels (5.9 µm)\n"
  f"Region {ri+1}: {reg['category']} reference ROI {reg['roi_id']}; per-round source coordinates (not warped R1 colocalization)",fontsize=12)
  fig.subplots_adjust(left=.01,right=.99,top=.81,bottom=.04,wspace=.06)
  path=O/f"STAGE1037_S500_REGION{ri+1}_{gene}_LARGE_NATIVE_FOOTPRINT.png"
  fig.savefig(path,dpi=150,bbox_inches="tight");plt.close(fig)
  rows.append({"section":500,"region":ri+1,"source_category":reg["category"],"reference_roi":int(reg["roi_id"]),
  "gene":gene,"source_native_Z_center":zi,"native_FOV_pixels":600,"FOV_um":138.,"z_slab_native_planes":15,
  "n_spots_in_cropped_slab":int(len(pts)),"n_bright_image_supported_footprints":len(masks),
  "pct_with_image_footprint":100*len(masks)/max(1,len(pts)),
  "median_estimated_support_radius_native_px":float(np.median(radii)) if radii else np.nan,
  "median_local_SNR_footprint":float(np.median(snrs)) if snrs else np.nan,
  "filename":path.name})
  print("STAGE1037_READY",path.name,"spots",len(pts),"supports",len(masks),flush=True)
pd.DataFrame(rows).to_csv(D/"stage1037_native_large_fov_footprint_qc.csv",index=False)
(D/"stage1037_native_footprints_authority.json").write_text(json.dumps({
 "stage":1037,"status":"DONE","section":500,"fields":3,"genes":list(inputs),
 "source":"Original N5 s0 fluorescence (0.23um XY, 0.42um Z), three independent imaging rounds, local image-supported footprints generated from original spot detections",
 "punctum_support_method":"For each detected xyz within 15-plane slab: native 3-slice MIP around spot Z; 23x23 neighborhood; robust annular background (8-11 px); 2.2 MAD or 22% peak excess threshold; connected footprint near spot center; keep 2-160 px area; trace actual raw-bright connected pixels. White 1px dot marks *all* centroids; colored filled outlines mark supported spot footprints. NOT RS-FISH fitted point spread function, NOT measured biological cell boundary.",
 "per_gene_Z_planes":"Source round Z center adaptively chosen to contain max nearby detections within ±60 native Z. Independent per gene so panels must not be interpreted as identical physical z planes or neuron colocalization.",
 "controls":"Raw and locally background-subtracted enhanced show same exact source voxel field. Diffuse VGAT foci are visible for focus/pointspread review. S500 only for this first batch.",
 "caveat":"Source round R77/R99 not transformed to R1 reference geometry, so spatial coincidence with a segmented R1 soma is NOT established. Annotated ROI identifies approximate region only. Source_R9 corrected spot catalog overlaid on uncorrected source N5."
},indent=2),encoding="utf8")
print("STAGE1037_ALL_FINISHED",len(rows),flush=True)

