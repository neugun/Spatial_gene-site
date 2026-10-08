"""Stage1036 entire 1119x1111 XY field from original acquisition TIFF.
Source TIFF is 8x XY /4x Z reduced from RS-FISH N5 s0.
Visually expanded detection disk denotes a detection (NOT footprint at native scale).
Use Stage1037 N5 s0 for real spot-sized image-supported footprint.
"""
from pathlib import Path
import json,numpy as np,pandas as pd,tifffile
from scipy.ndimage import gaussian_filter
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
R=Path(r"Z:\sternsonlab\Zhenggang\2acq")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");O=P/"assets"/"stage1036_fullXY_spot_context";O.mkdir(parents=True,exist_ok=True);D=P/"data"
genes={
"Snap25":(R/"outputs/M5L_500_4channels_5X5tile_1/R1_C0123.tif",R/"map6_allsections_slurm_20260926/baseline/S500/R1_Snap25_c1/spots_c1.txt",1),
"VGAT":(R/"outputs/M5L_500_4channels_5X5tile_R77/R77_C0-3.tif",R/"map6_allsections_slurm_20260926/baseline/S500/R7_Slc32a1_c0/spots_c0.txt",0),
"VGLUT2":(R/"outputs/M5L_500_4channels_5X5tile_R99/R99_C0-3.tif",R/"map6_allsections_slurm_20260926/corrected_spots/S500/R9_Slc17a6_c3/spots_c3.txt",3)}
# 3 TIFF Z center planes (intervals 15.12µm); 7 TIFF planes ~11.76µm MIP
centers=(150,175,200); half=3
summ=[]
for gene,(imfile,spfile,ch) in genes.items():
 data=[]
 for chunk in pd.read_csv(spfile,header=None,usecols=[0,1,2],dtype=np.float32,chunksize=350000):
  X=chunk.iloc[:,0].to_numpy()/1.84;Y=chunk.iloc[:,1].to_numpy()/1.84;Z=chunk.iloc[:,2].to_numpy()/1.68
  use=np.zeros(len(chunk),bool)
  for z in centers:use|=(abs(Z-z)<=half+.6)
  if use.any():data.append(np.column_stack([X[use],Y[use],Z[use]]))
 spots=np.vstack(data) if data else np.zeros((0,3))
 with tifffile.TiffFile(imfile) as tif:
  dims=tif.series[0].shape; assert len(dims)==4,dims
  print("INPUT",gene,"TIFF",dims,"selected spots",len(spots),flush=True)
  for z in centers:
   planes=[tif.pages[(zz*4)+ch].asarray().astype(np.float32) for zz in range(z-half,z+half+1)]
   raw=np.max(np.stack(planes),axis=0)
   bg=gaussian_filter(raw,2.8)
   residual=np.maximum(raw-bg,0)
   low=float(np.quantile(raw,.03));hi=float(np.quantile(raw,.998))
   enh=np.clip((raw-low)/max(hi-low,1),0,1)
   e99=float(np.quantile(residual,.997))
   hp=np.clip(residual/max(e99,1),0,1)
   use=np.abs(spots[:,2]-z)<=half+.6
   pts=spots[use]
   fig,axs=plt.subplots(1,3,figsize=(15.4,5.6))
   axs[0].imshow(raw,cmap="gray",vmin=low,vmax=hi,interpolation="nearest")
   axs[1].imshow(hp,cmap="gray",vmin=0,vmax=1,interpolation="nearest")
   axs[2].imshow(enh,cmap="gray",vmin=0,vmax=1,interpolation="nearest")
   if len(pts):
    # 2px visually expanded circles at 8x downsample (NOT real footprint)
    axs[2].scatter(pts[:,0],pts[:,1],marker="o",s=5.2,linewidths=.40,edgecolor="#eb5683",facecolor="none",alpha=.82,rasterized=True)
   for ax in axs:
    ax.set(xlim=(0,dims[-1]),ylim=(dims[-2],0))
    ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
    for spine in ax.spines.values():spine.set_visible(False)
   axs[0].set_title("Full XY raw acquisition · N5 downsampled preview")
   axs[1].set_title("Full XY high-pass signal (display only)")
   axs[2].set_title(f"Full XY raw + all RS-FISH detections, n={len(pts):,}")
   fig.suptitle(f"S500 FULL acquisition XY · {gene} · TIFF z={z}±{half} (11.8µm slab) · Z depends on round\n"
    "Ring radius expanded for screen visibility, NOT an actual measured spot footprint; Stage1037 uses native N5",fontsize=12)
   fig.subplots_adjust(left=.01,right=.99,top=.83,bottom=.015,wspace=.07)
   out=O/f"STAGE1036_S500_{gene}_FULLXY_TIFFz{z:03d}_RAW_HP_SPOTOVERLAY.png"
   fig.savefig(out,dpi=175,bbox_inches="tight");plt.close(fig)
   summ.append({"gene":gene,"section":500,"tiff_z_center":z,"z_halfwidth_tiff_planes":half,
    "n_detected_spots_overlay":len(pts),"full_XY_pixels":f"{dims[-1]}x{dims[-2]}", "source":"original acquisition TIFF 8xXY/4xZ downsample from N5 s0"})
   print("FINISHED",gene,z,len(pts),flush=True)
pd.DataFrame(summ).to_csv(D/"stage1036_fullXY_spots_by_gene_and_Z.csv",index=False)
(D/"stage1036_fullXY_provenance.json").write_text(json.dumps({"stage":1036,"scope":"S500 entire XY acquisition field in each of R1/R77/R99 raw channel coordinate systems",
"raw_image":"original acquisition TIFF is downsampled XY 8x and Z 4x relative to N5 s0. Z slabs z±3 TIFF frames; spots from own gene source coordinates /1.84um and /1.68um",
"overlay":"color rings are enlarged for display; NOT physical spot radius; scientifically meaningful footprints in high-res Stage1037",
"registration":"Rounds are shown each in OWN acquisition image coordinates. Their nominal XY fields are not transformed into one-to-one cell alignment. Do not conclude same-cell co-localization from these full fields.",
"expected_outputs":len(summ)},indent=2),encoding="utf8")
print("STAGE1036_FULLXY_COMPLETE",len(summ),flush=True)
