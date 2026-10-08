"""Stage1041 independent candidate bright-puncta screen in source raw N5.
Compares unsupervised 3D high-pass raw peak candidates to existing RSFISH centroids,
avoids treating image-derived peaks as validated transcripts or recall ground truth.
"""
from pathlib import Path
import json,gzip,struct,functools,numpy as np,pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from skimage.feature import peak_local_max
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r"Z:\sternsonlab\Zhenggang\2acq")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data";O=P/"assets"/"stage1041_peak_recall_screen";O.mkdir(exist_ok=True)
meta=pd.read_csv(D/"stage1037_native_large_fov_footprint_qc.csv")
info={
"Snap25":(R/"outputs/M5L_500_4channels_5X5tile_1/stitching/export.n5/c1/s0",R/"map6_allsections_slurm_20260926/baseline/S500/R1_Snap25_c1/spots_c1.txt"),
"VGAT":(R/"outputs/M5L_500_4channels_5X5tile_R77/stitching/export.n5/c0/s0",R/"map6_allsections_slurm_20260926/baseline/S500/R7_Slc32a1_c0/spots_c0.txt"),
"VGLUT2":(R/"outputs/M5L_500_4channels_5X5tile_R99/stitching/export.n5/c3/s0",R/"map6_allsections_slurm_20260926/corrected_spots/S500/R9_Slc17a6_c3/spots_c3.txt")}
sel=json.loads((D/"stage1033_nativeN5_overlay_provenance.json").read_text(encoding="utf8"))["samples"][0]
@functools.lru_cache(maxsize=48)
def block(base,ix,iy,iz):
 p=Path(base)/str(ix)/str(iy)/str(iz)
 if not p.exists():return None
 b=p.read_bytes();mode,nd=struct.unpack(">HH",b[:4]);dims=struct.unpack(">"+str(nd)+"I",b[4:4+4*nd])
 return np.frombuffer(gzip.decompress(b[4+4*nd:]),dtype=">u2").reshape((dims[2],dims[1],dims[0]))
def roi(base,midX,midY,midZ,half=300,zh=7):
 ix=int(midX);iy=int(midY);iz=int(round(midZ))
 x0=ix-half;x1=ix+half;y0=iy-half;y1=iy+half;z0=iz-zh;z1=iz+zh+1
 vol=np.zeros((2*zh+1,2*half,2*half),dtype=np.uint16)
 for kk in range(z0//70,(z1-1)//70+1):
  for jj in range(y0//128,(y1-1)//128+1):
   for ii in range(x0//128,(x1-1)//128+1):
    a=block(str(base),ii,jj,kk)
    if a is None:continue
    xx=max(x0,ii*128);xx1=min(x1,ii*128+a.shape[2]);yy=max(y0,jj*128);yy1=min(y1,jj*128+a.shape[1]);zz=max(z0,kk*70);zz1=min(z1,kk*70+a.shape[0])
    if min(xx1-xx,yy1-yy,zz1-zz)<=0:continue
    vol[zz-z0:zz1-z0,yy-y0:yy1-y0,xx-x0:xx1-x0]=a[zz-kk*70:zz1-kk*70,yy-jj*128:yy1-jj*128,xx-ii*128:xx1-ii*128]
 return vol.astype(np.float32),(x0,y0,z0)
records=[]
for gene,(img,spot) in info.items():
 row=meta[(meta.region==1)&(meta.gene==gene)].iloc[0]
 X=sel["cx"]*8;Y=sel["cy"]*8;Z=float(row.source_native_Z_center)
 im,(x0,y0,z0)=roi(img,X,Y,Z)
 # High-pass local 3D, no suppression of weak signal in visualization.
 smooth=gaussian_filter(im,(.8,1.0,1.0))
 base=gaussian_filter(im,(2.0,5,5))
 hp=smooth-base
 usable=hp[2:-2,12:-12,12:-12]
 # Unsupervised local peak candidates in actual raw volume, no overlap with detector list used for segmentation.
 targetspot=[]
 for c in pd.read_csv(spot,header=None,usecols=[0,1,2],dtype=np.float32,chunksize=330000):
  x=c.iloc[:,0].to_numpy()/.23;y=c.iloc[:,1].to_numpy()/.23;z=c.iloc[:,2].to_numpy()/.42
  good=(x>x0+12)&(x<x0+588)&(y>y0+12)&(y<y0+588)&(z>z0+2)&(z<z0+13)
  if good.any():targetspot.append(np.column_stack((z[good]-z0,y[good]-y0,x[good]-x0)))
 target=np.vstack(targetspot) if targetspot else np.empty((0,3))
 ttree=cKDTree(target*np.array([1.83,1,1])) if len(target) else None
 for q in [.995,.9975,.999]:
  thr=float(np.quantile(usable,q))
  peaks=peak_local_max(hp,min_distance=2,threshold_abs=thr,exclude_border=(2,12,12))
  # radius <= 3 native pix in XY, <=2 original z vox equivalently 3.66 pix
  dist=ttree.query(peaks*np.array([1.83,1,1]),k=1)[0] if len(peaks) and ttree else np.array([])
  nnear=int(np.sum(dist<=4.5))
  matches=None
  if q==.9975:
   mraw=im.max(axis=0);imh=hp.max(axis=0)
   fig3,ap=plt.subplots(1,3,figsize=(15.5,5.55))
   lo=float(np.quantile(mraw,.05));hi=float(np.quantile(mraw,.998))
   ap[0].imshow(mraw,cmap="gray",vmin=lo,vmax=hi)
   ap[1].imshow(imh,cmap="gray",vmin=max(0,float(np.quantile(imh,.03))),vmax=float(np.quantile(imh,.999)))
   ap[2].imshow(mraw,cmap="gray",vmin=lo,vmax=hi)
   match=(dist<=4.5)
   if np.any(match):ap[2].scatter(peaks[match,2],peaks[match,1],s=30,edgecolor="#4ef1b6",facecolor="none",linewidths=1.2,alpha=.95)
   if np.any(~match):ap[2].scatter(peaks[~match,2],peaks[~match,1],s=17,marker="x",c="#fa668c",linewidths=1.1,alpha=.93)
   for ax in ap:
    ax.set(aspect="equal",xlim=(0,600),ylim=(600,0))
    ax.set_xticks([]);ax.set_yticks([])
    for spine in ax.spines.values():spine.set_visible(False)
   ap[0].set_title("Native raw 3D MIP · no candidate markers")
   ap[1].set_title("3D source high-pass",fontsize=11)
   ap[2].set_title(f"Raw peaks: {np.sum(match)} matched / {len(peaks)} candidates",fontsize=11)
   fig3.suptitle(f"S500 native raw {gene} · 3D peaks above 99.75th percentile\nPink crosses are unverified candidate maxima, not proven missed molecules",fontsize=12)
   fig3.subplots_adjust(left=.02,right=.98,bottom=.02,top=.84,wspace=.18)
   fig3.savefig(O/f"STAGE1041_REGION1_{gene}_BRIGHT_UNMATCHED_OVERLAY.png",dpi=175,bbox_inches="tight");plt.close(fig3)
  records.append({"gene":gene,"region":1,"candidate_threshold_full_volume_quantile":q,
    "candidate_threshold_filtered_intensity":thr,"n_3D_raw_bright_peaks":len(peaks),"n_raw_peaks_near_RSFISH":nnear,
    "pct_raw_candidates_near_RSFISH":100*nnear/max(1,len(peaks)),
    "official_spots_in_comparable_interior":len(target),
    "other_bright_peaks_not_verified":len(peaks)-nnear})
  print("STAGE1041",gene,q,"peaks",len(peaks),"matched",nnear,"source_spots",len(target),flush=True)
summary=pd.DataFrame(records);summary.to_csv(D/"stage1041_source_image_candidate_peaks_vs_RSFISH.csv",index=False)
fig,axs=plt.subplots(1,3,figsize=(14.5,4.7),layout="constrained")
for ax,g in zip(axs,info):
 t=summary[summary.gene==g].sort_values("candidate_threshold_full_volume_quantile")
 ax.plot(t.candidate_threshold_full_volume_quantile*100,t.pct_raw_candidates_near_RSFISH,marker="o",color="#b75e99")
 for q in t.itertuples(index=False):
  ax.annotate(str(q.n_3D_raw_bright_peaks),(q.candidate_threshold_full_volume_quantile*100,q.pct_raw_candidates_near_RSFISH),fontsize=8,xytext=(0,8),textcoords="offset points",ha="center")
 ax.set(xlabel="Raw 3D high-pass threshold percentile",ylabel="% raw high-pass peaks near RS-FISH centroid",ylim=(0,105),
 title=f"{g} · independent bright-peak review")
 for sp in ("top","right"):ax.spines[sp].set_visible(False)
fig.suptitle("Source image high-pass peaks vs existing spot detection; unmatched maxima are CANDIDATES, not proven missed transcripts",fontsize=11)
fig.savefig(O/"STAGE1041_RAW_IMAGE_PEAKS_VS_RSFISH.png",dpi=195);plt.close(fig)
(D/"stage1041_candidate_recall_authority.json").write_text(json.dumps({"stage":1041,"status":"DONE",
 "source":"Region1 ~138um field, 3 original N5 gene channels S500. Gene-specific Z slab independently optimized in Stage1037.",
 "method":"3D Gaussian raw high-pass bandpass (sigma .8/1 XY ~1px vs 2Z/5XY background). peak_local_max with >=2 original pixels distance; three intensity thresholds 99.5/99.75/99.9% interior filtered voxels. Candidate considered close when distance in native 0.23um XY weighted Z×1.83 is <=4.5 pixels.",
 "interpretation":"Unmatched image high-pass local maxima may be genuine missed RNA, spurious noise, autofluorescence, physical spot double-counting, or other-round background. This screen is not a calibrated true molecular recall, precision or false-positive estimate.",
 "limitations":"One source region, three separate acquisition rounds, different selected z. No uncalled spot ground truth. S500 only."},indent=2))
print("STAGE1041_COMPLETE",flush=True)


