"""Stage1038: native 3D fluorescence spot sharpness on 120 source detections/gene.
Empirical XY and Z full-width at half maximum (FWHM), not microscope PSF fit.
Detect only in-chunk puncta to avoid huge image fetches, no inference across animals.
"""
from pathlib import Path
import gzip,struct,functools,json
import numpy as np,pandas as pd
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r"Z:\sternsonlab\Zhenggang\2acq")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");O=P/"assets"/"stage1038_threegene_focusQC";O.mkdir(exist_ok=True);D=P/"data"
src={
"Snap25":(R/"outputs/M5L_500_4channels_5X5tile_1/stitching/export.n5/c1/s0",R/"map6_allsections_slurm_20260926/baseline/S500/R1_Snap25_c1/spots_c1.txt",5307181),
"VGAT":(R/"outputs/M5L_500_4channels_5X5tile_R77/stitching/export.n5/c0/s0",R/"map6_allsections_slurm_20260926/baseline/S500/R7_Slc32a1_c0/spots_c0.txt",541321),
"VGLUT2":(R/"outputs/M5L_500_4channels_5X5tile_R99/stitching/export.n5/c3/s0",R/"map6_allsections_slurm_20260926/corrected_spots/S500/R9_Slc17a6_c3/spots_c3.txt",2981450)}
@functools.lru_cache(maxsize=80)
def chunk(p,x,y,z):
 f=Path(p)/str(x)/str(y)/str(z)
 if not f.is_file():return None
 b=f.read_bytes();mode,nd=struct.unpack(">HH",b[:4]);dims=struct.unpack(">"+str(nd)+"I",b[4:4+4*nd])
 return np.frombuffer(gzip.decompress(b[4+4*nd:]),dtype=">u2").reshape((dims[2],dims[1],dims[0]))
def w50(a,center):
 p=np.asarray(a,dtype=float)
 maximum=max(p[center],0.)
 if maximum<=0:return np.nan
 threshold=maximum*.5
 lo=center;hi=center
 while lo>0 and p[lo-1]>=threshold:lo-=1
 while hi<len(p)-1 and p[hi+1]>=threshold:hi+=1
 if lo==0 or hi==len(p)-1:return np.nan
 return float(hi-lo+1)
rng=np.random.default_rng(1038);rows=[]
for gene,(path,spfile,N) in src.items():
 reservoir=[]
 for c in pd.read_csv(spfile,header=None,usecols=[0,1,2],dtype=np.float32,chunksize=275000):
  k=max(1,int(np.ceil(len(c)*400/N)))
  j=rng.choice(len(c),size=min(len(c),k),replace=False)
  reservoir.append(c.iloc[j].to_numpy())
 samples=np.vstack(reservoir);rng.shuffle(samples)
 n=0
 for xyz in samples:
  x=int(round(xyz[0]/.23));y=int(round(xyz[1]/.23));z=int(round(xyz[2]/.42))
  ix=x%128;iy=y%128;iz=z%70
  if not(20<ix<106 and 20<iy<106 and 14<iz<55):continue
  blk=chunk(str(path),x//128,y//128,z//70)
  if blk is None:continue
  vol=blk[iz-8:iz+9,iy-12:iy+13,ix-12:ix+13].astype(np.float32)
  yy,xx=np.indices((25,25));ring=(xx-12)**2+(yy-12)**2>=64
  bg=float(np.median(vol[8,ring]))
  xy=np.max(vol[7:10],axis=0)-bg
  p_y,p_x=np.unravel_index(np.argmax(xy[10:15,10:15]),(5,5));p_y+=10;p_x+=10
  peak=float(xy[p_y,p_x])
  sig=1.4826*float(np.median(np.abs(vol[8,ring]-bg)))+1
  if peak<4*sig:continue
  profx=np.maximum(0,xy[p_y,:]);profy=np.maximum(0,xy[:,p_x])
  profz=np.maximum(0,np.max(vol[:,p_y-1:p_y+2,p_x-1:p_x+2],axis=(1,2))-bg)
  xw=w50(profx,p_x);yw=w50(profy,p_y);zw=w50(profz,int(np.argmax(profz[6:11]))+6)
  if not np.isfinite(xw) or not np.isfinite(yw) or not np.isfinite(zw):continue
  rows.append(dict(gene=gene,xy_fwhm_um=float(.23*np.sqrt(xw*yw)),z_fwhm_um=float(.42*zw),
   xy_FWHM_native_px=float(np.sqrt(xw*yw)),snr=peak/sig))
  n+=1
  if n>=120:break
 print("FOCUS_SAMPLE",gene,n,flush=True)
X=pd.DataFrame(rows)
S=X.groupby("gene").agg(n=("xy_fwhm_um","size"),
 xy_FWHM_median_um=("xy_fwhm_um","median"),z_FWHM_median_um=("z_fwhm_um","median"),
 xy_FWHM_p90_um=("xy_fwhm_um",lambda q:float(q.quantile(.9))),
 frac_XY_width_gt_1um=("xy_fwhm_um",lambda q:float((q>1).mean())),
 snr_median=("snr","median")).reset_index()
S.to_csv(D/"stage1038_n5_native_spot_sharpness_fwhm_summary.csv",index=False)
fig,axes=plt.subplots(1,2,figsize=(11.5,4.8),layout="constrained")
colors=["#22a897","#e3a038","#c45797"]
for ax,field,name in [(axes[0],"xy_fwhm_um","XY punctum FWHM (µm)"),(axes[1],"z_fwhm_um","Axial Z FWHM (µm)")]:
 vals=[X.loc[X.gene==g,field] for g in ("Snap25","VGAT","VGLUT2")]
 plots=ax.boxplot(vals,patch_artist=True,showfliers=False,medianprops=dict(color="#252525"))
 for patch,col in zip(plots["boxes"],colors):patch.set_facecolor(col);patch.set_alpha(.7)
 ax.set_xticks([1,2,3],["Snap25","VGAT","VGLUT2"])
 ax.set(ylabel=name,title="Native 3D source puncta (random sample)")
 for sp in ("top","right"):ax.spines[sp].set_visible(False)
fig.suptitle("Comparative intrinsic spot sharpness; source rounds differ and only detectable high-SNR spots included",fontsize=11.5)
fig.savefig(O/"STAGE1038_THREE_GENE_NATIVE_XY_Z_SHARPNESS.png",dpi=185);plt.close(fig)
(D/"stage1038_native_focusQC_authority.json").write_text(json.dumps({"stage":1038,"status":"DONE","quality_sampling":"random full-S500 original RSFISH catalog rows, retaining valid N5 source chunks and peak SNR>4; at most 120/gene",
"estimated_FWHM":"3D native fluorescence 25x25x17 voxel around detected centroid; ring-median local baseline; halfmax width contiguous in native XY and Z. Not a PSF model or calibrated optical-resolution estimate.",
"limits":"Different imaging rounds, z-dependent blur; brightness selected by detector, not universal optical focus across all attempted spots. Pixel XY .23µm and Z .42µm. One animal/section."},indent=2),encoding="utf8")
print("STAGE1038_DONE",S.to_string(index=False),flush=True)
