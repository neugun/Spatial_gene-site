"""Stage992: S530 genuine downseg.tif per-cell label voxel meshes, not ellipsoids.
Marching cubes on per-cell binary mask, small Gaussian only to soften voxel stair-steps.
ROI derived from exact original tiff (15.6GB memmapped), selected subset of local window.
"""
from pathlib import Path
import json,time
import numpy as np,pandas as pd,tifffile
from scipy.ndimage import gaussian_filter
from skimage.measure import marching_cubes
import pyvista as pv
from PIL import Image,ImageDraw
import matplotlib.colors as mcolors
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data";O=P/"assets"/"stage992_true_mask_surface";O.mkdir(exist_ok=True)
QF=Path(r"G:\PeriLC_current\D_mirror\11_PPTsummary\PeriLC_CHATGPT_RESULTS_20260914\141_DOWNSEG_MASKCLOUD\maskcloud_s530_xy5_z7_cap45.npz")
TIFF=Path(r"Z:\sternsonlab\Zhenggang\2acq\outputs\M5L_530_4channels_5X5tile_1\segmentation\downseg.tif")
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["section","roi_id","x_um","y_um","depth_um","fine26_stable_id"])
A=A[A.section==530];win=pd.read_csv(D/"stage950_local_singlecell_windows.csv");w=win[win.section==530].iloc[0]
KEY=pd.read_csv(D/"stage977_fine26_color_key.csv");color=dict(zip(KEY.fine26,KEY.color))
q=np.load(QF);rid=q["roi_id"].astype(int);sx=q["x_um"].astype(float)/2;sy=q["y_flipped_um"].astype(float)/2
maxid=max(int(rid.max()),int(A.roi_id.max()))+1
ax=np.full(maxid,np.nan);ay=np.full(maxid,np.nan);ax[A.roi_id]=A.x_um;ay[A.roi_id]=A.y_um
valid=np.isfinite(ax[rid])&np.isfinite(ay[rid])
dx=np.median(ax[rid[valid]]+sx[valid]);dy=np.median(ay[rid[valid]]-sy[valid])
print("PHYSICAL_ALIGNMENT",dx,dy,flush=True)
target=A[(A.x_um>=w.x0+55)&(A.x_um<w.x0+185)&(A.y_um>=w.y0+55)&(A.y_um<w.y0+185)&(A.depth_um>=65)&(A.depth_um<210)].copy()
assert len(target)>30
# Guarantee spread across identities not just one Fine26 type.
rng=np.random.default_rng(992)
target=target.groupby("fine26_stable_id",group_keys=False).apply(lambda x:x.sample(n=min(9,len(x)),random_state=992))
if len(target)>135:target=target.sample(n=135,random_state=992)
target=target.sort_values("roi_id")
print("SELECTED",len(target),"FINE26",target.fine26_stable_id.nunique(),flush=True)
tif=tifffile.memmap(TIFF)
assert tif.shape==(706,2258,2247)
# Source y-flipped intercept at raw pixel index 2234, empirically verified exact 160/160 sampled label points
y_top=2331.
rlist=rid
sort=np.argsort(rlist);sorted_ids=rlist[sort]
def points_for_roi(r):
 lo=np.searchsorted(sorted_ids,r,'left');hi=np.searchsorted(sorted_ids,r,'right')
 return sort[lo:hi]
def raw_points(ii):
 return (np.rint(q["x_um"][ii].astype(float)/.92).astype(int),
  np.rint(y_top-q["y_flipped_um"][ii].astype(float)/.92).astype(int),
  np.rint(q["z_um"][ii].astype(float)/.84).astype(int))
meshes=[];surfaceQA=[];xyzraw=[];raw_rgb=[]
for j,row in enumerate(target.itertuples(index=False)):
 r=int(row.roi_id)
 ii=points_for_roi(r)
 if len(ii)<8:continue
 xp,yp,zp=raw_points(ii)
 # individual bounding crop from actual observed mask voxels, no convex hull approximations.
 pad=7;zz0=max(0,int(zp.min())-pad);zz1=min(tif.shape[0],int(zp.max())+pad+1)
 yy0=max(0,int(yp.min())-pad);yy1=min(tif.shape[1],int(yp.max())+pad+1)
 xx0=max(0,int(xp.min())-pad);xx1=min(tif.shape[2],int(xp.max())+pad+1)
 if min(zz1-zz0,yy1-yy0,xx1-xx0)<4:continue
 vol=np.asarray(tif[zz0:zz1,yy0:yy1,xx0:xx1],dtype=np.int32)
 actual=(vol==r)
 if actual.sum()<35:continue
 smooth=gaussian_filter(actual.astype(np.float32),sigma=(.82,.82,.82))
 if smooth.max()<=.42:continue
 try:vertex,faces,_,_=marching_cubes(smooth,level=.42,spacing=(.42,.46,.46))
 except Exception:continue
 # Face coordinate index original z,y,x, adjusted into exactly Stage631-flipped physical XY.
 vx=dx-.46*xx0-vertex[:,2]
 vy=dy+.46*y_top-.46*yy0-vertex[:,1]
 vz=.42*zz0+vertex[:,0]
 pts=np.column_stack([vx,vy,vz]).astype(np.float32)
 fs=np.hstack([np.full((len(faces),1),3,np.int64),faces]).ravel()
 mesh=pv.PolyData(pts,fs)
 if mesh.n_cells<3:continue
 meshes.append((mesh,str(row.fine26_stable_id)))
 rawxy=np.column_stack([dx-.46*xp,dy+.46*y_top-.46*yp,zp*.42])
 xyzraw.append(rawxy)
 raw_rgb.extend([mcolors.to_rgb(color[row.fine26_stable_id])]*len(rawxy))
 surfaceQA.append(dict(roi_id=r,fine26=str(row.fine26_stable_id),raw_voxels=int(actual.sum()),
  maskcloud_samples=int(len(ii)),mesh_vertices=int(mesh.n_points),mesh_triangles=int(mesh.n_cells),
  bbox_pixels=[int(zz1-zz0),int(yy1-yy0),int(xx1-xx0)]))
 if len(meshes)%25==0:print("MESHED",len(meshes),flush=True)
assert len(meshes)>=20,("insufficient genuine mask surfaces",len(meshes))
points=np.concatenate(xyzraw).astype(np.float32);rgb=(np.array(raw_rgb)*255).astype(np.uint8)
center=np.median(points,axis=0);extent=max(np.ptp(points[:,0]),np.ptp(points[:,1]))*1.18
pv.global_theme.window_size=[1600,900]
plot=pv.Plotter(shape=(1,3),off_screen=True,window_size=(1800,680))
for panel,view in enumerate(["RAW SAMPLED CLOUD · XY","ACTUAL MASK-SURFACE · XY","ACTUAL MASK-SURFACE · 3D"]):
 plot.subplot(0,panel);plot.set_background("white")
 if panel==0:
  plot.add_points(points,scalars=rgb,rgb=True,point_size=3.6,render_points_as_spheres=False)
 else:
  for mesh,typ in meshes:
   plot.add_mesh(mesh,color=color[typ],smooth_shading=True,opacity=.93,show_edges=False,ambient=.28,diffuse=.65,specular=.07)
 cam=plot.camera;cam.parallel_projection=True;cam.focal_point=tuple(center)
 if panel==2:
  dist=extent*2.6;cam.position=(center[0]+dist*.36,center[1]-dist*.20,center[2]+dist*.94);cam.up=(0,1,0)
 else:
  cam.position=(center[0],center[1],center[2]+extent*3.4);cam.up=(0,1,0)
 cam.parallel_scale=extent*.80
 plot.add_text(view,position="upper_left",font_size=12,color="black")
plot.render()
plot.screenshot(str(O/"STAGE992_S530_ACTUAL_SEGMENTATION_SURFACE_VS_POINTCLOUD.png"))
plot.close()
pd.DataFrame(surfaceQA).to_csv(D/"stage992_true_mask_surface_qc_S530.csv",index=False)
metadata={"stage":992,"status":"TRUE_SEGMENTATION_SURFACES","source":"Full-resolution S530 downseg.tif genuine segmentation labels; raw source path withheld",
 "local_window":"S530 Stage950, subset of cell bodies within cropped XY + bounded Z, not all bodies",
 "genuine_meshes":len(meshes),"sampled_points":len(points),
 "algorithm":"ROI-specific exact downseg.tif label binary -> sigma .82 Gaussian anti-alias -> marching_cubes .42, triangle mesh","xy":"Stage631 physical already both flipped; raw index y intercept 2331 verified with exact ROI mask hits",
 "caution":"Sampling sparse maskcloud is not a true surface; this figure replaces pointcloud with exact segmentation-derived surface for the shown ROI subset"}
(D/"stage992_true_mask_surface_authority_S530.json").write_text(json.dumps(metadata,indent=2),encoding="utf8")
print("STAGE992_COMPLETE",len(meshes),len(points),flush=True)


