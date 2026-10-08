"""Stage980 3D Fine26 camera movie, Stage191/214/215 12s/25fps specification."""
from pathlib import Path
import numpy as np,pandas as pd,sys,json,math,time,cv2,subprocess
import matplotlib.colors as mc
import pyvista as pv
from PIL import Image
ROOT=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
MC=Path(r"G:\PeriLC_current\D_mirror\11_PPTsummary\PeriLC_CHATGPT_RESULTS_20260914\141_DOWNSEG_MASKCLOUD")
PUB=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
O=PUB/"assets"/"stage980_true3d_video";O.mkdir(parents=True,exist_ok=True)
QA=Path(r"G:\Map6_recover_all\stage980_movie_QC");QA.mkdir(parents=True,exist_ok=True)
K=pd.read_csv(PUB/"data"/"stage977_fine26_color_key.csv")
A=pd.read_csv(ROOT/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["section","roi_id","fine26","fine26_stable_id","x_um","y_um"])
mp=A[["fine26","fine26_stable_id"]].drop_duplicates().set_index("fine26").fine26_stable_id.to_dict()
colors=(255*np.array([mc.to_rgb(dict(zip(K.fine26,K.color))[mp[i]]) for i in range(26)])).astype(np.uint8)
P=pv.Plotter(shape=(1,3),off_screen=True,window_size=(1800,650),border=False)
cams=[];rg=np.random.default_rng(980)
for j,s in enumerate((500,530,560)):
 P.subplot(0,j);B=A[A.section==s];r=B.roi_id.to_numpy(int)
 Q=np.load(MC/f"maskcloud_s{s}_xy5_z7_cap45.npz");rid=Q["roi_id"].astype(int)
 m=max(rid.max(),r.max())+1;ax=np.full(m,np.nan);ay=np.full(m,np.nan);li=np.full(m,-1)
 ax[r]=B.x_um;ay[r]=B.y_um;li[r]=B.fine26
 x=Q["x_um"]/2;y=Q["y_flipped_um"]/2;z=Q["z_um"]/2
 ok=np.flatnonzero(np.isfinite(ax[rid])&np.isfinite(ay[rid])&(li[rid]>=0))
 ox=np.median(ax[rid[ok]]+x[ok]);oy=np.median(ay[rid[ok]]-y[ok])
 if len(ok)>70000:ok=np.sort(rg.choice(ok,70000,replace=False))
 xyz=np.stack([ox-x[ok],oy+y[ok],z[ok]],axis=1).astype(np.float32)
 P.set_background("black")
 P.add_points(xyz,scalars=colors[li[rid[ok]].astype(int)],rgb=True,point_size=2.15,opacity=.83,render_points_as_spheres=False)
 P.add_text(f"S{s}  Fine26 true segmentation bodies",position="upper_left",font_size=11,color="white")
 center=np.median(xyz,axis=0)
 span=max(np.ptp(np.quantile(xyz[:,0],[.003,.997])),np.ptp(np.quantile(xyz[:,1],[.003,.997])))
 cam=P.camera;cam.parallel_projection=True;cam.parallel_scale=span*.66
 cam.position=(center[0],center[1],center[2]+span*3.2);cam.focal_point=tuple(center);cam.up=(0,1,0)
 cams.append((cam,center,span));print("SECTION_SETUP",s,len(ok),round(span),flush=True)
def smooth(v):
 v=min(max(v,0),1);return v*v*(3-2*v)
def makeframe(k):
 t=k/25;zoom=smooth((t-2)/3);turn=smooth((t-5)/7)
 for cam,c,scale in cams:
  theta=math.radians(68*turn);d=scale*3.2
  cam.position=(c[0]+d*math.sin(theta),c[1]+d*.08*turn,c[2]+d*math.cos(theta))
  cam.focal_point=tuple(c);cam.up=(0,1,0);cam.parallel_scale=scale*.66*(1-.35*zoom)
 P.render();return P.screenshot(return_img=True)
ts=time.time()
if "--full" not in sys.argv:
 frames=[]
 for k in (0,50,100,175,299):
  img=makeframe(k);name=QA/f"QC_FRAME_{k:03d}.png";Image.fromarray(img).save(name);frames.append(img)
  print("FRAME_QC",k,img.shape,round(img.mean(),2),flush=True)
 contact=Image.new("RGB",(1800,975),"black")
 for i,f in enumerate(frames):
  small=Image.fromarray(f);small.thumbnail((900,325));contact.paste(small,((i%2)*900,(i//2)*325))
 contact.save(O/"STAGE980_FIVE_FRAME_QC.jpg",quality=86)
 print("PREVIEW_COMPLETE",round(time.time()-ts,1),flush=True)
else:
 temp=O/"intermediate_mp4v.mp4";w=cv2.VideoWriter(str(temp),cv2.VideoWriter_fourcc(*"mp4v"),25,(1800,650))
 assert w.isOpened()
 for k in range(300):
  w.write(cv2.cvtColor(makeframe(k),cv2.COLOR_RGB2BGR))
  if k%50==0:print("RENDER_FRAME",k,round(time.time()-ts,1),flush=True)
 w.release();out=O/"STAGE980_FINE26_MASKBODY_ZOOM_ORBIT_12S.mp4"
 res=subprocess.run(["ffmpeg","-y","-loglevel","error","-i",str(temp),"-c:v","libx264","-pix_fmt","yuv420p","-crf","24","-movflags","+faststart",str(out)],capture_output=True,text=True)
 assert res.returncode==0,res.stderr;temp.unlink()
 (PUB/"data"/"stage980_video_current_authority.json").write_text(json.dumps({"stage":980,"duration_seconds":12,"fps":25,"frames":300,"dimensions":[1800,650],"source":"Stage141 true segmentation maskcloud","frame":"Stage631 preflipped physical XY, verified ROI alignment","palette":"Stage977 26 distinct Fine26 colors","reference":["Stage191","Stage214","Stage215"],"file":out.name},indent=2))
 print("VIDEO_COMPLETE",str(out),out.stat().st_size,round(time.time()-ts,1),flush=True)
P.close()
