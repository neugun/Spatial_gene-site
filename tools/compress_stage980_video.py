from pathlib import Path
import subprocess,re,shutil
p=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review\assets\stage980_true3d_video")
src=p/"STAGE980_FINE26_MASKBODY_ZOOM_ORBIT_12S.mp4"
cand=p/"compressed_try.mp4"
cmd=["ffmpeg","-y","-loglevel","error","-i",str(src),"-an","-c:v","libx264","-preset","medium","-crf","31","-pix_fmt","yuv420p","-movflags","+faststart",str(cand)]
q=subprocess.run(cmd,capture_output=True,text=True)
assert q.returncode==0,q.stderr
a=subprocess.run(["ffmpeg","-hide_banner","-i",str(src),"-i",str(cand),"-lavfi","ssim","-f","null","NUL"],capture_output=True,text=True)
matches=re.findall(r"All:\s*([0-9.]+)",a.stderr)
score=float(matches[-1]) if matches else -1
print("COMPRESS_QC","orig",src.stat().st_size,"cand",cand.stat().st_size,"SSIM",score,flush=True)
if score>=.93 and cand.stat().st_size<src.stat().st_size*.65:
 archive=Path(r"G:\Map6_recover_all\stage980_movie_HQ");archive.mkdir(exist_ok=True)
 shutil.copy2(src,archive/src.name)
 cand.replace(src)
 print("MOBILE_MP4_PROMOTED",src.stat().st_size,flush=True)
else:
 cand.unlink()
 print("KEEP_HIGH_QUALITY",flush=True)
