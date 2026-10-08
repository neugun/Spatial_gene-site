from pathlib import Path
import subprocess,re,shutil
p=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review\assets\stage983_true3d_360_video")
src=p/"STAGE983_FINE26_MASKBODY_ZOOM_360_ZOOMOUT_27S.mp4"
assert src.is_file()
best=None
for crf in (32,30):
 c=p/f"_crf{crf}_test.mp4"
 q=subprocess.run(["ffmpeg","-y","-loglevel","error","-i",str(src),"-an","-c:v","libx264","-preset","medium","-crf",str(crf),"-pix_fmt","yuv420p","-movflags","+faststart",str(c)],capture_output=True,text=True)
 assert q.returncode==0,q.stderr
 a=subprocess.run(["ffmpeg","-hide_banner","-i",str(src),"-i",str(c),"-lavfi","ssim","-f","null","NUL"],capture_output=True,text=True)
 ms=re.findall(r"All:\s*([0-9.]+)",a.stderr)
 score=float(ms[-1]) if ms else -1
 print("VIDEO_CRF",crf,"bytes",c.stat().st_size,"ssim",score,flush=True)
 if score>=.94 and c.stat().st_size<65_000_000:
  best=c;break
 if best is None and score>=.94:best=c
 if crf==32 and best==c:break
 if best!=c:c.unlink()
if best:
 backup=Path(r"G:\Map6_recover_all\stage983_movie_HQ");backup.mkdir(exist_ok=True)
 shutil.copy2(src,backup/src.name)
 best.replace(src)
 print("VIDEO_COMPRESSED",src.stat().st_size,flush=True)
else:print("KEEP_INITIAL_MASTER",src.stat().st_size,flush=True)
