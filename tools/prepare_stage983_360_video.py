"""Stage983 film: full 360° orbit + zoom-in and zoom-out, Stage631 orientation unchanged."""
from pathlib import Path
src=Path(r"G:\Spatial_gene_site_publish\tools\build_stage980_true3d_video.py")
dst=src.with_name("build_stage983_true3d_360_video.py")
s=src.read_text(encoding="utf8")
replaces=[
('stage980_true3d_video','stage983_true3d_360_video'),
('stage980_movie_QC','stage983_movie_QC'),
('window_size=(1800,650)','window_size=(1650,600)'),
('t=k/25;zoom=smooth((t-2)/3);turn=smooth((t-5)/7)','t=k/25;zoom=smooth((t-2.5)/3.5)*(1-smooth((t-22.5)/3.5));turn=smooth((t-5.8)/16.8)'),
('theta=math.radians(68*turn)','theta=2*math.pi*turn'),
('for k in (0,50,100,175,299):','for k in (0,100,205,310,415,535,674):'),
('(1800,975)','(1650,1200)'),
('small.thumbnail((900,325))','small.thumbnail((825,300))'),
('((i%2)*900,(i//2)*325)','((i%2)*825,(i//2)*300)'),
('STAGE980_FIVE_FRAME_QC.jpg','STAGE983_SEVEN_FRAME_QC.jpg'),
('(1800,650))','(1650,600))'),
('for k in range(300):','for k in range(675):'),
('STAGE980_FINE26_MASKBODY_ZOOM_ORBIT_12S.mp4','STAGE983_FINE26_MASKBODY_ZOOM_360_ZOOMOUT_27S.mp4'),
('"-crf","24"','"-crf","28"'),
('stage980_video_current_authority.json','stage983_video_current_authority.json'),
('"stage":980,"duration_seconds":12,"fps":25,"frames":300,"dimensions":[1800,650]','"stage":983,"duration_seconds":27,"fps":25,"frames":675,"dimensions":[1650,600]'),
('"file":out.name}', '"file":out.name,"camera_path":"hold full view -> ease zoom in -> complete 360 degree orbit -> ease zoom out -> hold full view"}')
]
for old,new in replaces:
 assert old in s,old
 s=s.replace(old,new)
dst.write_text(s,encoding="utf8")
print("STAGE983_SCRIPT_READY",dst,flush=True)
