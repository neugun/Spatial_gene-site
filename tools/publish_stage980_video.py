"""Publish Stage980 current 3D anatomical video without changing orientation."""
from pathlib import Path
import json
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
mp4=P/"assets"/"stage980_true3d_video"/"STAGE980_FINE26_MASKBODY_ZOOM_ORBIT_12S.mp4"
poster=P/"assets"/"stage980_true3d_video"/"STAGE980_FIVE_FRAME_QC.jpg"
assert mp4.is_file() and mp4.stat().st_size>1000000
assert poster.is_file()
page=P/"single-cell.html"
s=page.read_text(encoding="utf8")
if 'id="video980"' not in s:
    video='''<section id="video980">
<h2>Fine26 true mask-body: coronal context → smooth zoom → 3D orbit</h2>
<p class="lead">12 s · 25 frames/s · S500/S530/S560 shown simultaneously. This video renders sampled points from <b>real downseg segmentation bodies</b>, not cell centers. Stage631 physical X and Y are already flipped and were verified against ROI coordinates; there is no second flip. The same 26-color Fine26 palette appears in the UMAP, spatial panels and true-body movie.</p>
<video controls playsinline preload="metadata" poster="assets/stage980_true3d_video/STAGE980_FIVE_FRAME_QC.jpg" style="display:block;width:100%;max-height:640px;background:#000;border-radius:8px">
<source src="assets/stage980_true3d_video/STAGE980_FINE26_MASKBODY_ZOOM_ORBIT_12S.mp4" type="video/mp4">
Your browser cannot play the video.</video>
<p class="lead">Camera QC: <a href="assets/stage980_true3d_video/STAGE980_FIVE_FRAME_QC.jpg">five reference frames</a> · historical choreography from Stages 191, 214 and 215. Model and cell identity are unchanged.</p>
</section>
'''
    s=s.replace('<section id="integration">',video+'<section id="integration">',1)
    s=s.replace('<a href="#maskbody">True maskbody</a>','<a href="#maskbody">True maskbody</a><a href="#video980">3D video</a>',1)
page.write_text(s,encoding="utf8")
print("STAGE980_VIDEO_PAGE_PUBLISHED",mp4.stat().st_size,flush=True)
