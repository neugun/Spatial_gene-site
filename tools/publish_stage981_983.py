from pathlib import Path
import re,json
p=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
page=p/"single-cell.html";s=page.read_text(encoding="utf8")
assert (p/"assets"/"stage981_vgat_vglut2"/"VGAT_VGLUT2_THRESHOLD_SCAN.png").is_file()
assert (p/"assets"/"stage983_true3d_360_video"/"STAGE983_FINE26_MASKBODY_ZOOM_360_ZOOMOUT_27S.mp4").is_file()
section='''<section id="threshold981">
<h2>VGAT vs VGLUT2: independent threshold and source-data audit</h2>
<p class="lead"><strong>The data do not yet establish two clean E/I marker groups.</strong> Nonzero detection badly overestimates positivity. The Yuhan Wang eLife 2023 >10-spots/cell reference shows Slc17a6 positive in 14.8% of frozen excitatory and 15.7% of frozen inhibitory cells, whereas Slc32a1 is enriched in inhibitory cells. A separate S500→S530 threshold screen is exploratory; S560 is untouched held-out evaluation. Raising both cutoffs to make co-positive cells disappear can turn real E cells into "neither", so we explicitly show coverage and false-negative tradeoffs. S500/S530/S560 are different sections of <em>one animal</em>, not biological replicates.</p>
<div class="cards">
<article class="card"><h3>Threshold parameter sweep</h3><a href="assets/stage981_vgat_vglut2/VGAT_VGLUT2_THRESHOLD_SCAN.png"><img loading="lazy" src="assets/stage981_vgat_vglut2/VGAT_VGLUT2_THRESHOLD_SCAN.png"></a><p>Validation-selected normalized thresholds, not ground truth. <a href="data/stage981_normalized_threshold_sweep.csv">All tested values and objectives</a></p></article>
<article class="card"><h3>Four-class composition across sections</h3><a href="assets/stage981_vgat_vglut2/VGAT_VGLUT2_FOUR_CLASS_SECTIONS.png"><img loading="lazy" src="assets/stage981_vgat_vglut2/VGAT_VGLUT2_FOUR_CLASS_SECTIONS.png"></a><p>VGAT-only, VGLUT2-only, double-positive, and neither, with section-by-class <a href="data/stage981_vgat_vglut2_four_class_by_section.csv">source table</a></p></article>
<article class="card"><h3>Same UMAP, fixed >10 spots vs normalized rule</h3><a href="assets/stage981_vgat_vglut2/VGAT_VGLUT2_FOUR_CLASS_ON_UMAP.png"><img loading="lazy" src="assets/stage981_vgat_vglut2/VGAT_VGLUT2_FOUR_CLASS_ON_UMAP.png"></a><p>Co-positive fractions must be interpreted alongside the number of unassigned cells.</p></article>
</div>
<p class="lead">Source audit: historical H5AD versus current corrected Route A yields near-identical class-specific Slc17a6 positivity and cell-level Pearson correlations around 0.9, so this anomaly is not introduced solely by the latest correction. <a href="data/stage984_h5ad_vs_routeA_source_audit.csv">Inspect the cross-version comparison</a>. A better t-SNE cannot correct inadequate marker specificity.</p>
</section>
'''
if 'id="threshold981"' not in s:
 s=s.replace('<section id="local">',section+'<section id="local">',1)
 assert s.count('id="threshold981"')==1
s=s.replace('<a href="#transmitters">Fine26 + markers</a>','<a href="#transmitters">Fine26 + markers</a><a href="#threshold981">VGAT/VGLUT2 QC</a>')
poster="assets/stage983_true3d_360_video/STAGE983_SEVEN_FRAME_QC.jpg"
movie="assets/stage983_true3d_360_video/STAGE983_FINE26_MASKBODY_ZOOM_360_ZOOMOUT_27S.mp4"
replacement='''<section id="video980">
<h2>Fine26 true-maskbody full 360° orbit · zoom in and out</h2>
<p class="lead">New Stage983: <strong>27 seconds, 25 fps, 675 frames, continuous 360° camera orbit</strong>. Start with complete coronal view, smoothly zoom into the 3D cell bodies, rotate all the way around, zoom back out and hold. S500, S530 and S560 are shown simultaneously. The Stage631 physical X/Y convention and Fine26 26-color key are unchanged; 3D bodies are sampled points from actual downseg masks rather than cell centers.</p>
<video controls playsinline preload="metadata" poster="'''+poster+'''" style="display:block;width:100%;max-height:640px;background:#000;border-radius:8px">
<source src="'''+movie+'''" type="video/mp4">Your browser cannot play the video.</video>
<p class="lead">Inspect <a href="'''+poster+'''">seven camera-angle quality-control frames</a>, or <a href="assets/stage980_true3d_video/STAGE980_FINE26_MASKBODY_ZOOM_ORBIT_12S.mp4">compare against the older 12-second video</a>. Full 360° refers to the physical camera orbit, not rotation of the fixed source images.</p>
</section>'''
s,n=re.subn(r'<section id="video980">.*?</section>',lambda m:replacement,s,count=1,flags=re.S)
assert n==1
page.write_text(s,encoding="utf8")
print("STAGE981_983_PAGE_UPDATED",page.stat().st_size,flush=True)
