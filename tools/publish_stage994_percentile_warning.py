from pathlib import Path
import re,json
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
f=P/"single-cell.html";h=f.read_text(encoding="utf8")
add="""<article class="card"><h3>Shared P70: exploratory candidate on full-data t-SNE</h3><a href="assets/stage988_gene_percentile/STAGE988_P70_SHARED_PERCENTILE_TSNE.png"><img loading="lazy" src="assets/stage988_gene_percentile/STAGE988_P70_SHARED_PERCENTILE_TSNE.png" alt="VGAT and VGLUT2 equal P70 thresholds with both-positive counts"/></a><p>Raw-count targets 15/5 are closest to equal percentiles around P69.5; practical P70 has VGAT norm=672.1 and VGLUT2 norm=396.6, 9.1% double-positive and 49.1% neither among all 71,950 cells. This is not evidence for independently validated cell classes. <a href="data/stage988_shared_p70_classification_by_section.csv">Every section and broad class</a>.</p></article>
<article class="card"><h3>Crucial control: normalized class separation is driven by unequal 27-gene totals</h3><a href="assets/stage988_gene_percentile/STAGE988_P70_DEPTH_DENOMINATOR_CONTROL.png"><img loading="lazy" src="assets/stage988_gene_percentile/STAGE988_P70_DEPTH_DENOMINATOR_CONTROL.png" alt="VGLUT2 E and I positive fractions after raw and 27 gene normalization"/></a><p>At raw-count P70, E and I VGLUT2+ are BOTH 30.05%. After dividing by 27 panel counts, E=60.94% and I=7.57% because E median panel total=48.3 versus I=312.1. Earlier H5AD shows the same confounding. Thus panel-total normalization must not be marketed as proven biological specificity. <a href="data/stage988_h5ad_routeA_depth_normalization_sensitivity.csv">Source-version comparisons</a>.</p></article>
"""
assert h.count('id="stage988-percentile"')==1
pattern=re.compile(r'(<section id="stage988-percentile">.*?<div class="cards">)(.*?)(</div>)',flags=re.S)
h,n=pattern.subn(lambda m:m.group(1)+m.group(2)+add+m.group(3),h,count=1)
assert n==1
for name in ['STAGE988_P70_SHARED_PERCENTILE_TSNE.png','STAGE988_P70_DEPTH_DENOMINATOR_CONTROL.png']:
 assert (P/'assets'/'stage988_gene_percentile'/name).is_file(),name
f.write_text(h,encoding="utf8")
print("STAGE994_DEPTH_CONTROL_PUBLISHED",len(h),flush=True)
