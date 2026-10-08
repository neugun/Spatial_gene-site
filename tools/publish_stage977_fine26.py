"""Publish Stage977 biologically qualified Fine26/marker UMAP figures."""
from pathlib import Path
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
page=P/"single-cell.html"
s=page.read_text(encoding="utf8")
anchor='<section id="local">'
if 'id="transmitters"' not in s:
    fragment='''<section id="transmitters">
<h2>Fine26 UMAP geometry and neurotransmitter-marker specificity</h2>
<p class="lead">All panels use the same Stage704 molecular embedding and frozen Fine26 identities (71,950 cells). The updated 26-color categorical palette separates adjacent Fine26 populations without redefining clusters. Raw detected spots are broad across cells; marker detection alone is <b>not</b> cell-type assignment. The top-5%-expression panels correct for each cell's total 27-gene count and highlight relative enrichment.</p>
<div class="cards">
<article class="card"><h3>Fine26 · distinct categorical colors</h3><a href="assets/stage977_fine26_umap/FINE26_CATEGORICAL_GLASBEY.png"><img loading="lazy" src="assets/stage977_fine26_umap/FINE26_CATEGORICAL_GLASBEY.png" alt="Fine26 molecular UMAP with 26 distinguishable colors"></a><p>The historical 26 labels are unchanged. New color key is available as <a href="data/stage977_fine26_color_key.csv">CSV</a>.</p></article>
<article class="card"><h3>VGAT / VGLUT2 / NE / ChAT · observed counts</h3><a href="assets/stage977_fine26_umap/VGAT_VGLUT2_NE_CHAT_MARKER_UMAP.png"><img loading="lazy" src="assets/stage977_fine26_umap/VGAT_VGLUT2_NE_CHAT_MARKER_UMAP.png"></a><p>Direct gene detections, with their own numeric log-count colorbars. High detection frequency must not be interpreted as a high fraction of corresponding neuronal types.</p></article>
<article class="card"><h3>Top 5% normalized marker expression</h3><a href="assets/stage977_fine26_umap/VGAT_VGLUT2_NE_CHAT_TOP5PCT.png"><img loading="lazy" src="assets/stage977_fine26_umap/VGAT_VGLUT2_NE_CHAT_TOP5PCT.png"></a><p>Gene-specific 95th-percentile thresholds on cell-depth-normalized counts; same UMAP axes and cell universe.</p></article>
<article class="card"><h3>Marker enrichment × frozen broad class</h3><a href="assets/stage977_fine26_umap/TRANSMITTER_TOP5_BROAD_ENRICHMENT.png"><img loading="lazy" src="assets/stage977_fine26_umap/TRANSMITTER_TOP5_BROAD_ENRICHMENT.png"></a><p>Enrichment factors are descriptive across cells, not animal-replicate inference. <a href="data/stage977_marker_top5_broad_enrichment.csv">Underlying class-by-marker table</a>.</p></article>
<article class="card"><h3>Broad-class control on the UMAP</h3><a href="assets/stage977_fine26_umap/BROAD_CLASSES_UMAP.png"><img loading="lazy" src="assets/stage977_fine26_umap/BROAD_CLASSES_UMAP.png"></a></article>
</div>
</section>
'''
    assert anchor in s
    s=s.replace(anchor,fragment+anchor,1)
    s=s.replace('<a href="#taxonomy">Molecular state</a>','<a href="#taxonomy">Molecular state</a><a href="#transmitters">Fine26 + markers</a>',1)
# Prior results maintained as historical comparison, new release entries directly visible.
for name in ["FINE26_CATEGORICAL_GLASBEY.png","VGAT_VGLUT2_NE_CHAT_MARKER_UMAP.png","VGAT_VGLUT2_NE_CHAT_TOP5PCT.png","TRANSMITTER_TOP5_BROAD_ENRICHMENT.png","BROAD_CLASSES_UMAP.png"]:
    assert (P/"assets"/"stage977_fine26_umap"/name).is_file(),name
page.write_text(s,encoding="utf8")
print("STAGE977_SINGLE_CELL_PAGE_UPDATED",s.count('id="transmitters"'))
