from pathlib import Path
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
f=P/"single-cell.html";h=f.read_text(encoding="utf8")
base="assets/stage982_tsne_comparison/"
hist="assets/stage985_historical_tsne/"
for path in [base+"multiscale_30_300_COMPARISON.png",base+"STAGE982_PAPER_GT10_ALLCELL_FOURCLASS_UMAP_TSNE.png",hist+"HISTORICAL_TSNE_VS_CURRENT_EMBEDDINGS_MATCHED.png"]:
 assert (P/path).is_file(),path
block='''<section id="tsne982">
<h2>Full-data t-SNE versus UMAP · matched historical reference</h2>
<p class="lead">A real multiscale t-SNE was fitted to <b>all 71,950 cells</b> using the same current 27-gene PCA15 input as UMAP (PCA initialization, perplexities 30 and 300, high learning rate). On an independent fixed 4,000-cell evaluation subset, local Fine26 purity improved from 0.182 (Stage704 UMAP) to <b>0.221</b> (new t-SNE), while high-dimensional 15-neighbor retention rose from 0.118 to <b>0.212</b>; broad-class purity remained around 0.72. These are descriptive structure checks, not proof of new clusters.</p>
<div class="cards">
<article class="card"><h3>Primary: VGAT / VGLUT2 on UMAP and true t-SNE</h3><a href="assets/stage982_tsne_comparison/STAGE982_PAPER_GT10_ALLCELL_FOURCLASS_UMAP_TSNE.png"><img loading="lazy" src="assets/stage982_tsne_comparison/STAGE982_PAPER_GT10_ALLCELL_FOURCLASS_UMAP_TSNE.png"></a><p>Same full 71,950 cells; Yuhan-style raw <b>&gt;10 spots/cell</b> calls: VGAT-only 20,709; VGLUT2-only 4,992; co-positive 6,116 (8.5%); neither 40,133 (55.8%). A small double-positive fraction does not demonstrate perfect exclusion when most cells remain below both thresholds. <a href="data/stage982_paper_gt10_fourclass_counts.csv">Counts</a></p></article>
<article class="card"><h3>Frozen Fine26: old versus new geometry</h3><a href="assets/stage982_tsne_comparison/multiscale_30_300_COMPARISON.png"><img loading="lazy" src="assets/stage982_tsne_comparison/multiscale_30_300_COMPARISON.png"></a><p>Upper row: same frozen 26 identities. Lower row: independent exploratory <em>normalized</em> cutoff sensitivity, not the primary >10 raw-count classification. <a href="data/stage982_multiscale_30_300_embedding_quality.csv">Full-dataset comparison metrics</a></p></article>
<article class="card"><h3>Revisited Stage127 historical t-SNE · exactly matched cells</h3><a href="assets/stage985_historical_tsne/HISTORICAL_TSNE_VS_CURRENT_EMBEDDINGS_MATCHED.png"><img loading="lazy" src="assets/stage985_historical_tsne/HISTORICAL_TSNE_VS_CURRENT_EMBEDDINGS_MATCHED.png"></a><p>35,249 historical cells matched by ROI identity. Fine26 15-NN purity: historical t-SNE <b>0.543</b>; current UMAP <b>0.253</b>; current multiscale t-SNE <b>0.347</b>. Crucially, high-D neighbor retention favors current t-SNE (0.253) over historical t-SNE (0.095), so the visually separated historical islands can distort present-day molecular geometry. <a href="data/stage985_historical_tsne_quality.json">Matched-sample audit</a></p></article>
</div>
<p class="lead">Method anchors: Kobak &amp; Berens, <a href="https://doi.org/10.1038/s41467-019-13056-x">Nature Communications 2019</a> for PCA-init/multiscale t-SNE; Yuhan Wang et al., <a href="https://doi.org/10.7554/eLife.84262">eLife 2023</a> for EASI-FISH cell-type mapping and marker-coexpression design. The older Stage127 t-SNE used an earlier EITv2 cell-classification stage; it is not an independent validation set.</p>
</section>
'''
if 'id="tsne982"' not in h:
 assert '<section id="local">' in h
 h=h.replace('<section id="local">',block+'<section id="local">',1)
 h=h.replace('<a href="#threshold981">VGAT/VGLUT2 QC</a>','<a href="#threshold981">VGAT/VGLUT2 QC</a><a href="#tsne982">t-SNE comparison</a>')
f.write_text(h,encoding="utf8")
print("STAGE982_985_PAGE_UPDATED",h.count('id="tsne982"'),flush=True)
