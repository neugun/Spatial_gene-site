"""Stage997 promote independently audited unsupervised gene-standardized t-SNE and genuine
S500/S530/S560 segmentation meshes across entire page."""
from pathlib import Path
import re,json
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
F=P/"single-cell.html";h=F.read_text(encoding="utf8")
h=h.replace("an unsupervised 27-gene PCA15 -> multiscale t-SNE","an unsupervised gene-wise log1p/standardized 27-gene PCA20 -> multiscale t-SNE")
h=h.replace("All t-SNE panels use one frozen 71,950-cell embedding.","All t-SNE panels use one 71,950-cell embedding (Stage995: log1p corrected counts, per-gene scaled, PCA20, multiscale 30/100, no label supervision).")
cards='''<article class="card"><h3>Fine26 upgraded t-SNE · new main geometry</h3><a href="assets/stage989_fine26_tsne/FINE26_CATEGORICAL_TSNE.png"><img loading="lazy" src="assets/stage989_fine26_tsne/FINE26_CATEGORICAL_TSNE.png" alt="Stage995 improved t-SNE Fine26 all cells"></a><p>Gene-wise log1p 27-count preprocessing, clip at p99.7, gene-wise z-score, PCA20, perplexity 30/100 t-SNE; no class labels in fitting.</p></article>
<article class="card"><h3>Stage982 versus improved Stage995 · frozen label check</h3><a href="assets/stage995_tsne_preprocess/STAGE995_ORIGINAL_VS_LOGZ_TSNE.png"><img loading="lazy" src="assets/stage995_tsne_preprocess/STAGE995_ORIGINAL_VS_LOGZ_TSNE.png" alt="Both label-free t-SNE geometries side by side"></a><p>On the <em>same</em> 4,000 randomly selected cells, Fine26 15-NN purity <b>0.227 → 0.435</b>; Broad purity <b>0.723 → 0.826</b>. Each embedding preserves its own high-dimensional preprocessing neighborhoods (0.213 versus 0.233); different gene preprocessing prevents treating these two last figures as a matched-data control. <a href="data/stage996_tsne_preprocess_quality_comparison.json">All baseline/cross-input scores</a>.</p></article>
<article class="card"><h3>Broad classes on improved t-SNE</h3><a href="assets/stage989_fine26_tsne/BROAD_CLASSES_TSNE.png"><img loading="lazy" src="assets/stage989_fine26_tsne/BROAD_CLASSES_TSNE.png"></a><p>The broad/fine labels remain frozen. More distinct islands are improved visualization of molecular neighborhoods, not proof of different neurotransmitter release.</p></article>'''
new='''<section id="tsne982"><h2>Stage995 · improved unsupervised full-data t-SNE</h2>
<p class="lead">The original PCA15 embedding retained substantial molecular overlap. A fully independent, still-unsupervised fit uses log1p(corrected 27-gene counts), gene-wise p99.7 winsorization, gene-wise standardization, PCA20 (92.2% input variance), PCA-initialized t-SNE with 30/100 multiscale perplexities and high learning rate. Exactly the same 71,950 cells and frozen Fine26 labels are displayed. We compared the new embedding on the same sampled cells without ever feeding Fine26 or Broad labels to optimization.</p>
<div class="cards">'''+cards+'''</div>
<p class="lead">Controls: Stage982 PCA15/30–300 had frozen-label local Fine26 purity 0.227, Broad 0.723. Stage995 has Fine26 0.435, Broad 0.826. Fifteen-neighbor retention relative to each method's <em>own</em> PCA input improves 0.213 → 0.233; however cross-input fidelity differs, and the two are not a controlled parameter-only t-SNE comparison. The Stage127 historical embedding can look more separated but retains less of its comparison high-dimensional neighborhood. Methods: <a href="https://doi.org/10.1038/s41467-019-13056-x">Kobak &amp; Berens 2019</a>. <a href="data/stage995_preprocessing_qa.json">Training preprocessing</a> · <a href="data/stage996_tsne_preprocess_quality_comparison.json">Metrics</a>. These are molecular visualization results from one specimen, not animal-level replications.</p></section>'''
h,n=re.subn(r'<section id="tsne982">.*?</section>',lambda m:new,h,count=1,flags=re.S);assert n==1
old='''<article class="card"><h3>S500 LOCAL: real segmentation surface vs sparse maskcloud</h3>'''
assert old in h
# Don't duplicate if script is rerun.
append='''<article class="card"><h3>S530 LOCAL: exact segmented masks converted to smooth surfaces</h3><a href="assets/stage992_true_mask_surface/STAGE992_S530_ACTUAL_SEGMENTATION_SURFACE_VS_POINTCLOUD.png"><img loading="lazy" src="assets/stage992_true_mask_surface/STAGE992_S530_ACTUAL_SEGMENTATION_SURFACE_VS_POINTCLOUD.png"></a><p>131 real downseg labels, 21 Fine26 types. Match original sampled voxel cloud to direct per-ROI marching cubes. <a href="data/stage992_true_mask_surface_qc_S530.csv">Per-mask voxel and mesh QA</a>.</p></article>
<article class="card"><h3>S560 LOCAL: exact segmented masks converted to smooth surfaces</h3><a href="assets/stage992_true_mask_surface/STAGE992_S560_ACTUAL_SEGMENTATION_SURFACE_VS_POINTCLOUD.png"><img loading="lazy" src="assets/stage992_true_mask_surface/STAGE992_S560_ACTUAL_SEGMENTATION_SURFACE_VS_POINTCLOUD.png"></a><p>116 real downseg labels, 20 Fine26 types, matched true-surface and pointcloud views. <a href="data/stage992_true_mask_surface_qc_S560.csv">Per-mask voxel and mesh QA</a>.</p></article>'''
m=re.search(r'(<section id="maskbody">.*?<div class="cards">)(.*?)(</div>)',h,flags=re.S)
assert m and 'STAGE992_S530' not in m.group(2)
h=h[:m.start(2)]+m.group(2)+append+h[m.end(2):]
h=h.replace("S500 surface QC: <a href='data/stage992_true_mask_surface_qc.csv'>per-ROI voxel and triangle counts</a>.","S500/S530/S560 surface QC: <a href='data/stage992_true_mask_surface_qc.csv'>S500</a>, <a href='data/stage992_true_mask_surface_qc_S530.csv'>S530</a>, <a href='data/stage992_true_mask_surface_qc_S560.csv'>S560</a> per-ROI voxel and triangle counts.")
h=h.replace("S500 surface meshes originate from full original downseg.tif segmentation labels","Local surface meshes originate from full original downseg.tif segmentation labels")
# Method provenance updated, local spot distributions and target remain unchanged.
assert 'stage989_fine26_tsne/' in h and 'stage995_tsne_preprocess/' in h
bad=[p for p in re.findall(r'<img[^>]+src="([^"]+)"',h) if not (P/p).exists()]
assert not bad,bad
F.write_text(h,encoding="utf8")
print("STAGE997_PRIMARY_TSNE_AND_ALL_LOCAL_SURFACES",len(h),flush=True)
