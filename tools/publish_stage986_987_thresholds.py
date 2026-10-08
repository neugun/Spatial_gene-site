"""Publish Stage986/987 without erasing prior reference or faking co-expression exclusion."""
from pathlib import Path
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
F=P/"single-cell.html"
h=F.read_text(encoding="utf8")
assert 'id="threshold981"' in h
new='''<section id="stage986-threshold">
<h2>Stage986 · Per-probe VGAT/VGLUT2 thresholds: >10 versus >2/3 corrected counts</h2>
<p class="lead"><strong>Marker-specific sensitivity test rather than an imposed E/I separation.</strong> VGAT (Slc32a1) stays at >10 count-equivalents per cell. VGLUT2 (Slc17a6) is compared at >2, >3, and >10. All numbers below come from 71,950 real cells in S500/S530/S560, from the same single specimen, keeping original Fine26 and broad labels frozen. Current post-correction Route A counts are numerical, not independently validated physical puncta counts; hence thresholds refer to corrected count-equivalents.</p>
<table style="border-collapse:collapse;width:100%;max-width:950px;margin:14px auto;text-align:center"><thead><tr><th>VGLUT2 count cutoff</th><th>VGLUT2+ within broad E</th><th>VGLUT2+ within broad I</th><th>All cells: VGAT+/VGLUT2+ co-positive</th><th>All cells: neither</th></tr></thead><tbody>
<tr><td><strong>&gt;2</strong></td><td>62.8%</td><td>61.6%</td><td>27.2%</td><td>28.6%</td></tr>
<tr><td><strong>&gt;3</strong></td><td>52.4%</td><td>51.3%</td><td>23.7%</td><td>35.1%</td></tr>
<tr><td>Historical comparator &gt;10</td><td>14.8%</td><td>15.6%</td><td>8.5%</td><td>55.8%</td></tr>
</tbody></table>
<p class="lead"><strong>Key interpretation:</strong> dropping VGLUT2 to 2–3 raises the detectable cell fraction substantially, but VGLUT2 by itself does <strong>not</strong> distinguish broad excitatory from inhibitory cells (E/I AUC 0.504). VGAT is far more discriminating (inverse VGAT AUC 0.844). If VGAT-positive cells are assigned first, VGLUT2&gt;3 and VGAT≤10 picks out 19,448 neuronal cells, 72.9% of which have frozen broad-E labels; this is a <em>VGAT-priority gating rule</em>, not evidence of absent dual-marker detection. The observed co-positive cells remain counted (23.0% of the E/I universe at &gt;3). Ratios and broad labels are partly non-independent, so no validation of neurotransmitter release is claimed.</p>
<div class="cards">
<article class="card"><h3>Figure A · thresholds versus real positive percentages</h3><a href="assets/stage986_marker_threshold/STAGE986_VGLUT2_THRESHOLD_SENSITIVITY.png"><img loading="lazy" alt="VGAT fixed 10, VGLUT2 threshold sweep; sensitivity and coexpression" src="assets/stage986_marker_threshold/STAGE986_VGLUT2_THRESHOLD_SENSITIVITY.png"></a><p>Positive, exclusive, double-positive, and neither as a function of VGLUT2 count cutoff. <a href="data/stage986_vgat_vglut_threshold_all_sections_broad.csv">Complete data table</a></p></article>
<article class="card"><h3>Figure B · four-class breakdown for S500, S530, S560</h3><a href="assets/stage986_marker_threshold/STAGE986_FOUR_CLASS_10_2_3_10_PER_SECTION.png"><img loading="lazy" alt="Four-class percentages across sections" src="assets/stage986_marker_threshold/STAGE986_FOUR_CLASS_10_2_3_10_PER_SECTION.png"></a><p>Same VGAT cutoff of 10, VGLUT2 2/3/10; section-specific denominators displayed.</p></article>
<article class="card"><h3>Figure C · UMAP vs multiscale t-SNE, three cutoff options</h3><a href="assets/stage986_marker_threshold/STAGE986_VGAT10_VGLUT2_2_3_10_UMAP_TSNE.png"><img loading="lazy" alt="VGAT VGLUT2 thresholds 2 3 10 on UMAP and tSNE" src="assets/stage986_marker_threshold/STAGE986_VGAT10_VGLUT2_2_3_10_UMAP_TSNE.png"></a><p>Identical 71,950 cells, fixed embeddings; nothing reclustered, so shape cannot be influenced by threshold choice.</p></article>
<article class="card"><h3>Figure D · cross-section consistency</h3><a href="assets/stage986_marker_threshold/STAGE986_VGLUT_SECTION_ROBUSTNESS.png"><img loading="lazy" alt="Section-specific VGLUT2 positivity and thresholds" src="assets/stage986_marker_threshold/STAGE986_VGLUT_SECTION_ROBUSTNESS.png"></a><p>Different sections show different detection levels; none show a strong E-specific Slc17a6 signal at 2/3. <a href="data/stage986_vgat10_vglut_threshold_h5ad_routeA_reproducibility.csv">Historic H5AD versus current correction</a>.</p></article>
<article class="card"><h3>Figure E · gene ratio and discriminative value</h3><a href="assets/stage986_marker_threshold/STAGE987_GENE_RATIO_AND_VGLUT_SPECIFICITY.png"><img loading="lazy" alt="Gene ratio and E I ROC AUC" src="assets/stage986_marker_threshold/STAGE987_GENE_RATIO_AND_VGLUT_SPECIFICITY.png"></a><p>VGLUT-only AUC 0.504; VGAT alone 0.844; VGLUT/VGAT log-ratio 0.817. Combination primarily benefits from VGAT. <a href="data/stage987_conditional_gating_and_marker_ratio.csv">VGAT-first gating coverage and precision</a>.</p></article>
</div>
<p class="lead"><b>Publication defensibility:</b> use Slc32a1&gt;10 / Slc17a6&gt;3 as a transparently marked lower-sensitivity-threshold <em>exploratory overlay</em>, retain &gt;10/&gt;10 comparator, raw detection histograms, and separate dual-positive fractions. Require probe-specific negative-control and true punctum/background inspection before calling Slc17a6&gt;3 an experimentally validated positive threshold. Cross-sections are technical/spatial consistency checks, not animal-level replication. Method context: <a href="https://doi.org/10.1016/j.cell.2021.11.024">Wang et al. Cell 2021</a>; <a href="https://doi.org/10.7554/eLife.84262">Wang et al. eLife 2023</a>.</p>
</section>
'''
if 'id="stage986-threshold"' not in h:
 h=h.replace('<section id="threshold981">',new+'<section id="threshold981">',1)
 if '<a href="#threshold981">' in h:
  h=h.replace('<a href="#threshold981">','<a href="#stage986-threshold">VGAT/VGLUT2 2–3 counts</a><a href="#threshold981">',1)
assert h.count('id="stage986-threshold"')==1
for relative in ["assets/stage986_marker_threshold/STAGE986_VGLUT2_THRESHOLD_SENSITIVITY.png","assets/stage986_marker_threshold/STAGE986_FOUR_CLASS_10_2_3_10_PER_SECTION.png","assets/stage986_marker_threshold/STAGE986_VGAT10_VGLUT2_2_3_10_UMAP_TSNE.png","assets/stage986_marker_threshold/STAGE986_VGLUT_SECTION_ROBUSTNESS.png","assets/stage986_marker_threshold/STAGE987_GENE_RATIO_AND_VGLUT_SPECIFICITY.png","data/stage987_conditional_gating_and_marker_ratio.csv"]:
 assert (P/relative).is_file(),relative
F.write_text(h,encoding="utf8")
print("STAGE986_PAGE_READY",len(h),flush=True)
