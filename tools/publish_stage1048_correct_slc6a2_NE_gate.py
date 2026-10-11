"""Add explicit error correction to current Map10 page without deleting prior quantitative audits.
As of Stage1048, NO validated NE gene cutoff. Original NE01/Stage1044 is sensitivity only.
"""
from pathlib import Path
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review')
h=P/'single-cell.html';s=h.read_text(encoding='utf8')
if 'id="stage1048-ne-marker-correction"' in s:
 print('SITE_ALREADY_PATCHED');raise SystemExit(0)
section="""<section id="stage1048-ne-marker-correction"><h2>Stage1048: Important correction — Slc6a2 expression is not an independently validated NE-cell identity</h2>
<p class="lead"><b>WITHDRAW previous 'Slc6a2&gt;200 is an accurate NE classification' inference.</b> Slc6a2 codes for the noradrenaline transporter NET, an informative biological marker. However, the NE01 labels were themselves made using the current 27-gene FISH panel, including Slc6a2. Comparing a threshold back to NE01 is circular. Existing 97.6% concordance cannot be reported as independent specificity. NE-positive neuronal identity is not proven by this gene alone. The 27-gene FISH panel lacks independent Dbh/Th/Slc18a2 measurements, so no current numerical raw-corrected-count threshold is validated for a definitive LC-NE class.</p>
<p class="lead">Direct check: among 50,630 Snap25&gt;5 cells, Slc6a2&gt;5 selects 11,892, only 15.0% of which match frozen NE01. Slc6a2&gt;200 selects 1,715 and 97.6% match NE01, but this is intrinsically non-independent. The high-count cells are anatomically localized in three sections and deserve independent LC marker confirmation; the much more widespread low count Slc6a2 could reflect weak true signal, extracellular/nearby cell effects or detection noise. Do not delete high marker-count cells or force NE/GABA/glutamate mutual exclusivity.</p>
<div class="cards">
<article class="card"><h3>Slc6a2 entire-section distribution vs original NE01 cell map</h3><a href="assets/stage1048_Slc6a2_marker_reaudit/STAGE1048_SLC6A2_GT5_GT200_VS_FROZEN_NE01_FULL3SECTION_XY.png"><img loading="lazy" src="assets/stage1048_Slc6a2_marker_reaudit/STAGE1048_SLC6A2_GT5_GT200_VS_FROZEN_NE01_FULL3SECTION_XY.png"></a><p>Entire S500/S530/S560 physical XY. Pale cyan = Slc6a2&gt;5; red = &gt;200; purple = old NE01. Spatial agreement supports informative anatomy, <b>not independent classification or spot accuracy</b>. <a href="data/stage1048_Slc6a2_NE01_continuous_not_truth.csv">Cell-class by section controls</a>.</p></article>
<article class="card"><h3>High threshold creates apparent NE01 concordance</h3><a href="assets/stage1048_Slc6a2_marker_reaudit/STAGE1048_SLC6A2_THRESHOLD_NOT_AN_INDEPENDENT_NE_TEST.png"><img loading="lazy" src="assets/stage1048_Slc6a2_marker_reaudit/STAGE1048_SLC6A2_THRESHOLD_NOT_AN_INDEPENDENT_NE_TEST.png"></a><p>Threshold-response curves demonstrate the classification circularity. This is not a NE ROC validated against orthogonal Dbh+Th labels.</p></article>
</div>
<p class="lead"><b>Independent MapMyCells warning:</b> the old NE01 group (1,833 cells) has modal Allen class "27 MY GABA" only 27.6%, followed by "24 MY Glut" 21.3%, "20 MB GABA" 13.9%. Class bootstrap median 0.18 and 0% of cells with class bootstrap &ge;0.5. Because Allen coarse taxonomy is not an independent LC-NE subcluster label and mapping used the same restricted panel, this does NOT conclusively reject LC localization; it strongly rejects treating NE01 as a validated neurotransmitter type reference. <a href="data/stage1048_NE01_OFFICIAL_MAPMYCELLS_CLASS_MIX.csv">Full reference mix</a> and <a href="data/stage1048_NE01_gate_correction_authority.json">withdrawal authority</a>.</p>
<p class="lead"><b>What follows:</b> Report Slc6a2 fluorescence and corrected count-equivalents as a continuous marker; show NE01 only as <em>existing provisional class</em>. Freeze VGAT/VGLUT2 thresholds independent of this NE01 cutoff; conduct separate double-positive and neither-positive cell audits without forcing an expected NE overlap. Test true Dbh/Th or Dbh/Th/Slc18a2 on independent transcriptomic reference and (ideally) additional multiplex FISH or immunostaining, with original source image quality controls. The earlier Stage1044 Slc6a2&gt;200 sensitivity calculation is retained for provenance but is <b>not a final classification</b>.</p></section>"""
target='<section id="stage1044-ne-calibration">'
assert target in s
s=s.replace(target,section+target,1)
s=s.replace('<section id="stage1044-ne-calibration"><h2>',
            '<section id="stage1044-ne-calibration"><p class="lead"><b>Superseded NE interpretation — see <a href="#stage1048-ne-marker-correction">Stage1048 correction</a>. All Slc6a2 cutoffs below are exploratory, not an independently validated NE gate.</b></p><h2>',1)
nav='<nav><a href="#stage1044-ne-calibration">'
if nav in s:
 s=s.replace(nav,'<nav><a href="#stage1048-ne-marker-correction">Correction: NE marker not validated</a><a href="#stage1044-ne-calibration">',1)
h.write_text(s,encoding='utf8')
print('STAGE1048_SITE_READY',len(s),flush=True)
