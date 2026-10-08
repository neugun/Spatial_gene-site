from pathlib import Path
import re
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
F=P/"single-cell.html";h=F.read_text(encoding="utf8")
root="assets/stage1014_doublepositive_spatial/"
def card(title,file,caption):
 p=root+file
 assert (P/p).is_file(),p
 return f'<article class="card"><h3>{title}</h3><a href="{p}"><img loading="lazy" src="{p}" alt="{title}"></a><p>{caption}</p></article>'
cards=[
card("Three COMPLETE sections · VGAT>10 VGLUT2>3, four-class","STAGE1014_raw_VGAT10_VGLUT2_3_FULL_SECTION_DUAL_SPATIAL.png","True cell centers from Stage631, fully x/y flipped already; purple = double-positive. Dashed tile grid is a 206 µm PITCH PROXY, not registered actual acquisition seams."),
card("Where are the double-positive cells? No grid overlay","STAGE1014_DUAL_ONLY_FULL_SECTION_NO_GRID.png","S500/S530/S560 entire measured fields; purple cells over neutral whole-tissue context."),
card("Cell-density and Fine26-adjusted dual spatial fractions","STAGE1014B_DUAL_FRACTION_AND_FINE26_ADJUSTED_SPATIAL.png","65 µm spatial bins; >=20 previously classified neurons per bin. Top = actual double-positive fraction; bottom = observed minus Fine26-class-composition-expected fraction. Section-wide color ranges."),
card("Same XY, VGAT>10 and VGLUT2>10 stringent comparator","STAGE1014_raw_VGAT10_VGLUT2_10_FULL_SECTION_DUAL_SPATIAL.png","Retains 6,040 double-positive reference neurons, 8.81% of 68,572."),
card("Same XY, VGAT>15 and VGLUT2>5 comparator","STAGE1014_raw_VGAT15_VGLUT2_5_FULL_SECTION_DUAL_SPATIAL.png","9,883 double-positive reference neurons, 14.41%."),
card("Near vs far from proposed borders","STAGE1014_EDGE_AND_NEIGHBOR_ENRICHMENT.png","Compare provisional tile-grid edges, tissue bounding rectangle and closeness to two exclusive-marker cell classes. Section-specific enrichment ratios; no animal-level p values."),
card("CRITICAL: assumed tile-grid phase sensitivity","STAGE1014_UNCERTAIN_TILE_PHASE_SENSITIVITY.png","Moving the proxy 5×5 grid changes seam enrichment. This analysis does not establish the real seam position; original stitch transform needed before concluding seam artifacts are absent."),
card("Normalized gene-specific P70 artifact control","STAGE1014_normalized_geneP70_FULL_SECTION_DUAL_SPATIAL.png","Not a validated E/I classifier: dividing by a 27-gene marker panel may create class-specific denominator artifacts.") ]
new='<section id="stage1014-dual-spatial"><h2>Stage1014 · VGAT/VGLUT2 co-detection in physical XY, tile and biological borders</h2><p class="lead"><b>First full-section maps of real dual-marker cells.</b> Of 68,572 previously defined reference neurons, 16,768 (24.45%) are double-positive with corrected VGAT&gt;10 / VGLUT2&gt;3; 6,040 (8.81%) remain double-positive using &gt;10/&gt;10, and 9,883 (14.41%) using &gt;15/&gt;5. Current criterion is exploratory and is <em>not</em> proof of dual neurotransmitter release. The 3 sections S500/S530/S560 originate from one mouse.</p><p class="lead">At the preselected assumed 206 µm grid phase, reference-neuron double-positive rates within 15 µm of a supposed tile seam versus farther away are S500 24.17% vs 25.97%, S530 23.05% vs 26.15%, S560 21.47% vs 21.89%. This is <b>not an observed consistent seam enrichment at that assumed phase</b>. However, phase is not validated against native acquisition tile boundaries. Results include the sensitivity sweep; no robust no-seam claim can yet be made. Co-detection could still arise from segmentation or true biological mixtures, and neighbor/tissue-border analyses are only descriptive.</p><div class="cards">'+''.join(cards)+'</div><p class="lead">Download QC summary: <a href="data/stage1014_double_positive_counts_all_definitions_gates.csv">class counts per section and neuron gate</a> · <a href="data/stage1014_double_positive_tile_and_cell_border_proxies.csv">edge/neighbor exposure-normalized prevalence</a> · <a href="data/stage1014_tile_phase_sensitivity_scan.csv">unknown-grid-phase sensitivity</a> · <a href="data/stage1014_dual_spatial_authority.json">definitions and caveats</a>. All per-cell source positions, source gene counts and class labels remain <b>private</b> on the workstation.</p></section>'
if 'id="stage1014-dual-spatial"' not in h:
 assert '<section id="stage988-percentile">' in h
 h=h.replace('<section id="stage988-percentile">',new+'<section id="stage988-percentile">',1)
 # navigation insert
 h=h.replace('<a href="#stage988-percentile">','<a href="#stage1014-dual-spatial">VGAT/VGLUT2 spatial double+</a><a href="#stage988-percentile">',1)
imgs=re.findall(r'<img[^>]+src="([^"]+)"',h)
assert all((P/x).is_file() for x in imgs if x.startswith("assets/"))
F.write_text(h,encoding="utf8")
print("STAGE1014_SITE_READY",h.count('id="stage1014-dual-spatial"'),len(imgs),flush=True)
