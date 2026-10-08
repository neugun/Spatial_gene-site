"""Stage993: promote full t-SNE, whole XY section, orthogonal 3D, genuine S500 mesh to public page.
Keep scientific boundaries; hide older UMAP pages and cropped "full" maps without rewriting data authority.
"""
from pathlib import Path
import re,json
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
F=P/"single-cell.html";h=F.read_text(encoding="utf8")
A=P/"assets";D=P/"data"
def img(rel,alt):
 assert (P/rel).exists(),rel
 return f'<a href="{rel}" target="_blank"><img loading="lazy" src="{rel}" alt="{alt}" /></a>'
def card(title,rel,description):
 return f'<article class="card"><h3>{title}</h3>{img(rel,title)}<p>{description}</p></article>'
def section(sid,head,intro,cards,footer=""):
 return f'<section id="{sid}"><h2>{head}</h2><p class="lead">{intro}</p><div class="cards">{"".join(cards)}</div>{footer}</section>'
def replace(sid,content):
 global h
 pat=re.compile(r'<section id="'+re.escape(sid)+r'">.*?</section>',flags=re.S)
 h,n=pat.subn(lambda _:content,h,count=1)
 assert n==1,sid
rt="assets/stage989_fine26_tsne/"
ss="assets/stage990_fullsection/"
th="assets/stage988_gene_percentile/"
orth="assets/stage991_orthogonal_maskbody/"
mesh="assets/stage992_true_mask_surface/"
replace("taxonomy",section("taxonomy","Fine26 molecular atlas — full-data t-SNE",
 "All 71,950 cells are displayed on an unsupervised 27-gene PCA15 -> multiscale t-SNE. Frozen Fine26 identity and gene expression have not been reclustered or relabeled. The former UMAP plots are no longer the primary presentation.",[
 card("Fine26 t-SNE · 26 consistent categorical colors",rt+"FINE26_CATEGORICAL_TSNE.png","Identical palette between cells in t-SNE, complete XY section and 3D true mask bodies."),
 card("All 27-gene t-SNE feature plots",rt+"ALL27_GENE_TSNE_MONTAGE.png","Every measured gene. Same t-SNE and 71,950-cell universe; each panel has its own numeric count colorbar."),
 card("Fine26 marker dotplot · frozen cell classes","assets/stage950_singlecell/FINE26_DOTPLOT_CURRENT.png","27-gene molecular marker summary, independent of 2D visualization."),
 card("Leave-one-section-out identity transfer","assets/stage950_singlecell/FINE26_CROSSSECTION_TRANSFER.png","Identity transportability between S500/S530/S560, not independent animal validation.")
 ]))
replace("transmitters",section("transmitters","Fine26 t-SNE, broad classes and four transmitter markers",
 "All t-SNE panels use one frozen 71,950-cell embedding. Slc32a1=VGAT, Slc17a6=VGLUT2, Slc6a2=NE transporter, Slc5a7=cholinergic transporter (not ChAT enzyme itself). Direct count overlays, depth-normalized signals and top-5% enrichment are separate biological questions. Gene detection alone does not prove cell identity.",[
 card("Fine26 · 26 high-contrast t-SNE colors",rt+"FINE26_CATEGORICAL_TSNE.png","Frozen Fine26; <a href='data/stage977_fine26_color_key.csv'>26-color key</a>."),
 card("Broad E / I / NE / cholinergic / other · t-SNE",rt+"BROAD_CLASSES_TSNE.png","Same cells, no UMAP overlays or reclassification."),
 card("VGAT/VGLUT2/NE/cholinergic · original corrected count",rt+"VGAT_VGLUT2_NE_CHAT_MARKER_TSNE_RAW.png","Log(1+corrected gene spot-count equivalents); high background positivity does not imply transmitter class."),
 card("VGAT/VGLUT2/NE/cholinergic · normalized expression",rt+"VGAT_VGLUT2_NE_CHAT_MARKER_TSNE_NORM.png","Log(1+10k depth-normalized expression), same embeddings and independent gene-specific numeric colorbars."),
 card("Top 5% normalized gene signal · t-SNE",rt+"VGAT_VGLUT2_NE_CHAT_TOP5PCT_TSNE.png","Per-gene 95th percentile, descriptive enrichment only."),
 card("Top 5% gene enrichment across Broad classes","assets/stage989_fine26_tsne/TRANSMITTER_TOP5_BROAD_ENRICHMENT.png","Molecular enrichment is quantified on original data and independent of t-SNE projection.")
 ]))
audit=json.loads((D/"stage988_gene_percentile_authority.json").read_text(encoding="utf8"))
audit_raw=pd=None
replace("threshold981",section("threshold981","Original threshold / source-data audits (archived controls)",
 "Older normalized-threshold and >10 raw-count screens are retained as <em>data tables</em>, not mixed into the promoted percentile or t-SNE primary results. Low co-expression fractions from high cutoffs alone do not prove molecular exclusivity.",[
 card("Threshold sensitivity by gene and E/I group","assets/stage986_marker_threshold/STAGE986_VGLUT2_THRESHOLD_SENSITIVITY.png","VGAT >10; VGLUT2 thresholds 0–20, all cells and section-specific prevalence.")
 ],"<p class='lead'>Data: <a href='data/stage981_normalized_threshold_sweep.csv'>exploratory normalized scan</a>; <a href='data/stage984_h5ad_vs_routeA_source_audit.csv'>original H5AD versus Route A</a>.</p>"))
replace("tsne982",section("tsne982","Unsupervised t-SNE geometry and historical comparison",
 "New 71,950-cell t-SNE: PCA initialization, multiscale perplexities 30/300, high learning rate; no broad/Fine26 labels supplied to the embedding. On 4,000 fixed sampled cells, Fine26 15-NN purity is 0.221 versus the previous embedding baseline 0.182; high-dimensional-neighbor overlap improves from 0.118 to 0.212. A second unsupervised run (30/100 perplexities) yields Fine26 0.219, Broad 0.725, high-dimensional overlap 0.212. We select 30/300 for Fine26 presentation, but the classes still overlap in real expression space.",[
 card("Current Fine26 full-data t-SNE",rt+"FINE26_CATEGORICAL_TSNE.png","No label-supervised embedding: cluster geometry is inferred solely from expression."),
 card("Broad identity mapping on the same t-SNE",rt+"BROAD_CLASSES_TSNE.png","Display labels stay frozen; brighter islands are not newly discovered molecular types."),
 card("Historical EITv2 and current t-SNE — quantitative audit",rt+"FINE26_CATEGORICAL_TSNE.png","Previously recovered Stage127 t-SNE has higher historical identity local-purity but considerably lower high-D neighbor fidelity (0.095 vs 0.253 for current t-SNE on matched 35,249 cells). <a href='data/stage985_historical_tsne_quality.json'>Complete matched-cell scores</a>.")
 ],"<p class='lead'>Methods: <a href='https://doi.org/10.1038/s41467-019-13056-x'>Kobak &amp; Berens 2019</a> and <a href='https://doi.org/10.7554/eLife.84262'>Yuhan Wang et al. 2023</a>. High-dimensional neighborhood metrics use a fixed 4,000-cell set; 3 sections are from one mouse. <a href='data/stage982_multiscale_30_300_embedding_quality.csv'>30/300 metrics</a>; <a href='data/stage982_multiscale_30_100_embedding_quality.csv'>30/100 metrics</a>.</p>"))
# Promote actual spatial FULL sections and 3D with unambiguous metadata; local windows are separate.
replace("local",section("local","Entire-section Fine26 / Broad identities and clearly marked local fields",
 "The two main XY atlases now contain <strong>all 31,461 S500 + 20,513 S530 + 19,976 S560 cells</strong>. They are not restricted to LC/periLC, polygon, or medial-subregion masks. Both physical x and y are already flipped in the Stage631 coordinate authority; no second axis transform is applied. Fine26 colors match the promoted t-SNE exactly.",[
 card("COMPLETE section Fine26 identities · 26-color t-SNE palette",ss+"STAGE990_FINE26_FULL_SECTION_XY.png","All 71,950 cells across all three measured XY sections; every Fine26 class included."),
 card("COMPLETE section Broad molecular classes",ss+"STAGE990_BROAD_FULL_SECTION_XY.png","Exact same complete fields and coordinates, with E/I/NE/cholinergic/other labels."),
 card("S500 · local 27-gene multiplex, explicitly cropped","assets/stage974_preflipped_reference/S500_LOCAL_MULTIPLEX_STAGE631_PREFLIPPED.png","This supplementary field is a deliberate local crop, not the full-section atlas."),
 card("S530 · local 27-gene multiplex, explicitly cropped","assets/stage974_preflipped_reference/S530_LOCAL_MULTIPLEX_STAGE631_PREFLIPPED.png","Deliberate local crop."),
 card("S560 · local 27-gene multiplex, explicitly cropped","assets/stage974_preflipped_reference/S560_LOCAL_MULTIPLEX_STAGE631_PREFLIPPED.png","Deliberate local crop.")
 ]))
replace("maskbody",section("maskbody","Genuine 3D segmentations — XY principal plus XZ / YZ anatomy",
 "Primary view is the complete coronal XY tissue plane. XZ and YZ orthogonal projections give the physical depth of the segmented bodies, with no odd camera angle or geometric stretching. Each section uses 200,000 subsampled <em>true downseg.tif label voxels</em> registered to Stage631 by ROI (XY agreement r > 0.99988); whole-section spatial scope is preserved.",[
 card("Fine26 3D true-body voxel atlas · XY / XZ / YZ",orth+"STAGE991_FINE26_XY_XZ_YZ_TRUE_MASKBODY.png","Whole S500/S530/S560 measured sections, 26-color palette identical to Fine26 t-SNE."),
 card("Broad-class 3D body atlas · XY / XZ / YZ",orth+"STAGE991_BROAD_XY_XZ_YZ_TRUE_MASKBODY.png","Whole-section physical XY, with XZ/YZ depth views."),
 card("S500 LOCAL: real segmentation surface vs sparse maskcloud",mesh+"STAGE992_S500_ACTUAL_SEGMENTATION_SURFACE_VS_POINTCLOUD.png","110 actual S500 downseg label surfaces, spanning 18 frozen Fine26 groups; gentle 0.82-voxel Gaussian anti-alias and marching cubes. No fake ellipsoid or synthetic sphere. Deliberate local subset shown and labelled, NOT all S500 cells.")
 ],"<p class='lead'>S500 surface QC: <a href='data/stage992_true_mask_surface_qc.csv'>per-ROI voxel and triangle counts</a>. Surface meshes originate from full original downseg.tif segmentation labels (not sparse points). Orthogonal whole-section panels use sampled true voxel cloud and therefore should not be described as full-resolution mesh surfaces.</p>"))
# Prior Fig C comparison displayed old UMAP; replace with tSNE-only three cutoffs.
h=h.replace("Figure C · UMAP vs multiscale t-SNE, three cutoff options","Figure C · full-data t-SNE, VGAT&gt;10 / VGLUT2&gt;2,3,10")
h=h.replace("assets/stage986_marker_threshold/STAGE986_VGAT10_VGLUT2_2_3_10_UMAP_TSNE.png",rt+"FOURCLASS_VGAT10_VGLUT2_2_3_10_TSNE.png")
h=h.replace("VGAT VGLUT2 thresholds 2 3 10 on UMAP and tSNE","VGAT VGLUT2 thresholds 2 3 10 on full-data t-SNE")
# Stage988 new front-page decision section placed before Stage986; units distinguished.
pct=section("stage988-percentile","Stage988 · gene-specific percentiles aligned across cell depth",
 "For each gene, calculate <code>10,000 × corrected gene count / summed 27-gene counts</code>. The request 'VGAT 15 and VGLUT2 5' is ambiguous unless the unit is explicit. If RAW corrected count equivalents, VGAT >15 lies at the 72.0th raw-count percentile and VGLUT2 >5 at the 64.1st; at <em>median total-27-gene depth</em> those represent normalized cutoffs 898 and 299, at the 80.3rd and 64.8th normalized percentiles. If 15 / 5 instead mean normalized counts, their percentiles are only 13.6 and 19.4. Individual cells have different total counts, so raw threshold cannot map to one fixed normalized cutoff. Shared gene-specific percentile thresholds are descriptive; they do not calibrate probe-specific false-positive rates or make VGAT/VGLUT2 biologically mutually exclusive.",[
 card("Gene-specific normalized and corrected-count CDF curves",th+"STAGE988_PERCENTILE_NORMALIZED_AND_RAW_CURVES.png","Two paired curves on all 71,950 cells; both zeros and percentiles are displayed explicitly."),
 card("RAW count VGAT=15 vs VGLUT2=5: normalized percentile conversion",th+"STAGE988_RAW15_RAW5_PERCENTILE_MAPPING.png","Median-depth conversions and empirical count CDF. <a href='data/stage988_raw15_raw5_normalized_percentile_mapping.csv'>Exact table</a>."),
 card("Shared percentile -> two distinct gene thresholds + coexpression",th+"STAGE988_COMMON_PERCENTILE_TRADEOFF.png","Percentiles 65-99 jointly evaluated; double-positive and neither fractions retained. <a href='data/stage988_shared_percentile_coexpression_tradeoff.csv'>All cutoff values</a>.")
 ],"<p class='lead'>Raw normalized mapping <a href='data/stage988_gene_specific_percentile_alignment.csv'>all sections and positive-only percentiles</a>; analysis authority <a href='data/stage988_gene_percentile_authority.json'>JSON</a>. A shared percentile controls each gene's marginal selected fraction but is <strong>not</strong> an experimental measure of transmitter identity; orthogonal negative-control spot analysis remains necessary.</p>")
if 'id="stage988-percentile"' not in h:
 marker='<section id="stage986-threshold">'
 assert marker in h
 h=h.replace(marker,pct+marker,1)
# Navigation link wiring — each existing section must be present and unique
h=h.replace('<a href="#transmitters">Fine26 + markers</a>','<a href="#transmitters">t-SNE + markers</a><a href="#stage988-percentile">Gene percentile QC</a>')
h=h.replace('Same UMAP, fixed >10 spots vs normalized rule','Frozen t-SNE and fixed >10 versus normalized rule (see current plots)')
# All former UMAP thumbnails as visible images should be gone (historical CSV/description may mention baseline).
assert h.count('id="stage988-percentile"')==1
import re
imgr=re.findall(r'<img[^>]*src="([^"]+)"',h,flags=re.I)
missing=[i for i in imgr if not (P/i).exists() and not i.startswith(('http','data:'))]
assert not missing,missing
bad=[i for i in imgr if 'umap' in i.lower()]
assert not bad,("Legacy UMAP shown in page",bad)
F.write_text(h,encoding="utf8")
(D/"stage993_public_display_authority.json").write_text(json.dumps({"stage":993,"status":"UPDATED","primary_embedding":"Stage982 full 71,950-cell t-SNE 30/300","umap_img_tags_remaining":len(bad),"all_img_src_checked":len(imgr),"section_999scope":"full S500/S530/S560","3d":"XY/XZ/YZ actual Stage141 segmented voxels","s500_surfaces":"Stage992 raw downseg label surfaces, 110 genuinely meshed ROI subset"},indent=2),encoding="utf8")
print("STAGE993_SITE_REPLACED",len(h),"IMAGES",len(imgr),"NO_UMAP_IMAGE",flush=True)

