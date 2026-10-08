from pathlib import Path
import re
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");F=P/"single-cell.html";h=F.read_text(encoding="utf8")
file="assets/stage1020_rong_rounds/STAGE1020_RONG_ROUND1_ROUND2_GABA_GLU_DOUBLE.png"
assert (P/file).exists()
if "STAGE1020_RONG_ROUND1_ROUND2_GABA_GLU_DOUBLE.png" not in h:
 m=re.search(r'(<section id="stage1018-real-tile-and-reference">.*?<div class="cards">)(.*?)(</div>)',h,flags=re.S)
 assert m
 card='<article class="card"><h3>Rong first versus second collection (true Cell Ranger files)</h3><a href="'+file+'"><img loading="lazy" src="'+file+'"></a><p>Original first-round neuron-enriched: female 2,137 cells, male 3,024; second-round mixed-sample female 28,864, male 30,382. Different cell compositions and preprocessing: ratios cannot be combined into one 25-cluster scRNA cohort. <a href="data/stage1020_rong_round1_round2_cellranger_direct_GABA_GLU.csv">All raw UMI detection counts</a>.</p></article>'
 h=h[:m.start(2)]+m.group(2)+card+h[m.end(2):]
F.write_text(h,encoding="utf8")
print("STAGE1022_PAGE_UPDATED",len(h),flush=True)
