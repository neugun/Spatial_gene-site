from pathlib import Path
import pandas as pd
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review');F=P/'single-cell.html';s=F.read_text(encoding='utf8')
path='assets/stage1043_strict_transmitter_gate/STAGE1046_NEW_TSNE_USER_SELECTED_CELLS.png'
assert (P/path).is_file()
t=pd.read_csv(P/'data'/'stage1046_strict_gate_old_vs_new_embedding_knn.csv')
print('PURITY_QC',t.to_string(index=False),flush=True)
k='<section id="stage1043-user-strict-gate">'
assert k in s
detail=('<section id="stage1046-gated-new-tsne"><h2>New t-SNE refitted on ONLY selected 31,115 cells: apparent identities and sectional effects</h2>'
'<p class="lead">In addition to applying the exact gate as a post hoc color overlay on the full-data Stage995 t-SNE, we independently refitted label-free unsupervised PCA20 / t-SNE using only the 31,115 selected rows (27 gene panel). No broad/Fine26/class labels supplied to optimization. <b>15-nearest-neighbor transmitter category concordance</b> changes from 0.660 (filtered old map) to 0.686 (gated refit), while <b>section-neighbor agreement</b> rises from 0.753 to 0.772. The slight increase in three-class concordance is not a biological validation because selection and displayed transmitter labels depend on VGAT/VGLUT2 expression; the stronger section association is an important potential batch/section geography confound to check before treating new islands as independent subtypes.</p>'
'<div class="cards"><article class="card"><h3>Newly refitted three transmitter classes and section-dependence</h3><a href="'+path+'"><img loading="lazy" src="'+path+'"></a>'
'<p>31,115 selected cells, matched identical 27-gene measurement panel; same-scale point plots by marker class and S500/S530/S560. Cluster positioning is not an independent GABA-vs-Glu identity discovery. <a href="data/stage1046_strict_gate_old_vs_new_embedding_knn.csv">Same-population 15NN comparison</a> · <a href="data/stage1046_strict_gate_tsne_authority.json">Exact t-SNE preprocessing</a>.</p>'
'</article></div></section>')
if 'stage1046-gated-new-tsne' not in s:s=s.replace(k,detail+k,1)
F.write_text(s,encoding='utf8');print('STAGE1046_PAGE_UPDATED',len(s),flush=True)
