"""Append factual NE01 cluster GABA/Glu spot count plots to existing Map10 without revising any source classification."""
from pathlib import Path
from bs4 import BeautifulSoup
import pandas as pd,json
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review')
h=P/'single-cell.html'
text=h.read_text(encoding='utf8')
secid='stage1049-ne01-cell-spots'
if f'id="{secid}"' in text:
 print('STAGE1049_ALREADY_PRESENT');raise SystemExit(0)
J=pd.read_csv(P/'data'/'stage1049_NE01_joint_VGAT_VGLUT2_threshold_counts.csv')
a=J[(J.source=='Stage943_corrected')&(J.VGAT_cutoff==10)&(J.VGLUT2_cutoff==5)].set_index('section')
old=J[(J.source=='Stage631_integer')&(J.VGAT_cutoff==10)&(J.VGLUT2_cutoff==5)].set_index('section')
S=pd.read_csv(P/'data'/'stage1049_NE01_VGAT_VGLUT2_by_section_spot_count_statistics.csv')
def stats(sec,gene):
 return S[(S.section.astype(str)==str(sec))&(S.gene==gene)&(S.source=='Stage943_RouteA_corrected')].iloc[0]
a0=a.loc['ALL'];b0=old.loc['ALL']
assert int(a0.n)==1833 and int(a0.both)==760 and int(b0.both)==690
parts=[
f'<section id="{secid}">',
'<h2>Stage1049 — VGAT and VGLUT2 spot counts INSIDE the visibly separated NE01 t-SNE cluster</h2>',
'<p class="lead">Start from the existing unsupervised full-data Stage995 t-SNE and its <b>frozen 1,833 NE01 cells</b>; neither Slc6a2 cutoff nor new transcript labels are used to choose cells. We directly extracted the paired per-cell Slc32a1 (VGAT) and Slc17a6 (VGLUT2) values from original earlier Stage631 integer per-cell count assignment and latest Stage943 RouteA corrected count-equivalents. NE01 is a visually distinct molecular island, not an independently proven catecholaminergic phenotype. These source gene spot counts are the thing to validate, not a target to suppress by thresholds.</p>',
'<div class="cards">',
'<article class="card"><h3>Where the NE01 island sits — per-cell VGAT and VGLUT2 counts</h3><a href="assets/stage1049_NE01_GABA_GLU_cell_spots/STAGE1049_NE01_TSNE_VGAT_VGLUT2_OLD_VS_CURRENT_PER_CELL.png"><img loading="lazy" src="assets/stage1049_NE01_GABA_GLU_cell_spots/STAGE1049_NE01_TSNE_VGAT_VGLUT2_OLD_VS_CURRENT_PER_CELL.png" alt="Separate NE01 tSNE group per-cell VGAT and VGLUT2 spot count from earlier integer and RouteA corrected"></a><p>Top: Slc32a1 (VGAT), bottom: Slc17a6 (VGLUT2); each point is one actual NE01 cell. Low count pale cyan; high count dark red; zeros grey. Same frozen embedding and cell selection for both genes.</p></article>',
'<article class="card"><h3>Distributions and threshold sensitivity within 1,833 NE01 cells</h3><a href="assets/stage1049_NE01_GABA_GLU_cell_spots/STAGE1049_NE01_VGAT_VGLUT2_CELL_COUNT_DISTRIBUTIONS.png"><img loading="lazy" src="assets/stage1049_NE01_GABA_GLU_cell_spots/STAGE1049_NE01_VGAT_VGLUT2_CELL_COUNT_DISTRIBUTIONS.png"></a><p>Both pre-correction Stage631 integer counts and current RouteA values are shown. Section-specific empirical CDFs and threshold sensitivity highlight that VGAT/VGLUT2 within NE01 are NOT one clean uniformly zero distribution.</p></article>',
'<article class="card"><h3>Same NE01 cells, jointly measured VGAT and VGLUT2 per section</h3><a href="assets/stage1049_NE01_GABA_GLU_cell_spots/STAGE1049_NE01_VGAT_VGLUT2_PER_CELL_JOINT_COUNTS.png"><img loading="lazy" src="assets/stage1049_NE01_GABA_GLU_cell_spots/STAGE1049_NE01_VGAT_VGLUT2_PER_CELL_JOINT_COUNTS.png"></a><p>Actual paired cell counts, not independent marginal overlap approximation. Dashed threshold lines VGAT&gt;10 / VGLUT2&gt;5; all three sections are separate anatomical samples, not independent animals.</p></article>',
'</div>',
'<h3>NE01 cell statistics (corrected RouteA per-cell count-equivalents)</h3>',
'<table><thead><tr><th>Section</th><th>NE01 n</th><th>VGAT median</th><th>VGLUT2 median</th><th>VGAT&gt;10</th><th>VGLUT2&gt;5</th><th>Both positive</th></tr></thead><tbody>'
]
for sec in [500,530,560,'ALL']:
 r=a.loc[str(sec)]
 vg=stats(sec,'VGAT');ve=stats(sec,'VGLUT2')
 parts.append(f'<tr><td>{str(sec) if sec=="ALL" else "S"+str(sec)}</td><td>{int(r.n):,}</td>'
              f'<td>{vg["median"]:.2f}</td><td>{ve["median"]:.2f}</td>'
              f'<td>{r.pct_vgat:.1f}%</td><td>{r.pct_vglut2:.1f}%</td>'
              f'<td>{int(r.both):,} ({r.pct_both:.1f}%)</td></tr>')
parts+=['</tbody></table>',
'<p class="lead">At VGAT&gt;10/VGLUT2&gt;5, paired NE01 calls are: <b>335 VGAT-only, 226 VGLUT2-only, 760 double-positive, 512 below both cutoffs</b>. Importantly, previous Stage631 integer count assignments ALREADY give <b>690/1,833 (37.6%)</b> double positive at the same cutoffs; current RouteA yields <b>760/1,833 (41.5%)</b>. Thus corrected matrix processing is NOT the sole origin of the overlap, although old integer counts also need raw RS-FISH source and cross-round cell assignment validation.</p>',
'<p class="lead"><b>Interpretation:</b> The NE01 island is strikingly separate in molecular t-SNE, largely because of Slc6a2 expression, yet VGAT/VGLUT2 spots frequently appear on cells assigned to it, with strong differences S500/S530/S560. This discrepancy is not resolved by choosing an arbitrary high transmitter cutoff, and t-SNE visual separation cannot certify true cell-level GABA/Glu co-transmission. Audit 3D original same-cell masks, registered R7/R9 spot centroids, local cell volume and boundaries, and NE-marker source fluorescence in each section before making any final neurotransmitter claim. A separate Dbh+Th/Slc18a2 experiment is required for independent NE ground truth.</p>',
'<p class="lead"><a href="data/stage1049_NE01_all1833_cells_VGAT_VGLUT2_original_baseline_and_corrected.csv">All 1,833 individual NE01 cell values (earlier integer and current corrected, section/ROI/t-SNE coordinates)</a> · <a href="data/stage1049_NE01_VGAT_VGLUT2_by_section_spot_count_statistics.csv">Gene spot count summaries</a> · <a href="data/stage1049_NE01_joint_VGAT_VGLUT2_threshold_counts.csv">Paired threshold combinations</a> · <a href="data/stage1049_NE01_GABA_GLU_spot_count_audit.json">Source and analysis limitations</a>.</p>',
'</section>'
]
chunk='\n'.join(parts)+'\n'
target='<section id="stage1048-ne-marker-correction">'
assert target in text
text=text.replace(target,chunk+target,1)
text=text.replace('<nav><a href="#stage1048-ne-marker-correction">','<nav><a href="#stage1049-ne01-cell-spots">New: NE cluster VGAT/Glu counts</a><a href="#stage1048-ne-marker-correction">',1)
h.write_text(text,encoding='utf8')
soup=BeautifulSoup(text,'html.parser')
section=soup.find(id=secid)
broken=[]
for a in section.find_all(['img','a']):
 href=a.get('src') or a.get('href')
 if href and not href.startswith(('http://','https://','#')) and not (P/href).is_file():
  broken.append(href)
assert not broken,broken
print('STAGE1049_PUBLISHED_LOCAL',len(section.find_all('img')),'FIGURES',len(section.find_all('a')),'LINKS','NO_MISSING')
