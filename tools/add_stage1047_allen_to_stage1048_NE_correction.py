from pathlib import Path
import json
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review')
h=P/'single-cell.html';s=h.read_text(encoding='utf8')
start=s.index('<section id="stage1048-ne-marker-correction">')
end=s.index('</section>',start)
section=s[start:end]
if 'STAGE1047_ALLEN_PONS_DBH_TH_INDEPENDENT_NE_VS_SLC6A2' not in section:
 mark='<p class="lead"><b>What follows:</b>'
 assert mark in section
 source=json.loads((P/'data'/'stage1047_allen_NE_independent_summary.json').read_text())
 b='<article class="card"><h3>Allen Pons direct independent DBH+TH cross-check (raw snRNA-seq)</h3>'
 b+='<a href="assets/stage1047_slc6a2_NE_audit/STAGE1047_ALLEN_PONS_DBH_TH_INDEPENDENT_NE_VS_SLC6A2.png">'
 b+='<img loading="lazy" src="assets/stage1047_slc6a2_NE_audit/STAGE1047_ALLEN_PONS_DBH_TH_INDEPENDENT_NE_VS_SLC6A2.png"></a>'
 b+=('<p>Allen Pons WMB 10xv3 whole dataset 143,661 cells: independent Dbh+Th co-detected in '
 f'{source["NE_DbhANDTh_n"]} cells; Slc6a2 detected in '
 f'{source["NE_DbhANDTh_NET_detected"]}/{source["NE_DbhANDTh_n"]} ({source["NET_detect_in_DbhTh_pct"]:.1f}%) versus '
 f'{source["NET_detect_without_joint_DbhTh_pct"]:.2f}% of other cells. In the 20 triple-marker-positive examples, '
 f'Slc17a6 is detected in {source["VGLUT2_detect_with_DbhTh_and_NET_pct"]:.0f}% and Slc32a1 in '
 f'{source["VGAT_detect_with_DbhTh_and_NET_pct"]:.0f}% <em>(n=20, preliminary)</em>. ')
 b+='This supports the biological relevance of Slc6a2 as a marker, but shows that its absolute FISH count cutoff cannot be obtained from different-platform raw UMI, and NE/VGLUT2 biological overlap must not be artificially forced to zero. Dbh+Th in single-cell UMI has severe dropout and is not all LC NE. '
 b+='<a href="data/stage1047_allen_pons_NE_multigene_vs_Slc6a2_UMI.csv">All thresholds and negative comparisons</a></p></article>'
 extra='<div class="cards">'+b+'</div>'
 section=section.replace(mark,extra+mark,1)
 s=s[:start]+section+s[end:]
 h.write_text(s,encoding='utf8')
print('ALLEN_INDEPENDENT_APPEND_READY')
