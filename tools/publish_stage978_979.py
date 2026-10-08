"""Publish Stage978 comparison and unify Stage979 Fine26 palette across spatial/3D."""
from pathlib import Path
import json
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
h=P/"single-cell.html"
s=h.read_text(encoding="utf8")
new="stage979_consistent_fine26/"
old="stage974_preflipped_reference/"
for f in ["FINE26_STAGE631_PREFLIPPED.jpg","FINE26_TRUE_MASKBODY_3D.png","S500_LOCAL_TRUE_MASKBODY.png","S530_LOCAL_TRUE_MASKBODY.png","S560_LOCAL_TRUE_MASKBODY.png"]:
    assert (P/"assets"/new/f).is_file(),f
    assert (old+f) in s,f
    s=s.replace(old+f,new+f)
metrics=json.loads((P/"data"/"stage978_geometry_audit.json").read_text())
assert metrics["status"]=="BUILT_AWAITING_VISUAL_REVIEW"
pic="assets/stage978_umap_geometry/FINE26_STAGE704_VS_978.png"
assert (P/pic).is_file()
if 'stage978_umap_geometry/' not in s:
    note='''<article class="card"><h3>Fine26 UMAP geometry · old versus new</h3>
<a href="assets/stage978_umap_geometry/FINE26_STAGE704_VS_978.png"><img loading="lazy" src="assets/stage978_umap_geometry/FINE26_STAGE704_VS_978.png" alt="Same frozen Fine26 colors and cell labels, two molecular embedding parameterizations"></a>
<p>Stage978 alternative (cosine / 55 neighbors / min distance 0.09) modestly improves descriptive 15-neighbor Fine26 purity from 0.176 to 0.193, but the negative silhouettes (−0.277 / −0.269) show considerable molecular overlap. Both remain visualization candidates, not evidence for distinct new cell types.</p></article>
'''
    s=s.replace('<article class="card"><h3>Broad-class control on the UMAP</h3>',note+'<article class="card"><h3>Broad-class control on the UMAP</h3>',1)
s=s.replace('Fine26 spatial identities · corrected XY','Fine26 spatial identities · matching 26-color UMAP palette')
s=s.replace('Fine26 true-maskbody 3D atlas','Fine26 true-maskbody 3D atlas · same 26-color palette')
h.write_text(s,encoding="utf8")
print("STAGE978_979_PAGE_LINKS_UPDATED",s.count("stage979_consistent_fine26/"),"files",flush=True)
