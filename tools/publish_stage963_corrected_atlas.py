from pathlib import Path
import re
p=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
f=p/"index.html"
s=f.read_text(encoding="utf8")
old="assets/stage947_final_atlas/FINAL_LC_PERILC_LOCATOR.png"
new="assets/stage963_orientation_smooth/FINAL_LC_PERILC_LOCATOR_XY_REVERSED.png"
assert s.count(old)==2
s=s.replace(old,new)
old="assets/stage956_final_atlas_v2/FINAL_MARKER_GUIDED_REGION10.png"
new="assets/stage963_orientation_smooth/FINAL_REGION10_SMOOTH_XY_REVERSED.png"
assert s.count(old)==2
s=s.replace(old,new)
pattern=r'assets/stage956_final_atlas_v2/([A-Za-z0-9]+_S(?:500|530|560)_expression_preview\.jpg)'
assert len(re.findall(pattern,s))==162
s=re.sub(pattern,r'assets/stage963_orientation_smooth/\1',s)
assert "Stage961" in s and "Stage962" in s
s=s.replace("S560 uses the stronger smoothing search selected by the marker-preservation gate.","S560 uses the stronger marker-field smoothing search selected by the marker-preservation gate. All displayed region outlines receive an additional sigma=3.0 grid-unit contour smoothing (Stage963), without changing any cell-level region identity.")
s=s.replace("All maps use reversed x/y display orientation.","All maps and the anatomical LC/periLC locator use exactly the same reversed X/Y display orientation; raw source coordinates and regional cell labels remain unchanged.")
f.write_text(s,encoding="utf8")
# Patch only legacy generator references so it cannot silently reintroduce mismatched axes.
b=p.parent/"tools"/"build_stage957_pages.py"
bs=b.read_text(encoding="utf8")
bs=bs.replace('MAN=pd.read_csv(DATA/"stage956_gene_section_manifest.csv")','MAN=pd.read_csv(DATA/"stage963_gene_section_manifest.csv")')
bs=bs.replace('assets/stage956_final_atlas_v2/{r.preview}','assets/stage963_orientation_smooth/{r.file}')
bs=bs.replace("assets/stage947_final_atlas/FINAL_LC_PERILC_LOCATOR.png","assets/stage963_orientation_smooth/FINAL_LC_PERILC_LOCATOR_XY_REVERSED.png")
bs=bs.replace("assets/stage956_final_atlas_v2/FINAL_MARKER_GUIDED_REGION10.png","assets/stage963_orientation_smooth/FINAL_REGION10_SMOOTH_XY_REVERSED.png")
b.write_text(bs,encoding="utf8")
print("STAGE963_PAGE_PATCHED 162 expression image refs + 2 locator + 2 region")
