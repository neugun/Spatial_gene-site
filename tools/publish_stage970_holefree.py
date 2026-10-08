"""Publish hole-free x flipped / y flipped spatial atlas with matching composites."""
from pathlib import Path
from bs4 import BeautifulSoup
import json
pub=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
p=pub/"index.html"
s=p.read_text(encoding="utf8")
from_old="assets/stage963_orientation_smooth/"
old_reg=from_old+"FINAL_REGION10_SMOOTH_X_FLIPPED_Y_FLIPPED.png"
new_reg="assets/stage970_holefree_xy/FINAL_REGION10_HOLEFREE_X_FLIPPED_Y_FLIPPED.png"
assert s.count(old_reg)==2
s=s.replace(old_reg,new_reg)
old_count=s.count(from_old)
# Preserve raw-coordinate-derived LC locator/Fine26, update only 81 gene panels.
import re
pattern=r"assets/stage963_orientation_smooth/([A-Za-z0-9]+_S(?:500|530|560)_expression_preview\.jpg)"
matches=re.findall(pattern,s)
assert len(matches)==162, len(matches)
s=re.sub(pattern,r"assets/stage970_holefree_xy/\1",s)
assert s.count(from_old)==old_count-162
s=s.replace("All displayed region outlines receive an additional sigma=3.0 grid-unit contour smoothing (Stage963), without changing any cell-level region identity.",
 "Stage970 refines region display topology from the marker-derived section grids: all enclosed holes are filled, minor disconnected label islands merged into a neighboring region, and boundaries smoothed. Every section has ten connected, hole-free domains; original per-cell region assignments and Stage943 expression remain unchanged.")
anchor=s.find("<section id=\"genes\">")
assert anchor>=0
s=s[:anchor]+"""<p class="lead" id="stage970-display-qa">Region display QA: S500 / S530 / S560 each has 10 connected, hole-free regions. Relative to the Stage956 reference, only 1.15%, 1.12%, and 1.68% of originally assigned tissue-grid pixels changed, respectively; no cell labels were changed. <a href="data/stage970_holefree_region_audit.json">Full geometry audit</a>.</p>
"""+s[anchor:]
p.write_text(s,encoding="utf8")
print("STAGE970_LINKS",len(matches),"REGION_MAIN",2)
