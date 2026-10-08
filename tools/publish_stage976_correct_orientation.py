from pathlib import Path
import re
pub=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
index=pub/"index.html"
s=index.read_text(encoding="utf8")
# gene panel references, all 81 duplicate image+link
src="assets/stage970_holefree_xy/"
matches=re.findall(r"assets/stage970_holefree_xy/[A-Za-z0-9]+_S(?:500|530|560)_expression_preview\.jpg",s)
assert len(matches)==162,len(matches)
s=re.sub(r"assets/stage970_holefree_xy/([A-Za-z0-9]+_S(?:500|530|560)_expression_preview\.jpg)",
         r"assets/stage973_correct_orientation_roundedge/\1",s)
old_reg="assets/stage970_holefree_xy/FINAL_REGION10_HOLEFREE_X_FLIPPED_Y_FLIPPED.png"
new_reg="assets/stage973_correct_orientation_roundedge/FINAL_REGION10_STAGE631_PREFLIPPED_SMOOTH_OUTER.png"
assert s.count(old_reg)>=2
s=s.replace(old_reg,new_reg)
s=s.replace("assets/stage963_orientation_smooth/FINAL_LC_PERILC_LOCATOR_X_FLIPPED_Y_FLIPPED.png",
           "assets/stage974_preflipped_reference/FINAL_LC_PERILC_STAGE631_PREFLIPPED.png")
s=s.replace("assets/stage971_all3_holefree_xy/","assets/stage975_all3_correct_orientation/")
s=s.replace("data/stage971_all3_holefree_manifest.csv","data/stage975_all3_correct_orientation_manifest.csv")
s=s.replace("assets/stage967_maskbody_xy/","assets/stage974_preflipped_reference/")
s=s.replace("data/stage970_holefree_region_audit.json","data/stage973_correct_orientation_region_audit.json")
s=s.replace("Stage970 refines region display topology","Stage973 uses already x-flipped and y-flipped Stage631 coordinates (no second flip) and refines region display topology")
s=s.replace("All maps and the anatomical LC/periLC locator use exactly the same x flipped, y flipped display orientation; raw source coordinates and regional cell labels remain unchanged.",
"""All current maps use the XY-flipped Stage631 physical coordinates directly, the same as historical Stage838; an additional flip would undo the desired orientation. Expression/count and per-cell region identity are unchanged.""")
s=s.replace("Geometry was rebuilt from section-specific, x-flipped and y-flipped cell coordinates and the marker-derived region grid, not by flipping an existing image.",
"""The Stage631 x/y coordinates were already flipped once relative to the acquisition array. Stage973 removes the previously incorrect second flip, and smooths both outer anatomy and internal boundaries. All panels are regenerated from cell coordinates and numerical region assignments, not mirrored raster pictures.""")
s=s.replace("Previous region shape","Prior double-flip shape (superseded)")
s=s.replace("Updated hole-free shape","Current original-orientation smooth-outline shape")
s=s.replace("assets/stage963_orientation_smooth/FINAL_REGION10_SMOOTH_X_FLIPPED_Y_FLIPPED.png","assets/stage970_holefree_xy/FINAL_REGION10_HOLEFREE_X_FLIPPED_Y_FLIPPED.png")
index.write_text(s,encoding="utf8")
sc=pub/"single-cell.html"
t=sc.read_text(encoding="utf8")
t=t.replace("assets/stage967_maskbody_xy/","assets/stage974_preflipped_reference/")
t=t.replace("assets/stage963_orientation_smooth/FINE26_XY_CURRENT.jpg","assets/stage974_preflipped_reference/FINE26_STAGE631_PREFLIPPED.jpg")
t=t.replace("assets/stage963_orientation_smooth/BROAD_XY_CURRENT.jpg","assets/stage974_preflipped_reference/BROAD_STAGE631_PREFLIPPED.jpg")
for sec in (500,530,560):
 t=t.replace(f"assets/stage963_orientation_smooth/S{sec}_LOCAL_MULTIPLEX_XY.png",
             f"assets/stage974_preflipped_reference/S{sec}_LOCAL_MULTIPLEX_STAGE631_PREFLIPPED.png")
t=t.replace("These three-section Fine26/broad spatial maps show selected medial periLC cells (not full-section class coverage). The XY display follows the same coordinate transform as the main gene atlas.",
"""These spatial panels use the already-flipped Stage631 physical coordinates (no double-flip). Fine26/Broad maps show selected medial periLC cells, not full-section class coverage.""")
t=t.replace("Stage967 displays Stage951 true segmentation bodies with x flipped, y flipped coordinates",
"""Stage974 displays Stage951 true segmentation bodies using preflipped Stage631 coordinates without an extra flip""")
sc.write_text(t,encoding="utf8")
print("STAGE976_PAGE_UPDATED",len(matches),"gene img+link",s.count("stage975_all3_correct_orientation/"),"composite links")
