from pathlib import Path
import json,re,pandas as pd
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
F=P/"single-cell.html";h=F.read_text(encoding="utf8")
D=P/"data";rows=pd.read_csv(D/"stage1037_native_large_fov_footprint_qc.csv")
def card(title,img,desc):
 assert (P/img).is_file(),img
 return '<article class="card"><h3>'+title+'</h3><a href="'+img+'"><img loading="lazy" src="'+img+'" alt="'+title+'"></a><p>'+desc+'</p></article>'
wide="assets/stage1036_fullXY_spot_context/"
native="assets/stage1037_large_native_footprints/"
cards=[]
for gene in ("Snap25","VGAT","VGLUT2"):
 r=wide+f"STAGE1036_S500_{gene}_FULLXY_TIFFz175_RAW_HP_SPOTOVERLAY.png"
 links=' · '.join(f'<a href="{wide}STAGE1036_S500_{gene}_FULLXY_TIFFz{z:03d}_RAW_HP_SPOTOVERLAY.png">Z{z}</a>' for z in (150,175,200))
 cards.append(card(f"{gene} · COMPLETE section XY raw/high-pass/RS-FISH detection",
 r,f"Whole S500 field 1119×~1116 pixels; 7-plane TIFF z-slab. All official detected spots shown (rings deliberately enlarged for full-field legibility, NOT real spot footprint). Alternative full-plane slices: {links}. This overview TIFF is downsampled 8× in XY and 4× in Z versus actual RS-FISH N5."))
for region in (1,2,3):
 for gene in ("Snap25","VGAT","VGLUT2"):
  t=rows[(rows.region==region)&(rows.gene==gene)].iloc[0]
  img=native+t.filename
  cards.append(card(f"Region {region} · {gene} · native N5 image-derived spot footprints",
 img,f"138 × 138 µm at original 0.23 µm XY pixel size; 15 native Z planes (5.9µm). "
 f"<b>{t.n_spots_in_cropped_slab:,} official RS-FISH detections</b>; {t.n_bright_image_supported_footprints:,} "
 f"have connected locally bright image-supported footprints ({t.pct_with_image_footprint:.1f}%). "
 f"Median estimated footprint radius {t.median_estimated_support_radius_native_px:.1f} native pixels. "
 f"White center = listed spot centroid, colored outlined area = local RAW image-derived connected punctum support, not a reported RS-FISH PSF. "
 f"Z center={t.source_native_Z_center:.0f} native voxels, selected within its own acquisition round to show puncta clearly."))
section='<section id="stage1036-full-raw-spots"><h2>Complete XY fields and true raw-N5 punctum footprints (Snap25, VGAT, VGLUT2)</h2>'
section+='<p class="lead">Following review of the prior six small-window overlays, we rebuilt each gene from its own original fluorescence volume and official RS-FISH spots. <b>Whole XY:</b> all ~1119×1116 TIFF-preview pixels, three Z windows per gene (150, 175, 200), with RAW, local background high-pass, and ALL spot overlays. This is the full acquisition XY and therefore avoids hidden arbitrary cropped fields. Because TIFF is 8× XY /4× Z downsampled, full-field rings are display markers, not footprints.</p>'
section+='<p class="lead"><b>Actually expanded punctum support:</b> nine additional high-resolution native N5 panels at ~138×138µm (3 fields × 3 genes), each with raw, high-pass and detected spots outlined over its local bright connected footprint. The source N5 resolution is 0.23µm XY, 0.42µm Z; only the visualization applies contrast enhancement. At region 1, VGAT has 58 listed detections and 51 native image-supported punctum footprints; VGLUT2 has 221 detections and 182 image-supported footprints. The previous small ROI plotting clearly underrepresented these detections.</p>'
section+='<p class="lead"><b>Scientific restrictions:</b> brightness-derived connected footprints are new VISUAL QC representations built around official centroid detections; they are not physical PSF estimates, new RS-FISH calls, adjusted count values, or a validated per-cell spot mask. Each gene comes from a different acquisition round (R1, R77, R99), and its Z slab is chosen in that round to contain native detections: the displayed columns are <b>not automatically registered to a common R1 3D cell mask</b>. This cannot yet distinguish true dual-transmitter identity from adjacent-cell spot assignment. Separate native spot sharpness measurements compare the potential VGAT defocus directly.</p>'
section+='<div class="cards">'+''.join(cards)+'</div>'
section+='<p class="lead">Sources and count tables: <a href="data/stage1036_fullXY_spots_by_gene_and_Z.csv">full XY spot counts per gene/Z</a> · <a href="data/stage1037_native_large_fov_footprint_qc.csv">9-panel image support, spot radii and Z centers</a> · <a href="data/stage1037_native_footprints_authority.json">native image footprint extraction specifications and limits</a>. Original N5 image chunks and single-cell coordinates remain in private workspace.</p></section>'
if 'id="stage1036-full-raw-spots"' not in h:
 assert '<section id="stage1033-native-rsfish-overlay">' in h
 h=h.replace('<section id="stage1033-native-rsfish-overlay">',section+'<section id="stage1033-native-rsfish-overlay">',1)
 h=h.replace('<a href="#stage1033-native-rsfish-overlay">','<a href="#stage1036-full-raw-spots">Full XY + real spot footprints</a><a href="#stage1033-native-rsfish-overlay">',1)
F.write_text(h,encoding="utf8")
print("STAGE1039_PUBLICATION_READY",len(cards),"new image cards",len(h),flush=True)
