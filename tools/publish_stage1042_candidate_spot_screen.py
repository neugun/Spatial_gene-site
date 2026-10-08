from pathlib import Path
import pandas as pd
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
F=P/"single-cell.html";h=F.read_text(encoding="utf8")
D=P/"data";T=pd.read_csv(D/"stage1041_source_image_candidate_peaks_vs_RSFISH.csv")
t=T[T.candidate_threshold_full_volume_quantile==.9975].set_index("gene")
def card(title,path,desc):
 assert (P/path).is_file(),path
 return f'<article class="card"><h3>{title}</h3><a href="{path}"><img loading="lazy" src="{path}"></a><p>{desc}</p></article>'
root="assets/stage1041_peak_recall_screen/"
cards=[]
for gene in ("VGAT","VGLUT2","Snap25"):
 v=t.loc[gene]
 cards.append(card(f"{gene}: independently bright RAW image puncta vs RS-FISH",
 root+f"STAGE1041_REGION1_{gene}_BRIGHT_UNMATCHED_OVERLAY.png",
 f"Source N5 native 3D region 1, 15 original Z planes, spot-independent high-pass local maxima threshold 99.75th percentile. "
 f"<b>{int(v.n_raw_peaks_near_RSFISH)} of {int(v.n_3D_raw_bright_peaks)} image candidates</b> match a nearby official detection (≤4.5 native-equivalent pixels). "
 "Green circles are matched image candidates; pink X marks are other bright maxima, NOT certified missed RNA. Raw image and 3D high-pass alongside."))
cards.append(card("Across thresholds: image-derived bright maxima supported by actual RS-FISH detections",
 root+"STAGE1041_RAW_IMAGE_PEAKS_VS_RSFISH.png",
 "Source N5 intensity 99.5%, 99.75% and 99.9% candidate cutoffs; relative comparison across three rounds, not calibrated detector sensitivity or transcript recall."))
html='<section id="stage1041-candidate-recall"><h2>Possible VGAT under-detection: image-derived 3D bright-puncta candidates not assigned to RS-FISH</h2>'
html+='<p class="lead"><b>More important than successful detection centroids is whether RS-FISH missed other real image puncta.</b> In the same 138×138µm S500 region 1 used above, we independently detected strong 3D local high-pass maxima in native raw source images, then checked whether a listed RS-FISH centroid was nearby. At the predeclared 99.75th percentile of 3D local filtered intensity: <b>VGAT 39/86 (45.3%)</b> bright candidate peaks had a nearby existing spot; <b>VGLUT2 128/159 (80.5%)</b>; <b>Snap25 184/193 (95.3%)</b>. At the stricter 99.9th percentile, VGAT is 39/61 (63.9%), VGLUT2 89/98 (90.8%), Snap25 101/105 (96.2%).</p>'
html+='<p class="lead">This identifies a <b>meaningful candidate VGAT spot-detection recall issue</b> worth direct image review: many bright VGAT source maxima are not paired to a formal RS-FISH centroid. However, the high-pass candidate detector is not ground truth; unmatched peaks may be noise, autofluorescence, multi-Z overlap, or duplicate peaks. It would be wrong to call the 47 unmatched VGAT candidates a measured number of undetected transcripts. Our earlier random QC mainly tested precision/spot centering among detected spots, not detection recall, so these findings are compatible.</p>'
html+='<div class="cards">'+''.join(cards)+'</div><p class="lead"><a href="data/stage1041_source_image_candidate_peaks_vs_RSFISH.csv">Threshold/count comparison table</a> · <a href="data/stage1041_candidate_recall_authority.json">3D candidate definitions and uncertainty</a>. This one-area result remains independent of same-cell cross-round segmentation validation.</p></section>'
if 'stage1041-candidate-recall' not in h:
 anchor='<section id="stage1038-vgat-focus">'
 assert anchor in h
 h=h.replace(anchor,html+anchor,1)
F.write_text(h,encoding="utf8")
print("STAGE1042_WEB_READY",len(cards),len(h),flush=True)
