from pathlib import Path
import json
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
q=P/"data"/"stage993_public_display_authority.json"
a=json.loads(q.read_text(encoding="utf8"))
a["primary_embedding"]="Stage995 corrected-log counts gene-wise standardized PCA20 -> unsupervised tSNE perplexities 30/100"
a["supervised_tsne"]=False
a["stage995_broad_and_fine_purity"]="Stage996 fixed sample Fine26 0.4349 Broad 0.8256, compare Stage982 0.2271/0.7229"
a["full_local_true_label_surfaces"]={"500":110,"530":131,"560":116}
a["stage996_metrics"]="stage996_tsne_preprocess_quality_comparison.json"
q.write_text(json.dumps(a,indent=2),encoding="utf8")
for sec,n in [(500,110),(530,131),(560,116)]:
 file=P/"data"/("stage992_true_mask_surface_authority.json" if sec==500 else f"stage992_true_mask_surface_authority_S{sec}.json")
 au=json.loads(file.read_text(encoding="utf8"))
 assert au["genuine_meshes"]==n,(sec,au["genuine_meshes"])
 assert "sternsonlab" not in file.read_text(encoding="utf8")
 plot=P/"assets"/"stage992_true_mask_surface"/f"STAGE992_S{sec}_ACTUAL_SEGMENTATION_SURFACE_VS_POINTCLOUD.png"
 assert plot.stat().st_size>20000
 print("VALID_TRUE_SEGMENTATION",sec,n,plot.stat().st_size,flush=True)
assert (P/"assets"/"stage995_tsne_preprocess"/"STAGE995_ORIGINAL_VS_LOGZ_TSNE.png").stat().st_size>100000
print("STAGE997_RELEASE_INTEGRITY_OK",flush=True)
