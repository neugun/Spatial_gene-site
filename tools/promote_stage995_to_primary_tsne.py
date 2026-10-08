from pathlib import Path
base=Path(r"G:\Spatial_gene_site_publish\tools")
newroot=r"G:\Map6_recover_all\stage995_tsne\STAGE995_CURRENT_LOGZ_PCA20_TSNE_30_100.npz"
for name in ["build_stage989_tsne_full_suite.py","build_stage989_all27_tsne_montage.py","build_stage988c_shared_percentile_selection.py"]:
 p=base/name;s=p.read_text(encoding="utf8")
 old=r"G:\Map6_recover_all\stage982_tsne\multiscale_30_300.npz"
 assert old in s,(name,"missing prior")
 s=s.replace(old,newroot)
 if name=="build_stage989_tsne_full_suite.py":
  s=s.replace('"input":"Stage704 PCA15, Stage982 unsupervised multiscale 30/300 t-SNE"','"input":"Stage995 log-count gene-wise StandardScaler PCA20, unsupervised multiscale 30/100 t-SNE"')
 p.write_text(s,encoding="utf8")
 print("PROMOTED_STAGE995_COORDINATES",name,flush=True)
