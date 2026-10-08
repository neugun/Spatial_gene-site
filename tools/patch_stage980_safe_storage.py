from pathlib import Path
root=Path(r"G:\Spatial_gene_site_publish\tools")
p=root/"build_stage978_umap_geometry.py"
s=p.read_text(encoding="utf8")
old='np.savez_compressed(out/"STAGE978_ALTERNATE_GEOMETRY.npz",embedding=new,labels=lab)'
new='private=Path(r"G:\\Map6_recover_all\\stage978_geometry");private.mkdir(parents=True,exist_ok=True)\nnp.savez_compressed(private/"STAGE978_ALTERNATE_GEOMETRY.npz",embedding=new,labels=lab)'
assert s.count(old)==1
p.write_text(s.replace(old,new),encoding="utf8")
p=root/"build_stage980_true3d_video.py";s=p.read_text(encoding="utf8")
needle='O=PUB/"assets"/"stage980_true3d_video";O.mkdir(parents=True,exist_ok=True)'
replacement=needle+'\nQA=Path(r"G:\\Map6_recover_all\\stage980_movie_QC");QA.mkdir(parents=True,exist_ok=True)'
assert needle in s
s=s.replace(needle,replacement).replace('name=O/f"QC_FRAME_{k:03d}.png"','name=QA/f"QC_FRAME_{k:03d}.png"')
p.write_text(s,encoding="utf8")
print("FUTURE_RAW_AND_QC_OUTPUTS_MOVED_OFF_PUBLIC",flush=True)
