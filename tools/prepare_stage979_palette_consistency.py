"""Re-use passed Stage974 XY geometry; update Fine26 colors only for Stage979."""
from pathlib import Path
tool=Path(r"G:\Spatial_gene_site_publish\tools")
base=(tool/"build_stage974_maskbody_preflipped.py").read_text(encoding="utf8")
assert base.count('color=cmap(c/25)')==2
base=base.replace("stage974_preflipped_reference","stage979_consistent_fine26")
base=base.replace("stage974_local_maskbody_summary.csv","stage979_local_maskbody_summary.csv")
base=base.replace("stage974_maskbody_authority.json","stage979_maskbody_authority.json")
base=base.replace("stage974_maskbody_coordinate_audit.json","stage979_maskbody_coordinate_audit.json")
base=base.replace("stage=974","stage=979")
base=base.replace('cmap=plt.get_cmap("turbo")', '''color_key=pd.read_csv(DATA/"stage977_fine26_color_key.csv")
num_to_stable=A[["fine26","fine26_stable_id"]].drop_duplicates().set_index("fine26").fine26_stable_id.to_dict()
pal=dict(zip(color_key.fine26,color_key.color))
palette_by_num={int(num):pal[stable] for num,stable in num_to_stable.items()}
assert len(palette_by_num)==26''')
base=base.replace("color=cmap(c/25)","color=palette_by_num[int(c)]")
base=base.replace("# atlas y_um is non-flipped; match ROI identities rather than coordinates for selection","# Stage631 XY are both already flipped; select local ROIs using the same cell-index authority")
(tool/"build_stage979_maskbody_colors.py").write_text(base,encoding="utf8")
s=(tool/"build_stage974_fine26_preflipped.py").read_text(encoding="utf8")
assert 'stage833_fine26_color_key.csv' in s
s=s.replace('stage833_fine26_color_key.csv','stage977_fine26_color_key.csv')
s=s.replace('stage974_preflipped_reference','stage979_consistent_fine26')
(tool/"build_stage979_spatial_colors.py").write_text(s,encoding="utf8")
print("STAGE979_RENDER_CODE_CREATED, identical XY methods, palette only",flush=True)
