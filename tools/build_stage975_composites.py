"""Stage971: 27 three-section composites assembled from the authoritative flipped panels."""
from pathlib import Path
from PIL import Image
import pandas as pd
pub=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
src=pub/"assets"/"stage973_correct_orientation_roundedge"
out=pub/"assets"/"stage975_all3_correct_orientation"
out.mkdir(exist_ok=True)
man=pd.read_csv(pub/"data"/"stage973_gene_section_manifest.csv")
rows=[]
for gene in man.gene.drop_duplicates():
 panels=[]
 for sec in (500,530,560):
  row=man[(man.gene==gene)&(man.section==sec)].iloc[0]
  assert row.x_flipped and row.y_flipped
  im=Image.open(src/row.file).convert("RGB")
  panels.append(im)
 h=max(im.height for im in panels);w=sum(im.width for im in panels)
 comp=Image.new("RGB",(w,h),"white")
 x=0
 for im in panels:
  comp.paste(im,(x,0));x+=im.width
 name=f"{gene}_S500_S530_S560_x_flipped_y_flipped.jpg"
 comp.save(out/name,quality=88,subsampling=0,optimize=True)
 rows.append(dict(gene=gene,file=name,sections="500;530;560",x_flipped=True,y_flipped=True,
                  shared_colorbar_vmax=float(man[man.gene==gene].vmax.iloc[0])))
pd.DataFrame(rows).to_csv(pub/"data"/"stage975_all3_correct_orientation_manifest.csv",index=False)
print("STAGE975_PASSED",len(rows),"composites",flush=True)
