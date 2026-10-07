from pathlib import Path
from bs4 import BeautifulSoup
p=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review\single-cell.html")
s=BeautifulSoup(p.read_text(encoding="utf8"),"html.parser")
local=s.find("section",id="local")
assert local
for k in (500,530,560):
 old=f"assets/stage950_singlecell/S{k}_LOCAL_SINGLECELL_MULTIPLEX.png"
 new=f"assets/stage963_orientation_smooth/S{k}_LOCAL_MULTIPLEX_XY.png"
 for e in local.find_all(["a","img"]):
  for att in ("href","src"):
   if e.get(att)==old:e[att]=new
for e in local.select(".stage964-spatial"):e.decompose()
extra=s.new_tag("div",attrs={"class":"cards stage964-spatial"})
for title,file in (("Fine26 spatial identities · corrected XY","FINE26_XY_CURRENT.jpg"),("Broad molecular classes · corrected XY","BROAD_XY_CURRENT.jpg")):
 art=s.new_tag("article",attrs={"class":"card"})
 h=s.new_tag("h3");h.string=title;art.append(h)
 href="assets/stage963_orientation_smooth/"+file
 a=s.new_tag("a",href=href,target="_blank")
 im=s.new_tag("img",src=href,alt=title)
 a.append(im);art.append(a);extra.append(art)
para=s.new_tag("p",attrs={"class":"lead stage964-spatial"})
para.string="These three-section Fine26/broad spatial maps show selected medial periLC cells (not full-section class coverage). The XY display follows the same coordinate transform as the main gene atlas. Local multiplex scatter panels below have also been regenerated from source coordinates."
local.find("p",class_="lead").insert_after(extra)
extra.insert_before(para)
p.write_text(str(s),encoding="utf8")
print("SINGLE_CELL_XY_PATCHED")
