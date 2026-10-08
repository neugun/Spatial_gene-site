from pathlib import Path
import pandas as pd
from bs4 import BeautifulSoup
pub=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
p=pub/"index.html"
s=BeautifulSoup(p.read_text(encoding="utf8"),"html.parser")
manifest=pd.read_csv(pub/"data"/"stage969_all3_flipped_manifest.csv").set_index("gene")
count=0
for card in s.select(".gene-card"):
 name=card.select_one("h3").get_text(strip=True)
 assert name in manifest.index
 t=card.select_one(".gene-title span")
 assert t
 for old in t.select(".allthree-link"):old.decompose()
 link=s.new_tag("a",href="assets/stage969_all3_flipped/"+manifest.loc[name,"file"],
  target="_blank",attrs={"class":"allthree-link"})
 link.string=" View three-section composite"
 t.append(" · ");t.append(link)
 count+=1
assert count==27
details=s.find("details")
assert details
link=s.new_tag("a",href="data/stage969_all3_flipped_manifest.csv")
link.string="27 x flipped, y flipped composites manifest"
details.append(" · ");details.append(link)
p.write_text(str(s),encoding="utf8")
print("STAGE969_PAGE_LINKS",count)
