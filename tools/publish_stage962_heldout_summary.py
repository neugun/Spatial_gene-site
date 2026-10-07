from pathlib import Path
from bs4 import BeautifulSoup
p=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review\index.html")
html=BeautifulSoup(p.read_text(encoding="utf8"),"html.parser")
sec=html.find("section",id="downstream")
assert sec
for e in sec.select(".stage962-summary"):e.decompose()
lead=sec.find("p",class_="lead")
note=html.new_tag("p",attrs={"class":"lead stage962-summary"})
note.append("Non-fitting-gene check (Stage962): 15 genes per section were not used in the region fitting. Their median within-Fine26 regional expression eta-squared is S500 0.0625, S530 0.0846, S560 0.0639. These values are descriptive within the same tissue, not an independent animal replication or a p-value. ")
a=html.new_tag("a",href="data/stage962_fit_marker_vs_heldout.csv")
a.string="Per-gene audit table"
note.append(a)
lead.insert_after(note)
p.write_text(str(html),encoding="utf8")
print("UPDATED",p)
