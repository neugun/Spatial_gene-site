"""Stage961: resolve false region fragmentation warnings without changing expression or biological labels."""
from pathlib import Path
from bs4 import BeautifulSoup
import json
root=Path(r"G:\Spatial_gene_site_publish")
pub=root/"perilc-map6-review"
data=pub/"data"
r=json.loads((data/"stage960_region10_topology_audit.json").read_text(encoding="utf8"))
old=json.loads((data/"stage959_atlas_integrity_audit.json").read_text(encoding="utf8"))
assert old["critical_issue_count"]==0
results={}
for ss in (500,530,560):
    a=r["sections"][str(ss)]
    islands=a["islands"]
    is_mask_island=(a["mask_components"] == len(islands)+1
      and a["component_count"] == len(islands)+10
      and all((x["mask_component"]>1 and x["neighbor_region"] is None)
              for x in islands)
      and a["small_island_fraction_of_tissue"]<0.001)
    assert is_mask_island, f"S{ss}: topology requires manual review; do not force relabel"
    results[str(ss)]={"state":"MAIN_TISSUE_10_CONNECTED_REGIONS",
        "mask_components":a["mask_components"],"minor_mask_pixels":sum(x["pixels"] for x in islands),
        "minor_mask_fraction":a["small_island_fraction_of_tissue"],
        "n_small_detached_mask_patches":len(islands),
        "original_region_labels_preserved":True}
summary={"stage":961,"status":"PASS_TOPOLOGY_QUALIFIED",
    "meaning":"Stage959 disconnected components are tiny disconnected components of the density mask; not splits of regions within the principal tissue mask.",
    "note":"Cannot infer whether detached mask patches are real tissue or low-density support; no region/cell reassignment.",
    "expression_authority":"Stage943 Route A","region_authority":"Stage956 unchanged",
    "sections":results,"image":"assets/stage960_region_topology/REGION10_TOPOLOGY_QA.png"}
(data/"stage961_region10_topology_resolution.json").write_text(json.dumps(summary,indent=2),encoding="utf8")
p=pub/"index.html"
s=BeautifulSoup(p.read_text(encoding="utf8"),"html.parser")
sec=s.find("section",id="regions")
assert sec is not None
for el in sec.select(".topology-audit-note"):el.decompose()
note=s.new_tag("p",attrs={"class":"lead topology-audit-note","id":"region-topology-qa"})
note.append("Topology QA (Stage961): Each main tissue mask contains 10 connected marker-guided regions. The 12 total components reported for S530/S560 arise from two tiny, detached density-mask patches per section (S530: 2 + 13 pixels; S560: 5 + 7 pixels), each under 0.1% of its section mask. Original expression and cell-region assignments are unchanged. ")
a=s.new_tag("a",href="assets/stage960_region_topology/REGION10_TOPOLOGY_QA.png",target="_blank")
a.string="View three-section topology QA."
note.append(a)
notes=sec.find("div",class_="region-notes")
assert notes is not None
notes.insert_after(note)
p.write_text(str(s),encoding="utf8")
print(json.dumps(summary,indent=2))
