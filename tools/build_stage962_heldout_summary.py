from pathlib import Path
import pandas as pd, json
p=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
d=p/"data"
w=pd.read_csv(d/"stage952_within_fine26_spatial_modulation.csv")
m=pd.read_csv(d/"stage956_marker_guided_region10_manifest.csv").set_index("section")
assert len(w)==81 and w.groupby("section").size().eq(27).all()
w["used_in_region_fit"]=[g in str(m.loc[s,"all_markers"]).split(";") for s,g in zip(w.section,w.gene)]
w.to_csv(d/"stage962_fit_marker_vs_heldout.csv",index=False)
v=w.groupby(["section","used_in_region_fit"]).weighted_within_fine26_eta2.agg(["count","median","mean"]).reset_index()
out={"stage":962,"status":"DESCRIPTIVE_WITHIN_FINE26_HOLDOUT_AUDIT","results":v.to_dict("records"),"caveat":"Within-section non-fitting genes, same tissue: descriptive only, no independent mouse replication or p-values","expression_authority":"Stage943 Route A","region_authority":"Stage956","statistics_source":"Stage952"}
(d/"stage962_holdout_authority.json").write_text(json.dumps(out,indent=2),encoding="utf8")
print(json.dumps(out,indent=2))
