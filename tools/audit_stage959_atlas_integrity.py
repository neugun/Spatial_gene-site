"""Current atlas integrity gate. PASS means display/data consistency, not biological validation."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[1] / "perilc-map6-review"
D = ROOT / "data"
AS = ROOT / "assets"
issues, warnings = [], []
def require(ok, msg):
    if not ok: issues.append(msg)

manifest = pd.read_csv(D / "stage956_gene_section_manifest.csv")
policy = pd.read_csv(D / "stage956_gene_colorbar_policy.csv")
regions = pd.read_csv(D / "stage956_cell_region10_labels.csv.gz")
summary = pd.read_csv(D / "stage956_marker_guided_region10_manifest.csv")
authority = json.loads((D / "stage956_final_atlas_authority.json").read_text(encoding="utf-8"))
single = json.loads((D / "stage958_singlecell_page_authority.json").read_text(encoding="utf-8"))
projection = json.loads((D / "stage954_projection_latent_authority.json").read_text(encoding="utf-8"))

genes = manifest.gene.drop_duplicates().tolist()
require(len(manifest) == 81, f"Expected 81 gene×section maps; got {len(manifest)}")
require(len(genes) == 27, f"Expected 27 genes; got {len(genes)}")
require(len(regions) == 71950, f"Expected 71,950 cell region assignments; got {len(regions)}")
require(manifest.section.isin([500,530,560]).all(), "Unexpected section in manifest")
require(manifest.groupby("gene").section.nunique().eq(3).all(), "A gene is missing a section")
require(manifest.x_flipped.astype(bool).all() and manifest.y_flipped.astype(bool).all(), "x flipped, y flipped flags invalid for all maps")
require(manifest.groupby("gene").vmax.nunique().eq(1).all(), "Unequal display vmax across sections within gene")
require(manifest.groupby("gene").gamma.nunique().eq(1).all(), "Unequal gamma across sections within gene")
require(np.isfinite(manifest.vmax).all() and (manifest.vmax > 0).all(), "Invalid colorbar maximum")
require(len(policy) == 27, "Colorbar policy incomplete")
snap = policy[policy.gene == "Snap25"]
require(len(snap) == 1 and abs(float(snap.vmax_quantile.iloc[0]) - .85)<1e-6, "Snap25 vmax must be 85th percentile")
require(authority["expression_authority"] == "Stage943 Route A", "Unrecognized expression authority")
require(single["spatial_authority"].startswith("Stage956"), "Single-cell bridge references stale spatial authority")
require(projection["status"] == "CIRCULARITY_AUDIT_NOT_FUNCTIONAL_EVIDENCE", "Projection proxy falsely promoted as functional evidence")
required={"Hcrtr1","Bcl11b","Slc5a7","Lmx1a","Piezo2","Slc6a2","Ghr"}
for _, row in summary.iterrows():
    ss=int(row.section)
    markers=set(str(row.all_markers).split(";"))
    require(required <= markers, f"S{ss} marker anchors incomplete")
    L=regions.loc[regions.section==ss, "region10"].astype(int)
    active=set(L.unique())-{0}
    require(active==set(range(1,11)), f"S{ss} not ten region identities: {sorted(active)}")
    require(float(row.retention) >= .88, f"S{ss} marker separation retention below 0.88")
    if int(row.components)>10:
        warnings.append(f"S{ss}: 10 labels but {int(row.components)} disconnected components; review islands/fragmentation")
    if (L==0).mean()>.05:
        warnings.append(f"S{ss}: unassigned region fraction {(L==0).mean():.2%}")
for _, row in manifest.iterrows():
    f=AS/"stage956_final_atlas_v2"/str(row.preview)
    require(f.exists() and f.stat().st_size>10000, f"Missing/empty map: {f.name}")

page_info={}
for fn in ("index.html","single-cell.html"):
    doc=BeautifulSoup((ROOT/fn).read_text(encoding="utf-8"),"html.parser")
    ids={e.get("id") for e in doc.select("[id]")}
    missing=[]
    for tag,attr in (("img","src"),("a","href")):
        for e in doc.find_all(tag):
            v=e.get(attr)
            if not v or v.startswith(("http://","https://","mailto:","data:")):continue
            if v.startswith("#"):
                if v[1:] and v[1:] not in ids:missing.append(v)
            else:
                target=ROOT/v.split("#")[0].split("?")[0]
                if not target.exists():missing.append(v)
    require(not missing, f"{fn}: {len(missing)} broken local refs: {missing[:5]}")
    page_info[fn]={"sections":len(doc.find_all("section")),"missing_links":len(missing)}
main=BeautifulSoup((ROOT/"index.html").read_text(encoding="utf-8"),"html.parser")
require(len(main.select(".gene-card"))==27, "Main page does not show 27 gene cards")
require(len(main.select(".gene-card .pane"))==81, "Main page does not show 81 gene×section panes")
result={
    "stage":959,
    "status":"PASS_WITH_REVIEW" if (not issues and warnings) else ("PASS" if not issues else "FAIL"),
    "critical_issue_count":len(issues),
    "critical_issues":issues,
    "review_warnings":warnings,
    "gene_count":len(genes),
    "gene_section_maps":len(manifest),
    "cell_count":len(regions),
    "section_regions":{str(int(s)):{'n_labels':int(regions.loc[regions.section==s,'region10'].replace(0,np.nan).nunique()),
                                   'connected_components':int(summary.loc[summary.section==s,'components'].iloc[0]),
                                   'unassigned_cells':int((regions.loc[regions.section==s,'region10']==0).sum())}
                       for s in (500,530,560)},
    "page_validation":page_info,
    "circularity_warning_preserved":True,
}
(D/"stage959_atlas_integrity_audit.json").write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
print(json.dumps(result,indent=2,ensure_ascii=False))
if issues: raise SystemExit(1)
