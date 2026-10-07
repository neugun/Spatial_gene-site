# Historical Stage957 page generator. It predates the Stage958–966 additions.
# Disabled by default so old page regeneration cannot silently erase newer validated content.
import sys
if "--allow-historical-stage957-overwrite" not in sys.argv:
    raise SystemExit("Historical Stage957 builder disabled; current site uses Stage958–966. Explicit --allow-historical-stage957-overwrite required.")
from pathlib import Path
import json, html
import pandas as pd

PUB=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
DATA=PUB/"data"
MAN=pd.read_csv(DATA/"stage963_gene_section_manifest.csv")
REG=pd.read_csv(DATA/"stage956_marker_guided_region10_manifest.csv")
POL=pd.read_csv(DATA/"stage956_gene_colorbar_policy.csv")
genes=MAN.gene.drop_duplicates().astype(str).tolist()
secs=[500,530,560]

range_map={r.gene:float(r.vmax) for _,r in POL.iterrows()}
class_map={r.gene:str(r.display_class) for _,r in POL.iterrows()}
cards=[]
for g in genes:
    panes=[]
    for ss in secs:
        r=MAN[(MAN.gene==g)&(MAN.section==ss)].iloc[0]
        panes.append(f'''<a class="pane" href="assets/stage963_orientation_smooth/{r.file}" target="_blank">
          <div class="pane-head"><span>S{ss}</span><span>0–{float(r.vmax):.1f}</span></div>
          <img src="assets/stage963_orientation_smooth/{r.file}" loading="lazy" alt="{html.escape(g)} S{ss} expression">
        </a>''')
    note="broad-expression scale" if class_map[g].startswith("broad") else "shared per-gene scale"
    cards.append(f'''<article class="gene-card" data-gene="{g.lower()}">
      <div class="gene-title"><h3>{html.escape(g)}</h3><span>{note} · vmax {range_map[g]:.1f}</span></div>
      <div class="gene-grid">{''.join(panes)}</div>
    </article>''')

region_notes=[]
for _,r in REG.iterrows():
    region_notes.append(f'''<div class="mini"><b>S{int(r.section)}</b><br>
      anchors + {html.escape(str(r.auto_specific_markers).replace(";", ", "))}<br>
      smoothing σ={float(r.smooth_sigma):.1f}</div>''')

css=r'''
*{box-sizing:border-box}
:root{--ink:#151515;--muted:#666;--line:#dedede;--soft:#f7f7f7;--max:1500px}
html{scroll-behavior:smooth}
body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif;color:var(--ink);background:#fff;line-height:1.45}
header{max-width:var(--max);margin:0 auto;padding:36px 24px 18px}
h1{font-size:34px;line-height:1.08;margin:0 0 10px;font-weight:720}
.subtitle{font-size:16px;color:var(--muted);margin:0}
nav{position:sticky;top:0;z-index:5;background:rgba(255,255,255,.96);backdrop-filter:blur(8px);border-top:1px solid var(--line);border-bottom:1px solid var(--line);padding:10px 24px}
nav .inner{max-width:var(--max);margin:auto;display:flex;gap:18px;flex-wrap:wrap}
nav a{color:#222;text-decoration:none;font-size:14px;font-weight:650}
main{max-width:var(--max);margin:auto;padding:12px 24px 70px}
section{padding:25px 0 36px;border-bottom:1px solid #ececec}
h2{font-size:25px;margin:0 0 8px}
.lead{color:var(--muted);max-width:1120px;margin:0 0 18px}
.hero{display:block;width:100%;max-width:1450px;height:auto;border:1px solid #eee;background:white}
.region-notes{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;margin:12px 0 0}
.mini{padding:10px 12px;background:#fafafa;border:1px solid #eee;border-radius:8px;color:#555;font-size:12px}
#geneSearch{width:min(520px,100%);font-size:16px;padding:11px 13px;border:1px solid #bbb;border-radius:8px;margin:4px 0 18px}
.gene-card{padding:18px 0 22px;border-top:1px solid #ececec}
.gene-title{display:flex;align-items:baseline;justify-content:space-between;gap:12px}
.gene-title h3{font-size:22px;margin:0 0 11px}.gene-title span{font-size:12px;color:#777}
.gene-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px}
.pane{display:block;text-decoration:none;color:#222;border:1px solid #e3e3e3;border-radius:8px;overflow:hidden;background:white}
.pane-head{display:flex;justify-content:space-between;padding:7px 10px 3px;font-size:12px;font-weight:700;color:#444}
.pane img{display:block;width:100%;height:auto}
.cards{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}
.card{border:1px solid #e5e5e5;border-radius:9px;overflow:hidden;background:white}
.card h3{font-size:17px;margin:12px 14px 5px}.card p{font-size:13px;color:#666;margin:0 14px 14px}
.card img{display:block;width:100%;height:auto}
.cta{display:inline-block;padding:9px 13px;border:1px solid #333;border-radius:7px;text-decoration:none;color:#111;font-weight:650;font-size:13px}
details{margin-top:22px;color:var(--muted);font-size:13px}details a{color:#444}
footer{max-width:var(--max);margin:auto;padding:0 24px 40px;color:#888;font-size:12px}
@media(max-width:760px){h1{font-size:28px}.gene-grid,.cards,.region-notes{grid-template-columns:1fr}main,header{padding-left:14px;padding-right:14px}nav{padding-left:14px;padding-right:14px}.gene-title{display:block}}
'''
js=r'''const q=document.getElementById('geneSearch');q.addEventListener('input',()=>{const s=q.value.trim().toLowerCase();document.querySelectorAll('.gene-card').forEach(c=>c.style.display=(!s||c.dataset.gene.includes(s))?'block':'none')});'''

index=f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Three-section spatial gene atlas</title><style>{css}</style></head><body>
<header><h1>Three-section spatial gene atlas</h1><p class="subtitle">S500 / S530 / S560 · Stage943 current expression · marker-guided section-local anatomy</p></header>
<nav><div class="inner"><a href="#locator">LC / periLC</a><a href="#regions">10 subregions</a><a href="#genes">Gene atlas</a><a href="#downstream">Downstream spatial</a><a href="single-cell.html">Single-cell bridge</a></div></nav>
<main>
<section id="locator"><h2>LC / periLC locator</h2>
<p class="lead">LC and the operational periLC field are shown only as anatomical locators. Gene-expression panels below show the complete section and are not cropped to these regions.</p>
<a href="assets/stage963_orientation_smooth/FINAL_LC_PERILC_LOCATOR_XY_REVERSED.png" target="_blank"><img class="hero" src="assets/stage963_orientation_smooth/FINAL_LC_PERILC_LOCATOR_XY_REVERSED.png" alt="LC periLC locator"></a></section>
<section id="regions"><h2>Marker-guided section-specific 10-subregion reference</h2>
<p class="lead">Each section is partitioned independently. The required anchors are <b>Hcrtr1, Bcl11b, Slc5a7, Lmx1a, Piezo2, Slc6a2 and Ghr</b>, supplemented by the most spatially informative genes in that section. Region numbers are local to each section; no cross-section homology is imposed. S560 uses the stronger smoothing search selected by the marker-preservation gate.</p>
<a href="assets/stage963_orientation_smooth/FINAL_REGION10_SMOOTH_XY_REVERSED.png" target="_blank"><img class="hero" src="assets/stage963_orientation_smooth/FINAL_REGION10_SMOOTH_XY_REVERSED.png" alt="marker guided section specific ten regions"></a>
<div class="region-notes">{''.join(region_notes)}</div></section>
<section id="genes"><h2>Final gene × section expression atlas</h2>
<p class="lead">All maps use reversed x/y display orientation. Every panel shows an explicit <b>spot-count colorbar</b> and numeric display range. The three sections share one range per gene. Broad genes use a lower saturation ceiling; Snap25 is intentionally shown with an 85th-percentile positive-count maximum so widespread expression remains visible.</p>
<input id="geneSearch" type="search" placeholder="Search gene…">{''.join(cards)}</section>
<section id="downstream"><h2>Downstream spatial analyses</h2>
<p class="lead">Only final biological summaries are retained here; correction sweeps and engineering intermediates remain outside the main atlas.</p>
<div class="cards">
<article class="card"><h3>Region10 × 27-gene molecular profiles</h3><a href="assets/stage956_final_atlas_v2/REGION10_GENE_PROFILE.png" target="_blank"><img src="assets/stage956_final_atlas_v2/REGION10_GENE_PROFILE.png" loading="lazy"></a><p>Section-local region identities are summarized independently; matching region numbers across sections are not assumed to be homologous.</p></article>
<article class="card"><h3>Fine26 populations in tissue space</h3><a href="assets/stage963_orientation_smooth/FINE26_XY_CURRENT.jpg" target="_blank"><img src="assets/stage963_orientation_smooth/FINE26_XY_CURRENT.jpg" loading="lazy"></a><p>Frozen single-cell molecular identities across the three-section atlas.</p></article>
<article class="card"><h3>Soma geometry by molecular population</h3><a href="assets/current_morph_size_distributions.jpg" target="_blank"><img src="assets/current_morph_size_distributions.jpg" loading="lazy"></a><p>Distribution-level morphology, not only effect-size summaries.</p></article>
<article class="card"><h3>Projection-associated molecular organization</h3><a href="assets/current_projection_reproducible.jpg" target="_blank"><img src="assets/current_projection_reproducible.jpg" loading="lazy"></a><p>Reproducible projection-linked enrichments are kept downstream of the direct molecular atlas.</p></article>
</div></section>
<section id="singlecell"><h2>Single-cell bridge</h2><p class="lead">The spatial atlas is linked to the frozen Fine26 single-cell state space without re-clustering for versioning. The dedicated page collects taxonomy, section reproducibility, morphology and projection analyses.</p>
<a class="cta" href="single-cell.html">Open current single-cell bridge →</a></section>
<details><summary>Data provenance</summary><p>
<a href="data/stage943_routeA_authority.json">Stage943 expression authority</a> ·
<a href="data/stage956_final_atlas_authority.json">Stage956 display authority</a> ·
<a href="data/stage956_gene_section_manifest.csv">81-map manifest</a> ·
<a href="data/stage956_marker_guided_region10_manifest.csv">region marker manifest</a> ·
<a href="data/stage956_gene_colorbar_policy.csv">colorbar policy</a>
</p></details>
</main><footer>Final biological atlas page. Intermediate correction/RS-FISH engineering results remain in repository provenance but are intentionally omitted here.</footer>
<script>{js}</script></body></html>'''
(PUB/"index.html").write_text(index,encoding="utf-8")

sc_css=css.replace(".gene-card{padding:18px 0 22px;border-top:1px solid #ececec}","")
single=f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Single-cell bridge · Three-section atlas</title><style>{sc_css}</style></head><body>
<header><h1>Single-cell bridge</h1><p class="subtitle">Frozen Fine26 identities · Stage943 expression authority · three-section integration</p></header>
<nav><div class="inner"><a href="index.html">← Spatial atlas</a><a href="#taxonomy">Taxonomy</a><a href="#spatial">Spatial bridge</a><a href="#morphology">Morphology</a><a href="#projection">Projection</a></div></nav><main>
<section id="taxonomy"><h2>Fine26 molecular taxonomy</h2><p class="lead">Fine26 labels remain frozen; the current workflow does not re-cluster merely to create a new version.</p>
<div class="cards"><article class="card"><h3>Current Fine26 marker dotplot</h3><a href="assets/current_fine26_dotplot.jpg" target="_blank"><img src="assets/current_fine26_dotplot.jpg"></a></article>
<article class="card"><h3>Cross-section reproducibility</h3><a href="assets/current_fine26_section_repro.jpg" target="_blank"><img src="assets/current_fine26_section_repro.jpg"></a></article></div></section>
<section id="spatial"><h2>Single-cell ↔ spatial atlas</h2><div class="cards">
<article class="card"><h3>Fine26 in tissue space</h3><a href="assets/stage963_orientation_smooth/FINE26_XY_CURRENT.jpg" target="_blank"><img src="assets/stage963_orientation_smooth/FINE26_XY_CURRENT.jpg"></a></article>
<article class="card"><h3>Marker-guided region10 molecular profiles</h3><a href="assets/stage956_final_atlas_v2/REGION10_GENE_PROFILE.png" target="_blank"><img src="assets/stage956_final_atlas_v2/REGION10_GENE_PROFILE.png"></a><p>Links section-local spatial domains to the 27-gene molecular state.</p></article></div></section>
<section id="morphology"><h2>Morphology</h2><div class="cards">
<article class="card"><h3>Population distributions</h3><a href="assets/current_morph_size_distributions.jpg" target="_blank"><img src="assets/current_morph_size_distributions.jpg"></a></article>
<article class="card"><h3>Section consistency</h3><a href="assets/current_morph_section_consistency.jpg" target="_blank"><img src="assets/current_morph_section_consistency.jpg"></a></article></div></section>
<section id="projection"><h2>Projection-associated organization</h2><div class="cards">
<article class="card"><h3>Projection cohort composition</h3><a href="assets/current_projection_composition.jpg" target="_blank"><img src="assets/current_projection_composition.jpg"></a></article>
<article class="card"><h3>Reproducible enrichments</h3><a href="assets/current_projection_reproducible.jpg" target="_blank"><img src="assets/current_projection_reproducible.jpg"></a></article></div>
<p class="lead">Projection correspondence remains downstream of direct molecular and spatial evidence.</p></section>
<details><summary>Integration plan</summary><p><a href="data/stage949_singlecell_integration_plan.md">Stage949 single-cell integration plan</a></p></details>
</main><footer><a href="index.html">Back to Three-section spatial gene atlas</a></footer></body></html>'''
(PUB/"single-cell.html").write_text(single,encoding="utf-8")

page_auth=dict(stage=957,status="CURRENT_PAGE",main_atlas="Stage956",expression="Stage943 Route A",
               page_sections=["LC/periLC locator","marker-guided section-local region10","81 gene×section maps","downstream spatial analyses","single-cell bridge"],
               intermediate_engineering_hidden=True,single_cell_page="single-cell.html")
(DATA/"stage957_page_authority.json").write_text(json.dumps(page_auth,indent=2),encoding="utf-8")
print(json.dumps(page_auth,indent=2))