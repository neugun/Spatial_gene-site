#!/usr/bin/env python3
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit, unquote
import re, sys

ROOT=Path(__file__).resolve().parents[1]
GAP=ROOT/"GAP"
TEXT_EXT={".html",".md",".json",".csv",".txt",".css",".yml",".yaml"}
DRIVE_RE=re.compile(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/][^\s\"'<>;,|]+")
PRIVATE_TOKENS=("sternsonlab","GAP_release_20261001","7_Grantwriting","analysis_workspace")
errors=[]

for p in ROOT.rglob("*"):
    rel=p.relative_to(ROOT)
    if not p.is_file() or ".git" in p.parts or (rel.parts and rel.parts[0] in {"tools",".github"}) or p.suffix.lower() not in TEXT_EXT:
        continue
    try: s=p.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        try: s=p.read_text(encoding="utf-8-sig")
        except: continue
    for m in DRIVE_RE.finditer(s):
        raw=m.group(0)
        if raw.lower().startswith("x:\\path\\to\\") or raw.lower().startswith("x:/path/to/"):
            continue
        errors.append(f"path leak: {p.relative_to(ROOT)} :: {raw[:120]}")
    low=s.lower()
    for tok in PRIVATE_TOKENS:
        if tok.lower() in low:
            errors.append(f"private token: {p.relative_to(ROOT)} :: {tok}")

main=GAP/"index.html"
if not main.exists():
    errors.append("missing GAP/index.html")
else:
    n=main.stat().st_size
    if n>180_000:
        errors.append(f"GAP/index.html too large: {n} bytes > 180000")
css=GAP/"docs/site_assets/20261004/site.css"
if css.exists():
    s=css.read_text(encoding="utf-8")
    if re.search(r"\.figure\s+img[^}]*min-width\s*:\s*520px",s,re.I):
        errors.append("mobile regression: figure img min-width:520px")

class P(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]; self.ids=set(); self.base=None
    def handle_starttag(self,tag,attrs):
        d=dict(attrs)
        if tag=="base" and d.get("href"): self.base=d["href"]
        if d.get("id"): self.ids.add(d["id"])
        for k in ("href","src"):
            if d.get(k): self.links.append((k,d[k]))

main_parser=P(); main_parser.feed(main.read_text(encoding="utf-8"))
main_ids=main_parser.ids

def resolve_local(html_path, parser, raw):
    if raw.startswith(("http://","https://","mailto:","data:","javascript:")):
        return
    u=urlsplit(raw)
    if not u.path:
        if u.fragment and html_path==main and u.fragment not in main_ids:
            errors.append(f"missing anchor in GAP/index.html: #{u.fragment}")
        return
    base_dir=html_path.parent
    if parser.base:
        base_dir=(html_path.parent/urlsplit(parser.base).path).resolve()
    target=(base_dir/unquote(u.path)).resolve()
    if str(u.path).endswith("/") or target.is_dir():
        target=target/"index.html"
    if not target.exists():
        errors.append(f"broken local link: {html_path.relative_to(ROOT)} -> {raw}")
    if u.fragment and target==main.resolve() and u.fragment not in main_ids:
        errors.append(f"missing target anchor: {raw}")

for rel in ["GAP/index.html","GAP/docs/index.html"]:
    p=ROOT/rel
    parser=P(); parser.feed(p.read_text(encoding="utf-8"))
    for _,u in parser.links: resolve_local(p,parser,u)

for p in sorted((GAP/"details").glob("*.html")):
    s=p.read_text(encoding="utf-8")
    m=re.search(r'url=\.\./#([A-Za-z0-9_-]+)',s)
    if not m:
        errors.append(f"detail page is not a redirect: {p.relative_to(ROOT)}")
    elif m.group(1) not in main_ids:
        errors.append(f"detail redirect missing main anchor: {p.name} -> {m.group(1)}")

blocked_ext={".py",".ipynb",".mat",".h5",".hdf5",".h5ad",".parquet",".pkl",".pickle",".npy",".npz",".pt",".pth",".ckpt",".ps1",".bat",".sh"}
for p in GAP.rglob("*"):
    if p.is_file() and p.suffix.lower() in blocked_ext:
        errors.append(f"non-results artifact in GAP public tree: {p.relative_to(ROOT)}")

# Map10 release gate: every CURRENT spatial output must use x flipped, y flipped.
# Lightweight stdlib check so GitHub Actions can prevent stale figure regressions.
MAP10=ROOT/"perilc-map6-review"
if MAP10.exists():
    import csv
    for name in ("index.html","single-cell.html"):
        page=MAP10/name
        html=page.read_text(encoding="utf-8")
        if "reversed" in html.lower() or "xy_reversed" in html.lower():
            errors.append(f"map10 obsolete orientation in {name}")
        parser=P();parser.feed(html)
        for _,src in parser.links:
            if src.startswith("assets/") and not (MAP10/src.split("?")[0].split("#")[0]).is_file():
                errors.append(f"map10 missing asset {name}: {src}")
        if "stage951_maskbody/" in html:
            errors.append(f"map10 obsolete maskbody panel in {name}")
    current=(MAP10/"index.html").read_text(encoding="utf-8")
    bridge=(MAP10/"single-cell.html").read_text(encoding="utf-8")
    if current.count("stage970_holefree_xy/") < 164:
        errors.append("map10 some current 81 gene panels are stale")
    for name in ("FINAL_LC_PERILC_LOCATOR_X_FLIPPED_Y_FLIPPED.png",
                 "FINAL_REGION10_HOLEFREE_X_FLIPPED_Y_FLIPPED.png"):
        if name not in current:errors.append(f"map10 missing current orientation: {name}")
    if bridge.count("stage967_maskbody_xy/") < 10:
        errors.append("map10 current 3D maskbody links are stale")
    for name in ("stage956_gene_section_manifest.csv","stage970_gene_section_manifest.csv"):
        with (MAP10/"data"/name).open(encoding="utf-8",newline="") as f:
            rows=list(csv.DictReader(f))
        if len(rows)!=81 or any(r.get("x_flipped")!="True" or r.get("y_flipped")!="True" for r in rows):
            errors.append(f"map10 flip metadata incomplete: {name}")
    with (MAP10/"data"/"stage971_all3_holefree_manifest.csv").open(encoding="utf-8",newline="") as f:
        composites=list(csv.DictReader(f))
    if len(composites)!=27 or any(not (MAP10/"assets"/"stage971_all3_holefree_xy"/r["file"]).is_file() for r in composites):
        errors.append("map10 27 corrected three-section composites missing")
    # Region topology is a testable release contract, not just an image caption.
    import json
    qa=json.loads((MAP10/"data"/"stage970_holefree_region_audit.json").read_text(encoding="utf-8"))
    if qa.get("stage")!=970 or qa.get("gene_maps")!=81:
        errors.append("map10 Stage970 geometry authority missing")
    for sec in ("500","530","560"):
        row=qa.get("region_QA",{}).get(sec,{})
        topo=row.get("after",{})
        if topo.get("white_holes")!=0 or topo.get("white_pixels")!=0:
            errors.append(f"map10 S{sec} tissue holes remain")
        if len(topo.get("region_components",{}))!=10 or any(x!=1 for x in topo.get("region_components",{}).values()):
            errors.append(f"map10 S{sec} disconnected region remains")
        if len(topo.get("region_holes",{}))!=10 or any(x!=0 for x in topo.get("region_holes",{}).values()):
            errors.append(f"map10 S{sec} internal region hole remains")
        if row.get("changed_existing_fraction",1)>.03:
            errors.append(f"map10 S{sec} boundary moved too far")
    if current.count("stage971_all3_holefree_xy/")!=27:
        errors.append("map10 missing current section-composite links")

if errors:
    print("PUBLIC SITE AUDIT: FAIL")
    for e in errors[:300]: print(" -",e)
    if len(errors)>300: print(f" ... {len(errors)-300} more")
    sys.exit(2)
print("PUBLIC SITE AUDIT: PASS")
print(f"main_bytes={main.stat().st_size}")
print(f"main_anchors={len(main_ids)}")
print(f"details_redirects={len(list((GAP/'details').glob('*.html')))}")
