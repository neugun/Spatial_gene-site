from pathlib import Path
import re, json
d=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review\data")
changed=[]
for p in sorted(d.glob("stage899_*_same_detector_result.json")):
    old=p.read_text(encoding="utf8")
    if "sternsonlab" not in old.lower():continue
    clean,n=re.subn(r'("(?:baseline_path|fresh_path)"\s*:\s*)"(?:\\.|[^"\\])*"',r'\1"[REDACTED_LOCAL_PATH]"',old)
    assert n==2 and "sternsonlab" not in clean.lower(),p.name
    a=json.loads(old);b=json.loads(clean)
    assert sorted(a)==sorted(b)
    for k in a:
        if k not in ("baseline_path","fresh_path"):assert a[k]==b[k],(p.name,k)
    p.write_text(clean,encoding="utf8")
    changed.append(p.name)
assert len(changed)==6,changed
print("REDACTED_STAGE899_FILES",len(changed),"SCIENTIFIC_FIELDS_UNCHANGED",True)
