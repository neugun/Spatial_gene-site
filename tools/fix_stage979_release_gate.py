from pathlib import Path
p=Path(r"G:\Spatial_gene_site_publish\tools\public_site_audit.py")
s=p.read_text(encoding="utf-8")
old='if bridge.count("stage974_preflipped_reference/") < 15:'
new='if bridge.count("stage974_preflipped_reference/") + bridge.count("stage979_consistent_fine26/") < 15:'
assert s.count(old)==1
p.write_text(s.replace(old,new),encoding="utf-8")
print("STAGE979_PUBLIC_GATE_FIXED",flush=True)
