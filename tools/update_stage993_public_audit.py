from pathlib import Path
import json
P=Path(r"G:\Spatial_gene_site_publish")
a=P/"perilc-map6-review"/"data"/"stage992_true_mask_surface_authority.json"
q=json.loads(a.read_text(encoding="utf8"))
q["source"]="Full-resolution S500 downseg.tif genuine segmentation labels; source disk path withheld"
a.write_text(json.dumps(q,indent=2),encoding="utf8")
s=P/"tools"/"build_stage992_true_s500_surface.py"
x=s.read_text(encoding="utf8")
x=x.replace('"source":str(TIFF),','"source":"Full-resolution S500 downseg.tif genuine segmentation labels; raw source path withheld",')
s.write_text(x,encoding="utf8")
f=P/"tools"/"public_site_audit.py";h=f.read_text(encoding="utf8")
old='''    if bridge.count("stage974_preflipped_reference/") + bridge.count("stage979_consistent_fine26/") < 15:
        errors.append("map10 current 3D maskbody links are stale")'''
new='''    # Stage993: current main page must now contain complete-section Fine26/Broad
    # and orthogonal physical segmentation 3D. Historical Stage974/979 counts
    # were a legacy view and may no longer be shown.
    if bridge.count("stage990_fullsection/") < 4:
        errors.append("map10 full-section (not periLC crop) XY atlas missing")
    if bridge.count("stage991_orthogonal_maskbody/") < 4:
        errors.append("map10 whole-section XY XZ YZ segmentation missing")
    if bridge.count("stage992_true_mask_surface/") < 2:
        errors.append("map10 S500 genuine local segmentation mesh missing")
    import re
    if re.search(r'<img[^>]+src=["\\\'][^"\\\']*[Uu][Mm][Aa][Pp]',bridge):
        errors.append("map10 current single-cell page still shows deprecated UMAP image")'''
assert old in h
h=h.replace(old,new)
f.write_text(h,encoding="utf8")
print("STAGE993_AUDIT_GATE_PATCHED, no source path leak",flush=True)
