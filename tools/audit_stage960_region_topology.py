"""Stage960: topology QA for the current section-local 10-region marker reference.
Research QA only; does not overwrite Stage943/Stage956 authorities.
"""
from pathlib import Path
import json, csv
import numpy as np
from scipy.ndimage import label, binary_dilation, binary_erosion
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root=Path(r"G:\Spatial_gene_site_publish")
src=(root/"tools"/"build_stage956_final_atlas.py").read_text(encoding="utf-8")
prefix=src.split('pd.DataFrame(sweep_rows).to_csv')[0]
assert 'region_info[ss]=dict' in prefix
scope={"__name__":"stage960_diagnostic"}
exec(compile(prefix,"stage956_unchanged_reconstruction","exec"),scope)
out=root/"perilc-map6-review"/"assets"/"stage960_region_topology"
out.mkdir(exist_ok=True)
dat=root/"perilc-map6-review"/"data"
regions=scope["region_info"]
report={"stage":960,"parent_stage":956,"authority":"diagnostic only","sections":{}}
rows=[]
fig,axs=plt.subplots(2,3,figsize=(15,9))
for col,ss in enumerate((500,530,560)):
    ri=regions[ss]
    G=ri["local"]; M=ri["mask"]
    compmask,nmask=label(M)
    ax=axs[0,col]
    ax.imshow(np.ma.masked_where(~M,G),origin="lower",interpolation="nearest",cmap="tab10",vmin=1,vmax=10)
    ax.set_title(f"S{ss}: current local 10 regions | mask parts={nmask}")
    all_extra=np.zeros(M.shape,bool)
    chunks=[]
    for k in range(1,11):
        cc,n=label((G==k)&M)
        sizes=np.bincount(cc.ravel())[1:]
        if len(sizes)==0:continue
        big=np.argmax(sizes)+1
        for cid,size in enumerate(sizes,1):
            island=(cc==cid)
            touching=binary_dilation(island)&M&(~island)
            touching_classes=G[touching]
            touching_classes=touching_classes[(touching_classes>0)&(touching_classes!=k)]
            neighbour=int(np.bincount(touching_classes,minlength=11).argmax()) if len(touching_classes) else None
            is_minor=cid!=big
            if is_minor:all_extra |= island
            ys,xs=np.where(island)
            rec={"section":ss,"region":k,"component":cid,"pixels":int(size),
              "fraction_of_region":float(size/sizes.sum()),"minor":bool(is_minor),
              "bbox_x0":int(xs.min()),"bbox_x1":int(xs.max()),
              "bbox_y0":int(ys.min()),"bbox_y1":int(ys.max()),
              "neighbor_region":neighbour,
              "mask_component":int(np.bincount(compmask[island],minlength=nmask+1).argmax())}
            rows.append(rec)
            if is_minor:chunks.append(rec)
    axs[1,col].imshow(np.ma.masked_where(~M,G),origin="lower",interpolation="nearest",cmap="tab10",vmin=1,vmax=10,alpha=.4)
    if all_extra.any():
        axs[1,col].contour(all_extra,levels=[.5],colors="red",linewidths=1.2)
        ys,xs=np.where(all_extra)
        axs[1,col].scatter(xs,ys,s=.2,c="red")
    axs[1,col].set_title(f"S{ss}: minor components marked red (n={len(chunks)})")
    for row in axs[:,col]:row.set_xticks([]);row.set_yticks([])
    report["sections"][str(ss)]={"mask_components":int(nmask),
     "component_count":len([r for r in rows if r["section"]==ss]),
     "islands":chunks,
     "small_island_fraction_of_tissue":float(all_extra.sum()/max(M.sum(),1)),
     "mask_coverage_grid_pixels":int(M.sum())}
fig.tight_layout()
fig.savefig(out/"REGION10_TOPOLOGY_QA.png",dpi=180)
plt.close(fig)
with (dat/"stage960_region10_components.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
    w.writeheader();w.writerows(rows)
(dat/"stage960_region10_topology_audit.json").write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding="utf-8")
print(json.dumps(report,ensure_ascii=False,indent=2))
