"""Stage990: entire measured section (not LC/periLC crop) 26-type and broad maps."""
from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(r"Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");D=P/"data";O=P/"assets"/"stage990_fullsection";O.mkdir(exist_ok=True)
A=pd.read_csv(R/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",usecols=["section","x_um","y_um","fine26_stable_id","roi_id"])
K=pd.read_csv(D/"stage977_fine26_color_key.csv")
colors=dict(zip(K.fine26,K.color));bmap=dict(zip(K.fine26,K.broad))
BC={"E":"#ec862f","I":"#397cbd","NE":"#2f9b75","ChAT":"#b75fb3","Other":"#999999"}
A["broad"]=A.fine26_stable_id.map(bmap).fillna("Other")
assert len(A)==71950 and set(A.fine26_stable_id.unique())==set(K.fine26)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})
for mode in ("Fine26","Broad"):
 fig,axs=plt.subplots(1,3,figsize=(16.8,6),facecolor="white")
 counts={}
 for ax,s in zip(axs,[500,530,560]):
  q=A[A.section==s];counts[str(s)]=int(len(q));assert len(q)>18000
  gr="fine26_stable_id" if mode=="Fine26" else "broad"
  palette=colors if mode=="Fine26" else BC
  for typ,rr in q.groupby(gr,sort=True):
   ax.scatter(rr.x_um,rr.y_um,s=1.15 if mode=="Fine26" else 1.25,
    color=palette.get(typ,"#888888"),alpha=.82,linewidths=0,rasterized=True)
  ax.set_aspect("equal",adjustable="box")
  ax.set(xlim=(-25,1050),ylim=(-25,1050),xlabel="Stage631 x (um; X flipped)",ylabel="Stage631 y (um; Y flipped)")
  ax.set_title(f"S{s} · all {len(q):,} cells",fontsize=12)
  ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
 fig.suptitle(f"{mode} spatial atlas · COMPLETE measured XY section, not periLC-only | 71,950 cells",fontsize=15)
 palette=colors if mode=="Fine26" else BC
 items=K.fine26.tolist() if mode=="Fine26" else ["E","I","NE","ChAT","Other"]
 handles=[Line2D([0],[0],marker="o",linestyle="",color=palette[q],label=q,markersize=6) for q in items]
 fig.legend(handles=handles,loc="lower center",ncol=13 if mode=="Fine26" else 5,frameon=False,fontsize=8,bbox_to_anchor=(.5,.015))
 fig.subplots_adjust(bottom=.16,top=.86,left=.07,right=.985,wspace=.24)
 fig.savefig(O/f"STAGE990_{mode.upper()}_FULL_SECTION_XY.png",dpi=205)
 fig.savefig(O/f"STAGE990_{mode.upper()}_FULL_SECTION_XY.pdf")
 plt.close(fig)
 print("STAGE990_SAVED",mode,counts,flush=True)
(D/"stage990_fullsection_authority.json").write_text(json.dumps({"stage":990,"all_cells":71950,"section_counts":counts,"source":"Stage631 full x/y cell atlas, without medial/periLC/LC ROI gate","orientation":"Stage631 x and y already flipped; do not flip again","visuals":["Fine26 all-section XY","Broad all-section XY"]},indent=2))
