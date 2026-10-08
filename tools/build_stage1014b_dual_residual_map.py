"""Stage1014b: real neuronal dual-pos spatial fractions and section×Fine26
expected residuals; do not mistake higher cell density for enrichment.
"""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
R=Path(r"G:\Map6_recover_all\stage1014_spatial_doublepositive_private")
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
D=P/"data";O=P/"assets"/"stage1014_doublepositive_spatial"
binwidth=65;edges=np.arange(0,1040+binwidth,binwidth)
fig,axs=plt.subplots(2,3,figsize=(16.5,10.5),facecolor="white")
dat=[];matrices={}
for col,ss in enumerate((500,530,560)):
 f=R/f"raw_VGAT10_VGLUT2_3_S{ss}_cells.csv.gz";df=pd.read_csv(f)
 oldflags=pd.read_csv(Path(r"G:\Map6_recover_all\stage1001_neuron_gate_cell_flags_PRIVATE.csv.gz"),usecols=["section","roi_id","reference_neurons"])
 old=oldflags[oldflags.section==ss]
 df=df.merge(old,on=["section","roi_id"],how="left",validate="one_to_one")
 df=df[df.reference_neurons].copy()
 df["dual"]=(df.class4==2).astype(int)
 # expected per Fine26 in each section, no E/I typology circularity in outcome
 q=df.groupby("Fine26").dual.agg(["mean","count"])
 df["fine26_expect"]=df.Fine26.map(q["mean"]).astype(float)
 ix=np.searchsorted(edges,df.x_um,side="right")-1
 iy=np.searchsorted(edges,df.y_um,side="right")-1
 ix=np.clip(ix,0,len(edges)-2);iy=np.clip(iy,0,len(edges)-2)
 nx=len(edges)-1
 sample=np.zeros((nx,nx),int);obs=np.zeros((nx,nx),float);exp=np.zeros((nx,nx),float)
 np.add.at(sample,(iy,ix),1)
 np.add.at(obs,(iy,ix),df.dual)
 np.add.at(exp,(iy,ix),df.fine26_expect)
 frac=np.divide(obs,sample,out=np.full(obs.shape,np.nan),where=sample>=20)
 resid=np.divide(obs-exp,sample,out=np.full(obs.shape,np.nan),where=sample>=20)
 frac[~(sample>=20)]=np.nan;resid[~(sample>=20)]=np.nan
 matrices[ss]=(frac,resid,sample)
 dat.append(dict(section=ss,n=int(len(df)),bin_um=binwidth,bin_min_n=20,sections_one_mouse=True,
  frac_min=float(np.nanmin(frac)),frac_max=float(np.nanmax(frac)),
  residual_min=float(np.nanmin(resid)),residual_max=float(np.nanmax(resid))))
# maps calibrated across all sections
allres=np.concatenate([matrices[s][1].ravel() for s in [500,530,560]])
lim=max(.08,float(np.nanquantile(np.abs(allres),.985)))
for col,ss in enumerate((500,530,560)):
 frac,res,sample=matrices[ss]
 ax=axs[0,col]
 im=ax.imshow(frac,origin="lower",extent=[0,edges[-1],0,edges[-1]],vmin=0,vmax=.65,cmap="PuRd",interpolation="nearest")
 ax.set_title(f"S{ss} · double-positive fraction\nVGAT>10 / VGLUT2>3, neuron-only")
 ax.set_aspect("equal")
 ax.set(xlim=(0,1030),ylim=(0,1030),xlabel="Stage631 X (µm)",ylabel="Stage631 Y (µm)")
 ax=axs[1,col]
 z=ax.imshow(res,origin="lower",extent=[0,edges[-1],0,edges[-1]],norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),cmap="RdBu_r",interpolation="nearest")
 ax.set_title(f"S{ss} · observed − Fine26-expected double+ fraction\nSame spatial bins, >=20 neurons/bin")
 ax.set_aspect("equal")
 ax.set(xlim=(0,1030),ylim=(0,1030),xlabel="Stage631 X (µm)",ylabel="Stage631 Y (µm)")
fig.subplots_adjust(left=.055,right=.92,top=.92,bottom=.07,wspace=.14,hspace=.30)
cax1=fig.add_axes([.94,.58,.012,.30]);cb1=fig.colorbar(im,cax=cax1);cb1.set_label("Double+ / neurons")
cax2=fig.add_axes([.94,.15,.012,.30]);cb2=fig.colorbar(z,cax=cax2);cb2.set_label("Observed − expected")
fig.suptitle("Full-section dual-positive geography, with cell-type-conditional spatial residual (one specimen)",fontsize=15)
fig.savefig(O/"STAGE1014B_DUAL_FRACTION_AND_FINE26_ADJUSTED_SPATIAL.png",dpi=195,bbox_inches="tight")
plt.close(fig)
pd.DataFrame(dat).to_csv(D/"stage1014b_fine26_adjusted_spatial_bins_metadata.csv",index=False)
print("STAGE1014B_SPATIAL_RESIDUAL_COMPLETE",dat,flush=True)
