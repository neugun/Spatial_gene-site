"""Spatial and anatomical boundary checks on candidate Snap25>5 (NE/ChAT excluded) VGAT/Glu classifications.
All original section XY; raw RouteA genes; source labels preserved, no new tSNE retraining.
"""
from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.ndimage import binary_closing,binary_fill_holes,distance_transform_edt,label
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
R=Path(r'Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926')
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review');O=P/'assets'/'stage999_transmitter_spatial_strictNE';O.mkdir(exist_ok=True)
z=np.load(R/'stage704_current_umap27'/'CURRENT_UMAP27_STAGE704.npz')
genes=list(z['genes'].astype(str))
C=np.concatenate([np.load(Path(r'G:\Map6_recover_all\CURRENT_ROUTEA')/f'S{s}_cell_gene27_current.npy') for s in (500,530,560)],axis=0)
df=pd.read_csv(R/'stage631_current_3d_atlas'/'CURRENT_3D_CELL_ATLAS.csv.gz',usecols=['x_um','y_um','section','fine26_stable_id'])
B=df.fine26_stable_id.to_numpy()
get=lambda name:C[:,genes.index(name)]
snap=get('Snap25')>5
ne=(B=='NE01')|(get('Slc6a2')>200)
chat=(B=='ChAT01')|(get('Slc5a7')>50)
allow=snap&~ne&~chat
v=get('Slc32a1');g=get('Slc17a6')
allrows=[]
for sec in (500,530,560):
 mask=df.section.to_numpy()==sec
 x=df.x_um.to_numpy()[mask];y=df.y_um.to_numpy()[mask]
 binsize=40
 bx=np.clip(np.floor(x/binsize).astype(int),0,32);by=np.clip(np.floor(y/binsize).astype(int),0,32)
 shape=(33,33);occ=np.bincount(by*33+bx,minlength=33*33).reshape(shape)
 cover=binary_closing(occ>=3,structure=np.ones((3,3),bool),iterations=1)
 lab,n=label(cover)
 if n:
  c=np.bincount(lab.ravel());c[0]=0
  cover=lab==np.argmax(c)
 dist=distance_transform_edt(binary_fill_holes(cover))*binsize
 edge=dist[by,bx]<=80
 for tv,te in [(10,5),(15,5),(20,5),(25,8),(30,10),(40,15)]:
  valid=allow[mask]
  both=valid&(v[mask]>tv)&(g[mask]>te)
  med=valid&~edge;ed=valid&edge
  allrows.append(dict(section=sec,t_vgat=tv,t_vglut2=te,eligible=int(valid.sum()),
   n_both=int(both.sum()),both_edge=int((both&edge).sum()),both_interior=int((both&~edge).sum()),
   pct_double_edge=float(100*np.mean(both[ed])) if ed.any() else np.nan,
   pct_double_interior=float(100*np.mean(both[med])) if med.any() else np.nan,
   ratio_edge_to_interior=float(np.mean(both[ed])/np.mean(both[med])) if ed.any() and med.any() and np.mean(both[med])>0 else np.nan))
S=pd.DataFrame(allrows);S.to_csv(P/'data'/'stage999_double_positive_true_tissue_edge_screen_STRICT_NE.csv',index=False)
colors=['#65aec0','#df945e','#9a75aa','#d7dce1']
for tv,te in [(10,5),(15,5),(20,5),(25,8),(40,15)]:
 fig,axs=plt.subplots(1,3,figsize=(14.3,5.2))
 for ax,sec in zip(axs,[500,530,560]):
  m=df.section.to_numpy()==sec
  vloc=v[m];gloc=g[m];a=allow[m]
  cls=np.full(int(m.sum()),4)
  cls[a&~(vloc>tv)&~(gloc>te)]=3
  cls[a&(vloc>tv)&~(gloc>te)]=0
  cls[a&~(vloc>tv)&(gloc>te)]=1
  cls[a&(vloc>tv)&(gloc>te)]=2
  x=df.x_um.to_numpy()[m];y=df.y_um.to_numpy()[m]
  ax.scatter(x[cls==4],y[cls==4],s=.7,color='#e7e7e8',alpha=.12,linewidths=0,rasterized=True)
  for k in [3,0,1,2]:
   flag=cls==k
   ax.scatter(x[flag],y[flag],s=1.2 if k!=2 else 3.8,color=colors[k],alpha=.33 if k==3 else .79,
       linewidths=0,rasterized=True)
  ax.set(xlabel='Original flipped XY X (µm)',ylabel='Original flipped XY Y (µm)',
          title=f'S{sec}: {int((cls==2).sum()):,} double+ / {a.sum():,} eligible')
  ax.set_aspect('equal');ax.invert_yaxis()
  for sp in ['top','right']:ax.spines[sp].set_visible(False)
 fig.suptitle(f'Full sections: corrected counts VGAT>{tv}, VGLUT2>{te}; Snap25>5, excluding NE/ChAT\n'
              'Blue VGAT-only; orange Glu-only; purple double-positive; gray both-negative; light background excluded',fontsize=11)
 fig.tight_layout()
 fig.savefig(O/f'STAGE999_FULL_SECTION_GABA_GLU_COLOCALIZATION_VGAT{tv}_VGLUT{te}.png',dpi=190);plt.close(fig)
fig,ax=plt.subplots(figsize=(8.7,5.1))
for (tv,te),color in zip([(10,5),(20,5),(25,8),(40,15)],['#68aebe','#de9368','#9a75aa','#7188a7']):
 q=S[(S.t_vgat==tv)&(S.t_vglut2==te)]
 ax.plot([500,530,560],q.ratio_edge_to_interior,marker='o',label=f'VGAT>{tv}, VGLUT2>{te}',lw=1.7,color=color)
ax.axhline(1,lw=.8,ls='--',color='#858585')
ax.set(xlabel='Section',ylabel='Double+ rate near tissue edge / deep interior',
 title='Does GABA/Glu colocalization concentrate near tissue boundary?')
ax.set_xticks([500,530,560])
ax.legend(frameon=False,fontsize=8)
for sp in ['top','right']:ax.spines[sp].set_visible(False)
fig.tight_layout();fig.savefig(O/'STAGE999_DOUBLE_POSITIVE_TISSUE_EDGE_ENRICHMENT.png',dpi=190);plt.close(fig)
print('STAGE999_DOUBLE_BOUNDARY',S.to_string(index=False,float_format=lambda x:f'{x:.3f}'),flush=True)
