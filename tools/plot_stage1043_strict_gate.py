from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review');D=P/'data';O=P/'assets'/'stage1043_strict_transmitter_gate'
R=Path(r'Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926')
A=pd.read_csv(R/'stage631_current_3d_atlas'/'CURRENT_3D_CELL_ATLAS.csv.gz',usecols=['section','roi_id','x_um','y_um','fine26_stable_id'])
mask=np.load(Path(r'G:\Map6_recover_all')/'stage1043_private_filter_masks.npz')
u=mask['selected'];cl=mask['class3'];ref=mask['ref']
colors=['#3f7fad','#dd8b52','#9559a6']
fig,axs=plt.subplots(1,3,figsize=(15.8,6.3))
for ax,sec in zip(axs,[500,530,560]):
 m=A.section.to_numpy(int)==sec
 ax.scatter(A.x_um[m],A.y_um[m],s=.4,alpha=.23,color='#d5d9dc',linewidth=0,rasterized=True)
 for k in [0,1,2]:
  q=m&(cl==k)
  ax.scatter(A.x_um[q],A.y_um[q],s=1.1 if k!=2 else 2.0,c=colors[k],alpha=.80,linewidth=0,rasterized=True)
 ax.set(xlim=(0,1040),ylim=(0,1040),xlabel='Stage631 X (already flipped), µm',ylabel='Stage631 Y (already flipped), µm',
        title=f'S{sec}: retained {np.sum(m&u):,}, dual {np.sum(m&(cl==2)):,}')
 ax.set_aspect('equal')
 for side in ('top','right'):ax.spines[side].set_visible(False)
fig.suptitle('Snap25>5, exclude frozen NE/ChAT, VGAT>10 OR VGLUT2>5: true full-section XY',fontsize=12)
fig.subplots_adjust(left=.06,right=.985,bottom=.10,top=.89,wspace=.25)
fig.savefig(O/'STAGE1043_USER_GATE_ALL_SECTIONS_SPATIAL.png',dpi=180,bbox_inches='tight');plt.close(fig)
pos=Path(r'G:\Map6_recover_all\stage995_tsne\STAGE995_CURRENT_LOGZ_PCA20_TSNE_30_100.npz')
z=np.load(pos);print('STAGE995_KEYS',z.files,flush=True)
xy=next((np.asarray(z[k]) for k in z.files if np.asarray(z[k]).shape==(len(A),2)),None)
assert xy is not None,'no known tSNE coordinates'
fig,axs=plt.subplots(1,2,figsize=(12,6.4))
for ax,keep,title in zip(axs,[u,u&ref],['USER Snap25>5 + OR gate','Additional independent ref-neuron gate']):
 ax.scatter(xy[:,0],xy[:,1],c='#e1e4e7',s=.2,alpha=.13,linewidth=0,rasterized=True)
 for k in (0,1,2):
  sel=keep&(cl==k)
  ax.scatter(xy[sel,0],xy[sel,1],s=1.0,c=colors[k],alpha=.74,linewidth=0,rasterized=True)
 ax.set_title(f'{title}\nretained {keep.sum():,}')
 ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect('equal')
 for side in ax.spines.values():side.set_visible(False)
fig.suptitle('Filtered cells on original unsupervised 71,950-cell Stage995 t-SNE; not refitted')
fig.subplots_adjust(top=.87,bottom=.07,wspace=.15)
fig.savefig(O/'STAGE1043_USER_GATE_ORIGINAL_STAGE995_TSNE.png',dpi=190,bbox_inches='tight');plt.close(fig)
# Exact overlap split by frozen E/I and subtype, expose circularity
M=pd.read_csv(D/'stage977_fine26_color_key.csv');b=A.fine26_stable_id.map(dict(zip(M.fine26,M.broad))).to_numpy()
rows=[]
for sec in [0,500,530,560]:
 q=(A.section.to_numpy(int)==sec) if sec else np.ones(len(A),bool)
 for broad in ['E','I','Other']:
  a=q&(b==broad)&u
  if not a.any():continue
  rows.append(dict(section=sec,broad=broad,retained=int(a.sum()),onlyVGAT=int(np.sum(a&(cl==0))),
   onlyVGLUT2=int(np.sum(a&(cl==1))),double=int(np.sum(a&(cl==2))),dual_pct=100*np.mean(cl[a]==2)))
pd.DataFrame(rows).to_csv(D/'stage1043_frozen_broad_versus_actual_marker_call.csv',index=False)
print('STAGE1043_PLOTS_DONE',len(rows),flush=True)
