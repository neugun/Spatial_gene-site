"""Examine all NE01 cluster cells with per-cell VGAT/Slc32a1 and Glu/Slc17a6 counts.
Compare Stage631 original catalog per-cell integer assignment (OLD), to CURRENT
RouteA corrected count equivalents. No forced labels, no threshold tuning.
"""
from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
NE_CMAP=LinearSegmentedColormap.from_list('pale_cyan_red',['#d9f6f4','#a9dce7','#f0dcce','#e58c71','#a32545','#68122e'])
from scipy.stats import spearmanr
R=Path(r'Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926')
Q=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review')
D=Q/'data';F=Q/'assets'/'stage1049_NE01_GABA_GLU_cell_spots';F.mkdir(exist_ok=True)
cols=['cell_id','roi_id','section','x_um','y_um','fine26_stable_id','Slc32a1','Slc17a6','Slc6a2']
A=pd.read_csv(R/'stage631_current_3d_atlas'/'CURRENT_3D_CELL_ATLAS.csv.gz',usecols=cols)
C=np.load(r'G:\Map6_recover_all\CURRENT_ROUTEA\current_routeA_all27.npy')
gn=list(np.load(R/'stage704_current_umap27'/'CURRENT_UMAP27_STAGE704.npz')['genes'].astype(str))
E=np.load(r'G:\Map6_recover_all\stage995_tsne\STAGE995_CURRENT_LOGZ_PCA20_TSNE_30_100.npz')['embedding']
assert len(A)==len(C)==len(E)==71950
label=A.fine26_stable_id.eq('NE01').to_numpy()
M=A.loc[label].copy()
idx=np.flatnonzero(label)
for g in ['Slc32a1','Slc17a6','Slc6a2','Snap25']:
 M[f'{g}_current']=C[label,gn.index(g)]
M['tsne_x']=E[label,0];M['tsne_y']=E[label,1]
M['stage631_VGAT_integer']=M.Slc32a1
M['stage631_VGLUT2_integer']=M.Slc17a6
M['vgat_baseline_diff']=M.Slc32a1_current-M.Slc32a1
M['vglut_baseline_diff']=M.Slc17a6_current-M.Slc17a6
M=M.rename(columns={'Slc32a1':'VGAT_Stage631', 'Slc17a6':'VGLUT2_Stage631', 'Slc6a2':'Slc6a2_Stage631'})
M.to_csv(D/'stage1049_NE01_all1833_cells_VGAT_VGLUT2_original_baseline_and_corrected.csv',index=False)
gene_pairs=[('VGAT','VGAT_Stage631','Slc32a1_current'),('VGLUT2','VGLUT2_Stage631','Slc17a6_current')]
qs=[]
for sec in ['ALL',500,530,560]:
 X=M if sec=='ALL' else M[M.section==sec]
 for g,old,cur in gene_pairs:
  for source,col in [('Stage631_previous_integer',old),('Stage943_RouteA_corrected',cur)]:
   v=X[col].to_numpy(float)
   q=dict(section=sec,gene=g,source=source,n=len(v),mean=float(v.mean()),median=float(np.median(v)),
       q25=float(np.quantile(v,.25)),q75=float(np.quantile(v,.75)),q90=float(np.quantile(v,.9)),
       q95=float(np.quantile(v,.95)),q99=float(np.quantile(v,.99)),max=float(v.max()),
       pct_eq0=float(np.mean(v==0)*100),pct_gt0=float(np.mean(v>0)*100))
   for t in [1,2,3,5,10,15,20,30,40,50]:
    q[f'pct_gt{t}']=float(np.mean(v>t)*100)
   qs.append(q)
S=pd.DataFrame(qs);S.to_csv(D/'stage1049_NE01_VGAT_VGLUT2_by_section_spot_count_statistics.csv',index=False)
joint=[]
for sec in ['ALL',500,530,560]:
 X=M if sec=='ALL' else M[M.section==sec]
 for source,gg,ee in [('Stage631_integer','VGAT_Stage631','VGLUT2_Stage631'),
                       ('Stage943_corrected','Slc32a1_current','Slc17a6_current')]:
  g=X[gg].to_numpy(float);e=X[ee].to_numpy(float)
  for a,b in [(0,0),(5,2),(10,5),(15,5),(20,5),(20,10),(30,10),(40,15)]:
   bg=g>a;be=e>b;n=len(X)
   joint.append(dict(section=sec,source=source,VGAT_cutoff=a,VGLUT2_cutoff=b,n=n,
                     vg_only=int(np.sum(bg&~be)),glut_only=int(np.sum(~bg&be)),
                     both=int(np.sum(bg&be)),neither=int(np.sum(~bg&~be)),
                     pct_both=float(100*np.mean(bg&be)),pct_vgat=float(100*np.mean(bg)),
                     pct_vglut2=float(100*np.mean(be)),both_given_either=float(np.sum(bg&be)/max(np.sum(bg|be),1)*100),
                     spearman_counts=float(spearmanr(g,e).statistic)))
J=pd.DataFrame(joint);J.to_csv(D/'stage1049_NE01_joint_VGAT_VGLUT2_threshold_counts.csv',index=False)
# Histograms all data: 0 etc, fractions from same source; pale warm colors
fig,axs=plt.subplots(2,3,figsize=(14.1,8.2))
for j,(g,old,cur) in enumerate(gene_pairs):
 vals=M[[old,cur]].to_numpy(float)
 axs[j,0].hist(vals[:,0],bins=np.r_[-.5,np.arange(.5,80.5,2)],density=True,alpha=.55,color='#61afbe',label='Stage631 integer')
 axs[j,0].hist(vals[:,1],bins=np.r_[-.5,np.arange(.5,80.5,2)],density=True,alpha=.46,color='#c87663',label='RouteA corrected')
 axs[j,0].set(xlabel=f'{g} count (0-80 range)',ylabel='Density',title=f'{g} NE01 count density, n={len(M)}')
 axs[j,0].legend(frameon=False,fontsize=7.8)
 for sec,c in zip([500,530,560],['#63aabb','#d29364','#a36caa']):
  v=M.loc[M.section==sec,cur].to_numpy(float)
  x=np.sort(v)
  y=100*np.arange(1,len(x)+1)/len(x)
  axs[j,1].plot(x,y,color=c,lw=1.6,label=f'S{sec} n={len(x)}')
 axs[j,1].set(xlabel=f'RouteA {g} corrected count',ylabel='Cumulative % of NE01 cells',title=f'{g} continuous distribution')
 axs[j,1].set_xscale('symlog',linthresh=8)
 axs[j,1].legend(frameon=False,fontsize=8)
 for sec,c in zip([500,530,560],['#63aabb','#d29364','#a36caa']):
  v=M.loc[M.section==sec,cur].to_numpy(float)
  cuts=[0,2,5,10,15,20,30,40]
  axs[j,2].plot(cuts,[100*np.mean(v>t) for t in cuts],marker='o',ms=3.5,lw=1.5,color=c,label=f'S{sec}')
 axs[j,2].set(xlabel=f'{g} cutoff (corrected)',ylabel='% NE01 exceeding cutoff',title=f'{g} threshold sensitivity',ylim=(0,105))
 axs[j,2].legend(frameon=False,fontsize=8)
 for ax in axs[j]:
  ax.grid(False)
  for sp in ['top','right']:ax.spines[sp].set_visible(False)
fig.suptitle('All frozen NE01 cells: previous per-cell integer assignment vs corrected RouteA count equivalents\nNeither is independently verified true RNA per cell; NE01 labels partly defined by same 27-gene panel',fontsize=11)
fig.tight_layout();fig.savefig(F/'STAGE1049_NE01_VGAT_VGLUT2_CELL_COUNT_DISTRIBUTIONS.png',dpi=185);plt.close(fig)
# actual NE01 cluster, scatter one panel broad t-SNE gray and zoom NE01 scatter colored
fig,axs=plt.subplots(2,3,figsize=(15.2,9.2))
for j,(g,old,cur) in enumerate(gene_pairs):
 color='#d9f2f4'
 ax=axs[j,0]
 ax.scatter(E[:,0],E[:,1],s=.34,alpha=.11,color='#d5dadd',linewidths=0,rasterized=True)
 ax.scatter(M.tsne_x,M.tsne_y,s=2.5,alpha=.85,color='#b75478',linewidths=0,rasterized=True)
 ax.set_title(f'Frozen NE01 on all 71,950 t-SNE\nn={len(M)}, gene={g}')
 for k,(col,title) in enumerate([(old,'Earlier Stage631 integer count'),(cur,'Corrected RouteA count-equivalents')],start=1):
  ax=axs[j,k];v=M[col].to_numpy(float)
  ax.scatter(M.tsne_x,M.tsne_y,s=1.6,c='#dfe7e9',alpha=.37,linewidths=0,rasterized=True)
  positive=v>0
  q=np.quantile(v[positive],.98) if positive.any() else 1
  p=ax.scatter(M.tsne_x[positive],M.tsne_y[positive],s=3.2,c=np.clip(v[positive],0,q),
          cmap=NE_CMAP,vmin=0,vmax=q,alpha=.94,linewidths=0,rasterized=True)
  fig.colorbar(p,ax=ax,shrink=.70,aspect=20,label=f'{g} count')
  ax.set_title(f'{g}: {title}\nmedian={np.median(v):.1f}; >0 {100*positive.mean():.1f}%')
 for ax in axs[j]:
  ax.set_xticks([]);ax.set_yticks([])
  for sp in ax.spines.values():sp.set_visible(False)
fig.suptitle('Actual NE01 cluster on frozen unsupervised Stage995 t-SNE: GABA and Glu spots per cell\nValue and spatial pattern displayed WITHOUT re-clustering or redefining NE identities',fontsize=11)
fig.tight_layout();fig.savefig(F/'STAGE1049_NE01_TSNE_VGAT_VGLUT2_OLD_VS_CURRENT_PER_CELL.png',dpi=190);plt.close(fig)
# joint RNA counts in NE cells
fig,axs=plt.subplots(1,3,figsize=(14.3,4.8))
for ax,sec in zip(axs,[500,530,560]):
 X=M[M.section==sec]
 g=X.Slc32a1_current.to_numpy(float);e=X.Slc17a6_current.to_numpy(float)
 ax.hexbin(g,e,gridsize=26,mincnt=1,bins='log',cmap='YlGnBu',extent=(0,100,0,60))
 ax.axvline(10,c='#ba5a79',ls='--',lw=1)
 ax.axhline(5,c='#ba5a79',ls='--',lw=1)
 ax.set(xlim=(0,100),ylim=(0,60),xlabel='VGAT RouteA count (0-100)',ylabel='VGLUT2 RouteA count (0-60)',
     title=f'S{sec} NE01 cells n={len(X)}\nVGAT>10, VGLUT2>5 joint {100*np.mean((g>10)&(e>5)):.1f}%')
 for sp in ['top','right']:ax.spines[sp].set_visible(False)
fig.suptitle('Distribution of VGAT and VGLUT2 corrected count-equivalents within existing NE01, by section\nNo artificial threshold optimization and counts above display ranges not shown',fontsize=11)
fig.tight_layout();fig.savefig(F/'STAGE1049_NE01_VGAT_VGLUT2_PER_CELL_JOINT_COUNTS.png',dpi=190);plt.close(fig)
# summary original vs corrected diagnostic
out={'status':'COMPUTED_NE01_SOURCE_COUNT_DISTRIBUTIONS','n_NE01':len(M),'NE01_per_section':M.section.value_counts().sort_index().to_dict(),
'original_baseline':'Stage631 per-cell integer count array, predates RouteA corrective stage 943; not verified pixel-level original RSFISH truth',
'current':'Stage943 RouteA corrected count-equivalents (fractional), cannot equate to literal physical number of puncta',
'cluster_identity':'Frozen Fine26 NE01, not independent confirmation of LC noradrenergic identity',
'tsne':'Stage995 unsupervised 71,950 full cells, fixed embedding, only selected NE01 for display',
'no_reclassification':True,'no_force_double_pos':True}
(D/'stage1049_NE01_GABA_GLU_spot_count_audit.json').write_text(json.dumps(out,indent=2))
print('NE01_N',len(M),'SECTIONS',out['NE01_per_section'],flush=True)
print('NE01_SUMMARY',S.to_string(index=False,float_format=lambda x:f'{x:.3f}'),flush=True)
print('NE01_JOINT',J[(J.section=='ALL')].to_string(index=False,float_format=lambda x:f'{x:.3f}'),flush=True)
