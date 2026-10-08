"""Stage1044: for every *unmatched VGAT bright raw candidate* check listed 
RS-FISH detections at same XY through +/- 2,4,8,16,32,64 source Z planes.
Avoid counting unrelated x/y nearby without matching actual Z.
Reconstruct exact Stage1041 raw candidates, no cherry-picked maxima.
"""
from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from skimage.feature import peak_local_max
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review');D=P/'data';O=P/'assets'/'stage1044_axial_spot_rescue';O.mkdir(exist_ok=True)
src=(Path(r'G:\Spatial_gene_site_publish\tools\build_stage1041_native_peak_recall_screen.py').read_text(encoding='utf-8-sig')).split('records=[]')[0]
scope={};exec(compile(src,'stage1041_funcs','exec'),scope)
metadata=scope['meta'];info=scope['info'];sel=scope['sel'];roi=scope['roi']
rows=[];long=[]
for gene in ['VGAT','VGLUT2','Snap25']:
 base,sp=info[gene]
 t=metadata[(metadata.region==1)&(metadata.gene==gene)].iloc[0]
 cz=float(t.source_native_Z_center)
 img,(x0,y0,z0)=roi(base,sel['cx']*8,sel['cy']*8,cz)
 hp=gaussian_filter(img,(.8,1,1))-gaussian_filter(img,(2,5,5))
 thr=float(np.quantile(hp[2:-2,12:-12,12:-12],.9975))
 peaks=peak_local_max(hp,min_distance=2,threshold_abs=thr,exclude_border=(2,12,12))
 origin=peaks.astype(float)+np.array([z0,y0,x0])[None,:]
 # all official source RSFISH spots within XY ROI, any z within ±80 native planes:
 matches=[]
 for ch in pd.read_csv(sp,header=None,usecols=[0,1,2],dtype=np.float32,chunksize=280000):
  xx=ch.iloc[:,0].to_numpy()/.23;yy=ch.iloc[:,1].to_numpy()/.23;zz=ch.iloc[:,2].to_numpy()/.42
  m=(xx>x0-6)&(xx<x0+606)&(yy>y0-6)&(yy<y0+606)&(abs(zz-cz)<82)
  if m.any():matches.append(np.column_stack([zz[m],yy[m],xx[m]]))
 spots=np.concatenate(matches) if matches else np.zeros((0,3))
 tree=cKDTree(spots[:,1:]) if len(spots) else None
 for jj,p in enumerate(origin):
  cand=tree.query_ball_point(p[1:],r=3.) if tree else []
  zdelta=np.abs(spots[cand,0]-p[0]) if len(cand) else np.array([])
  zmin=float(np.min(zdelta)) if len(zdelta) else np.inf
  # Original 1041 matching uses 3D metric X,Y 1px, Z*1.83, threshold4.5.
  c3=spots[cand]-p
  d3=np.sqrt((c3[:,0]*1.83)**2+(c3[:,1]**2)+(c3[:,2]**2)) if len(cand) else np.array([])
  original=bool(len(d3) and np.min(d3)<=4.5)
  long.append(dict(gene=gene,candidate_idx=jj,original_3d_match=original,nearest_xy3_zdiff=zmin,
      zbin=('matched_3d' if original else 'Zdiff_2to4' if zmin<=4 else 'Zdiff_4to8' if zmin<=8 else 'Zdiff_8to16' if zmin<=16 else 'Zdiff_16to32' if zmin<=32 else 'Zdiff_32to64' if zmin<=64 else 'not_matched_within_64')))
 M=pd.DataFrame([v for v in long if v['gene']==gene])
 a=M[~M.original_3d_match]
 row=dict(gene=gene,n_native_raw_candidates=len(M),n_matched_original=int(M.original_3d_match.sum()),
   n_notmatched_original=int(len(a)))
 for n in (4,8,16,32,64):
  row[f'n_unmatched_with_nearXY_and_abs_dZ_le_{n}']=int((a.nearest_xy3_zdiff<=n).sum())
  row[f'n_unmatched_noXYz_rescue_by_{n}']=int((a.nearest_xy3_zdiff>n).sum())
 rows.append(row)
 print('STAGE1044_AXIAL',row,flush=True)
T=pd.DataFrame(rows);T.to_csv(D/'stage1044_near_xy_upperlower_plane_rescue_counts.csv',index=False)
# full per candidate data kept private, not exposed on public
pd.DataFrame(long).to_csv(Path(r'G:\Map6_recover_all')/'stage1044_axial_candidates_private.csv',index=False)
fig,axs=plt.subplots(1,3,figsize=(14.8,4.8),layout='constrained')
for ax,gene in zip(axs,['VGAT','VGLUT2','Snap25']):
 row=T.set_index('gene').loc[gene]
 axis=[4,8,16,32,64]
 n=[row[f'n_unmatched_with_nearXY_and_abs_dZ_le_{n}'] for n in axis]
 ax.plot(axis,n,marker='o',color='#985aa4',lw=1.6)
 ax.axhline(row.n_notmatched_original,ls='--',color='#a4a4a4',lw=1)
 ax.set(xlabel='Upper/lower matching Z depth (native planes)',ylabel='Previously unmatched candidates with\nnearby XY official spot',title=f'{gene}: unmatched {row.n_notmatched_original} / {row.n_native_raw_candidates}')
 for sp in ('top','right'):ax.spines[sp].set_visible(False)
fig.suptitle('Does a missed-looking raw fluorescence peak actually have an RS-FISH detection in other Z planes?',fontsize=12)
fig.savefig(O/'STAGE1044_UNMATCHED_VGAT_UP_DOWN_Z_RESCUE.png',dpi=195)
plt.close(fig)
(D/'stage1044_axial_rescue_authority.json').write_text(json.dumps({
'stage':1044,'method':'Reproduces Stage1041 native raw 3D local high-pass candidate threshold 99.75%, same region1; for originally nonmatching candidates search complete official spot catalog at XY radius <=3 original 0.23µm pixels, across absolute up/down z difference 4/8/16/32/64 native 0.42µm planes. Distinct from original 3D match <=4.5 equivalent source pixels.',
'important':'A listed detected spot at similar XY but a distant Z could be another true molecule. XY-only rescues at large Z must not be treated as corrected 3D matches. Missing if no official spot nearby does not prove physical RNA.',
'origin':'One S500 region, sampled raw source N5 vs each corresponding RS-FISH official catalog. No biological replication.'},indent=2))
print('STAGE1044_COMPLETE',flush=True)

