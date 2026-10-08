from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review');D=P/'data';O=P/'assets'/'stage1043_strict_transmitter_gate'
R=Path(r'Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926')
z=np.load(R/'stage704_current_umap27'/'CURRENT_UMAP27_STAGE704.npz');g=list(z['genes'].astype(str))
X=np.concatenate([np.load(Path(r'G:\Map6_recover_all\CURRENT_ROUTEA')/f'S{s}_cell_gene27_current.npy') for s in (500,530,560)])
a=pd.read_csv(R/'stage631_current_3d_atlas'/'CURRENT_3D_CELL_ATLAS.csv.gz',usecols=['fine26_stable_id'])
key=pd.read_csv(D/'stage977_fine26_color_key.csv');b=a.fine26_stable_id.map(dict(zip(key.fine26,key.broad))).to_numpy()
v=X[:,g.index('Slc32a1')];w=X[:,g.index('Slc17a6')];s=X[:,g.index('Snap25')]
base=(s>5)&~np.isin(b,['NE','ChAT'])
rows=[]
for vt in (5,10,15,20,30,50):
 for wt in (2,3,5,8,10,15,20):
  V=base&(v>vt);W=base&(w>wt);any_=V|W
  rows.append(dict(VGAT_cutoff=vt,VGLUT2_cutoff=wt,n=int(any_.sum()),double=int((V&W).sum()),pct_double=100*(V&W).sum()/max(any_.sum(),1)))
q=pd.DataFrame(rows);q.to_csv(D/'stage1045_exact_gate_threshold_sensitivity.csv',index=False)
print(q.query('VGAT_cutoff in [10,20,30] and VGLUT2_cutoff in [5,10,15]').to_string(index=False),flush=True)
m=base&(v>10)&(w>5)
ratio=(v[m]+1)/(w[m]+1)
out=dict(n_dual=int(m.sum()),n_both_strong_20_10=int(np.sum((v[m]>20)&(w[m]>10))),n_both_strong_30_15=int(np.sum((v[m]>30)&(w[m]>15))),
 n_raw_vgat_dominant_3x=int((ratio>=3).sum()),n_raw_vglut_dominant_3x=int((ratio<=1/3).sum()),n_ambiguous_ratio=int(((ratio>1/3)&(ratio<3)).sum()),probe_efficiency_warning='Raw count ratios not comparable across genes')
(D/'stage1045_dual_strength.json').write_text(json.dumps(out,indent=2));print('DUAL_STRENGTH',out,flush=True)
fig,axs=plt.subplots(1,2,figsize=(11.6,4.7),layout='constrained')
for ax,f,title,cm in [(axs[0],'pct_double','Double-positive among retained (%)','PuRd'),(axs[1],'n','Number of selected cells','Blues')]:
 t=q.pivot(index='VGAT_cutoff',columns='VGLUT2_cutoff',values=f)
 im=ax.imshow(t,cmap=cm,aspect='auto')
 ax.set_xticks(range(len(t.columns)),t.columns);ax.set_yticks(range(len(t.index)),t.index)
 ax.set(xlabel='VGLUT2 count >',ylabel='VGAT count >',title=title)
 fig.colorbar(im,ax=ax,fraction=.045)
fig.savefig(O/'STAGE1045_GATE_THRESHOLD_SENSITIVITY.png',dpi=195)
print('STAGE1045_COMPLETE',flush=True)
