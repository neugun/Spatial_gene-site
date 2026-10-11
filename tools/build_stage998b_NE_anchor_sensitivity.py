from pathlib import Path
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(r'Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926')
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review')
Z=np.load(R/'stage704_current_umap27'/'CURRENT_UMAP27_STAGE704.npz')
gen=list(Z['genes'].astype(str))
C=np.concatenate([np.load(Path(r'G:\Map6_recover_all\CURRENT_ROUTEA')/f'S{s}_cell_gene27_current.npy') for s in (500,530,560)],axis=0)
A=pd.read_csv(R/'stage631_current_3d_atlas'/'CURRENT_3D_CELL_ATLAS.csv.gz',usecols=['section','fine26_stable_id'])
get=lambda g:C[:,gen.index(g)]
snap=get('Snap25')>5
ne0=A.fine26_stable_id.eq('NE01').to_numpy()
ch0=A.fine26_stable_id.eq('ChAT01').to_numpy()
n=get('Slc6a2');h=get('Slc5a7');vg=get('Slc32a1');gl=get('Slc17a6')
rows=[]
for s in [0,500,530,560]:
 m=(A.section.to_numpy()==s) if s else np.ones(len(C),bool)
 for t in [0,5,10,20,40,80,120,200,400,800,1200]:
  mask=m&snap&(n>t)
  rows.append(dict(section=s,marker='Slc6a2',t=t,n=int(mask.sum()),fraction_frozen_NE=float(np.mean(ne0[mask])) if mask.any() else np.nan,
      sensitivity_to_frozen_NE=float(np.sum(mask&ne0)/max(1,np.sum(m&snap&ne0))),
      count_frozen_NE=int(np.sum(mask&ne0)),total_frozen_NE=int(np.sum(m&snap&ne0))))
 for t in [5,20,50,100,200,400]:
  mask=m&snap&(h>t)
  rows.append(dict(section=s,marker='Slc5a7',t=t,n=int(mask.sum()),fraction_frozen_NE=np.nan,
   sensitivity_to_frozen_NE=np.nan,count_frozen_NE=0,total_frozen_NE=0))
D=pd.DataFrame(rows)
D.to_csv(P/'data'/'stage998b_NE_anchor_sensitivity_cutoff_scan.csv',index=False)
options=[]
for rule in ['frozen_only','Slc6a2>5','Slc6a2>20','Slc6a2>80','Slc6a2>200']:
 if rule=='frozen_only':NE=ne0
 else:NE=ne0|(n>float(rule.split('>')[1]))
 CH=ch0|(h>50)
 elig=snap&~NE&~CH
 for ss in [0,500,530,560]:
  m=elig&((A.section.to_numpy()==ss) if ss else np.ones(len(C),bool))
  for tg,te in [(10,5),(15,5),(20,5),(25,8),(30,10),(40,15)]:
   z=(vg>tg)&(gl>te)&m
   either=((vg>tg)|(gl>te))&m
   options.append(dict(rule_NE=rule,section=ss,t_vgat=tg,t_vglut2=te,eligible=int(m.sum()),
      both=int(z.sum()),both_pct=float(100*z.sum()/max(1,m.sum())),
      both_pct_of_union=float(100*z.sum()/max(1,either.sum())),
      neither_pct=float(100*np.mean(~((vg[m]>tg)|(gl[m]>te))) if m.any() else np.nan),
      ne_excluded=int(np.sum(snap&NE))))
T=pd.DataFrame(options);T.to_csv(P/'data'/'stage998b_ne_anchor_selection_robustness.csv',index=False)
print('NE_ROC',D[(D.section==0)&(D.marker=='Slc6a2')].to_string(index=False,float_format=lambda q:f'{q:.3f}'),flush=True)
print('THRESHOLD_ANCHOR',T[(T.section==0)&(T.t_vgat==15)&(T.t_vglut2==5)].to_string(index=False),flush=True)
fig,axs=plt.subplots(1,2,figsize=(10.2,4.3))
j=D[(D.section==0)&(D.marker=='Slc6a2')&(D.t>0)]
axs[0].plot(j.t,j.fraction_frozen_NE*100,marker='o',lw=1.8,color='#62a6bd')
axs[0].set_xscale('log')
axs[0].set(ylabel='% high-Slc6a2 cells belonging to frozen NE01',xlabel='Slc6a2 corrected count cutoff',title='Slc6a2 >5 is not NE-specific')
for rule,co in zip(['frozen_only','Slc6a2>5','Slc6a2>20','Slc6a2>80','Slc6a2>200'],
                   ['#555555','#ca765c','#cfac73','#559bb7','#7c72b5']):
 r=T[(T.section==0)&(T.rule_NE==rule)]
 axs[1].plot(r.t_vgat,r.both_pct,marker='o',lw=1.4,label=rule,color=co)
axs[1].set(xlabel='VGAT cutoff (VGLUT2 cutoff co-varies by predeclared candidates)',ylabel='Double+ fraction (%)',title='NE anchor choice changes eligible denominator')
axs[1].legend(frameon=False,fontsize=7)
for a in axs:
 for sp in ['top','right']:a.spines[sp].set_visible(False)
fig.tight_layout()
out=P/'assets'/'stage998_ne_reference_threshold'/'STAGE998b_NE_ANCHOR_SENSITIVITY.png'
fig.savefig(out,dpi=190);plt.close(fig)
