"""Compare users' independent scRNA25 NE 36-cell sparse transcript detections to
Map10 FISH frozen NE labels. Prevent target leakage: explicitly show what quantile
forcing scRNA rate would cost in transmitter coverage, do NOT recommend forcing.
"""
from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportion_confint
R=Path(r'Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926')
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review')
g=list(np.load(R/'stage704_current_umap27'/'CURRENT_UMAP27_STAGE704.npz')['genes'].astype(str))
C=np.concatenate([np.load(Path(r'G:\Map6_recover_all\CURRENT_ROUTEA')/f'S{s}_cell_gene27_current.npy') for s in (500,530,560)],axis=0)
lab=pd.read_csv(R/'stage631_current_3d_atlas'/'CURRENT_3D_CELL_ATLAS.csv.gz',usecols=['fine26_stable_id','section'])
SCR=Path(r'D:\cluster\Desktop\EASIFISH\PeriLC_VTA_work\projection_scrna_overlay_v136\scRNA25_projection_gene_expression_v136.csv')
X=pd.read_csv(SCR)
ne=X[(X.broad_class=='NE')&(X.gene.isin(['Slc17a6','Slc32a1']))].copy()
assert len(ne)==2 and ne.n.eq(36).all()
refmask=lab.fine26_stable_id.eq('NE01').to_numpy()
neuron=C[:,g.index('Snap25')]>5
refmask=refmask&neuron
v=C[:,g.index('Slc32a1')]
e=C[:,g.index('Slc17a6')]
h=C[:,g.index('Slc5a7')]
# Additional marker NE>200 sensitivity
NE=refmask|(C[:,g.index('Slc6a2')]>200)
CHAT=lab.fine26_stable_id.eq('ChAT01').to_numpy()|(h>50)
eligible=neuron&~NE&~CHAT
rows=[]
for rr in ne.itertuples():
 gene=rr.gene;Y=v if gene=='Slc32a1' else e
 det=int(round(rr.detect_fraction*rr.n))
 ci=proportion_confint(det,int(rr.n),alpha=.05,method='wilson')
 # forcing calibrated to 1/36 looks like extreme FISH threshold because different assay dropout
 t=float(np.quantile(Y[refmask],1-rr.detect_fraction))
 for tt in [3,5,10,15,20,30,40,60,80,100,150,200,t]:
  n_retain=int(np.sum(eligible&(Y>tt)))
  rows.append(dict(gene=gene,threshold_raw_corrected=float(tt),
   scRNA_NE_n=int(rr.n),scRNA_NE_detected=det,scRNA_NE_detection_pct=float(100*rr.detect_fraction),
   scRNA_NE_Wilson95_lo=float(ci[0]*100),scRNA_NE_Wilson95_hi=float(ci[1]*100),
   observed_FISH_frozen_NE_positivity_pct=float(100*np.mean(Y[refmask]>tt)),
   eligible_non_NE_neuron_retained_pct=float(100*n_retain/max(1,eligible.sum())),
   n_eligible_non_NE_neuron_retained=n_retain,
   threshold_that_would_force_scRNA_pct=t,
   label_source='NE01 class based partly on same 27 FISH gene panel; not independent gold standard'))
T=pd.DataFrame(rows)
T.to_csv(P/'data'/'stage998c_independent_NE_scRNA36_vs_Map10_FISH_threshold_QC.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(10.2,4.4))
for ax,(gene,clr) in zip(axs,[('Slc32a1','#648fb0'),('Slc17a6','#cf8468')]):
 d=T[(T.gene==gene)].sort_values('threshold_raw_corrected')
 a=d.iloc[0]
 ax.plot(d.threshold_raw_corrected,d.observed_FISH_frozen_NE_positivity_pct,marker='o',label='NE01 FISH positivity',color=clr)
 ax.plot(d.threshold_raw_corrected,d.eligible_non_NE_neuron_retained_pct,marker='s',label='Candidate non-NE retained',color='#777d85')
 ax.axhspan(a.scRNA_NE_Wilson95_lo,a.scRNA_NE_Wilson95_hi,color='#93b4bf',alpha=.17,label='scRNA NE 95% CI (1/36)')
 ax.axhline(a.scRNA_NE_detection_pct,color='#487b94',ls='--',lw=.9)
 ax.set(xlabel='Corrected spot-count threshold',ylabel='% of respective cells',title=f'{gene}: no assay-independent threshold identity')
 ax.set_xscale('symlog',linthresh=10);ax.set_ylim(0,105);ax.legend(frameon=False,fontsize=7.2)
 for side in ['right','top']:ax.spines[side].set_visible(False)
fig.suptitle('Independent scRNA NE cell n=36 is too small and modality mismatched to FORCE an FISH cutoff')
fig.tight_layout();fig.savefig(P/'assets'/'stage998_ne_reference_threshold'/'STAGE998c_NE_SCSEQ36_REFERENCE_CAVEAT.png',dpi=180);plt.close(fig)
print('NE_SCRNA_36_REFERENCE',T.groupby('gene').first()[['scRNA_NE_detection_pct','scRNA_NE_Wilson95_lo','scRNA_NE_Wilson95_hi','threshold_that_would_force_scRNA_pct']].to_string(),flush=True)
