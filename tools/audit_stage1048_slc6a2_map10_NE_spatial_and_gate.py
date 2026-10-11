"""Audit Slc6a2 as NET gene, not an independently verified LC-NE classification.
Explicit NE01 double positive, section-dependent confounds, NE high cell spatial distribution.
Do not alter cell IDs or fit threshold to expected GABA/Glu proportions.
"""
from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(r'Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926')
ROOT=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review')
D=ROOT/'data';O=ROOT/'assets'/'stage1048_Slc6a2_marker_reaudit';O.mkdir(exist_ok=True,parents=True)
g=list(np.load(R/'stage704_current_umap27'/'CURRENT_UMAP27_STAGE704.npz')['genes'].astype(str))
C=np.concatenate([np.load(Path(r'G:\Map6_recover_all\CURRENT_ROUTEA')/f'S{s}_cell_gene27_current.npy') for s in (500,530,560)],axis=0)
A=pd.read_csv(R/'stage631_current_3d_atlas'/'CURRENT_3D_CELL_ATLAS.csv.gz',usecols=['section','fine26_stable_id','x_um','y_um','roi_id'])
assert len(A)==len(C)==71950
N=C[:,g.index('Slc6a2')];S=C[:,g.index('Snap25')]
V=C[:,g.index('Slc32a1')];E=C[:,g.index('Slc17a6')]
label=A.fine26_stable_id.to_numpy()=='NE01'
neuron=S>5
rows=[]
for ss in [0,500,530,560]:
 section=(A.section.to_numpy()==ss) if ss else np.ones(len(A),bool)
 for pop,mask in [('ALL',section),('SNAP25_GT5',section&neuron),
                  ('Frozen_NE01',section&label),('Frozen_NE01_snap_pos',section&label&neuron),
                  ('Outside_NE01_snap_pos',section&~label&neuron)]:
  vals=N[mask]
  rows.append(dict(section=int(ss),pop=pop,n=int(mask.sum()),slc6a2_median=float(np.median(vals)) if len(vals) else None,
   slc6a2_q90=float(np.quantile(vals,.9)) if len(vals) else None,
   pct_gt5=float((vals>5).mean()*100) if len(vals) else None,
   pct_gt80=float((vals>80).mean()*100) if len(vals) else None,
   pct_gt200=float((vals>200).mean()*100) if len(vals) else None,
   vgat_gt10_pct=float((V[mask]>10).mean()*100) if len(vals) else None,
   vglut_gt5_pct=float((E[mask]>5).mean()*100) if len(vals) else None))
 pd.DataFrame(rows).to_csv(D/'stage1048_Slc6a2_NE01_continuous_not_truth.csv',index=False)
T=pd.DataFrame(rows)
# full physical 3 section reference; not source LC anatomical probability
fig,axs=plt.subplots(2,3,figsize=(14.2,9.4))
pal=['#a7d7e1','#edba92','#a13746']
for i,ss in enumerate([500,530,560]):
 mask=A.section.to_numpy()==ss
 x=A.x_um.to_numpy()[mask];y=A.y_um.to_numpy()[mask]
 count=N[mask];labels=label[mask];snap=neuron[mask]
 axs[0,i].scatter(x,y,s=.6,color='#dee5e8',alpha=.22,linewidths=0,rasterized=True)
 sel=snap&(count>5)
 axs[0,i].scatter(x[sel],y[sel],s=.95,color='#97d4df',alpha=.48,linewidths=0,rasterized=True)
 sel2=snap&(count>200)
 axs[0,i].scatter(x[sel2],y[sel2],s=2.5,color='#b54755',alpha=.78,linewidths=0,rasterized=True)
 axs[0,i].set(title=f'S{ss} Slc6a2>5 cyan, >200 red',xlabel='Current physical X (um)',ylabel='Current physical Y (um)')
 axs[1,i].scatter(x,y,s=.7,color='#dedfe2',alpha=.2,linewidths=0,rasterized=True)
 axs[1,i].scatter(x[labels],y[labels],s=2.7,color='#9a65a9',alpha=.85,linewidths=0,rasterized=True)
 axs[1,i].set(title=f'S{ss} existing frozen NE01 labels, n={int(labels.sum())}',xlabel='Current physical X (um)',ylabel='Current physical Y (um)')
 for ax in [axs[0,i],axs[1,i]]:
  ax.set_aspect('equal');ax.invert_yaxis()
  for sp in ['top','right']:ax.spines[sp].set_visible(False)
fig.suptitle('SLC6A2 is a biological NE-associated transporter but its FISH threshold is not an independently verified NE cell label\n'
 'Current 27-gene NE01 reference used Slc6a2; source counts are corrected count-equivalents not integer molecules',fontsize=10.8)
fig.tight_layout();fig.savefig(O/'STAGE1048_SLC6A2_GT5_GT200_VS_FROZEN_NE01_FULL3SECTION_XY.png',dpi=175);plt.close(fig)
# distributions show clearly high off-target at low cutoff; uses labels only as internal circularity evidence
fig,axs=plt.subplots(1,2,figsize=(11.3,4.2))
bins=np.linspace(0,4.2,85)
for cat,mm,clr in [('Frozen NE01',label&neuron,'#a56096'),('Outside NE01',~label&neuron,'#64a8ba')]:
 vals=np.log10(1+N[mm])
 axs[0].hist(vals,bins=bins,histtype='step',density=True,lw=1.7,color=clr,label=f'{cat} (n={mm.sum():,})')
axs[0].set(xlabel='log10(1 + corrected Slc6a2 counts)',ylabel='Cell density',title='SLC6A2 gradient and 27-gene class reuse')
axs[0].legend(frameon=False,fontsize=8)
cuts=[0,5,10,20,40,80,120,200,400,800]
y=[int(np.sum(neuron&(N>t))) for t in cuts]
pref=[float(np.mean(label[neuron&(N>t)])) for t in cuts]
axs[1].plot(cuts,y,marker='o',color='#65aaba')
axs[1].set(xlabel='Corrected Slc6a2 cutoff',ylabel='Number of Snap25>5 cells retained',title='Marker threshold excludes non-NE01 neurons at low cutoff')
ax2=axs[1].twinx();ax2.plot(cuts,np.array(pref)*100,marker='s',color='#a96686')
ax2.set_ylabel('% agreeing with pre-existing (circular) NE01',color='#a96686')
for ax in axs:
 for sp in ['top','right']:ax.spines[sp].set_visible(False)
fig.tight_layout();fig.savefig(O/'STAGE1048_SLC6A2_THRESHOLD_NOT_AN_INDEPENDENT_NE_TEST.png',dpi=185);plt.close(fig)
K=pd.read_csv(D/'stage1009_mapmycells_fine26_class_fraction_crosswalk.csv').set_index('fine26_stable_id')
k=K.loc['NE01'].sort_values(ascending=False).head(12)
k.to_csv(D/'stage1048_NE01_OFFICIAL_MAPMYCELLS_CLASS_MIX.csv',header=['fraction'])
B=pd.read_csv(D/'stage1009_official_fine26_confidence_and_reportability.csv')
b=B[B.Fine26=='NE01'].iloc[0]
auth={'status':'RETRACT_SLC6A2_AS_VALIDATED_NE_GATE','marker':'Slc6a2 encodes NET. Marker positivity not proven cell identity',
 'provisional_gene_only_threshold_200':'NOT RECOMMENDED as NE ground-truth or GABA/Glu calibration target',
 'frozen_NE01_source':'same targeted 27 FISH genes, especially Slc6a2, circular reference',
 'mapmycells_NE01_modal':str(k.index[0]),'mapmycells_NE01_modal_fraction':float(k.iloc[0]),
 'mapmycells_NE01_class_bootstrap_median':float(b['median_class_bootstrap']),
 'mapmycells_NE01_ge50_pct':float(b['pct_class_bootstrap_ge_50']),
 'corrected_RouteA_counts_are_not_physical_spot_integer_threshold':True,
 'independent_confirmation_needed':['Dbh','Th','Slc18a2','original Slc6a2 native N5 raw puncta','LC atlas registration'],
 'never_force_NEVGAT_or_NEVGLUT_to_scRNA_1of36':True}
(D/'stage1048_NE01_gate_correction_authority.json').write_text(json.dumps(auth,indent=2))
print('STAGE1048_OUTPUT',auth,flush=True)
print('SLC6A2_POP',T[(T.section==0)].to_string(index=False),flush=True)
