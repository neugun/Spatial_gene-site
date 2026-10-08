from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(r'Z:\sternsonlab\Zhenggang\2acq\map6_allsections_slurm_20260926')
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review');D=P/'data';O=P/'assets'/'stage1043_strict_transmitter_gate';O.mkdir(exist_ok=True)
A=pd.read_csv(R/'stage631_current_3d_atlas'/'CURRENT_3D_CELL_ATLAS.csv.gz',usecols=['section','roi_id','fine26_stable_id','x_um','y_um'])
G=np.load(R/'stage704_current_umap27'/'CURRENT_UMAP27_STAGE704.npz');genes=list(G['genes'].astype(str))
C=np.concatenate([np.load(Path(r'G:\Map6_recover_all\CURRENT_ROUTEA')/f'S{s}_cell_gene27_current.npy') for s in (500,530,560)])
F=pd.read_csv(Path(r'G:\Map6_recover_all\stage1001_neuron_gate_cell_flags_PRIVATE.csv.gz'),usecols=['section','roi_id','reference_neurons'])
assert len(A)==len(C)==len(F)==71950 and np.array_equal(A.roi_id.to_numpy(),F.roi_id.to_numpy())
v=lambda name:C[:,genes.index(name)]
vg=v('Slc32a1');gl=v('Slc17a6');snap=v('Snap25')
mapping=pd.read_csv(D/'stage977_fine26_color_key.csv')
broad=A.fine26_stable_id.map(dict(zip(mapping.fine26,mapping.broad))).fillna('Other').to_numpy()
nosp=~np.isin(broad,['NE','ChAT'])
vgp=vg>10;glp=gl>5;sp=snap>5
u=sp&nosp&(vgp|glp);ref=F.reference_neurons.to_numpy(bool)
cl=np.full(len(C),-1,np.int8)
cl[u&vgp&~glp]=0;cl[u&~vgp&glp]=1;cl[u&vgp&glp]=2
sec=A.section.to_numpy(int)
records=[]
for s in (0,500,530,560):
 m=(sec==s) if s else np.ones(len(C),bool)
 for nm,gate in [('all',m),('Snap25>5',m&sp),('Snap25>5_noNEChat',m&sp&nosp),('user_exact',m&u),('user_exact_plus_refNeuron',m&u&ref)]:
  nn=gate.sum();numbers=[int(np.sum(gate&(cl==k))) for k in (0,1,2)]
  records.append(dict(section=s,gate=nm,n=int(nn),VGATonly=numbers[0],GLUonly=numbers[1],double=numbers[2],dual_pct=100*numbers[2]/nn if nn else np.nan,
   retained_pct_of_section=100*nn/m.sum(),n_ref_neurons=int((gate&ref).sum()),Fine_E=int((gate&(broad=='E')).sum()),Fine_I=int((gate&(broad=='I')).sum())))
S=pd.DataFrame(records);S.to_csv(D/'stage1043_exact_gate_counts.csv',index=False)
print('EXACT_GATE_RESULTS',S.to_string(index=False,float_format=lambda x:f'{x:.2f}'),flush=True)
extras=[]
for cut in (0,1,3,5,10):
 alt=sp&(vgp|glp)&(v('Slc6a2')<=cut)&(v('Slc5a7')<=cut)
 extras.append(dict(NE_ChAT_marker_max=cut,n=int(alt.sum()),double=int((alt&vgp&glp).sum()),dual_pct=100*(alt&vgp&glp).sum()/max(1,alt.sum()),frozen_NE_ChAT_remaining=int((alt&~nosp).sum())))
pd.DataFrame(extras).to_csv(D/'stage1043_special_marker_exclusion_sensitivity.csv',index=False)
t=pd.crosstab(A.fine26_stable_id,cl);t.columns=[f'class_{a}' for a in t.columns];t.to_csv(D/'stage1043_fine26_by_3classes.csv')
(D/'stage1043_strict_gate_authority.json').write_text(json.dumps({'criteria':'Corrected Snap25>5 AND frozen Broad not NE/ChAT AND (VGAT Slc32a1>10 OR VGLUT2 Slc17a6>5)','universe':71950,'n_selected':int(u.sum()),'n_dual':int((cl==2).sum()),
'warning':'Double-positive remains a third category. A gate cannot eliminate GABA/glutamate overlap and can increase percent double among retained cells. Broad NE/ChAT labels are from same 27-gene panel; Slc6a2 and Slc5a7 sensitivity table provided. No Chat gene in panel. Reuse existing Stage995 full-data tSNE, NOT refit.'},indent=2))
np.savez_compressed(Path(r'G:\Map6_recover_all')/'stage1043_private_filter_masks.npz',selected=u,ref=ref,class3=cl)
print('STAGE1043_COMPUTED',flush=True)
