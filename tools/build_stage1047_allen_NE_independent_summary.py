"""Stage1047 independent Dbh/Th Allen Pons NE-marker benchmark; no FISH cutoff calibration."""
from pathlib import Path
import pandas as pd,numpy as np,json
from statsmodels.stats.proportion import proportion_confint
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review')
T=pd.read_csv(P/'data'/'stage1047_allen_pons_NE_multigene_vs_Slc6a2_UMI.csv')
r=T[(T.slc6a2_raw_UMI_cutoff==0)&(T.population.isin(['Dbh_AND_Th','Not_Dbh_Th']))]
assert len(r)==2
r=r.set_index('population')
fig,ax=plt.subplots(figsize=(6,5.1))
xx=np.arange(2)
v=[r.loc['Dbh_AND_Th','slc6a2_pos_pct'],r.loc['Not_Dbh_Th','slc6a2_pos_pct']]
n=[int(r.loc['Dbh_AND_Th','n']),int(r.loc['Not_Dbh_Th','n'])]
k=[int(r.loc['Dbh_AND_Th','slc6a2_pos']),int(r.loc['Not_Dbh_Th','slc6a2_pos'])]
ax.bar(xx,v,width=.40,color=['#c96a6b','#8ac9d3'],edgecolor='none')
for i,(k0,n0) in enumerate(zip(k,n)):
 ci=proportion_confint(k0,n0,method='wilson')
 ax.errorbar(i,100*k0/n0,yerr=np.array([[100*k0/n0-100*ci[0]],[100*ci[1]-100*k0/n0]]),fmt='none',elinewidth=1.0,capsize=3,color='#444444')
 ax.text(i,v[i]+4.8,f'{k0}/{n0:,}\n{v[i]:.2f}%',ha='center',fontsize=10)
ax.set(xlim=(-.55,1.55),ylim=(0,95),ylabel='Raw Slc6a2 UMI detected in reference (%)',
 title='Allen Pons: Dbh+Th reference independently predicts NET expression')
ax.set_xticks(xx,['Dbh + Th\n(26 cells)','Not Dbh+Th\n(143,635 cells)'])
for side in ['top','right']:ax.spines[side].set_visible(False)
fig.text(.05,.016,'Different modality from Map10 FISH. Dropout and rare joint-Dhb/Th positives; small n. No LC-only anatomy or per-mouse replication.',fontsize=8)
fig.subplots_adjust(bottom=.20,top=.89,left=.13,right=.96)
out=P/'assets'/'stage1047_slc6a2_NE_audit';out.mkdir(exist_ok=True)
fig.savefig(out/'STAGE1047_ALLEN_PONS_DBH_TH_INDEPENDENT_NE_VS_SLC6A2.png',dpi=200);plt.close(fig)
s0=r.loc['Dbh_AND_Th']
result={'allen_pons_cells':143661,
 'NE_DbhANDTh_n':int(s0.n),'NE_DbhANDTh_NET_detected':int(s0.slc6a2_pos),
 'NET_detect_in_DbhTh_pct':float(s0.slc6a2_pos_pct),
 'NET_detect_without_joint_DbhTh_pct':float(r.loc['Not_Dbh_Th','slc6a2_pos_pct']),
 'VGAT_detect_with_DbhTh_and_NET_pct':float(100*s0.vgat_pos_among_selected),
 'VGLUT2_detect_with_DbhTh_and_NET_pct':float(100*s0.vglut2_pos_among_selected),
 'do_not_force_target_rate':True,'single_cell_UMI_vs_FISH_not_comparable':True,
 'important':'Dbh+Th raw UMI selected only 26 cells of 143661 total; absence of joint call includes NE dropout. 60% VGLUT2 among 20 positive joint DbhTh NET shows zero VGLUT2 overlap is unjustified but not a stable LC rate.'}
(P/'data'/'stage1047_allen_NE_independent_summary.json').write_text(json.dumps(result,indent=2))
print('ALLEN_INDEPENDENT_NE_SUMMARY',result,flush=True)
