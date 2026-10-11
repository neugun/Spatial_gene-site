"""Allen pons independent NE-gene validation. True marker reference Dbh+Th,
optional Slc18a2, never Slc6a2/NE01/Fine26 used to select reference.
Raw scRNAseq UMI by Pons 10xv3, NOT calibrated against FISH.
"""
from pathlib import Path
import pandas as pd,numpy as np,json
import anndata as ad
from scipy import sparse
ROOT=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review')
D=ROOT/'data';O=ROOT/'assets'/'stage1047_slc6a2_NE_audit';O.mkdir(exist_ok=True)
path=Path(r'Z:\sternsonlab\ABC_Atlas_Downloads\Pons\abc_pipeline_outputs\expression_matrices\WMB-10Xv3\20230630\WMB-10Xv3-P-raw.h5ad')
a=ad.read_h5ad(path,backed='r')
genes=['Slc6a2','Dbh','Th','Slc18a2','Slc32a1','Slc17a6','Snap25','Gad1','Gad2','Slc6a4']
sym=a.var['gene_symbol'].astype(str).to_numpy()
miss=[s for s in genes if s not in sym]
print('ALLEN_SHAPE',a.shape,'MISSING',miss,flush=True)
genes=[s for s in genes if s in sym]
ix=[int(np.flatnonzero(sym==s)[0]) for s in genes]
print('SYMBOL_TO_ENSEMBL',list(zip(genes,a.var_names[ix])),flush=True)
x=a[:,ix].X
if sparse.issparse(x):x=x.toarray()
M=np.asarray(x,dtype=np.float32)
a.file.close()
v={g:M[:,i] for i,g in enumerate(genes)}
n=len(M)
necore=(v['Dbh']>0)&(v['Th']>0)
ne_loose=(v['Dbh']>0)|(v['Th']>0)
ne_strict=necore&(v['Slc18a2']>0)
assert necore.sum()>0,(necore.sum(),n)
rows=[]
for source_name,m in [('ALL_PONS',np.ones(n,bool)),('Dbh_AND_Th',necore),('Dbh_OR_Th',ne_loose),('Dbh_Th_Slc18a2',ne_strict),
 ('Not_Dbh_Th',~necore),('Snap25_pos',v['Snap25']>0)]:
 for t in [0,1,2,3,5,10]:
  slc=v['Slc6a2']>t
  z=m&slc
  rows.append(dict(population=source_name,n=int(m.sum()),slc6a2_raw_UMI_cutoff=int(t),
     slc6a2_pos=int(z.sum()),slc6a2_pos_pct=float(100*z.sum()/max(1,m.sum())),
     n_joint_Th_Dbh=int((m&necore).sum()),
     frac_Slc6a2_pos_Dbh_Th=float(np.mean(necore[z])) if z.any() else np.nan,
     vgat_pos_among_selected=float(np.mean(v['Slc32a1'][z]>0)) if z.any() else np.nan,
     vglut2_pos_among_selected=float(np.mean(v['Slc17a6'][z]>0)) if z.any() else np.nan))
T=pd.DataFrame(rows);T.to_csv(D/'stage1047_allen_pons_NE_multigene_vs_Slc6a2_UMI.csv',index=False)
meta={'source':'Allen Pons WMB-10Xv3 20230630 RAW 143661 cells', 'independent_NE_reference':'Dbh+Th both >0 UMI; supplemental Slc18a2>0; independent of Slc6a2/Map10',
'raw_cells':n,'DbhTh_n':int(necore.sum()),'DbhThSlc18a2_n':int(ne_strict.sum()),'DbhOrTh_n':int(ne_loose.sum()),
'warning':'Dbh+Th scRNA UMI is itself dropout-prone; not all independently NE cells have both. Whole pons includes other catecholaminergic and multi-marker cells; no anatomical LC guarantee. Slc6a2 expression differs by modality, cannot directly apply UMI cutoff to FISH corrected spots.'}
(D/'stage1047_allen_pons_NE_authority.json').write_text(json.dumps(meta,indent=2))
print('INDEPENDENT_NE',meta,flush=True)
print(T[T.population.isin(['Dbh_AND_Th','Dbh_Th_Slc18a2','Not_Dbh_Th']) & T.slc6a2_raw_UMI_cutoff.isin([0,1])].to_string(index=False),flush=True)
