from pathlib import Path
import json, math, shutil
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter, label as cc_label, binary_dilation
from sklearn.cluster import KMeans
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from PIL import Image
from map10_display_xy import flip_xy,verify_flip_xy
from map10_display_xy import flip_xy,verify_flip_xy

PUB=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
DATA=PUB/"data"
OLD=PUB/"assets"/"stage947_final_atlas"
OUT=PUB/"assets"/"stage956_final_atlas_v2"
OUT.mkdir(parents=True,exist_ok=True)

X=np.load(DATA/"stage943_routeA_all27.npy").astype(float)
cells=pd.read_csv(DATA/"stage947_sectionwise_region10_display_labels.csv.gz")
oldman=pd.read_csv(DATA/"stage947_final_gene_section_manifest.csv")
genes=oldman["gene"].drop_duplicates().astype(str).tolist()
assert X.shape==(len(cells),len(genes))
secs=[500,530,560]
x0=cells["x"].to_numpy(float); y0=cells["y"].to_numpy(float); section=cells["section"].to_numpy(int)

# Canonical display convention: x flipped, y flipped within each section.
xd=np.empty_like(x0); yd=np.empty_like(y0)
for ss in secs:
    m=section==ss
    xd[m],yd[m]=flip_xy(x0[m],y0[m])
    verify_flip_xy(x0[m],y0[m],xd[m],yd[m])

required=["Hcrtr1","Bcl11b","Slc5a7","Lmx1a","Piezo2","Slc6a2","Ghr"]
# previously useful spatially selective genes; broad housekeeping-like genes are not allowed to dominate region fitting.
specific_pool=["Npy1r","Pdyn","Tacr1","Slc6a4","Foxp2","Fstl4","Epb41l4a","Pde11a","Ebf3","Cep112","Otx2os1","Dscaml1","Mc4r","Meis2","Zfhx3","Npy5r"]
gidx={g:i for i,g in enumerate(genes)}

def grid_fields(ss, marker_genes, n=165, smooth=3.0):
    m=section==ss
    xs,ys=xd[m],yd[m]
    xe=np.linspace(xs.min(),xs.max(),n+1); ye=np.linspace(ys.min(),ys.max(),n+1)
    count,_,_=np.histogram2d(ys,xs,bins=[ye,xe])
    dens=gaussian_filter(count,1.5)
    mask=dens>max(0.5,0.012*dens.max())
    feats=[]
    for g in marker_genes:
        v=np.log1p(X[m,gidx[g]])
        sm,_,_=np.histogram2d(ys,xs,bins=[ye,xe],weights=v)
        num=gaussian_filter(sm,smooth)
        den=gaussian_filter(count,smooth)
        field=num/np.maximum(den,1e-6)
        vals=field[mask]
        field=(field-np.nanmedian(vals))/(np.nanstd(vals)+1e-6)
        feats.append(field)
    xc=(xe[:-1]+xe[1:])/2; yc=(ye[:-1]+ye[1:])/2
    return np.stack(feats), mask, xe, ye, xc, yc, count

def spatial_score_for_gene(ss,g,n=90):
    # Score on a smoothed grid: spatial field variance × occupied coverage, normalized by cell-level variance.
    m=section==ss
    xs,ys=xd[m],yd[m];v=np.log1p(X[m,gidx[g]])
    xe=np.linspace(xs.min(),xs.max(),n+1);ye=np.linspace(ys.min(),ys.max(),n+1)
    cnt,_,_=np.histogram2d(ys,xs,bins=[ye,xe])
    sm,_,_=np.histogram2d(ys,xs,bins=[ye,xe],weights=v)
    field=gaussian_filter(sm,2.2)/np.maximum(gaussian_filter(cnt,2.2),1e-6)
    mask=gaussian_filter(cnt,1.4)>max(.5,.015*gaussian_filter(cnt,1.4).max())
    spatial=np.var(field[mask]) if mask.any() else 0
    cell=np.var(v)+1e-9
    pos=(X[m,gidx[g]]>0).mean()
    return float(spatial/cell*np.sqrt(max(pos,.02)))

def separation_score(F, lab, mask):
    vals=F[:,mask].T
    overall=np.var(vals,axis=0).sum()+1e-9
    between=0.0
    mu=vals.mean(0)
    for k in np.unique(lab[mask]):
        q=lab[mask]==k
        if q.sum():
            d=vals[q].mean(0)-mu
            between+=q.mean()*float(np.dot(d,d))
    return float(between/overall)

def boundary_complexity(lab,mask):
    # normalized edge disagreement; lower = smoother
    a=(lab[:,1:]!=lab[:,:-1]) & mask[:,1:] & mask[:,:-1]
    b=(lab[1:,:]!=lab[:-1,:]) & mask[1:,:] & mask[:-1,:]
    denom=(mask[:,1:]&mask[:,:-1]).sum()+(mask[1:,:]&mask[:-1,:]).sum()
    return float((a.sum()+b.sum())/max(denom,1))


def enforce_contiguous(lab, mask, prob_stack, max_iter=8):
    """Keep each label as one dominant component; reassign islands to adjacent labels.
    This is display cleanup only and is constrained by the marker-derived probability field."""
    out=lab.copy()
    for _ in range(max_iter):
        changed=False
        for k in range(10):
            comp,n=cc_label((out==k)&mask)
            if n<=1:
                continue
            sizes=np.bincount(comp.ravel())[1:]
            keep=1+int(np.argmax(sizes))
            for ci in range(1,n+1):
                if ci==keep:
                    continue
                island=comp==ci
                border=binary_dilation(island,iterations=1)&mask&(~island)
                neigh=out[border]
                neigh=neigh[(neigh>=0)&(neigh!=k)]
                if len(neigh)==0:
                    continue
                cand=np.unique(neigh)
                # boundary contact first, marker probability second
                best=None
                for j in cand:
                    contact=float(np.mean(neigh==j))
                    pscore=float(np.mean(prob_stack[int(j)][island]))
                    score=contact+0.35*pscore
                    if best is None or score>best[0]:
                        best=(score,int(j))
                if best is not None:
                    out[island]=best[1]
                    changed=True
        if not changed:
            break
    return out

def components_count(lab,mask):
    total=0
    for k in range(10):
        _,n=cc_label((lab==k)&mask)
        total+=n
    return int(total)

region_info={}
region_manifest=[]
sweep_rows=[]
cell_region10=np.zeros(len(cells),dtype=int)

for ss in secs:
    ranked=sorted([(spatial_score_for_gene(ss,g),g) for g in specific_pool], reverse=True)
    auto=[g for _,g in ranked[:5] if g not in required]
    markers=required+auto
    # marker fields are deliberately smooth before clustering: regions reflect broad anatomical domains, not cell-level speckle.
    F,mask,xe,ye,xc,yc,cnt=grid_fields(ss,markers,n=165,smooth=3.4 if ss!=560 else 4.2)
    YY,XX=np.meshgrid(yc,xc,indexing="ij")
    # Geometry is a supporting prior only; markers dominate.
    xz=(XX-XX[mask].mean())/(XX[mask].std()+1e-6)
    yz=(YY-YY[mask].mean())/(YY[mask].std()+1e-6)
    V=np.concatenate([F[:,mask].T,0.42*xz[mask,None],0.42*yz[mask,None]],axis=1)
    km=KMeans(n_clusters=10,n_init=30,random_state=956+ss).fit(V)
    raw=np.full(mask.shape,-1,int); raw[mask]=km.labels_
    sep0=separation_score(F,raw,mask)
    # Smooth categorical probability fields; S560 gets a wider candidate range.
    sigmas=[1.4,2.0,2.8,3.6,4.5] + ([5.5,6.5] if ss==560 else [])
    best=None
    for sig in sigmas:
        P=[]
        for k in range(10):
            P.append(gaussian_filter(((raw==k)&mask).astype(float),sig))
        sm=np.argmax(np.stack(P),axis=0)
        sm[~mask]=-1
        sep=separation_score(F,sm,mask)
        bc=boundary_complexity(sm,mask)
        comp=components_count(sm,mask)
        frac=np.array([np.mean(sm[mask]==k) for k in range(10)])
        minf=float(frac.min())
        retention=sep/(sep0+1e-9)
        # Require biology retention. S560 explicitly rewards smoother boundary.
        valid=(retention>=0.88 if ss==560 else retention>=0.90) and minf>=0.008
        obj=sep - (0.90 if ss==560 else 0.65)*bc - 0.006*max(comp-10,0)
        sweep_rows.append(dict(section=ss,sigma=sig,separation=sep,retention=retention,boundary_complexity=bc,
                               components=comp,min_region_fraction=minf,valid=valid,objective=obj))
        if valid and (best is None or obj>best[0]):
            best=(obj,sig,sm,sep,retention,bc,comp,minf)
    if best is None:
        # fallback to the smooth candidate with best objective
        q=[r for r in sweep_rows if r["section"]==ss]
        br=max(q,key=lambda z:z["objective"])
        sig=br["sigma"]
        P=[gaussian_filter(((raw==k)&mask).astype(float),sig) for k in range(10)]
        sm=np.argmax(np.stack(P),axis=0);sm[~mask]=-1
        best=(br["objective"],sig,sm,br["separation"],br["retention"],br["boundary_complexity"],br["components"],br["min_region_fraction"])
    _,sig,sm,sep,ret,bc,comp,minf=best
    if ss==560:
        valid560=[r for r in sweep_rows if r["section"]==ss and r["valid"]]
        if valid560:
            br=max(valid560,key=lambda z:z["sigma"])
            sig=br["sigma"]
            P=[gaussian_filter(((raw==k)&mask).astype(float),sig) for k in range(10)]
            sm=np.argmax(np.stack(P),axis=0); sm[~mask]=-1
            sep=separation_score(F,sm,mask);ret=sep/(sep0+1e-9)
    Pstack=np.stack([gaussian_filter(((raw==k)&mask).astype(float),sig) for k in range(10)])
    sm=enforce_contiguous(sm,mask,Pstack)
    sep=separation_score(F,sm,mask)
    ret=sep/(sep0+1e-9)
    bc=boundary_complexity(sm,mask)
    comp=components_count(sm,mask)
    frac=np.array([np.mean(sm[mask]==k) for k in range(10)])
    minf=float(frac.min())
    # Section-local numbering: left-to-right then top-to-bottom in DISPLAY coordinates. No homology across sections.
    cent=[]
    for k in range(10):
        yy,xx=np.where((sm==k)&mask)
        cent.append((k,float(xc[xx].mean()) if len(xx) else 1e9,float(yc[yy].mean()) if len(yy) else 1e9))
    order=sorted(cent,key=lambda z:(z[1],z[2]));mp={old:i+1 for i,(old,_,__) in enumerate(order)}
    local=np.zeros_like(sm)
    for old,new in mp.items(): local[sm==old]=new
    local[~mask]=0
    region_info[ss]=dict(local=local,mask=mask,xe=xe,ye=ye,xc=xc,yc=yc,markers=markers,sigma=sig)
    # Map every cell to its section-local region on the final smoothed grid.
    mm=section==ss
    ix=np.clip(np.digitize(xd[mm],xe)-1,0,len(xc)-1)
    iy=np.clip(np.digitize(yd[mm],ye)-1,0,len(yc)-1)
    cell_region10[mm]=local[iy,ix]
    region_manifest.append(dict(section=ss,required_markers=";".join(required),auto_specific_markers=";".join(auto),
                                all_markers=";".join(markers),smooth_sigma=sig,separation=sep,retention=ret,
                                boundary_complexity=bc,components=comp,min_region_fraction=minf))

pd.DataFrame(sweep_rows).to_csv(DATA/"stage956_marker_guided_region10_sweep.csv",index=False)
pd.DataFrame(region_manifest).to_csv(DATA/"stage956_marker_guided_region10_manifest.csv",index=False)
pd.DataFrame(dict(section=section,x_display=xd,y_display=yd,region10=cell_region10)).to_csv(DATA/"stage956_cell_region10_labels.csv.gz",index=False,compression="gzip")

# Downstream region × gene profile (final analysis, not a fitting intermediate).
rows=[]
for ss in secs:
    for rr in range(1,11):
        q=(section==ss)&(cell_region10==rr)
        if not q.any(): continue
        mu=np.log1p(X[q]).mean(0)
        for gi,g in enumerate(genes): rows.append(dict(section=ss,region=rr,gene=g,mean_log1p=float(mu[gi]),n_cells=int(q.sum())))
prof=pd.DataFrame(rows); prof.to_csv(DATA/"stage956_region10_gene_profiles.csv",index=False)
mat=prof.pivot_table(index=["section","region"],columns="gene",values="mean_log1p").reindex(columns=genes)
Z=(mat-mat.mean(0))/(mat.std(0)+1e-9)
fig,ax=plt.subplots(figsize=(13.8,9.2))
im=ax.imshow(Z.values,aspect="auto",cmap="coolwarm",vmin=-2.5,vmax=2.5)
ax.set_xticks(range(len(genes)));ax.set_xticklabels(genes,rotation=60,ha="right",fontsize=7)
ax.set_yticks(range(len(mat.index)));ax.set_yticklabels([f"S{s} R{r}" for s,r in mat.index],fontsize=7)
for y in [9.5,19.5]: ax.axhline(y,color="white",lw=1.5)
ax.set_title("Marker-guided region10 × 27-gene profiles",fontsize=12)
cb=fig.colorbar(im,ax=ax,fraction=.025,pad=.02);cb.set_label("gene-wise z score",fontsize=8);cb.ax.tick_params(labelsize=7)
fig.tight_layout();fig.savefig(OUT/"REGION10_GENE_PROFILE.png",dpi=180,bbox_inches="tight");plt.close(fig)

# Region reference figure
fig,axs=plt.subplots(1,3,figsize=(13.2,4.5))
tab=plt.get_cmap("tab10")
for ax,ss in zip(axs,secs):
    info=region_info[ss]; G=info["local"]; mask=info["mask"]; xc=info["xc"];yc=info["yc"]
    for L in range(1,11):
        z=((G==L)&mask).astype(float)
        if z.max()>0:
            ax.contourf(xc,yc,z,levels=[.5,1.5],colors=[tab(L-1)],alpha=.38)
            ax.contour(xc,yc,z,levels=[.5],colors=["#343434"],linewidths=.8,alpha=.85)
            yy,xx=np.where((G==L)&mask)
            if len(xx):
                ax.text(xc[xx].mean(),yc[yy].mean(),str(L),ha="center",va="center",fontsize=8,
                        bbox=dict(boxstyle="circle,pad=.18",fc="white",ec="#555",lw=.5,alpha=.9))
    ax.set_title(f"S{ss} · local regions 1–10",fontsize=10.5)
    ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
    for sp in ax.spines.values():sp.set_visible(False)
fig.suptitle("Marker-guided section-specific 10-subregion reference",fontsize=13.5,y=.99)
fig.text(.5,.015,"Required anchors: Hcrtr1 · Bcl11b · Slc5a7 · Lmx1a · Piezo2 · Slc6a2 · Ghr  |  + section-specific spatial markers",
         ha="center",fontsize=8.5,color="#555")
fig.tight_layout(rect=[0,.04,1,.94])
regpng=OUT/"FINAL_MARKER_GUIDED_REGION10.png";fig.savefig(regpng,dpi=190,bbox_inches="tight");plt.close(fig)

# colorbar policy
policy=[]
manifest=[]
for gi,g in enumerate(genes):
    v=X[:,gi];pos=v[v>0];posfrac=float((v>0).mean())
    if len(pos)==0:
        q=.995;vmax=1.0;gamma=.7;kind="empty"
    elif g=="Snap25":
        q=.85;vmax=float(np.quantile(pos,q));gamma=.50;kind="broad-special"
    elif posfrac>=.80:
        q=.95;vmax=float(np.quantile(pos,q));gamma=.58;kind="broad"
    elif posfrac>=.50:
        q=.98;vmax=float(np.quantile(pos,q));gamma=.64;kind="intermediate"
    else:
        q=.995;vmax=float(np.quantile(pos,q));gamma=.72;kind="sparse"
    vmax=max(vmax,1e-6)
    policy.append(dict(gene=g,positive_fraction=posfrac,vmax_quantile=q,vmax=vmax,gamma=gamma,display_class=kind))
    norm=PowerNorm(gamma=gamma,vmin=0,vmax=vmax)
    for ss in secs:
        m=section==ss;vv=X[m,gi];xx=xd[m];yy=yd[m]
        fig,ax=plt.subplots(figsize=(5.0,4.7))
        ax.scatter(xx,yy,s=.45,c="#e8e8e8",alpha=.58,linewidths=0,rasterized=True)
        qp=vv>0
        sc=ax.scatter(xx[qp],yy[qp],s=1.55,c=np.clip(vv[qp],0,vmax),cmap="Reds",norm=norm,
                      alpha=.90,linewidths=0,rasterized=True)
        info=region_info[ss];G=info["local"];mask=info["mask"];xc=info["xc"];yc=info["yc"]
        for L in range(1,11):
            z=((G==L)&mask).astype(float)
            if z.max()>0: ax.contour(xc,yc,z,levels=[.5],colors=["#777"],linewidths=.42,alpha=.32)
        ax.set_title(f"{g} · S{ss}",fontsize=11)
        ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
        for sp in ax.spines.values():sp.set_visible(False)
        cb=fig.colorbar(sc,ax=ax,fraction=.047,pad=.025)
        # explicit range; tick midpoint only when useful
        mid=vmax/2
        cb.set_ticks([0,mid,vmax]);cb.set_ticklabels([f"0",f"{mid:.0f}" if vmax>=10 else f"{mid:.1f}",f"{vmax:.0f}" if vmax>=10 else f"{vmax:.1f}"])
        cb.ax.tick_params(labelsize=7);cb.set_label("spot count",fontsize=8)
        fig.text(.5,.015,f"display range: 0–{vmax:.1f} spot count · x flipped, y flipped",ha="center",fontsize=7.5,color="#666")
        fig.tight_layout(rect=[0,.035,1,1])
        fp=OUT/f"{g}_S{ss}_expression.png";fig.savefig(fp,dpi=185,bbox_inches="tight");plt.close(fig)
        im=Image.open(fp).convert("RGB");im.thumbnail((820,780));im.save(OUT/f"{g}_S{ss}_expression_preview.jpg",quality=86)
        try: fp.unlink()
        except Exception: pass
        manifest.append(dict(gene=g,section=ss,file=f"{g}_S{ss}_expression_preview.jpg",preview=f"{g}_S{ss}_expression_preview.jpg",
                             positive_fraction=posfrac,vmax=vmax,vmax_quantile=q,gamma=gamma,display_class=kind,
                             x_flipped=True,y_flipped=True,region_reference="stage956_marker_guided_region10"))
pd.DataFrame(policy).to_csv(DATA/"stage956_gene_colorbar_policy.csv",index=False)
pd.DataFrame(manifest).to_csv(DATA/"stage956_gene_section_manifest.csv",index=False)

auth=dict(stage=956,status="CURRENT_FINAL_DISPLAY_ATLAS",expression_authority="Stage943 Route A",
          orientation="x flipped, y flipped for display in every gene map and region reference",
          colorbar="explicit on every gene×section map; shared range across three sections for each gene",
          snap25="special broad-expression display: 85th percentile positive-count vmax + gamma 0.50",
          region_reference=dict(n_regions=10,section_specific=True,cross_section_homology_assumed=False,
              required_marker_anchors=required,auto_specific_markers_per_section={str(r["section"]):r["auto_specific_markers"].split(";") for r in region_manifest},
              note="marker-guided smoothed fields + geometry support; S560 receives stronger smoothing search"),
          main_page_rule="final atlas only; intermediate correction/validation plots omitted from main page")
(DATA/"stage956_final_atlas_authority.json").write_text(json.dumps(auth,indent=2),encoding="utf-8")
print(json.dumps(auth,indent=2))