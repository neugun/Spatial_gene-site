from pathlib import Path
import sys
import json, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PUB=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review")
AS=PUB/"assets"/"stage967_maskbody_xy";AS.mkdir(parents=True,exist_ok=True);DATA=PUB/"data"
ROOT=Path(sys.argv[1])
MC=Path(sys.argv[2])
A=pd.read_csv(ROOT/"stage631_current_3d_atlas"/"CURRENT_3D_CELL_ATLAS.csv.gz",
              usecols=["section","roi_id","fine26","fine26_stable_id","fine26_name","broad_x","x_um","y_um","depth_um"])
windows=pd.read_csv(DATA/"stage950_local_singlecell_windows.csv")
cmap=plt.get_cmap("turbo")
broad_cols={"E":"#E69F00","I":"#56B4E9","NE":"#00BFAE","ChAT":"#F0E442","Other":"#8A8A8A"}

def lut_for(sec,rid,field):
    q=A[A.section==sec]
    ids=q.roi_id.astype(int).to_numpy()
    maxid=max(int(ids.max()),int(rid.max()))
    if field=="fine26":
        lut=np.full(maxid+1,-1,int);lut[ids]=q.fine26.astype(int).to_numpy();return lut[rid]
    mp={k:i for i,k in enumerate(["E","I","NE","ChAT","Other"])}
    lut=np.full(maxid+1,-1,int);lut[ids]=[mp.get(str(x),4) for x in q.broad_x];return lut[rid]


def cloud_to_display(sec,rid,cloud):
    """Match maskcloud coordinates to atlas x-flipped/y-flipped convention by ROI IDs.
    Source maskcloud X has negative correlation with raw atlas X and is already
    x-flipped; despite its name y_flipped_um correlates positively with atlas raw Y.
    Fit only additive origins from robust matched-ROI medians, never distort bodies.
    """
    q=A[A.section==sec]
    ids=q.roi_id.astype(int).to_numpy()
    hi=max(int(ids.max()),int(rid.max()))
    rx=np.full(hi+1,np.nan); ry=np.full(hi+1,np.nan)
    rx[ids]=q.x_um.to_numpy(float);ry[ids]=q.y_um.to_numpy(float)
    src_x=cloud["x_um"].astype(float)/2.
    src_y=cloud["y_flipped_um"].astype(float)/2.
    ok=np.isfinite(rx[rid])&np.isfinite(ry[rid])
    assert ok.sum()>1000
    xs=q.x_um.to_numpy(float);ys=q.y_um.to_numpy(float)
    desired_x=xs.min()+xs.max()-rx[rid[ok]]
    desired_y=ys.min()+ys.max()-ry[rid[ok]]
    shiftx=float(np.median(desired_x-src_x[ok]))
    shifty=float(np.median(desired_y+src_y[ok]))
    transformed_x=src_x+shiftx
    transformed_y=-src_y+shifty
    corrx=float(np.corrcoef(transformed_x[ok],desired_x)[0,1])
    corry=float(np.corrcoef(transformed_y[ok],desired_y)[0,1])
    assert corrx>.995 and corry>.995,(sec,corrx,corry)
    coordinate_audit[str(sec)]=dict(x_flipped=True,y_flipped=True,
         roi_points_matched=int(ok.sum()), x_agreement_correlation=corrx,
         y_agreement_correlation=corry,
         median_abs_x_um=float(np.median(np.abs(transformed_x[ok]-desired_x))),
         median_abs_y_um=float(np.median(np.abs(transformed_y[ok]-desired_y))))
    return transformed_x,transformed_y

coordinate_audit={}

# Whole three-section true mask-body Fine26 atlas
fig=plt.figure(figsize=(18,6),facecolor="black")
for j,sec in enumerate([500,530,560]):
    C=np.load(MC/f"maskcloud_s{sec}_xy5_z7_cap45.npz");rid=C["roi_id"]
    xx,yy=cloud_to_display(sec,rid,C); zz=C["z_um"]/2.
    lab=lut_for(sec,rid,"fine26")
    ax=fig.add_subplot(1,3,j+1,projection="3d");ax.set_facecolor("black")
    for c in range(26):
        m=lab==c
        if m.any():ax.scatter(xx[m],zz[m],yy[m],s=.16,color=cmap(c/25),alpha=.58,linewidths=0,depthshade=False,rasterized=True)
    ax.view_init(elev=16,azim=-66);ax.set_title(f"S{sec}",color="white",fontsize=10)
    ax.set_xlabel("X display (flipped) µm",color="white");ax.set_ylabel("Z µm",color="white");ax.set_zlabel("Y display (flipped) µm",color="white")
    ax.tick_params(colors="white",labelsize=5);ax.grid(False)
    for pane in (ax.xaxis.pane,ax.yaxis.pane,ax.zaxis.pane):
        pane.set_facecolor((0,0,0,0));pane.set_edgecolor((.6,.6,.6,.25))
fig.suptitle("Stage943 true segmentation bodies | x flipped, y flipped",color="white",fontsize=14)
fig.tight_layout(rect=[0,0,1,.94]);fig.savefig(AS/"FINE26_TRUE_MASKBODY_3D.png",dpi=230,bbox_inches="tight",facecolor="black");plt.close(fig)

# Broad classes on true mask bodies
fig=plt.figure(figsize=(18,6),facecolor="black")
keys=["E","I","NE","ChAT","Other"]
for j,sec in enumerate([500,530,560]):
    C=np.load(MC/f"maskcloud_s{sec}_xy5_z7_cap45.npz");rid=C["roi_id"]
    xx,yy=cloud_to_display(sec,rid,C); zz=C["z_um"]/2.
    code=lut_for(sec,rid,"broad")
    ax=fig.add_subplot(1,3,j+1,projection="3d");ax.set_facecolor("black")
    for ci,k in enumerate(keys):
        m=code==ci
        if m.any():ax.scatter(xx[m],zz[m],yy[m],s=.17,color=broad_cols[k],alpha=.62,linewidths=0,depthshade=False,rasterized=True)
    ax.view_init(elev=16,azim=-66);ax.set_title(f"S{sec}",color="white",fontsize=10);ax.tick_params(colors="white",labelsize=5);ax.grid(False)
    for pane in (ax.xaxis.pane,ax.yaxis.pane,ax.zaxis.pane):
        pane.set_facecolor((0,0,0,0));pane.set_edgecolor((.6,.6,.6,.25))
fig.legend(handles=[Line2D([0],[0],marker="o",ls="",color=broad_cols[k],label=k,markersize=6) for k in keys],
           loc="upper center",ncol=5,frameon=False,labelcolor="white",bbox_to_anchor=(.5,.95))
fig.suptitle("Broad classes on true segmentation bodies | x flipped, y flipped",color="white",fontsize=14)
fig.tight_layout(rect=[0,0,1,.91]);fig.savefig(AS/"BROAD_TRUE_MASKBODY_3D.png",dpi=230,bbox_inches="tight",facecolor="black");plt.close(fig)

# local maskbody fields corresponding to Stage950 windows
local=[]
for _,r in windows.iterrows():
    sec=int(r.section);x0=float(r.x0);y0=float(r.y0);win=float(r.window_um)
    C=np.load(MC/f"maskcloud_s{sec}_xy5_z7_cap45.npz");rid=C["roi_id"];xx,yy=cloud_to_display(sec,rid,C);zz=C["z_um"]/2.
    # atlas y_um is non-flipped; match ROI identities rather than coordinates for selection
    qA=A[(A.section==sec)&(A.x_um>=x0)&(A.x_um<x0+win)&(A.y_um>=y0)&(A.y_um<y0+win)]
    ids=set(qA.roi_id.astype(int)); m=np.isin(rid,list(ids)); lab=lut_for(sec,rid,"fine26")
    fig=plt.figure(figsize=(7.2,6),facecolor="black");ax=fig.add_subplot(111,projection="3d");ax.set_facecolor("black")
    for c in sorted(qA.fine26.unique()):
        mm=m&(lab==c)
        if mm.any():ax.scatter(xx[mm],zz[mm],yy[mm],s=.8,color=cmap(c/25),alpha=.78,linewidths=0,depthshade=False,rasterized=True)
    ax.view_init(elev=18,azim=-64);ax.set_title(f"S{sec} local true segmentation bodies | x flipped, y flipped",color="white",fontsize=11);ax.tick_params(colors="white",labelsize=6);ax.grid(False)
    for pane in (ax.xaxis.pane,ax.yaxis.pane,ax.zaxis.pane):
        pane.set_facecolor((0,0,0,0));pane.set_edgecolor((.6,.6,.6,.25))
    fig.tight_layout();fig.savefig(AS/f"S{sec}_LOCAL_TRUE_MASKBODY.png",dpi=240,bbox_inches="tight",facecolor="black");plt.close(fig)
    local.append(dict(section=sec,n_cells=int(len(qA)),n_mask_points=int(m.sum()),n_fine26=int(qA.fine26.nunique())))
pd.DataFrame(local).to_csv(DATA/"stage967_local_maskbody_summary.csv",index=False)
auth=dict(stage=967,status="DISPLAY_FLIPS_X_AND_Y",expression_authority="Stage943 Route A",
          segmentation_source="Stage141 downseg-derived maskcloud; true segmentation-body samples, not cell-center points",
          coordinates="x flipped, y flipped relative to current 2D atlas; source 2x gel coordinates converted by /2 and aligned by matched ROI IDs",
          outputs=["FINE26_TRUE_MASKBODY_3D.png","BROAD_TRUE_MASKBODY_3D.png"]+[f"S{s}_LOCAL_TRUE_MASKBODY.png" for s in [500,530,560]])
(DATA/"stage967_maskbody_xy_authority.json").write_text(json.dumps(auth,indent=2),encoding="utf-8")
print(json.dumps(auth,indent=2))

(DATA/'stage967_maskbody_xy_coordinate_audit.json').write_text(json.dumps(coordinate_audit,indent=2),encoding='utf8')
print('STAGE967_COORDINATE_AUDIT',json.dumps(coordinate_audit),flush=True)
