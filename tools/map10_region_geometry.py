"""Topology-aware map10 display-region refinement, never changes cell assignments.
No post-hoc flipping of raster images. Geometry uses the already x-flipped/y-flipped
source-coordinate grid. Do not imply that smoothing validates anatomical identity.
"""
import numpy as np
from scipy.ndimage import (binary_fill_holes, binary_dilation, distance_transform_edt,
                           label, gaussian_filter)

def topology_counts(G, mask):
    from scipy.ndimage import label
    void=binary_fill_holes(mask)&~mask
    _,whole_holes=label(void)
    n_components={k:label(G==k)[1] for k in range(1,11)}
    n_holes={k:label(binary_fill_holes(G==k)&(G!=k))[1] for k in range(1,11)}
    return dict(white_holes=int(whole_holes),white_pixels=int(void.sum()),
                region_components=n_components,region_holes=n_holes)

def refine_region10(original, original_mask, *, sigma=2.0, max_change_fraction=.045):
    original=np.asarray(original,dtype=int)
    orig_mask=np.asarray(original_mask,bool)
    cc,n=label(orig_mask)
    counts=np.bincount(cc.ravel())
    main=cc==(int(np.argmax(counts[1:]))+1)
    # Fill ALL fully enclosed voids; do not fill channels that open to background.
    mask=binary_fill_holes(main)
    G=original.copy()
    G[~mask]=0
    # At most 3 grid pixels of smoothing at category-level; retain marker-like geometry.
    missing=mask & (G==0)
    _,idx=distance_transform_edt(G==0,return_indices=True)
    nearest=G[tuple(idx)]
    # Single enclosed void gets an adjacent regional identity; no mixed / dotted fill.
    islands,n=label(missing)
    for i in range(1,n+1):
        q=islands==i
        b=binary_dilation(q)&(G>0)
        neigh=G[b]
        if neigh.size:
            val=int(np.bincount(neigh,minlength=11)[1:].argmax()+1)
            G[q]=val
        else:G[q]=nearest[q]
    # Spatially smooth category probabilities; restrict categorical smoothing
    # with agreement to the frozen section-local region reference.
    best=None
    for sg in [0.8,1.1,1.4,1.7,2.0,2.4,2.8,3.0]:
        P=np.stack([gaussian_filter((G==k).astype(float),sg) for k in range(1,11)])
        L=P.argmax(0)+1
        L[~mask]=0
        agreement=float(np.mean(L[main]==original[main]))
        ious=[float(np.logical_and(L==k, original==k).sum()/max(1,np.logical_or((L==k)&main,(original==k)&main).sum()))
              for k in range(1,11)]
        if agreement>=.97 and min(ious)>=.88:best=(L,sg,agreement,min(ious))
    if best is None:raise ValueError("Smoothing fails shape preservation gate")
    L,sg,agreement,minio=best
    changed=[]
    # Region fragments are small disjoint pieces within a connected tissue mask;
    # reassign ONLY fragments, leaving largest component of each region intact.
    for it in range(10):
        alterations=0
        for k in range(1,11):
            c,n=label(L==k)
            if n<=1:continue
            sizes=np.bincount(c.ravel())[1:]
            major=int(np.argmax(sizes))+1
            for i in range(1,n+1):
                if i==major:continue
                q=c==i
                b=binary_dilation(q)&(L>0)&~q
                neighbors=L[b]
                neighbors=neighbors[neighbors!=k]
                if neighbors.size==0:continue
                replacement=int(np.bincount(neighbors,minlength=11).argmax())
                if replacement<1:continue
                changed.append(dict(from_region=k,to_region=replacement,pixels=int(q.sum())))
                L[q]=replacement;alterations+=int(q.sum())
        if alterations==0:break
    # Fix enclosed holes in individual regions (including internal tiny islands).
    for it in range(10):
        fixes=0
        for k in range(1,11):
            hole=binary_fill_holes(L==k)&(L!=k)
            if not hole.any():continue
            changed.append(dict(from_region="nested",to_region=k,pixels=int(hole.sum())))
            L[hole]=k;fixes+=int(hole.sum())
        if fixes==0:break
    after=topology_counts(L,mask)
    total_changes=int(((L!=original)&main).sum())
    fraction=total_changes/max(1,int(main.sum()))
    if fraction>max_change_fraction:raise ValueError(f"Large label change fraction {fraction:.3f}")
    if after["white_holes"] or any(v!=1 for v in after["region_components"].values()) or any(after["region_holes"].values()):
        raise ValueError(f"Failed topological gate {after}")
    return L,mask,dict(before=topology_counts(original,orig_mask),
                       after=after,orig_main_pixels=int(main.sum()),
                       inner_voids_filled=int((mask&~main).sum()),
                       detached_support_pixels=int((orig_mask&~main).sum()),
                       changed_existing_pixels=total_changes,
                       changed_existing_fraction=fraction,
                       smoothing_sigma=float(sg),agreement_original=agreement,
                       minimum_region_iou_original=minio,component_edits=changed)
