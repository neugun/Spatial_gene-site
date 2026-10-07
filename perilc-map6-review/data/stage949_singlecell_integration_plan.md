# Stage949 — single-cell ↔ final spatial atlas integration plan

## Current authority
- Expression: Stage943 Route A, 71,950 segmentation-defined cells × 27 genes.
- Final display atlas: Stage947/948.
- Main page remains minimal: LC/periLC locator + section-local 10-region reference + 81 gene×section maps.
- Do not reintroduce correction/Route B/intermediate panels into the main atlas page.

## What from the prior single-cell work should be connected
1. Frozen Fine26 molecular identities; do not re-cluster merely to create a new version.
2. UMAP / molecular evidence for the 27-gene single-cell state space.
3. Cross-section transfer of Fine26 signatures.
4. Local single-cell multiplex views that connect the global atlas back to individual cells.
5. Section-wise spatial domains, now using 10 local regions as the preferred display granularity.
6. True segmentation-body / maskbody views for selected single-cell examples and 3D atlas validation.
7. Allen / MapMyCells correspondence as secondary external annotation, not the organizing backbone.
8. High-dimensional latent information / atlas-lifting analyses remain separate downstream work; cell-type labels are not assumed to explain function.

## Remaining gaps to close
### A. Section-wise 10-region authority
Current Stage947 uses the historical global K10_w1 cell labels and then relabels/smooths them separately within S500/S530/S560.
This satisfies the display requirement but is not yet a true per-section independent 10-region fit.
Next:
- fit 10 spatially constrained regions independently within each section using current Stage943 expression + section geometry;
- compare to historical K10_w1 only as a reference;
- smooth display boundaries sectionwise;
- do not force region correspondence across sections;
- report contiguity, minimum region size, molecular separation, seed stability, and boundary sensitivity;
- promote only if it improves the current Stage947 display without making biologically implausible fragments.

### B. Re-run single-cell molecular figures from Stage943 current counts
The existing Stage704/718/725 scripts read CURRENT_ROUTEA, but historical output labels/titles still reference Stage934.
Rebuild under a new current bundle:
- 27-gene UMAP montage;
- Fine26 dotplot / molecular hierarchy;
- leave-one-section-out Fine26 transfer;
- one local multiplex field in each section (not only S530);
- section abundance / stability;
- selected-cell heatmaps.
Keep Fine26 frozen; no opportunistic re-clustering.

### C. True cell-body rendering
The final gene atlas currently uses point/scatter rendering.
For the single-cell branch:
- restore true segmentation-body / maskbody rendering for representative local fields and the 3D atlas;
- preserve physical coordinates and flipped-Y convention;
- use points only for overview/QC, not as the only single-cell anatomical rendering.

### D. Link 10-region spatial structure to single-cell molecular state
For each section independently:
- Fine26 × region10 enrichment;
- 27-gene region profiles;
- within-Fine26 spatial modulation;
- continuous gradients within regions;
- morphology / local-density context.
Do not assume region 1 in S500 equals region 1 in S530 or S560.

### E. External annotation
Update MapMyCells / Allen consensus using the current Stage943 matrix where applicable.
Treat it as correspondence, not definitive identity, because only 27 genes are measured.
Keep external atlas mapping downstream of the direct molecular/spatial evidence.

### F. Functional / atlas-lifting bridge
Separate from the final display page:
- test whether continuous 27-gene latent state predicts projection / neural / behavioral endpoints beyond Fine26 labels and beyond section/region;
- compare cell-type-only vs latent-state vs spatial+latent models;
- use held-out section / dataset gates;
- only retain spatial priors when real-vs-shuffle and held-out endpoint gains are positive.

## Deliverables
1. Stage949 section-wise independent region10 authority + smooth maps.
2. Stage950 current single-cell molecular bundle (UMAP, Fine26, transfer, local multiplex).
3. Stage951 true-maskbody single-cell / 3D validation.
4. Stage952 region10 × Fine26 / gene / morphology integration.
5. Stage953 optional Allen/MapMyCells current correspondence.
6. Stage954 latent-function / atlas-lifting bridge.
7. Main atlas page stays minimal; single-cell results go to a separate page/tab until explicitly promoted.
