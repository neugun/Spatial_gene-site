# Figure requirements v3 — current authority

**Status:** this document supersedes the older `Fig*_FULL_v2.png` layout as the figure-design authority.
The old PNGs remain in the repository for provenance but should not be treated as the final current figures.

The design rule is **biological question first, matched evidence second, model name third**. Every quantitative panel must state its target, held biological unit, metric and matched null/baseline. The current browser figures are rendered by [`scripts/render_requirement_aligned_figures.py`](../scripts/render_requirement_aligned_figures.py); the renderer follows this file rather than defining new scientific claims.

## Figure 1 — What part of function is explained by identity, within-identity molecular state, and position?

Required panels:

A. **Three-layer law:** discrete identity/type → continuous molecular state within identity → position/region context.

B. **Mouse representative hierarchy:** Xu/PVH and Bugeon/VISp, type vs type+genes, showing R², rank and AUROC rather than R² alone.
- Xu Ghrelin: R² 0.1259→0.2492; rho 0.3774→0.5299; AUROC 0.7465→0.8470.
- Bugeon running/state: R² 0.1668→0.1783; rho 0.4806→0.5170; AUROC 0.8027→0.8337.

C. **Zebrafish WARP uses actual cell types, not region labels.**
Across 8 axes: Δrho positive 8/8; ΔAUROC 8/8; ΔR² 5/8; null-significant R²/rho/AUROC = 4/8, 5/8, 6/8.

D. **C. elegans identity-dominant boundary.**
Pumping identity R² 0.1067; velocity identity R² 0.5328; within-identity molecular residual adds approximately zero.

E. **Cross-system frequency summary:** how often within-type continuous molecular state adds held-unit information.

F. **Position/context test, not anatomy-as-truth.**
Mouse same-region transfer is not automatically stronger (mean rho -0.0566) than different-region same-superclass pairs (+0.01217). WARP region-resolved transfer is heterogeneous and shows no global homologous-visual-region advantage.

G. **Cross-region / homologous-region comparison** remains a required panel. It must use the same source/target semantic axis, target cells, held unit and metric. No visual-vs-nonvisual conclusion is promoted from unmatched endpoints.

## Figure 2 — Why AtlasLift is needed beyond coarse type/markers

The figure must not imply that cell type is unimportant. The claim is: **type is a useful first coordinate, but sparse type/marker information can leave functionally relevant molecular state unresolved.**

Required panels:

A. Type / measured sparse markers → continuous AtlasLift state → optional anatomy/space expert.

B. **Within-type recovery:** Bugeon within-type ordering AUROC 0.654→0.705 (+0.051).

C. **Correspondence control:** real AtlasLift functional-prior gain +0.0146 R² versus shuffled correspondence -0.0271.

D. **Unified visual-system contract:** Bugeon and Allen VC2P/VBO under the same held-unit/multimetric ladder; anatomy-conditioned increments shown separately.

E. **PVH MC4R must be visible in the main figure.**
Atlas-inferred Mc4r raises Xu Ghrelin R² **0.2589→0.2711** (Δ +0.01223; random-gene p=0.00995), extreme AUROC +0.00349 (p=0.0149), with negative predicted direction. Exact-PVH restriction preserves the effect (ΔR² +0.00725).

F. Kirrel3 and Bugeon named-gene controls belong beside MC4R to separate set-level atlas gain from gene-specific evidence.

G. **Type vs type+AtlasLift vs type+AtlasLift+position/anatomy** must be shown wherever the required inputs exist. Space is optional and must beat a coordinate shuffle.

## Figure 3 — Which compact measurements should the experiment buy?

Required panels:

A. AtlasLift candidate state → nested single-gene test → compact panel → prospective assay.

B. Gene-specific null: Gabra4, Prkacb, Mc4r, Kirrel3; Cacna2d3 specificity control.

C. Exact-anatomy sensitivity for nominated genes.

D. **Fixed-budget selector comparison**, not a universal winner:
M1 Patch-seq K25/K50/K100 with FSO-joint, PERSIST, random and task-aware alternatives.

E. Panel Pareto axes should include functional information, operator/shape preservation and assay budget.

F. Receptor/effector-rich WARP panel remains a phenotype-specific example, not a generic marker claim.

G. Prospective PVH validation: Ghrelin vs saline, Fos/cFos + Mc4r/Kirrel3 + matched random genes + identity/anatomy controls.

## Figure 4 — BehLift plus gene-predictable functional manifold

Required panels:

A. Leave-one-behavior/state-out completion.

B. Xu/PVH conditional complementarity:
other 10 states R² ≈0.422; genes ≈0.276; joint ≈0.486; conditional molecular ΔR² +0.0639, p=0.00498.

C. WARP-D2 AtlasLift×BehLift:
OMR-backward conditional ΔR² +0.024776 and Δ extreme-AUROC +0.013643; p=0.00990.

D. **Held-unit neural-manifold reconstruction**, not decorative UMAP:
C. elegans latent R² 0.267 versus direct 0.274 and shuffled -0.024.

E. **Gene/program → manifold geometry** must be explicit.
WARP between-group temporal-subspace angle 0.743 vs within-group split-half floor 0.629; cross-animal geometry Spearman 0.643.
Neuromodulator/neuropeptide family-distance predicts manifold distance: pooled rho 0.219; all three fish 0.222–0.261; within-fish permutation p<0.005.
Transcription-factor, adhesion/guidance and synaptic-effector families do not show the same reproducible positive relation.

F. Boundary regimes: WARP/Shainer global latent R² near -0.006 but modestly above matched molecular shuffles; Neuroplex approximately null. A manifold is promoted only when geometry/reconstruction survives biological holdout and nulls.

## Figure 5 — Molecularly parameterized state-dependent response law

Required panels:

A. Theta stable response geometry + Gamma state/context susceptibility.

B. **Main operator comparison:** independent Ridge heads R² 0.14454 → shared low-rank operator 0.16356; parameters ~438→220.5 (~49.7% reduction); 3/4 held mice improve.

C. Supporting same-question stress tests:
shape 0.19822→0.21721; 16-target 0.03774→0.04280; strict 24-D 0.00927→0.01852. The 24-D result is not the overall operator headline.

D. Conserved vs animal-specific axes: conserved axes mean gene AUROC ≈0.657 vs ≈0.456 for the other axes.

E. Gamma susceptibility prediction:
type lookup ≈0.0732; type+BrainBeacon4+Atlas8+receptor6 ≈0.11546, 4/4 held animals positive.

F. Xu/PVH cross-state flexibility:
amplitude variability ΔR² +0.0580; state range +0.0990; 3/3 mice; matched-null p=0.00498.

G. Reverse shared geometry:
NuCLR bACC 0.4429 → +shape 0.4824 → +shape+gain 0.5173, with matched-dimensional and wrong-neuron controls.

## Figure 6 — Foundation / evolutionary transfer, all regimes in one figure

This figure must contain the full transfer story while keeping regimes separate.

A. Contract schematic: strict zero-shot / matched representation / mechanism-ablation / few-shot / mammalian-primate extension.

B. **Strict zero-shot six-direction TF64 conserved-axis map**
mouse→fish 0.09766; fish→mouse 0.14599; mouse→worm 0.09049; worm→mouse 0.29161; fish→worm 0.11236; worm→fish 0.16082 (Spearman).

C. Semantic-axis permutation results and directional asymmetry; do not average these into one global transfer score.

D. **Matched Zhao→WARP model comparison** must use the same target cells, semantic axis, held-fish split and metric for TF64 / BrainBeacon / LangPatch / GAP. Rows that are not exactly matched stay out of the winner comparison.

E. Mechanism/ablation:
HCR/reference, ortholog mapping, gene-program gate, contextual TF geometry, region boundary.
Do not claim ion-channel transfer is globally strongest; current strongest interpretable program evidence is axis-specific, with neuromodulator/neuropeptide emphasized where supported.

F. Few-shot:
strict WARP→Shainer curves and Bugeon V2.2 sample-efficiency curves, visibly separated from zero-shot.

G. Foundation/hybrid boundary:
Z4 versus BB/TF additions under the same Bugeon contract; do not call BB/TF generic additive wins.

H. Mammalian/primate extension:
mouse→macaque→human evolutionary gene geometry and, where available, thalamus/visual-cortex functional/ephys/program tests. TF64 vs BrainBeacon must be framed as complementary capabilities unless a matched functional endpoint says otherwise.

I. **Pharaoh ant control is required**, but no frozen local authority was located in the current repository audit. It remains MISSING/PENDING and must not be fabricated.

## Figure 7 — Frozen GAP → perturbation nomination → real causal test

Figure 7 must visually depend on Figures 1–6.

A. Freeze the observational GAP/operator from training animals.

B. Nominate a molecule/controller and a specific Theta/Gamma coordinate from the frozen model.

C. Prospective intervention prediction: sign, magnitude and state dependence are specified before experiment.

D. Perturb in new biological units and measure the predicted coordinate plus orthogonal coordinates.

E. Causal support requires coordinate-selective movement versus matched perturbation/anatomy/state controls.

F. Synthetic six-animal counterfactual remains a **sanity check only**:
Gamma-only true 1.547, predicted 1.545, orthogonal ΔTheta 0.014; counterfactual R² 0.997; cosine 0.99999. Theta-only R² 0.998.

G. The panel must distinguish **nomination / prospective prediction / causal proof**. Observational prediction alone is never labelled causal.

## Global completion gate

Do not mark the figure set COMPLETE until:
1. per-dataset authority table is complete: native metric → matched held-unit definition → baseline → GAP → provenance;
2. ranking/extreme/retrieval feasibility matrix is PASS / N.A. / MISSING, with retrieval coverage explicitly audited;
3. unused FSO papers are mapped to Molecular→Function claim, region/cell type, molecular and functional variables, supporting panel, observational/causal evidence level and provenance.
