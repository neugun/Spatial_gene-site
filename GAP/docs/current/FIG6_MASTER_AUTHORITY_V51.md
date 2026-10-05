# Figure 6 Master Authority v5.1 — local finalization, 2026-10-02

## Figure question

How far can one routed GAP framework generalize without replacing it with a different final model for every dataset?

## Public model identity

The public model identity throughout this figure is **GAP**. TranscriptFormer, BrainBeacon, AtlasLift, V2.2/function heads, residual adapters, transport gates, and Compact/Large backbones are internal modules, routing branches, or ablations. They are not competing final model identities.

Across five pre-specified capability blocks, GAP improves **0.874** of matched metric views after equal weighting the blocks. This is an audited capability summary, not a pooled leaderboard across unrelated endpoints. Block-level definitions remain in Source Data; the structured Operator result itself remains in Figure 5.

## Main figure contract

**A. One GAP framework.** Ortholog-aware molecular input enters a shared GAP encoder; context routing selects reference/function/transport modules before held-unit functional prediction. Audited matched-capability win fraction = **0.874**. GAP-Compact uses 1,007,639 versus 5,957,243 trainable parameters for GAP-Large (**83.1% fewer**).

**B. Sample efficiency.** Bugeon held-mouse R2 improves from 0.0753 to 0.1289 at 20% target labels (delta +0.0536) and from 0.0767 to 0.1188 at 50% labels (delta +0.0421).

**C. Same-species reuse.** Zhao→Bugeon temporal variability at 10% labels gives delta R2 +0.1546, delta Spearman +0.0781 and delta AUROC +0.0632, with 4/4 held mice improving in R2. Gains remain positive at 20% and 50%.

**D. Direct mouse→fish champion.** Across five matched motor/OMR axes, GAP-Transport improves mean R2 by +0.0567, mean Spearman by +0.0687 and mean AUROC by +0.0503 relative to the frozen TF-only branch. Held-fish delta R2 is positive in 15/15 task×fish comparisons. Routing falls back to the simpler branch when residualization is harmful.

**E. Six-direction coverage.** Under one frozen routing contract, held-target Spearman is positive for all six exact mouse/fish/worm directions: mouse→fish 0.1402, fish→mouse 0.1551, mouse→worm 0.0946, worm→mouse 0.3698, fish→worm 0.0857 and worm→fish 0.1816. This is routing coverage, not a claim that every internal branch improves every edge.

**F. Direct mouse→human function.** Under the same source-PCA8/Ridge-alpha=1 contract, TF cell-context exceeds TF gene-pooled by delta rho **+0.1038** (15/18 endpoints improve; paired p=0.002373). GAP FullContext exceeds matched pooled GAP by **+0.1103** (14/18 improve; p=0.006016). FullContext human mean rho = **0.0650**, with 15/18 endpoints positive. Human target-donor bootstrap remains positive, 95% CI **[0.0250, 0.0973]**, P(mean<=0)=0.0018. A strict mouse-source donor bootstrap that refits scaling, PCA and Ridge has 95% CI **[-0.0191, 0.1098]**, P(mean<=0)=0.099. All 76 strict source-donor leave-one-out estimates remain positive (minimum rho 0.0227), as do all 56 target-donor exclusions (minimum rho 0.0586). Cell-level TF is the dominant transferable branch; FullContext does not significantly exceed TF-cell alone.

**G. Compactness.** GAP-Compact improves 13/21 matched metric cells across Condylis, Xu and Zhao while using 83.1% fewer trainable parameters than GAP-Large. Report as an efficiency tradeoff, not universal superiority.

**H. Broad held-group utility.** XC64 evaluation is complete for 60/60 held-group runs and H05S5 mouse evaluation for 42/42 runs. Compact GAP retains positive rank/ROC utility across WARP, Xu, Condylis, Shainer and Zhao even where absolute R2 is small or negative.

## Extended Data contract

**ED1 — audited GAP winner + generalization.** Show the single audited capability-level summary (0.874), Compact efficiency (13/21 matched metric cells; 83.1% fewer parameters), and H05/XC held-group Spearman/AUROC across independent datasets. The full classical/native model zoo is retained in Source Data and is not the visual headline.

**ED2 — full directionality.** Show all predeclared positive and failed/reverse cross-species directions. Main Fig. 6E shows the frozen six-direction routed coverage; ED2 prevents that from being interpreted as universal transfer.

**ED3 — human mechanism and biological-unit robustness.** Show all 18 human ephys endpoints, same-contract cell-context versus pooled tests, endpoint/target-donor/strict-source-donor uncertainty, leave-one-donor ranges, and TF/BB branch decomposition. The strict source-donor interval crosses zero and must remain visible.

## Boundaries

1. GAP is the selected reusable framework across audited capabilities; do not claim one scalar architecture wins every endpoint.
2. Do not mix unrelated endpoint definitions into a raw leaderboard.
3. Structured Operator results remain in Figure 5; Figure 6 only references the audited capability summary.
4. Mouse→human means intrinsic electrophysiology, not behavior, disease or clinical prediction.
5. Do not claim BrainBeacon is required for human transfer; cell-level TranscriptFormer is the dominant transferable component.
6. Do not claim FullContext significantly exceeds TF-cell alone in mouse→human Patch-seq.
7. Human target-donor resampling is positive; strict source-donor bootstrap crosses zero, so mouse source-cohort composition is the main remaining sampling uncertainty.
8. Shainer→WARP common-reference zero-shot and Allen strict zero-shot remain boundary/negative evidence.
9. Six-direction coverage reflects source-validated routing/fallback, not universal improvement by one internal branch.

## Reproducibility

Main renderer: scripts/895_render_fig6_gap_generalization_v51_local.py
ED1 renderer: scripts/896_render_fig6_v51_ed1_winner_local.py
ED2 renderer: scripts/885_render_fig6_crossspecies_ed2_local.py, copied into the v5.1 canonical output
ED3 renderer: scripts/897_render_fig6_v51_human_ed3_local.py
Human authority: results/mouse_human_gap_20261002/HUMAN_TRANSFER_AUTHORITY_V3.json
Current figure directory: analysis_workspace/results/manuscript\Fig6_GAP_Generalization_v51_20261002
Network/Titan is not required to reproduce the current v5.1 package from mirrored local assets.
