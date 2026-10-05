# Current Figure 1–3 authority — 2026-10-05

This reviewer-audited authority supersedes the 2026-10-04 narrative for interpretation while retaining the frozen 136/136 runtime validation and all prior files for provenance.

## Global tool state
- Runtime: 136/136 PASS.
- Capability-aware execution: 14/14 audited datasets, 91/91 requested Fig1–3 capability checks, zero manifest errors in the frozen integration audit.
- Integration/stress benchmarks: McLachlan/Jerry Chen Cell Reports 11/11 PASS; Allen Visual Learning 2P+HCR 12/12 PASS; Condylis/Jerry Chen Science 8/8 PASS; MICRONS structural-only 13/13 PASS plus 8/8 selection-bias guard.
- EASI-PASS three-level intake safe-stops incomplete public inputs and executes the full applicable contract when same-cell molecular+functional multi-animal input is supplied.

## Figure 1 — reviewer-level inference
- Shainer 2025: space vs coarse type, held-fish mean ΔSpearman = +0.24206, 6/6 positive, exact one-sided sign-flip P=0.015625.
- Xu 2020 PVH: type+genes vs type, held-mouse mean ΔSpearman = +0.17292, 3/3 positive; n=3 limits the minimum unit-level P to 0.125.
- Bugeon 2022: genes have strong absolute readability (mean held-mouse rho≈0.396), but genes-vs-type mean ΔSpearman≈+0.02295; 3/4 positive; exact P=0.1875. This is not a strong incremental-gene claim.
- Zhao 2026: connectivity vs type gives relative Δrho=+0.11932 and exceeds a 1000× within-animal response-shuffle null (P=0.04496), but absolute connectivity rho=0.05437 does not exceed the same null (P=0.17183) and held-unit R2 is negative. Zhao is a relative-structure boundary, not an absolute-prediction anchor.
- Identity-derived molecular proxies are excluded from independent continuous molecular replication: C. elegans held-identity gene prediction collapses to median rho≈0.00055; legacy Allen VC2P gene__ columns are Cre driver-line one-hot features, not measured transcripts.

## Figure 2 — reference-qualified AtlasLift
- McLachlan exact PERI/ECT is now compared to isocortex-GABA and whole-brain ABC at identical K, selector, held-animal geometry and matched-random controls.
- Exact PERI/ECT function-guided geometry is higher than whole-brain ABC at every tested K. Formal matched-random support survives at K=20 (rho 0.20833→0.30833, Δ=+0.10000, global q=0.02994) and K=23 (rho 0.20833→0.33452, Δ=+0.12619, q=0.01996).
- Isocortex-GABA function-guided geometry is not supported after correction. Broad references can still win under atlas-variance selection, so anatomy match and selector choice remain distinct questions.
- Xu/PVH external PVN bridge is now a resolved boundary, not an unfinished optimization target. Target-free reference selection chooses external Berkhout/HypoMap C185 by measured-gene recovery; masked-gene median rho≈0.244 across 9 measured genes. Strict held-mouse residual-guided external-gene addition is non-positive; raw-guided best K=1 gives only ΔR2≈+0.00375 and does not survive matched-random/K-search control. Prespecified Mc4r gives ΔR2≈+0.0020, P≈0.129.
- Therefore measured molecular information (Fig4 conditional completion) and external atlas imputation (Fig2) remain separate claims.

## Figure 3 — fixed-budget assay design
- M1 Patch-seq remains the primary true-measured fixed-budget benchmark: K=25 joint R2≈0.3686, 84.5% of full1000; +0.042 vs PERSIST and +0.123 vs random-p95; essentially tied with sBNN at K=25.
- M1 provenance is now resolved: Scala et al., Nature 598, 144–150 (2021), DOI 10.1038/s41586-020-2907-3; GEO GSE163764; BioProject PRJNA687490; SRA SRP299088; main ephys DANDI 000008; source repository berenslab/mini-atlas.
- Allen HCR K=5 remains the strongest dataset-specific compact-panel result: 6/6 held mice improve, one-sided Wilcoxon P=0.015625, BH q=0.046875 across the fixed-budget metric family.
- WARP K=5 is directional in 3/3 held fish but unit-level P=0.125; measured-panel compression is retained, not relabeled as hidden-gene replacement.
- Wang 2023 remains the true hidden-gene replacement validation. Xu 9→new9 and Bugeon 72→new72 remain prospective when replacement genes were not measured in the original functional cells.

## Remaining boundaries, not failures
- Bugeon exact VISp has strong reference recovery and endpoint-specific utility, but the formal global geometry null remains NOT_SUPPORTED.
- EASI-PASS matching-only public material cannot support function↔molecule inference; this is a data-availability boundary.
- MICRONS is structural-only and must not be converted into a molecular claim.
- Small-N unit-level P resolution and conditional permutation P values must never be conflated.

## Current reviewer-audit sources
- results/manuscript/FIG123_REVIEWER_AUDIT_AUTHORITY_20261005.json
- results/manuscript/M1_PATCHSEQ_PROVENANCE_AUTHORITY_20261005.json
- reviewer_reanalysis_20261005/FIG1_UNIT_LEVEL_EXACT_REANALYSIS_20261005.csv
- reviewer_reanalysis_20261005/IDENTITY_DERIVED_MOLECULAR_AUDIT_20261005.json
- reviewer_reanalysis_20261005/MCLACHLAN_SAMEK_ATLAS_AUTHORITY_20261005.csv
- reviewer_reanalysis_20261005/ZHAO_ABSOLUTE_NULL_AUTHORITY_20261005.csv
