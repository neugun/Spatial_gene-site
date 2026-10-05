# Current Figure 1-3 authority — 2026-10-04

Runtime validation: 136/136 PASS (results/validate_runtime_fig123_v1paper_v7_20261004.txt).

## Figure 1
- Authority: FIG1_MULTIAXIS_AUTHORITY_V3_20261003.csv (120 evidence rows).
- R2, Pearson/Spearman, AUROC/AUPRC, relation index, delta-AIC and rank-biserial remain separate evidence axes.
- Zhao connectivity vs type: delta R2 ~+0.0011 while delta Spearman ~+0.1269, delta AUROC ~+0.0660 and delta AUPRC ~+0.1072.
- Relation-index target-validity audit now supersedes the unfiltered table: 31 old inferential rows -> 29 valid; 13 old FDR<0.05 -> 11 valid. Two MERGE standardized rows that used QC metadata rather than functional/projection response targets are excluded from inference.

## Figure 2
- Authority: FIG2_ATLASFLOW_MULTIMETRIC_AUTHORITY_20261003.csv (6 datasets; measured/atlas/measured+atlas).
- Reference QC and downstream functional recovery are separate from held-unit global geometry.
- Bugeon exact VISp: strong reference recovery, endpoint-specific multi-metric gains, but formal 1000-null global geometry remains NOT_SUPPORTED.
- McLachlan exact PERI/ECT: matched-reference global geometry remains SUPPORTED_RECOVERY.
- Zhao exact VISp: checkerboard response shows strong rank/classification lift despite modest reference-recovery rho.
- Current geometry/headroom authority separates positive and boundary cases: WARP K10 rho 0.1649->0.1991 (delta +0.0342, 3/3 fish, q=0.0053); C. elegans has significant relative recovery but remains below zero in absolute geometry; McLachlan exact PERI/ECT 0.2083->0.3202 (delta +0.1119, q=0.0400); Bugeon/Shainer/Condylis global geometry remain boundary/not supported.

## Figure 3
- Authority: FIG3_CURRENT_AUTHORITY_V3_20261003.csv.
- M1 is the primary fixed-budget benchmark; Wang provides true measured hidden-gene replacement validation.
- WARP supports measured-panel compression but not strict fixed-budget replacement; Zhao marker compression remains prospective.
- Broad fixed-budget table currently has 0/8 BH-q<0.05 rows; do not convert directional effects into a universal win claim.
- Directional replication and universal significance are kept separate: Chevee and MERGE provide held-group directional/p95 evidence under their own contracts, but this does not override the broad BH-FDR boundary.
- Public `gap fixed-budget` and high-level `gap run --goals fixed_budget` are now available; K equal to the full measured panel is explicitly non-inferential.
- McLachlan / Jerry Chen Cell Reports integration test passes Fig1, Fig2 and Fig3 end-to-end (11/11 checks), including one-command execution of all seven core Fig1-3 capabilities with zero errors.
- Allen Visual Learning 2P + HCR same-cell integration test passes 12/12 checks: 1321 unique HCR-linked neurons after outcome-blind duplicate-ROI removal / 6 mice / 21 common HCR genes / 13 primary functional features, with the same seven-capability one-command contract and zero errors.
- Allen held-mouse Fig1 after strict HCR dedup: genes vs coarse type delta Spearman ~+0.0885 (6/6 mice). Fig2 exact VISp reference is HIGH (19/21 overlap; recovery rho~0.206, p~0.0476). Fig3 K=5 beats matched-random for Spearman/AUROC/AUPRC/NRMSE at BH q~0.0469; treat this as dataset-specific evidence, not a universal selector claim.
- Condylis 2022 / Jerry Chen lab Science is the third integration benchmark (8/8 PASS): exact SSp-bfd AtlasLift works, 9 held mice execute correctly, continuous quartile effect-size is explicitly blocked for the binary marker panel, and K=6 full panel remains non-inferential.
- Integrity-safe runtime now uses fold-aware target support, scope-aware biological/session holdouts, session-scoped cell IDs, stable-ID duplicate checks (including MICrONS root_id), sparse-predictor exclusion, deterministic target ordering, and no held-out outcome imputation.
- Global integrity regression passes 14/14 datasets with 0 blocking errors; remaining WARN states are real sparse/partial endpoints or predictor coverage, not hidden execution failures.
- Five one-command boundary tests pass 5/5 with zero hidden errors: MICrONS (no transcriptome/sparse connectivity), MERGE-seq (projection-only), Condylis (binary marker panel), Allen HCR (partial event support), and Sorensen (low-N/many biological units/dual domains).
- Six ordinary core datasets also pass 6/6 one-command regression with zero manifest errors: Bugeon, Allen VC2P, Atanas C. elegans, Shainer, Xu and Zhao. Numeric sanity locks dataset-specific structure rather than forcing a universal positive result.
- Across all 14 currently audited Fig1-3 datasets, capability-aware execution coverage is 14/14 PASS with 91/91 requested capability checks and zero manifest errors.
- WARP scale contract: 306,108 input rows; coordinate quick mode uses 60,000 balanced rows (20k/fish), fixed-budget uses 150,000 balanced rows (50k/fish). K=5 Spearman gain over matched-random is +0.0461 with 3/3 fish wins but p=0.125; K=41 full panel is non-inferential.
- EASI-PASS 2026 is now the first explicit incomplete-public-data external intake test. The authors' golden JS078 matching output is independently reproduced at 1,199/1,439 IoU matches (83.3%), 1,149/1,439 confident Soma-print calls (79.85%), and median z-IoU 0.286. Matching-only input is explicitly blocked from Fig1/3 inference; a dedicated adapter, multi-animal join contract and pairwise correlation-distance stage are available for the full PBN/V1 data when supplied.
- EASI-PASS reference readiness is now closed for both branches: awake V1 uses exact VISp; PBN now uses mouse_pbn_paulichen_exact (21 neuronal subclusters x 14,437 genes). PBN reference QC is HIGH and sparse self-recovery is strong, but these are reference-headroom results, not EASI-PASS response-biology reproduction.
- EASI-PASS now exercises a three-level intake contract. Official candidate matching/QC: all seven Fig1-3 requests safe-stop with zero errors. Raw aligned molecular feature schema: automatically canonicalized (PV->Pvalb, SST->Sst), exact VISp AtlasLift runs, function-dependent stages remain blocked. Multi-animal molecular+functional same-cell input: full applicable Fig1-3 contract executes with zero errors.
- Fig1-3 now have an explicit downstream composability contract. Allen Visual Learning HCR reaches Fig1-Fig7 in one public tool chain while preserving a Fig4 saturation boundary and a negative Allen-HCR->Zhao Fig6 transfer. McLachlan/Condylis/MICRONS correct stops are frozen as PASS states rather than treated as tool failures.
- MICrONS V1 EM is the structural-only boundary benchmark (13/13 PASS + 8/8 selection-bias guard): 24,713 unique neurons, one animal/13 sessions, no genes, full-cohort connectivity correctly excluded at 4.66% coverage, and a separate 1,151-neuron connectivity-complete analysis. Space strongly recovers retinotopy (RFx rho~0.497, RFy rho~0.536). Connectivity shows local gains for several orientation/direction endpoints and a modest Pearson gain over space, but no robust broad increment after nesting on type+space(+nucleus volume). Matched per-session random controls show the connectivity-complete cohort is selection-shifted, so connectivity claims are restricted to the EM-resolved subpopulation.
- Cross-dataset integration matrix: FIG123_INDEPENDENT_INTEGRATION_BENCHMARKS_20261003.csv; append future independent datasets here rather than creating incompatible one-off summaries.

## Supersession rules
- Old AtlasFlow runs without atlas/measured+atlas consumption are obsolete.
- Old paired-delta outputs limited to R2/Spearman are obsolete.
- Whole-ABC Bugeon is retained as sensitivity/history, not current anatomy-matched Fig2 authority.

SHA256 manifest: CURRENT_FIG123_AUTHORITY_SHA256_20261003.csv.