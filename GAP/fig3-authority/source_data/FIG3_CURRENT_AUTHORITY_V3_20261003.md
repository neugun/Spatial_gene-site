# Figure 3 current authority — 2026-10-03

Question: under a fixed molecular-assay budget, which coordinates should be measured, and what kind of validation supports the choice?

## Frozen interpretation

- Primary fixed-budget benchmark (M1): FSO-functional K25 joint R2=0.369, retains 84.5% of full1000, delta vs PERSIST=+0.042, delta vs random-p95=+0.123. K50 retains 87.2%.
- True hidden-gene validation (Wang2023): selected hidden genes beat matched random in AUROC advantage for 5/5 targets; advantage range +0.011 to +0.035. This is stronger than prospective AtlasLift nomination because the hidden genes are actually measured.
- WARP measured compression: K25 AUROC 0.556 vs full31 0.555; within-fish Spearman 0.104 vs 0.104. Compression can preserve the assay.
- WARP strict replacement: fixed-budget bridge/functional replacements do not beat full31 AUROC; adding functional genes reaches AUROC 0.559. Compression and replacement must not be conflated.
- Bugeon selector dependence: no universal selector is assumed; operator and shape objectives can prefer different panels and the preferred method changes with K.
- Zhao marker budget: no registered marker-compression q-value is <0.05; retain as prospective design evidence, not true hidden-gene replacement proof.
- Broad replication: 0/8 rows pass BH q<0.05 in the current broad table. Chevee is 6/6 held folds positive but broad-table q remains 0.125; report effect/direction without overstating FDR support.

## Guardrails

- Never rank selectors with one aggregate score; report objective, K, held-unit split, random comparator and native metric.
- Distinguish measured-panel compression, true measured hidden-gene replacement, inferred-gene nomination and prospective marker compression.
- Foundation low-label gains are auxiliary capability evidence unless they replicate across datasets under the same label-fraction contract.

Authority CSV: results/manuscript/FIG3_CURRENT_AUTHORITY_V3_20261003.csv