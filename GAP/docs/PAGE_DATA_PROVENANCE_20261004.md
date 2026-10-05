# Browser-page provenance and claim ledger — 2026-10-04

**Current authority:** `FIG123-current-20261004-v18-external-paper-native`
**Runtime validation:** 136/136 PASS, 0 failed.

Every browser-facing quantitative claim is backed by a compact file under `docs/source_data/20261004/`. Large raw datasets are not redistributed; original public accessions remain in `docs/DATA_CATALOG.md`.

| Browser claim | Displayed value | Source file | Inference / boundary |
|---|---:|---|---|
| Bugeon gene rank signal | aggregate OOF rho = 0.398 | `BUGEON_COORDINATE_MODELS.csv` + `BUGEON_COORDINATE_FOLDS.csv` | 4 held mice; thin lines are mice, diamond is frozen aggregate OOF |
| Xu genes vs coarse type | aggregate OOF 0.338 vs 0.201 | `XU_COORDINATE_MODELS.csv` + `XU_COORDINATE_FOLDS.csv` | 3 held mice |
| Shainer space vs genes | aggregate OOF 0.261 vs 0.032 | `SHAINER_COORDINATE_MODELS.csv` + `SHAINER_COORDINATE_FOLDS.csv` | 6 held fish; spatial regime |
| Zhao metric dissociation | aggregate OOF delta R2 ~ +0.0011; delta Spearman ~ +0.1269 | `FIG1_MULTIAXIS_AUTHORITY_V3.csv` + `ZHAO_COORDINATE_FOLDS.csv` | 6 held animals shown individually; aggregate OOF and mean-fold R2 are deliberately not conflated |
| McLachlan exact-reference geometry | delta rho +0.112; 2/3 folds improve; q = 0.03996 | `FIG2_ATLASLIFT_GEOMETRY_HEADROOM.csv` | exact PERI/ECT; n=3 fold limit retained |
| Allen HCR K=5 fixed budget | delta Spearman +0.0967; 6/6 mice; p=0.015625; q=0.046875 | `ALLEN_HCR_FIXED_BUDGET_FOLDS.csv` + `ALLEN_HCR_FIXED_BUDGET_SUMMARY.csv` | paired selected-vs-matched-random inference across 6 held mice |
| M1 K=25 fixed budget | R2 = 0.369; 84.5% of full; mean±SD across n=10 folds | `M1_FIXED_BUDGET_MULTIMETRIC.csv` | **upstream publication/accession not yet frozen**; result visible, release-provenance gate remains open |
| WARP measured compression | K25 AUROC 0.556 vs full31 0.555 | `WARP_MEASURED_PANEL_COMPRESSION.csv` | descriptive measured-panel compression, not hidden-gene replacement |
| WARP K=5 fixed-budget companion | 3/3 fish improve; one-sided p=0.125 | `WARP_FIXED_BUDGET_FOLDS.csv` | paired held-fish evidence; directional, not significant with n=3 |
| Wang 2023 true hidden validation | selected hidden genes beat matched random for 5/5 projection targets | `WANG_TRUE_HIDDEN_GENE_REPLACEMENT.csv` | target-wise measured-hidden validation; measured EASI-FISH source Figshare 10.25378/janelia.21171373; independent reference GSE213828 |
| EASI-PASS method check | 1199/1439 = 83.3% | `EASIPASS_GOLDEN_MATCHING_SUMMARY.json` | direct method reproduction only |
| V1 pairwise fixture | within rho -0.812; across rho -0.828 | `EASIPASS_V1_PAIRWISE_SUMMARY.csv` | engineering fixture; not paper biology |

## Claim-level labels used on the page

- **Direct reproduction** = authors' released data + matching metric definition.
- **Held-unit evidence** = primary inference is animal/fish/worm/donor, not cells.
- **Reference headroom** = self-consistency of an external atlas/reference, never functional-biology proof.
- **Engineering contract** = route and guardrail validation on fixtures.
- **Boundary / safe stop** = scientifically insufficient or negative result retained rather than hidden.

## Superseded presentation

The old `docs/figures/v3/` graphics remain in Git history/provenance. They are **not** the current browser-facing authority after this refresh.


## Figure 4-7 browser claims

| Browser claim | Displayed value | Public compact source | Boundary |
|---|---:|---|---|
| Fig4 C. elegans pumping completion | delta R2 +0.117; delta Spearman +0.214; delta AUROC +0.115 | `FIG4_7_STAGE_METRIC_LEDGER.csv` | held biological units; not universal |
| Fig4 Sorensen completion | delta R2 +0.0777; delta Spearman +0.451; delta AUROC +0.273 | `FIG4_7_STAGE_METRIC_LEDGER.csv` | independent low-dimensional dataset |
| Fig4 Xu Ghrelin completion | delta R2 +0.0639 | `FIG4_7_STAGE_METRIC_LEDGER.csv` | 3 held mice |
| Fig5 Bugeon exact 24D operator | R2 0.0093 -> 0.0185 | `FIG4_7_STAGE_METRIC_LEDGER.csv` | exact dataset-specific contract |
| Fig5 Sorensen low-rank replication | delta R2 +0.0116; delta Spearman +0.103 | `FIG4_7_STAGE_METRIC_LEDGER.csv` | shared law not universal |
| Fig5 structural branch | Chevee R2 0.677; worm connectivity R2 0.781; MERGE R2 0.222 | `FIG4_7_ARTICLE_ROLE_MATRIX.csv` | structural predictability != shared-operator superiority |
| Fig6 Zhao->Bugeon 10% labels | delta R2 +0.1546; delta Spearman +0.0781; delta AUROC +0.0632 | `FIG6_MASTER_AUTHORITY_V51.md` | same-species target-calibrated transfer |
| Fig6 Allen HCR->Zhao semantic axis | source OOF rho ~+0.305; target rho ~-0.192 | `FIG4_7_END2END_TOOL_V2_AUTHORITY.md` | **no transfer claim** |
| Fig7 interaction law | molecular 75.81%; projected functional 85.48% | `FIG7_PERTURBGAP_AUTHORITY.md` | prospective design, not causal proof |
| Allen HCR full chain | 44,133 conditions; 7,781 eligible; 2,046 perturbations; 23 contexts | `FIG4_7_END2END_TOOL_V2_AUTHORITY.md` | Fig7 runs only after held-mouse gate |
