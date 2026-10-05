# Fig4-7 end-to-end tool V2 authority - 2026-10-04

## Core result
Fig4-7 now operates as a continuation of the Fig1-3 tool contract rather than as separate analysis scripts. The public gap run can optionally accept both a cross-dataset generalization target and a perturbation-DE atlas. It reports stage readiness, scientific interpretation, and explicit stop reasons.

### Allen Visual Learning HCR full-chain
- Fig4: saturation/mixed boundary; no inflated completion claim.
- Fig5: perturbation-compatible 15-gene operator remains held-mouse predictive (R2 about 0.156, Spearman about 0.357, AUROC about 0.742); deployment is full-rank, but superiority over independent Ridge is not claimed.
- Fig6 within-dataset: raw+GAP-Large adds mean R2 about +0.00689 and is positive on 8/9 targets.
- Fig6 cross-dataset HCR->Zhao visual-evoked axis: source OOF rho about 0.305 but target rho about -0.192, so status is EXECUTED_NO_TARGET_TRANSFER and claim_supported=false.
- Fig7: 44,133 perturbation-context conditions, 7,781 eligible, 2,046 perturbations, 23 contexts.

### Sorensen
- Fig4 completion support: best delta R2 about +0.0777.
- Fig5 deployment selection correctly chooses rank-1 RRR: R2 about 0.0548, Spearman about 0.148, AUROC about 0.567; held target-fold R2/Spearman win fractions are both 0.70.
- Fig7 does not run because the five molecular coordinates are txPC01-05 latent transcriptomic PCs rather than perturbation-atlas gene symbols; overlap is 0/5. This is a namespace boundary, not an operator failure.

### Tool rules now frozen
- Fig4 biology-aware profiles prevent same-axis derivative leakage.
- Fig5 deployment is selected source-only between full-rank and selected-RRR by pooled held-unit OOF R2.
- Fig7 requires absolute operator support plus held-unit stability, and separately checks perturbation gene namespace/coverage.
- Fig6 route/hyperparameter selection is source-only. Target labels are evaluation-only, but when target labels exist they must support the final transfer/conservation claim; source success alone is insufficient.
- Fig4/5 validated-stage outer CV remains frozen at exact LOGO for <=10 groups and deterministic 5-fold GroupKFold otherwise.

## QA
- Full frozen Fig4/5 regression: 21/21 PASS.
- End-to-end contract hardening: 21/21 PASS.

Machine-readable stage authority: analysis_workspace/results/manuscript\FIG4_7_END2END_TOOL_V2_AUTHORITY_20261004.csv
QA authority: analysis_workspace/results/manuscript\FIG4_7_END2END_TOOL_V2_QA_20261004.csv