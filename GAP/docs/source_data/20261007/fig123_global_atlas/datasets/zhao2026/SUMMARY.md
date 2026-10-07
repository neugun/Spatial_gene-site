# GAP analysis bundle

Generated: 2026-10-07 12:32:24

## Dataset

- rows: 141
- genes: 30
- functional targets: 68
- species: mouse

## Execution policy

- local-first
- Slurm is never auto-submitted

## Steps

- function_map: done (5 summary rows); backend=local_cpu; held-unit Ridge/Lasso/PLS model zoo
- coordinate_mediation: done (16 coordinate models; 6 held groups); backend=local_cpu; same-split held-unit comparison of type, genes, space/anatomy and coordinate combinations
- atlaslift: done (mouse_visp_exact); backend=local_cpu; reference scoring and AtlasLift
- fixed_budget: done (function:25 rows/5 targets; projection:25 rows/16 targets); backend=local_cpu; Fig3 held-unit fixed-budget panel design with matched-K random comparators
- projection_map: done (independent held-unit projection Ridge; r2_vw=-0.0267287); backend=local_cpu; held-unit molecular-to-projection benchmark
- relation_index: done (2 eligible identity-domain levels; function:relation=1,AIC=1,effects=150; projection:relation=1,AIC=1,effects=480); backend=local_cpu; Zhao-style discrete identity-response graph with within-biological-unit matched null
- information_criterion: done (2 identity-domain summaries; function:relation=1,AIC=1,effects=150; projection:relation=1,AIC=1,effects=480); backend=local_cpu; same-sample type-vs-gene Gaussian AIC with explicit complexity penalty
- effect_size: done (630 gene-target-domain effect rows; function:relation=1,AIC=1,effects=150; projection:relation=1,AIC=1,effects=480); backend=local_cpu; within-biological-unit high-vs-low molecular rank-biserial effect sizes

## Data integrity

- severity: WARN
- predictor_columns_excluded:morphology:2
- partial_targets:12
- invalid_targets:30
- outer split: group_kfold_whole_biological_units (6 folds from 9 biological units)

## Fig4-to-Fig7 continuation readiness

- Fig4: NOT_RUN
- Fig5: PROJECTION_OPERATOR_BOUNDARY; target_domain=projection, fig7_ready=False, reason=projection-domain response law is not auto-routed into functional Fig7
- Fig6: NOT_RUN
- Fig7: NOT_READY; reason=no qualified Fig5 operator

## Scientific interpretation

### Fig1-3 evidence
- Question: How strongly are molecular identity and function related under complementary Fig1-3 statistical views (relation graph, complexity-penalized AIC, and effect size)?
- Allowed claim: function domain: relation index: 1 eligible identity levels, 0 with BH q<=0.05; strongest=0.0 (type_l1); function domain: same-sample AIC: strongest median deltaAIC(type-gene)=-22.427549086286035 at type_l1; gene-better fraction=0.0; function domain: rank-biserial: 150 gene-target pairs; median absolute effect=0.0984898589065255; >=0.2 with >=0.75 sign consistency=25; projection domain: relation index: 1 eligible identity levels, 0 with BH q<=0.05; strongest=0.6666666666666666 (type_l1); projection domain: same-sample AIC: strongest median deltaAIC(type-gene)=-23.450759232397395 at type_l1; gene-better fraction=0.0625; projection domain: rank-biserial: 480 gene-target pairs; median absolute effect=0.1324404761904762; >=0.2 with >=0.75 sign consistency=146
- Guardrail: Relation index, AIC, and rank-biserial answer different questions. Positive deltaAIC(type-gene) favors genes only for same-sample models with explicit parameter counts; rank-biserial magnitude is not by itself a multiple-testing-significant gene claim.

### Fig4
- Question: Do molecular measurements add unique information about a missing functional coordinate after the other measured functions are known?
- Allowed claim: Fig4 completion was not executed in this bundle (status=NOT_RUN).
- Guardrail: Do not interpret an unexecuted or preflight-blocked stage as a negative biological result.

### Fig5
- Question: Can molecular state predict a multivariate functional response law across held biological units, and does shared structure improve on independent prediction?
- Allowed claim: The held-unit operator does not meet the prespecified prospective-use gate.
- Guardrail: Stop before Fig7 rather than generating virtual perturbation rankings from an unsupported law.

### Fig6
- Question: Does the frozen GAP representation add information beyond raw measured genes, or support a prespecified cross-dataset transfer?
- Allowed claim: Foundation/generalization stage is NOT_RUN.
- Guardrail: Do not lower gene/region preflight thresholds or invent semantic endpoint pairs merely to force a Fig6 result.

### Fig7
- Question: Is there a held-unit-qualified molecular-to-function law that can be used for prospective perturbation ranking across cellular contexts?
- Allowed claim: Fig7 status: NOT_READY.
- Guardrail: Prospective nomination requires a held-unit-qualified molecular operator and compatible perturbation gene coverage.


## Provenance

- Python: 3.14.6
- platform: Windows-11-10.0.26200-SP0