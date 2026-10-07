# GAP analysis bundle

Generated: 2026-10-07 12:42:13

## Dataset

- rows: 1283
- genes: 23
- functional targets: 2
- species: mouse

## Execution policy

- local-first
- Slurm is never auto-submitted

## Steps

- function_map: done (2 summary rows); backend=local_cpu; held-unit Ridge/Lasso/PLS model zoo
- coordinate_mediation: not_applicable (coordinate mediation requires >=3 independent holdout groups); backend=local_cpu; same-split held-unit comparison of type, genes, space/anatomy and coordinate combinations
- atlaslift: not_applicable (no anatomically eligible registered reference could be scored); backend=local_cpu; reference scoring and AtlasLift
- fixed_budget: done (function:20 rows/2 targets; projection:20 rows/11 targets); backend=local_cpu; Fig3 held-unit fixed-budget panel design with matched-K random comparators
- projection_map: done (independent held-unit projection Ridge; r2_vw=-0.121102); backend=local_cpu; held-unit molecular-to-projection benchmark
- relation_index: blocked (input requirements not met); backend=blocked; input requirements not met
- information_criterion: blocked (input requirements not met); backend=blocked; input requirements not met
- effect_size: done (117 gene-target-domain effect rows; function:relation=0,AIC=0,effects=18; projection:relation=0,AIC=0,effects=99); backend=local_cpu; within-biological-unit high-vs-low molecular rank-biserial effect sizes

## Data integrity

- severity: WARN
- invalid_targets:1
- outer split: leave_one_biological_unit_out (2 folds from 2 biological units)

## Fig4-to-Fig7 continuation readiness

- Fig4: NOT_RUN
- Fig5: PROJECTION_OPERATOR_BOUNDARY; target_domain=projection, fig7_ready=False, reason=projection-domain response law is not auto-routed into functional Fig7
- Fig6: NOT_RUN
- Fig7: NOT_READY; reason=no qualified Fig5 operator

## Scientific interpretation

### Fig1-3 evidence
- Question: How strongly are molecular identity and function related under complementary Fig1-3 statistical views (relation graph, complexity-penalized AIC, and effect size)?
- Allowed claim: function domain: rank-biserial: 18 gene-target pairs; median absolute effect=0.0984454327589376; >=0.2 with >=0.75 sign consistency=4; projection domain: rank-biserial: 99 gene-target pairs; median absolute effect=0.0332025576513024; >=0.2 with >=0.75 sign consistency=8
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