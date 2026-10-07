# GAP analysis bundle

Generated: 2026-10-07 12:28:31

## Dataset

- rows: 673
- genes: 1500
- functional targets: 1
- species: mouse

## Execution policy

- local-first
- Slurm is never auto-submitted

## Steps

- function_map: done (1 summary rows); backend=local_cpu; held-unit Ridge/Lasso/PLS model zoo
- coordinate_mediation: done (4 coordinate models; 6 held groups); backend=local_cpu; same-split held-unit comparison of type, genes, space/anatomy and coordinate combinations
- atlaslift: not_applicable (no anatomically eligible registered reference could be scored); backend=local_cpu; reference scoring and AtlasLift
- fixed_budget: done (function:30 rows/1 targets); backend=local_cpu; Fig3 held-unit fixed-budget panel design with matched-K random comparators
- projection_map: blocked (input requirements not met); backend=blocked; input requirements not met
- relation_index: blocked (input requirements not met); backend=blocked; input requirements not met
- information_criterion: done (0 identity-domain summaries; function:relation=0,AIC=0,effects=128); backend=local_cpu; same-sample type-vs-gene Gaussian AIC with explicit complexity penalty
- effect_size: done (128 gene-target-domain effect rows; function:relation=0,AIC=0,effects=128); backend=local_cpu; within-biological-unit high-vs-low molecular rank-biserial effect sizes

## Data integrity

- severity: PASS
- outer split: group_kfold_whole_biological_units (6 folds from 17 biological units)

## Fig4-to-Fig7 continuation readiness

- Fig4: NOT_RUN
- Fig5: NOT_RUN
- Fig6: NOT_RUN
- Fig7: NOT_APPLICABLE_FUNCTION_OPERATOR; reason=Fig7 default contract requires a validated function-domain operator; projection-only datasets are kept structural

## Scientific interpretation

### Fig1-3 evidence
- Question: How strongly are molecular identity and function related under complementary Fig1-3 statistical views (relation graph, complexity-penalized AIC, and effect size)?
- Allowed claim: rank-biserial: 128 gene-target pairs; median absolute effect=0.0979506604506605; >=0.2 with >=0.75 sign consistency=5
- Guardrail: Relation index, AIC, and rank-biserial answer different questions. Positive deltaAIC(type-gene) favors genes only for same-sample models with explicit parameter counts; rank-biserial magnitude is not by itself a multiple-testing-significant gene claim.

### Fig4
- Question: Do molecular measurements add unique information about a missing functional coordinate after the other measured functions are known?
- Allowed claim: Fig4 completion was not executed in this bundle (status=NOT_RUN).
- Guardrail: Do not interpret an unexecuted or preflight-blocked stage as a negative biological result.

### Fig5
- Question: Can molecular state predict a multivariate functional response law across held biological units, and does shared structure improve on independent prediction?
- Allowed claim: Fig5 operator was not executed in this bundle (status=NOT_RUN).
- Guardrail: Do not interpret an unexecuted operator stage as failure of the response-law hypothesis.

### Fig6
- Question: Does the frozen GAP representation add information beyond raw measured genes, or support a prespecified cross-dataset transfer?
- Allowed claim: Foundation/generalization stage is NOT_RUN.
- Guardrail: Do not lower gene/region preflight thresholds or invent semantic endpoint pairs merely to force a Fig6 result.

### Fig7
- Question: Is there a held-unit-qualified molecular-to-function law that can be used for prospective perturbation ranking across cellular contexts?
- Allowed claim: Fig7 correctly stops: Fig7 default contract requires a validated function-domain operator; projection-only datasets are kept structural.
- Guardrail: A scientifically correct stop is a tool PASS and should remain visible as a boundary.


## Provenance

- Python: 3.14.6
- platform: Windows-11-10.0.26200-SP0