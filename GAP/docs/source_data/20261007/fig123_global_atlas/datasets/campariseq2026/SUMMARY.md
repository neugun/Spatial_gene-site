# GAP analysis bundle

Generated: 2026-10-07 12:24:35

## Dataset

- rows: 15559
- genes: 1500
- functional targets: 1
- species: zebrafish

## Execution policy

- local-first
- Slurm is never auto-submitted

## Steps

- function_map: blocked (input requirements not met); backend=blocked; input requirements not met
- coordinate_mediation: blocked (input requirements not met); backend=blocked; input requirements not met
- atlaslift: done (zebrafish_gse269_k120); backend=local_cpu; reference scoring and AtlasLift
- fixed_budget: blocked (input requirements not met); backend=blocked; input requirements not met
- projection_map: blocked (input requirements not met); backend=blocked; input requirements not met
- relation_index: blocked (input requirements not met); backend=blocked; input requirements not met
- information_criterion: blocked (input requirements not met); backend=blocked; input requirements not met
- effect_size: blocked (input requirements not met); backend=blocked; input requirements not met

## Data integrity

- severity: WARN
- invalid_targets:1
- outer split: leave_one_biological_unit_out (4 folds from 4 biological units)

## Fig4-to-Fig7 continuation readiness

- Fig4: NOT_RUN
- Fig5: NOT_RUN
- Fig6: NOT_RUN
- Fig7: NOT_APPLICABLE_FUNCTION_OPERATOR; reason=Fig7 default contract requires a validated function-domain operator; projection-only datasets are kept structural

## Scientific interpretation

### Fig1-3 evidence
- Question: How strongly are molecular identity and function related under complementary Fig1-3 statistical views (relation graph, complexity-penalized AIC, and effect size)?
- Allowed claim: Fig1-3 evidence modules were not run in this bundle.
- Guardrail: Do not interpret an unexecuted evidence module as a null result.

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