# GAP analysis bundle

Generated: 2026-10-07 12:30:40

## Dataset

- rows: 33139
- genes: 29
- functional targets: 0
- species: mouse

## Execution policy

- local-first
- Slurm is never auto-submitted

## Steps

- function_map: blocked (input requirements not met); backend=blocked; input requirements not met
- coordinate_mediation: blocked (input requirements not met); backend=blocked; input requirements not met
- atlaslift: not_applicable (no anatomically eligible registered reference could be scored); backend=local_cpu; reference scoring and AtlasLift
- fixed_budget: done (projection:25 rows/5 targets); backend=local_cpu; Fig3 held-unit fixed-budget panel design with matched-K random comparators
- projection_map: done (independent held-unit projection Ridge; r2_vw=-0.0977188); backend=local_cpu; held-unit molecular-to-projection benchmark
- relation_index: done (1 eligible identity-domain levels; projection:relation=1,AIC=1,effects=145); backend=local_cpu; Zhao-style discrete identity-response graph with within-biological-unit matched null
- information_criterion: done (1 identity-domain summaries; projection:relation=1,AIC=1,effects=145); backend=local_cpu; same-sample type-vs-gene Gaussian AIC with explicit complexity penalty
- effect_size: done (145 gene-target-domain effect rows; projection:relation=1,AIC=1,effects=145); backend=local_cpu; within-biological-unit high-vs-low molecular rank-biserial effect sizes

## Data integrity

- severity: PASS
- outer split: leave_one_biological_unit_out (2 folds from 2 biological units)

## Fig4-to-Fig7 continuation readiness

- Fig4: NOT_RUN
- Fig5: PROJECTION_OPERATOR_BOUNDARY; target_domain=projection, fig7_ready=False, reason=projection-domain response law is not auto-routed into functional Fig7
- Fig6: NOT_RUN
- Fig7: NOT_APPLICABLE_FUNCTION_OPERATOR; reason=Fig7 default contract requires a validated function-domain operator; projection-only datasets are kept structural

## Scientific interpretation

### Fig1-3 evidence
- Question: How strongly are molecular identity and function related under complementary Fig1-3 statistical views (relation graph, complexity-penalized AIC, and effect size)?
- Allowed claim: relation index: 1 eligible identity levels, 0 with BH q<=0.05; strongest=0.004122092452645 (type_l2); same-sample AIC: strongest median deltaAIC(type-gene)=19.60167791316053 at type_l2; gene-better fraction=0.6; rank-biserial: 145 gene-target pairs; median absolute effect=0.0048378702557493; >=0.2 with >=0.75 sign consistency=0
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
- Allowed claim: Fig7 correctly stops: Fig7 default contract requires a validated function-domain operator; projection-only datasets are kept structural.
- Guardrail: A scientifically correct stop is a tool PASS and should remain visible as a boundary.


## Provenance

- Python: 3.14.6
- platform: Windows-11-10.0.26200-SP0