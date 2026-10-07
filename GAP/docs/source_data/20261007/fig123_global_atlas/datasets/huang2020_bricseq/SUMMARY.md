# GAP analysis bundle

Generated: 2026-10-07 12:39:36

## Dataset

- rows: 2156
- genes: 153
- functional targets: 0
- species: mouse

## Execution policy

- local-first
- Slurm is never auto-submitted

## Steps

- function_map: blocked (input requirements not met); backend=blocked; input requirements not met
- coordinate_mediation: blocked (input requirements not met); backend=blocked; input requirements not met
- atlaslift: not_applicable (no anatomically eligible registered reference could be scored); backend=local_cpu; reference scoring and AtlasLift
- fixed_budget: done (projection:30 rows/16 targets); backend=local_cpu; Fig3 held-unit fixed-budget panel design with matched-K random comparators
- projection_map: done (independent held-unit projection Ridge; r2_vw=-0.229062); backend=local_cpu; held-unit molecular-to-projection benchmark
- relation_index: done (1 eligible identity-domain levels; projection:relation=1,AIC=1,effects=2048); backend=local_cpu; Zhao-style discrete identity-response graph with within-biological-unit matched null
- information_criterion: done (1 identity-domain summaries; projection:relation=1,AIC=1,effects=2048); backend=local_cpu; same-sample type-vs-gene Gaussian AIC with explicit complexity penalty
- effect_size: done (2048 gene-target-domain effect rows; projection:relation=1,AIC=1,effects=2048); backend=local_cpu; within-biological-unit high-vs-low molecular rank-biserial effect sizes

## Data integrity

- severity: WARN
- sparse_or_degenerate_predictor_block:genes
- predictor_columns_excluded:connectivity:8
- partial_gene_panel:153
- partial_targets:5
- invalid_targets:10
- outer split: leave_one_biological_unit_out (8 folds from 8 biological units)

## Fig4-to-Fig7 continuation readiness

- Fig4: NOT_RUN
- Fig5: PROJECTION_OPERATOR_BOUNDARY; target_domain=projection, fig7_ready=False, reason=projection-domain response law is not auto-routed into functional Fig7
- Fig6: NOT_RUN
- Fig7: NOT_APPLICABLE_FUNCTION_OPERATOR; reason=Fig7 default contract requires a validated function-domain operator; projection-only datasets are kept structural

## Scientific interpretation

### Fig1-3 evidence
- Question: How strongly are molecular identity and function related under complementary Fig1-3 statistical views (relation graph, complexity-penalized AIC, and effect size)?
- Allowed claim: relation index: 1 eligible identity levels, 0 with BH q<=0.05; strongest=1.0 (type_l1); same-sample AIC: strongest median deltaAIC(type-gene)=386.6949347638874 at type_l1; gene-better fraction=1.0; rank-biserial: 2048 gene-target pairs; median absolute effect=0.0607998747625264; >=0.2 with >=0.75 sign consistency=30
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