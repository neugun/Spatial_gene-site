# GAP analysis bundle

Generated: 2026-10-07 12:36:55

## Dataset

- rows: 306108
- genes: 41
- functional targets: 33
- species: zebrafish

## Execution policy

- local-first
- Slurm is never auto-submitted

## Steps

- function_map: done (8 summary rows); backend=local_cpu; held-unit Ridge/Lasso/PLS model zoo
- coordinate_mediation: done (7 coordinate models; 3 held groups); backend=local_cpu; same-split held-unit comparison of type, genes, space/anatomy and coordinate combinations
- atlaslift: done (zebrafish_gse158_5dpf); backend=local_cpu; reference scoring and AtlasLift
- fixed_budget: done (function:25 rows/8 targets); backend=local_cpu; Fig3 held-unit fixed-budget panel design with matched-K random comparators
- projection_map: blocked (input requirements not met); backend=blocked; input requirements not met
- relation_index: done (2 eligible identity-domain levels; function:relation=2,AIC=2,effects=24); backend=local_cpu; Zhao-style discrete identity-response graph with within-biological-unit matched null
- information_criterion: done (2 identity-domain summaries; function:relation=2,AIC=2,effects=24); backend=local_cpu; same-sample type-vs-gene Gaussian AIC with explicit complexity penalty
- effect_size: done (24 gene-target-domain effect rows; function:relation=2,AIC=2,effects=24); backend=local_cpu; within-biological-unit high-vs-low molecular rank-biserial effect sizes

## Data integrity

- severity: PASS
- outer split: leave_one_biological_unit_out (3 folds from 3 biological units)

## Fig4-to-Fig7 continuation readiness

- Fig4: NOT_RUN
- Fig5: NOT_RUN
- Fig6: NOT_RUN
- Fig7: NOT_READY; reason=no qualified Fig5 operator

## Scientific interpretation

### Fig1-3 evidence
- Question: How strongly are molecular identity and function related under complementary Fig1-3 statistical views (relation graph, complexity-penalized AIC, and effect size)?
- Allowed claim: relation index: 2 eligible identity levels, 2 with BH q<=0.05; strongest=0.0032622543657585 (type_l2); same-sample AIC: strongest median deltaAIC(type-gene)=-970.3929461390944 at type_l1; gene-better fraction=0.0; rank-biserial: 24 gene-target pairs; median absolute effect=0.0198227153488785; >=0.2 with >=0.75 sign consistency=0
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
- Allowed claim: Fig7 status: NOT_READY.
- Guardrail: Prospective nomination requires a held-unit-qualified molecular operator and compatible perturbation gene coverage.


## Provenance

- Python: 3.14.6
- platform: Windows-11-10.0.26200-SP0