# GAP analysis bundle

Generated: 2026-10-07 12:25:11

## Dataset

- rows: 740
- genes: 0
- functional targets: 48
- species: fly

## Execution policy

- local-first
- Slurm is never auto-submitted

## Steps

- function_map: blocked (input requirements not met); backend=blocked; input requirements not met
- coordinate_mediation: done (8 coordinate models; 6 held groups); backend=local_cpu; same-split held-unit comparison of type, genes, space/anatomy and coordinate combinations
- atlaslift: blocked (input requirements not met); backend=blocked; input requirements not met
- fixed_budget: blocked (input requirements not met); backend=blocked; input requirements not met
- projection_map: blocked (input requirements not met); backend=blocked; input requirements not met
- relation_index: done (4 eligible identity-domain levels; function:relation=2,AIC=0,effects=0; projection:relation=2,AIC=0,effects=0); backend=local_cpu; Zhao-style discrete identity-response graph with within-biological-unit matched null
- information_criterion: blocked (input requirements not met); backend=blocked; input requirements not met
- effect_size: blocked (input requirements not met); backend=blocked; input requirements not met

## Data integrity

- severity: PASS
- outer split: group_kfold_whole_biological_units (6 folds from 20 biological units)

## Fig4-to-Fig7 continuation readiness

- Fig4: STRUCTURAL_ONLY; reason=no transcriptomic gene block
- Fig5: NOT_APPLICABLE; reason=no molecular operator can be defined
- Fig6: NOT_RUN
- Fig7: NOT_APPLICABLE; reason=no validated molecular operator

## Scientific interpretation

### Fig1-3 evidence
- Question: How strongly are molecular identity and function related under complementary Fig1-3 statistical views (relation graph, complexity-penalized AIC, and effect size)?
- Allowed claim: function domain: relation index: 2 eligible identity levels, 2 with BH q<=0.05; strongest=0.0958900784113805 (type_l1); function domain: rank-biserial: N/A (no_gene_columns; 0 valid within-unit quartile pairs); projection domain: relation index: 2 eligible identity levels, 2 with BH q<=0.05; strongest=0.603658536585366 (type_l2); projection domain: rank-biserial: N/A (no_gene_columns; 0 valid within-unit quartile pairs)
- Guardrail: Relation index, AIC, and rank-biserial answer different questions. Positive deltaAIC(type-gene) favors genes only for same-sample models with explicit parameter counts; rank-biserial magnitude is not by itself a multiple-testing-significant gene claim.

### Fig4
- Question: Do molecular measurements add unique information about a missing functional coordinate after the other measured functions are known?
- Allowed claim: This dataset supports structural/functional geometry analysis but not molecular functional completion.
- Guardrail: Do not convert structural coordinates into a gene-to-function claim.

### Fig5
- Question: Can molecular state predict a multivariate functional response law across held biological units, and does shared structure improve on independent prediction?
- Allowed claim: No molecular multivariate operator is scientifically defined for this dataset.
- Guardrail: Do not fabricate a gene-to-response operator from non-molecular structural labels.

### Fig6
- Question: Does the frozen GAP representation add information beyond raw measured genes, or support a prespecified cross-dataset transfer?
- Allowed claim: Foundation/generalization stage is NOT_RUN.
- Guardrail: Do not lower gene/region preflight thresholds or invent semantic endpoint pairs merely to force a Fig6 result.

### Fig7
- Question: Is there a held-unit-qualified molecular-to-function law that can be used for prospective perturbation ranking across cellular contexts?
- Allowed claim: Fig7 correctly stops: no validated molecular operator.
- Guardrail: A scientifically correct stop is a tool PASS and should remain visible as a boundary.


## Provenance

- Python: 3.14.6
- platform: Windows-11-10.0.26200-SP0