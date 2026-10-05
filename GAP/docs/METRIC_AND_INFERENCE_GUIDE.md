# Metric and inference guide: what each number means biologically

GAP reports several metrics because “does molecular identity predict function?” is not one statistical question. The same cells can preserve response order while failing absolute cross-animal calibration, or separate extreme responders while explaining little continuous variance.

| Readout | Biological question | Interpretation |
|---|---|---|
| **R²** | Can the model reproduce the calibrated magnitude of the held-out phenotype? | Sensitive to cross-animal offsets and scale. A negative value means the predictor is worse than the held-out mean for calibrated amplitude; it does not automatically imply absence of rank/order information. |
| **Pearson r** | Does predicted magnitude vary linearly with the measured phenotype? | Less strict than calibrated R² because a linear scale/offset mismatch can leave correlation intact. |
| **Spearman ρ** | Are cells placed in approximately the correct functional order? | Useful when absolute amplitude is noisy but relative position is reproducible. |
| **AUROC** | Can molecular information separate two biologically defined response groups? | Used for sign, responsive/nonresponsive, or top-vs-bottom extremes. It asks an ordering/discrimination question, not amplitude calibration. |
| **AUPRC** | How well are positive cells recovered when classes are imbalanced? | Reported beside AUROC when prevalence differs strongly. |
| **Balanced accuracy** | Are all classes predicted rather than only the majority class? | Important for molecular identity decoding and imbalanced response categories. |
| **Normalized RMSE/MAE** | How large are prediction errors relative to the phenotype scale? | Complements R² with an error-size view. |
| **Operator/kernel score** | Does molecular identity improve prediction of a state-dependent response law? | Used only when the response-law coordinate itself passes a reliability/reproducibility gate. |

## The biological unit is the unit of inference

Thousands of cells from one animal do not equal thousands of independent biological replicates. Main analyses therefore hold out mouse, fish, worm, donor, or another acquisition-level biological unit whenever the dataset permits it.

Hyperparameters and feature selection are learned only inside the outer-training biological units. Cell-level bootstraps can describe uncertainty conditional on the sampled animals, but they are not relabeled as population-of-animal confidence intervals.

## Association is not held-unit prediction

A gene-function association estimated after pooling cells answers a different question from prediction in a held-out animal. The oPhys/Xenium audit is a useful external example: many gene-metric associations can be detected within the observed cohort while leave-one-mouse prediction of scalar tuning metrics remains poorly calibrated. GAP therefore does not promote pooled association as evidence of biological-unit generalization.

Likewise, a significant within-cohort correlation is not substituted for a held-unit R2, rank metric, or classification result. Association, calibration, ordering and discrimination are reported as distinct evidence layers.

## Matched nulls preserve everything except the hypothesized correspondence

A useful null destroys the biological link being tested while preserving obvious nuisance structure. Examples include:

- shuffle expression among cells within the same animal × hard type to test within-type molecular information;
- shuffle cell-to-atlas posterior rows while preserving atlas dimensionality and marginal distributions;
- shuffle spatial coordinates while preserving molecular measurements and sample composition;
- shuffle projection labels within the appropriate biological strata;
- permute behavior/condition labels only at the level allowed by the experimental design.

A generic global shuffle is often too easy and can exaggerate evidence.

## Within-type residual analysis

Hard type is first treated as an informative baseline, not as a nuisance to remove blindly. Training-fold type means are estimated without the held animal. Molecular and functional residuals are then defined relative to those training-derived type means. The question becomes: **among cells already assigned the same hard type, does continuous molecular variation predict continuous functional variation in a new animal?**

This is the cleanest test of information that a categorical taxonomy cannot contain by construction.

## Reliability and headroom

A functional phenotype cannot be predicted more reproducibly than it is itself measured. Split-half/repeat stability is therefore reported when available. A small raw R² can be meaningful near a low reliability ceiling, while the same R² can be weak evidence for a highly reproducible phenotype.

Similarly, a rich measured molecular panel can leave little headroom for an external atlas. A failed lift after a strong measured baseline is not equivalent to a failed molecular bridge. GAP keeps bridge quality, molecular headroom, and downstream functional gain as separate quantities.

## Reporting rule

Every main result should state the target, biological holdout unit, molecular-evidence tier, immediately simpler baseline, matched null/control, and metric family. The interpretation should say explicitly whether the result concerns **amplitude, order, classification, response geometry, or state susceptibility**.

