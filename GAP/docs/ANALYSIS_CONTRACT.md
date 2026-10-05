# Common GAP analysis contract

The repository uses a common evidentiary contract so that a result does not become “replicated” merely because different datasets contribute different favorable metrics.

## A. Define the biological target before the model

For every compatible scalar functional phenotype, score several question forms under the same biological split:

- continuous magnitude: R², Pearson, Spearman, normalized RMSE/MAE, calibration slope/intercept;
- direction/sign: AUROC, AUPRC, balanced accuracy, F1, Brier score;
- responsive or large-|response| status;
- top-versus-bottom functional extremes;
- ordinal/tertile agreement.

A poor amplitude R² does not erase reproducible rank/order information, and an easy binary endpoint does not replace continuous prediction.

## B. Use a nested biological representation ladder

When available, compare with the identical outer split and downstream head:

1. context/intercept only;
2. hard class/type identity;
3. continuous measured molecular state;
4. type + continuous molecular state;
5. anatomy/space or projection, separately gated;
6. molecular + contextual coordinates;
7. state-conditioned interaction/operator.

## C. Test information below taxonomy

Whenever a meaningful hard type exists, remove the training-fold type mean from the functional target and ask whether continuous measured genes or molecular state predict the residual. The outer split remains animal/fish/worm held out.

## D. Validate molecular completion independently of function

An external atlas/reference is not accepted because it improves a downstream score. Before using lifted molecular information, audit the bridge without held-out activity whenever possible:

- masked measured-gene recovery;
- posterior/mapping entropy and support;
- reference/anatomy sensitivity;
- panel-size/headroom sweep.

External-reference values remain **inferred**. Datasets without direct expression ground truth receive `N/A`, not zero, for masked-gene recovery.

## E. Use representation-matched nulls

Examples include atlas-row shuffle within animal × hard identity, within-type expression shuffle, coordinate shuffle, projection shuffle, and behavior/condition permutation. The null destroys the correspondence being claimed while preserving dimensionality and obvious marginal structure.

## F. Biological replication is the inferential unit

Report biological-group mean/median, paired held-unit deltas, fraction of held-out units with positive gain, and bootstrap/confidence intervals where supported. A large number of cells, trials, or sessions never substitutes for independent biological replicates.

## G. Modality-specific analyses are extensions, not loopholes

Response shape/gain, operator dynamics, BehLift, projection completion, spatial experts, and transfer/sample-efficiency analyses enter only when the dataset supports those modalities. Missing modalities are marked `NOT_APPLICABLE` instead of being replaced by a convenient endpoint.

## H. Report identifiability with the result

For every main target, report the number of independent biological units, subgroup size, target reliability/headroom where available, and the effective floor imposed by the split. Repeated cells, sessions or trials improve precision but do not create new biological replicates.

A small-subgroup null is therefore distinguished from a well-powered negative result. A negative calibrated R2 on a noisy target is not interpreted the same way as a negative result on a reliable target with adequate biological replication.

## I. Give every claim an evidence status

Public results use four states: `PROMOTED`, `SUPPORTIVE`, `BOUNDARY`, or `PENDING`. The status is attached to the biological claim, not to a model name. A boundary result remains visible when it constrains the scope of a positive result.

The current status map is maintained in `CLAIM_EVIDENCE_LEDGER.md`.

## J. Do not compare across abstraction levels as one leaderboard

Population-generative realism, molecular embedding quality, behavior-aligned activity geometry, same-cell response-law inference, and transfer/sample efficiency are different endpoints. External methods are compared only when the biological target, held-unit split, dimensionality/capacity control and evaluation head make the comparison interpretable.

A specialist is allowed to win a narrow scalar endpoint. Such a result defines the boundary of the broader response-law claim rather than being hidden.
