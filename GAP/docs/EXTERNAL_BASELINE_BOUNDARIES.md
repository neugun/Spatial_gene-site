# External baseline boundaries: what each large model contributes

The external models are not arranged as a leaderboard. They operate at different biological abstraction levels, and the useful comparison is whether each closes a specific missing-information gap in the same held-biological-unit endpoint.

## Allen / Ito / Arkhipov mechanistic V1

A full bio-trained GLIF simulation contains 66,658 neurons, with 16,711 neurons in the central-core sensitivity readout. Under the standardized population-distribution benchmark, mean 1-KS is approximately 0.769 for the Allen mechanistic model versus 0.326 for XC64.

The same model remains weak for held-donor individual-cell OSI/DSI decoding: the best tested GLIF signatures retain negative R². This is not a contradiction. Wiring + physiology can reproduce a population response distribution without identifying which held-out cell will occupy which functional coordinate.

**Role in GAP/FSO:** mechanistic population-dynamics prior. It defines a complementary success criterion to same-cell molecular→functional inference.

## BrainBeacon molecular foundation encoder

The official Stage-1 checkpoint is a 1024-D, 16-layer, 16-head molecular encoder with 92,076 gene tokens. In Bugeon, 71/72 measured genes map to its mouse token space and 1,065 cells can be embedded; the reliability-qualified operator benchmark uses 376 cells across four held animals.

The raw 1024-D representation is not automatically a better same-cell decoder than the measured genes. For response-shape PC1, BrainBeacon1024 gives R² ≈ 0.111 versus ≈0.169 for raw71. The informative result appears after the same dimensionality and held-animal protocol are enforced.

## Matched low-dimensional BrainBeacon benchmark

BrainBeacon PCA4 vs raw71 PCA4 for response-shape PC1:

- 10% target labels: R² 0.0297 vs -0.0159; +0.0456, 4/4 held animals improve.
- 20%: 0.1213 vs 0.0836; +0.0377, 3/4 improve.
- 50%: 0.1750 vs 0.1511; +0.0240, 2/4 improve.
- 100%: 0.1845 vs 0.1676; +0.0169, 2/4 improve.
The conclusion is therefore narrower and more useful than “foundation model wins”: BrainBeacon contains a compact molecular prior for the dominant Theta-like response-shape axis, with the clearest advantage in low-label conditions after dimensionality is matched.

## BrainBeacon + function-aligned V2.2 core

The hybrid combines BrainBeacon PCA4 with a frozen six-seed V2.2 Z4 function-aligned core.

- 10% labels: hybrid is worse than BrainBeacon alone (-0.0022 vs 0.0297).
- 20%: 0.1331 vs 0.1213.
- 50%: 0.1859 vs 0.1750.
- 100%: 0.2066 vs 0.1845; 3/4 animals improve; rho ≈ 0.481.

The functional core therefore adds value only once enough target calibration is available. It should not be described as universally sample-efficient.

## The unresolved Gamma layer

For log-gain / Gamma-like susceptibility, absolute held-animal R² remains negative. At full labels:

- BrainBeacon PCA4: -0.1121.
- BrainBeacon + V2.2 Z4: -0.0867, 4/4 animals improve.
- BrainBeacon + V2.2 Z4 + Atlas10: -0.0745, 4/4 improve.

Adding generic molecular representation and a function-aligned core moves the result in the correct direction but does not solve Gamma **under that representation-only contract**. A later receptor-aware exact 633-cell LOAO screen reaches R² ≈0.11546 with type + BrainBeacon4 + Atlas8 + receptor6 and 4/4 held animals positive. The updated interpretation is therefore that Gamma requires a susceptibility-specific information layer; generic embedding capacity alone is insufficient.

## Architecture implied by the current evidence

The current data support a modular design:

`general molecular encoder → function-aligned core → context/receptor susceptibility layer → held-unit functional operator`

This is a scientific decomposition, not a software preference. Each block is retained only if it adds information to the biological endpoint that the preceding block does not already contain.

## Public interpretation rule

Allen GLIF, BrainBeacon, CEBRA, molecular foundation models, and GAP/FSO should be compared at their own abstraction level. Population-distribution realism, molecular embedding quality, activity/behavior geometry, and same-cell operator prediction are complementary objects. A model is not promoted as globally superior because it wins one of them.

## oPhys/Xenium as an analysis-logic stress test

The frozen `oPhys_Transcriptomics_Analysis` repository is not treated as another model in the leaderboard. It is useful because it cleanly exposes the difference between pooled molecular association and held-mouse prediction. Its native gene-to-tuning notebook reports many significant associations, while leave-one-mouse LASSO remains poorly calibrated for several scalar tuning metrics.

**Current role in GAP/FSO:** methodological stress test and independent analysis template. DANDI 001532 is not promoted as a matched GAP benchmark until persistent assets, same-cell joins, native target reproduction and identical held-mouse baseline/GAP evaluation all pass.

## Where the external baselines belong

- **Q5 / Figure 5:** Allen and BrainBeacon define abstraction-level boundaries around the response-law/operator question; they help identify whether the missing information is population mechanism, generic molecular representation, function alignment, or context/receptor susceptibility.
- **Q6 / Figure 6:** BrainBeacon, V2.2 and cross-system experts are evaluated by target-label efficiency and biological alignment.
- **Analysis discipline:** oPhys demonstrates why association, scalar calibration, rank/order and held-mouse generalization must not be conflated.

These roles should remain separate even when the same model appears in more than one diagnostic analysis.
