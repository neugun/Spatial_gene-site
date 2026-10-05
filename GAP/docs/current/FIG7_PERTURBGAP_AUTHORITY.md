# Figure 7 Perturb-GAP / Virtual Neuron authority — 2026-10-01

## Biological question
Can a frozen molecular-to-functional GAP model transform realistic molecular perturbation responses into preregisterable predictions about selective changes in neural response laws, thereby choosing the next causal experiment?

## Boundary
The computation does NOT establish causality. It nominates experiments. Causal closure requires applying the nominated perturbation in the specified cell context and measuring the same functional coordinates prospectively.

## Pipeline
observed molecule + function -> learn/freeze GAP -> perturb gene in a specified cellular context -> obtain or predict the distributed transcriptomic response (Delta G) -> project Delta G through frozen GAP -> predict Delta Theta / Delta Gamma / other response-law coordinates -> rank the wet-lab experiment -> validate by intervention.

## Why target-gene-only perturbation is insufficient
Whole-brain in-vivo CRISPR Perturb-seq responses were aligned to the 70/72 GAP-compatible Bugeon genes.
For perturbations whose target gene itself is in the GAP panel, the cosine between:
1. the direct target-gene GAP Jacobian and
2. the full Perturb-seq-informed molecular-network response projected through GAP
has a median around -0.34.
Therefore changing only the target-gene coordinate is generally not an adequate model of a real perturbation response.

## Cell-context dependence
Observed whole-brain perturbation tensor:
- 2046 perturbations
- 23 cellular contexts
- 70 GAP-compatible genes
- 44,133 valid perturbation x context conditions.

Balanced two-way decomposition using the 439 perturbations observed in all 23 contexts:
Molecular response variance:
- perturbation main effect 4.23%
- context main effect 19.96%
- perturbation x context interaction 75.81%.

After projection through frozen GAP, functional-response-law variance:
- perturbation main effect 4.31%
- context main effect 10.21%
- perturbation x context interaction 85.48%.

Per functional coordinate, interaction explains:
- state-running 87.46%
- speed 72.46%
- pupil 86.43%
- network-total 89.74%
- network-instant 89.04%
- network-lagged 91.14%
- self-memory 79.83%
- behavior-gain 82.97%.

This is the main mechanistic design law for Fig7: a perturbation does not have a single context-free virtual neural effect.

## Missing-edge prediction
Five-fold held perturbation x context edges:
- global functional R2 ~ 0
- perturbation mean across other contexts R2 -0.0133
- context mean R2 0.0776
- two-way additive R2 0.0723
- low-rank residual ranks 1/2/4/8 do not improve and become progressively worse.

Therefore a simple low-rank interaction model is not currently supported. Do not promote it. The interaction is large but not trivially low-dimensional.

## Three kinds of prospective experiment nomination
1. Target-coordinate test:
Choose a perturbation x cell-context condition predicted to move one GAP coordinate while minimizing orthogonal functional effects.
Examples from the current scorecard include Syn2/Nlgn2/Sec23a for state-running, Myo6/Apc/Cntnap2 for speed, Ift43/Slc7a14/Klf13 for pupil, and Rest/Atp6v1a/Ppargc1b for self-memory. These are nominations, not causal-gene claims.

2. Model-falsification test:
Choose cases where target-gene-only and network-aware predictions disagree.
Examples:
- Htr3a in L4-5 IT CTX Glut: cosine -0.911; network-total direct effect positive but Perturb-seq-informed effect negative.
- Cck in Pvalb Gaba: cosine -0.784 with sign disagreement on network-total.
These experiments distinguish whether distributed transcriptomic response is required for functional prediction.

3. Context-dependence test:
Apply the same perturbation in two cellular contexts with highly discordant predicted functional effects.
Examples include Notch4 Pvalb vs MB Glut (cosine -0.974), Gprc5b CNU-HYa GABA vs TH/Prkcd/Grin2c Glut (-0.962), and Pax7 Pvalb vs TH/Prkcd/Grin2c Glut (-0.957).

## Link to earlier figures
Fig1 establishes discrete identity plus continuous molecular state as levels of molecular-functional organization.
Fig2 asks whether missing molecular state can be recovered.
Fig3 asks which molecular measurements are most useful.
Fig4 asks which functional manifold is molecularly readable and which gene programs organize it.
Fig5 defines structured response-law coordinates.
Fig6 tests whether those coordinates transfer.
Fig7 intervenes on the molecular side and asks which experiment should be performed next to test a prespecified change in those response-law coordinates.

## Current manuscript claim
GAP converts observed or predicted molecular perturbation responses into experimentally testable predictions about selective changes in neural response laws, enabling prospective design of the next causal experiment.

Live figure:
results/manuscript/figures_refresh_20261001_live/Fig7_PerturbGAP_hierarchical_live_20261001.{png,pdf}
