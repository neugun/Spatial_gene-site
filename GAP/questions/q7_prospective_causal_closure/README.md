# Prospective Figure 7 — Can a molecular perturbation move the predicted functional coordinate?

> **Status: synthetic/prospective only.** This page defines the causal falsification test. It is not experimental evidence that GAP/FSO has already established causality.

![Prospective Figure 7](../../docs/figures/v3/Fig7_REQUIREMENT_V3.png)

## Biological question

The observational model nominates a molecular/controller coordinate and predicts how perturbing it should change the stable response geometry (`Theta`) or state susceptibility (`Gamma`). Figure 7 asks whether that prediction is specific enough to be falsified experimentally.

## Preregistered causal test

A real Perturb-FSO experiment should:

1. nominate a molecular handle from training data only;
2. specify before perturbation which operator coordinate should move;
3. perturb that handle in a new biological replicate;
4. measure `Delta Theta`, `Delta Gamma`, and the state-dependent activity response;
5. require the predicted coordinate to move while orthogonal coordinates remain comparatively preserved;
6. compare with matched perturbation, expression, state and anatomical controls.

## Synthetic sanity check

The current six-animal synthetic system has known causal ground truth.

For a Gamma-only intervention:
- true `Delta Gamma = 1.547`;
- predicted `Delta Gamma = 1.545`;
- predicted orthogonal `Delta Theta = 0.014`;
- counterfactual response `R2 = 0.997`;
- Gamma-change cosine `= 0.99999`.

A complementary Theta-only intervention gives counterfactual response `R2 = 0.998`.

These numbers validate the **logic of the falsification test**, not a biological mechanism.

## What would count as support?

Support requires coordinate-selective movement in held biological units, with the direction and state dependence preregistered from the observational operator.

A generic activity change, a broad cFos increase, or a post-hoc axis chosen after the perturbation would not close the causal loop.

## What would falsify the model?

The causal interpretation would be weakened if:
- the nominated perturbation changes the wrong coordinate;
- the effect does not generalize across biological replicates;
- orthogonal response-law coordinates change as much as the predicted coordinate;
- a matched molecular/anatomical control produces the same effect;
- the observed perturbation response lies outside the preregistered operator prediction.

## Relationship to Figures 1–6

Figures 1–6 establish what is reproducible, molecularly readable, recoverable, worth measuring, operator-structured and reusable. Prospective Figure 7 asks the final question: **does perturbing the nominated molecular/controller variable move the exact functional coordinate predicted by the model?**

Until real perturbation data pass this test, the public claim remains predictive rather than causal.


## Current 2026-10-04 empirical design engine

Figure 7 now also has a real perturbation-DE design layer in addition to the synthetic causal sanity check. A held-unit-qualified Fig5 operator projects whole-brain perturbation programs into prospective functional coordinates; 85.48% of projected functional variance is assigned to context x perturbation interaction. This remains experiment nomination/falsification, not experimental causal closure.
