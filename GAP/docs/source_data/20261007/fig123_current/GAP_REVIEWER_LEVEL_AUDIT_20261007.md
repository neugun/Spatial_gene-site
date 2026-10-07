# GAP reviewer-level manuscript audit — 2026-10-07

## Executive judgment

The strongest version of the GAP paper is **not** a model-zoo paper, a gene-panel paper, or a universal gene→activity paper. Its distinctive scientific claim is:

> Stable biological identity does not specify all neural activity. It specifies **particular functional coordinates and response-law parameters**, at a resolution that differs across circuits. Those coordinates can then be used to decide what information is missing, what should be measured, and when transfer should or should not be attempted.

This gives the manuscript one continuous causal/logical ladder:

**readable coordinate → missing molecular information → compact measurement → missing functional information → response law → transfer/routing → falsification experiment**

The POCO paper is useful because it demonstrates how to make one architectural idea carry an entire paper without turning the Results into an inventory. GAP should copy that discipline.

## Recommended title

### Preferred
**Molecular and anatomical identity constrain transferable coordinates of neuronal function**

Why this is safer and stronger than “Molecular identity specifies…”:
- Fig.1 already shows that space can dominate sparse molecular identity in Shainer.
- structural/projection information becomes important in Fig.5.
- “constrain” avoids implying deterministic reconstruction of moment-to-moment activity.
- “transferable coordinates” still points directly to the Fig.5–6 conceptual payoff.

### Stronger alternate if Fig.5/6 are further closed
**GAP reveals biologically anchored response laws across neural states and species**

Use this only if the response-law and transfer sections become the unquestioned center of the final figure set.

## The paper's real novelty

### Novelty 1 — target formulation before model choice
Most multimodal work asks whether one modality predicts another. GAP first asks **what functional object is actually stable/readable**: exact amplitude, rank/extreme state, normalized response geometry, state susceptibility, or projection/structural coordinate.

This is scientifically stronger than simply improving R2.

### Novelty 2 — typed information rather than “multimodal integration”
GAP keeps measured transcript, hard identity, atlas-inferred expression, space, projection and function semantically distinct. That prevents several common overclaims:
- Cre line ≠ measured gene expression.
- CeNGEN class lookup ≠ same-cell transcriptomics.
- atlas inference ≠ observation.
- projection/connectivity ≠ behavior.
- response gain ≠ Gamma.

### Novelty 3 — missing-variable logic
Figures 2 and 4 are dual:
- Fig.2 asks whether missing **molecular** information is worth adding.
- Fig.4 asks whether missing **functional** information is worth measuring.

The conceptual variable is **headroom/nonredundancy**, not imputation accuracy.

### Novelty 4 — response law rather than static phenotype
Fig.5 is the conceptual center. The paper becomes more than a benchmark only if the reader understands:
- theta = stable response geometry / response-law coordinate;
- Gamma = state/context susceptibility;
- shared law = reusable structure across targets/conditions;
- direct predictability and shared-law gain are distinct.

### Novelty 5 — transfer as routing/fallback
Fig.6 is strongest when it explicitly rejects naïve universality. GAP does not say “one representation works everywhere.” It asks whether a source-supported route should be shared, adapted, or rejected.

POCO's negative multispecies pooling result is an excellent external conceptual precedent for this logic.

### Novelty 6 — end in a falsifiable experiment
Fig.7 is valuable only if it is framed as **prospective experimental design**, not causal simulation. This makes the paper end with a concrete way to prove GAP wrong.

---

# Main reviewer vulnerabilities and required fixes

## Risk 1 — “You have many modules but no single method”
This is the largest narrative risk.

### Current problem
AtlasLift, FSO-Select, BehLift/ManifoldLift, Operator, Transport and Perturb-GAP can look like separate projects.

### Fix
Present GAP as one **typed decision framework**:
1. identify a stable coordinate;
2. estimate what information is missing;
3. decide whether to measure, infer or transfer it;
4. stop when a scientific gate fails.

Every module should be introduced as an operation on that same object.

The paper schematic should therefore show one DAG, not seven boxes of equal status.

## Risk 2 — Fig.1 could be misread as “genes beat types”
That is not the current evidence.

### Current authority
- Xu: type+genes > type directionally in 3/3 mice, but n=3 gives an exact P floor of .125.
- Bugeon: strong absolute gene readability but only a small genes-over-type increment, P=.1875.
- Shainer: space is the strongest current held-unit result.

### Fix
Headline:
**The biological resolution required to read out function is circuit dependent.**

Do not headline:
**Continuous genes outperform cell types.**

## Risk 3 — Fig.2 can be attacked as imputation dressed up as function
The corrected Xu result is actually the defense.

### Fix
Make the negative Xu PVH result central:
- external reference recovers measured molecular structure;
- downstream Ghrelin lift is not robust.

This proves that GAP distinguishes “good molecular mapping” from “useful functional information.”

The negative is therefore mechanistic evidence for the liftability rule.

## Risk 4 — Fig.3 is close to PERSIST/sBNN prior art
At K25, FSO-functional is tied with sBNN.

### Fix
Do not claim a novel generic supervised selector.
The novelty is the **experimental-design contract**:
- in-vivo functional object;
- fixed gene budget;
- held biological unit;
- measured compression versus true hidden replacement kept separate;
- independent same-cell HCR replication;
- explicit objective-specific Pareto tradeoff.

PERSIST and sBNN should be treated as serious baselines, not straw men.

## Risk 5 — Fig.4 can become a vague “BehLift” story
The useful scientific result is target dependence.

### Fix
Build the main figure around:
**headroom × orthogonality → molecular gain**

The most important panel should be an across-target scatter/phase diagram, not a list of positive datasets.

## Risk 6 — Fig.5 absolute R2 can look small
The exact24D value is intentionally hard and should not be the headline.

### Fix
Headline the six-task shared-law result and parameter efficiency.
Use exact24D as a stress test/inset.

Make the comparator structure explicit:
same X, same Y, same folds, direct predictor versus structured/shared predictor.

## Risk 7 — structural and functional operators can be conflated
Chevee is strongly structurally predictable but has little extra shared-law gain.

### Fix
Use two questions:
A. Is structure/function predictable?
B. Does a **shared law** improve over matched independent mapping?

This prevents a reviewer from correctly pointing out that direct predictability alone does not prove shared operator structure.

## Risk 8 — Fig.6 can look like cherry-picked transfer
The defense must be built into the main figure.

### Fix
Show:
- source-only routing;
- negative/reverse routes;
- fallback;
- biological n separate from task count;
- source-donor uncertainty for human.

Mouse→fish: 15 task×fish contrasts = 3 biological fish.
Human: target-donor CI positive, strict source-donor CI crosses zero.

This should be visible, not buried in Methods.

## Risk 9 — Zhao→Bugeon is still post-freeze qualified
Do not silently treat it as fully closed conservation evidence.

### Fix
Keep the frozen v5.1 result but visually mark:
**post-freeze transfer-null qualification pending**

If the target-only/ortholog-permutation audit closes positively, remove the qualifier through a new authority rather than silently editing the frozen package.

## Risk 10 — Fig.7 provenance and causal language
The Oct1 Perturb-GAP authority uses a frozen 8-coordinate functional map. The exact upstream compute generator is not yet fully recovered, and that frozen map should not automatically be called identical to the current Fig.5 main operator.

### Fix
In the current manuscript:
- call it a **frozen GAP functional map** or **prospective functional projection**;
- do not say “the Fig.5 operator predicts…” unless the artifact/hash correspondence is explicitly recovered;
- retain Fig.7 as prospective experiment design.

Until the lineage is fully closed, Fig.7 should be conceptually strong but numerically conservative.

## Risk 11 — POCO could distract from the main paper
POCO is valuable, but there are no GAP-conditioned-POCO results yet.

### Fix
For now:
- Introduction: related dynamic forecasting problem.
- Discussion: complementary distinction between learned unit identity and biological identity.
- Fig.6/ED: conceptual negative-pooling precedent.
- Do **not** create a main Results subsection around POCO until E3/E4 passes held-animal tests.

## Risk 12 — “stable” may be stronger than what every dataset measures
Held-animal readability is not always the same as longitudinal within-cell stability.

### Fix
Prefer:
- “held-unit generalizable”
- “reproducible functional coordinate”
- “biologically readable coordinate”

Reserve “stable” for datasets/targets with direct temporal/session stability evidence.

---

# What should be in the main paper versus ED

## Main Fig.1
Xu, Bugeon, Shainer + one WARP within-type panel.
Use Zhao as a boundary, not a positive anchor.

## Main Fig.2
WARP + McLachlan positive.
Xu PVH negative decoupling is essential.
Bugeon exact VISp can be compact boundary/ED.

## Main Fig.3
M1 K curves + K25 paired result.
Allen HCR K5 replication.
Wang hidden replacement.
Bugeon objective Pareto.
Move large selector inventories to ED.

## Main Fig.4
Bugeon/Xu/WARP/Sorensen selected completion.
Condylis saturation.
Headroom/orthogonality summary should be the mechanistic centerpiece.

## Main Fig.5
Bugeon six-task shared law.
Xu and Sorensen functional replication.
One compact structural-law contrast.
Exact24D as stress-test inset.
Full molecular-program zoo in ED.

## Main Fig.6
Sample efficiency.
Mouse→fish.
Six-direction route/fallback.
Human with both CIs.
Zhao→Bugeon only if qualification is visible.
Move broad foundation/model zoo to ED.

## Main Fig.7
Interaction decomposition.
Missing-edge negative.
2–3 falsification candidates.
One prospective experiment card.
Everything else ED.

---

# Manuscript-specific sentence rules

Use:
- “reads out”, “constrains”, “is informative about”, “improves held-unit prediction”, “supports a shared law”, “nominates”.
Avoid unless directly demonstrated:
- “determines”, “encodes” for transcriptome→function causal language;
- “generalizes” without naming the held unit;
- “causal” in Fig.7;
- “foundation model” as the paper identity;
- “significant” when the inferential unit is a conditional shuffle rather than biological-unit exchangeability.

For every important result, write the sentence in this order:
**biological question → exact comparator → held biological unit → effect → inference scope → meaning.**

---

# Recommended abstract logic

1. Problem: stable identity versus dynamic activity.
2. Observation: biological identity reads specific functional coordinates, with circuit-dependent resolution.
3. Missing-molecule result: reference recovery and functional lift decouple.
4. Measurement result: compact panels preserve function under fixed budgets.
5. Mechanism: response-law sharing separates coordinate from susceptibility.
6. Generalization: source-validated routing transfers selectively and exposes boundaries.
7. End: GAP is a calibrated experimental-design framework, not universal activity reconstruction.

---

# Critical experiments that would most increase paper strength now

1. **Close Zhao→Bugeon transfer-null qualification.**
2. **Recover and hash-lock Fig.7 upstream compute lineage.**
3. **Run POCO-strict on Bugeon first, then test molecular/GAP initialization of unit embeddings.**
   The most valuable endpoint is zero-shot/adaptation speed, not asymptotic forecast accuracy.
4. **Make Fig.4 headroom×orthogonality summary a real cross-target analysis if not already frozen.**
5. **Convert Fig.5 six-task/shared-law result into the dominant visual rather than letting exact24D define reader perception.**
6. If feasible, obtain one truly prospective new measurement validating a Fig.3 panel or Fig.7 perturbation; this would materially change the paper's impact ceiling.

## Current impact ceiling

Without new wet-lab validation, the manuscript is strongest as a rigorous multimodal computational/experimental-design framework with unusually explicit boundaries.

With one prospective validation — compact panel or perturbation — the paper can credibly claim that GAP not only organizes existing multimodal data but **improves what experiment is performed next**.
