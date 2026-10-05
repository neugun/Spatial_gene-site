# Start here: the biology before the code

> **2026-10-04 current entry point:** [browser-first scientific page](../index.html) · [v18 Figure 1–3 authority](current/CURRENT_FIG123_AUTHORITY_V18_20261004.md) · [page provenance](PAGE_DATA_PROVENANCE_20261004.md).
> Current runtime: **136/136 PASS**. The older v3 figure layouts remain available for historical provenance but are not the current numerical authority.

GAP/FSO starts from one biological question: **which parts of a neuron's in-vivo functional response are stable enough across biological replicates to be specified by molecular identity?**

The repository is organized around that question rather than around model classes or software history. A reader should be able to move from biological claim → dataset → matched baseline → held-biological-unit test → canonical script without learning internal project names first.

## The three principles

**1. Predict a response law, not every activity sample.** Molecular state is comparatively slow; activity is fast and context dependent. GAP therefore first asks whether molecular identity predicts response order, normalized response shape, within-type functional position, or state susceptibility.

**2. Treat context as biology when it carries independent information.** Space, projection, state, stimulus, history, and network context are added only when they improve a held-out biological endpoint under a matched control.

**3. Expect multiple coding regimes.** Molecular readability is a region × function × coding-architecture property. Some functions are type-anchored, some depend on continuous molecular programs, some on spatial/projection coordinates, and some are distributed or history dependent.

## Five guardrails from the current audit

1. **Association is not held-biological-unit prediction.** The external oPhys/Xenium notebooks are a useful stress test: many molecular associations can coexist with poor leave-one-mouse calibrated prediction of a scalar tuning metric.
2. **Conserved versus animal-specific is the operator gate.** The earlier shorthand “shape stable, gain unstable” is too coarse; both Theta- and Gamma-like axes can be molecularly readable when the axis itself is conserved across animals.
3. **Gamma is a susceptibility family.** Bugeon tests a state-contrast operator; Xu/PVH adds a cross-state amplitude envelope/range across 11 conditions. These are related susceptibility objects, not evidence that every dynamic coefficient is molecularly specified.
4. **Architecture blocks are biological operations.** Generic molecular representation, function-aligned structure, context/receptor susceptibility, atlas completion and transfer are separate information layers. A larger model is not automatically a better biological explanation.
5. **Negative and mixed results stay visible.** A result that fails calibration, a null, a transfer gate or a bridge-validity test remains part of the public evidence boundary rather than disappearing from the story.

For the exact current wording, read `CURRENT_PUBLIC_AUTHORITY.md`; for claim status, read `CLAIM_EVIDENCE_LEDGER.md`.

## The six questions

The six directories under `questions/` are the public entry points. They cover molecularly specified function, molecular completion, compact measurement design, functional/state completion, state-dependent response laws, and transfer to new animals or systems.

If you are reading the paper, follow Q1→Q6 in order. For the browser-first result view, open `FIGURES.md`. Prospective Figure 7 is synthetic/causal-design only and is intentionally separated from the six empirical questions. If you are planning an experiment, begin with Q1 and then jump to the missing-information question that applies to your dataset.
## What every main analysis must show

A main result should specify: the biological unit held out (mouse/fish/worm/donor/session), the functional target, the molecular-evidence tier, the matched baseline, the null or correspondence control, and at least one metric that answers the stated biological question.

This means that a higher model score is not itself the endpoint. For example, calibrated R² asks whether amplitude transfers; Spearman asks whether functional order transfers; AUROC asks whether biologically meaningful response groups can be separated. These can agree or dissociate, and the dissociation is itself informative.

## Where the code lives

- `gxa/` contains reusable data adapters and model/analysis components.
- `questions/` explains each scientific question in biological language and points to canonical scripts.
- `scripts/` contains canonical figure/reproduction scripts selected from the larger development workspace.
- `docs/DATA_CATALOG.md` lists the data sources and their roles.
- `docs/MODEL_AND_BASELINE_GUIDE.md` explains what each comparator actually tests.
- `docs/FIGURE_TO_CODE_MAP.md` maps manuscript figures to canonical code.
- `docs/ANALYSIS_CONTRACT.md` defines the held-unit and null-control rules.

Large datasets and internal experiment outputs are deliberately not vendored into the repository. Paths are configured through environment variables described in `docs/DATA_POLICY.md` and `gxa/paths.py`.

## Reading order for reviewers and collaborators

Read `README.md` → this page → `docs/CURRENT_PUBLIC_AUTHORITY.md` → `docs/CLAIM_EVIDENCE_LEDGER.md` → the relevant `questions/q*/README.md` → `docs/ANALYSIS_CONTRACT.md` → `docs/FIGURE_TO_CODE_MAP.md`. Only then is it useful to inspect individual model classes. This keeps the biological question ahead of engineering implementation.


## Continue beyond Figure 3

The current Page now follows the same audited contract through Figures 4-7:

- **Figure 4:** biology-aware functional completion; same-axis technical derivatives are collapsed before estimating molecular headroom.
- **Figure 5:** identical X/Y/folds for direct Ridge versus shared full-rank / reduced-rank operators; functional and structural endpoints remain distinct.
- **Figure 6:** source-only route selection with target evaluation required for a final transfer claim when target labels exist.
- **Figure 7:** only a held-unit-qualified molecular->function operator may enter prospective perturbation ranking.

Compact public authorities: `docs/current/FIG4_7_END2END_TOOL_V2_AUTHORITY.md`, `docs/current/FIG6_MASTER_AUTHORITY_V51.md`, and `docs/current/FIG7_PERTURBGAP_AUTHORITY.md`.
