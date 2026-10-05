# GAP/FSO - Molecular coordinates of in-vivo neural function

**GAP asks which reproducible components of a neuron's in-vivo response law are specified by molecular identity, and how that information should change the next functional-spatial-omics experiment.**

> **Start with the [browser-first scientific results page](index.html).** It is the current reader-facing entry point: biology question → public data source → molecular-evidence tier → held biological unit → matched baseline/null → quantitative result → allowed interpretation/boundary → source data/code.

**Current Figure 1–3 authority:** [v18, 2026-10-04](docs/current/CURRENT_FIG123_AUTHORITY_V18_20261004.md) • **136/136 runtime checks PASS** • **14/14 audited datasets** • **91/91 capability checks**. The 2026-10-04 Page refresh is built directly from compact authority/source-data files under [docs/source_data/20261004](docs/source_data/20261004/).

The project does not assume that genes reconstruct every moment of neural activity. Instead, it separates stable functional coordinates from context-dependent expression of those coordinates, evaluates them under held-biological-unit generalization, and asks when molecular state, anatomy, projection, behavior, or network context adds nonredundant information.

**Claim discipline:** boundary/negative results remain visible; Atlas/reference self-recovery is not relabeled as functional proof; compression is not relabeled as hidden-gene replacement; engineering fixtures are not relabeled as biological reproduction. See the [page provenance ledger](docs/PAGE_DATA_PROVENANCE_20261004.md) and [Page QA](docs/PAGE_QA_20261004.md).

Older requirement-aligned v3 graphics remain in [docs/FIGURES.md](docs/FIGURES.md) for provenance only. They are no longer the browser-facing authority after the 2026-10-04 refresh.

## Quick start

```bash
python -m venv .venv
# activate the environment, then:
pip install -e .
python -m pytest -q tests
```

Data are configured separately through `GAP_DATA_ROOT`; see [`docs/DATA_CATALOG.md`](docs/DATA_CATALOG.md) for the original public source of each core resource and [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for path setup.

## Repository layout

- `questions/`: biological entry points (Q1-Q6); start here for a paper claim.
- `gxa/`: reusable adapters, representations, operators, transfer and reporting code.
- `carma/`: paired Bugeon/Xu loaders and core same-cell response/molecular helpers.
- `scripts/`: canonical paper/reproduction entry points selected by the authority map.
- `docs/`: data provenance, baseline principles, inference contract and figure-to-code mapping.
- `tests/`: release smoke tests that do not require vendored raw data.


## Three biological principles

1. **Response law, not activity snapshot.** Molecular state changes more slowly than moment-to-moment activity. We therefore test response order, normalized response shape, selectivity, and state susceptibility before demanding full activity reconstruction.
2. **Molecular priors + biological context.** Space, projection, state, and history are candidate coordinates rather than nuisance variables. We ask what molecular information remains after the relevant context is represented.
3. **No universal gene-to-behavior law.** Molecular readability is a `region x function x coding architecture` property. Type-anchored, continuous-program, spatial, projection, state/history, and distributed regimes require different analyses.

## Six empirical questions organize the repository

| Question | Biological object | Main analysis | Repository section |
|---|---|---|---|
| **Q1. Which components of in-vivo function are molecularly specified?** | amplitude, order, shape, within-type position | held-unit target hierarchy; type vs continuous molecular state | `questions/q1_*` |
| **Q2. When can missing molecular state be recovered?** | unmeasured molecular dimensions | reference bridge validation; molecular completion; correspondence null | `questions/q2_*` |
| **Q3. What should be measured in a compact assay?** | genes, space, projection | marginal information gain; panel selection; gated spatial/projection experts | `questions/q3_*` |
| **Q4. Which state or behavior should be measured next?** | omitted functional dimensions | leave-one-behavior/state-out completion | `questions/q4_*` |
| **Q5. How does molecular identity parameterize a state-dependent response law?** | stable response geometry + susceptibility | nested response-law/operator model | `questions/q5_*` |
| **Q6. Are the recovered coordinates reusable?** | shared molecular-functional representation | held-animal low-label and ontology-matched transfer | `questions/q6_*` |

**Prospective Figure 7 is a separate synthetic causal-closure design**, not a seventh empirical claim. See `questions/q7_prospective_causal_closure/`.

## One analysis contract

Whenever the modalities permit, the same nested comparison is used:

`context/intercept -> hard identity -> + continuous molecular state -> + space/projection if independently supported -> + state interaction/operator`

The unit of inference is the **biological replicate** (animal/fish/worm), not the number of cells. Main analyses pair calibrated continuous metrics with rank/order or biologically meaningful classification metrics and a representation-matched null.

See [`docs/ANALYSIS_CONTRACT.md`](docs/ANALYSIS_CONTRACT.md) for the full contract.

## What counts as molecular evidence?

The repository keeps molecular resolution explicit rather than treating every identity label as transcriptomics:

- **same-cell measured expression:** molecular measurements and function come from the same neuron;
- **hard molecular identity:** type/driver/class labels without continuous same-cell expression;
- **external-reference molecular state:** continuous expression or posterior inferred from an independent atlas/reference;
- **class-linked external expression:** class-level lookup joined to identified neurons, not same-cell measurement.

Claims and nulls depend on the evidence tier. Atlas-derived expression is always labeled **inferred**, never measured.

## Baselines are matched to the question

The project includes simple linear specialists, supervised low-rank models, nonlinear multimodal models, representation learners, molecular foundation models, and gene-panel selectors. They solve different objects and are not collapsed into one leaderboard.

See [`docs/MODEL_AND_BASELINE_GUIDE.md`](docs/MODEL_AND_BASELINE_GUIDE.md) for the biological purpose of each comparator and the matched endpoint on which it is meaningful.

## Claim discipline

The repository keeps scientific boundaries visible instead of optimizing the story after seeing model scores. In particular:

- **Q5 is an operator/susceptibility question; Q6 is a reuse/sample-efficiency question.** A transfer model does not replace the state-dependent operator, and an operator result is not evidence of cross-system transfer.
- **Gamma is a susceptibility family, not a single universal gain.** Bugeon provides contrast-specific operator evidence; Xu/PVH independently supports a cross-state amplitude envelope/range across 11 conditions.
- **External models live at different abstraction levels.** Allen mechanistic V1, BrainBeacon, CEBRA and target-specific specialists are compared only on matched endpoints that they actually solve.
- **Boundary results remain public.** Negative held-unit R2, weak zero-shot transfer, unsupported oPhys benchmark status, and failed/null endpoints are retained because they define what the framework does not yet explain.

See [the claim/evidence ledger](docs/CLAIM_EVIDENCE_LEDGER.md).

## Data

GAP standardizes paired molecule/identity x neural-function resources across mouse, zebrafish, *C. elegans*, and additional systems. Public datasets remain at their original repositories/accessions; this code repository stores adapters, compact manifests, and reproducibility metadata rather than duplicating large raw datasets.

See [`docs/DATA_CATALOG.md`](docs/DATA_CATALOG.md) and [`docs/DATA_POLICY.md`](docs/DATA_POLICY.md).

## Reproducing a paper analysis

The public interface is figure/question first. Each question directory points to canonical scripts; reusable logic lives in `gxa/`. Internal historical sweeps are not the entry point.

See [`docs/FIGURE_TO_CODE_MAP.md`](docs/FIGURE_TO_CODE_MAP.md) and [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

For the current release boundary and consolidation status, see [docs/CODE_INVENTORY.md](docs/CODE_INVENTORY.md) and [docs/RELEASE_STATUS.md](docs/RELEASE_STATUS.md).
## Fast navigation

- [Figure requirements v3](docs/FIGURE_REQUIREMENTS_V3.md): current panel-by-panel design authority and completion gates.
- [Current v3 figures](docs/FIGURES.md): requirement-aligned Figure 1–7 graphics, with unresolved matched panels explicitly marked PENDING.
- [Latest integrated update](docs/LATEST_UPDATES_20260930.md): newest figure/result corrections and superseded values.
- [Cross-species extension](docs/CROSS_SPECIES_EXTENSION_20260930.md): latest mouse/fish/worm transport-space results kept separate from frozen Fig.6.
- [Model evolution / negative results](docs/MODEL_EVOLUTION_20260930.md): attention/flow/diffusion negatives, Gamma fusion status, and continual-learning replay policy.

- [Results at a glance](docs/RESULTS_AT_A_GLANCE.md): six biological questions, matched comparisons, representative held-unit results.
- [Baseline primer](docs/BASELINE_PRIMER.md): plain-language explanation of Ridge, RRR, PLS, CEBRA, NEUROeSTIMator, molecular foundation models, panel selectors, and the GAP response-law model.
- [Metric and inference guide](docs/METRIC_AND_INFERENCE_GUIDE.md): explains what R2, rank, AUROC, held-biological-unit splits, reliability, and matched nulls mean biologically.
- [Minimal examples](examples/): runnable synthetic examples that expose the split, residualization, and matched-null logic without large datasets.
- [Narrative crosswalk](docs/NARRATIVE_CROSSWALK.md): keeps the presentation, manuscript, and repository in the same biological order.
- [Current public authority](docs/CURRENT_PUBLIC_AUTHORITY.md): latest wording locks and superseded shorthand.
- [Claim/evidence ledger](docs/CLAIM_EVIDENCE_LEDGER.md): promoted, supportive, boundary, and pending claims.
- [External baseline boundaries](docs/EXTERNAL_BASELINE_BOUNDARIES.md): Allen mechanistic V1, BrainBeacon, and the modular molecular→function→context architecture boundary.
- [Model architecture by biological operation](docs/MODEL_ARCHITECTURE_BY_BIOLOGY.md): maps functional object, molecular encoder, Theta/Gamma operator, atlas bridge, and transfer heads to public code and biological claims.
