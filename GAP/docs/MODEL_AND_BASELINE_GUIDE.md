# Model and baseline guide

This document explains **what biological/statistical object each comparator solves**. It is intentionally not a single model leaderboard.

| Family | Core idea | Appropriate GAP role | What it does *not* establish |
|---|---|---|---|
| Ridge / ElasticNet | regularized supervised linear prediction | simple target-specific specialist and capacity floor | reusable functional geometry or mechanism |
| RRR | supervised low-rank multivariate regression | test whether a multidimensional phenotype is predictably low dimensional | nonlinear or state-conditioned dynamics |
| PLS | supervised latent directions maximizing cross-view covariance | compact supervised molecular↔functional coordinates | independent evidence of biological causality |
| sBNN | sparse nonlinear transcriptome→physiology predictor | nonlinear physiology specialist | a general in-vivo response-law representation |
| UnitedNet / coupled AE | coupled latent spaces across modalities | multimodal latent/prediction baseline | biological meaning of an arbitrary latent axis |
| NEUROeSTIMator | transcriptome→scalar activation/excitability score | molecular activation prior | multidimensional or state-conditioned response law |
| CEBRA | behavior/time-aligned representation learning | functional embedding baseline and optional late fusion | gene/receptor mechanism or molecular completion |
| BrainBeacon | large pretrained spatial-molecular encoder | external molecular prior; ask whether a compact pretrained coordinate improves the same held-cell response-law endpoint | same-cell functional operator or state susceptibility by itself |
| Allen / Ito / Arkhipov GLIF V1 | wiring + electrophysiology constrained mechanistic population simulation | population-dynamics/generative prior and mechanistic comparator | accurate held-cell molecular→function identity mapping |
| PERSIST / geneBasis / scGeneFit / Spapros | targeted feature/panel selection under different objectives | matched-K panel-design comparators | universal best panel for every functional objective |
| scGPT / Geneformer / UCE / scFoundation | large pretrained molecular representations | alternative molecular encoders under a frozen downstream GAP task | neural functional superiority by pretraining alone |
| GAP/FSO | hierarchical molecular/context representation → recoverable functional coordinates/operator | identify what is molecularly specified, what information is missing, and what transfers | universal superiority on every fully supervised scalar target |

## Comparator rule

A comparison is interpretable only when the models receive the same biological training split, predict the same endpoint, and are evaluated with the same metric family. A target-specific specialist is allowed to win a narrow scalar regression; that result defines the boundary of the broader GAP claim rather than being hidden.

## Abstraction-level rule

Population realism, molecular embedding quality, behavior-aligned activity geometry, and same-cell response-law inference are different objects. External models should therefore be used to locate which information layer they solve, not collapsed into a global rank. See `EXTERNAL_BASELINE_BOUNDARIES.md` for the matched Allen/BrainBeacon examples.
