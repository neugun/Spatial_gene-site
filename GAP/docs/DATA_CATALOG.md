# Data catalog

GAP does not redistribute large public datasets. Dataset adapters standardize molecular/identity variables, biological-unit identifiers, functional targets, spatial/projection metadata, and provenance while preserving the original source's molecular-resolution tier.

The table below answers four questions for each core resource: **what biological system is measured, what functional object is available, what molecular evidence is available, and where the original public data/code come from.** The repository stores adapters and compact provenance, not a second copy of the raw archive.

## Core paired / identity-function resources

| Dataset | System and functional object | Molecular / identity evidence | Public source / accession | GAP role |
|---|---|---|---|---|
| Bugeon 2022 | mouse VISp; visual responses and running/state modulation | same-cell 72-gene coppaFISH + hierarchy | Figshare `10.6084/m9.figshare.19448531` | continuous molecular state, within-type signal, operator, AtlasLift positive regime |
| Xu 2020 / CaRMA | mouse PVH; responses across 11 survival-related states | same-cell multiplex FISH markers | Science `10.1126/science.abb2494`; GEO `GSE148568`; CaRMA code/example data: `github.com/neugun/CaRMA-imaging` (also maintained/forked by Saintgene-Xu-lab) | state susceptibility, within-type molecular information, headroom boundary |
| WARP | larval zebrafish; brain-wide visual/behavior responses | same-cell gene panels + registered anatomy | Figshare dataset `29962931`, author `postprocessed/` arrays | within-type, spatial expert, BehLift, transfer |
| Shainer 2025 | zebrafish optic tectum; visual response/topography | HCR marker identity + spatial/morphological context | Zenodo `10.5281/zenodo.14146655` | topographic boundary / spatial coding regime |
| Zhao 2026 | mouse VISp; visual activity + morphology/projection | same-neuron activity, morphology/projection and 2cEASI-FISH | `github.com/WKLabION/Trimodal-Data-and-Analysis` | multimodal complementarity / projection endpoint |
| Condylis 2022 | mouse S1 L2/3; tactile task/passive responses | same-cell molecular barcode | G-Node GIN `10.12751/g-node.7q0lz0` | task-selectivity boundary under current extracted panel |
| Atanas 2023 + CeNGEN + Cook connectome | *C. elegans*; brain-wide behavior tuning | identified-neuron/class join to external class transcriptomes + connectome | Cell `10.1016/j.cell.2023.07.035`; code `github.com/flavell-lab/AtanasKim-Cell2023`; processed-data Zenodo record `8185377` | class-linked molecular + connectome systems validation |
| Allen Visual Coding 2P | mouse visual cortex; standardized visual/state features | Cre identity + external exact-CCF/MERFISH state | Allen Brain Observatory Visual Coding API / `ApiCamCellMetric` | external continuous molecular-state test |
| Allen Visual Behavior 2P | mouse visual cortex; task/change/omission and state susceptibility | Cre identity + external atlas state | Allen Brain Observatory Visual Behavior release | calibration-versus-order dissociation |
## Evidence tiers are part of the data model

A row with same-cell measured expression can support different claims from a row with only a Cre line, transcriptomic class lookup, or atlas-inferred state. The adapters therefore retain this distinction explicitly. A class-linked or atlas-inferred molecular feature is never silently relabeled as direct same-cell transcriptomics.

## Additional resources

Additional adapters cover projection-rich datasets, activity tagging, Patch-seq/electrophysiology, structural/connectomic controls, spatial references, and cross-species transfer resources. They enter the paper only when they answer a prespecified biological question and have an appropriate biological-unit split and null/control.

Before a public tag, each question-specific reproduction note should freeze the exact upstream file/version required from these sources so that the repository can be reproduced without relying on a lab-local directory layout.

## 2026-10-04 tool-generalization and external stress resources

| Dataset | System / role | Public source | Current GAP interpretation |
|---|---|---|---|
| McLachlan et al. 2025 | mouse parahippocampal / PERI-ECT; independent Fig1–3 integration test | Cell Reports 44(1):115175; DOI 10.1016/j.celrep.2024.115175 | exact PERI/ECT reference; Fig2 geometry supported; Fig3 directional only; later operator gate can stop |
| M1 Patch-seq fixed-budget benchmark | mouse M1 Patch-seq; measured transcriptome → 16 electrophysiology outputs; n=10 grouped folds in the frozen benchmark | **PROVENANCE PENDING:** exact upstream publication/accession has not yet been frozen; current browser-facing compact extraction is `docs/source_data/20261004/M1_FIXED_BUDGET_MULTIMETRIC.csv` | Fig3 primary fixed-budget benchmark may be inspected, but public-release provenance remains open until the exact upstream file/version is frozen |
| Wang et al. 2023 CEA EASI-FISH | mouse central amygdala; 29 measured same-cell EASI-FISH genes + five retrograde projection channels (BNST, SNlat, vlPAG, PBN, PCRt) | measured table: Janelia Figshare DOI `10.25378/janelia.21171373`; independent CEA scRNA reference: GEO `GSE213828` | Fig3 true measured hidden-gene replacement; hidden genes are evaluated because they are actually measured, not merely atlas-nominated |
| MICrONS minnie65 | mouse V1 EM + functional coreg; structural boundary benchmark | MICrONS public CAVE / static repositories (minnie65_public) | one animal; session/acquisition-held only; no transcriptome; sparse connectivity evaluated only in the measured 1,151-neuron cohort |
| Allen Visual Learning 2P + HCR same-cell extraction | mouse visual learning; 6 mice; 21 shared HCR genes; same-cell functional features | exact source/session contract is frozen in the public extraction/source-data authority rather than inferred from a generic Allen release label | 1,321 outcome-blind deduplicated same-cell rows; hit-sparse endpoints excluded from cross-mouse primary benchmark |
| EASI-PASS 2026 JS078 public demo | method-level cross-modal matching stress test; one 2P plane + 5-channel HCR | github.com/orena1/easi-pass; preprint DOI 10.64898/2026.08.21.746328 | direct matching benchmark can be reproduced; full PBN/V1 functional biology is not present in the demo and is not claimed reproduced |
| Pauli/Chen PBN reference | mouse parabrachial nucleus scRNA-seq reference | eLife 2022, DOI 10.7554/eLife.81868; GEO GSE207708 | exact-region PBN reference: 21 neuronal subclusters × 14,437 genes; reference-headroom evidence only |

### Source-status rule

A dataset name is not treated as sufficient provenance. For each promoted result, the Page reports the molecular-evidence tier, biological unit, exact public accession/DOI when available, and the compact GAP extraction/authority file used for analysis. If an exact upstream public file/version has not been frozen, that limitation stays visible rather than being replaced by a generic repository label.