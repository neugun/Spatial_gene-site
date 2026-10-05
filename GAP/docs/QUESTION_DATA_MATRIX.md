# Question-centric data matrix

The same dataset can answer several biological questions, but only when it contains the required modality and biological replication. This matrix complements `DATA_CATALOG.md` by organizing resources around Q1-Q6 rather than by publication.

| Resource | Public source | Molecular evidence tier | Functional object | Primary GAP questions |
|---|---|---|---|---|
| **Bugeon 2022, mouse VISp** | Figshare DOI `10.6084/m9.figshare.19448531` | same-cell 72-gene coppaFISH + hierarchy | visual response, running/state modulation, repeated response kernels | Q1 within-type/shape; Q2 atlas completion; Q3 panel design; Q5 operator; Q6 sample efficiency |
| **Xu 2020, mouse PVH** | Science `10.1126/science.abb2494`; GEO `GSE148568`; CaRMA code `github.com/neugun/CaRMA-imaging` | same-cell multiplex FISH markers | survival/internal-state responses including Ghrelin | Q1 molecular continuum; Q2 headroom/reference boundary; Q4 state complementarity; **Q5 cross-state flexibility envelope** |
| **WARP zebrafish** | Figshare dataset `29962931` | same-cell gene panels + registered 3-D anatomy | brain-wide visual/behavior response battery | Q1 within-type; Q3 spatial/receptor-panel design; Q4 omitted behavior; Q6 aligned fish transfer |
| **Shainer 2025 zebrafish** | Zenodo `10.5281/zenodo.14146655` | HCR marker identity + spatial/morphological context | optic-tectum visual response/topography | Q1/Q3 coding-regime boundary; Q4 response geometry; Q6 WARP→Shainer looming transfer |
| **Zhao 2026 mouse VISp** | `github.com/WKLabION/Trimodal-Data-and-Analysis` | activity + morphology/projection + 2cEASI-FISH in matched neurons | visual activity and projection phenotype | Q1 multimodal complementarity; Q2 phenotype-specific lift; Q3 projection measurement |
| **Condylis 2022 mouse S1** | G-Node GIN `10.12751/g-node.7q0lz0` | same-cell molecular barcode | tactile task/passive response | Q1 boundary/coding regime |
| **Atanas 2023 + CeNGEN/connectome** | Cell `10.1016/j.cell.2023.07.035`; `github.com/flavell-lab/AtanasKim-Cell2023`; Zenodo `8185377` | identified neuron/class linked to external class transcriptome + connectome | multi-behavior whole-nervous-system tuning | Q1 class-linked systems validation; Q4 omitted behavior/manifold; reverse molecular geometry |
| **Allen Visual Coding 2P** | Allen Brain Observatory Visual Coding release/API | Cre identity + external exact-CCF/MERFISH continuous state | visual tuning/reliability and running modulation | Q1 coarse vs continuous external state; Q2 independent visual-system replication |
| **Allen Visual Behavior 2P** | Allen Brain Observatory Visual Behavior release | Cre identity + external atlas state | task/change/omission and state susceptibility | Q2 calibration-vs-order dissociation; Q5 independent state-operator test |
| **OpenScope GLO** | Allen/OpenScope GLO release | optotagged/coarse identity + ephys features | electrophysiological state/operator response | Q5 independent dynamic validation |
| **DANDI 001532 / oPhys Xenium** | public DANDI + `oPhys_Transcriptomics_Analysis` repository | same-cell 299-gene Xenium with 2P visual physiology | tuning, state/context, temporal response features | external method stress test; **not yet promoted as a GAP benchmark** |

## Why evidence tier is shown next to the dataset

Same-cell expression, class-linked transcriptomics, Cre identity, projection identity, and atlas-inferred molecular state are not interchangeable. The same numerical performance supports different claims depending on what was physically measured in the functional cells.

## Dataset eligibility by question

Q1 requires a function target and an identity/molecular variable plus biological-unit replication. Q2 additionally requires a defensible reference bridge and molecular headroom. Q3 requires a candidate measurement block and a matched way to destroy only that block's correspondence. Q4 requires multiple functional/state dimensions measured in the same system. Q5 requires repeat/state structure sufficient to identify a reliable response-law coordinate. Q6 requires an explicitly aligned source-target biological variable and a target-label budget curve.

## External oPhys gate

The frozen oPhys repository is useful now as a methodological reference because it makes the hierarchy from transcriptomic identity to tuning/state/context analyses explicit. Its published notebook outputs are not imported as GAP benchmark evidence. DANDI 001532 enters the matched benchmark only after persistent assets, cell-level joins, native target reproduction, and identical held-mouse baseline/GAP evaluation have passed the intake gate.
