# Biological question / figure -> canonical code map

This is the public-facing map. The **scientific panel contract** is [`FIGURE_REQUIREMENTS_V3.md`](FIGURE_REQUIREMENTS_V3.md); the **current rendered figures and numerical reading notes** are in [`FIGURES.md`](FIGURES.md). The browser figures themselves are reproduced by [`scripts/render_requirement_aligned_figures.py`](../scripts/render_requirement_aligned_figures.py).

| Question / figure | Scientific role | Canonical entry points |
|---|---|---|
| Q1 / Fig. 1 | identity → within-identity molecular state → position/context; held-unit target hierarchy | `scripts/114_target_formulation_scorecard.py`; `scripts/115_type_gene_function_decomposition.py`; `scripts/102_bugeon_lowdim_state_controller.py`; `scripts/106_state_controller_authority_v2.py` |
| Q2 / Fig. 2 | reference bridge, within-type AtlasLift, correspondence null, named-gene recovery, optional anatomy/space | `scripts/315_exact_ccf_prior_atlaslift.py`; `scripts/506_unified_visual_state_benchmark.py`; `scripts/507_summarize_unified_visual_metrics.py`; `scripts/509_atlaslift_model_gain_multimetric.py`; `scripts/513_aggregate_core9_atlaslift_multimetric.py`; `scripts/516_directlift_sensitivity.py` |
| Q3 / Fig. 3 | gene-specific nomination, exact-anatomy sensitivity, fixed-budget panel design, receptor/effector panel, prospective assay | `scripts/110_panel_pareto_unified.py` plus the frozen AtlasLift single-gene and panel-selector entry points documented in the question README |
| Q4 / Fig. 4 | BehLift, conditional molecular information, held-unit functional manifold and gene-family→geometry analysis | `scripts/108_behlift.py`; `scripts/110_panel_pareto_unified.py` plus the manifold/FNP analyses listed in the manuscript authority |
| Q5 / Fig. 5 | low-dimensional shared operator, conservation gate, Gamma susceptibility, reverse molecular-functional geometry | `scripts/102_bugeon_lowdim_state_controller.py`; `scripts/103_bugeon_scalar_controller_validation.py`; `scripts/104_bugeon_scalar_controller_compact.py`; `scripts/105_bugeon_ranked_controller.py`; `scripts/613_xu_functional_flexibility_20260929.py`; `scripts/614_xu_functional_flexibility_null_20260929.py` |
| Q6 / Fig. 6 | strict zero-shot evolutionary transfer, matched representations, mechanism/ablation, few-shot and primate geometry | `scripts/520_foundation_v20_arch_sweep.py`; `scripts/521_eval_foundation_v20_sweep.py`; `scripts/522_foundation_v21_parallel.py`; `scripts/523_eval_foundation_v21_parallel.py`; `scripts/526_cv_gap_v21_factorized.py`; `scripts/527_cv_v22_expanded_compact.py`; `scripts/531_cv_gap_factor_compact_hybrid.py`; `scripts/535_cv_gap_crossspecies_refexpert.py`; see `docs/CROSS_SPECIES_EXTENSION_20260930.md` |
| Prospective Fig. 7 | frozen observational model → perturbation nomination → coordinate-specific causal falsification | `questions/q7_prospective_causal_closure/README.md`; current schematic `docs/figures/v3/Fig7_REQUIREMENT_V3.png` |

## Public-interface rule

Numeric script IDs are provenance identifiers, not the scientific story. Readers should enter through the biological question pages and figure guide; scripts remain stable so each manuscript result can be traced to an implementation.

The figure renderer is intentionally separate from analysis scripts: it reads the **frozen public numerical authority** and produces the browser figures. A new analysis result should first change the relevant result/claim authority before it changes a plotted panel.
