# Model architecture by biological operation

The public architecture is organized by **what biological information is being represented**, not by the historical sequence of Transformer versions. Different figures activate different blocks; no figure should be read as “the biggest model wins.”

| Biological block | Input | Output | Public code | Scientific purpose |
|---|---|---|---|---|
| Functional-object definition | calcium-derived responses + labeled stimulus/state/task/event + biological-unit IDs | standardized response matrix, reliability, groups | `gxa/function_objects.py` | define the phenotype before molecular prediction |
| Function-aligned core | support context/activity from training data | low-D functional barcode `zF` and held-context prediction | `gxa/fso_function_core.py` | learn functional structure without allowing gene labels to shape it |
| Molecular/context encoder | measured molecular state + species/region/panel descriptors + optional reference memory | compact cell representation / operator coordinate | `gxa/factorized_gap_router_v21.py` and canonical scripts | map static identity/context to the frozen functional object |
| State-dependent operator | dynamic inputs `U(t)` + state/context proxy `m(t)` | stable coefficients `Theta` and susceptibility `Gamma` | `gxa/controller.py` | separate stable response organization from context-dependent modulation |
| Molecular completion | measured sparse panel + independent reference atlas | posterior/reference state + expected unmeasured molecular dimensions | `gxa/atlaslift.py` | add missing molecular information only after a valid bridge |
| Measurement / transfer heads | frozen coordinates + candidate gene/space/projection/state or source-domain prior | held-unit utility / low-label prediction | question-specific canonical scripts | ask what to measure next and what transfers |

## 1. Freeze the functional object before asking genes to explain it

`FunctionObjectBundle` standardizes paired cell-level responses, response names, biological groups and optional reliability. Raw movie pixels and image embeddings are excluded from this interface. The purpose is to prevent the molecular model from quietly changing the phenotype definition.

The first gate is therefore functional: a coordinate must reproduce across animals/fish/worms before its molecular readability is interpreted.

## 2. Learn functional structure without molecular labels

`FSOFunctionCore` infers a cell-specific functional barcode `zF` from support context/activity and predicts held-out activity under new contexts. Dataset-specific adapters absorb assay differences while the latent/operator logic is shared.

Genes, space and projection are intentionally absent from this stage. They are mapped to the frozen functional coordinate only later. This separation makes reverse decoding and matched controls interpretable instead of circular.
## 3. Encode molecular state and biological context at the resolution the task needs

The factorized GAP router uses a shared molecular encoder plus residual experts conditioned by species, region family, panel descriptors and functional-task semantics. The public interpretation is not “mixture of experts” as an engineering claim. It is a test of whether the same molecular program is sufficient everywhere or whether region/panel/task-specific residual structure is required.

The default information ladder remains:

`context/intercept → hard identity → continuous molecular state → optional anatomy/projection → state interaction/operator`

Every extra block must improve a held-biological-unit endpoint over the immediately simpler block and an appropriate matched null.

## 4. Separate stable response geometry from susceptibility

The controller implements

`A_i(t) = Theta_i^T U(t) + Gamma_i^T [phi(m(t)) * U_C(t)]`.

`Theta` describes stable response organization; `Gamma` describes how a state/context variable changes selected response components. The controller code **does not assume Gamma is molecularly predictable**. That is tested downstream.

This distinction matters because a state-dependent parameter can be well identified from activity but still lack the receptor/state/history information required for cross-animal molecular prediction.

Gamma should now be read as a **susceptibility family**, not one scalar coefficient. Bugeon supplies contrast-specific operator coordinates, whereas Xu/PVH shows that continuous molecular state predicts cross-state amplitude variability and total state range across 11 conditions beyond hard type. The promoted common object is the range or direction in which context can move a stable response law, not generic temporal complexity.

## 5. Add atlas information only as a separately validated bridge

AtlasLift maps cells to an external reference using only measured genes, retains the posterior, and derives expected expression of genes that were not assayed in vivo. Masked measured-gene recovery audits whether the bridge carries real molecular information before downstream function is scored.

Atlas-derived expression remains labeled inferred. It is never substituted for same-cell measurement in the evidence tier.

## 6. Latest external baselines refine the architecture

BrainBeacon shows that a generic pretrained molecular encoder can provide a useful compact prior for the dominant response-shape/Theta-like coordinate once dimensionality and held-animal splits are matched. The frozen V2.2 function-aligned core adds further information after sufficient target calibration. Neither solves log-gain/Gamma.

Allen/Ito/Arkhipov provides the orthogonal boundary: wiring + electrophysiology can reproduce population functional distributions while held-donor individual-cell OSI/DSI calibration remains weak.

The current evidence therefore motivates:

`general molecular encoder → function-aligned core → context/receptor susceptibility layer → held-unit functional operator`

The third block is a **missing-information hypothesis**, not a cosmetic architectural addition. Receptor, neuromodulator, state, history and projection variables should enter only if they close Gamma headroom in held biological units.

The latest receptor-aware screen partially closes this gap: under the exact 633-cell LOAO contract, type lookup is about R² 0.0732, while type + BrainBeacon4 + Atlas8 + receptor6 Ridge reaches about **0.11546 with 4/4 held animals positive**. This is evidence for a receptor/context susceptibility layer, not evidence that Gamma is completely solved or causally identified.
## 7. How to read historical model files

The repository retains several Transformer/router versions because exact manuscript results must remain traceable. These filenames are provenance, not the conceptual hierarchy.

For a new reader:

1. start with `function_objects.py` to understand the target;
2. read `fso_function_core.py` to see how functional structure is frozen independently of genes;
3. read `controller.py` for the Theta/Gamma decomposition;
4. use `atlaslift.py` only for the molecular-completion question;
5. inspect `factorized_gap_router_v21.py` when studying shared versus species/region/panel/task residual structure;
6. follow `FIGURE_TO_CODE_MAP.md` for the exact canonical script that produced a paper result.

## 8. What the model is expected to output

The scientific outputs are not arbitrary embeddings. Depending on the question they are:

- calibrated functional amplitude;
- relative response order / extremes;
- normalized response geometry or shape PCs;
- within-hard-type functional residual;
- stable operator coordinate `Theta`;
- context/state susceptibility `Gamma`;
- marginal information gain from an additional molecular/anatomical/functional measurement;
- target-label efficiency or ontology-matched transfer.

A new architecture is useful only when it improves one of these prespecified objects under the corresponding held-unit split and null/control.

## 9. Keep operator identification separate from transfer

The same representation block may be used in both analyses, but the biological endpoints differ. Q5/Figure 5 asks whether a stable response-law coordinate or susceptibility is molecularly specified. Q6/Figure 6 asks whether a recovered coordinate reduces paired-label demand or transfers under an explicitly aligned ontology.

A model can pass one question and fail the other. Public text should therefore report the question-specific endpoint rather than assigning one global architecture rank.

## 10. Architecture changes are gated by frozen ablations

Generic attention, flow, diffusion and a distribution-preserving residual have been tested and do not replace the current stack. Replay + L2-SP is currently preferred for continual adaptation because it preserves the shared prior while learning a new region.

See [MODEL_EVOLUTION_20260930.md](MODEL_EVOLUTION_20260930.md) for the negative-result and continual-learning ledger.
