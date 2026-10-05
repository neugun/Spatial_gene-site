# Results at a glance — current browser authority (2026-10-04)

> **Start with the [browser-first scientific page](../index.html).** Current Figure 1–3 authority is [v18](current/CURRENT_FIG123_AUTHORITY_V18_20261004.md): **136/136 runtime checks PASS**, **14/14 audited datasets**, **91/91 capability checks**. Compact machine-readable sources are under [source_data/20261004](source_data/20261004/).

The table below is deliberately not a cross-dataset leaderboard. Every row is interpreted under its own molecular-evidence tier, biological holdout and matched comparator/null.

| Figure / question | System | Held unit | Current audited anchor | What it supports | What it does not support |
|---|---|---|---|---|---|
| **Fig1: which functional coordinate is readable?** | Bugeon VISp | 4 mice | genes rho 0.398; type 0.386; space 0.040 | gene/type-rich regime | genes universally dominate |
| Fig1 | Xu/PVH | 3 mice | genes rho 0.338 vs type 0.201 | continuous genes add rank information | every dynamic feature is calibrated |
| Fig1 | Shainer tectum | 6 fish | space rho 0.261 vs genes 0.032 | spatial/topographic regime | “gene model failed” as a biological conclusion |
| Fig1 | Zhao VISp | held animals | connectivity vs type: delta R2 ~+0.0011, delta Spearman ~+0.1269, delta AUROC ~+0.0660, delta AUPRC ~+0.1072 | useful rank/classification structure can coexist with tiny variance gain | one metric should replace the others |
| Fig1 boundary | MICRONS V1 EM | 13 acquisition sessions in one animal | RFx rho ~0.497, RFy ~0.536 from space; connectivity adds no broad robust increment after type+space | structural/retinotopic sanity and capability boundary | cross-animal transcriptomic inference |
| **Fig2: can missing molecular state be recovered?** | WARP | 3 fish | geometry rho 0.165→0.199; 3/3 improve; q=0.0053 | held-fish functional-geometry recovery | all zebrafish datasets should lift |
| Fig2 | McLachlan PERI/ECT | 3 held folds | rho 0.208→0.320; delta +0.112; matched-random q=0.0400 | exact anatomy matters | whole-brain reference is equivalent |
| Fig2 relative-only | C. elegans | 40 worms | delta rho +0.00236; q=0.0032 but absolute geometry remains <0 | significant relative recovery | positive absolute geometry |
| Fig2 boundary | Bugeon / Shainer / Condylis | held biological units | formal global geometry remains unsupported/directional | bridge/endpoint and geometry claims must stay separate | every successful map is a Fig2 win |
| **Fig3: what should a compact assay measure?** | M1 Patch-seq | n=10 grouped folds in frozen benchmark | K25 joint R2 0.369; 84.5% of full1000; +0.042 vs PERSIST; +0.123 vs random-p95 | primary fixed-budget benchmark; mean±SD and matched comparators shown | **exact upstream publication/accession not yet frozen; provenance gate remains open** |
| Fig3 | Allen Visual Learning 2P+HCR | 6 mice | K5 delta Spearman +0.0967; 6/6 wins; BH q=0.046875 | dataset-specific small-panel support | broad universal significance |
| Fig3 true hidden validation | Wang 2023 CEA EASI-FISH | 5 target-wise measured hidden-gene tests | selected hidden genes beat matched random for 5/5 targets | true measured hidden-gene replacement evidence; measured source Figshare 10.25378/janelia.21171373; reference GSE213828 | atlas nomination equals measurement |
| Fig3 compression | WARP | 3 fish | K25 AUROC 0.556 vs full31 0.555 | measured-panel compression can preserve function | strict hidden-gene replacement |
| Fig3 broad boundary | eight broad rows | dataset-specific held units | 0/8 BH q<0.05 | direction/effect remains reportable | universal selector claim |

## External new-data stress test: EASI-PASS 2026

The public JS078 golden matching output provides a **direct method-level check**: 1,199/1,439 accepted IoU matches (83.3%), 1,149/1,439 confident soma-print calls (79.85%), median z-restricted IoU 0.286.

More importantly, the public GAP intake behaves differently at three evidence levels:

1. matching/QC-only table → Fig1–3 biological stages **safe-stop**;
2. aligned molecular feature table → channels auto-canonicalize and exact-reference Fig2 can run;
3. multi-animal molecular+functional table → applicable Fig1–3 stages open under held animals.

The released EASI-PASS demo does **not** contain the full PBN/V1 functional biology data, so those biological findings are not labeled reproduced. See [external-readiness authority](current/EASIPASS_2026_EXTERNAL_READINESS.json) and [paper-native replication contract](current/EASIPASS_PAPER_NATIVE_REPLICATION_CONTRACT.json).

## Fig1–3 now compose into Fig4–7

Allen Visual Learning HCR reaches Fig7 through one public workflow, but later-stage gates remain scientifically active: Fig4 is a saturation boundary, the prespecified Allen-HCR→Zhao Fig6 target transfer is negative, and Fig7 proceeds only after a held-mouse operator gate. McLachlan, Condylis, MICRONS and Sorensen retain their own correct stop conditions.

See [end-to-end authority](current/FIG123_TO_FIG47_END2END_AUTHORITY_V2.md).

## How to audit a number

Use [PAGE_DATA_PROVENANCE_20261004.md](PAGE_DATA_PROVENANCE_20261004.md) to map a browser claim to its compact source file; use [DATA_CATALOG.md](DATA_CATALOG.md) for original public dataset accessions; use [ANALYSIS_CONTRACT.md](ANALYSIS_CONTRACT.md) for held-unit/null definitions.


## Figures 4-7 current extension

### Figure 4 - conditional functional completion
The strongest current multimetric examples are C. elegans pumping (delta R2 ~ +0.117; delta Spearman ~ +0.214; delta AUROC ~ +0.115), Sorensen ephys latent (delta R2 ~ +0.0777; delta Spearman ~ +0.451; delta AUROC ~ +0.273), Xu Ghrelin (delta R2 ~ +0.0639), Bugeon running/state, Allen VC2P reliability and WARP OMR-backward. Zhao native visual completion and Condylis saturation remain explicit boundaries.

### Figure 5 - response laws
Bugeon exact 24D improves direct Ridge R2 0.00927 -> 0.01852 with concordant rank/AUROC/AUPRC gains. Xu raw multitarget and Sorensen low-rank provide independent functional replications. A separate structural branch shows strong molecular readout of projection/connectivity in Chevee, C. elegans, MERGE-seq and Sorensen; shared-operator superiority is not required for that structural claim.

### Figure 6 - selective reuse
The main authority remains selective biologically aligned reuse: Bugeon sample efficiency, Zhao->Bugeon same-species transfer, mouse->fish, six-direction species transfer and mouse->human intrinsic ephys. A new prespecified Allen-HCR->Zhao visual-evoked transfer is negative on the target despite positive source OOF, and is therefore explicitly **not** promoted.

### Figure 7 - prospective intervention design
Whole-brain perturbation responses are strongly context dependent. Molecular context x perturbation interaction is 75.81%; after the frozen functional projection it is 85.48%. Allen HCR now exercises the full public chain through 44,133 perturbation-context conditions, but the output remains prospective ranking rather than causal evidence.
