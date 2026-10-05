# Reproducibility

## Public interface

Start from a biological question in `questions/`, then follow the linked canonical entry point. Reusable functions live in `gxa/`; numeric script IDs are stable provenance identifiers, not the intended narrative entry point.

## Data paths

Large external datasets are not stored in the repository. Configure:

```bash
export GAP_DATA_ROOT=/path/to/gap-data
export GAP_RESULTS_ROOT=/path/to/gap-results
export GAP_CACHE_ROOT=/path/to/gap-cache
```

The public package uses `gxa.paths` to resolve these locations. Dataset-specific preparation instructions will specify the expected subdirectory below `GAP_DATA_ROOT`.

On Windows PowerShell the equivalent is:

```powershell
$env:GAP_DATA_ROOT = 'internal_authority/gap-data'
$env:GAP_RESULTS_ROOT = 'internal_authority/gap-results'
$env:GAP_CACHE_ROOT = 'internal_authority/gap-cache'
```


## Biological split

Paper-level evaluation is grouped by animal/fish/worm whenever the acquisition supports that split. Inner hyperparameter selection is restricted to the outer-training biological units. Any exception is explicitly labeled descriptive or sensitivity analysis.

## Provenance

`docs/PUBLIC_RELEASE_MANIFEST.csv` records SHA256 hashes for the staged package and canonical scripts. Canonical scripts are checked against the frozen manuscript code map before release.

## Results policy

Generated `results/`, large caches, checkpoints, logs and temporary files are ignored. Source-data-sized tables may be released separately with a figure/question provenance record.

## Release gate

Before a public tag:
1. portability/security audit has no unresolved private paths or credentials;
2. `python -m compileall gxa scripts` passes;
3. smoke tests pass on synthetic/small public examples;
4. every main figure/question has a public entry point and data note;
5. script hashes are reconciled with the manuscript authority.
