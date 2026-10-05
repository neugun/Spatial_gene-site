# Public-release data policy

## Included in the code repository

- reusable source code and dataset adapters;
- canonical paper-analysis entry points;
- compact configuration files, schemas, and small manifests needed to identify public inputs;
- source-data-sized result tables when release permissions allow;
- synthetic/small example data for testing;
- provenance hashes and question/figure-to-code maps.

## Not included by default

- raw or large derived datasets;
- model checkpoints and large cached embeddings;
- workstation `results/`, logs, temporary caches, or exploratory sweeps;
- unpublished/private experimental data;
- internal handoff documents, application/presentation materials, or private manuscript drafts;
- credentials, tokens, machine-specific paths, or software licenses.

## External public datasets

Public resources remain at their authoritative repository/accession. The GAP repository records the exact source/version and provides an adapter or preparation instruction rather than silently repackaging the underlying data.

## Inferred versus measured molecular data

Atlas/reference-derived molecular values are always labeled `inferred`. A continuous reference posterior, a Cre/driver label, a class-level external lookup, and same-cell measured expression are distinct evidence tiers and must not be merged into one “transcriptomics” label.
