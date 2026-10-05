# Fig4-7 full dataset acceptance matrix — 2026-10-04

This matrix covers every dataset in the adapter-universe audit plus the independent Fig1-3 integration datasets. Applicability is separated from execution and from evidence strength. A blocked or N/A row is not a failed biological result.

## Summary

figure / applicability / n_datasets
figure           applicability  n_datasets
  Fig4              APPLICABLE          13
  Fig4          NOT_APPLICABLE          13
  Fig4         BLOCKED_ADAPTER          10
  Fig4                EXECUTED           2
  Fig4     EXECUTED_STRUCTURAL           1
  Fig5          NOT_APPLICABLE          14
  Fig5              APPLICABLE          13
  Fig5         BLOCKED_ADAPTER          10
  Fig5                EXECUTED           2
  Fig6     APPLICABLE_PRECHECK          14
  Fig6         BLOCKED_ADAPTER          10
  Fig6        BLOCKED_PRECHECK           8
  Fig6          NOT_APPLICABLE           6
  Fig6                EXECUTED           1
  Fig7 NOT_APPLICABLE_UPSTREAM          13
  Fig7             CONDITIONAL          12
  Fig7         BLOCKED_ADAPTER          10
  Fig7     NOT_APPLICABLE_GATE           2
  Fig7                EXECUTED           1
  Fig7          NOT_APPLICABLE           1

## Current end-to-end acceptance datasets

- Allen Visual Learning HCR: full Fig4→Fig7 one-command continuation with perturbation atlas.
- McLachlan: Fig4/5 execute; Fig6 coverage preflight and Fig7 held-unit gate stop further claims.
- Condylis: Fig4/5 execute; sparse marker panel blocks Fig6 and fails Fig7 operator gate.
- MICrONS: structural/session-held route only; molecular stages are explicitly N/A.

Machine-readable matrix: internal_authority/FIG4_7_FULL_DATASET_ACCEPTANCE_MATRIX_20261004.csv