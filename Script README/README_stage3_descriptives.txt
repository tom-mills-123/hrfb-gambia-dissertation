STAGE 3 — Descriptive Analysis (raw inputs for all primary results)

Inputs:
- outputs/step02_births_qc_enhanced.rds  (working dataset: raw + QA variables)

Policy:
- All Stage 3 prevalence estimates, subgroup comparisons, chi-square tests, and plots use raw variables only.
- The primary outcome is hrfb_binary (any HRFB).
- Component descriptives use the raw Stage 1 HRFB component indicators.
- Stage 2 banded and smoothed variables are retained in the dataset for sensitivity analyses later in the pipeline, but are not used here.

Tables:
- tables_stage3/stage3_national_prevalence_raw.csv
- tables_stage3/stage3_subgroup_prevalence_raw.csv
- tables_stage3/stage3_subgroup_raoscott_tests_raw.csv
- tables_stage3/stage3_components_by_residence_RAW.csv
- tables_stage3/stage3_region_prevalence_RAW.csv

Figures:
- PDF and PNG versions are written for all main plotted outputs in figs_stage3/
- figs_stage3/stage3_prev_by_residence_RAW.*
- figs_stage3/stage3_prev_by_region_RAW.*
- figs_stage3/stage3_prev_by_education_RAW.*
- figs_stage3/stage3_prev_by_wealth_RAW.*
- figs_stage3/stage3_prev_by_ethnicity_RAW.*
- figs_stage3/stage3_components_by_residence_RAW.*

Notes:
- Estimates are design-adjusted using weights, PSU identifiers (v021), and the carried-forward strata variable.
- Using the Stage 2 dataset preserves quality-control provenance without changing the primary descriptive definitions used in Stage 3.

Saved on: 2026-03-28 22:47:25 GMT
