STAGE 6 — Scenario modelling (exploratory; policy-focused)

Inputs:
- outputs/step02_births_qc_enhanced.rds (working dataset: raw + QA variables)

Method summary:
- Refit the Stage 4 primary survey-adjusted logistic model (svyglm with quasibinomial family).
- Generate counterfactual predictions using prediction-based g-computation.
- Average predicted probabilities using the survey design to obtain design-adjusted prevalence under each scenario.

Scenarios:
- Baseline: observed data with no intervention.
- S1 (+20pp contextual education): share_sec_plus := pmin(1, share_sec_plus + 0.20).
- S2 (+15pp contraceptive use): reassign approximately 15 percentage points of the weighted population from 'not using' to the weighted modal non-'not using' contraceptive category.
- S3 (Both): apply S1 and S2 simultaneously.

Outputs (tables):
- tables_stage6/stage6_refit_primary_logistic_OR.csv
- tables_stage6/stage6_scenario_prevalence_overall.csv
- tables_stage6/stage6_scenario_prevalence_changes.csv
- tables_stage6/stage6_scenario_prevalence_by_group.csv

Outputs (figures):
- PDF and PNG versions are written for all main plotted outputs in figs_stage6/
- figs_stage6/stage6_overall_scenario_prevalence.*
- figs_stage6/stage6_by_residence_scenario_prevalence.*
- figs_stage6/stage6_by_region_scenario_prevalence.*

Notes:
- The Stage 6 model specification matches Stage 4 to maintain consistency and avoid definitional leakage.
- These scenario estimates are exploratory model-based counterfactual predictions rather than identified causal effects.
- Predictions are averaged using survey weights and the complex survey design.

Saved on: 2026-03-29 10:47:06 BST
