STAGE 6A — Scenario Modelling (Exploratory sensitivity analysis)

Method summary:
- Refit the Stage 4 primary design-adjusted logistic model (svyglm quasibinomial).
- Evaluate alternative counterfactual scenarios using design-adjusted g-computation.
- This script differs from Stage 6 by replacing the contextual education scenario
  with an individual-level education scenario.

Scenarios:
- S1 (+20pp individual secondary+): reassign approximately 20 percentage points of the
  weighted population from 'no education' first, then 'primary', into 'secondary'.
- S2 (+15pp contraceptive use): reassign approximately 15 percentage points of the
  weighted population from 'not using' to the weighted modal non-'not using' method.
- S3 (Both): apply S1 and S2 simultaneously.

Interpretation note:
- The education scenario is a long-run structural counterfactual, not a short-run policy effect.
- Results remain associational and should not be interpreted causally.

Key files written:
- tables_stage6a/stage6a_refit_primary_logistic_OR.csv
- tables_stage6a/stage6a_scenario_prevalence_overall.csv
- tables_stage6a/stage6a_scenario_prevalence_by_group.csv
- tables_stage6a/stage6_vs_stage6a_overall_comparison.csv
- figs_stage6a/stage6a_overall_scenario_prevalence.png
- figs_stage6a/stage6a_by_residence_scenario_prevalence.png
- figs_stage6a/stage6a_by_region_scenario_prevalence.png
