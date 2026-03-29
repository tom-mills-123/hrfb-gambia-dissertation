STAGE 4 — Regression modelling (raw outcomes; design-adjusted)

Inputs:
- outputs/step02_births_qc_enhanced.rds (working dataset: raw + QA variables)

Primary model:
- svyglm with quasibinomial family, outcome: hrfb_binary.
- Predictors: residence, region, education, wealth, ethnicity, FP exposure (radio/TV/print),
  contraceptive_use, fert_pref, and share_sec_plus.
- Maternal age at birth, parity, and short birth interval are excluded because they define the HRFB outcome itself.

Secondary analyses:
- Model comparison summaries (null, core, full).
- Residual and functional-form diagnostics.
- Adjusted prevalence plots from the primary binary model.
- Optional multinomial modelling for hrfb_multicat.
- Sensitivity models for <25 and >=35 age strata.
- One-birth-per-woman sensitivity model.

Tables:
- tables_stage4/stage4_primary_logistic_OR.csv
- tables_stage4/appendix_fp_exposure_terms.csv
- tables_stage4/rq2_pseudoR2.csv
- tables_stage4/rq2_joint_wald_tests.csv
- tables_stage4/rq2_adjprev_by_residence.csv
- tables_stage4/rq2_adjprev_by_education.csv
- tables_stage4/rq2_adjprev_by_wealth.csv
- tables_stage4/rq2_adjprev_by_region.csv
- tables_stage4/rq2_adjprev_by_ethnicity.csv
- tables_stage4/rq2_adjprev_by_share_sec_plus.csv
- tables_stage4/stage4_primary_logistic_lt25_OR.csv (if fitted)
- tables_stage4/stage4_primary_logistic_ge35_OR.csv (if fitted)
- tables_stage4/stage4_primary_logistic_onebirth_per_woman_OR.csv
- tables_stage4/stage4_compare_main_vs_onebirth_OR.csv
- tables_stage4/stage4_multinomial_RRR.csv (if fitted)

Figures:
- PDF and PNG versions are written for all main plotted outputs in figs_stage4/
- figs_stage4/stage4_primary_logistic_forest.*
- figs_stage4/rq2_forest_fertpref_contraception.*
- figs_stage4/rq2_adjprev_by_residence.*
- figs_stage4/rq2_adjprev_by_education.*
- figs_stage4/rq2_adjprev_by_wealth.*
- figs_stage4/rq2_adjprev_by_region.*
- figs_stage4/rq2_adjprev_by_ethnicity.*
- figs_stage4/rq2_adjprev_by_share_sec_plus.*
- figs_stage4/rq2_resid_pearson_vs_fitted.*
- figs_stage4/rq2_resid_deviance_vs_fitted.*
- figs_stage4/rq2_resid_pearson_vs_share_sec_plus.*

Notes:
- All main inferences are design-adjusted using weights, PSU identifiers, and strata.
- QA variables carried forward from Stage 2 remain in the working dataset but are not used as primary outcomes.
- The multinomial model is treated as exploratory and may be omitted if instability prevents reliable estimation.

Saved on: 2026-03-28 23:16:14 GMT
