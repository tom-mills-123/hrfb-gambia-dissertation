STAGE 5 — Machine-learning classification (Random Forest and XGBoost)

Inputs:
- outputs/step02_births_qc_enhanced.rds (working dataset; raw HRFB outcome carried forward from Stage 1 and Stage 2).

Outcome and predictors:
- Outcome: hrfb_binary, stored here as factor hrfb_f with 'yes' as the event class.
- Predictors: residence, region, education, wealth, ethnicity, FP message exposure (radio/TV/print),
  contraceptive_use, fertility preference, and share_sec_plus.
- Exclusion rule: age-at-birth, parity, and short birth interval are not used as predictors because they define HRFB.
- v021 (PSU) is used only for grouped splitting and validation, not as a predictor.

Methods:
- Train/test split: 70/30 PSU-held-out split using rsample::group_initial_split(group = v021).
- Cross-validation: 10-fold PSU-grouped CV on the training set.
- Models: Random Forest (ranger) and XGBoost.
- Metrics: Accuracy, ROC-AUC, Sensitivity, Specificity, Precision, Recall, F1, Brier score, and confusion matrices.
- Interpretation: XGBoost SHAP importance and XGBoost partial dependence plots.

Outputs (tables):
- tables_stage5/stage5_missingness_snapshot.csv
- tables_stage5/stage5_groupcv_fold_balance.csv
- tables_stage5/stage5_rf_cv_metrics.csv
- tables_stage5/stage5_xgb_cv_metrics.csv
- tables_stage5/stage5_rf_test_metrics.csv
- tables_stage5/stage5_xgb_test_metrics.csv
- tables_stage5/stage5_baseline_benchmark.csv
- tables_stage5/stage5_additional_metrics.csv
- tables_stage5/stage5_brier_scores.csv
- tables_stage5/stage5_rf_confusion_matrix_test.csv
- tables_stage5/stage5_xgb_confusion_matrix_test.csv
- tables_stage5/stage5_xgb_shap_importance_top15.csv
- tables_stage5/stage5_test_performance_summary.csv

Outputs (figures):
- PDF and PNG versions are written for all main plotted outputs in figs_stage5/
- figs_stage5/stage5_rf_roc.*
- figs_stage5/stage5_xgb_roc.*
- figs_stage5/stage5_roc_combined_rf_vs_xgb.*
- figs_stage5/stage5_rf_confusion_matrix.*
- figs_stage5/stage5_xgb_confusion_matrix.*
- figs_stage5/stage5_rf_calibration.*
- figs_stage5/stage5_xgb_calibration.*
- figs_stage5/stage5_xgb_shap_importance_top15.*
- figs_stage5/stage5_xgb_pdp_*.png/pdf

Key interpretation note:
- Performance estimates reflect out-of-cluster generalisation because entire PSUs are held out at test time.
- This reduces the risk of overly optimistic performance that could arise if births from the same PSU appeared in both training and test data.

Saved on: 2026-03-29 10:18:39 BST
