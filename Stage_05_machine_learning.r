############################################################
# STAGE 5 — Machine-learning classification (RQ3)
# Project: High-risk fertility behaviour (HRFB) — Gambia DHS 2019–20
#
# What this script does:
#   1) Load the Stage 2 enriched working dataset.
#   2) Build a machine-learning analysis dataset using the same
#      non-definitional predictors carried forward from Stage 4.
#   3) Split the data using PSU-held-out train/test partitioning
#      and PSU-grouped cross-validation to reduce overly optimistic
#      performance estimates.
#   4) Fit two classifiers:
#        - Random Forest
#        - XGBoost
#   5) Evaluate predictive performance using grouped cross-validation
#      and held-out PSU test data.
#   6) Produce ROC curves, calibration plots, confusion matrices,
#      Brier scores, and summary performance tables.
#   7) Export model-interpretation outputs for XGBoost:
#        - SHAP importance
#        - partial dependence plots
#   8) Write a Stage 5 README documenting the predictive workflow.
#
# Key principle for this stage:
#   This is a predictive extension of the dissertation rather than
#   a replacement for the design-based inferential models. The same
#   leakage-avoidance rule used in Stage 4 is retained here:
#   variables that directly define HRFB are not included as predictors.
############################################################

# --------------------------
# 0) Packages + global options
# --------------------------
AUTO_INSTALL <- FALSE

# These packages support data preparation, machine-learning workflows,
# performance evaluation, interpretation, and figure export.
pkgs <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "forcats",
  "readr",
  "janitor",
  "ggplot2",
  "recipes",
  "rsample",
  "yardstick",
  "parsnip",
  "workflows",
  "tune",
  "dials",
  "ranger",
  "xgboost",
  "vip",
  "pdp",
  "purrr",
  "scales",
  "Cairo"
)

missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))

# The script stops early if required packages are unavailable,
# unless AUTO_INSTALL has been switched on.
if (length(missing_pkgs) > 0) {
  if (isTRUE(AUTO_INSTALL)) {
    install.packages(missing_pkgs)
  } else {
    stop(
      "Missing packages: ", paste(missing_pkgs, collapse = ", "),
      "\nInstall them first, or set AUTO_INSTALL <- TRUE at the top of the script."
    )
  }
}

# Startup messages are suppressed to keep the log output compact.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(readr)
  library(janitor)
  library(ggplot2)
  library(recipes)
  library(rsample)
  library(yardstick)
  library(parsnip)
  library(workflows)
  library(tune)
  library(dials)
  library(ranger)
  library(xgboost)
  library(vip)
  library(pdp)
  library(purrr)
  library(scales)
  library(Cairo)
})

# --------------------------
# Helper functions
# --------------------------

# Define a consistent plotting theme used across dissertation figures.
theme_dissertation <- function(base_size = 11, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major.y = element_line(linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 10, 8, 10)
    )
}

# Save figures in both PDF and PNG formats.
save_fig <- function(p, filename, fig_dir, w = 8, h = 5) {
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  
  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".pdf")),
    plot = p,
    width = w,
    height = h,
    units = "in",
    device = Cairo::CairoPDF
  )
  
  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".png")),
    plot = p,
    width = w,
    height = h,
    units = "in",
    dpi = 450
  )
}

# A common caption is used across Stage 5 plots.
caption_rq3 <- paste0(
  "Data: The Gambia DHS 2019–20. Predictive models evaluated using PSU-held-out validation."
)

# This helper computes a few additional confusion-matrix metrics that are
# useful for reporting alongside the standard outputs.
calc_extra_metrics <- function(cm_obj) {
  tab <- as.data.frame(cm_obj$table)
  names(tab)[names(tab) == "Freq"] <- "n"
  
  TP <- tab$n[tab$Prediction == "yes" & tab$Truth == "yes"]
  TN <- tab$n[tab$Prediction == "no"  & tab$Truth == "no"]
  FP <- tab$n[tab$Prediction == "yes" & tab$Truth == "no"]
  FN <- tab$n[tab$Prediction == "no"  & tab$Truth == "yes"]
  
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  npv  <- TN / (TN + FN)
  bal_acc <- (sens + spec) / 2
  
  tibble::tibble(
    sensitivity = round(sens, 3),
    specificity = round(spec, 3),
    balanced_accuracy = round(bal_acc, 3),
    npv = round(npv, 3)
  )
}

# This helper calculates the Brier score for predicted probabilities.
brier_score <- function(df) {
  mean((as.numeric(df$hrfb_f == "yes") - df$.pred)^2)
}

# This helper draws a single ROC curve in the dissertation plot style.
plot_roc_curve <- function(roc_df, title, auc_val) {
  ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_abline(linetype = 2, linewidth = 0.4) +
    geom_path(linewidth = 0.9) +
    coord_equal() +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    labs(
      title = title,
      subtitle = paste0("ROC–AUC = ", round(auc_val, 3)),
      x = "1 − Specificity",
      y = "Sensitivity",
      caption = caption_rq3
    ) +
    theme_dissertation()
}

# This helper produces a decile-based calibration plot using held-out predictions.
calibration_plot <- function(df, model_name) {
  df_cal <- df %>%
    mutate(prob_bin = ntile(.pred, 10)) %>%
    group_by(prob_bin) %>%
    summarise(
      mean_pred = mean(.pred),
      obs_rate = mean(hrfb_f == "yes"),
      .groups = "drop"
    )
  
  ggplot(df_cal, aes(x = mean_pred, y = obs_rate)) +
    geom_abline(linetype = 2, linewidth = 0.4) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.9) +
    coord_equal() +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      title = paste("Calibration plot —", model_name),
      subtitle = "Held-out PSU test set",
      x = "Mean predicted probability",
      y = "Observed HRFB rate",
      caption = caption_rq3
    ) +
    theme_dissertation() +
    theme(
      plot.margin = margin(10, 10, 32, 10),
      plot.caption = element_text(margin = margin(t = 8))
    )
}

# This helper draws a confusion-matrix heatmap from a conf_mat object.
plot_confusion_matrix <- function(cm_obj, title) {
  df_cm <- as.data.frame(cm_obj$table)
  
  if ("Freq" %in% names(df_cm)) {
    df_cm <- df_cm %>% rename(n = Freq)
  } else if (!"n" %in% names(df_cm)) {
    stop(
      "Could not find confusion-matrix count column. Columns are: ",
      paste(names(df_cm), collapse = ", ")
    )
  }
  
  ggplot(df_cm, aes(x = Truth, y = Prediction, fill = n)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = n), size = 3.6) +
    labs(
      title = title,
      subtitle = "Counts on held-out PSU test set",
      x = "True class",
      y = "Predicted class",
      caption = caption_rq3
    ) +
    theme_dissertation() +
    theme(legend.position = "none")
}

# This helper gives more readable labels for machine-learning feature names
# after recipe-based dummy coding.
label_feature_ml <- function(x) {
  x %>%
    stringr::str_replace("^educ_no\\.education$", "Education: No education") %>%
    stringr::str_replace("^educ_primary$", "Education: Primary") %>%
    stringr::str_replace("^educ_secondary$", "Education: Secondary") %>%
    stringr::str_replace("^educ_higher$", "Education: Higher") %>%
    stringr::str_replace("^wealth_poorest$", "Wealth: Poorest") %>%
    stringr::str_replace("^wealth_poorer$", "Wealth: Poorer") %>%
    stringr::str_replace("^wealth_middle$", "Wealth: Middle") %>%
    stringr::str_replace("^wealth_richer$", "Wealth: Richer") %>%
    stringr::str_replace("^wealth_richest$", "Wealth: Richest") %>%
    stringr::str_replace("^fert_pref_have\\.another$", "Fertility preference: Have another") %>%
    stringr::str_replace("^fert_pref_no\\.more$", "Fertility preference: No more") %>%
    stringr::str_replace("^fert_pref_undecided$", "Fertility preference: Undecided") %>%
    stringr::str_replace("^fert_pref_missing$", "Fertility preference: Missing") %>%
    stringr::str_replace("^fert_pref_sterilized\\.infecund$", "Fertility preference: Sterilised/infecund") %>%
    stringr::str_replace("^fert_pref_other\\.infecund\\.missing$", "Fertility preference: Other/infecund/missing") %>%
    stringr::str_replace("^fp_radio_1$", "FP message exposure (radio): Not at all") %>%
    stringr::str_replace("^fp_radio_2$", "FP message exposure (radio): < once a week") %>%
    stringr::str_replace("^fp_radio_3$", "FP message exposure (radio): ≥ once a week") %>%
    stringr::str_replace("^fp_tv_1$", "FP message exposure (TV): Not at all") %>%
    stringr::str_replace("^fp_tv_2$", "FP message exposure (TV): < once a week") %>%
    stringr::str_replace("^fp_tv_3$", "FP message exposure (TV): ≥ once a week") %>%
    stringr::str_replace("^fp_print_1$", "FP message exposure (print): Not at all") %>%
    stringr::str_replace("^fp_print_2$", "FP message exposure (print): < once a week") %>%
    stringr::str_replace("^fp_print_3$", "FP message exposure (print): ≥ once a week") %>%
    stringr::str_replace("^contraceptive_use_not\\.using$", "Contraceptive use: Not using") %>%
    stringr::str_replace("^contraceptive_use_pill$", "Contraceptive use: Pill") %>%
    stringr::str_replace("^contraceptive_use_iud$", "Contraceptive use: IUD") %>%
    stringr::str_replace("^contraceptive_use_injections$", "Contraceptive use: Injections") %>%
    stringr::str_replace("^contraceptive_use_male\\.condom$", "Contraceptive use: Male condom") %>%
    stringr::str_replace("^contraceptive_use_implants\\.norplant$", "Contraceptive use: Implants/Norplant") %>%
    stringr::str_replace("^contraceptive_use_female\\.sterilization$", "Contraceptive use: Female sterilisation") %>%
    stringr::str_replace("^contraceptive_use_periodic\\.abstinence$", "Contraceptive use: Periodic abstinence") %>%
    stringr::str_replace("^contraceptive_use_withdrawal$", "Contraceptive use: Withdrawal") %>%
    stringr::str_replace("^contraceptive_use_other\\.traditional$", "Contraceptive use: Other traditional") %>%
    stringr::str_replace("^contraceptive_use_missing$", "Contraceptive use: Missing") %>%
    stringr::str_replace("^share_sec_plus$", "Cluster share with secondary+")
}

# This prediction wrapper returns XGBoost probabilities for use in PDPs.
xgb_pred_fun <- function(object, newdata) {
  as.numeric(predict(object, newdata = as.matrix(newdata)))
}

# --------------------------
# 1) Inputs and outputs
# --------------------------

# Stage 5 uses the Stage 2 enriched dataset as its working input.
# This preserves the quality-control trail while keeping the predictive
# outcome aligned with the raw HRFB definition.
in_births <- "outputs/step02_births_qc_enhanced.rds"
stopifnot(file.exists(in_births))

df <- readRDS(in_births)

# Output directories for Stage 5 tables, figures, and README.
out_dir <- "outputs"
tab_dir <- file.path(out_dir, "tables_stage5")
fig_dir <- file.path(out_dir, "figs_stage5")

if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# --------------------------
# 2) Outcome + predictors
# --------------------------

# These variables mirror the non-definitional predictors used in Stage 4.
# Variables that directly define HRFB are excluded to avoid leakage.
vars_keep <- c(
  "hrfb_binary",
  "v021",
  "residence", "region", "educ", "wealth", "ethnicity",
  "fp_radio", "fp_tv", "fp_print",
  "contraceptive_use", "fert_pref",
  "share_sec_plus"
)

miss <- setdiff(vars_keep, names(df))
if (length(miss) > 0) {
  stop("Missing required Stage 5 variables: ", paste(miss, collapse = ", "))
}

dat <- df %>%
  select(any_of(vars_keep)) %>%
  filter(!is.na(hrfb_binary), !is.na(v021)) %>%
  mutate(
    # Reference categories are aligned as closely as possible with Stage 4.
    residence = forcats::fct_relevel(residence, "urban"),
    educ = forcats::fct_relevel(educ, "higher", "secondary", "primary", "no education"),
    wealth = forcats::fct_relevel(wealth, "richest", "richer", "middle", "poorer", "poorest"),
    region = factor(stringr::str_to_title(as.character(region))),
    ethnicity = factor(stringr::str_to_title(as.character(ethnicity))),
    fp_radio = forcats::fct_explicit_na(fp_radio, na_level = "Missing"),
    fp_tv = forcats::fct_explicit_na(fp_tv, na_level = "Missing"),
    fp_print = forcats::fct_explicit_na(fp_print, na_level = "Missing"),
    contraceptive_use = forcats::fct_explicit_na(contraceptive_use, na_level = "Missing"),
    fert_pref = forcats::fct_explicit_na(fert_pref, na_level = "Missing"),
    
    # The event level is set to "yes" and placed first so yardstick metrics
    # interpret the positive class consistently.
    hrfb_f = factor(hrfb_binary, levels = c(1, 0), labels = c("yes", "no"))
  ) %>%
  mutate(
    # FP exposure levels are explicitly ordered for consistency.
    fp_radio = factor(
      fp_radio,
      levels = c("not at all", "less than once a week", "at least once a week", "Missing")
    ),
    fp_tv = factor(
      fp_tv,
      levels = c("not at all", "less than once a week", "at least once a week", "Missing")
    ),
    fp_print = factor(
      fp_print,
      levels = c("not at all", "less than once a week", "at least once a week", "Missing")
    )
  ) %>%
  mutate(
    across(where(is.factor), forcats::fct_drop)
  )

# A simple missingness snapshot is saved before recipe-based preprocessing.
message("\n[Stage 5] Missingness snapshot (pre-recipe):")
miss_tbl <- dat %>%
  summarise(across(everything(), ~ sum(is.na(.))))
print(miss_tbl, n = Inf, width = Inf)
readr::write_csv(miss_tbl, file.path(tab_dir, "stage5_missingness_snapshot.csv"))

# --------------------------
# 3) Group-aware train/test split + group cross-validation
# --------------------------
set.seed(2025)

# Entire PSUs are held out at test time so the evaluation reflects
# out-of-cluster performance rather than within-cluster memorisation.
split <- rsample::group_initial_split(dat, group = v021, prop = 0.7)

train <- rsample::training(split)
test  <- rsample::testing(split)

message("\n[Stage 5] Train/Test sizes (group split by PSU):")
print(
  tibble::tibble(
    set = c("train", "test"),
    n = c(nrow(train), nrow(test)),
    pos_rate = c(mean(train$hrfb_f == "yes"), mean(test$hrfb_f == "yes"))
  )
)

# This overlap check should return zero if the PSU split worked correctly.
psu_overlap <- length(intersect(train$v021, test$v021))
message("[Stage 5] PSU overlap between train and test: ", psu_overlap)

message("[Stage 5] Unique PSUs in train: ", dplyr::n_distinct(train$v021))
message("[Stage 5] Unique PSUs in test:  ", dplyr::n_distinct(test$v021))

set.seed(2025)

# Group-based v-fold cross-validation is created on the training set only.
cv10 <- rsample::group_vfold_cv(train, group = v021, v = 10)
message("[Stage 5] 10-fold group CV created on training set.")

# Fold-level balance is saved because it helps document whether outcome
# prevalence remains broadly similar across the held-out folds.
fold_balance <- cv10 %>%
  mutate(
    train_pos = purrr::map_dbl(splits, ~ mean(rsample::analysis(.)$hrfb_f == "yes")),
    test_pos  = purrr::map_dbl(splits, ~ mean(rsample::assessment(.)$hrfb_f == "yes")),
    train_n   = purrr::map_int(splits, ~ nrow(rsample::analysis(.))),
    test_n    = purrr::map_int(splits, ~ nrow(rsample::assessment(.)))
  ) %>%
  select(id, train_n, test_n, train_pos, test_pos)

print(fold_balance)
readr::write_csv(fold_balance, file.path(tab_dir, "stage5_groupcv_fold_balance.csv"))

# --------------------------
# 4) Preprocessing recipe
# --------------------------

# The recipe handles missing values, dummy coding, zero-variance removal,
# and numeric standardisation. It is fit on the training data only.
rec <- recipes::recipe(
  hrfb_f ~ residence + region + educ + wealth + ethnicity +
    fp_radio + fp_tv + fp_print +
    contraceptive_use + fert_pref + share_sec_plus,
  data = train
) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

rec_prep <- prep(rec)
train_bake <- bake(rec_prep, new_data = train)
test_bake  <- bake(rec_prep, new_data = test)

message("\n[Stage 5] Feature matrix dimensions after preprocessing:")
print(
  tibble::tibble(
    dataset = c("train_bake", "test_bake"),
    rows = c(nrow(train_bake), nrow(test_bake)),
    cols = c(ncol(train_bake), ncol(test_bake))
  )
)

# --------------------------
# 5) Model specifications
# --------------------------

# A random forest and an XGBoost model are fit as complementary learners.
# The specifications here are deliberately moderate rather than heavily tuned.
rf_spec <- parsnip::rand_forest(
  mode = "classification",
  trees = 1000,
  mtry = min(8, ncol(train_bake) - 1),
  min_n = 5
) %>%
  set_engine("ranger", importance = "impurity")

xgb_spec <- parsnip::boost_tree(
  mode = "classification",
  trees = 1000,
  learn_rate = 0.05,
  tree_depth = 4,
  loss_reduction = 0,
  min_n = 10
) %>%
  set_engine("xgboost", nthread = 0)

wf_rf <- workflows::workflow() %>%
  add_model(rf_spec) %>%
  add_recipe(rec)

wf_xgb <- workflows::workflow() %>%
  add_model(xgb_spec) %>%
  add_recipe(rec)

# --------------------------
# 6) Cross-validated performance
# --------------------------

# The same metric set is used for both CV and final held-out evaluation.
metric_set_all <- yardstick::metric_set(
  accuracy,
  yardstick::roc_auc,
  sens,
  spec,
  precision,
  recall,
  f_meas
)

set.seed(2025)
rf_cv <- tune::fit_resamples(
  wf_rf,
  resamples = cv10,
  metrics = metric_set_all,
  control = control_resamples(save_pred = TRUE)
)

set.seed(2025)
xgb_cv <- tune::fit_resamples(
  wf_xgb,
  resamples = cv10,
  metrics = metric_set_all,
  control = control_resamples(save_pred = TRUE)
)

rf_cv_metrics <- collect_metrics(rf_cv)
xgb_cv_metrics <- collect_metrics(xgb_cv)

readr::write_csv(rf_cv_metrics, file.path(tab_dir, "stage5_rf_cv_metrics.csv"))
readr::write_csv(xgb_cv_metrics, file.path(tab_dir, "stage5_xgb_cv_metrics.csv"))

message("\n[Stage 5] Cross-validated metrics — Random Forest:")
print(rf_cv_metrics %>% arrange(.metric))

message("\n[Stage 5] Cross-validated metrics — XGBoost:")
print(xgb_cv_metrics %>% arrange(.metric))

# --------------------------
# 7) Final fit on training data + held-out test evaluation
# --------------------------

# Each workflow is refit on the full training set and then evaluated on
# the PSU-held-out test set.
rf_fit <- fit(wf_rf, data = train)
xgb_fit <- fit(wf_xgb, data = train)

rf_pred_class <- predict(rf_fit, new_data = test, type = "class")
rf_pred_prob  <- predict(rf_fit, new_data = test, type = "prob") %>%
  select(.pred_yes)

xgb_pred_class <- predict(xgb_fit, new_data = test, type = "class")
xgb_pred_prob  <- predict(xgb_fit, new_data = test, type = "prob") %>%
  select(.pred_yes)

rf_res <- bind_cols(test, rf_pred_class, rf_pred_prob) %>%
  rename(.pred = .pred_yes)

xgb_res <- bind_cols(test, xgb_pred_class, xgb_pred_prob) %>%
  rename(.pred = .pred_yes)

rf_test_metrics <- metric_set_all(
  rf_res,
  truth = hrfb_f,
  estimate = .pred_class,
  .pred,
  event_level = "first"
)

xgb_test_metrics <- metric_set_all(
  xgb_res,
  truth = hrfb_f,
  estimate = .pred_class,
  .pred,
  event_level = "first"
)

readr::write_csv(rf_test_metrics, file.path(tab_dir, "stage5_rf_test_metrics.csv"))
readr::write_csv(xgb_test_metrics, file.path(tab_dir, "stage5_xgb_test_metrics.csv"))

message("\n[Stage 5] Test-set metrics — Random Forest:")
print(rf_test_metrics %>% arrange(.metric))

message("\n[Stage 5] Test-set metrics — XGBoost:")
print(xgb_test_metrics %>% arrange(.metric))

# --------------------------
# 8) ROC curves
# --------------------------

# ROC curves are drawn separately and jointly for the held-out PSU test set.
roc_rf <- yardstick::roc_curve(rf_res, truth = hrfb_f, .pred)
roc_xgb <- yardstick::roc_curve(xgb_res, truth = hrfb_f, .pred)

auc_rf <- yardstick::roc_auc(rf_res, truth = hrfb_f, .pred, event_level = "first")$.estimate
auc_xgb <- yardstick::roc_auc(xgb_res, truth = hrfb_f, .pred, event_level = "first")$.estimate

p_rf_roc <- plot_roc_curve(roc_rf, "ROC curve — Random Forest", auc_rf)
p_xgb_roc <- plot_roc_curve(roc_xgb, "ROC curve — XGBoost", auc_xgb)

print(p_rf_roc)
save_fig(p_rf_roc, "stage5_rf_roc", fig_dir = fig_dir, w = 6.5, h = 5.5)

print(p_xgb_roc)
save_fig(p_xgb_roc, "stage5_xgb_roc", fig_dir = fig_dir, w = 6.5, h = 5.5)

roc_both <- bind_rows(
  roc_rf %>% mutate(model = paste0("Random Forest (AUC = ", round(auc_rf, 3), ")")),
  roc_xgb %>% mutate(model = paste0("XGBoost (AUC = ", round(auc_xgb, 3), ")"))
)

p_roc_combined <- ggplot(roc_both, aes(x = 1 - specificity, y = sensitivity, linetype = model)) +
  geom_abline(linetype = 2, linewidth = 0.4) +
  geom_path(linewidth = 0.95) +
  coord_equal() +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "ROC curves — Random Forest vs XGBoost (RQ3)",
    subtitle = "Performance evaluated on held-out PSUs",
    x = "1 − Specificity",
    y = "Sensitivity",
    linetype = NULL,
    caption = caption_rq3
  ) +
  theme_dissertation() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(10, 10, 28, 10),
    plot.caption = element_text(margin = margin(t = 8))
  )

print(p_roc_combined)
save_fig(p_roc_combined, "stage5_roc_combined_rf_vs_xgb", fig_dir = fig_dir, w = 7.2, h = 6.5)

# --------------------------
# 9) Baseline benchmark
# --------------------------

# A naive benchmark is included so the model accuracy can be interpreted
# against the simplest possible classifier.
test_prev <- mean(test$hrfb_f == "yes")

baseline_tbl <- tibble::tibble(
  metric = c(
    "Test prevalence (HRFB = yes)",
    "Naive accuracy (always predict yes)"
  ),
  value = c(round(test_prev, 3), round(test_prev, 3))
)

print(baseline_tbl)
readr::write_csv(baseline_tbl, file.path(tab_dir, "stage5_baseline_benchmark.csv"))

# --------------------------
# 10) Confusion matrices + additional metrics
# --------------------------

# Confusion matrices help show the balance between positive and negative
# classification performance more directly than a single summary metric.
cm_rf <- yardstick::conf_mat(rf_res, truth = hrfb_f, estimate = .pred_class)
cm_xgb <- yardstick::conf_mat(xgb_res, truth = hrfb_f, estimate = .pred_class)

message("\n[Stage 5] Confusion matrix — Random Forest:")
print(cm_rf)

message("\n[Stage 5] Confusion matrix — XGBoost:")
print(cm_xgb)

readr::write_csv(as.data.frame(cm_rf$table), file.path(tab_dir, "stage5_rf_confusion_matrix_test.csv"))
readr::write_csv(as.data.frame(cm_xgb$table), file.path(tab_dir, "stage5_xgb_confusion_matrix_test.csv"))

extra_rf <- calc_extra_metrics(cm_rf) %>%
  mutate(model = "Random Forest")

extra_xgb <- calc_extra_metrics(cm_xgb) %>%
  mutate(model = "XGBoost")

extra_metrics_tbl <- bind_rows(extra_rf, extra_xgb)

print(extra_metrics_tbl)
readr::write_csv(extra_metrics_tbl, file.path(tab_dir, "stage5_additional_metrics.csv"))

p_cm_rf <- plot_confusion_matrix(cm_rf, "Confusion matrix — Random Forest")
p_cm_xgb <- plot_confusion_matrix(cm_xgb, "Confusion matrix — XGBoost")

print(p_cm_rf)
save_fig(p_cm_rf, "stage5_rf_confusion_matrix", fig_dir = fig_dir, w = 6, h = 5)

print(p_cm_xgb)
save_fig(p_cm_xgb, "stage5_xgb_confusion_matrix", fig_dir = fig_dir, w = 6, h = 5)

# --------------------------
# 11) Brier scores + calibration plots
# --------------------------

# Brier scores summarise overall probability calibration.
brier_tbl <- tibble::tibble(
  model = c("Random Forest", "XGBoost"),
  brier_score = c(
    round(brier_score(rf_res), 4),
    round(brier_score(xgb_res), 4)
  )
)

print(brier_tbl)
readr::write_csv(brier_tbl, file.path(tab_dir, "stage5_brier_scores.csv"))

p_cal_rf <- calibration_plot(rf_res, "Random Forest")
p_cal_xgb <- calibration_plot(xgb_res, "XGBoost")

print(p_cal_rf)
save_fig(p_cal_rf, "stage5_rf_calibration", fig_dir = fig_dir, w = 6.8, h = 6.0)

print(p_cal_xgb)
save_fig(p_cal_xgb, "stage5_xgb_calibration", fig_dir = fig_dir, w = 6.8, h = 6.0)


# --------------------------
# 12) XGBoost SHAP importance
# --------------------------

# SHAP values are extracted from the fitted XGBoost booster and summarised
# as mean absolute contribution on the held-out test set.
xgb_fit_extracted <- workflows::extract_fit_parsnip(xgb_fit)
xgb_booster <- xgb_fit_extracted$fit

X_test <- bake(rec_prep, new_data = test) %>%
  dplyr::select(-hrfb_f)

shap <- predict(xgb_booster, newdata = as.matrix(X_test), predcontrib = TRUE)
shap_df <- as.data.frame(shap)

if (".BIAS" %in% names(shap_df)) shap_df <- shap_df %>% dplyr::select(-.BIAS)
if ("BIAS"  %in% names(shap_df)) shap_df <- shap_df %>% dplyr::select(-BIAS)

shap_imp <- data.frame(
  feature = names(shap_df),
  mean_abs_shap = apply(shap_df, 2, function(x) mean(abs(x), na.rm = TRUE))
) %>%
  tibble::as_tibble() %>%
  dplyr::arrange(dplyr::desc(mean_abs_shap)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::mutate(feature_label = label_feature_ml(feature))

message("\n[Stage 5] XGBoost SHAP importance (top rows):")
print(head(shap_imp, 10))
readr::write_csv(shap_imp, file.path(tab_dir, "stage5_xgb_shap_importance_top15.csv"))

p_shap <- ggplot(shap_imp, aes(x = reorder(feature_label, mean_abs_shap), y = mean_abs_shap)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.001),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "XGBoost — SHAP importance (top 15)",
    subtitle = "Mean absolute SHAP contribution on held-out test set",
    x = NULL,
    y = "Mean |SHAP|",
    caption = caption_rq3
  ) +
  theme_dissertation() +
  theme(
    axis.text.y = element_text(size = 9, margin = margin(r = 8)),
    plot.margin = margin(10, 10, 10, 20)
  )

print(p_shap)
save_fig(p_shap, "stage5_xgb_shap_importance_top15", fig_dir = fig_dir, w = 8.5, h = 6.2)

# --------------------------
# 13) XGBoost partial dependence plots
# --------------------------

# PDPs are used here as simple model-interpretation graphics for a small
# set of selected features.
key_feats <- c("educ_secondary", "educ_primary", "wealth_richest", "share_sec_plus")

X_train_df <- train_bake %>%
  dplyr::select(-hrfb_f)

for (f in key_feats) {
  if (f %in% colnames(X_train_df)) {
    pd <- pdp::partial(
      object = xgb_booster,
      pred.var = f,
      train = X_train_df,
      pred.fun = xgb_pred_fun,
      grid.resolution = 20,
      ice = FALSE
    )
    
    p <- ggplot(pd, aes(x = .data[[f]], y = .data[["yhat"]])) +
      geom_line(linewidth = 0.9) +
      {
        if (is.numeric(X_train_df[[f]])) {
          geom_rug(
            data = X_train_df,
            aes(x = .data[[f]]),
            inherit.aes = FALSE,
            alpha = 0.15
          )
        } else {
          geom_rug(
            data = X_train_df %>% mutate(`..tmp..` = as.numeric(.data[[f]])),
            aes(x = `..tmp..`),
            inherit.aes = FALSE,
            alpha = 0.15
          )
        }
      } +
      scale_y_continuous(
        labels = scales::number_format(accuracy = 0.01),
        expand = expansion(mult = c(0.02, 0.08))
      ) +
      labs(
        title = paste("Partial dependence —", f),
        subtitle = "XGBoost; marginal effect with other features at observed distribution",
        x = f,
        y = "Predicted P(HRFB = yes)",
        caption = caption_rq3
      ) +
      theme_dissertation()
    
    print(p)
    save_fig(p, paste0("stage5_xgb_pdp_", f), fig_dir = fig_dir, w = 6.5, h = 5.5)
  }
}

# --------------------------
# 14) Concise held-out performance summary
# --------------------------

# A wide summary table is exported to make model comparison easy in the write-up.
perf_tbl <- bind_rows(
  rf_test_metrics %>% mutate(model = "Random Forest"),
  xgb_test_metrics %>% mutate(model = "XGBoost")
) %>%
  select(model, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate) %>%
  mutate(
    accuracy = round(accuracy, 3),
    roc_auc = round(roc_auc, 3),
    sens = round(sens, 3),
    spec = round(spec, 3),
    precision = round(precision, 3),
    recall = round(recall, 3),
    f_meas = round(f_meas, 3)
  )

message("\n[Stage 5] Test performance summary:")
print(perf_tbl, n = Inf, width = Inf)
readr::write_csv(perf_tbl, file.path(tab_dir, "stage5_test_performance_summary.csv"))

# --------------------------
# 15) README
# --------------------------

# The README records the predictive design, leakage-avoidance rule,
# and the main exported outputs from this stage.
readme_path <- file.path(out_dir, "README_stage5_ml.txt")
readme_txt <- c(
  "STAGE 5 — Machine-learning classification (Random Forest and XGBoost)",
  "",
  "Inputs:",
  "- outputs/step02_births_qc_enhanced.rds (working dataset; raw HRFB outcome carried forward from Stage 1 and Stage 2).",
  "",
  "Outcome and predictors:",
  "- Outcome: hrfb_binary, stored here as factor hrfb_f with 'yes' as the event class.",
  "- Predictors: residence, region, education, wealth, ethnicity, FP message exposure (radio/TV/print),",
  "  contraceptive_use, fertility preference, and share_sec_plus.",
  "- Exclusion rule: age-at-birth, parity, and short birth interval are not used as predictors because they define HRFB.",
  "- v021 (PSU) is used only for grouped splitting and validation, not as a predictor.",
  "",
  "Methods:",
  "- Train/test split: 70/30 PSU-held-out split using rsample::group_initial_split(group = v021).",
  "- Cross-validation: 10-fold PSU-grouped CV on the training set.",
  "- Models: Random Forest (ranger) and XGBoost.",
  "- Metrics: Accuracy, ROC-AUC, Sensitivity, Specificity, Precision, Recall, F1, Brier score, and confusion matrices.",
  "- Interpretation: XGBoost SHAP importance and XGBoost partial dependence plots.",
  "",
  "Outputs (tables):",
  paste0("- ", file.path("tables_stage5", "stage5_missingness_snapshot.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_groupcv_fold_balance.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_rf_cv_metrics.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_xgb_cv_metrics.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_rf_test_metrics.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_xgb_test_metrics.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_baseline_benchmark.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_additional_metrics.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_brier_scores.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_rf_confusion_matrix_test.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_xgb_confusion_matrix_test.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_xgb_shap_importance_top15.csv")),
  paste0("- ", file.path("tables_stage5", "stage5_test_performance_summary.csv")),
  "",
  "Outputs (figures):",
  "- PDF and PNG versions are written for all main plotted outputs in figs_stage5/",
  paste0("- ", file.path("figs_stage5", "stage5_rf_roc.*")),
  paste0("- ", file.path("figs_stage5", "stage5_xgb_roc.*")),
  paste0("- ", file.path("figs_stage5", "stage5_roc_combined_rf_vs_xgb.*")),
  paste0("- ", file.path("figs_stage5", "stage5_rf_confusion_matrix.*")),
  paste0("- ", file.path("figs_stage5", "stage5_xgb_confusion_matrix.*")),
  paste0("- ", file.path("figs_stage5", "stage5_rf_calibration.*")),
  paste0("- ", file.path("figs_stage5", "stage5_xgb_calibration.*")),
  paste0("- ", file.path("figs_stage5", "stage5_xgb_shap_importance_top15.*")),
  paste0("- ", file.path("figs_stage5", "stage5_xgb_pdp_*.png/pdf")),
  "",
  "Key interpretation note:",
  "- Performance estimates reflect out-of-cluster generalisation because entire PSUs are held out at test time.",
  "- This reduces the risk of overly optimistic performance that could arise if births from the same PSU appeared in both training and test data.",
  "",
  paste0("Saved on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)

writeLines(readme_txt, readme_path)

message("Wrote README: ", normalizePath(readme_path))
message("Stage 5 complete. Outputs written to: ", normalizePath(out_dir))
############################################################
