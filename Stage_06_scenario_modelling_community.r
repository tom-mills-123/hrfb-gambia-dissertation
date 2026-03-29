############################################################
# STAGE 6 — Scenario modelling (exploratory; policy-focused)
# Project: High-risk fertility behaviour (HRFB) — Gambia DHS 2019–20
#
# What this script does:
#   1) Load the Stage 2 enriched working dataset.
#   2) Rebuild the survey design directly on that dataset so that
#      all scenario estimates remain aligned to the current analysis file.
#   3) Refit the Stage 4 primary survey-adjusted logistic model.
#   4) Define three policy-relevant counterfactual scenarios:
#        - S1: +20 percentage points in cluster share_sec_plus
#        - S2: +15 percentage points in contraceptive use
#        - S3: S1 and S2 combined
#   5) Use design-aware g-computation to estimate predicted HRFB
#      prevalence under each scenario.
#   6) Export overall and subgroup scenario tables and produce
#      publication-ready figures.
#   7) Write a Stage 6 README documenting the scenarios and outputs.
#
# Key principle for this stage:
#   This is an exploratory prediction-based extension of the main
#   design-adjusted regression analysis from Stage 4. The model
#   specification remains aligned with Stage 4, and predictors that
#   directly define HRFB are not introduced into the counterfactual
#   scenarios in order to avoid definitional leakage.
############################################################

# --------------------------
# 0) Packages + global options
# --------------------------
AUTO_INSTALL <- FALSE

# These packages support data preparation, survey-weighted modelling,
# prediction-based g-computation, and figure export.
pkgs <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "forcats",
  "readr",
  "janitor",
  "survey",
  "srvyr",
  "broom",
  "ggplot2",
  "scales",
  "purrr",
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
  library(survey)
  library(srvyr)
  library(broom)
  library(ggplot2)
  library(scales)
  library(purrr)
  library(Cairo)
})

# This matches the survey option used in earlier survey-based stages.
options(survey.lonely.psu = "adjust")

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
      plot.margin = margin(10, 12, 12, 12)
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

# Common caption used across Stage 6 plots.
caption_rq_policy <- paste0(
  "Data: The Gambia DHS 2019–20. Design-adjusted g-computation using predictions from the Stage 4 primary model.\n",
  "Error bars are 95% Wald confidence intervals."
)

# limit numeric values to the unit interval so scenario shifts do not
# produce impossible proportions below 0 or above 1.
cap01 <- function(x) {
  pmax(0, pmin(1, x))
}

# Align the survey design columns from common DHS-style names.
# This makes the design rebuilding step more robust if naming varies.
detect_design_columns <- function(dat) {
  candidates <- list(
    ids = c("v021", "psu", "cluster", "hv021"),
    strata = c("strata", "v022", "v023", "hv022", "stratum"),
    weights = c("weight", "v005", "wt", "wgt", "perweight")
  )
  
  pick_first_match <- function(x) {
    found <- x[x %in% names(dat)]
    if (length(found) > 0) found[1] else NA_character_
  }
  
  ids_col <- pick_first_match(candidates$ids)
  strata_col <- pick_first_match(candidates$strata)
  weight_col <- pick_first_match(candidates$weights)
  
  if (anyNA(c(ids_col, strata_col, weight_col))) {
    stop(
      "Could not detect survey design columns. Found ids = ", ids_col,
      ", strata = ", strata_col,
      ", weights = ", weight_col,
      ". Available columns are: ", paste(names(dat), collapse = ", ")
    )
  }
  
  list(ids = ids_col, strata = strata_col, weights = weight_col)
}

# Build a fresh survey design object directly from the supplied dataset.
make_design <- function(dat) {
  dc <- detect_design_columns(dat)
  
  dat[[dc$weights]] <- suppressWarnings(as.numeric(dat[[dc$weights]]))
  if (any(!is.finite(dat[[dc$weights]]))) {
    stop("Non-finite values detected in survey weights column: ", dc$weights)
  }
  
  # DHS weights are stored on a 1e6 scale and need rescaling.
  if (max(dat[[dc$weights]], na.rm = TRUE) > 10) {
    dat[[dc$weights]] <- dat[[dc$weights]] / 1e6
  }
  
  dat[[dc$ids]] <- droplevels(as.factor(dat[[dc$ids]]))
  dat[[dc$strata]] <- droplevels(as.factor(dat[[dc$strata]]))
  
  survey::svydesign(
    ids = stats::as.formula(paste0("~", dc$ids)),
    strata = stats::as.formula(paste0("~", dc$strata)),
    weights = stats::as.formula(paste0("~", dc$weights)),
    data = dat,
    nest = TRUE
  )
}

# Generate predicted probabilities from a fitted svyglm model
safe_predict_prob <- function(model, newdata) {
  tt <- stats::terms(model)
  
  mf <- stats::model.frame(
    stats::delete.response(tt),
    newdata,
    na.action = stats::na.pass,
    xlev = model$xlevels
  )
  
  mm <- stats::model.matrix(
    stats::delete.response(tt),
    mf,
    contrasts.arg = model$contrasts
  )
  
  cf <- stats::coef(model)
  
  # Add any coefficient columns missing from the new model matrix.
  missing_cols <- setdiff(names(cf), colnames(mm))
  if (length(missing_cols) > 0) {
    add <- matrix(0, nrow = nrow(mm), ncol = length(missing_cols))
    colnames(add) <- missing_cols
    mm <- cbind(mm, add)
  }
  
  # Remove columns in the model matrix that are not in the fitted model.
  extra_cols <- setdiff(colnames(mm), names(cf))
  if (length(extra_cols) > 0) {
    mm <- mm[, setdiff(colnames(mm), extra_cols), drop = FALSE]
  }
  
  # Force exact coefficient order before matrix multiplication.
  mm <- mm[, names(cf), drop = FALSE]
  
  eta <- as.numeric(mm %*% cf)
  p <- model$family$linkinv(eta)
  p[!is.finite(p)] <- NA_real_
  
  p
}

# Compute the design-weighted mean predicted prevalence from a scenario
# dataset using g-computation. Predicted probabilities are averaged using
# the survey design weights rather than simple arithmetic means.
predict_prevalence <- function(d_new, model) {
  p <- safe_predict_prob(model, d_new)
  
  des_new <- make_design(d_new)
  des_new$variables$p <- p
  
  est_obj <- survey::svymean(~p, design = des_new, na.rm = TRUE)
  
  estimate <- as.numeric(stats::coef(est_obj))
  se <- as.numeric(survey::SE(est_obj))
  low <- estimate - 1.96 * se
  upp <- estimate + 1.96 * se
  
  tibble::tibble(
    estimate = estimate,
    se = se,
    low = low,
    upp = upp
  )
}

# Compute design-weighted predicted prevalence by subgroup. This is used
# for the residence and region scenario plots.
predict_prevalence_by <- function(d_new, model, group_var) {
  p <- safe_predict_prob(model, d_new)
  
  des <- make_design(d_new)
  des$variables$p <- p
  des$variables$grp <- droplevels(as.factor(d_new[[group_var]]))
  
  by_obj <- survey::svyby(
    ~p,
    ~grp,
    design = des,
    FUN = survey::svymean,
    na.rm = TRUE
  )
  
  estimate <- as.numeric(by_obj$p)
  se <- as.numeric(survey::SE(by_obj))
  low <- estimate - 1.96 * se
  upp <- estimate + 1.96 * se
  
  tibble::tibble(
    group = by_obj$grp,
    estimate = estimate,
    se = se,
    low = low,
    upp = upp
  )
}

# Identify the weighted modal contraceptive category excluding "not using"
# and missing. This provides the destination category for the uptake shift
# scenario.
modal_non_no_contra_w <- function(v, w, no_label = "not using") {
  v_chr <- as.character(v)
  
  is_no <- tolower(v_chr) == tolower(no_label)
  is_missing <- tolower(v_chr) == "missing" | is.na(v_chr)
  
  keep <- !(is_no | is_missing)
  if (!any(keep)) return(NA_character_)
  
  tmp <- tibble::tibble(cat = v_chr[keep], w = w[keep]) %>%
    dplyr::group_by(cat) %>%
    dplyr::summarise(w_sum = sum(w, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(w_sum))
  
  tmp$cat[1]
}

# Shift a specified weighted proportion of women from "not using" to the
# weighted modal non-"not using" contraceptive category. 
shift_contraception_weighted <- function(v, w, increase_pp = 0.15, no_label = "not using") {
  stopifnot(length(v) == length(w))
  
  v_chr <- as.character(v)
  target <- modal_non_no_contra_w(v, w, no_label = no_label)
  
  if (is.na(target)) return(v)
  
  is_no <- tolower(v_chr) == tolower(no_label)
  idx_no <- which(is_no)
  
  if (length(idx_no) == 0) return(v)
  
  w_total <- sum(w, na.rm = TRUE)
  w_move_target <- increase_pp * w_total
  
  w_no <- sum(w[idx_no], na.rm = TRUE)
  w_move_target <- min(w_move_target, w_no)
  
  idx_sorted <- idx_no[order(w[idx_no], decreasing = TRUE)]
  cum_w <- cumsum(w[idx_sorted])
  
  move_idx <- idx_sorted[cum_w <= w_move_target]
  if (length(move_idx) == 0 && length(idx_sorted) > 0) {
    move_idx <- idx_sorted[1]
  }
  
  v2_chr <- v_chr
  v2_chr[move_idx] <- target
  
  factor(v2_chr, levels = levels(v))
}

# Helper to plot subgroup scenario results in a consistent style.
plot_subgroup <- function(tbl, var_name, file_stub, fig_dir) {
  df_plot <- tbl %>%
    filter(var == var_name) %>%
    mutate(
      scenario = dplyr::recode(
        scenario,
        "baseline" = "Baseline",
        "S1_educ_ctx_plus20pp" = "S1: +20pp contextual secondary+",
        "S2_contra_plus15pp" = "S2: +15pp contraceptive use",
        "S3_both" = "S3: Both"
      ),
      scenario = factor(
        scenario,
        levels = c(
          "Baseline",
          "S1: +20pp contextual secondary+",
          "S2: +15pp contraceptive use",
          "S3: Both"
        )
      ),
      group = droplevels(as.factor(group))
    )
  
  p <- ggplot(df_plot, aes(x = group, y = estimate, fill = scenario)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_errorbar(
      aes(ymin = low, ymax = upp),
      position = position_dodge(width = 0.8),
      width = 0.18,
      linewidth = 0.5
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = paste0("Predicted HRFB prevalence by ", str_to_title(var_name), " under counterfactual scenarios"),
      subtitle = "Design-adjusted g-computation from Stage 4 model",
      x = str_to_title(var_name),
      y = "Predicted prevalence (survey-weighted)",
      fill = "Scenario",
      caption = caption_rq_policy
    ) +
    theme_dissertation() +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = margin(10, 12, 18, 12),
      plot.caption = element_text(margin = margin(t = 8))
    )
  
  print(p)
  save_fig(
    p,
    paste0("stage6_by_", file_stub, "_scenario_prevalence"),
    fig_dir = fig_dir,
    w = 9.2,
    h = 5.8
  )
}

# --------------------------
# 1) Inputs and outputs
# --------------------------

# Stage 6 uses the Stage 2 enriched dataset as its working input.
# The survey design is rebuilt directly from this file so that scenario
# predictions remain bound to the current analysis dataset.
in_births <- "outputs/step02_births_qc_enhanced.rds"
stopifnot(file.exists(in_births))

df <- readRDS(in_births)

out_dir <- "outputs"
tab_dir <- file.path(out_dir, "tables_stage6")
fig_dir <- file.path(out_dir, "figs_stage6")

if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# --------------------------
# 2) Guardrails and factor housekeeping
# --------------------------

# These variables are required both for refitting the Stage 4-style model
# and for computing design-adjusted scenario predictions.
raw_required <- c(
  "hrfb_binary", "hrfb_multicat",
  "residence", "region", "educ", "wealth", "ethnicity",
  "fp_radio", "fp_tv", "fp_print",
  "contraceptive_use", "fert_pref", "share_sec_plus",
  "v021", "strata", "weight"
)

miss <- setdiff(raw_required, names(df))
if (length(miss) > 0) {
  stop("Missing required Stage 6 variables: ", paste(miss, collapse = ", "))
}

# Factor coding is aligned to Stage 4 so that the refitted logistic model
# uses the same reference categories and predictor interpretation.
df <- df %>%
  mutate(
    residence = fct_relevel(residence, "urban"),
    educ = fct_relevel(educ, "higher", "secondary", "primary", "no education"),
    wealth = fct_relevel(wealth, "richest", "richer", "middle", "poorer", "poorest"),
    region = factor(stringr::str_to_title(as.character(region))),
    ethnicity = factor(stringr::str_to_title(as.character(ethnicity))),
    ethnicity = fct_drop(ethnicity),
    fp_radio = fct_explicit_na(fp_radio, na_level = "Missing"),
    fp_tv = fct_explicit_na(fp_tv, na_level = "Missing"),
    fp_print = fct_explicit_na(fp_print, na_level = "Missing"),
    contraceptive_use = fct_explicit_na(contraceptive_use, na_level = "Missing"),
    fert_pref = fct_explicit_na(fert_pref, na_level = "Missing")
  )

# --------------------------
# 3) Refit Stage 4 primary model
# --------------------------

# The Stage 6 scenario analysis is built on the same model structure as
# Stage 4 so that counterfactual predictions remain aligned with the main model.
f_bin <- hrfb_binary ~
  residence + region + educ + wealth + ethnicity +
  fp_radio + fp_tv + fp_print +
  contraceptive_use + fert_pref +
  share_sec_plus

des0 <- make_design(df)

m_bin <- svyglm(
  formula = f_bin,
  design = des0,
  family = quasibinomial()
)

# The refitted coefficient table is exported.
tidy_bin <- broom::tidy(
  m_bin,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
) %>%
  rename(OR = estimate, LCL = conf.low, UCL = conf.high) %>%
  select(term, OR, LCL, UCL, std.error, statistic, p.value)

readr::write_csv(tidy_bin, file.path(tab_dir, "stage6_refit_primary_logistic_OR.csv"))
print(head(tidy_bin, 12))

# --------------------------
# 4) Define counterfactual scenarios
# --------------------------

scenarios <- list(
  baseline = function(d) {
    d
  },
  
  S1_educ_ctx_plus20pp = function(d) {
    d %>%
      mutate(share_sec_plus = cap01(share_sec_plus + 0.20))
  },
  
  S2_contra_plus15pp = function(d) {
    d %>%
      mutate(
        contraceptive_use = shift_contraception_weighted(
          contraceptive_use,
          w = weight,
          increase_pp = 0.15
        )
      )
  },
  
  S3_both = function(d) {
    d %>%
      mutate(
        share_sec_plus = cap01(share_sec_plus + 0.20),
        contraceptive_use = shift_contraception_weighted(
          contraceptive_use,
          w = weight,
          increase_pp = 0.15
        )
      )
  }
)

# A diagnostic message is printed to show which contraceptive
# category receives the weighted shift in Scenario 2.
target_method <- modal_non_no_contra_w(
  df$contraceptive_use,
  df$weight,
  no_label = "not using"
)

message("[Stage 6] Scenario 2 target contraceptive category: ", target_method)

# --------------------------
# 5) Overall scenario evaluation via g-computation
# --------------------------

# Each scenario is applied to the full dataset, predicted probabilities are
# generated from the refitted model, and the design-weighted mean predicted
# prevalence is then estimated.
scn_results_raw <- purrr::imap_dfr(
  scenarios,
  ~{
    d_cf <- .x(df)
    res <- predict_prevalence(d_cf, m_bin)
    res$scenario <- .y
    res
  }
)

# Reshape to compute absolute and relative differences from baseline.
scn_results_wide <- scn_results_raw %>%
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = c(estimate, se, low, upp)
  ) %>%
  mutate(
    abs_change_S1 = estimate_S1_educ_ctx_plus20pp - estimate_baseline,
    abs_change_S2 = estimate_S2_contra_plus15pp - estimate_baseline,
    abs_change_S3 = estimate_S3_both - estimate_baseline,
    rel_change_S1 = abs_change_S1 / estimate_baseline,
    rel_change_S2 = abs_change_S2 / estimate_baseline,
    rel_change_S3 = abs_change_S3 / estimate_baseline
  )

# Build a tidy long table for output and plotting.
scn_results_long <- tibble::tibble(
  scenario = c(
    "Baseline",
    "S1: +20pp contextual secondary+",
    "S2: +15pp contraceptive use",
    "S3: Both"
  ),
  estimate = c(
    scn_results_wide$estimate_baseline,
    scn_results_wide$estimate_S1_educ_ctx_plus20pp,
    scn_results_wide$estimate_S2_contra_plus15pp,
    scn_results_wide$estimate_S3_both
  ),
  se = c(
    scn_results_wide$se_baseline,
    scn_results_wide$se_S1_educ_ctx_plus20pp,
    scn_results_wide$se_S2_contra_plus15pp,
    scn_results_wide$se_S3_both
  ),
  low = c(
    scn_results_wide$low_baseline,
    scn_results_wide$low_S1_educ_ctx_plus20pp,
    scn_results_wide$low_S2_contra_plus15pp,
    scn_results_wide$low_S3_both
  ),
  upp = c(
    scn_results_wide$upp_baseline,
    scn_results_wide$upp_S1_educ_ctx_plus20pp,
    scn_results_wide$upp_S2_contra_plus15pp,
    scn_results_wide$upp_S3_both
  )
) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "Baseline",
        "S1: +20pp contextual secondary+",
        "S2: +15pp contraceptive use",
        "S3: Both"
      )
    )
  )

readr::write_csv(scn_results_long, file.path(tab_dir, "stage6_scenario_prevalence_overall.csv"))
readr::write_csv(scn_results_wide, file.path(tab_dir, "stage6_scenario_prevalence_changes.csv"))
print(scn_results_long)

# This diagnostic confirms whether Scenario 2 achieves approximately the
# intended weighted shift away from "not using."
wmean <- function(x, w) {
  sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

is_not_using <- function(v) {
  tolower(as.character(v)) == "not using"
}

base_not <- wmean(is_not_using(df$contraceptive_use), df$weight)
s2_df <- scenarios$S2_contra_plus15pp(df)
s2_not <- wmean(is_not_using(s2_df$contraceptive_use), s2_df$weight)

message(sprintf(
  "[Stage 6] Weighted 'not using' prevalence: baseline = %.3f; S2 = %.3f; change = %.3f",
  base_not, s2_not, (s2_not - base_not)
))

message(sprintf(
  "[Stage 6] Implied weighted uptake increase (percentage points): %.1f",
  100 * (base_not - s2_not)
))

# --------------------------
# 6) Subgroup scenario evaluation
# --------------------------

# Subgroup estimates are produced for residence and region
subgroup_vars <- c("residence", "region")

subgroup_tables <- purrr::imap(
  scenarios,
  ~{
    d_cf <- .x(df)
    
    purrr::map(subgroup_vars, function(gv) {
      predict_prevalence_by(d_cf, m_bin, gv) %>%
        mutate(scenario = .y, var = gv)
    }) %>%
      dplyr::bind_rows()
  }
) %>%
  dplyr::bind_rows()

readr::write_csv(subgroup_tables, file.path(tab_dir, "stage6_scenario_prevalence_by_group.csv"))
print(head(subgroup_tables, 12))

# --------------------------
# 7) Plots
# --------------------------

# Overall scenario prevalence plot.
p_overall <- ggplot(scn_results_long, aes(x = scenario, y = estimate)) +
  geom_col(width = 0.75) +
  geom_errorbar(aes(ymin = low, ymax = upp), width = 0.15, linewidth = 0.5) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Predicted prevalence of any HRFB under counterfactual scenarios",
    subtitle = "Design-adjusted g-computation from the Stage 4 primary model",
    x = NULL,
    y = "Predicted prevalence (survey-weighted)",
    caption = caption_rq_policy
  ) +
  theme_dissertation() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.margin = margin(10, 12, 18, 12),
    plot.caption = element_text(margin = margin(t = 8))
  )

print(p_overall)
save_fig(
  p_overall,
  "stage6_overall_scenario_prevalence",
  fig_dir = fig_dir,
  w = 8.8,
  h = 5.6
)

# Residence and region scenario plots.
plot_subgroup(subgroup_tables, "residence", "residence", fig_dir = fig_dir)
plot_subgroup(subgroup_tables, "region", "region", fig_dir = fig_dir)

# --------------------------
# 8) Write Stage 6 README
# --------------------------

# The README documents the scenario definitions, model basis, and key
# files written by this stage.
readme_path <- file.path(out_dir, "README_stage6_scenarios.txt")

readme_txt <- c(
  "STAGE 6 — Scenario modelling (exploratory; policy-focused)",
  "",
  "Inputs:",
  "- outputs/step02_births_qc_enhanced.rds (working dataset: raw + QA variables)",
  "",
  "Method summary:",
  "- Refit the Stage 4 primary survey-adjusted logistic model (svyglm with quasibinomial family).",
  "- Generate counterfactual predictions using prediction-based g-computation.",
  "- Average predicted probabilities using the survey design to obtain design-adjusted prevalence under each scenario.",
  "",
  "Scenarios:",
  "- Baseline: observed data with no intervention.",
  "- S1 (+20pp contextual education): share_sec_plus := pmin(1, share_sec_plus + 0.20).",
  "- S2 (+15pp contraceptive use): reassign approximately 15 percentage points of the weighted population from 'not using' to the weighted modal non-'not using' contraceptive category.",
  "- S3 (Both): apply S1 and S2 simultaneously.",
  "",
  "Outputs (tables):",
  paste0("- ", file.path("tables_stage6", "stage6_refit_primary_logistic_OR.csv")),
  paste0("- ", file.path("tables_stage6", "stage6_scenario_prevalence_overall.csv")),
  paste0("- ", file.path("tables_stage6", "stage6_scenario_prevalence_changes.csv")),
  paste0("- ", file.path("tables_stage6", "stage6_scenario_prevalence_by_group.csv")),
  "",
  "Outputs (figures):",
  "- PDF and PNG versions are written for all main plotted outputs in figs_stage6/",
  paste0("- ", file.path("figs_stage6", "stage6_overall_scenario_prevalence.*")),
  paste0("- ", file.path("figs_stage6", "stage6_by_residence_scenario_prevalence.*")),
  paste0("- ", file.path("figs_stage6", "stage6_by_region_scenario_prevalence.*")),
  "",
  "Notes:",
  "- The Stage 6 model specification matches Stage 4 to maintain consistency and avoid definitional leakage.",
  "- These scenario estimates are exploratory model-based counterfactual predictions rather than identified causal effects.",
  "- Predictions are averaged using survey weights and the complex survey design.",
  "",
  paste0("Saved on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)

writeLines(readme_txt, readme_path)

message("Wrote README: ", normalizePath(readme_path))
message("Stage 6 complete. Outputs written to: ", normalizePath(out_dir))
############################################################