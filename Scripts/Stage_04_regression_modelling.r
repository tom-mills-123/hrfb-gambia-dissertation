############################################################
# STAGE 4 — Regression modelling (RQ2)
# Project: High-risk fertility behaviours (HRFB) — Gambia DHS 2019–20
#
# What this script does:
#   1) Load the Stage 2 enriched working dataset and rebuild the
#      survey design on this dataset.
#   2) Enforce a raw-outcome policy for the primary regression models.
#   3) Fit the primary survey-adjusted binary logistic regression model
#      for any HRFB using socio-demographic, FP exposure, reproductive,
#      and contextual predictors.
#   4) Export odds-ratio tables, model comparison summaries, and
#      publication-ready forest plots.
#   5) Produce adjusted prevalence plots from the primary binary model.
#   6) Explore an optional survey-adjusted multinomial model for
#      none / single / multiple HRFB.
#   7) Run sensitivity checks using age-stratified models and a
#      one-birth-per-woman restriction.
#   8) Write a Stage 4 README documenting model choices and outputs.
#
# Key principle for this stage:
#   The primary regression analysis uses the Stage 2 dataset as the
#   working container, but the main outcome remains the raw Stage 1
#   binary HRFB indicator. Variables that directly define HRFB
#   (maternal age at birth, short birth interval, parity) are not
#   included as predictors in the primary model in order to avoid
#   circularity between predictors and outcome.
############################################################

# --------------------------
# 0) Packages + global options
# --------------------------
AUTO_INSTALL <- FALSE

# These packages support data handling, survey-weighted modelling,
# table export, and figure production for the regression stage.
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
  "VGAM",
  "svyVGAM",
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
  library(VGAM)
  library(svyVGAM)
  library(purrr)
  library(Cairo)
})

# This matches the survey option used in earlier stages.
options(survey.lonely.psu = "adjust")

# --------------------------
# Helper functions
# --------------------------

# Define a consistent plotting theme used across all dissertation figures.
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

# Save figures in both PDF and PNG formats
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

# This helper standardises predictions from a fitted binary survey model
# by rebuilding the model matrix carefully and returning probabilities on
# the response scale.
predict_prob_safe <- function(model, newdata) {
  tt <- stats::delete.response(stats::terms(model))
  
  mf <- stats::model.frame(
    tt,
    data = newdata,
    xlev = model$xlevels,
    na.action = stats::na.pass
  )
  
  mm <- stats::model.matrix(
    tt,
    data = mf,
    contrasts.arg = model$contrasts
  )
  
  beta <- stats::coef(model)
  
  X <- matrix(
    0,
    nrow = nrow(mm),
    ncol = length(beta),
    dimnames = list(NULL, names(beta))
  )
  
  common <- intersect(colnames(mm), names(beta))
  if (length(common) == 0) {
    stop(
      "No overlapping columns between model matrix and coefficients. ",
      "This suggests a formula or data mismatch."
    )
  }
  
  X[, common] <- mm[, common, drop = FALSE]
  eta <- as.vector(X %*% beta)
  
  model$family$linkinv(eta)
}

# Design-based standardised prevalence for a categorical predictor.
# Each category is set in turn while the remaining covariates are left
# at their observed distribution.
std_prev_by_cat <- function(model, var, label_map = NULL) {
  des <- model$survey.design
  dat <- des$variables
  
  if (!(var %in% names(dat))) stop("Variable not found in model design data: ", var)
  
  if (!is.null(model$xlevels) && length(model$xlevels) > 0) {
    for (v in names(model$xlevels)) {
      if (v %in% names(dat)) {
        dat[[v]] <- factor(dat[[v]], levels = model$xlevels[[v]])
      }
    }
  }
  
  if (!is.null(model$xlevels) && var %in% names(model$xlevels)) {
    levs <- model$xlevels[[var]]
  } else {
    dat[[var]] <- factor(dat[[var]])
    levs <- levels(dat[[var]])
  }
  
  res <- purrr::map_dfr(levs, function(L) {
    newdat <- dat
    newdat[[var]] <- factor(L, levels = levs)
    
    p <- predict_prob_safe(model, newdata = newdat)
    des_p <- update(des, .p = p)
    est <- survey::svymean(~.p, des_p, na.rm = TRUE)
    ci <- suppressMessages(stats::confint(est))
    
    tibble::tibble(
      var = var,
      level = L,
      prev = unname(coef(est)[1]),
      se = unname(SE(est)[1]),
      low = unname(ci[1]),
      upp = unname(ci[2])
    )
  })
  
  if (!is.null(label_map)) {
    res <- res %>% mutate(level = dplyr::recode(level, !!!label_map))
  }
  
  res
}

# Design-based standardised prevalence for a continuous predictor at a
# specified grid of values.
std_prev_by_cont_grid <- function(model, var, grid_vals) {
  des <- model$survey.design
  dat <- des$variables
  
  if (!(var %in% names(dat))) stop("Variable not found in model design data: ", var)
  
  res <- purrr::map_dfr(grid_vals, function(g) {
    newdat <- dat
    
    if (!is.null(model$xlevels) && length(model$xlevels) > 0) {
      for (v in names(model$xlevels)) {
        if (v %in% names(newdat)) {
          newdat[[v]] <- factor(newdat[[v]], levels = model$xlevels[[v]])
        }
      }
    }
    
    newdat[[var]] <- g
    
    p <- predict_prob_safe(model, newdata = newdat)
    des_p <- update(des, .p = p)
    
    est <- survey::svymean(~.p, des_p, na.rm = TRUE)
    ci <- suppressMessages(stats::confint(est))
    
    tibble::tibble(
      var = var,
      value = g,
      prev = unname(coef(est)[1]),
      se = unname(SE(est)[1]),
      low = unname(ci[1]),
      upp = unname(ci[2])
    )
  })
  
  res
}

# A simple null-coalescing helper is used in the adjusted prevalence plot function.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Plot adjusted prevalence results in the same visual style used in Stage 3.
plot_std_prev <- function(df_prev, title, file_stub, fig_dir, level_order = NULL, xlab = NULL) {
  df_plot <- df_prev %>%
    mutate(level = as.character(level)) %>%
    filter(!is.na(level))
  
  if (!is.null(level_order)) {
    df_plot <- df_plot %>%
      mutate(level = factor(level, levels = level_order)) %>%
      arrange(level)
  } else {
    df_plot <- df_plot %>%
      mutate(level = forcats::fct_inorder(level))
  }
  
  p <- ggplot(df_plot, aes(x = level, y = prev)) +
    geom_col() +
    geom_errorbar(aes(ymin = low, ymax = upp), width = 0.2) +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, 1),
      oob = scales::squish,
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = title,
      subtitle = "Adjusted predicted prevalence (standardised over the sample); 95% confidence intervals",
      x = xlab %||% NULL,
      y = "Adjusted prevalence of any HRFB",
      caption = "Data: The Gambia DHS 2019–20. Predictions from survey-adjusted logistic regression; 95% confidence intervals."
    ) +
    theme_dissertation()
  
  print(p)
  save_fig(p, file_stub, fig_dir = fig_dir, w = 8, h = 5)
}

# This helper assigns a stable conceptual order to coefficient groups
# so forest plots remain interpretable and consistent across runs.
term_block_order <- function(term) {
  case_when(
    str_detect(term, "^residence") ~ 1,
    str_detect(term, "^region") ~ 2,
    str_detect(term, "^educ") ~ 3,
    str_detect(term, "^wealth") ~ 4,
    str_detect(term, "^ethnicity") ~ 5,
    str_detect(term, "^fp_") ~ 6,
    str_detect(term, "^contraceptive_use") ~ 7,
    str_detect(term, "^fert_pref") ~ 8,
    term == "share_sec_plus" ~ 9,
    TRUE ~ 99
  )
}

# Label model terms for export tables and plots.
label_term <- function(term) {
  term %>%
    str_replace("^residencerural$", "Residence: Rural (ref: Urban)") %>%
    str_replace("^educsecondary$", "Education: Secondary (ref: Higher)") %>%
    str_replace("^educprimary$", "Education: Primary (ref: Higher)") %>%
    str_replace("^educno education$", "Education: No education (ref: Higher)") %>%
    str_replace("^wealthricher$", "Wealth: Richer (ref: Richest)") %>%
    str_replace("^wealthmiddle$", "Wealth: Middle (ref: Richest)") %>%
    str_replace("^wealthpoorer$", "Wealth: Poorer (ref: Richest)") %>%
    str_replace("^wealthpoorest$", "Wealth: Poorest (ref: Richest)") %>%
    str_replace("^fp_radio", "FP message exposure (radio): ") %>%
    str_replace("^fp_tv", "FP message exposure (TV): ") %>%
    str_replace("^fp_print", "FP message exposure (print): ") %>%
    str_replace("^contraceptive_use", "Contraceptive use: ") %>%
    str_replace("^fert_pref", "Fertility preference: ") %>%
    str_replace("^share_sec_plus$", "Cluster share with ≥ secondary education (per 1.0 increase)") %>%
    str_replace("^region", "Region: ") %>%
    str_replace("^ethnicity", "Ethnicity: ")
}

# A helper for fitting the multinomial model safely so that failures
# are trapped cleanly rather than stopping the whole script.
fit_mn_safe <- function(formula_obj, design_obj, label) {
  message("\n--- Trying ", label, " ---")
  
  out <- tryCatch(
    svyVGAM::svy_vglm(
      formula = formula_obj,
      design = design_obj,
      family = VGAM::multinomial(refLevel = "none")
    ),
    error = function(e) e
  )
  
  if (inherits(out, "error")) {
    message("FAILED: ", conditionMessage(out))
    return(NULL)
  } else {
    message("SUCCESS")
    return(out)
  }
}

# Manual tidying for svy_vglm output.
tidy_svy_vglm_manual <- function(model, outcome_levels = c("single", "multiple")) {
  cf <- stats::coef(model)
  V <- stats::vcov(model)
  
  est <- as.numeric(cf)
  se <- sqrt(diag(V))
  z <- est / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  low <- est - 1.96 * se
  upp <- est + 1.96 * se
  
  nm <- names(cf)
  parsed <- stringr::str_match(nm, "^(.*?):([0-9]+)$")
  
  term <- parsed[, 2]
  out_idx <- parsed[, 3]
  
  bad <- is.na(term) | is.na(out_idx)
  if (any(bad)) {
    term[bad] <- nm[bad]
  }
  
  outcome_level <- rep(NA_character_, length(nm))
  outcome_level[out_idx == "1"] <- outcome_levels[1]
  outcome_level[out_idx == "2"] <- outcome_levels[2]
  
  term[term == "(Intercept)"] <- "(Intercept)"
  
  tibble::tibble(
    outcome_level = outcome_level,
    term = term,
    estimate = est,
    conf.low = low,
    conf.high = upp,
    std.error = se,
    statistic = z,
    p.value = p
  ) %>%
    mutate(
      RRR = exp(estimate),
      LCL = exp(conf.low),
      UCL = exp(conf.high)
    ) %>%
    select(outcome_level, term, RRR, LCL, UCL, std.error, statistic, p.value)
}

# Label multinomial terms for clearer export tables and plots.
label_term_mn <- function(term) {
  term %>%
    str_replace("^residencerural$", "Residence: Rural (ref: Urban)") %>%
    str_replace("^educsecondary$", "Education: Secondary (ref: Higher)") %>%
    str_replace("^educprimary$", "Education: Primary (ref: Higher)") %>%
    str_replace("^educno education$", "Education: No education (ref: Higher)") %>%
    str_replace("^wealthricher$", "Wealth: Richer (ref: Richest)") %>%
    str_replace("^wealthmiddle$", "Wealth: Middle (ref: Richest)") %>%
    str_replace("^wealthpoorer$", "Wealth: Poorer (ref: Richest)") %>%
    str_replace("^wealthpoorest$", "Wealth: Poorest (ref: Richest)") %>%
    str_replace("^contraceptive_use_mnmodern method$", "Contraceptive use: Modern method (ref: Not using)") %>%
    str_replace("^contraceptive_use_mntraditional/other$", "Contraceptive use: Traditional/other (ref: Not using)") %>%
    str_replace("^contraceptive_use_mnMissing$", "Contraceptive use: Missing (ref: Not using)") %>%
    str_replace("^fert_pref_mnundecided$", "Fertility preference: Undecided (ref: Have another)") %>%
    str_replace("^fert_pref_mnno more$", "Fertility preference: No more (ref: Have another)") %>%
    str_replace("^fert_pref_mnother/infecund/missing$", "Fertility preference: Other/infecund/missing (ref: Have another)") %>%
    str_replace("^share_sec_plus$", "Cluster share with ≥ secondary education") %>%
    str_replace("^region", "Region: ") %>%
    str_replace("^ethnicity", "Ethnicity: ") %>%
    str_replace("^fp_radio", "FP message exposure (radio): ") %>%
    str_replace("^fp_tv", "FP message exposure (TV): ") %>%
    str_replace("^fp_print", "FP message exposure (print): ")
}

term_block_order_mn <- function(term) {
  dplyr::case_when(
    stringr::str_detect(term, "^residence") ~ 1,
    stringr::str_detect(term, "^region") ~ 2,
    stringr::str_detect(term, "^educ") ~ 3,
    stringr::str_detect(term, "^wealth") ~ 4,
    stringr::str_detect(term, "^ethnicity") ~ 5,
    stringr::str_detect(term, "^fp_") ~ 6,
    stringr::str_detect(term, "^contraceptive_use_mn") ~ 7,
    stringr::str_detect(term, "^fert_pref_mn") ~ 8,
    term == "share_sec_plus" ~ 9,
    TRUE ~ 99
  )
}

# Export odds-ratio tables from fitted stratified binary models.
export_or <- function(mod, out_path) {
  if (inherits(mod, "try-error")) return(invisible(NULL))
  
  tb <- broom::tidy(
    mod,
    conf.int = TRUE,
    conf.level = 0.95,
    exponentiate = TRUE
  ) %>%
    rename(OR = estimate, LCL = conf.low, UCL = conf.high) %>%
    select(term, OR, LCL, UCL, std.error, statistic, p.value)
  
  readr::write_csv(tb, out_path)
  tb
}

# --------------------------
# 1) Inputs and outputs
# --------------------------

# Stage 4 uses the Stage 2 enriched dataset as its working input.
# This preserves the QA trail while keeping the primary regression
# analysis anchored to the raw Stage 1 outcome definition.
in_births <- "outputs/step02_births_qc_enhanced.rds"
stopifnot(file.exists(in_births))

df <- readRDS(in_births)

# Output directories for Stage 4 tables, figures, and README files.
out_dir <- "outputs"
tab_dir <- file.path(out_dir, "tables_stage4")
fig_dir <- file.path(out_dir, "figs_stage4")

if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# --------------------------
# 2) Guardrails and factor housekeeping
# --------------------------

# These variables are required for the main modelling stage.
raw_required <- c(
  "hrfb_binary", "hrfb_multicat",
  "residence", "region", "educ", "wealth", "ethnicity",
  "fp_radio", "fp_tv", "fp_print",
  "contraceptive_use", "fert_pref",
  "share_sec_plus",
  "woman_id", "child_age_months",
  "weight", "strata", "v021",
  "age_at_birth_years"
)

miss <- setdiff(raw_required, names(df))
if (length(miss) > 0) {
  stop("Missing required variables: ", paste(miss, collapse = ", "))
}

# Factor levels are standardised here so that model reference categories
# are explicit and interpretation remains stable across runs.
df <- df %>%
  mutate(
    residence = fct_relevel(residence, "urban"),
    educ = fct_relevel(educ, "higher", "secondary", "primary", "no education"),
    wealth = fct_relevel(wealth, "richest", "richer", "middle", "poorer", "poorest"),
    ethnicity = fct_drop(ethnicity),
    region = factor(stringr::str_to_title(as.character(region))),
    ethnicity = factor(stringr::str_to_title(as.character(ethnicity))),
    fp_radio = fct_explicit_na(fp_radio, na_level = "Missing"),
    fp_tv = fct_explicit_na(fp_tv, na_level = "Missing"),
    fp_print = fct_explicit_na(fp_print, na_level = "Missing"),
    contraceptive_use = fct_explicit_na(contraceptive_use, na_level = "Missing"),
    fert_pref = fct_explicit_na(fert_pref, na_level = "Missing"),
    hrfb_multicat = fct_relevel(hrfb_multicat, "none", "single", "multiple")
  )

# FP exposure variables are forced to unordered factors so that the model
# uses category-specific contrasts rather than polynomial contrasts.
df <- df %>%
  mutate(
    fp_radio = factor(fp_radio, ordered = FALSE),
    fp_tv = factor(fp_tv, ordered = FALSE),
    fp_print = factor(fp_print, ordered = FALSE)
  )

# Explicit level ordering is applied to keep the reference category clear.
df <- df %>%
  mutate(
    fp_radio = factor(fp_radio, levels = c("not at all", "less than once a week", "at least once a week", "Missing")),
    fp_tv = factor(fp_tv, levels = c("not at all", "less than once a week", "at least once a week", "Missing")),
    fp_print = factor(fp_print, levels = c("not at all", "less than once a week", "at least once a week", "Missing"))
  )

# --------------------------
# 3) Build the survey design used for modelling
# --------------------------

# The Stage 4 design object is rebuilt directly on the current working
# dataset so that all model estimates are aligned to the same analysis frame.
df_design <- df %>%
  dplyr::filter(!is.na(weight), !is.na(v021), !is.na(strata))

svy_base <- survey::svydesign(
  ids = ~v021,
  strata = ~strata,
  weights = ~weight,
  data = df_design,
  nest = TRUE
)

message(
  "Design built on df_design: ",
  nrow(df_design), " rows; ",
  dplyr::n_distinct(df_design$v021), " PSUs; ",
  dplyr::n_distinct(df_design$strata), " strata."
)

# --------------------------
# 4) Primary binary model — survey-adjusted logistic regression
# --------------------------

# The primary model predicts any HRFB using socio-demographic, FP exposure,
# reproductive preference, contraceptive use, and contextual predictors.
# Definitional HRFB components are deliberately excluded from the predictor set.
f_bin <- hrfb_binary ~
  residence + region + educ + wealth + ethnicity +
  fp_radio + fp_tv + fp_print +
  contraceptive_use + fert_pref +
  share_sec_plus

# Quasibinomial estimation is used because it is a practical and standard
# choice for survey-weighted logistic models under possible mild dispersion issues.
m_bin <- svyglm(
  formula = f_bin,
  design = svy_base,
  family = quasibinomial()
)

# Two comparison models are fitted to support interpretation of overall fit.
m0 <- svyglm(
  hrfb_binary ~ 1,
  design = svy_base,
  family = quasibinomial()
)

m_core <- svyglm(
  hrfb_binary ~ residence + region + educ + wealth + ethnicity,
  design = svy_base,
  family = quasibinomial()
)

# --------------------------
# 4A) Model comparison summaries
# --------------------------

# McFadden-style pseudo-R2 is used here as a simple descriptive fit summary.
pseudo_r2_mcfadden <- function(m_full, m_null) {
  1 - (deviance(m_full) / deviance(m_null))
}

r2_core <- pseudo_r2_mcfadden(m_core, m0)
r2_full <- pseudo_r2_mcfadden(m_bin, m0)

pseudo_r2_tbl <- tibble::tibble(
  model = c("core", "full"),
  pseudo_r2 = c(r2_core, r2_full)
)

print(pseudo_r2_tbl)
readr::write_csv(pseudo_r2_tbl, file.path(tab_dir, "rq2_pseudoR2.csv"))

anova_core_full <- anova(m_core, m_bin)
print(anova_core_full)
capture.output(anova_core_full, file = file.path(tab_dir, "rq2_anova_core_vs_full.txt"))

# Joint Wald tests are used to assess the overall contribution of the
# main multi-level categorical predictors.
jt_region <- survey::regTermTest(m_bin, ~region)
jt_educ <- survey::regTermTest(m_bin, ~educ)
jt_wealth <- survey::regTermTest(m_bin, ~wealth)
jt_eth <- survey::regTermTest(m_bin, ~ethnicity)

print(jt_region)
print(jt_educ)
print(jt_wealth)
print(jt_eth)

joint_tests <- tibble::tibble(
  term = c("region", "educ", "wealth", "ethnicity"),
  F = c(
    unname(jt_region$Ftest[1]),
    unname(jt_educ$Ftest[1]),
    unname(jt_wealth$Ftest[1]),
    unname(jt_eth$Ftest[1])
  ),
  df_num = c(
    unname(jt_region$df[1]),
    unname(jt_educ$df[1]),
    unname(jt_wealth$df[1]),
    unname(jt_eth$df[1])
  ),
  df_den = c(
    unname(jt_region$df[2]),
    unname(jt_educ$df[2]),
    unname(jt_wealth$df[2]),
    unname(jt_eth$df[2])
  ),
  p_value = c(jt_region$p, jt_educ$p, jt_wealth$p, jt_eth$p)
)

readr::write_csv(joint_tests, file.path(tab_dir, "rq2_joint_wald_tests.csv"))

# --------------------------
# 4B) Residual and functional-form diagnostics
# --------------------------

# These residual plots are informal checks rather than formal tests,
# but they are useful for spotting obvious outliers or patterns.
df_diag <- svy_base$variables %>%
  dplyr::mutate(
    fitted = stats::fitted(m_bin),
    resid_pearson = stats::residuals(m_bin, type = "pearson"),
    resid_deviance = stats::residuals(m_bin, type = "deviance")
  )

p_resid1 <- ggplot(df_diag, aes(x = fitted, y = resid_pearson)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    title = "Pearson residuals vs fitted values",
    x = "Fitted probability",
    y = "Pearson residual"
  ) +
  theme_dissertation()

print(p_resid1)
save_fig(p_resid1, "rq2_resid_pearson_vs_fitted", fig_dir = fig_dir, w = 8, h = 5)

p_resid2 <- ggplot(df_diag, aes(x = fitted, y = resid_deviance)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    title = "Deviance residuals vs fitted values",
    x = "Fitted probability",
    y = "Deviance residual"
  ) +
  theme_dissertation()

print(p_resid2)
save_fig(p_resid2, "rq2_resid_deviance_vs_fitted", fig_dir = fig_dir, w = 8, h = 5)

p_resid3 <- ggplot(df_diag, aes(x = share_sec_plus, y = resid_pearson)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    title = "Pearson residuals vs share_sec_plus",
    x = "share_sec_plus",
    y = "Pearson residual"
  ) +
  theme_dissertation()

print(p_resid3)
save_fig(p_resid3, "rq2_resid_pearson_vs_share_sec_plus", fig_dir = fig_dir, w = 8, h = 5)

# A quadratic version of share_sec_plus is tested as a simple check on whether
# the continuous contextual effect is approximately linear on the logit scale.
f_bin_quad <- update(f_bin, . ~ . + I(share_sec_plus^2))

m_bin_quad <- svyglm(
  formula = f_bin_quad,
  design = svy_base,
  family = quasibinomial()
)

jt_quad <- survey::regTermTest(m_bin_quad, ~I(share_sec_plus^2))
print(jt_quad)

cmp_lin_quad <- anova(m_bin, m_bin_quad)
print(cmp_lin_quad)

capture.output(jt_quad, file = file.path(tab_dir, "rq2_share_sec_plus_quadratic_test.txt"))
capture.output(cmp_lin_quad, file = file.path(tab_dir, "rq2_anova_linear_vs_quadratic.txt"))

# --------------------------
# 4C) Primary odds-ratio table and forest plots
# --------------------------

# The main export table reports odds ratios with 95% confidence intervals
# and p-values from the primary survey-adjusted binary model.
tidy_bin <- broom::tidy(
  m_bin,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
) %>%
  rename(OR = estimate, LCL = conf.low, UCL = conf.high) %>%
  select(term, OR, LCL, UCL, std.error, statistic, p.value)

readr::write_csv(tidy_bin, file.path(tab_dir, "stage4_primary_logistic_OR.csv"))
print(head(tidy_bin, 12))

# A focused appendix table is also exported for the FP exposure terms.
fp_terms <- tidy_bin %>%
  dplyr::filter(stringr::str_detect(term, "^fp_radio|^fp_tv|^fp_print")) %>%
  dplyr::mutate(term_label = label_term(term)) %>%
  dplyr::select(term_label, OR, LCL, UCL, p.value)

readr::write_csv(fp_terms, file.path(tab_dir, "appendix_fp_exposure_terms.csv"))
print(fp_terms)

caption_rq2 <- "Data: The Gambia DHS 2019–20. Survey-weighted logistic regression; 95% confidence intervals via Taylor linearisation."

forest_df <- tidy_bin %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_label = label_term(term),
    block = term_block_order(term)
  ) %>%
  arrange(block, term_label) %>%
  mutate(term_label = factor(term_label, levels = rev(unique(term_label))))

p_forest <- ggplot(forest_df, aes(x = term_label, y = OR)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.15, linewidth = 0.5) +
  geom_point(size = 1.8) +
  coord_flip() +
  scale_y_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    labels = scales::number_format(accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Adjusted associations with any high-risk fertility behaviour (RQ2)",
    subtitle = "Odds ratios from survey-weighted logistic regression",
    x = NULL,
    y = "Odds ratio (log scale)",
    caption = caption_rq2
  ) +
  theme_dissertation(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.margin = margin(10, 12, 10, 10)
  )

print(p_forest)
save_fig(p_forest, "stage4_primary_logistic_forest", fig_dir = fig_dir, w = 8.5, h = 8.5)

# A smaller forest plot is also produced for fertility preference and
# contraceptive-use terms, which may be easier to read in the write-up.
forest_behav_df <- tidy_bin %>%
  filter(term != "(Intercept)") %>%
  filter(str_detect(term, "^fert_pref") | str_detect(term, "^contraceptive_use")) %>%
  mutate(
    term_label = label_term(term),
    block = case_when(
      str_detect(term, "^fert_pref") ~ 1,
      str_detect(term, "^contraceptive_use") ~ 2,
      TRUE ~ 99
    )
  ) %>%
  arrange(block, term_label) %>%
  mutate(term_label = factor(term_label, levels = rev(unique(term_label))))

p_forest_behav <- ggplot(forest_behav_df, aes(x = term_label, y = OR)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.15, linewidth = 0.5) +
  geom_point(size = 1.8) +
  coord_flip() +
  scale_y_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
    labels = scales::number_format(accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Adjusted associations with any HRFB: reproductive intentions and contraception",
    subtitle = "Subset of full model coefficients",
    x = NULL,
    y = "Odds ratio (log scale)",
    caption = caption_rq2
  ) +
  theme_dissertation(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.margin = margin(10, 12, 10, 10)
  )

print(p_forest_behav)
save_fig(p_forest_behav, "rq2_forest_fertpref_contraception", fig_dir = fig_dir, w = 8.5, h = 6.0)

# --------------------------
# 4D) Adjusted prevalence from the primary binary model
# --------------------------

# These adjusted prevalence plots complement the odds-ratio results by
# showing predicted prevalence on the outcome scale for key predictors.
ap_residence <- std_prev_by_cat(model = m_bin, var = "residence")

plot_std_prev(
  ap_residence,
  title = "Adjusted prevalence of any HRFB by residence",
  file_stub = "rq2_adjprev_by_residence",
  fig_dir = fig_dir,
  level_order = c("urban", "rural")
)

ap_educ <- std_prev_by_cat(model = m_bin, var = "educ")

plot_std_prev(
  ap_educ,
  title = "Adjusted prevalence of any HRFB by maternal education",
  file_stub = "rq2_adjprev_by_education",
  fig_dir = fig_dir,
  level_order = c("higher", "secondary", "primary", "no education")
)

ap_wealth <- std_prev_by_cat(model = m_bin, var = "wealth")

plot_std_prev(
  ap_wealth,
  title = "Adjusted prevalence of any HRFB by household wealth quintile",
  file_stub = "rq2_adjprev_by_wealth",
  fig_dir = fig_dir,
  level_order = c("richest", "richer", "middle", "poorer", "poorest")
)

ap_region <- std_prev_by_cat(model = m_bin, var = "region")

plot_std_prev(
  ap_region,
  title = "Adjusted prevalence of any HRFB by region",
  file_stub = "rq2_adjprev_by_region",
  fig_dir = fig_dir
)

ap_ethnicity <- std_prev_by_cat(model = m_bin, var = "ethnicity")

plot_std_prev(
  ap_ethnicity,
  title = "Adjusted prevalence of any HRFB by ethnicity",
  file_stub = "rq2_adjprev_by_ethnicity",
  fig_dir = fig_dir
)

# A small grid of values is used to show adjusted prevalence across the
# continuous contextual predictor share_sec_plus.
des_dat <- m_bin$survey.design$variables
q <- stats::quantile(des_dat$share_sec_plus, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
grid_vals <- as.numeric(q)

ap_share <- std_prev_by_cont_grid(m_bin, "share_sec_plus", grid_vals)

p_share <- ggplot(ap_share, aes(x = value, y = prev)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.2) +
  scale_y_continuous(
    limits = c(0, 1),
    oob = scales::squish,
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Adjusted prevalence of any HRFB across cluster education context",
    subtitle = "Predicted prevalence at selected percentiles of share_sec_plus; 95% confidence intervals",
    x = "Cluster share with secondary education or higher (share_sec_plus)",
    y = "Adjusted prevalence of any HRFB",
    caption = "Data: The Gambia DHS 2019–20. Predictions from survey-adjusted logistic regression; 95% confidence intervals."
  ) +
  theme_dissertation()

print(p_share)
save_fig(p_share, "rq2_adjprev_by_share_sec_plus", fig_dir = fig_dir, w = 8, h = 5)

# Adjusted prevalence tables are exported for use in the write-up.
readr::write_csv(ap_residence, file.path(tab_dir, "rq2_adjprev_by_residence.csv"))
readr::write_csv(ap_educ, file.path(tab_dir, "rq2_adjprev_by_education.csv"))
readr::write_csv(ap_wealth, file.path(tab_dir, "rq2_adjprev_by_wealth.csv"))
readr::write_csv(ap_region, file.path(tab_dir, "rq2_adjprev_by_region.csv"))
readr::write_csv(ap_ethnicity, file.path(tab_dir, "rq2_adjprev_by_ethnicity.csv"))
readr::write_csv(ap_share, file.path(tab_dir, "rq2_adjprev_by_share_sec_plus.csv"))

# --------------------------
# 5) OPTIONAL — Multinomial model (exploratory extension of RQ2)
# --------------------------

# PURPOSE OF THIS SECTION
# --------------------------
# The primary analysis models HRFB as a binary outcome (any vs none).
# This section extends that by modelling a 3-category outcome:
#   - none
#   - single HRFB
#   - multiple HRFB
#
# This allows us to assess whether predictors operate differently
# for increasing levels of risk.
#
# IMPORTANT:
# - This section is exploratory rather than primary inference.
# - Multinomial models are more fragile than binary models under
#   survey weighting, particularly when some categories are sparse.
# - To reduce instability, this section uses multinomial-specific
#   recodes and fits models in stages.

message("\n============================================================")
message("SECTION 5 — MULTINOMIAL MODEL EXPLORATION")
message("============================================================")

# --------------------------
# 5a) Basic diagnostics of outcome + key predictors
# --------------------------

# Before fitting a multinomial model, the first step is to inspect the
# raw cross-tabulations. This helps identify sparse cells that may lead to:
#   - convergence failure
#   - non-full-rank model matrices
#   - unstable or extreme estimates
message("\nOutcome distribution:")
print(table(df_design$hrfb_multicat, useNA = "ifany"))

message("\nCross-tab: outcome by contraceptive_use")
print(table(df_design$hrfb_multicat, df_design$contraceptive_use, useNA = "ifany"))

message("\nCross-tab: outcome by fert_pref")
print(table(df_design$hrfb_multicat, df_design$fert_pref, useNA = "ifany"))

# --------------------------
# 5b) Multinomial-specific recoding to reduce sparsity
# --------------------------

# The original DHS contraceptive-use variable contains many detailed
# categories. That level of detail is useful descriptively, but in a
# multinomial regression it can create tiny cells and unstable estimates.
#
# To make the multinomial model more robust, detailed categories are
# collapsed into broader groups that remain substantively meaningful.
#
# IMPORTANT:
# - These recodes are used only for the multinomial model.
# - They do not replace the original variables used in the primary binary model.
df_mn <- df_design %>%
  mutate(
    
    # Contraceptive-use categories are collapsed into:
    #   - not using
    #   - modern method
    #   - traditional/other
    #   - Missing
    contraceptive_use_mn = case_when(
      contraceptive_use == "not using" ~ "not using",
      contraceptive_use %in% c(
        "pill", "iud", "injections", "male condom",
        "female sterilization", "implants/norplant",
        "lactational amenorrhea (lam)", "standard days method (sdm)"
      ) ~ "modern method",
      contraceptive_use %in% c(
        "periodic abstinence", "withdrawal", "other traditional"
      ) ~ "traditional/other",
      TRUE ~ "Missing"
    ),
    
    # Fertility-preference categories are also simplified to reduce
    # sparse cells while keeping the main substantive distinctions.
    fert_pref_mn = case_when(
      fert_pref %in% c("have another", "undecided", "no more") ~ as.character(fert_pref),
      TRUE ~ "other/infecund/missing"
    )
  ) %>%
  mutate(
    
    # Explicit factor levels are set so that reference categories remain clear.
    contraceptive_use_mn = factor(
      contraceptive_use_mn,
      levels = c("not using", "modern method", "traditional/other", "Missing")
    ),
    fert_pref_mn = factor(
      fert_pref_mn,
      levels = c("have another", "undecided", "no more", "other/infecund/missing")
    ),
    
    # The multinomial outcome is explicitly ordered so that "none" is
    # the base outcome for later relative risk ratio interpretation.
    hrfb_multicat = factor(hrfb_multicat, levels = c("none", "single", "multiple"))
  ) %>%
  droplevels()

# --------------------------
# 5c) Post-recode diagnostics
# --------------------------

# After collapsing categories, the same cross-tabs are inspected again.
# This is a critical check because the recoding should leave the model
# with fewer sparse combinations and a more stable design matrix.
message("\nOutcome by simplified contraceptive_use_mn")
print(table(df_mn$hrfb_multicat, df_mn$contraceptive_use_mn, useNA = "ifany"))

message("\nOutcome by simplified fert_pref_mn")
print(table(df_mn$hrfb_multicat, df_mn$fert_pref_mn, useNA = "ifany"))

# --------------------------
# 5d) Build survey design for multinomial modelling
# --------------------------

# The multinomial model is fit on the recoded dataset, so a separate
# survey design object is built here using the same weights, PSU, and strata.
svy_base_mn <- survey::svydesign(
  ids = ~v021,
  strata = ~strata,
  weights = ~weight,
  data = df_mn,
  nest = TRUE
)

# --------------------------
# 5e) Define staged model specifications
# --------------------------

# Multinomial models can fail for technical reasons once too many sparse
# predictors are added. To diagnose where instability appears, models are
# fitted in stages from simplest to most complex.

# Model A includes only socio-demographic predictors.
mn_formula_a <- hrfb_multicat ~
  residence + region + educ + wealth + ethnicity

# Model B adds fertility preference and the contextual cluster measure.
mn_formula_b <- hrfb_multicat ~
  residence + region + educ + wealth + ethnicity +
  fert_pref_mn + share_sec_plus

# Model C adds simplified contraceptive use.
mn_formula_c <- hrfb_multicat ~
  residence + region + educ + wealth + ethnicity +
  fert_pref_mn + contraceptive_use_mn + share_sec_plus

# Model D is the most complete specification and also adds FP media exposure.
mn_formula_d <- hrfb_multicat ~
  residence + region + educ + wealth + ethnicity +
  fp_radio + fp_tv + fp_print +
  fert_pref_mn + contraceptive_use_mn + share_sec_plus

# --------------------------
# 5f) Fit models sequentially using a safe wrapper
# --------------------------

# fit_mn_safe() catches estimation failures and allows the script to continue,
# rather than stopping entirely if one multinomial specification breaks down.
m_mn_a <- fit_mn_safe(mn_formula_a, svy_base_mn, "Model A: socio-demographic only")
m_mn_b <- fit_mn_safe(mn_formula_b, svy_base_mn, "Model B: + fertility preference + context")
m_mn_c <- fit_mn_safe(mn_formula_c, svy_base_mn, "Model C: + simplified contraceptive use")
m_mn_d <- fit_mn_safe(mn_formula_d, svy_base_mn, "Model D: + FP media exposures")

# --------------------------
# 5g) Select the most complete successful model
# --------------------------

# The final multinomial model is chosen hierarchically:
#   D > C > B > A
# In other words, the most complete model that converges successfully
# is kept for export and plotting.
m_mn_final <- NULL
m_mn_label <- NULL

if (!is.null(m_mn_d)) {
  m_mn_final <- m_mn_d
  m_mn_label <- "Model D"
} else if (!is.null(m_mn_c)) {
  m_mn_final <- m_mn_c
  m_mn_label <- "Model C"
} else if (!is.null(m_mn_b)) {
  m_mn_final <- m_mn_b
  m_mn_label <- "Model B"
} else if (!is.null(m_mn_a)) {
  m_mn_final <- m_mn_a
  m_mn_label <- "Model A"
}

# A clear message is printed so it is obvious whether the multinomial
# extension succeeded and, if so, which specification was retained.
if (is.null(m_mn_final)) {
  message("\nNo multinomial model could be fitted successfully.")
} else {
  message("\nSelected multinomial model for export: ", m_mn_label)
  print(summary(m_mn_final))
}

# --------------------------
# 5h) Extract and export multinomial results
# --------------------------

# svy_vglm output does not tidy neatly with broom, so a custom helper is
# used to extract coefficients and convert them into relative risk ratios.
if (!is.null(m_mn_final)) {
  tidy_mn <- tidy_svy_vglm_manual(m_mn_final)
  
  # The full multinomial results table is exported first.
  readr::write_csv(
    tidy_mn,
    file.path(tab_dir, "stage4_multinomial_RRR.csv")
  )
  
  message("\nHead of multinomial RRR table:")
  print(head(tidy_mn, 20))
  
  # A smaller key-terms table is also exported to make interpretation easier.
  # This subset keeps the predictors most likely to be discussed in the write-up.
  tidy_mn_key <- tidy_mn %>%
    filter(term != "(Intercept)") %>%
    filter(
      stringr::str_detect(term, "^educ") |
        stringr::str_detect(term, "^wealth") |
        stringr::str_detect(term, "^contraceptive_use_mn") |
        stringr::str_detect(term, "^fert_pref_mn") |
        stringr::str_detect(term, "^fp_radio") |
        stringr::str_detect(term, "^fp_tv") |
        stringr::str_detect(term, "^fp_print") |
        term == "share_sec_plus"
    )
  
  readr::write_csv(
    tidy_mn_key,
    file.path(tab_dir, "stage4_multinomial_RRR_key_terms.csv")
  )
  
  message("\nHead of multinomial key-terms table:")
  print(head(tidy_mn_key, 30))
}

# --------------------------
# 5i) Export a labelled multinomial table
# --------------------------

# Raw model term names are not very readable, so a labelled version
# is created for easier inspection and possible appendix use.
if (!is.null(m_mn_final)) {
  tidy_mn_labeled <- tidy_mn %>%
    mutate(term_label = label_term_mn(term)) %>%
    select(outcome_level, term_label, RRR, LCL, UCL, p.value)
  
  readr::write_csv(
    tidy_mn_labeled,
    file.path(tab_dir, "stage4_multinomial_RRR_labeled.csv")
  )
  
  message("\nHead of labeled multinomial table:")
  print(head(tidy_mn_labeled, 30))
}

# --------------------------
# 5j) Multinomial plots
# --------------------------

# Two forest-style plots are produced if a multinomial model succeeds:
#   1) a broader selected-terms plot
#   2) a more focused plot for education, fertility preference,
#      and contraceptive use
if (!is.null(m_mn_final)) {
  caption_mn <- paste0(
    "Data: The Gambia DHS 2019–20. Survey-adjusted multinomial regression; ",
    "95% confidence intervals via Taylor linearisation.\n",
    "Reference outcome = no HRFB."
  )
  
  # The first plot keeps a selected set of substantively important predictors.
  mn_plot_df <- tidy_mn_labeled %>%
    dplyr::filter(term_label != "(Intercept)") %>%
    dplyr::filter(
      stringr::str_detect(term_label, "^Education:") |
        stringr::str_detect(term_label, "^Wealth:") |
        stringr::str_detect(term_label, "^Fertility preference:") |
        stringr::str_detect(term_label, "^Contraceptive use:") |
        stringr::str_detect(term_label, "^FP message exposure \\(TV\\):") |
        term_label == "Cluster share with ≥ secondary education"
    ) %>%
    dplyr::mutate(
      block = term_block_order_mn(
        dplyr::case_when(
          stringr::str_detect(term_label, "^Education:") ~ "educ",
          stringr::str_detect(term_label, "^Wealth:") ~ "wealth",
          stringr::str_detect(term_label, "^Fertility preference:") ~ "fert_pref_mn",
          stringr::str_detect(term_label, "^Contraceptive use:") ~ "contraceptive_use_mn",
          stringr::str_detect(term_label, "^FP message exposure \\(TV\\):") ~ "fp_tv",
          term_label == "Cluster share with ≥ secondary education" ~ "share_sec_plus",
          TRUE ~ "other"
        )
      )
    ) %>%
    dplyr::arrange(outcome_level, block, term_label) %>%
    dplyr::mutate(
      term_label = factor(term_label, levels = rev(unique(term_label))),
      outcome_level = factor(
        outcome_level,
        levels = c("single", "multiple"),
        labels = c("Single HRFB vs none", "Multiple HRFB vs none")
      )
    )
  
  p_mn_forest <- ggplot(mn_plot_df, aes(x = term_label, y = RRR)) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.15, linewidth = 0.5) +
    geom_point(size = 1.8) +
    coord_flip() +
    facet_wrap(~ outcome_level, ncol = 2) +
    scale_y_log10(
      breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32),
      labels = scales::number_format(accuracy = 0.01),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    labs(
      title = "Selected multinomial associations with high-risk fertility behaviour",
      subtitle = "Relative risk ratios comparing single and multiple HRFB with no HRFB",
      x = NULL,
      y = "Relative risk ratio (log scale)",
      caption = caption_mn
    ) +
    theme_dissertation(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(10, 12, 10, 10)
    )
  
  print(p_mn_forest)
  save_fig(p_mn_forest, "stage4_multinomial_forest_selected", fig_dir = fig_dir, w = 10, h = 7.5)
  
  # The second plot narrows attention to education, fertility preference,
  # and contraceptive use, which may be easier to discuss in the dissertation text.
  mn_focus_df <- tidy_mn_labeled %>%
    dplyr::filter(term_label != "(Intercept)") %>%
    dplyr::filter(
      stringr::str_detect(term_label, "^Education:") |
        stringr::str_detect(term_label, "^Fertility preference:") |
        stringr::str_detect(term_label, "^Contraceptive use:")
    ) %>%
    dplyr::mutate(
      block = dplyr::case_when(
        stringr::str_detect(term_label, "^Education:") ~ 1,
        stringr::str_detect(term_label, "^Fertility preference:") ~ 2,
        stringr::str_detect(term_label, "^Contraceptive use:") ~ 3,
        TRUE ~ 99
      )
    ) %>%
    dplyr::arrange(outcome_level, block, term_label) %>%
    dplyr::mutate(
      term_label = factor(term_label, levels = rev(unique(term_label))),
      outcome_level = factor(
        outcome_level,
        levels = c("single", "multiple"),
        labels = c("Single HRFB vs none", "Multiple HRFB vs none")
      )
    )
  
  p_mn_focus <- ggplot(mn_focus_df, aes(x = term_label, y = RRR)) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.15, linewidth = 0.5) +
    geom_point(size = 1.8) +
    coord_flip() +
    facet_wrap(~ outcome_level, ncol = 2) +
    scale_y_log10(
      breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32),
      labels = scales::number_format(accuracy = 0.01),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    labs(
      title = "Multinomial associations for education, reproductive intentions, and contraceptive use",
      subtitle = "Relative risk ratios comparing single and multiple HRFB with no HRFB",
      x = NULL,
      y = "Relative risk ratio (log scale)",
      caption = caption_mn
    ) +
    theme_dissertation(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(10, 12, 10, 10)
    )
  
  print(p_mn_focus)
  save_fig(p_mn_focus, "stage4_multinomial_forest_focus", fig_dir = fig_dir, w = 10, h = 6.8)
}

# --------------------------
# 5k) Write a short multinomial note
# --------------------------

# A separate README note is written so that the exploratory status of the
# multinomial model is documented independently of the primary binary model.
mn_note_path <- file.path(out_dir, "README_stage4_multinomial_note.txt")

mn_note <- c(
  "STAGE 4 — Multinomial model note",
  "",
  "This stage explored survey-adjusted multinomial regression for hrfb_multicat (none/single/multiple).",
  "To reduce sparse-cell instability, multinomial-specific recodes were used:",
  "- contraceptive_use_mn: not using / modern method / traditional-other / Missing",
  "- fert_pref_mn: have another / undecided / no more / other-infecund-missing",
  "",
  paste0("Final selected model for export: ", ifelse(is.null(m_mn_label), "none (all failed)", m_mn_label)),
  "",
  "Files written if a model succeeded:",
  "- tables_stage4/stage4_multinomial_RRR.csv",
  "- tables_stage4/stage4_multinomial_RRR_key_terms.csv",
  "- tables_stage4/stage4_multinomial_RRR_labeled.csv",
  "- figs_stage4/stage4_multinomial_forest_selected.pdf/.png",
  "- figs_stage4/stage4_multinomial_forest_focus.pdf/.png"
)

writeLines(mn_note, mn_note_path)
message("\nMultinomial note written to: ", normalizePath(mn_note_path))

# --------------------------
# 6) Sensitivity analyses
# --------------------------

# Age-stratified models are fitted as a sensitivity check to reduce
# concern about age-related circularity around the HRFB outcome.
df <- df %>%
  mutate(
    age_group_strat = case_when(
      age_at_birth_years < 25 ~ "<25",
      age_at_birth_years >= 35 ~ ">=35",
      TRUE ~ "25-34"
    )
  )

fit_in_stratum <- function(level) {
  svy_strat <- subset(svy_base, age_group_strat == level)
  svyglm(f_bin, design = svy_strat, family = quasibinomial())
}

m_lt25 <- try(fit_in_stratum("<25"), silent = TRUE)
m_ge35 <- try(fit_in_stratum(">=35"), silent = TRUE)

tb_lt25 <- export_or(m_lt25, file.path(tab_dir, "stage4_primary_logistic_lt25_OR.csv"))
tb_ge35 <- export_or(m_ge35, file.path(tab_dir, "stage4_primary_logistic_ge35_OR.csv"))

# A one-birth-per-woman restriction is also fitted to remove within-woman
# dependence entirely and check whether the main binary model is being
# driven disproportionately by women contributing multiple births.
df_onebirth <- df %>%
  filter(!is.na(woman_id), !is.na(child_age_months)) %>%
  group_by(woman_id) %>%
  slice_min(order_by = child_age_months, n = 1, with_ties = FALSE) %>%
  ungroup()

qa_onebirth <- df %>%
  summarise(
    n_births_total = n(),
    n_women = n_distinct(woman_id),
    n_births_kept = nrow(df_onebirth),
    pct_removed = 100 * (1 - nrow(df_onebirth) / n())
  )

print(qa_onebirth)
readr::write_csv(qa_onebirth, file.path(tab_dir, "stage4_onebirth_QA.csv"))

svy_onebirth <- survey::svydesign(
  ids = ~v021,
  strata = ~strata,
  weights = ~weight,
  data = df_onebirth %>% filter(!is.na(weight), !is.na(v021), !is.na(strata)),
  nest = TRUE
)

m_bin_onebirth <- svyglm(
  formula = f_bin,
  design = svy_onebirth,
  family = quasibinomial()
)

tidy_onebirth <- broom::tidy(
  m_bin_onebirth,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
) %>%
  rename(OR = estimate, LCL = conf.low, UCL = conf.high) %>%
  select(term, OR, LCL, UCL, std.error, statistic, p.value)

readr::write_csv(
  tidy_onebirth,
  file.path(tab_dir, "stage4_primary_logistic_onebirth_per_woman_OR.csv")
)

tidy_main <- tidy_bin %>%
  select(term, OR_main = OR, LCL_main = LCL, UCL_main = UCL)

comp <- tidy_onebirth %>%
  select(term, OR_onebirth = OR, LCL_onebirth = LCL, UCL_onebirth = UCL) %>%
  left_join(tidy_main, by = "term") %>%
  mutate(pct_change_or = 100 * (OR_onebirth - OR_main) / OR_main)

readr::write_csv(
  comp,
  file.path(tab_dir, "stage4_compare_main_vs_onebirth_OR.csv")
)

print(head(comp, 12))

# --------------------------
# 7) Write Stage 4 README
# --------------------------

readme_path <- file.path(out_dir, "README_stage4_regression.txt")
readme_txt <- c(
  "STAGE 4 — Regression modelling (raw outcomes; design-adjusted)",
  "",
  "Inputs:",
  "- outputs/step02_births_qc_enhanced.rds (working dataset: raw + QA variables)",
  "",
  "Primary model:",
  "- svyglm with quasibinomial family, outcome: hrfb_binary.",
  "- Predictors: residence, region, education, wealth, ethnicity, FP exposure (radio/TV/print),",
  "  contraceptive_use, fert_pref, and share_sec_plus.",
  "- Maternal age at birth, parity, and short birth interval are excluded because they define the HRFB outcome itself.",
  "",
  "Secondary analyses:",
  "- Model comparison summaries (null, core, full).",
  "- Residual and functional-form diagnostics.",
  "- Adjusted prevalence plots from the primary binary model.",
  "- Optional multinomial modelling for hrfb_multicat.",
  "- Sensitivity models for <25 and >=35 age strata.",
  "- One-birth-per-woman sensitivity model.",
  "",
  "Tables:",
  paste0("- ", file.path("tables_stage4", "stage4_primary_logistic_OR.csv")),
  paste0("- ", file.path("tables_stage4", "appendix_fp_exposure_terms.csv")),
  paste0("- ", file.path("tables_stage4", "rq2_pseudoR2.csv")),
  paste0("- ", file.path("tables_stage4", "rq2_joint_wald_tests.csv")),
  paste0("- ", file.path("tables_stage4", "rq2_adjprev_by_residence.csv")),
  paste0("- ", file.path("tables_stage4", "rq2_adjprev_by_education.csv")),
  paste0("- ", file.path("tables_stage4", "rq2_adjprev_by_wealth.csv")),
  paste0("- ", file.path("tables_stage4", "rq2_adjprev_by_region.csv")),
  paste0("- ", file.path("tables_stage4", "rq2_adjprev_by_ethnicity.csv")),
  paste0("- ", file.path("tables_stage4", "rq2_adjprev_by_share_sec_plus.csv")),
  paste0("- ", file.path("tables_stage4", "stage4_primary_logistic_lt25_OR.csv"), " (if fitted)"),
  paste0("- ", file.path("tables_stage4", "stage4_primary_logistic_ge35_OR.csv"), " (if fitted)"),
  paste0("- ", file.path("tables_stage4", "stage4_primary_logistic_onebirth_per_woman_OR.csv")),
  paste0("- ", file.path("tables_stage4", "stage4_compare_main_vs_onebirth_OR.csv")),
  paste0("- ", file.path("tables_stage4", "stage4_multinomial_RRR.csv"), " (if fitted)"),
  "",
  "Figures:",
  "- PDF and PNG versions are written for all main plotted outputs in figs_stage4/",
  paste0("- ", file.path("figs_stage4", "stage4_primary_logistic_forest.*")),
  paste0("- ", file.path("figs_stage4", "rq2_forest_fertpref_contraception.*")),
  paste0("- ", file.path("figs_stage4", "rq2_adjprev_by_residence.*")),
  paste0("- ", file.path("figs_stage4", "rq2_adjprev_by_education.*")),
  paste0("- ", file.path("figs_stage4", "rq2_adjprev_by_wealth.*")),
  paste0("- ", file.path("figs_stage4", "rq2_adjprev_by_region.*")),
  paste0("- ", file.path("figs_stage4", "rq2_adjprev_by_ethnicity.*")),
  paste0("- ", file.path("figs_stage4", "rq2_adjprev_by_share_sec_plus.*")),
  paste0("- ", file.path("figs_stage4", "rq2_resid_pearson_vs_fitted.*")),
  paste0("- ", file.path("figs_stage4", "rq2_resid_deviance_vs_fitted.*")),
  paste0("- ", file.path("figs_stage4", "rq2_resid_pearson_vs_share_sec_plus.*")),
  "",
  "Notes:",
  "- All main inferences are design-adjusted using weights, PSU identifiers, and strata.",
  "- QA variables carried forward from Stage 2 remain in the working dataset but are not used as primary outcomes.",
  "- The multinomial model is treated as exploratory and may be omitted if instability prevents reliable estimation.",
  "",
  paste0("Saved on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)

writeLines(readme_txt, readme_path)

message("Wrote README: ", normalizePath(readme_path))
message("Stage 4 complete. Outputs written to: ", normalizePath(out_dir))

# ------------------------------------------------------------------
