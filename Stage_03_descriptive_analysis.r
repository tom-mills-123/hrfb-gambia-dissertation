############################################################
# STAGE 3 — Descriptive Analysis (RQ1)
# Project: High-risk fertility behaviours (HRFB) — Gambia DHS 2019–20
#
# What this script does:
#   1) Load the Stage 2 enriched working dataset and rebuild the
#      survey design on this dataset.
#   2) Enforce a raw-only policy for all primary descriptive results.
#   3) Estimate national, design-adjusted prevalence of any HRFB and
#      each component indicator with standard errors and 95% CIs.
#   4) Estimate subgroup prevalence of any HRFB by key covariates and
#      run Rao–Scott chi-square tests for group differences.
#   5) Export publication-ready tables and figures.
#   6) Plot subgroup prevalence of any HRFB.
#   7) Plot HRFB components by residence.
#   8) Export a regional prevalence table.
#   9) Write a Stage 3 README documenting outputs and the raw-only policy.
#
# Key principle for this stage:
#   Stage 3 uses the Stage 2 dataset as the working data container
#   because it carries forward the cleaned Stage 1 variables together
#   with the Stage 2 QA additions. However, all primary descriptive
#   outputs in this script are calculated strictly from the raw HRFB
#   outcome and component variables carried forward from Stage 1.
############################################################

# --------------------------
# 0) Packages + global options
# --------------------------
AUTO_INSTALL <- FALSE

# Import required packages
pkgs <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "forcats",
  "ggplot2",
  "readr",
  "janitor",
  "survey",
  "srvyr",
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
  library(ggplot2)
  library(readr)
  library(janitor)
  library(survey)
  library(srvyr)
  library(scales)
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

# Grouped design-adjusted prevalence for a binary outcome.
# This helper returns the point estimate, standard error, 95% CI,
# and a weighted denominator for each subgroup.
svy_group_prev <- function(design_obj, outcome, group_var) {
  design_obj %>%
    group_by(!!as.name(group_var)) %>%
    summarise(
      prev = survey_mean(!!as.name(outcome), vartype = c("se", "ci"), na.rm = TRUE),
      wt_n = survey_total(1, na.rm = TRUE)
    ) %>%
    rename(group = !!as.name(group_var)) %>%
    arrange(desc(prev))
}

# Rao–Scott chi-square test for differences in a binary outcome across groups.
# A standard survey design object is rebuilt internally from the srvyr design
# variables to keep the test specification explicit.
svy_rs_chisq <- function(design_obj, outcome, group_var) {
  dat <- design_obj$variables %>%
    dplyr::filter(!is.na(.data[[outcome]]), !is.na(.data[[group_var]]))
  
  svy_obj <- survey::svydesign(
    ids = ~v021,
    strata = ~strata,
    weights = ~weight,
    data = dat,
    nest = TRUE
  )
  
  f <- as.formula(paste0("~", outcome, " + ", group_var))
  rs <- survey::svychisq(f, design = svy_obj, statistic = "Chisq")
  
  tibble::tibble(
    outcome = outcome,
    group = group_var,
    statistic = unname(rs$statistic),
    df = unname(rs$parameter),
    p_value = unname(rs$p.value)
  )
}

# A common caption is used across the main descriptive plots.
caption_rq1 <- "Data: The Gambia DHS 2019–20. Weighted estimates accounting for complex survey design; 95% confidence intervals."

# --------------------------
# 1) Inputs and outputs
# --------------------------

# Stage 3 uses the Stage 2 enriched dataset as its working input.
# This preserves the full QA trail while allowing the primary
# descriptive analysis to remain based on the raw variables only.
in_births <- "outputs/step02_births_qc_enhanced.rds"
stopifnot(file.exists(in_births))

df <- readRDS(in_births)

# Output directories for Stage 3 tables, figures, and README.
out_dir <- "outputs"
fig_dir <- file.path(out_dir, "figs_stage3")
tab_dir <- file.path(out_dir, "tables_stage3")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

# --------------------------
# 2) Rebuild the survey design on the Stage 3 working dataset
# --------------------------

# The Stage 3 survey design is rebuilt directly on the current working
# dataset so that all subsequent weighted estimates align exactly to the
# data frame being analysed in this script.
stopifnot(all(c("weight", "v021", "strata") %in% names(df)))

df_design <- df %>%
  dplyr::filter(!is.na(weight), !is.na(v021), !is.na(strata))

design_svy <- srvyr::as_survey_design(
  .data = df_design,
  ids = v021,
  strata = strata,
  weights = weight,
  nest = TRUE
)

message(
  "Design built on df_design: ",
  nrow(df_design), " rows; ",
  dplyr::n_distinct(df_design$v021), " PSUs; ",
  dplyr::n_distinct(df_design$strata), " strata."
)

# --------------------------
# 3) Guardrails: enforce raw-only analysis
# --------------------------

# A backward-compatibility safeguard is included here in case the Stage 2
# dataset does not already contain age_at_birth_years as a named variable.
if (!("age_at_birth_years" %in% names(df)) && ("age_at_birth" %in% names(df))) {
  df <- df %>% mutate(age_at_birth_years = age_at_birth)
}

# These are the raw variables required for the primary Stage 3 results.
# If any are missing, the script stops immediately.
raw_required <- c(
  "age_at_birth_years",
  "hrfb_binary", "hrfb_multicat",
  "risk_age_young", "risk_age_old", "risk_short_int", "risk_high_parity",
  "residence", "region", "educ", "wealth", "ethnicity"
)

missing_raw <- setdiff(raw_required, names(df))
if (length(missing_raw) > 0) {
  stop(
    "Missing required raw variables in step02 dataset: ",
    paste(missing_raw, collapse = ", ")
  )
}

# Stage 2 sensitivity variables may still be present in the dataset.
# Their presence is expected, but they are not used in the primary Stage 3 outputs.
qa_variants_present <- intersect(
  c(
    "age_band",
    "hrfb_binary_banded", "hrfb_multicat_banded",
    "age_at_birth_smooth",
    "hrfb_binary_smooth", "hrfb_multicat_smooth"
  ),
  names(df)
)

if (length(qa_variants_present) > 0) {
  message(
    "QA variants detected (expected from Stage 2): ",
    paste(qa_variants_present, collapse = ", "),
    "\n→ These will NOT be used for Stage 3 primary descriptives."
  )
}

# The primary Stage 3 descriptive outcome is the binary HRFB indicator.
# Component prevalence and subgroup prevalence are both built around this.
outcome_primary <- "hrfb_binary"
subgroup_vars <- c("residence", "region", "educ", "wealth", "ethnicity")

# --------------------------
# 4) National prevalence (raw, design-adjusted)
# --------------------------

# This section reports the headline prevalence of any HRFB and each of
# the four component indicators, with standard errors and 95% confidence intervals.
national_prev <- design_svy %>%
  summarise(
    prev_hrfb_any    = survey_mean(hrfb_binary, vartype = c("se", "ci"), na.rm = TRUE),
    prev_age_young   = survey_mean(risk_age_young, vartype = c("se", "ci"), na.rm = TRUE),
    prev_age_old     = survey_mean(risk_age_old, vartype = c("se", "ci"), na.rm = TRUE),
    prev_short_int   = survey_mean(risk_short_int, vartype = c("se", "ci"), na.rm = TRUE),
    prev_high_parity = survey_mean(risk_high_parity, vartype = c("se", "ci"), na.rm = TRUE)
  ) %>%
  janitor::clean_names()

readr::write_csv(national_prev, file.path(tab_dir, "stage3_national_prevalence_raw.csv"))
print(national_prev)

# --------------------------
# 5) Subgroup prevalence of any HRFB + Rao–Scott chi-square tests
# --------------------------

# Subgroup prevalence tables are produced for the primary binary outcome.
subgroup_prev <- purrr::map_dfr(
  subgroup_vars,
  ~ svy_group_prev(design_svy, outcome = outcome_primary, group_var = .x) %>%
    mutate(var = .x) %>%
    select(var, group, prev, prev_se, prev_low, prev_upp, wt_n, wt_n_se)
)

readr::write_csv(subgroup_prev, file.path(tab_dir, "stage3_subgroup_prevalence_raw.csv"))
print(head(subgroup_prev, 10))

# Rao–Scott chi-square tests assess whether the prevalence of any HRFB
# differs significantly across each subgroup variable.
subgroup_tests <- purrr::map_dfr(
  subgroup_vars,
  ~ svy_rs_chisq(design_svy, outcome_primary, .x)
)

readr::write_csv(subgroup_tests, file.path(tab_dir, "stage3_subgroup_raoscott_tests_raw.csv"))
print(subgroup_tests)

# --------------------------
# 6) Plots — subgroup prevalence with 95% confidence intervals
# --------------------------

# This plotting helper draws a single bar chart for one subgroup variable,
# ordered by prevalence and annotated with 95% confidence intervals.
plot_subgroup <- function(df_prev, var_label, file_stub, title_pretty = NULL) {
  df_plot <- df_prev %>%
    filter(var == var_label) %>%
    mutate(group = as.character(group)) %>%
    filter(!is.na(group)) %>%
    arrange(prev) %>%
    mutate(group = forcats::fct_inorder(group))
  
  if (is.null(title_pretty)) title_pretty <- str_to_title(var_label)
  
  p <- ggplot(df_plot, aes(x = group, y = prev)) +
    geom_col() +
    geom_errorbar(aes(ymin = prev_low, ymax = prev_upp), width = 0.2) +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, 1),
      oob = scales::squish,
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = paste0("Prevalence of any HRFB by ", title_pretty),
      subtitle = "Weighted (DHS complex design); 95% confidence intervals",
      x = NULL,
      y = "Prevalence",
      caption = caption_rq1
    ) +
    theme_dissertation()
  
  print(p)
  save_fig(p, paste0("stage3_prev_by_", file_stub, "_RAW"), fig_dir = fig_dir, w = 8, h = 5)
}

plot_subgroup(subgroup_prev, "residence", "residence", "Residence")
plot_subgroup(subgroup_prev, "region", "region", "Region")
plot_subgroup(subgroup_prev, "educ", "education", "Education")
plot_subgroup(subgroup_prev, "wealth", "wealth", "Wealth")
plot_subgroup(subgroup_prev, "ethnicity", "ethnicity", "Ethnicity")

# --------------------------
# 7) Components by residence (illustrative component plot)
# --------------------------

# This component plot is included to show which dimensions of HRFB
# differ most clearly between urban and rural settings.
component_prev <- purrr::imap_dfr(
  c(
    age_young = "risk_age_young",
    age_old = "risk_age_old",
    short_int = "risk_short_int",
    parity4p = "risk_high_parity"
  ),
  ~ svy_group_prev(design_svy, outcome = .x, group_var = "residence") %>%
    mutate(component = .y) %>%
    select(component, group, prev, prev_se, prev_low, prev_upp)
)

readr::write_csv(component_prev, file.path(tab_dir, "stage3_components_by_residence_RAW.csv"))

component_prev <- component_prev %>%
  mutate(
    component = recode(
      component,
      "age_young" = "Age < 18",
      "age_old" = "Age ≥ 35",
      "short_int" = "Birth interval < 24m",
      "parity4p" = "Parity ≥ 4"
    ),
    component = factor(
      component,
      levels = c("Age < 18", "Age ≥ 35", "Birth interval < 24m", "Parity ≥ 4")
    )
  )

p_comp <- ggplot(component_prev, aes(x = component, y = prev, fill = group)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = prev_low, ymax = prev_upp),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Components of HRFB by residence",
    subtitle = "Weighted (DHS complex design); 95% confidence intervals",
    x = NULL,
    y = "Prevalence",
    fill = "Residence",
    caption = caption_rq1
  ) +
  theme_dissertation()

print(p_comp)
save_fig(p_comp, "stage3_components_by_residence_RAW", fig_dir = fig_dir, w = 8, h = 5)

# --------------------------
# 8) Regional prevalence table
# --------------------------

# Regional prevalence of any HRFB is exported as a table in all cases.
region_prev <- svy_group_prev(design_svy, outcome = outcome_primary, group_var = "region") %>%
  rename(region_label = group) %>%
  mutate(region_label = as.character(region_label))

readr::write_csv(region_prev, file.path(tab_dir, "stage3_region_prevalence_RAW.csv"))



# --------------------------
# 9) Write Stage 3 README
# --------------------------

# The Stage 3 README records the raw-only policy and documents the main
# tables and figures produced in this script.
readme_path <- file.path(out_dir, "README_stage3_descriptives.txt")
readme_txt <- c(
  "STAGE 3 — Descriptive Analysis (raw inputs for all primary results)",
  "",
  "Inputs:",
  "- outputs/step02_births_qc_enhanced.rds  (working dataset: raw + QA variables)",
  "",
  "Policy:",
  "- All Stage 3 prevalence estimates, subgroup comparisons, chi-square tests, and plots use raw variables only.",
  "- The primary outcome is hrfb_binary (any HRFB).",
  "- Component descriptives use the raw Stage 1 HRFB component indicators.",
  "- Stage 2 banded and smoothed variables are retained in the dataset for sensitivity analyses later in the pipeline, but are not used here.",
  "",
  "Tables:",
  paste0("- ", file.path("tables_stage3", "stage3_national_prevalence_raw.csv")),
  paste0("- ", file.path("tables_stage3", "stage3_subgroup_prevalence_raw.csv")),
  paste0("- ", file.path("tables_stage3", "stage3_subgroup_raoscott_tests_raw.csv")),
  paste0("- ", file.path("tables_stage3", "stage3_components_by_residence_RAW.csv")),
  paste0("- ", file.path("tables_stage3", "stage3_region_prevalence_RAW.csv")),
  "",
  "Figures:",
  "- PDF and PNG versions are written for all main plotted outputs in figs_stage3/",
  paste0("- ", file.path("figs_stage3", "stage3_prev_by_residence_RAW.*")),
  paste0("- ", file.path("figs_stage3", "stage3_prev_by_region_RAW.*")),
  paste0("- ", file.path("figs_stage3", "stage3_prev_by_education_RAW.*")),
  paste0("- ", file.path("figs_stage3", "stage3_prev_by_wealth_RAW.*")),
  paste0("- ", file.path("figs_stage3", "stage3_prev_by_ethnicity_RAW.*")),
  paste0("- ", file.path("figs_stage3", "stage3_components_by_residence_RAW.*")),
  "",
  "Notes:",
  "- Estimates are design-adjusted using weights, PSU identifiers (v021), and the carried-forward strata variable.",
  "- Using the Stage 2 dataset preserves quality-control provenance without changing the primary descriptive definitions used in Stage 3.",
  "",
  paste0("Saved on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)

writeLines(readme_txt, readme_path)

message("Wrote README: ", normalizePath(readme_path))
message("Stage 3 complete (raw-only). Outputs written to: ", normalizePath(out_dir))
############################################################