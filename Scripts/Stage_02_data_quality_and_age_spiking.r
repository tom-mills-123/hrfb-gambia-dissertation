############################################################
# STAGE 2 — Data Quality Assessment: Age Heaping & Spiking
# Project: High-risk fertility behaviours (HRFB) — Gambia DHS 2019–20
#
# What this script does (and why):
#   1) Load the Stage 1 birth-level checkpoint.
#   2) Create tidy age variables and basic histograms to visually detect heaping
#      for mother’s age-at-birth and child age (b19/derived).
#   3) Quantify heaping with standard indices:
#        - Whipple’s Index (preference for terminal digits 0/5, ages 23–62)
#        - Myers’ Blended Index (digit preference 0–9 across adult ages)
#      (computed for age-at-birth in whole years, and for child months as a
#       light-touch diagnostic on multiples of 6/12).
#   4) Produce digit-preference plots (0–9) and save all figures.
#   5) Mitigation strategies:
#        - Banding: group age-at-birth into 5-year bands (15–19, 20–24, …)
#        - Smoothing (only for sensitivity): conservative, reproducible jitter
#          applied to exact integer ages ending in 0 or 5, then reclassify.
#   6) Sensitivity analysis:
#        - Recompute HRFB outcomes with (a) raw ages, (b) banded, (c) smoothed
#        - Compare weighted prevalence (survey design) across versions.
#   7) Save:
#        - Figures (histograms, digit preference)
#        - Heaping indices table (appendix)
#        - Updated objects for downstream stages
#        - README for Stage 2
############################################################

# --------------------------
# 0) Packages + global options
# --------------------------
AUTO_INSTALL <- FALSE

# Install packages required
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

# Configure variance estimation when a stratum contains only one PSU after subsetting
options(survey.lonely.psu = "adjust")

# --------------------------
# Helper functions
# --------------------------

# Define a consistent plotting theme
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

# Save figures in both PDF and PNG
save_fig <- function(p, filename, fig_dir, w = 7, h = 4.5) {
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



# --------------------------
# 1) Load Stage 1 birth-level checkpoint
# --------------------------
in_births <- "outputs/step01_births_clean.rds"
stopifnot(file.exists(in_births))

df <- readRDS(in_births)

# Output dirs for this stage
out_dir       <- "outputs"
fig_dir       <- file.path(out_dir, "figs_stage2")
tab_dir       <- file.path(out_dir, "tables_stage2")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

# --------------------------
# 2) Prepare age variables
#    NOTE: age_at_birth in Stage 1 is derived from CMC dates (b3, v011),
#    so classic age-heaping indices (Whipple/Myers) are not strictly applicable.
#    Here we use integer ages only for plotting/QC convenience and threshold checks.
# --------------------------
df <- df %>%
  mutate(
    age_at_birth_years     = age_at_birth,                 # numeric, possibly fractional
    age_at_birth_years_int = floor(age_at_birth_years),    # completed years (avoid round() artifacts)
    child_age_months_int   = floor(child_age_months)       # completed months
  )


# --------------------------
# 2A) Diagnostic: Is age_at_birth_years effectively an integer?
# Why: Terminal-digit heaping is meaningful only if ages are recorded (or constructed) in completed years.
#      If age_at_birth_years is genuinely fractional (e.g., derived from months/dates), using floor()
#      can distort terminal-digit patterns by systematically shifting ages downward.
# What we check:
#   - prop_near_integer: share of ages extremely close to an integer (tolerance handles floating-point noise)
#   - frac_part summary: distribution of fractional parts (if spread across 0–1, ages are truly fractional)
# Output:
#   - Console message + CSV saved for audit trail / appendix methods.
# --------------------------
prop_near_integer <- mean(abs(df$age_at_birth_years - round(df$age_at_birth_years)) < 1e-8, na.rm = TRUE)

frac_part <- df$age_at_birth_years - floor(df$age_at_birth_years)

age_integer_diagnostic <- tibble::tibble(
  prop_near_integer = prop_near_integer,
  frac_p01 = as.numeric(quantile(frac_part, 0.01, na.rm = TRUE)),
  frac_p25 = as.numeric(quantile(frac_part, 0.25, na.rm = TRUE)),
  frac_p50 = as.numeric(quantile(frac_part, 0.50, na.rm = TRUE)),
  frac_p75 = as.numeric(quantile(frac_part, 0.75, na.rm = TRUE)),
  frac_p99 = as.numeric(quantile(frac_part, 0.99, na.rm = TRUE))
)

message("Age integer-likeness check: prop_near_integer = ",
        round(prop_near_integer, 3),
        " (rule of thumb: >0.95 suggests ages are effectively whole-year).")

readr::write_csv(age_integer_diagnostic,
                 file.path(tab_dir, "stage2_age_integer_diagnostic.csv"))
print(age_integer_diagnostic)


# Quick QA snapshot of ranges
message("Age-at-birth (yrs) range: ",
        round(min(df$age_at_birth_years, na.rm = TRUE),1),"–",
        round(max(df$age_at_birth_years, na.rm = TRUE),1))
message("Child age (months) range: ",
        min(df$child_age_months_int, na.rm = TRUE),"–",
        max(df$child_age_months_int, na.rm = TRUE))

# --------------------------
# 3) Visual diagnostics
# --------------------------

# The first set of checks is graphical. These figures are intended to make
# any obvious spikes or clustering visible before moving to summary measures.

# 3a) Histogram: age-at-birth (years, numeric)
p_age_hist <- ggplot(df, aes(x = age_at_birth_years)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left",
                 color = "white", linewidth = 0.2) +
  scale_x_continuous(breaks = seq(10, 60, by = 5)) +
  labs(
    title = "Mother’s age at birth",
    subtitle = "Unweighted diagnostic histogram (1-year bins) to visually assess heaping",
    x = "Age at birth (years)",
    y = "Number of births",
    caption = "Data: Gambia DHS 2019–20 (birth-level file)."
  ) +
  theme_dissertation()

print(p_age_hist)
save_fig(p_age_hist, "age_at_birth_histogram", fig_dir = fig_dir)

# NOTE: Terminal-digit heaping diagnostics are most interpretable when age_at_birth_years is effectively integer-valued.
# See Section 2A (stage2_age_integer_diagnostic.csv) for evidence supporting the integer conversion used here.

# 3b) Bar chart: terminal digit preference (0–9) for age-at-birth integer ages
digit_pref_df <- df %>%
  filter(!is.na(age_at_birth_years_int),
         age_at_birth_years_int >= 10,
         age_at_birth_years_int <= 60) %>%
  mutate(term_digit = age_at_birth_years_int %% 10) %>%
  count(term_digit, name = "n") %>%
  mutate(share = n / sum(n))

n_age <- sum(digit_pref_df$n)

p_digit <- ggplot(digit_pref_df, aes(x = factor(term_digit), y = share)) +
  geom_hline(yintercept = 0.10, linetype = "dashed", linewidth = 0.4) +
  geom_col(width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Terminal digit preference in reported maternal age",
    subtitle = "If ages are uniformly reported, each digit would be ~10% (dashed line)",
    x = "Terminal digit (0–9)",
    y = "Share of ages",
    caption = paste0("Integers only, ages 10–60. N = ", scales::comma(n_age), ".")
  ) +
  theme_dissertation()

print(p_digit)
save_fig(p_digit, "age_at_birth_terminal_digit_preference", fig_dir = fig_dir)


# 3c) Child age (months) quick check: multiples of 6/12
child_mod_df <- df %>%
  filter(!is.na(child_age_months_int)) %>%
  mutate(mod_6  = child_age_months_int %% 6,
         mod_12 = child_age_months_int %% 12) %>%
  summarise(
    pct_div6  = mean(mod_6  == 0) * 100,
    pct_div12 = mean(mod_12 == 0) * 100
  )
print(child_mod_df)

# visual: histogram of child age in months (shown + saved)
p_child_hist <- ggplot(df, aes(x = child_age_months_int)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left",
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = seq(0, 59, by = 12), linetype = "dotted", linewidth = 0.3) +
  scale_x_continuous(breaks = seq(0, 60, by = 6), limits = c(0, 59)) +
  labs(
    title = "Child age in months",
    subtitle = "Unweighted diagnostic histogram (dotted lines at 12-month intervals)",
    x = "Child age (months)",
    y = "Number of births",
    caption = "Used to assess month-of-age reporting spikes (e.g., 12, 24, 36)."
  ) +
  theme_dissertation() +
  theme(axis.text.x = element_text(size = 9))

print(p_child_hist)
save_fig(p_child_hist, "child_age_months_histogram", fig_dir = fig_dir)


# --------------------------
# 4) QC metrics (aligned with Stage 1 definitions)
#    Why: age_at_birth is CMC-derived, so Whipple/Myers (reported-age heaping indices)
#         are not strictly appropriate. Instead we quantify:
#         - digit preference on completed-year ages (descriptive only)
#         - boundary mass near HRFB thresholds (<18 and >34)
#         - child-month spikes at 12-month intervals (more relevant here)
# --------------------------

qc_summary <- df %>%
  summarise(
    n_births = n(),
    
    # How much mass sits close to key HRFB thresholds?
    pct_within_0_5yr_of_18 = mean(abs(age_at_birth_years - 18) <= 0.5, na.rm = TRUE) * 100,
    pct_within_0_5yr_of_35 = mean(abs(age_at_birth_years - 35) <= 0.5, na.rm = TRUE) * 100,
    
    # Month spikes for child age distribution
    pct_child_months_div12 = mean(child_age_months_int %% 12 == 0, na.rm = TRUE) * 100,
    pct_child_months_div6  = mean(child_age_months_int %% 6  == 0, na.rm = TRUE) * 100
  )

print(qc_summary)
readr::write_csv(qc_summary, file.path(tab_dir, "stage2_qc_summary_thresholds_and_spikes.csv"))

# 4a) Whipple’s Index
whipple_index <- function(age_int_vec) {
  a <- age_int_vec[!is.na(age_int_vec) & age_int_vec >= 23 & age_int_vec <= 62]
  if (length(a) == 0) return(NA_real_)
  is_05 <- (a %% 5 == 0)
  100 * (sum(is_05) / (length(a) / 5))
}


# 4b) Myers’ Blended Index
# Measures preference for each terminal digit 0–9; returns total heaping (0–90).
# Compute shares by terminal digit across a broad adult range (e.g., 10–80),
# then take the sum of absolute deviations from uniform (10% per digit), times 1/2.
myers_index <- function(age_int_vec, lower = 10, upper = 80) {
  a <- age_int_vec[!is.na(age_int_vec) & age_int_vec >= lower & age_int_vec <= upper]
  if (length(a) == 0) return(NA_real_)
  dtab <- table(a %% 10)
  # ensure digits 0..9 all present
  all_digits <- 0:9
  counts <- sapply(all_digits, function(d) ifelse(!is.na(dtab[as.character(d)]), as.numeric(dtab[as.character(d)]), 0))
  shares <- counts / sum(counts)
  50 * sum(abs(shares - 0.1)) # 0=no heaping, 90=max
}

# Compute indices for mother’s age at birth
W_mother <- whipple_index(df$age_at_birth_years_int)
M_mother <- myers_index(df$age_at_birth_years_int, lower = 10, upper = 60)

# For completeness, a light-touch index on child age months (divisible by 6/12)
# Not standard Whipple/Myers, but useful for spotting reporting spikes.
child_div6  <- mean(df$child_age_months_int %% 6  == 0, na.rm = TRUE) * 100
child_div12 <- mean(df$child_age_months_int %% 12 == 0, na.rm = TRUE) * 100

heaping_summary <- tibble::tibble(
  metric  = c("Whipple_age_at_birth_23_62", "Myers_age_at_birth_10_60",
              "Pct_child_age_months_divisible_by_6", "Pct_child_age_months_divisible_by_12"),
  value   = c(W_mother, M_mother, child_div6, child_div12),
  note    = c("100=no heaping; 500=all ages on 0/5",
              "0=no heaping; 90=max heaping",
              "Share (%) of reported months on 0,6,12,18,...",
              "Share (%) of reported months on 0,12,24,36,48")
)
readr::write_csv(heaping_summary, file.path(tab_dir, "stage2_heaping_indices.csv"))
print(heaping_summary)


# --------------------------
# 5) Mitigation Strategy A — Threshold-preserving banding
# --------------------------

# This banded version is included as a sensitivity exercise. The bands are chosen
# so that the Stage 1 HRFB thresholds are preserved as closely as possible.
age_breaks <- c(10,15,18,20,25,30,35,40,45,50,55,60) 
age_labels <- c("10-14","15-17","18-19","20-24","25-29","30-34",
                "35-39","40-44","45-49","50-54","55-59")

df <- df %>%
  mutate(
    age_band = cut(
      age_at_birth_years,
      breaks = age_breaks, labels = age_labels,
      right = FALSE, include.lowest = TRUE
    ),
    
    # Threshold-preserving banded flags
    risk_age_young_banded = if_else(!is.na(age_band) & age_band == "15-17", 1L, 0L, missing = 0L),
    risk_age_old_banded   = if_else(!is.na(age_band) & age_band %in% c("35-39","40-44","45-49","50-54","55-59"),
                                    1L, 0L, missing = 0L),
    
    # Composite outcomes using same non-age risks from Stage 1
    hrfb_count_banded    = risk_age_young_banded + risk_age_old_banded + risk_short_int + risk_high_parity,
    hrfb_binary_banded   = if_else(hrfb_count_banded >= 1, 1L, 0L),
    hrfb_multicat_banded = case_when(
      hrfb_count_banded == 0 ~ "none",
      hrfb_count_banded == 1 ~ "single",
      hrfb_count_banded >= 2 ~ "multiple",
      TRUE ~ NA_character_
    )
  )


# --------------------------
# 6) Sensitivity variant B — conservative smoothing
# --------------------------

# This smoothing step is included only as a secondary sensitivity check.
# It is applied only when age_at_birth behaves very similarly to an integer-valued
# variable, and only to exact integer ages ending in 0 or 5.
do_smoothing <- prop_near_integer > 0.95
message("Smoothing gate: do_smoothing = ", do_smoothing)

set.seed(2025)

if (do_smoothing) {
  df <- df %>%
    mutate(
      is_exact_int        = abs(age_at_birth_years - round(age_at_birth_years)) < 1e-8,
      ends_0_or_5         = (round(age_at_birth_years) %% 5 == 0),
      needs_jitter        = is_exact_int & ends_0_or_5,
      
      age_at_birth_smooth = if_else(
        needs_jitter,
        age_at_birth_years + runif(dplyr::n(), min = -0.5, max = 0.5),
        age_at_birth_years
      ),
      
      risk_age_young_smooth = if_else(age_at_birth_smooth < 18, 1L, 0L, missing = 0L),
      risk_age_old_smooth   = if_else(age_at_birth_smooth > 34, 1L, 0L, missing = 0L),
      
      hrfb_count_smooth     = risk_age_young_smooth + risk_age_old_smooth + risk_short_int + risk_high_parity,
      hrfb_binary_smooth    = if_else(hrfb_count_smooth >= 1, 1L, 0L),
      hrfb_multicat_smooth  = case_when(
        hrfb_count_smooth == 0 ~ "none",
        hrfb_count_smooth == 1 ~ "single",
        hrfb_count_smooth >= 2 ~ "multiple",
        TRUE ~ NA_character_
      )
    )
} else {
  # If smoothing not appropriate, define smooth vars identical to raw for downstream code consistency
  df <- df %>%
    mutate(
      age_at_birth_smooth     = age_at_birth_years,
      risk_age_young_smooth   = risk_age_young,
      risk_age_old_smooth     = risk_age_old,
      hrfb_count_smooth       = hrfb_risk_count,
      hrfb_binary_smooth      = hrfb_binary,
      hrfb_multicat_smooth    = hrfb_multicat
    )
}


# Visual comparison of the raw and smoothed age distributions
df_long_age <- df %>%
  transmute(raw = age_at_birth_years,
            smoothed = age_at_birth_smooth) %>%
  pivot_longer(everything(), names_to = "version", values_to = "age")

p_age_compare <- ggplot(df_long_age, aes(x = age, fill = version)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left",
                 position = "identity", alpha = 0.45,
                 color = "white", linewidth = 0.2) +
  scale_x_continuous(breaks = seq(10, 60, by = 5)) +
  labs(
    title = "Maternal age at birth: raw vs smoothed",
    subtitle = "Smoothing applied only to exact integer ages ending in 0 or 5 (sensitivity analysis)",
    x = "Age at birth (years)",
    y = "Number of births",
    fill = NULL
  ) +
  theme_dissertation()

print(p_age_compare)
save_fig(p_age_compare, "age_at_birth_histogram_raw_vs_smoothed", fig_dir = fig_dir)


# --------------------------
# 6A) Before vs after smoothing
# --------------------------

# These comparisons show whether smoothing materially changes the descriptive
# heaping diagnostics or whether the raw and smoothed results remain very similar.
df <- df %>%
  mutate(
    # Keep integer conversion consistent with earlier digit diagnostics
    # (completed years; matches Stage 2 Section 2)
    age_at_birth_years_int        = floor(age_at_birth_years),
    age_at_birth_years_int_smooth = floor(age_at_birth_smooth)
  )

W_mother_raw    <- whipple_index(df$age_at_birth_years_int)
M_mother_raw    <- myers_index(df$age_at_birth_years_int, lower = 10, upper = 60)

W_mother_smooth <- whipple_index(df$age_at_birth_years_int_smooth)
M_mother_smooth <- myers_index(df$age_at_birth_years_int_smooth, lower = 10, upper = 60)

heaping_compare <- tibble::tibble(
  metric  = c("Whipple_23_62","Myers_10_60"),
  raw     = c(W_mother_raw, M_mother_raw),
  smooth  = c(W_mother_smooth, M_mother_smooth),
  delta   = smooth - raw
)
print(heaping_compare)
readr::write_csv(heaping_compare,
                 file.path(tab_dir, "stage2_heaping_indices_before_after_smoothing.csv"))

# Terminal-digit preference AFTER smoothing (for transparency)
digit_compare <- bind_rows(
  df %>%
    filter(!is.na(age_at_birth_years_int),
           age_at_birth_years_int >= 10, age_at_birth_years_int <= 60) %>%
    mutate(term_digit = age_at_birth_years_int %% 10) %>%
    count(term_digit, name = "n") %>%
    mutate(share = n / sum(n), version = "Raw"),
  df %>%
    filter(!is.na(age_at_birth_years_int_smooth),
           age_at_birth_years_int_smooth >= 10, age_at_birth_years_int_smooth <= 60) %>%
    mutate(term_digit = age_at_birth_years_int_smooth %% 10) %>%
    count(term_digit, name = "n") %>%
    mutate(share = n / sum(n), version = "Smoothed")
)

p_digit_compare <- ggplot(digit_compare, aes(x = factor(term_digit), y = share)) +
  geom_hline(yintercept = 0.10, linetype = "dashed", linewidth = 0.4) +
  geom_col(width = 0.75) +
  facet_wrap(~ version, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Terminal digit preference: raw vs smoothed",
    subtitle = "Dashed line indicates 10% per digit under uniform reporting",
    x = "Terminal digit (0–9)",
    y = "Share of ages"
  ) +
  theme_dissertation() +
  theme(legend.position = "none")

print(p_digit_compare)
save_fig(p_digit_compare, "terminal_digit_preference_raw_vs_smoothed", fig_dir = fig_dir)


# --------------------------
# 6B) Agreement across versions
# --------------------------

# These tables quantify how often the sensitivity versions differ from the
# original Stage 1 HRFB classification.
agree_banded <- df %>%
  summarise(
    pct_agree        = mean(hrfb_binary == hrfb_binary_banded, na.rm = TRUE) * 100,
    pct_flip_to_risk = mean(hrfb_binary == 0 & hrfb_binary_banded == 1, na.rm = TRUE) * 100,
    pct_flip_to_no   = mean(hrfb_binary == 1 & hrfb_binary_banded == 0, na.rm = TRUE) * 100
  )
print(agree_banded)
readr::write_csv(agree_banded,
                 file.path(tab_dir, "stage2_agreement_raw_vs_banded.csv"))

agree_smooth <- df %>%
  summarise(
    pct_agree        = mean(hrfb_binary == hrfb_binary_smooth, na.rm = TRUE) * 100,
    pct_flip_to_risk = mean(hrfb_binary == 0 & hrfb_binary_smooth == 1, na.rm = TRUE) * 100,
    pct_flip_to_no   = mean(hrfb_binary == 1 & hrfb_binary_smooth == 0, na.rm = TRUE) * 100
  )
print(agree_smooth)
readr::write_csv(agree_smooth,
                 file.path(tab_dir, "stage2_agreement_raw_vs_smoothed.csv"))

# --------------------------
# 7) Weighted prevalence across raw, banded, and smoothed versions
# --------------------------

# The Stage 1 survey design is loaded above for continuity, but once the Stage 2
# variants have been added to the dataset, the design is rebuilt so the weighted
# summaries are aligned exactly to the enriched data frame.

design_svy2 <- srvyr::as_survey_design(
  .data   = df %>% dplyr::filter(!is.na(weight), !is.na(v021), !is.na(strata)),
  ids     = v021,
  strata  = strata,
  weights = weight,
  nest    = TRUE
)

prev_tbl <- design_svy2 %>%
  summarise(
    prev_hrfb_raw    = survey_mean(hrfb_binary,        vartype = c("ci","se"), na.rm = TRUE),
    prev_hrfb_banded = survey_mean(hrfb_binary_banded, vartype = c("ci","se"), na.rm = TRUE),
    prev_hrfb_smooth = survey_mean(hrfb_binary_smooth, vartype = c("ci","se"), na.rm = TRUE)
  ) %>%
  janitor::clean_names()

print(prev_tbl)
readr::write_csv(prev_tbl, file.path(tab_dir, "stage2_prev_sensitivity.csv"))

precision_tbl <- tibble::tibble(
  version  = c("raw","banded","smoothed"),
  estimate = c(prev_tbl$prev_hrfb_raw,
               prev_tbl$prev_hrfb_banded,
               prev_tbl$prev_hrfb_smooth),
  se       = c(prev_tbl$prev_hrfb_raw_se,
               prev_tbl$prev_hrfb_banded_se,
               prev_tbl$prev_hrfb_smooth_se),
  ci_low   = c(prev_tbl$prev_hrfb_raw_low,
               prev_tbl$prev_hrfb_banded_low,
               prev_tbl$prev_hrfb_smooth_low),
  ci_upp   = c(prev_tbl$prev_hrfb_raw_upp,
               prev_tbl$prev_hrfb_banded_upp,
               prev_tbl$prev_hrfb_smooth_upp)
) %>%
  mutate(ci_width = ci_upp - ci_low)

print(precision_tbl)
readr::write_csv(precision_tbl, file.path(tab_dir, "stage2_prev_precision_comparison.csv"))


# --------------------------
# 8) Threshold-focused sanity check
# --------------------------

# If the differences introduced by smoothing are concentrated very tightly around
# the age thresholds, that suggests the sensitivity exercise is mainly nudging
# boundary cases rather than materially changing the data structure.

near_thresh <- df %>%
  mutate(
    dist_18 = age_at_birth_years - 18,
    dist_35 = age_at_birth_years - 35
  ) %>%
  summarise(
    flips_young_within_0_5yr = mean(abs(dist_18) <= 0.5 &
                                      (risk_age_young_smooth != (age_at_birth_years < 18)), na.rm = TRUE) * 100,
    flips_old_within_0_5yr   = mean(abs(dist_35) <= 0.5 &
                                      (risk_age_old_smooth   != (age_at_birth_years > 34)), na.rm = TRUE) * 100
  )
print(near_thresh)
readr::write_csv(near_thresh, file.path(tab_dir, "stage2_threshold_flip_rates.csv"))

# --------------------------
# 9) Save enriched Stage 2 dataset
# --------------------------
stage2_path <- file.path(out_dir, "step02_births_qc_enhanced.rds")
saveRDS(df, stage2_path)

# --------------------------
# 10) Write Stage 2 README
# --------------------------

readme_path <- file.path(out_dir, "README_stage2_age_heaping.txt")
readme_txt <- c(
  "STEP 2 OUTPUTS — AGE HEAPING AND SPIKING DIAGNOSTICS",
  "",
  "Files:",
  "- figs_stage2/age_at_birth_histogram.pdf/.png",
  "- figs_stage2/age_at_birth_terminal_digit_preference.pdf/.png",
  "- figs_stage2/child_age_months_histogram.pdf/.png",
  "- figs_stage2/age_at_birth_histogram_raw_vs_smoothed.pdf/.png",
  "- figs_stage2/terminal_digit_preference_raw_vs_smoothed.pdf/.png",
  "- tables_stage2/stage2_age_integer_diagnostic.csv",
  "- tables_stage2/stage2_qc_summary_thresholds_and_spikes.csv",
  "- tables_stage2/stage2_heaping_indices.csv",
  "- tables_stage2/stage2_heaping_indices_before_after_smoothing.csv",
  "- tables_stage2/stage2_agreement_raw_vs_banded.csv",
  "- tables_stage2/stage2_agreement_raw_vs_smoothed.csv",
  "- tables_stage2/stage2_prev_sensitivity.csv",
  "- tables_stage2/stage2_prev_precision_comparison.csv",
  "- tables_stage2/stage2_threshold_flip_rates.csv",
  paste0("- ", basename(stage2_path)),
  "",
  "Notes:",
  "- age_at_birth is derived from DHS CMC dates, so the heaping indices are descriptive diagnostics only.",
  "- The banded and smoothed outcomes are used for sensitivity analysis rather than as replacement outcomes.",
  "- Weighted prevalence comparisons are based on a rebuilt survey design aligned to the enriched Stage 2 dataset.",
  "",
  paste0("Saved on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)

writeLines(readme_txt, readme_path)

message("Saved: ", normalizePath(stage2_path))
message("Wrote README: ", normalizePath(readme_path))
############################################################
