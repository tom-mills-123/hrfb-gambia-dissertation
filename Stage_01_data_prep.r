############################################################
# STAGE 1 — Data Preparation & Checkpoint
# Project: High-risk fertility behaviour (HRFB) — Gambia DHS 2019–20
#
# What this script does (and why):
#   1) Load BR (births) and IR (women) DHS files
#      - BR is the birth-level workhorse; IR contributes woman-level covariates and context.
#   2) Clean variable names and keep only key design + analysis fields
#      - reduces memory, speeds joins, and makes code readable.
#   3) Enforce unit of analysis: births in the 5 years before interview (0–59 months)
#      → improves temporal alignment of outcome and covariates, reduces recall bias.
#   4) Construct HRFB components per birth:
#      - maternal age at birth <18 or >34
#      - short interval (<24 months)
#      - high parity (birth order ≥4)
#   5) Create human-readable factors for core covariates and ordered FP message exposure.
#   6) Join selected IR covariates (v312 contraceptive use; v602 fertility preference).
#   7) Create a PSU contextual aggregate: share of women with ≥ secondary education (v001).
#   8) Build a survey design object (weights, PSU, strata) for design-based inference.
#   9) Quick QA snapshots: missingness and sanity checks.
#  10) Save analysis-ready objects (births df + survey design) as a checkpoint for later stages.
############################################################

# --------------------------
# 0) Packages + global options
# --------------------------

# Set to TRUE only if you want the script to install missing packages automatically.
AUTO_INSTALL <- FALSE

# List of packages required for import, cleaning, plotting, and complex survey inference
pkgs <- c(
  "haven","dplyr","tidyr","stringr","labelled","survey","srvyr",
  "janitor","forcats","ggplot2","purrr","Cairo"
)

# Identify which required packages are not installed
missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))

# Install missing packages only if explicitly requested; otherwise stop with instructions
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

# Load packages (suppress startup messages to keep logs readable)
suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(labelled)
  library(survey)
  library(srvyr)
  library(janitor)
  library(forcats)
  library(ggplot2)
  library(purrr)
  library(Cairo)
})

# Configure variance estimation when a stratum contains only one PSU after subsetting
options(survey.lonely.psu = "adjust")


# --------------------------
# Plot style + export helpers
# --------------------------

# Define a consistent publication-ready plotting theme
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

# Save figures in both PDF (vector) and PNG (high resolution)
save_fig <- function(p, filename, fig_dir, w = 7, h = 4.5) {
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  
  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".pdf")),
    plot = p, width = w, height = h, units = "in",
    device = Cairo::CairoPDF
  )
  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".png")),
    plot = p, width = w, height = h, units = "in", dpi = 450
  )
}

# --------------------------
# 1) Locate and import DHS data (BR & IR)
#    Why: BR has birth histories; IR has woman-level covariates/context.
# --------------------------


# --------------------------
# 1) Locate and import DHS data (BR & IR)
# --------------------------

# Set directory
data_dir <- "PATH_TO_YOUR_DHS_DATA"

if (!dir.exists(data_dir)) {
  stop("Please set `data_dir` to the local folder containing the DHS BR and IR .dta files.")
}

# Identify Births Recode (BR) and Individual Recode (IR) files
br_path <- list.files(data_dir, pattern = "(?i)BR.*\\.dta$", full.names = TRUE, recursive = TRUE)
ir_path <- list.files(data_dir, pattern = "(?i)IR.*\\.dta$", full.names = TRUE, recursive = TRUE)
if (length(br_path) == 0) stop("Births Recode (BR) .dta file not found.")
if (length(ir_path) == 0) stop("Individual Recode (IR) .dta file not found.")
if (length(br_path) > 1) warning("Multiple BR files found; using: ", br_path[1])
if (length(ir_path) > 1) warning("Multiple IR files found; using: ", ir_path[1])
br_path <- br_path[1]; ir_path <- ir_path[1]

# Import and standardise variable names
br_raw <- read_dta(br_path) %>% clean_names()
ir_raw <- read_dta(ir_path) %>% clean_names()
message("BR dims: ", paste(dim(br_raw), collapse = " x"))
message("IR dims: ", paste(dim(ir_raw), collapse = " x"))

# Output directories (used for figures + checkpoints)
out_dir <- "outputs"
fig_dir <- file.path(out_dir, "figs_stage1")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)


# --------------------------
# 2) Keep key variables & set initial design helpers
#    Why: Trim to essentials; create weight/strata fields for later.
# --------------------------
# Start from the raw Births Recode data (already read and name-cleaned as `br_raw`)
br <- br_raw %>%
  # Keep only the columns we need for design, core outcomes, and covariates
  select(
    caseid,                 # Case ID (unique within file; helps with tracking and QA)
    v001, v002, v003,       # IDs: v001 = cluster/PSU; v002 = household; v003 = respondent (woman)
    
    v021,                   # Survey design: Primary Sampling Unit (PSU) used by survey package
    v022,                   # Survey design: Strata (preferred) — may be missing in some releases
    v023,                   # Survey design: Alternate strata variable (fallback if v022 is NA)
    v005,                   # Survey weight; stored as an integer that must be scaled by 1e6
    v011,                   # Woman's date of birth in CMC (Century Month Code)
    v012,                   # Woman's age at interview in completed years (quick age check/QA)
    v024,                   # Region (admin unit used for disaggregation/maps)
    v025,                   # Place of residence (Urban/Rural)
    v106,                   # Highest education level attained
    v190,                   # Wealth quintile (1=poorest ... 5=richest)
    v131,                   # Ethnicity (country-specific categories)
    v157, v158, v159,       # Exposure to family planning messages: radio (v157), TV (v158), print (v159)
    v008,                   # Interview date in CMC (used to compute fallback child age if b19 missing)
    b3,                     # Child's date of birth in CMC (used for age-at-birth and age computations)
    b11,                    # Preceding birth interval in months (used to flag short intervals <24)
    bord,                   # Birth order (1,2,3,...) — used to flag high parity (≥4)
    b19                     # Child's age in months at time of interview (used to restrict to last 5 years)
  ) %>%
  # Create standardized design helpers used later by the survey design object
  mutate(
    strata = coalesce(v022, v023),  # If v022 exists use it; otherwise use v023 as the strata variable
    weight = v005 / 1e6             # DHS weight must be divided by 1,000,000 to get the correct scale
  )



# --------------------------
# 3) Restrict births to the last 5 years
#    Why: Align outcomes (e.g. birth timing, risk factors)
#         with covariates measured at interview.
#         DHS collects full birth histories, but we only
#         keep births within 5 years to reduce recall bias
#         and maintain temporal consistency.
# --------------------------

br <- br %>%
  mutate(
    # Create a new variable 'child_age_months' representing the child's age
    # at the time of the mother's interview.
    #
    # - Some DHS BR files already contain b19 (age in months).
    # - If b19 is missing (NA), we compute it manually using v008 (interview date)
    #   minus b3 (child's date of birth), both stored in Century Month Code (CMC).
    #   The result is age in months.
    #
    # coalesce(x, y) returns the first non-missing value between x and y.
    # So this expression says: "use b19 if available; otherwise use v008 - b3."
    child_age_months = coalesce(b19, v008 - b3)
  ) %>%
  
  # Filter the dataset to include only births that occurred
  # within the last 60 months (i.e., 5 years) prior to interview.
  #
  # - Exclude rows where 'child_age_months' is missing (NA),
  #   as those cases have insufficient timing data.
  # - Exclude births with child_age_months >= 60, i.e., older than 5 years.
  #
  # This ensures that all births analyzed are recent relative
  # to the mother's interview and therefore comparable across covariates.
  filter(!is.na(child_age_months) & child_age_months < 60)

# Print a message to the console showing how many rows remain
# after the 5-year restriction (useful QA to confirm the subset size).
message("Rows after 5-year restriction: ", nrow(br))

# There are now 8362 rows after restriction

# Double-check that the 5-year restriction worked by looking at the distribution 
# of child ages in months

br %>%
  summarise(
    min_age = min(child_age_months, na.rm = TRUE),
    max_age = max(child_age_months, na.rm = TRUE),
    mean_age = mean(child_age_months, na.rm = TRUE)
  )


# Visual QC: distribution of child ages (months) to confirm 0–59 restriction behaves as expected
p_child_age_stage1 <- ggplot(br, aes(x = child_age_months)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left",
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = seq(0, 59, by = 12), linetype = "dotted", linewidth = 0.3) +
  scale_x_continuous(breaks = seq(0, 60, by = 6), limits = c(0, 59)) +
  labs(
    title = "Child age in months (analytic 5-year window)",
    subtitle = "Unweighted QC check of 0–59 month restriction (dotted lines at 12-month intervals)",
    x = "Child age (months)",
    y = "Number of births",
    caption = "Data: Gambia DHS 2019–20 BR file; ages use b19 with v008-b3 fallback if missing."
  ) +
  theme_dissertation() +
  theme(axis.text.x = element_text(size = 9))

print(p_child_age_stage1)
save_fig(p_child_age_stage1, "stage1_child_age_months_distribution", fig_dir = fig_dir)

# QA: check for potential top-coding/heaping at the upper bound
br %>%
  count(child_age_months) %>%
  arrange(desc(n)) %>%
  dplyr::slice_head(n = 10) %>%
  print()


br %>%
  summarise(
    pct_at_59 = mean(child_age_months == 59, na.rm = TRUE) * 100,
    pct_at_58 = mean(child_age_months == 58, na.rm = TRUE) * 100
  ) %>%
  print()


# --------------------------
# 4) Maternal age at birth (years)
#    Why: Defines age-based HRFB components.
# --------------------------
br <- br %>%
  mutate(
    # Compute mother's age (in years) at the time of each birth.
    #
    # DHS stores all dates as "Century Month Code" (CMC) —
    #   the number of months since January 1900.
    # This makes it easy to compute time differences in months.
    #
    # b3   = child's date of birth (CMC)
    # v011 = mother's date of birth (CMC)
    #
    # Subtracting gives the mother's age in MONTHS at the child's birth.
    # Dividing by 12 converts months → years.
    #
    # Example:
    #   If v011 = 1400 (born in April 2016) and b3 = 1700 (April 1941)
    #   Real example: v011 = 1650 (born June 1987), b3 = 2100 (June 2012)
    #   → 2100 - 1650 = 450 months = 37.5 years
    # So age_at_birth = 37.5
    #
    # Result is a numeric variable (may include decimals).
    age_at_birth = (b3 - v011) / 12
  )


# --------------------------
# 5) HRFB components & summaries
#    Why: Construct the outcome variables used throughout the study.
# --------------------------
br <- br %>%
  mutate(
    # 1) Young maternal age at birth:
    #    Flag births where mother's age at delivery was strictly < 18 years.
    #    - age_at_birth is continuous (from Step 4).
    #    - Using 1L/0L stores integers (handy for modeling & memory).
    #    - missing = 0L ensures NA ages do not get flagged as "at risk" by accident.
    risk_age_young = if_else(age_at_birth < 18, 1L, 0L, missing = 0L),
    
    # 2) Older maternal age at birth:
    #    Flag births where mother's age at delivery was strictly > 34 years.
    #    - Exactly 35 qualifies as risk (because >34).
    risk_age_old   = if_else(age_at_birth > 34, 1L, 0L, missing = 0L),
    
    # 3) Short preceding birth interval:
    #    Flag births if the interval since the previous birth (b11) is < 24 months.
    #    - b11 is NA for first births by DHS design, so require !is.na(b11)
    #      to avoid incorrectly labeling first births as "short interval."
    risk_short_int = if_else(!is.na(b11) & b11 < 24, 1L, 0L, missing = 0L),
    
    # 4) High parity (birth order ≥ 4):
    #    Use 'bord' as a proxy for parity at the time of the birth.
    #    - bord is 1 for first births, 2 for second, etc.
    risk_high_parity = if_else(!is.na(bord) & bord >= 4, 1L, 0L, missing = 0L),
    
    # 5) Total number of risks present (0 to 4):
    #    Sum the four binary indicators to capture cumulative risk load.
    hrfb_risk_count = risk_age_young + risk_age_old + risk_short_int + risk_high_parity,
    
    # 6) Primary study outcome: any HRFB (binary):
    #    1 if at least one risk present, else 0.
    hrfb_binary = if_else(hrfb_risk_count >= 1, 1L, 0L),
    
    # 7) Secondary outcome: multinomial risk classification.
    #    - "none"    : 0 risks
    #    - "single"  : exactly 1 risk
    #    - "multiple": 2 or more risks (captures compounding effects)
    hrfb_multicat = case_when(
      hrfb_risk_count == 0 ~ "none",
      hrfb_risk_count == 1 ~ "single",
      hrfb_risk_count >= 2 ~ "multiple",
      TRUE ~ NA_character_   # defensive default (should rarely trigger)
    )
  )

# Check the variables behave as expected
# Summarize proportions (unweighted preview)
br %>%
  summarise(
    pct_young = mean(risk_age_young) * 100,
    pct_old = mean(risk_age_old) * 100,
    pct_short = mean(risk_short_int) * 100,
    pct_parity = mean(risk_high_parity) * 100,
    pct_any = mean(hrfb_binary) * 100
  )

# Cross-tab to see how many births fall into each risk category
table(br$hrfb_multicat, useNA = "ifany")

# --------------------------
# 6) Human-readable covariates (labels → factors)
#    Why: Easier QA/plots; FP as ordered to preserve frequency (dose-response).
# --------------------------

# Define the desired order of FP exposure categories from lowest → highest frequency.
# Using a fixed order ensures consistent plotting and interpretation across variables.
fp_levels <- c("not at all", "less than once a week", "at least once a week", "almost every day")

br <- br %>%
  mutate(
    # Convert DHS-coded integers into their labelled text equivalents.
    # haven::as_factor() reads the value labels embedded in the imported Stata file.
    residence     = haven::as_factor(v025),  # Urban / Rural
    region        = haven::as_factor(v024),  # Administrative regions used by DHS
    educ          = haven::as_factor(v106),  # Highest education level attained
    wealth        = haven::as_factor(v190),  # Wealth quintile (Poorest → Richest)
    
    # Keep an explicit copy of the raw-labelled ethnicity so we can lump later
    # but still have the pre-lumped version for descriptives if needed.
    ethnicity_raw = haven::as_factor(v131),
    
    # Family planning message exposure variables:
    # v157 = Radio, v158 = TV, v159 = Print (newspaper/magazine).
    # We convert to ordered factors so that the categories reflect increasing exposure.
    # intersect() guards against surveys that have slightly different label spellings,
    # preventing invalid-level warnings by using only the levels actually present.
    fp_radio = ordered(
      haven::as_factor(v157),
      levels = intersect(fp_levels, levels(haven::as_factor(v157)))
    ),
    fp_tv = ordered(
      haven::as_factor(v158),
      levels = intersect(fp_levels, levels(haven::as_factor(v158)))
    ),
    fp_print = ordered(
      haven::as_factor(v159),
      levels = intersect(fp_levels, levels(haven::as_factor(v159)))
    )
  ) %>%
  mutate(
    # Ethnicity can have many small groups. To avoid tiny cells in cross-tabs and models,
    # lump categories that constitute <2% of the sample into "Other".
    # Ethnicity: standardise case first, then lump small groups
    ethnicity = ethnicity_raw %>%
      as.character() %>%
      stringr::str_trim() %>%
      stringr::str_to_lower() %>%
      forcats::as_factor() %>%
      forcats::fct_recode(Other = "other") %>%   # if "other" exists, force it to "Other"
      forcats::fct_lump_prop(prop = 0.02, other_level = "Other") %>%
      forcats::fct_drop()
  )

# Run QA checks
# 1) Verify factor levels (and ordering for FP)
levels(br$residence); levels(br$region); levels(br$educ); levels(br$wealth); levels(br$ethnicity)
levels(br$fp_radio);  is.ordered(br$fp_radio)
levels(br$fp_tv);     is.ordered(br$fp_tv)
levels(br$fp_print);  is.ordered(br$fp_print)

# 2) Show FP exposure distributions (include zero-count levels to see what's missing)
br %>% count(fp_radio, .drop = FALSE)
br %>% count(fp_tv,    .drop = FALSE)
br %>% count(fp_print, .drop = FALSE)

# 3) Make missing explicit (good for scan of data quality)
br %>%
  mutate(fp_tv = forcats::fct_explicit_na(fp_tv, na_level = "Missing")) %>%
  count(fp_tv, .drop = FALSE)

# 4) Ethnicity before vs after lumping (sanity check)
br %>% count(ethnicity_raw, sort = TRUE) %>% print(n = 20)
br %>% count(ethnicity,      sort = TRUE) %>% print(n = 20)



# --------------------------
# 7) Join extra IR variables & build contextual aggregate
#    Why:
#      1. Add woman-level variables from the IR file that aren't in BR
#         (v312 = contraceptive use, v602 = fertility preference).
#      2. Create a cluster-level (% with secondary education or higher)
#         to capture contextual educational environment.
# --------------------------

# 7a) --- Bring selected IR covariates ---
ir_small <- ir_raw %>%
  select(
    v001, v002, v003,  # Keys: cluster, household, respondent (to join on)
    v312,              # Current contraceptive method
    v602,              # Fertility preference ("wants more children", etc.)
    v106               # Education level (kept for cluster-level aggregate)
  ) %>%
  distinct()            # Ensure one row per woman (avoid duplicates)

# Join these IR variables to the BR dataset.
# Each birth will inherit its mother's contraceptive and fertility preference info.
df <- br %>%
  left_join(ir_small, by = c("v001", "v002", "v003")) %>%
  mutate(
    # >>> ADDED: stable unique identifier for each woman <<<
    # Why: used later for within-woman sensitivity analysis (one birth per woman).
    woman_id = paste(v001, v002, v003, sep = "_"),
    
    # Convert DHS coded variables into labelled factors.
    # haven::as_factor() uses the embedded value labels to turn codes into text.
    contraceptive_use = haven::as_factor(v312),  # e.g. "No method", "Pill", "Injectable", etc.
    fert_pref         = haven::as_factor(v602)   # e.g. "Wants within 2 years", "Wants no more"
  )


# 7b) --- Contextual aggregate: share of women with ≥ secondary education per PSU (v001) ---
clus_educ <- ir_raw %>%
  # Create a binary indicator: 1 if education level >= "secondary" (v106 >= 2)
  # DHS coding for v106:
  # 0 = no education, 1 = primary, 2 = secondary, 3 = higher
  transmute(
    v001,
    sec_plus = as.integer(as.numeric(v106) >= 2)
  ) %>%
  # Group by cluster (PSU)
  group_by(v001) %>%
  # Compute the proportion of women with secondary+ education in each PSU
  summarise(
    share_sec_plus = mean(sec_plus, na.rm = TRUE),
    .groups = "drop"
  )

# Merge this cluster-level contextual variable back into the birth-level data.
# Every birth now carries information about the educational environment of its cluster.
df <- df %>% left_join(clus_educ, by = "v001")

# Step 7 QA checks
# Check row count (should not change)
nrow(br); nrow(ir_small); nrow(df)

# Confirm successful join (no lost or duplicated rows)
df %>% summarise(n_births = n(), n_missing_contra = sum(is.na(v312)))

# Check unique women in IR subset
ir_small %>% count(v001, v002, v003) %>% summarise(max_n = max(n))  # should return 1

# Sanity-check contextual variable: should range from 0 to 1
summary(df$share_sec_plus)

# PSU-level distribution (one row per PSU)
psu_dist <- df %>%
  distinct(v001, share_sec_plus)

p_share_sec_plus <- ggplot(psu_dist, aes(x = share_sec_plus)) +
  geom_histogram(bins = 20, boundary = 0, closed = "left",
                 color = "white", linewidth = 0.2) +
  labs(
    title = "Cluster educational context",
    subtitle = "Distribution of PSU-level share of women with ≥ secondary education",
    x = "Share with ≥ secondary education (0–1)",
    y = "Number of PSUs",
    caption = "Computed from IR file at PSU (v001) level; one observation per PSU."
  ) +
  theme_dissertation()

print(p_share_sec_plus)
save_fig(p_share_sec_plus, "stage1_psu_share_sec_plus_distribution", fig_dir = fig_dir)




# ---------------------------------------------------------
# STEP 8 — Build the survey design object (DHS complex design)
#
# Purpose:
#   Turn the flat data frame `df` into a design-aware object that
#   encodes DHS sampling features (weights, PSUs, strata). This is
#   required for correct weighted estimates, standard errors, and CIs.
#
# Inputs expected in `df`:
#   - weight : numeric = v005 / 1e6 (scaled DHS sampling weight)
#   - v021   : integer PSU/cluster identifier
#   - strata : stratification variable (coalesce(v022, v023))
#
# Output:
#   - design_srvyr : an srvyr survey design object
# ---------------------------------------------------------

# 1) Safety check — ensure required variables exist before building the design
stopifnot(all(c("weight", "strata", "v021") %in% names(df)))

# 2) Filter out any rows lacking design info
#    DHS files are usually complete, but this guard prevents NA-related errors later.
df_design <- df %>%
  dplyr::filter(!is.na(weight), !is.na(v021), !is.na(strata))

# 3) Build the survey design object
#    - ids    : PSU/cluster identifier (v021)
#    - strata : stratification variable (often region x urban/rural)
#    - weights: probability weight (v005 scaled to weight = v005 / 1e6)
#    - nest   : tells the survey package PSUs are nested within strata (DHS standard)
design_srvyr <- srvyr::as_survey_design(
  .data   = df_design,
  ids     = v021,
  strata  = strata,
  weights = weight,
  nest    = TRUE
)

# 4) QA Checks

# How many observations, PSUs, and strata made it into the design?
df_design %>%
  dplyr::summarise(
    n_rows   = dplyr::n(),
    n_psu    = dplyr::n_distinct(v021),
    n_strata = dplyr::n_distinct(strata),
    w_sum    = sum(weight, na.rm = TRUE),
    w_min    = min(weight, na.rm = TRUE),
    w_max    = max(weight, na.rm = TRUE)
  ) %>%
  print()

# Smoke test: compute a weighted prevalence with a 95% CI
# (Replace `hrfb_binary` with any 0/1 variable you’ve created.)
design_srvyr %>%
  srvyr::summarise(hrfb_prev = survey_mean(hrfb_binary, vartype = "ci", na.rm = TRUE)) %>%
  print()


# ---------------------------------------------------------
# STEP 9 — Core QA snapshots
#
# Purpose:
#   1) Verify completeness of key variables used to define outcomes
#      and the survey design (weights, PSU, strata).
#   2) Do a quick design-adjusted prevalence check for HRFB to confirm
#      the survey design object behaves as expected.
#
# Assumes:
#   - `df` exists (your cleaned birth-level dataset after Steps 1–8).
#   - `design_srvyr` exists (survey design object built in Step 8).
#
# Outputs:
#   - A tibble `qa_missing` with missing-value counts.
#   - A printed weighted prevalence of HRFB (with 95% CI).
#   - Optional sanity checks to confirm b11 missingness matches first births.
# ---------------------------------------------------------

# 1) Snapshot of missingness in critical fields
qa_missing <- df %>%
  dplyr::summarise(
    n                 = dplyr::n(),                 # total births in analytic sample
    miss_age_at_birth = sum(is.na(age_at_birth)),   # should be 0 (computed in Step 4)
    miss_b11          = sum(is.na(b11)),            # expected missing for first births (no preceding interval)
    miss_bord         = sum(is.na(bord)),           # should be 0
    miss_weight       = sum(is.na(weight)),         # should be 0 (else design/weights break)
    miss_psu          = sum(is.na(v021)),           # should be 0
    miss_strata       = sum(is.na(strata))          # should be 0
  )

# Print the missingness table
message("\n== Step 9: Missingness snapshot ==")
print(qa_missing)

# 2) Quick design-adjusted prevalence of any HRFB (binary)
#    Verifies the survey design object is working and result
#    is consistent with earlier checks (~0.54 for GMB 2019–20).
message("\n== Step 9: Weighted prevalence of any HRFB (binary) ==")
design_srvyr %>%
  srvyr::summarise(
    prev_hrfb = survey_mean(hrfb_binary, vartype = "ci", na.rm = TRUE)
  ) %>%
  print()

# 3) Sanity check: most missing b11 should be first births (bord == 1).
#    This confirms missing intervals are by design not data loss.
message("\n== Are missing b11 values basically first births? ==")
df %>%
  dplyr::summarise(
    total_b11_na = sum(is.na(b11)),
    first_births = sum(bord == 1, na.rm = TRUE),
    overlap      = sum(is.na(b11) & bord == 1, na.rm = TRUE)
  ) %>%
  print()

# 4) Bounds check on age_at_birth — catches data entry outliers.
message("\n== Logical bounds check for age_at_birth (years) ==")
df %>%
  dplyr::summarise(
    n_lt_10 = sum(age_at_birth < 10,  na.rm = TRUE),
    n_gt_55 = sum(age_at_birth > 55,  na.rm = TRUE)
  ) %>%
  print()

# Notes for interpretation:
# - miss_age_at_birth, miss_bord, miss_weight, miss_psu, miss_strata should be 0.
# - miss_b11 should roughly equal the number of first births (bord == 1).
# - prev_hrfb should be close to the smoke test (around 0.54 with a tight CI).
# - n_lt_10 and n_gt_55 should both be 0 (or extremely rare) in cleaned DHS data.


# ---------------------------------------------------------
# STEP 10 — Save analysis-ready objects (checkpoint)
#
# Purpose:
#   Persist the cleaned births dataset (`df`) and the survey design
#   object (`design_srvyr`) so downstream scripts can start from here
#   without re-running Steps 1–9. Also write a small README describing
#   what's inside and key assumptions.
# ---------------------------------------------------------


# 2) Construct clear, versionable file paths for the checkpoint objects
births_path <- file.path(out_dir, "step01_births_clean.rds")      # cleaned birth-level dataframe
design_path <- file.path(out_dir, "step01_svy_design.rds")        # srvyr survey design object

# 3) Save the R objects as .rds 
saveRDS(df,           births_path)                     # write the cleaned births dataset (5-year window)
saveRDS(design_srvyr, design_path)                     # write the survey design (weights/PSU/strata)

# 4) Write a README file
readme_path <- file.path(out_dir, "README_step1_checkpoint.txt")  
readme_txt <- c(
  "STEP 1 CHECKPOINT",
  "",
  "Files:",
  paste0("- ", basename(births_path), ": Cleaned birth-level dataset (5-year restriction). ",
         "Includes HRFB flags, socio-demographic factors, IR covariates (v312, v602), ",
         "and PSU contextual share_sec_plus."),
  paste0("- ", basename(design_path), ": srvyr survey design object (v021 PSU, strata coalesce(v022,v023), weight v005/1e6)."),
  "",
  "Notes:",
  "- HRFB: age<18, age>34, interval<24 months, parity>=4.",
  "- Unit of analysis: births with b19<60 (or v008-b3 fallback).",
  "- Contextual: share_sec_plus = proportion women with >= secondary education by PSU (IR).",
  "- Sensitivity: refit primary binary model using one birth per woman (most recent birth) to assess within-woman clustering.",
  "",
  paste0("Saved on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
writeLines(readme_txt, readme_path)                   

# 5) Console confirmations (nice for logs / debugging paths)
message("✓ Saved: ", normalizePath(births_path))
message("✓ Saved: ", normalizePath(design_path))
message("✓ Wrote README: ", normalizePath(readme_path))
# ---------------------------------------------------------
