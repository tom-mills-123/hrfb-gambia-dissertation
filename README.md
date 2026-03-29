# hrfb-gambia-dissertation

MSc dissertation code for analysing **high-risk fertility behaviour (HRFB)** in **The Gambia** using **DHS 2019–20** data.

## Project overview

This repository contains the full analysis pipeline used for my dissertation on high-risk fertility behaviour (HRFB) in The Gambia. The project uses Demographic and Health Survey (DHS) data and combines data preparation, survey-weighted descriptive analysis, regression modelling, machine learning, and exploratory scenario modelling.

The analysis is conducted at the **birth level**, with the analytic sample restricted to **births occurring in the 5 years before interview**.

## Research focus

The project examines patterns and predictors of high-risk fertility behaviour (HRFB) in The Gambia using DHS 2019–20 data.

HRFB is defined using standard demographic risk components:
- maternal age at birth **under 18 years**
- maternal age at birth **over 34 years**
- **short birth interval** (<24 months)
- **high parity** (birth order 4 or higher)

## Data source

This project uses data from:

- **The Gambia DHS 2019–20**
- DHS recode files used:
  - **BR file** (Births Recode)
  - **IR file** (Individual Recode)

## Important note on data access

The DHS datasets are **not included** in this repository.

To run the code, users must:
1. Register for access via the DHS Program
2. Download the required DHS data files
3. Update the local file path in `Stage_01_data_prep.r`

DHS data are available from the DHS Program website.

## Project structure

The analysis is organised into sequential stages. Each script builds on outputs from the previous stage.

### `Stage_01_data_prep.r`
Data preparation and checkpoint creation.

This script:
- loads the DHS BR and IR files
- cleans and standardises variables
- restricts the sample to births in the last 5 years
- constructs the HRFB indicators
- creates core covariates
- builds the survey design variables
- saves analysis-ready objects for later stages

### `Stage_02_data_quality_and_age_spiking.r`
Data quality assessment and sensitivity checks.

This script:
- evaluates potential age heaping and age spiking
- produces diagnostic plots
- calculates age-quality summary measures
- explores banded and smoothed age-based sensitivity checks
- saves enhanced working datasets for downstream analysis

### `Stage_03_descriptive_analysis.r`
Descriptive analysis for Research Question 1.

This script:
- estimates national HRFB prevalence
- estimates subgroup prevalence by key covariates
- applies survey-weighted methods
- produces tables and figures for descriptive results

### `Stage_04_regression_modelling.r`
Regression modelling for Research Question 2.

This script:
- fits survey-adjusted logistic regression models
- estimates adjusted associations with HRFB
- exports odds ratio tables
- produces forest plots and adjusted prevalence plots

### `Stage_05_machine_learning.r`
Machine-learning classification for Research Question 3.

This script:
- builds predictive models for HRFB classification
- fits Random Forest and XGBoost models
- evaluates model performance on held-out data
- produces ROC curves, calibration plots, confusion matrices, and feature-importance outputs

### `Stage_06_scenario_modelling_community.r`
Exploratory scenario modelling.

This script:
- refits the primary Stage 4 model
- simulates policy-relevant counterfactual scenarios
- estimates predicted HRFB prevalence under alternative scenarios
- produces overall and subgroup scenario outputs

### `Stage_06A_individual.R`
Exploratory sensitivity scenario modelling.

This script:
- applies an alternative individual-level scenario design
- evaluates sensitivity of the scenario results
- produces additional tables and figures for comparison with Stage 6

## How to run the analysis

Run the scripts in the following order:

1. `Stage_01_data_prep.r`
2. `Stage_02_data_quality_and_age_spiking.r`
3. `Stage_03_descriptive_analysis.r`
4. `Stage_04_regression_modelling.r`
5. `Stage_05_machine_learning.r`
6. `Stage_06_scenario_modelling_community.r`
7. `Stage_06A_individual.R`

## Required R packages

The scripts use a range of R packages for:
- data import and cleaning
- survey analysis
- regression modelling
- machine learning
- plotting and export

Each script includes its own package checks at the top. If a required package is missing, the script will stop and report which packages need to be installed.

## Outputs

Outputs are saved to the `outputs/` directory created during the workflow. These include:

- cleaned intermediate datasets (`.rds`)
- tables for each analysis stage
- figures in `.pdf` and `.png` formats
- model summaries and diagnostic outputs

Typical output folders include:
- `outputs/figs_stage1`
- `outputs/figs_stage2`
- `outputs/figs_stage3`
- `outputs/figs_stage4`
- `outputs/figs_stage5`
- `outputs/tables_stage2`
- `outputs/tables_stage3`
- `outputs/tables_stage4`
- `outputs/tables_stage5`

## Methodological notes

- All survey-based analyses account for **sampling weights**, **clustering (PSU)**, and **stratification**.
- The primary descriptive and regression analyses use the **raw HRFB outcome definition**.
- Variables that directly define HRFB are excluded from predictive models where appropriate to avoid circularity or leakage.
- The machine-learning stage uses grouped validation at PSU level to reduce overly optimistic performance estimates.
- The scenario-modelling stages are **exploratory** and should be interpreted as model-based simulations rather than causal effects.

## Repository purpose

This repository is intended to document and organise the code used for the dissertation analysis. It is structured to make the workflow transparent, reproducible, and easy to follow.

## Author

Tom

## Dissertation title

**High-risk fertility behaviour in The Gambia (DHS 2019–20)**

## License

This repository is for academic dissertation use. Reuse of DHS data remains subject to DHS Program data access conditions.