# Higher SARS-CoV-2 Transmission Burden Among Racialized Individuals: Evidence from Canadian Serology Data

**Authors:** Simran Mann, Natalie Wilson, Clara Lee, David N. Fisman

**Affiliation:** Dalla Lana School of Public Health, University of Toronto

**Corresponding author:** David Fisman, MD MPH (david.fisman@utoronto.ca)

## Overview

This repository contains the analysis code for our study examining racial disparities in SARS-CoV-2 force of infection (FOI) in Canada using seroprevalence data from the Canadian Blood Services serosurveillance program. We fit a dynamic susceptible-infected (SI) compartmental model with seroreversion to race-stratified seroprevalence data (June 2020 – April 2023) to estimate FOI separately for White and racialized populations in the pre-Omicron and Omicron periods.

**Key finding:** Racialized individuals experienced significantly higher force of infection throughout the pandemic. Although the racial disparity narrowed during the Omicron period (IRR from 2.21 to 1.24), this compression was driven by a larger proportional increase in FOI among White individuals—not by reduced transmission among racialized populations. Apparent seroprevalence convergence does not reflect epidemiologic equity.

## Data

Seroprevalence data used in this study are publicly available from the [COVID-19 Immunity Task Force (CITF) Data Portal](https://portal.citf.mcgill.ca).

Place the downloaded data file (`seroprevalence_by_social_determinant.csv`) in the `data/` directory before running the analysis.

**Note:** One erroneous observation is excluded by the analysis scripts (White group, December 28, 2021; reported seroprevalence 89.6%). The anti-N and anti-S values for this observation appear to have been transposed in the source data (surrounding anti-N values were approximately 5%).

## Repository Structure

```
covid19-seroprev-race/
├── R/
│   ├── 00_helpers.R              # Shared functions (SI model, delta method CIs)
│   ├── 01_saturated_model.R      # Fit saturated model (race-specific Omicron effects)
│   ├── 02_parsimonious_model.R   # Fit parsimonious model (shared Omicron multiplier)
│   ├── 03_compare_models.R       # Model comparison (AIC, BIC, LR test)
│   ├── 04_irr_with_cis.R        # Incidence rate ratios with 95% CIs
│   ├── 05_multipliers_with_cis.R # Omicron multipliers with 95% CIs
│   └── 06_figure1.R             # Publication figure (BMC format)
├── data/
│   └── (place seroprevalence_by_social_determinant.csv here)
├── output/
│   └── figures/
├── README.md
└── .gitignore
```

## How to Run

Scripts are numbered and should be run in order from the project root directory. All scripts source `R/00_helpers.R` for shared functions.

```r
# From the project root directory:
source("R/01_saturated_model.R")      # ~1 min
source("R/02_parsimonious_model.R")   # ~1 min
source("R/03_compare_models.R")       # seconds
source("R/04_irr_with_cis.R")        # seconds
source("R/05_multipliers_with_cis.R") # seconds
source("R/06_figure1.R")             # seconds
```

Model objects are saved as `.rds` files in `output/` after steps 1–2, and subsequent scripts load these. Results tables are saved as `.csv` files in `output/`.

## Dependencies

- R ≥ 4.4.3
- `bbmle` — maximum likelihood estimation
- `dplyr` — data manipulation
- `lubridate` — date parsing
- `readr` — CSV reading
- `ggplot2` — figure generation

Install all dependencies with:

```r
install.packages(c("bbmle", "dplyr", "lubridate", "readr", "ggplot2"))
```

## Model Description

We use a discrete-time SI model with seroreversion:

- S(t+1) = S(t) × exp(−λ(t)) + I(t) × (1 − exp(−ρ))
- I(t+1) = I(t) × exp(−ρ) + S(t) × (1 − exp(−λ(t)))

where λ(t) is the time-varying force of infection and ρ is the seroreversion rate. The model is initialized at March 2020 (S=1, I=0) and fit to observed seroprevalence data starting June 2020 using maximum likelihood estimation with binomial likelihood.

**Saturated model** (5 parameters): Race-specific FOI for pre-Omicron and Omicron periods (4 λ) plus a shared seroreversion rate (ρ).

**Parsimonious model** (4 parameters): Race-specific pre-Omicron FOI (2 λ), a shared Omicron multiplier, and ρ.

The saturated model was selected based on superior fit (ΔAIC = −551.2, LR test p = 2.53 × 10⁻¹²²).

## Citation

If you use this code, please cite:

> Mann S, Wilson N, Lee C, Fisman DN. Higher SARS-CoV-2 Transmission Burden Among Racialized Individuals: Evidence from Canadian Serology Data. *BMC Public Health*. [year];[vol]:[pages].

This work builds on:

> Hassan A, Nassrallah EI, Fisman DN. Seroprevalence Convergence Does Not Reflect Transmission Equity: Persistent Socioeconomic Disparities in COVID-19 Force of Infection in Canada. *Epidemics*. 2026.

## Funding

This project was supported by a Canadian Institutes for Health Research (CIHR) project grant (#518192).

## License

MIT License. See LICENSE file.
