# MODEL STRUCTURE:
# - 5 parameters total:
#   * 2 pre-Omicron λ (White & Racialized)
#   * 2 Omicron λ (White & Racialized)  
#   * 1 seroreversion rate ρ
# - SI compartmental model with seroreversion
# - Model starts March 1, 2020 (S=1, I=0)
# - Fits to data from June 14, 2020 onwards (when race data is available)

# Set working directory to project root (adjust as needed)
# setwd("path/to/covid19-seroprev-race")

library(readr)
library(bbmle)
library(dplyr)
library(lubridate)
library(ggplot2)

source("R/00_helpers.R")

# FITTING SATURATED MODEL: RACE-SPECIFIC OMICRON EFFECTS --------------
# 1. LOAD AND PREPARE DATA 
ses_data <- read.csv("data/seroprevalence_by_social_determinant.csv")
ses_data$date <- mdy(ses_data$samplingdate)

table(ses_data$population)
table(ses_data$ab_target)
table(ses_data$ab_target, ses_data$population, exclude=NULL) # No missing values

# Filter to CBS, anti-N, race data
cbs_data <- ses_data %>%
  filter(ab_target == "N",
         population %in% c("white", "racialized"),
         project == "CBS") %>%
  arrange(date)

# Exclude erroneous observation: Dec 28 2021, white, seroprev=89.6%
# This observation has anti-N and anti-S values that appear transposed
# (89.6% anti-N vs 10.4% anti-S; surrounding anti-N values are ~5%)
n_before <- nrow(cbs_data)
cbs_data <- cbs_data %>%
  filter(!(population == "white" & date == as.Date("2021-12-28") & seroprev_est > 80))
cat(sprintf("Excluded %d erroneous observation(s)\n", n_before - nrow(cbs_data)))

cat("Data loaded successfully\n")
cat(sprintf("Total observations: %d\n", nrow(cbs_data)))
cat(sprintf("Date range: %s to %s\n\n", min(cbs_data$date), max(cbs_data$date)))

# 2. DEFINE KEY DATES
model_start_date <- as.Date("2020-03-01")  # Model starts here (S=1, I=0)
omicron_date <- as.Date("2021-12-26")      # Omicron emergence

cat(sprintf("Model start date: %s (S=1, I=0)\n", model_start_date))
cat(sprintf("Omicron emergence: %s\n\n", omicron_date))

# 3. PREPARE DATA FOR LIKELIHOOD CALCULATION
race <- c("white", "racialized")
all_observed <- list()
all_dates <- list()
all_sample_sizes <- list()

for (i in 1:2) {
  rdata <- cbs_data %>%
    filter(population == race[i]) %>%
    arrange(date)
  
  all_observed[[i]] <- rdata$seroprev_est
  all_dates[[i]] <- rdata$date
  all_sample_sizes[[i]] <- rdata$sero_denom
}

# 4. DEFINE SI MODEL WITH SEROREVERSION
# SI model simulation from March 2020
simulate_si_saturated <- function(lambda_pre, lambda_omi, rho, 
                                  target_dates, model_start_date, omicron_date) {
  
  # Create complete date sequence from model start to last observation
  all_dates_seq <- seq(model_start_date, max(target_dates), by = "7 days")
  num_steps <- length(all_dates_seq)
  
  susceptible <- numeric(num_steps)
  infected <- numeric(num_steps)
  
  # Initialize at March 2020
  susceptible[1] <- 1
  infected[1] <- 0
  
  # Run model forward from March 2020
  for (t in 2:num_steps) {
    # Switch force of infection at Omicron emergence
    if (all_dates_seq[t] < omicron_date) {
      r1 <- lambda_pre
    } else {
      r1 <- lambda_omi
    }
    r2 <- rho
    
    # SI dynamics with seroreversion
    new_susceptible <- susceptible[t - 1] * exp(-r1) + infected[t - 1] * (1 - exp(-r2))
    new_infected <- infected[t - 1] * exp(-r2) + susceptible[t - 1] * (1 - exp(-r1))
    
    susceptible[t] <- new_susceptible
    infected[t] <- new_infected
  }
  
  # Extract predictions for target dates only
  predictions <- numeric(length(target_dates))
  for (i in 1:length(target_dates)) {
    idx <- which(all_dates_seq == target_dates[i])
    if (length(idx) > 0) {
      predictions[i] <- infected[idx]
    } else {
      # Interpolate if exact date not in sequence
      closest_idx <- which.min(abs(all_dates_seq - target_dates[i]))
      predictions[i] <- infected[closest_idx]
    }
  }
  
  return(predictions)
}

# 5. DEFINE LIKELIHOOD FUNCTION
# Joint log-likelihood across all racial groups
joint_loglik_saturated <- function(lambda_white_pre, lambda_racialized_pre,
                                   lambda_white_omi, lambda_racialized_omi, rho) {
  
  lambdas_pre <- c(lambda_white_pre, lambda_racialized_pre)
  lambdas_omi <- c(lambda_white_omi, lambda_racialized_omi)
  
  total_ll <- 0
  
  for (i in 1:2) {
    # Simulate model for this racial group
    pred <- simulate_si_saturated(
      lambdas_pre[i], lambdas_omi[i], rho,
      all_dates[[i]], model_start_date, omicron_date
    )
    
    obs_prop <- all_observed[[i]] / 100
    n <- all_sample_sizes[[i]]
    
    # Binomial likelihood with actual sample sizes
    obs_positive <- round(n * obs_prop)
    ll <- sum(dbinom(x = obs_positive, 
                     size = n, 
                     prob = pred, 
                     log = TRUE))
    total_ll <- total_ll + ll
  }
  
  return(total_ll)
}

# Objective function (negative log-likelihood)
f_saturated <- function(logit_lambda_white_pre, logit_lambda_racialized_pre,
                        logit_lambda_white_omi, logit_lambda_racialized_omi,
                        logit_rho) {
  
  # Transform parameters from logit scale
  lambda_white_pre <- plogis(logit_lambda_white_pre)
  lambda_racialized_pre <- plogis(logit_lambda_racialized_pre)
  
  lambda_white_omi <- plogis(logit_lambda_white_omi)
  lambda_racialized_omi <- plogis(logit_lambda_racialized_omi)
  
  # Transform rho to [0, 0.01] range (0-1% seroreversion per week)
  rho <- plogis(logit_rho) * 0.01
  
  # Return negative log-likelihood
  -joint_loglik_saturated(lambda_white_pre, lambda_racialized_pre,
                          lambda_white_omi, lambda_racialized_omi, rho)
}

# 6. SET INITIAL PARAMETER GUESSES
guess_sat <- list(
  logit_lambda_white_pre = qlogis(0.0005),
  logit_lambda_racialized_pre = qlogis(0.0008),
  logit_lambda_white_omi = qlogis(0.025),
  logit_lambda_racialized_omi = qlogis(0.028),
  logit_rho = qlogis(0.1)  # Starts at rho = 0.001
)

# 7. FIT MODEL USING MAXIMUM LIKELIHOOD ESTIMATION
cat("=== FITTING SATURATED MODEL ===\n")
cat("Parameters: 2 λ_pre + 2 λ_omi + 1 ρ = 5 parameters\n")
cat("Starting maximum likelihood estimation...\n\n")

fit_saturated <- mle2(
  function(logit_lambda_white_pre, logit_lambda_racialized_pre,
           logit_lambda_white_omi, logit_lambda_racialized_omi, logit_rho)
    f_saturated(logit_lambda_white_pre, logit_lambda_racialized_pre,
                logit_lambda_white_omi, logit_lambda_racialized_omi, logit_rho),
  start = guess_sat,
  method = "L-BFGS-B"
)

cat("✓ Fitting complete!\n\n")

# 8. EXTRACT AND DISPLAY RESULTS
coefs <- coef(fit_saturated)

# Transform back to natural scale
lambda_white_pre <- plogis(coefs["logit_lambda_white_pre"])
lambda_racialized_pre <- plogis(coefs["logit_lambda_racialized_pre"])

lambda_white_omi <- plogis(coefs["logit_lambda_white_omi"])
lambda_racialized_omi <- plogis(coefs["logit_lambda_racialized_omi"])

rho <- plogis(coefs["logit_rho"]) * 0.01

cat("=== ESTIMATED SEROREVERSION RATE ===\n\n")
cat(sprintf("ρ = %.5f per week (%.3f%% per week)\n", rho, rho * 100))
cat(sprintf("Half-life = %.1f weeks (%.1f months)\n\n", log(2)/rho, log(2)/rho/4.33))

cat("=== PRE-OMICRON FORCE OF INFECTION ===\n\n")
cat(sprintf("White: %.5f per week\n", lambda_white_pre))
cat(sprintf("Racialized: %.5f per week\n", lambda_racialized_pre))

cat("\n=== OMICRON FORCE OF INFECTION ===\n\n")
cat(sprintf("White: %.5f per week\n", lambda_white_omi))
cat(sprintf("Racialized: %.5f per week\n", lambda_racialized_omi))

cat("\n=== OMICRON MULTIPLIERS ===\n\n")
cat(sprintf("White: %.2fx\n", lambda_white_omi / lambda_white_pre))
cat(sprintf("Racialized: %.2fx\n", lambda_racialized_omi / lambda_racialized_pre))

cat("\n=== PRE-OMICRON IRRs (vs White) ===\n\n")
cat(sprintf("White: 1.000 (reference)\n"))
cat(sprintf("Racialized: %.3f\n", lambda_racialized_pre / lambda_white_pre))

cat("\n=== OMICRON IRRs (vs White) ===\n\n")
cat(sprintf("White: 1.000 (reference)\n"))
cat(sprintf("Racialized: %.3f\n", lambda_racialized_omi / lambda_white_omi))

results <- list(
  fit = fit_saturated,
  loglik = logLik(fit_saturated),
  k = 5,
  aic = AIC(fit_saturated),
  bic = BIC(fit_saturated),
  rho = rho,
  model_start_date = model_start_date,
  omicron_date = omicron_date
)
cat("Model fit statistics:\n")
cat(sprintf("  Log-likelihood: %.2f\n", results$loglik))
cat(sprintf("  AIC: %.2f\n", results$aic))
cat(sprintf("  Parameters: %d\n\n", results$k))

saveRDS(results, "output/saturated_model_estimate_rho.rds")

cat("\n✓ Results saved: saturated_model_estimate_rho.rds\n\n")
cat("Model fit statistics:\n")
cat(sprintf("  Log-likelihood: %.2f\n", results$loglik))
cat(sprintf("  AIC: %.2f\n", results$aic))
cat(sprintf("  Parameters: %d\n\n", results$k))