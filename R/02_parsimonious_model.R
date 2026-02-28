# Set working directory to project root (adjust as needed)
# setwd("path/to/covid19-seroprev-race")

library(readr)
library(bbmle)
library(dplyr)
library(lubridate)

source("R/00_helpers.R")

# FITTING PARSIMONIOUS MODEL: SHARED OMICRON MULTIPLIER ------------------- 
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
# (anti-N and anti-S values appear transposed; surrounding values are ~5%)
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

# 4. DEFINE SI MODEL (with shared multiplier)
simulate_si_parsimonious <- function(lambda_pre, multiplier, rho,
                                     target_dates, model_start_date, omicron_date) {
  
  all_dates_seq <- seq(model_start_date, max(target_dates), by = "7 days")
  num_steps <- length(all_dates_seq)
  
  susceptible <- numeric(num_steps)
  infected <- numeric(num_steps)
  
  susceptible[1] <- 1
  infected[1] <- 0
  
  for (t in 2:num_steps) {
    # Apply shared multiplier after Omicron
    if (all_dates_seq[t] < omicron_date) {
      r1 <- lambda_pre
    } else {
      r1 <- lambda_pre * multiplier  # Shared multiplier
    }
    r2 <- rho
    
    new_susceptible <- susceptible[t - 1] * exp(-r1) + infected[t - 1] * (1 - exp(-r2))
    new_infected <- infected[t - 1] * exp(-r2) + susceptible[t - 1] * (1 - exp(-r1))
    
    susceptible[t] <- new_susceptible
    infected[t] <- new_infected
  }
  
  predictions <- numeric(length(target_dates))
  for (i in 1:length(target_dates)) {
    idx <- which(all_dates_seq == target_dates[i])
    if (length(idx) > 0) {
      predictions[i] <- infected[idx]
    } else {
      closest_idx <- which.min(abs(all_dates_seq - target_dates[i]))
      predictions[i] <- infected[closest_idx]
    }
  }
  
  return(predictions)
}

# 5. DEFINE LIKELIHOOD FUNCTION
joint_loglik_parsimonious <- function(lambda_white_pre, lambda_racialized_pre,
                                      multiplier, rho) {
  
  lambdas_pre <- c(lambda_white_pre, lambda_racialized_pre)
  
  total_ll <- 0
  
  for (i in 1:2) {
    pred <- simulate_si_parsimonious(
      lambdas_pre[i], multiplier, rho,
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

# Objective function
f_parsimonious <- function(logit_lambda_white_pre, logit_lambda_racialized_pre,
                           log_multiplier, logit_rho) {
  
  lambda_white_pre <- plogis(logit_lambda_white_pre)
  lambda_racialized_pre <- plogis(logit_lambda_racialized_pre)
  
  multiplier <- exp(log_multiplier)  # Ensure positive
  rho <- plogis(logit_rho) * 0.01
  
  -joint_loglik_parsimonious(lambda_white_pre, lambda_racialized_pre,
                             multiplier, rho)
}

# 6. SET INITIAL GUESSES
guess_pars <- list(
  logit_lambda_white_pre = qlogis(0.0005),
  logit_lambda_racialized_pre = qlogis(0.0008),
  log_multiplier = log(35),  # Starts at ~35x multiplier
  logit_rho = qlogis(0.1)
)

# 7. FIT PARSIMONIOUS MODEL
cat("=== FITTING PARSIMONIOUS MODEL ===\n")
cat("Parameters: 2 λ_pre + 1 shared multiplier + 1 ρ = 4 parameters\n")
cat("Starting maximum likelihood estimation...\n\n")

fit_parsimonious <- mle2(
  function(logit_lambda_white_pre, logit_lambda_racialized_pre,
           log_multiplier, logit_rho)
    f_parsimonious(logit_lambda_white_pre, logit_lambda_racialized_pre, 
                   log_multiplier, logit_rho),
  start = guess_pars,
  method = "L-BFGS-B"
)

cat("✓ Fitting complete!\n\n")

# 8. EXTRACT AND DISPLAY RESULTS
coefs <- coef(fit_parsimonious)

lambda_white_pre <- plogis(coefs["logit_lambda_white_pre"])
lambda_racialized_pre <- plogis(coefs["logit_lambda_racialized_pre"])

multiplier <- exp(coefs["log_multiplier"])
rho <- plogis(coefs["logit_rho"]) * 0.01

cat("=== ESTIMATED SEROREVERSION RATE ===\n\n")
cat(sprintf("ρ = %.5f per week (%.3f%% per week)\n", rho, rho * 100))
cat(sprintf("Half-life = %.1f weeks (%.1f months)\n\n", log(2)/rho, log(2)/rho/4.33))

cat("=== SHARED OMICRON MULTIPLIER ===\n\n")
cat(sprintf("Multiplier = %.2fx (applied to all racial groups)\n\n", multiplier))

cat("=== PRE-OMICRON FORCE OF INFECTION ===\n\n")
cat(sprintf("White: %.5f per week\n", lambda_white_pre))
cat(sprintf("Racialized: %.5f per week\n", lambda_racialized_pre))

cat("\n=== OMICRON FORCE OF INFECTION (λ_pre × multiplier) ===\n\n")
cat(sprintf("White: %.5f per week\n", lambda_white_pre * multiplier))
cat(sprintf("Racialized: %.5f per week\n", lambda_racialized_pre * multiplier))

cat("\n=== PRE-OMICRON IRRs (vs White) ===\n\n")
cat(sprintf("White: 1.000 (reference)\n"))
cat(sprintf("Racialized: %.3f\n", lambda_racialized_pre / lambda_white_pre))

cat("\n=== OMICRON IRRs (vs White) - SAME AS PRE-OMICRON ===\n")
cat("(Because shared multiplier preserves gradient)\n\n")
cat(sprintf("White: 1.000 (reference)\n"))
cat(sprintf("Racialized: %.3f\n", lambda_racialized_pre / lambda_white_pre))

results <- list(
  fit = fit_parsimonious,
  loglik = logLik(fit_parsimonious),
  k = 4,
  aic = AIC(fit_parsimonious),
  bic = BIC(fit_parsimonious),
  rho = rho,
  multiplier = multiplier,
  model_start_date = model_start_date,
  omicron_date = omicron_date
)

cat("Model fit statistics:\n")
cat(sprintf("  Log-likelihood: %.2f\n", results$loglik))
cat(sprintf("  AIC: %.2f\n", results$aic))
cat(sprintf("  Parameters: %d\n\n", results$k))

saveRDS(results, "output/parsimonious_model_estimate_rho.rds")

cat("\n✓ Results saved: parsimonious_model_estimate_rho.rds\n\n")
cat("Model fit statistics:\n")
cat(sprintf("  Log-likelihood: %.2f\n", results$loglik))
cat(sprintf("  AIC: %.2f\n", results$aic))
cat(sprintf("  Parameters: %d\n\n", results$k))
