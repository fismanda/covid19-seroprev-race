# ==============================================================================
# Helper Functions for COVID-19 Seroprevalence by Race Analysis
# ==============================================================================
# Shared functions used across multiple analysis scripts.
# Source this file before running scripts that compute CIs.
# ==============================================================================

# Ensure bbmle is available (needed for mle2 objects)
if (!requireNamespace("bbmle", quietly = TRUE)) {
  stop("Package 'bbmle' is required. Install with: install.packages('bbmle')")
}
library(bbmle)

#' Calculate confidence interval for a ratio of transformed parameters
#' using the delta method on the log scale.
#'
#' @param fit A fitted mle2 model object
#' @param param_numerator Name of the logit-scale numerator parameter
#' @param param_denominator Name of the logit-scale denominator parameter
#' @param alpha Significance level (default 0.05 for 95% CI)
#' @return Named numeric vector with 'lower' and 'upper' bounds
calc_irr_ci <- function(fit, param_numerator, param_denominator, alpha = 0.05) {
  
  # Define the ratio function
  ratio_fn <- function(params) {
    num <- plogis(params[param_numerator])
    den <- plogis(params[param_denominator])
    return(num / den)
  }
  
  # Get point estimate
  point_est <- ratio_fn(coef(fit))
  
  # Get variance-covariance matrix
  vcov_mat <- vcov(fit)
  
  # Numerical gradient for delta method
  eps <- 1e-7
  coefs <- coef(fit)
  
  grad <- numeric(length(coefs))
  for (i in 1:length(coefs)) {
    coefs_plus <- coefs
    coefs_plus[i] <- coefs[i] + eps
    
    coefs_minus <- coefs
    coefs_minus[i] <- coefs[i] - eps
    
    grad[i] <- (ratio_fn(coefs_plus) - ratio_fn(coefs_minus)) / (2 * eps)
  }
  
  # Calculate standard error using delta method
  se <- sqrt(t(grad) %*% vcov_mat %*% grad)
  
  # Calculate CI on log scale (better for ratios)
  log_est <- log(point_est)
  log_se <- se / point_est
  
  z_crit <- qnorm(1 - alpha/2)
  
  ci_lower <- exp(log_est - z_crit * log_se)
  ci_upper <- exp(log_est + z_crit * log_se)
  
  return(c(lower = ci_lower, upper = ci_upper))
}

#' Calculate confidence interval for an Omicron multiplier (ratio of
#' Omicron to pre-Omicron FOI) using the delta method on the log scale.
#'
#' @param fit A fitted mle2 model object
#' @param param_pre Name of the logit-scale pre-Omicron parameter
#' @param param_omi Name of the logit-scale Omicron parameter
#' @param alpha Significance level (default 0.05 for 95% CI)
#' @return Named numeric vector with 'lower' and 'upper' bounds
calc_mult_ci <- function(fit, param_pre, param_omi, alpha = 0.05) {
  
  ratio_fn <- function(params) {
    omi <- plogis(params[param_omi])
    pre <- plogis(params[param_pre])
    return(omi / pre)
  }
  
  point_est <- ratio_fn(coef(fit))
  vcov_mat <- vcov(fit)
  
  eps <- 1e-7
  coefs <- coef(fit)
  
  grad <- numeric(length(coefs))
  for (i in 1:length(coefs)) {
    coefs_plus <- coefs
    coefs_plus[i] <- coefs[i] + eps
    
    coefs_minus <- coefs
    coefs_minus[i] <- coefs[i] - eps
    
    grad[i] <- (ratio_fn(coefs_plus) - ratio_fn(coefs_minus)) / (2 * eps)
  }
  
  se <- sqrt(t(grad) %*% vcov_mat %*% grad)
  
  log_est <- log(point_est)
  log_se <- se / point_est
  
  z_crit <- qnorm(1 - alpha/2)
  
  ci_lower <- exp(log_est - z_crit * log_se)
  ci_upper <- exp(log_est + z_crit * log_se)
  
  return(c(lower = ci_lower, upper = ci_upper))
}

#' SI model simulation with seroreversion (saturated specification)
#'
#' @param lambda_pre Pre-Omicron force of infection (per week)
#' @param lambda_omi Omicron force of infection (per week)
#' @param rho Seroreversion rate (per week)
#' @param target_dates Dates at which to extract predictions
#' @param model_start_date Date at which model is initialized (S=1, I=0)
#' @param omicron_date Date of Omicron emergence
#' @return Numeric vector of predicted seroprevalence (proportion) at target_dates
simulate_si_saturated <- function(lambda_pre, lambda_omi, rho, 
                                  target_dates, model_start_date, omicron_date) {
  
  all_dates_seq <- seq(model_start_date, max(target_dates), by = "7 days")
  num_steps <- length(all_dates_seq)
  
  susceptible <- numeric(num_steps)
  infected <- numeric(num_steps)
  
  susceptible[1] <- 1
  infected[1] <- 0
  
  for (t in 2:num_steps) {
    if (all_dates_seq[t] < omicron_date) {
      r1 <- lambda_pre
    } else {
      r1 <- lambda_omi
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
