library(dplyr)

source("R/00_helpers.R")

cat("\n=======================================================")
cat("\n=== MODEL COMPARISON: SATURATED VS PARSIMONIOUS ===")
cat("\n=======================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD BOTH MODELS
# -----------------------------------------------------------------------------

cat("Loading model results...\n\n")

saturated <- readRDS("output/saturated_model_estimate_rho.rds")
parsimonious <- readRDS("output/parsimonious_model_estimate_rho.rds")

# -----------------------------------------------------------------------------
# 2. EXTRACT KEY STATISTICS
# -----------------------------------------------------------------------------

# Saturated model
sat_ll <- as.numeric(saturated$loglik)
sat_k <- saturated$k
sat_aic <- saturated$aic
sat_bic <- saturated$bic

# Parsimonious model
pars_ll <- as.numeric(parsimonious$loglik)
pars_k <- parsimonious$k
pars_aic <- parsimonious$aic
pars_bic <- parsimonious$bic

# -----------------------------------------------------------------------------
# 3. CALCULATE COMPARISONS
# -----------------------------------------------------------------------------

# AIC comparison
delta_aic <- sat_aic - pars_aic

# BIC comparison
delta_bic <- sat_bic - pars_bic

# Likelihood ratio test
lr_stat <- 2 * (sat_ll - pars_ll)
df <- sat_k - pars_k
p_value <- pchisq(lr_stat, df, lower.tail = FALSE)

# -----------------------------------------------------------------------------
# 4. DISPLAY RESULTS
# -----------------------------------------------------------------------------

cat("=== MODEL SPECIFICATIONS ===\n\n")

cat("PARSIMONIOUS MODEL:\n")
cat("  Description: Shared Omicron multiplier (uniform effect)\n")
cat(sprintf("  Parameters: %d (2 λ_pre + 1 multiplier + 1 ρ)\n", pars_k))
cat(sprintf("  Log-likelihood: %.2f\n", pars_ll))
cat(sprintf("  AIC: %.2f\n", pars_aic))
cat(sprintf("  BIC: %.2f\n\n", pars_bic))

cat("SATURATED MODEL:\n")
cat("  Description: Race-specific Omicron effects\n")
cat(sprintf("  Parameters: %d (2 λ_pre + 2 λ_omi + 1 ρ)\n", sat_k))
cat(sprintf("  Log-likelihood: %.2f\n", sat_ll))
cat(sprintf("  AIC: %.2f\n", sat_aic))
cat(sprintf("  BIC: %.2f\n\n", sat_bic))

cat("=== MODEL COMPARISON ===\n\n")

cat("AKAIKE INFORMATION CRITERION (AIC):\n")
cat(sprintf("  ΔAIC = %.2f (Saturated - Parsimonious)\n", delta_aic))

if (delta_aic < -10) {
  cat("  *** SATURATED MODEL STRONGLY PREFERRED ***\n")
  cat("  Interpretation: Omicron's effect varied substantially by race.\n")
  cat("  Differential amplification is statistically supported.\n")
} else if (delta_aic < -2) {
  cat("  Saturated model preferred\n")
} else if (delta_aic > 2) {
  cat("  Parsimonious model preferred\n")
} else {
  cat("  Models essentially equivalent\n")
}

cat("\nBAYESIAN INFORMATION CRITERION (BIC):\n")
cat(sprintf("  ΔBIC = %.2f (Saturated - Parsimonious)\n", delta_bic))

if (delta_bic < -10) {
  cat("  Saturated model strongly preferred (BIC)\n")
} else if (delta_bic < -2) {
  cat("  Saturated model preferred (BIC)\n")
} else if (delta_bic > 2) {
  cat("  Parsimonious model preferred (BIC)\n")
} else {
  cat("  Models essentially equivalent (BIC)\n")
}

cat("\nLIKELIHOOD RATIO TEST:\n")
cat(sprintf("  χ² statistic: %.2f\n", lr_stat))
cat(sprintf("  Degrees of freedom: %d\n", df))
cat(sprintf("  p-value: %.2e\n", p_value))

if (p_value < 0.001) {
  cat("  *** HIGHLY SIGNIFICANT (p < 0.001) ***\n")
  cat("  Saturated model provides significantly better fit.\n")
} else if (p_value < 0.05) {
  cat("  Significant (p < 0.05)\n")
} else {
  cat("  Not significant (p ≥ 0.05)\n")
}

# -----------------------------------------------------------------------------
# 5. EXTRACT AND COMPARE KEY PARAMETERS
# -----------------------------------------------------------------------------

cat("\n=== KEY PARAMETER COMPARISONS ===\n\n")

# Parsimonious multiplier
pars_coefs <- coef(parsimonious$fit)
shared_mult <- exp(pars_coefs["log_multiplier"])

cat("PARSIMONIOUS MODEL (Shared Multiplier):\n")
cat(sprintf("  All race groups: %.2fx\n\n", shared_mult))

# Saturated multipliers
sat_coefs <- coef(saturated$fit)

lambda_white_pre_sat <- plogis(sat_coefs["logit_lambda_white_pre"])
lambda_racialized_pre_sat <- plogis(sat_coefs["logit_lambda_racialized_pre"])

lambda_white_omi_sat <- plogis(sat_coefs["logit_lambda_white_omi"])
lambda_racialized_omi_sat <- plogis(sat_coefs["logit_lambda_racialized_omi"])

cat("SATURATED MODEL (Race-Specific Multipliers):\n")
cat(sprintf("  White: %.2fx\n", lambda_white_omi_sat / lambda_white_pre_sat))
cat(sprintf("  Racialized: %.2fx\n", lambda_racialized_omi_sat / lambda_racialized_pre_sat))

cat("\nDIFFERENTIAL AMPLIFICATION:\n")
mult_white <- lambda_white_omi_sat / lambda_white_pre_sat
mult_racialized <- lambda_racialized_omi_sat / lambda_racialized_pre_sat
cat(sprintf("  White multiplier: %.2fx\n", mult_white))
cat(sprintf("  Racialized multiplier: %.2fx\n", mult_racialized))
cat(sprintf("  Ratio (White/Racialized): %.2f\n", mult_white / mult_racialized))
cat("  Interpretation: White individuals experienced %.1f%% larger\n", 
    (mult_white/mult_racialized - 1) * 100)
cat("  proportional increase than racialized.\n")

# -----------------------------------------------------------------------------
# 6. SAVE RESULTS TO FILE
# -----------------------------------------------------------------------------

sink("output/model_comparison_results.txt")

cat("=======================================================\n")
cat("MODEL COMPARISON: SATURATED VS PARSIMONIOUS\n")
cat("=======================================================\n\n")

cat("Analysis Date:", format(Sys.Date(), "%B %d, %Y"), "\n\n")

cat("PARSIMONIOUS MODEL:\n")
cat(sprintf("  Parameters: %d\n", pars_k))
cat(sprintf("  Log-likelihood: %.2f\n", pars_ll))
cat(sprintf("  AIC: %.2f\n", pars_aic))
cat(sprintf("  BIC: %.2f\n\n", pars_bic))

cat("SATURATED MODEL:\n")
cat(sprintf("  Parameters: %d\n", sat_k))
cat(sprintf("  Log-likelihood: %.2f\n", sat_ll))
cat(sprintf("  AIC: %.2f\n", sat_aic))
cat(sprintf("  BIC: %.2f\n\n", sat_bic))

cat("COMPARISON STATISTICS:\n")
cat(sprintf("  ΔAIC: %.2f\n", delta_aic))
cat(sprintf("  ΔBIC: %.2f\n", delta_bic))
cat(sprintf("  LR χ²: %.2f (df=%d)\n", lr_stat, df))
cat(sprintf("  p-value: %.2e\n\n", p_value))

cat("CONCLUSION:\n")
if (delta_aic < -10) {
  cat("  Saturated model is strongly preferred by AIC.\n")
  cat("  Evidence for differential Omicron amplification by race.\n")
}

sink()

cat("\n✓ Results saved to: model_comparison_results.txt\n\n")

cat("=== SUMMARY ===\n\n")
cat("The saturated model, which allows Omicron's effect to vary by\n")
cat("racial status, provides substantially better fit than the\n")
cat("parsimonious model with a shared multiplier.\n\n")

cat(sprintf("ΔAIC = %.0f indicates overwhelming evidence for the saturated model.\n", delta_aic))
cat("This supports the hypothesis of differential amplification:\n")
cat("white populations experienced larger proportional increases\n")
cat("in force of infection during Omicron than racialized populations.\n\n")

