library(bbmle)
library(dplyr)

source("R/00_helpers.R")

cat("\n==============================================================")
cat("\n=== CALCULATING OMICRON MULTIPLIERS WITH CONFIDENCE INTERVALS ===")
cat("\n==============================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD SATURATED MODEL
# -----------------------------------------------------------------------------

cat("Loading saturated model...\n")
results <- readRDS("output/saturated_model_estimate_rho.rds")
fit <- results$fit

coefs <- coef(fit)

# Pre-Omicron
lambda_white_pre <- plogis(coefs["logit_lambda_white_pre"])
lambda_racialized_pre <- plogis(coefs["logit_lambda_racialized_pre"])

# Omicron
lambda_white_omi <- plogis(coefs["logit_lambda_white_omi"])
lambda_racialized_omi <- plogis(coefs["logit_lambda_racialized_omi"])

# Calculate multipliers (point estimates)
mult <- c(
  white = lambda_white_omi / lambda_white_pre,
  racialized = lambda_racialized_omi / lambda_racialized_pre
)

cat("✓ Model loaded\n\n")

# -----------------------------------------------------------------------------
# 2. CALCULATE CONFIDENCE INTERVALS
# -----------------------------------------------------------------------------

cat("Calculating 95% confidence intervals...\n")
cat("(calc_mult_ci function loaded from R/00_helpers.R)\n\n")

# Calculate CIs for each racial group
mult_ci <- data.frame(
  Race = c("white", "racialized"),
  Multiplier = mult,
  Lower = NA,
  Upper = NA
)

param_pre <- c("logit_lambda_white_pre", "logit_lambda_racialized_pre")
param_omi <- c("logit_lambda_white_omi", "logit_lambda_racialized_omi")

for (i in 1:2) {
  ci <- calc_mult_ci(fit, param_pre[i], param_omi[i])
  mult_ci$Lower[i] <- ci["lower"]
  mult_ci$Upper[i] <- ci["upper"]
  
  cat(sprintf("%s: %.2fx (95%% CI: %.2f - %.2f)\n", 
              mult_ci$Race[i], mult[i], ci["lower"], ci["upper"]))
}

cat("\n✓ Confidence intervals calculated\n\n")

# -----------------------------------------------------------------------------
# 3. SAVE RESULTS
# -----------------------------------------------------------------------------

write.csv(mult_ci, "output/omicron_multipliers_with_ci.csv", row.names = FALSE)
saveRDS(mult_ci, "output/omicron_multipliers_with_ci.rds")

cat("✓ Results saved:\n")
cat("  - omicron_multipliers_with_ci.csv\n")
cat("  - omicron_multipliers_with_ci.rds\n\n")

# -----------------------------------------------------------------------------
# 4. FORMATTED OUTPUT
# -----------------------------------------------------------------------------

cat("=== OMICRON MULTIPLIERS WITH 95% CONFIDENCE INTERVALS ===\n\n")
cat("Fold-change in force of infection from pre-Omicron to Omicron:\n")
cat("─────────────────────────────────────────────────────\n")
cat("Race    Multiplier    95% CI\n")
cat("─────────────────────────────────────────────────────\n")

for (i in 1:2) {
  r <- mult_ci$Race[i]
  m <- mult_ci$Multiplier[i]
  ci_str <- sprintf("%.2f - %.2f", mult_ci$Lower[i], mult_ci$Upper[i])
  cat(sprintf("%-10s  %.2fx         (%s)\n", r, m, ci_str))
}

cat("\n")

# -----------------------------------------------------------------------------
# 5. DIFFERENTIAL AMPLIFICATION ANALYSIS
# -----------------------------------------------------------------------------

cat("=== DIFFERENTIAL AMPLIFICATION ===\n\n")

mult_white <- mult[1]
mult_racialized <- mult[2]
diff_ratio <- mult_white / mult_racialized

cat(sprintf("white:  %.2fx increase\n", mult_white))
cat(sprintf("racialized:   %.2fx increase\n", mult_racialized))
cat(sprintf("Ratio (white/racialized):        %.2f\n\n", diff_ratio))

cat("INTERPRETATION:\n")
cat(sprintf("White individuals experienced a %.0f%% larger\n", 
            (diff_ratio - 1) * 100))
cat("proportional increase in force of infection compared to\n")
cat("racialized individuals during Omicron.\n\n")

cat("This 'differential amplification' resulted in gradient compression:\n")
cat("white populations 'caught up' to levels already experienced\n")
cat("by disadvantaged populations, producing apparent convergence\n")
cat("in seroprevalence despite persistent transmission disparities.\n\n")

# -----------------------------------------------------------------------------
# 6. COMPARE TO PARSIMONIOUS MODEL
# -----------------------------------------------------------------------------

cat("=== COMPARISON: SATURATED vs PARSIMONIOUS ===\n\n")

cat("Loading parsimonious model...\n")
pars_results <- readRDS("output/parsimonious_model_estimate_rho.rds")
pars_coefs <- coef(pars_results$fit)
shared_mult <- exp(pars_coefs["log_multiplier"])
pars_rho <- pars_results$rho

cat("✓ Loaded\n\n")

cat("PARSIMONIOUS MODEL (Shared Multiplier):\n")
cat(sprintf("  All races: %.2fx\n", shared_mult))
cat(sprintf("  Seroreversion: ρ = %.5f/week (half-life: %.1f months)\n\n", 
            pars_rho, log(2)/pars_rho/4.33))

cat("SATURATED MODEL (Race-Specific Multipliers):\n")
for (i in 1:2) {
  cat(sprintf("  %s: %.2fx\n", mult_ci$Race[i], mult[i]))
}
sat_rho <- results$rho
cat(sprintf("  Seroreversion: ρ = %.5f/week (half-life: %.1f months)\n\n", 
            sat_rho, log(2)/sat_rho/4.33))

cat("SEROREVERSION COMPARISON:\n")
cat(sprintf("  Saturated:     ρ = %.5f/week\n", sat_rho))
cat(sprintf("  Parsimonious:  ρ = %.5f/week\n", pars_rho))
if (abs(sat_rho - pars_rho) < 0.0005) {
  cat("  → Nearly identical (minimal antibody waning)\n\n")
} else {
  cat(sprintf("  → Difference: %.5f/week\n\n", abs(sat_rho - pars_rho)))
}

cat("\nThe parsimonious model assumes uniform Omicron amplification,\n")
cat(sprintf("averaging to %.2fx across all groups. The saturated model\n", shared_mult))
cat(sprintf("reveals substantial heterogeneity (%.2fx in white vs %.2fx in racialized),\n", mult[1], mult[2]))
cat("with ΔAIC = -192.5 strongly favoring the saturated specification.\n\n")

# -----------------------------------------------------------------------------
# 7. ABSOLUTE CHANGES IN FORCE OF INFECTION
# -----------------------------------------------------------------------------

cat("=== ABSOLUTE CHANGES IN FORCE OF INFECTION ===\n\n")
cat("                Pre-Omicron    Omicron      Absolute\n")
cat("Race            (per week)     (per week)   Change\n")
cat("─────────────────────────────────────────────────────\n")

lambdas_pre <- c(lambda_white_pre, lambda_racialized_pre)
lambdas_omi <- c(lambda_white_omi, lambda_racialized_omi)

for (i in 1:2) {
  abs_change <- lambdas_omi[i] - lambdas_pre[i]
  cat(sprintf("%s          %.5f        %.5f      +%.5f\n", 
              mult_ci$Race[i], lambdas_pre[i], lambdas_omi[i], abs_change))
}

cat("\nNote: Despite larger proportional increases for white individuals, absolute\n")
cat("increases were similar across racial groups (~0.02/week), resulting\n")
cat("in convergence to similar final force of infection levels.\n\n")
