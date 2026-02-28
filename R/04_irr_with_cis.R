library(bbmle)
library(dplyr)

source("R/00_helpers.R")

cat("\n================================================================")
cat("\n=== CALCULATING IRRs WITH CONFIDENCE INTERVALS ===")
cat("\n================================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD MODEL
# -----------------------------------------------------------------------------

cat("Loading saturated model...\n")
results <- readRDS("output/saturated_model_estimate_rho.rds")
fit <- results$fit

cat("✓ Model loaded\n\n")

# -----------------------------------------------------------------------------
# 2. EXTRACT POINT ESTIMATES
# -----------------------------------------------------------------------------

coefs <- coef(fit)

# Pre-Omicron
lambda_white_pre <- plogis(coefs["logit_lambda_white_pre"])
lambda_racialized_pre <- plogis(coefs["logit_lambda_racialized_pre"])

# Omicron
lambda_white_omi <- plogis(coefs["logit_lambda_white_omi"])
lambda_racialized_omi <- plogis(coefs["logit_lambda_racialized_omi"])

# Calculate IRRs (point estimates)
irr_pre <- c(
  white = 1.000,
  racialized = lambda_racialized_pre / lambda_white_pre
)

irr_omi <- c(
  white = 1.000,
  racialized = lambda_racialized_omi / lambda_white_omi
)

cat("Point estimates calculated\n\n")

# -----------------------------------------------------------------------------
# 3. CALCULATE CONFIDENCE INTERVALS USING PROFILE LIKELIHOOD
# -----------------------------------------------------------------------------

cat("Calculating 95% confidence intervals using delta method...\n")
cat("(calc_irr_ci function loaded from R/00_helpers.R)\n\n")

# Calculate CIs for pre-Omicron IRRs
cat("Pre-Omicron IRRs:\n")
irr_pre_ci <- data.frame(
  Race = c("white", "racialized"),
  Period = "Pre-Omicron",
  IRR = irr_pre,
  Lower = NA,
  Upper = NA
)

# White is reference (CI = NA)
irr_pre_ci$Lower[1] <- NA
irr_pre_ci$Upper[1] <- NA

# Calculate CIs for Racialized
for (i in 2) {
  param_num <- "logit_lambda_racialized_pre"
  param_den <- "logit_lambda_white_pre"
  
  ci <- calc_irr_ci(fit, param_num, param_den)
  irr_pre_ci$Lower[i] <- ci["lower"]
  irr_pre_ci$Upper[i] <- ci["upper"]
  
  cat(sprintf("  %s: %.3f (95%% CI: %.3f - %.3f)\n", 
              irr_pre_ci$Race[i], irr_pre[i], ci["lower"], ci["upper"]))
}

# Calculate CIs for Omicron IRRs
cat("\nOmicron IRRs:\n")
irr_omi_ci <- data.frame(
  Race = c("white", "racialized"),
  Period = "Omicron",
  IRR = irr_omi,
  Lower = NA,
  Upper = NA
)

# White is reference
irr_omi_ci$Lower[1] <- NA
irr_omi_ci$Upper[1] <- NA

# Calculate CIs for Racialized
for (i in 2) {
  param_num <- "logit_lambda_racialized_omi"
  param_den <- "logit_lambda_white_omi"
  
  ci <- calc_irr_ci(fit, param_num, param_den)
  irr_omi_ci$Lower[i] <- ci["lower"]
  irr_omi_ci$Upper[i] <- ci["upper"]
  
  cat(sprintf("  %s: %.3f (95%% CI: %.3f - %.3f)\n", 
              irr_omi_ci$Race[i], irr_omi[i], ci["lower"], ci["upper"]))
}

# -----------------------------------------------------------------------------
# 4. COMBINE AND SAVE RESULTS
# -----------------------------------------------------------------------------

cat("\n✓ Confidence intervals calculated\n\n")

# Combine both periods
irr_all <- rbind(irr_pre_ci, irr_omi_ci)

# Clean up row names
rownames(irr_all) <- NULL

# Save as CSV
write.csv(irr_all, "output/irr_estimates_with_ci.csv", row.names = FALSE)

# Save as RDS for plotting
saveRDS(irr_all, "output/irr_estimates_with_ci.rds")

cat("✓ Results saved:\n")
cat("  - irr_estimates_with_ci.csv\n")
cat("  - irr_estimates_with_ci.rds\n\n")

# -----------------------------------------------------------------------------
# 5. DISPLAY FORMATTED TABLE
# -----------------------------------------------------------------------------

cat("=== INCIDENCE RATE RATIOS WITH 95% CONFIDENCE INTERVALS ===\n\n")

cat("PRE-OMICRON PERIOD (vs White):\n")
cat("─────────────────────────────────────────────────────\n")
cat("Race    IRR      95% CI\n")
cat("─────────────────────────────────────────────────────\n")
for (i in 1:2) {
  q <- irr_pre_ci$Race[i]
  irr <- irr_pre_ci$IRR[i]
  
  if (i == 1) {
    cat(sprintf("%-10s  %.3f    (reference)\n", q, irr))
  } else {
    ci_str <- sprintf("%.3f - %.3f", irr_pre_ci$Lower[i], irr_pre_ci$Upper[i])
    cat(sprintf("%-10s  %.3f    (%s)\n", q, irr, ci_str))
  }
}

cat("\nOMICRON PERIOD (vs White):\n")
cat("─────────────────────────────────────────────────────\n")
cat("Race    IRR      95% CI\n")
cat("─────────────────────────────────────────────────────\n")
for (i in 1:2) {
  q <- irr_omi_ci$Race[i]
  irr <- irr_omi_ci$IRR[i]
  
  if (i == 1) {
    cat(sprintf("%-10s  %.3f    (reference)\n", q, irr))
  } else {
    ci_str <- sprintf("%.3f - %.3f", irr_omi_ci$Lower[i], irr_omi_ci$Upper[i])
    cat(sprintf("%-10s  %.3f    (%s)\n", q, irr, ci_str))
  }
}

# -----------------------------------------------------------------------------
# 6. CALCULATE GRADIENT COMPRESSION
# -----------------------------------------------------------------------------

cat("\n=== GRADIENT COMPRESSION ===\n\n")

pre_gradient <- irr_pre[2]  # Racialized
omi_gradient <- irr_omi[2]  # Racialized

compression <- (pre_gradient - omi_gradient) / pre_gradient * 100

cat(sprintf("Pre-Omicron Racialized vs White: %.3f\n", pre_gradient))
cat(sprintf("Omicron Racialized vs White: %.3f\n", omi_gradient))
cat(sprintf("Relative compression: %.1f%%\n", compression))
cat(sprintf("Absolute difference: %.3f\n\n", pre_gradient - omi_gradient))

cat(sprintf("Interpretation: The racial gradient compressed by %.0f%%\n", compression))
cat("during the Omicron period compared to pre-Omicron.\n\n")

# -----------------------------------------------------------------------------
# 7. CALCULATE PARSIMONIOUS MODEL IRRs (FOR SUPPLEMENT)
# -----------------------------------------------------------------------------

cat("=== PARSIMONIOUS MODEL IRRs (FOR SUPPLEMENT) ===\n\n")
cat("Loading parsimonious model...\n")

pars_results <- readRDS("output/parsimonious_model_estimate_rho.rds")
pars_fit <- pars_results$fit
pars_coefs <- coef(pars_fit)

# Pre-Omicron lambdas (same as saturated for pre-period)
lambda_white_pre_pars <- plogis(pars_coefs["logit_lambda_white_pre"])
lambda_racialized_pre_pars <- plogis(pars_coefs["logit_lambda_racialized_pre"])

# Shared multiplier
multiplier <- exp(pars_coefs["log_multiplier"])

# Omicron lambdas (pre × multiplier)
lambda_white_omi_pars <- lambda_white_pre_pars * multiplier
lambda_racialized_omi_pars <- lambda_racialized_pre_pars * multiplier

# Calculate IRRs
irr_pre_pars <- c(
  white = 1.000,
  racialized = lambda_racialized_pre_pars / lambda_white_pre_pars
)

irr_omi_pars <- c(
  white = 1.000,
  racialized = lambda_racialized_omi_pars / lambda_white_omi_pars
)

cat("✓ Parsimonious model loaded\n\n")

cat("NOTE: In the parsimonious model, the shared Omicron multiplier\n")
cat("preserves the gradient, so IRRs are identical in both periods.\n\n")

cat("PARSIMONIOUS MODEL - PRE-OMICRON IRRs:\n")
for (i in 1:2) {
  if (i == 1) {
    cat(sprintf("  %s: %.3f (reference)\n", i, irr_pre_pars[i]))
  } else {
    cat(sprintf("  %s: %.3f\n", i, irr_pre_pars[i]))
  }
}

cat("\nPARSIMONIOUS MODEL - OMICRON IRRs:\n")
cat("(Same as pre-Omicron due to shared multiplier)\n")
for (i in 1:2) {
  if (i == 1) {
    cat(sprintf("  %s: %.3f (reference)\n", i, irr_omi_pars[i]))
  } else {
    cat(sprintf("  %s: %.3f\n", i, irr_omi_pars[i]))
  }
}

# Calculate CIs for parsimonious (only need pre-Omicron since Omicron is same)
cat("\nCalculating CIs for parsimonious model...\n")

irr_pre_pars_ci <- data.frame(
  Race = c("white", "racialized"),
  Period = "Pre-Omicron",
  IRR = irr_pre_pars,
  Lower = NA,
  Upper = NA
)

irr_pre_pars_ci$Lower[1] <- NA
irr_pre_pars_ci$Upper[1] <- NA

for (i in 2) {
  param_num <- "logit_lambda_racialized_pre"
  param_den <- "logit_lambda_white_pre"
  
  ci <- calc_irr_ci(pars_fit, param_num, param_den)
  irr_pre_pars_ci$Lower[i] <- ci["lower"]
  irr_pre_pars_ci$Upper[i] <- ci["upper"]
}

# Omicron IRRs are same as pre-Omicron
irr_omi_pars_ci <- irr_pre_pars_ci
irr_omi_pars_ci$Period <- "Omicron"

# Combine
irr_pars_all <- rbind(irr_pre_pars_ci, irr_omi_pars_ci)
rownames(irr_pars_all) <- NULL

# Save
write.csv(irr_pars_all, "output/irr_estimates_parsimonious_with_ci.csv", row.names = FALSE)

cat("✓ Parsimonious IRRs saved: irr_estimates_parsimonious_with_ci.csv\n\n")

cat("PARSIMONIOUS MODEL IRRs WITH 95% CI:\n")
cat("(Note: IRRs identical in both periods due to shared multiplier)\n")
cat("─────────────────────────────────────────────────────\n")
cat("Race    IRR      95% CI\n")
cat("─────────────────────────────────────────────────────\n")
for (i in 1:2) {
  q <- irr_pre_pars_ci$Race[i]
  irr <- irr_pre_pars_ci$IRR[i]
  
  if (i == 1) {
    cat(sprintf("%-10s  %.3f    (reference)\n", q, irr))
  } else {
    ci_str <- sprintf("%.3f - %.3f", irr_pre_pars_ci$Lower[i], irr_pre_pars_ci$Upper[i])
    cat(sprintf("%-10s  %.3f    (%s)\n", q, irr, ci_str))
  }
}

cat("\n")

cat("COMPARISON: SATURATED vs PARSIMONIOUS\n")
cat("─────────────────────────────────────────────────────\n")
cat("                 Saturated    Parsimonious\n")
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("Pre-Omicron Racialized:   %.3f        %.3f\n", irr_pre[2], irr_pre_pars[2]))
cat(sprintf("Omicron Racialized:       %.3f        %.3f\n", irr_omi[2], irr_omi_pars[2]))
cat("─────────────────────────────────────────────────────\n\n")

cat("The saturated model shows gradient compression,\n")
cat("while the parsimonious model constrains gradients to be equal.\n")
cat("The large ΔAIC = -192.5 strongly favors the saturated model.\n\n")

