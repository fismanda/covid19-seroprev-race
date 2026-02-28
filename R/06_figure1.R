# ==============================================================================
# Figure 1: Infection-Induced SARS-CoV-2 Seroprevalence by Racial Group
# ==============================================================================
# Generates publication-quality figure for BMC Public Health.
#
# BMC figure specs:
#   - Single column: 85 mm (~3.35 in); 1.5 column: 120 mm (~4.72 in)
#   - Full width: 170 mm (~6.69 in)
#   - Min resolution: 300 dpi (halftone/colour), 1000 dpi (line art)
#   - Formats: TIFF, PNG, PDF, EPS
#   - Font: Arial preferred, 8-12 pt
#
# This script uses full-width format with Arial font.
# ==============================================================================

# Set working directory to project root (adjust as needed)
# setwd("path/to/covid19-seroprev-race")

library(ggplot2)
library(dplyr)
library(lubridate)

source("R/00_helpers.R")

cat("\n=== GENERATING FIGURE 1: MODEL FIT (BMC FORMAT) ===\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA
# -----------------------------------------------------------------------------

ses_data <- read.csv("data/seroprevalence_by_social_determinant.csv")
ses_data$date <- mdy(ses_data$samplingdate)

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

# Create clean race labels
cbs_data$race_clean <- case_when(
  cbs_data$population == "white" ~ "White",
  cbs_data$population == "racialized" ~ "Racialized",
  TRUE ~ NA_character_
)

cat("Racial groups in data:\n")
print(table(cbs_data$race_clean, useNA = "ifany"))

# Filter to fitting period
fitting_start_date <- as.Date("2020-06-14")
plot_data_obs <- cbs_data %>%
  filter(date >= fitting_start_date,
         !is.na(race_clean)) %>%
  select(date, seroprev_est, race_clean)

cat(sprintf("Loaded %d observations from %s to %s\n",
            nrow(plot_data_obs), min(plot_data_obs$date), max(plot_data_obs$date)))

# -----------------------------------------------------------------------------
# 2. LOAD MODEL AND EXTRACT PARAMETERS
# -----------------------------------------------------------------------------

results <- readRDS("output/saturated_model_estimate_rho.rds")
coefs <- coef(results$fit)

lambda_white_pre <- plogis(coefs["logit_lambda_white_pre"])
lambda_racialized_pre <- plogis(coefs["logit_lambda_racialized_pre"])
lambda_white_omi <- plogis(coefs["logit_lambda_white_omi"])
lambda_racialized_omi <- plogis(coefs["logit_lambda_racialized_omi"])

rho <- results$rho
omicron_date <- results$omicron_date
model_start_date <- results$model_start_date

cat(sprintf("Parameters loaded: ρ = %.5f/week\n\n", rho))

# -----------------------------------------------------------------------------
# 3. SIMULATE MODEL PREDICTIONS
# -----------------------------------------------------------------------------

race <- c("White", "Racialized")
lambdas_pre <- c(lambda_white_pre, lambda_racialized_pre)
lambdas_omi <- c(lambda_white_omi, lambda_racialized_omi)

prediction_data <- data.frame()

for (i in 1:2) {
  obs_dates <- plot_data_obs %>%
    filter(race_clean == race[i]) %>%
    pull(date)

  if (length(obs_dates) > 0) {
    date_seq <- seq(as.Date("2020-06-14"), max(obs_dates), by = "7 days")

    # simulate_si_saturated returns proportions; multiply by 100 for %
    pred <- simulate_si_saturated(lambdas_pre[i], lambdas_omi[i], rho,
                                  date_seq, model_start_date, omicron_date) * 100

    pred_df <- data.frame(
      date = date_seq,
      predicted = pred,
      race_clean = race[i]
    )

    prediction_data <- rbind(prediction_data, pred_df)
    cat(sprintf("Generated predictions for %s: %d timepoints\n",
                race[i], length(date_seq)))
  }
}

# -----------------------------------------------------------------------------
# 4. SET UP FACTOR LEVELS
# -----------------------------------------------------------------------------

plot_data_obs$race <- factor(plot_data_obs$race_clean,
                             levels = c("White", "Racialized"))

prediction_data$race <- factor(prediction_data$race_clean,
                               levels = c("White", "Racialized"))

# -----------------------------------------------------------------------------
# 5. CREATE PUBLICATION FIGURE
# -----------------------------------------------------------------------------

# BMC-friendly colour palette (colourblind-safe)
race_colours <- c("White" = "#2166AC", "Racialized" = "#B2182B")

# Font: BMC prefers Arial. On systems with Arial installed, change to "Arial".
# "sans" maps to the system default sans-serif (usually Helvetica or DejaVu Sans).
fig_font <- "sans"

p <- ggplot() +
  # Model predictions (smooth lines)
  geom_line(data = prediction_data,
            aes(x = date, y = predicted, colour = race),
            linewidth = 0.8) +
  # Observed data (points)
  geom_point(data = plot_data_obs,
             aes(x = date, y = seroprev_est, colour = race),
             size = 1.5, alpha = 0.7) +
  # Omicron emergence line
  geom_vline(xintercept = as.numeric(omicron_date),
             linetype = "dashed", colour = "gray40", linewidth = 0.5) +
  annotate("text", x = omicron_date + 10, y = 82,
           label = "Omicron\nemergence",
           hjust = 0, size = 2.8, colour = "gray30",
           family = fig_font, lineheight = 0.9) +
  # Colour scale
  scale_colour_manual(values = race_colours) +
  # Axes
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y",
               limits = c(as.Date("2020-05-01"), as.Date("2023-05-01")),
               expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 15),
                     expand = c(0.02, 0)) +
  labs(
    x = "Date",
    y = "Seroprevalence (%)",
    colour = "Racial group"
  ) +
  theme_bw(base_size = 10, base_family = fig_font) +
  theme(
    # Panel
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.border = element_rect(colour = "black", linewidth = 0.5),
    # Legend
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.2, "cm"),
    legend.margin = margin(t = -5),
    # Axes
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),
    # Margins
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  ) +
  guides(colour = guide_legend(nrow = 1))

# -----------------------------------------------------------------------------
# 6. SAVE IN MULTIPLE FORMATS
# -----------------------------------------------------------------------------

# BMC full-width figure: 170 mm wide
fig_width_mm  <- 170
fig_height_mm <- 110
fig_width_in  <- fig_width_mm / 25.4
fig_height_in <- fig_height_mm / 25.4

# TIFF (BMC preferred for raster)
ggsave("output/figures/Figure1_Model_Fit.tiff", p,
       width = fig_width_in, height = fig_height_in,
       dpi = 300, compression = "lzw", bg = "white")

# PNG (alternative raster)
ggsave("output/figures/Figure1_Model_Fit.png", p,
       width = fig_width_in, height = fig_height_in,
       dpi = 300, bg = "white")

# PDF (vector, useful for review and proofs)
ggsave("output/figures/Figure1_Model_Fit.pdf", p,
       width = fig_width_in, height = fig_height_in,
       bg = "white")

cat("\n✓ Figure 1 saved to output/figures/:\n")
cat(sprintf("  - Figure1_Model_Fit.tiff  (%.0f × %.0f mm, 300 dpi)\n",
            fig_width_mm, fig_height_mm))
cat("  - Figure1_Model_Fit.png\n")
cat("  - Figure1_Model_Fit.pdf\n")

# -----------------------------------------------------------------------------
# 7. MODEL FIT STATISTICS
# -----------------------------------------------------------------------------

cat("\n=== MODEL FIT QUALITY ===\n\n")

for (i in 1:2) {
  obs <- plot_data_obs %>%
    filter(race_clean == race[i])

  if (nrow(obs) > 0) {
    pred <- simulate_si_saturated(lambdas_pre[i], lambdas_omi[i], rho,
                                  obs$date, model_start_date, omicron_date) * 100

    ss_res <- sum((obs$seroprev_est - pred)^2)
    ss_tot <- sum((obs$seroprev_est - mean(obs$seroprev_est))^2)
    r2 <- 1 - (ss_res / ss_tot)
    rmse <- sqrt(mean((obs$seroprev_est - pred)^2))

    cat(sprintf("%s: R² = %.4f, RMSE = %.2f%%, n = %d\n",
                race[i], r2, rmse, nrow(obs)))
  }
}

cat("\n✓ Done!\n\n")
