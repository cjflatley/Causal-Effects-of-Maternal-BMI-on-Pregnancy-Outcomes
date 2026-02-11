# ===============================================
# Mediation Analysis Template
# ===============================================
# Purpose: Calculate indirect, direct, and total effects for a mediator
# Author: Christopher Flatley
# Date: 2026-02-12
# Notes:
#   - Replace the effect sizes and standard errors with your own summary statistics.
#   - Works with any exposure -> mediator -> outcome combination.
# ===============================================

library(data.table)
library(dplyr)

# ----------------------------
# Input effect sizes and SEs
# ----------------------------
# Replace these with your own estimates
# Exposure -> Mediator
beta_exp_med <- -0.0404454695270945
se_exp_med <- 0.0184042763170435

# Mediator -> Outcome
beta_med_out <- -0.0273955543184189
se_med_out <- 0.0105728620453891

# Total effect of exposure -> outcome
beta_total <- 0.04372567
se_total <- 0.0156195

# ----------------------------
# Calculate indirect effect
# ----------------------------
beta_indirect <- beta_exp_med * beta_med_out

# Standard error for indirect effect
se_indirect <- sqrt(
  (beta_exp_med * se_med_out)^2 +
  (beta_med_out * se_exp_med)^2 +
  (se_exp_med * se_med_out)^2
)

# Confidence intervals
ci_lower_indirect <- beta_indirect - 1.96 * se_indirect
ci_upper_indirect <- beta_indirect + 1.96 * se_indirect

# Z-score and p-value
z_indirect <- beta_indirect / se_indirect
p_indirect <- 2 * pnorm(-abs(z_indirect))

# Create data frame
df_indirect <- data.frame(
  Beta = beta_indirect,
  SE = se_indirect,
  CI_lower = ci_lower_indirect,
  CI_upper = ci_upper_indirect,
  P = p_indirect,
  Effect = "Indirect Effect"
)

# ----------------------------
# Calculate direct effect
# ----------------------------
beta_direct <- beta_total - beta_indirect
se_direct <- sqrt(se_total^2 + se_indirect^2)

ci_lower_direct <- beta_direct - 1.96 * se_direct
ci_upper_direct <- beta_direct + 1.96 * se_direct

z_direct <- beta_direct / se_direct
p_direct <- 2 * pnorm(-abs(z_direct))

df_direct <- data.frame(
  Beta = beta_direct,
  SE = se_direct,
  CI_lower = ci_lower_direct,
  CI_upper = ci_upper_direct,
  P = p_direct,
  Effect = "Direct Effect"
)

# ----------------------------
# Total effect
# ----------------------------
ci_lower_total <- beta_total - 1.96 * se_total
ci_upper_total <- beta_total + 1.96 * se_total

z_total <- beta_total / se_total
p_total <- 2 * pnorm(-abs(z_total))

df_total <- data.frame(
  Beta = beta_total,
  SE = se_total,
  CI_lower = ci_lower_total,
  CI_upper = ci_upper_total,
  P = p_total,
  Effect = "Total Effect"
)

# ----------------------------
# Combine results
# ----------------------------
df_mediation <- bind_rows(df_indirect, df_direct, df_total)

# ----------------------------
# Save results
# ----------------------------
write.table(
  df_mediation,
  "results/mediation/mediation_effects.txt",  # Replace with your desired output path
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  na = ""
)

message("Mediation analysis complete. Results saved to 'results/mediation/mediation_effects.txt'.")
# ===============================================
# End of Template Script
# ===============================================
