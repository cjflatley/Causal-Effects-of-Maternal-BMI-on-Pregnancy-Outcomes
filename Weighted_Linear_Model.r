# ===============================================
# Maternal-Fetal Weighted Linear Model (WLM) Adjustment Template
# ===============================================
# Purpose: Adjust maternal and fetal SNP effects using the WLM method
# Author: Christopher Flatley
# Date: 2026-02-12
# Notes:
#   - Replace file paths with your own data.
#   - Works on summary statistics only (no individual-level data).
#   - Calculates adjusted maternal and fetal effects, standard errors, Z-scores, and p-values.
# ===============================================

# ----------------------------
# Load required packages
# ----------------------------
library(data.table)
library(dplyr)

# ----------------------------
# Load input summary statistics
# ----------------------------
# Maternal and fetal GWAS data
maternal_file <- "data/maternal_gwas.txt.gz"  # Replace with actual maternal GWAS summary stats
fetal_file <- "data/fetal_gwas.txt.gz"        # Replace with actual fetal GWAS summary stats

# SNPs from exposure of interest (e.g., BMI)
exposure_snps_file <- "data/exposure_snps.txt"

# Read data
maternal_gwas <- fread(maternal_file)
fetal_gwas <- fread(fetal_file)
exposure_snps <- fread(exposure_snps_file)

# ----------------------------
# Prepare exposure SNP list
# ----------------------------
exposure_snps <- exposure_snps %>% 
  select(SNP) %>%        # Keep only SNP column
  distinct(SNP, .keep_all = TRUE) %>%  # Remove duplicates
  mutate(check = "EXPOSURE")           # Generic identifier

# ----------------------------
# Prepare fetal GWAS
# ----------------------------
fetal_gwas <- fetal_gwas %>% 
  rename(
    SNP = rsID,
    CHR = Chr,
    POS = Pos,
    effect_allele = A1,
    other_allele = A0
  )

# Keep only SNPs in exposure
fetal_gwas <- left_join(exposure_snps, fetal_gwas, by = "SNP")

# ----------------------------
# Prepare maternal GWAS
# ----------------------------
maternal_gwas <- maternal_gwas %>% 
  rename(
    SNP = rsID,
    CHR = Chr,
    POS = Pos,
    effect_allele = A1,
    other_allele = A0
  )

# Keep only SNPs in exposure
maternal_gwas <- left_join(exposure_snps, maternal_gwas, by = "SNP")

# ----------------------------
# Calculate weighted allele frequency for maternal GWAS
# ----------------------------
n_maternal <- 210266  # replace with your sample size
n_fetal <- 59736      # replace with your sample size

maternal_gwas$eaf <- (maternal_gwas$`IS-frq` * n_fetal + maternal_gwas$`EGG-frq` * n_maternal) / (n_maternal + n_fetal)
# Handle missing EAF
maternal_gwas$eaf <- ifelse(is.na(maternal_gwas$eaf) & !is.na(maternal_gwas$`IS-frq`), maternal_gwas$`IS-frq`, maternal_gwas$eaf)
maternal_gwas$eaf <- maternal_gwas$eaf / 100  # Convert percentage to fraction

# ----------------------------
# Rename columns for clarity
# ----------------------------
maternal_gwas <- maternal_gwas %>% rename(
  beta_m = `Beta-A1`,
  p_m = P,
  eaf_m = eaf
) %>% select(SNP, effect_allele, other_allele, beta_m, p_m, eaf_m)

fetal_gwas <- fetal_gwas %>% rename(
  beta_f = `Beta-A1`,
  p_f = P
) %>% select(SNP, effect_allele, other_allele, beta_f, p_f)

# ----------------------------
# Merge maternal and fetal GWAS
# ----------------------------
wlm <- full_join(maternal_gwas, fetal_gwas, by = c("SNP", "effect_allele", "other_allele"))

# ----------------------------
# Calculate adjusted maternal and fetal betas
# ----------------------------
wlm$beta_maternal_adjusted <- ((4/3) * wlm$beta_m) - ((2/3) * wlm$beta_f)
wlm$beta_fetal_adjusted <- ((-2/3) * wlm$beta_m) + ((4/3) * wlm$beta_f)

# ----------------------------
# Calculate Z-scores and standard errors
# ----------------------------
wlm$z_m <- qnorm(1 - wlm$p_m / 2)
wlm$z_f <- qnorm(1 - wlm$p_f / 2)

wlm$se_m <- wlm$beta_m / wlm$z_m
wlm$se_f <- wlm$beta_f / wlm$z_f

# Covariance intercept (from LD Score regression)
cm <- 0.1722

# Calculate SE for adjusted maternal and fetal effects
wlm$wlm_se_m <- sqrt((16 * wlm$se_m^2) + 4 * wlm$se_f^2 - 16 * cm * sqrt(wlm$se_m^2 * wlm$se_f^2)) / 3
wlm$wlm_se_f <- sqrt((16 * wlm$se_f^2) + 4 * wlm$se_m^2 - 16 * cm * sqrt(wlm$se_m^2 * wlm$se_f^2)) / 3

# Calculate maternal Z-score and two-tailed p-value
wlm$wlm_z_score <- wlm$beta_maternal_adjusted / wlm$wlm_se_m
wlm$wlm_p_value <- 2 * (1 - pnorm(abs(wlm$wlm_z_score)))

# ----------------------------
# Select output columns
# ----------------------------
wlm_output <- wlm %>% select(
  SNP, effect_allele, other_allele,
  beta_m, p_m, eaf_m,
  beta_maternal_adjusted, wlm_se_m, wlm_p_value
)

# ----------------------------
# Save results
# ----------------------------
write.table(
  wlm_output,
  "results/wlm_adjusted_beta.txt",  # Replace with your desired output path
  row.names = FALSE,
  na = "",
  quote = FALSE,
  sep = "\t"
)

message("WLM maternal-fetal adjustment complete.")
# ===============================================
# End of Template Script
# ===============================================
