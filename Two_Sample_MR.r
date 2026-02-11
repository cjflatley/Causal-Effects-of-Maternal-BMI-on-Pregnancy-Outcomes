# ===============================================
# Two-sample MR + Mediation Template Script
# ===============================================
# Purpose: Perform two-sample MR for multiple outcomes using a generic exposure.
# Author: Christopher Flatley
# Date: 2026-02-12
# Notes:
#   - Replace file paths with your own data.
#   - Designed to work with summary statistics (no individual-level data included).
#   - Outputs MR results, heterogeneity, pleiotropy, leave-one-out, single SNP results, and scatter plots.
# ===============================================

# ----------------------------
# Load required packages
# ----------------------------
library(data.table)
library(dplyr)
library(TwoSampleMR)       # For MR analysis
library(MRInstruments)     # Optional: provides curated instruments
library(ggplot2)           # For plotting

# ----------------------------
# Load exposure summary statistics
# ----------------------------
# Replace with path to your exposure GWAS summary stats
exposure_file <- "data/exposure.txt.gz"

exposure <- fread(exposure_file)

# Rename columns to standard MR names
exposure <- exposure %>% rename(
  beta = BETA,                  # effect size
  se = SE,                      # standard error
  eaf = Freq_Tested_Allele_in_HRS,  # effect allele frequency
  pval = P,                     # p-value
  effect_allele = Tested_Allele,
  other_allele = Other_Allele,
  samplesize = N
)

# Format exposure data for TwoSampleMR
exposure$Phenotype <- "BMI"  # Replace with generic exposure 
exposure_data <- format_data(exposure)

# Optionally filter genome-wide significant SNPs
exposure_data <- exposure_data %>% filter(pval.exposure < 5e-08)

# ----------------------------
# Clump exposure data for independent SNPs
# ----------------------------
# Adjust clumping parameters as needed
exposure_data <- clump_data(
  exposure_data,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR",  # adjust for your population
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "data/1kg_ref/EUR"  # reference panel path
)

# Calculate F-statistics for exposure instruments
exposure_data$F.exposure <- (exposure_data$beta.exposure / exposure_data$se.exposure)^2

# ----------------------------
# Define outcomes
# ----------------------------
# Replace with your outcome GWAS files or names
outcomes <- c("Basophil", "Eosinophil", "Lymphocyte", "Monocyte", "Neutrophil", "Platelet")

# ----------------------------
# Loop over outcomes
# ----------------------------
for (outcome_name in outcomes) {

  # Load outcome data
  outcome_file <- paste0("data/", outcome_name, ".txt")  # Replace with your path
  outcome <- fread(outcome_file)
  
  # Standardize columns
  outcome <- outcome %>% rename(
    SNP = variant_id,
    beta = beta,
    se = standard_error,
    pos = base_pair_location,
    eaf = effect_allele_frequency,
    pval = p_value,
    samplesize = N
  )
  
  outcome$Phenotype <- outcome_name
  
  # Format for TwoSampleMR
  outcome_data <- outcome %>% rename(
    beta.outcome = beta,
    se.outcome = se,
    samplesize.outcome = samplesize,
    pval.outcome = pval,
    eaf.outcome = eaf,
    effect_allele.outcome = effect_allele,
    other_allele.outcome = other_allele
  )
  
  outcome_data$outcome <- outcome_name
  outcome_data$id.outcome <- "id"

  # ----------------------------
  # Harmonize exposure and outcome
  # ----------------------------
  harmonised_data <- harmonise_data(
    exposure_dat = exposure_data,
    outcome_dat = outcome_data,
    action = 1  # automatically align alleles
  )
  
  # ----------------------------
  # Define output directory
  # ----------------------------
  output_dir <- paste0("results/MR/", outcome_name, "/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Save F-statistics
  f_stat_summary <- capture.output(summary(harmonised_data$F.exposure))
  writeLines(f_stat_summary, paste0(output_dir, "F_stat_summary.txt"))

  # ----------------------------
  # Run MR analyses
  # ----------------------------
  mr_results <- mr(harmonised_data)
  heterogeneity <- mr_heterogeneity(harmonised_data)
  pleiotropy <- mr_pleiotropy_test(harmonised_data)
  leave_one_out <- mr_leaveoneout(harmonised_data)
  single_snp <- mr_singlesnp(harmonised_data)

  # ----------------------------
  # Save MR results
  # ----------------------------
  write.table(as.data.frame(mr_results), paste0(output_dir, "mr_results.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(as.data.frame(heterogeneity), paste0(output_dir, "heterogeneity.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(as.data.frame(pleiotropy), paste0(output_dir, "pleiotropy.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(as.data.frame(leave_one_out), paste0(output_dir, "leave_one_out.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(as.data.frame(single_snp), paste0(output_dir, "single_snp.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")

  # ----------------------------
  # Generate scatter plot
  # ----------------------------
  scatter_plot <- mr_scatter_plot(mr_results, harmonised_data)
  
  png(paste0(output_dir, "scatter_plot_", outcome_name, ".png"), width = 900, height = 900)
  print(scatter_plot)
  dev.off()
  
  message("MR analysis complete for: ", outcome_name)
}

# ===============================================
# End of Template Script
# ===============================================
