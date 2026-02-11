This repository contains R scripts used to investigate the causal effects of maternal body mass index (BMI) on pregnancy outcomes using publicly available summary statistics. 
The study applied Mendelian randomisation (MR) to assess the causal relationships between maternal BMI and pregnancy outcomes, including birth weight, placental weight, gestational duration, and pre-eclampsia.

To further explore potential mechanisms, two-step MR and mediation analyses were performed to evaluate whether immune-related blood cell counts, such as neutrophils, lymphocytes, and platelets, mediate these relationships. 
Summary statistics for maternal BMI and pregnancy outcomes were obtained from public genome-wide association studies (GWAS), with maternal genetic effects used to proxy intrauterine environmental influences.

Two-sample Mendelian Randomisation (MR): estimating causal effects of maternal BMI on blood cell traits.

Maternal-fetal Weighted Linear Model (WLM) adjustment: separating maternal and fetal genetic contributions.

Mediation analysis: decomposing total effects into direct and indirect pathways through intermediate traits.

The code is designed for reproducible research and accompanies a published manuscript:
"Causal Effects of Maternal BMI on Pregnancy Outcomes: A Mendelian Randomisation Study Investigating the Mediating Role of Blood Counts."

No individual-level data are included; only summary statistics are required. 
Researchers can adapt these scripts to their own datasets for similar causal analyses.

## Scripts

### 1. `Two_Sample_MR.r`
- Performs **two-sample MR** for a set of exposure SNPs on multiple blood cell outcomes.
- Uses the **TwoSampleMR** R package.
- Outputs:
  - MR results, heterogeneity tests, pleiotropy tests
  - Leave-one-out and single SNP analyses
  - Scatter plots for each outcome

**Usage:**  
Update the paths to your exposure and outcome GWAS summary statistics and run the script.

---

### 2. `Weighted_Linear_Model.r`
- Performs **maternal-fetal Weighted Linear Model (WLM)** adjustment.
- Inputs maternal and fetal summary statistics and SNPs from the exposure of interest.
- Calculates:
  - Adjusted maternal and fetal betas
  - Standard errors, Z-scores, p-values
- Outputs a table of adjusted maternal effects ready for downstream analyses.

**Usage:**  
Provide maternal/fetal GWAS summary statistics and exposure SNP list. 
Update sample sizes if different from defaults.

---

### 3. `Mediation_Analysis.r`
- Performs **mediation analysis** for an exposure → mediator → outcome pathway.
- Calculates:
  - Indirect effect (exposure → mediator → outcome)
  - Direct effect (remaining effect of exposure on outcome)
  - Total effect (exposure → outcome)
- Provides beta, standard error, confidence intervals, Z-scores, and p-values.

**Usage:**  
Replace the example effect sizes and standard errors with your own summary statistics.
