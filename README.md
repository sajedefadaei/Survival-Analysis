# Survival Analysis of TCN1 Expression in TCGA-LUAD

## Overview

This project investigates the prognostic relevance of *TCN1* (transcobalamin 1) gene expression in lung adenocarcinoma (LUAD) using publicly available data from The Cancer Genome Atlas (TCGA). RNA-sequencing-derived gene expression and clinical outcome data for primary tumor samples are integrated, patients are stratified by *TCN1* expression level, and survival differences between strata are assessed using Kaplan–Meier estimation and Cox proportional hazards regression.

## Data Source

All data are retrieved programmatically from the NCI Genomic Data Commons (GDC) via the `TCGAbiolinks` package, using the **TCGA-LUAD** (Lung Adenocarcinoma) project:

- **Clinical data**: Clinical Supplement, parsed at the patient level (vital status, days to death, days to last follow-up).
- **Gene expression data**: RNA-Seq, STAR-Counts workflow, Gene Expression Quantification, restricted to primary tumor samples (open access). A cohort of 200 primary tumor cases is used.

## Methodology

**1. Clinical data retrieval and survival variable construction.** Clinical metadata are queried and parsed using `GDCprepare_clinic`. Vital status is reconciled across multiple fields (`vital_status`, `days_to_death`, `days_to_last_known_alive`, `days_to_last_followup`), and patients with missing vital status are removed. A binary `deceased` indicator and a combined `overall_survival` time variable (days to death for deceased patients, days to last follow-up for living patients) are derived for use in survival modeling.

**2. Gene expression data retrieval and preprocessing.** RNA-Seq count data for 200 TCGA-LUAD primary tumor samples are downloaded and prepared as a `SummarizedExperiment` object. Raw counts are extracted along with gene-level and sample-level metadata.

**3. Variance-stabilizing transformation (VST).** Genes with fewer than 10 total reads across all samples are filtered out. The remaining count matrix is normalized using `DESeq2`'s variance-stabilizing transformation (`vst`), producing an expression matrix suitable for downstream comparison across samples.

**4. Gene-of-interest extraction and stratification.** Expression values for *TCN1* are extracted from the VST-normalized matrix and reshaped to long format. Patients are stratified into **HIGH** and **LOW** expression groups based on whether their *TCN1* expression is above or below the cohort median.

**5. Integration with clinical outcomes.** Expression-based strata are merged with the clinical survival data using harmonized case identifiers (TCGA barcodes truncated to the patient-level identifier). Records with missing survival time, event status, or strata are excluded from the final analysis dataset.

**6. Survival modeling.**
- A Kaplan–Meier survival curve is fitted via `survfit()`, stratified by *TCN1* expression group (HIGH vs. LOW).
- A Cox proportional hazards model is fitted via `coxph()` to estimate the hazard ratio (HR), 95% confidence interval, and Wald test p-value for the association between *TCN1* expression strata and overall survival.
- A log-rank test (`survdiff()`) is additionally computed to test for a difference in survival distributions between strata.

**7. Visualization.** The Kaplan–Meier curve is plotted using `survminer::ggsurvplot()`, including a risk table and confidence intervals, with the Cox model's hazard ratio, confidence interval, and p-value annotated directly on the plot.

## Required Packages

```r
# Bioconductor packages
BiocManager::install(c("TCGAbiolinks", "EDASeq", "SummarizedExperiment", "DESeq2"))

# CRAN packages
install.packages(c("survival", "survminer", "tidyverse"))
```

## Usage

1. Ensure all required packages (above) are installed.
2. Run the script sequentially in R or RStudio. Note that the `GDCdownload()` step retrieves RNA-Seq data for 200 samples from the GDC and may take a significant amount of time and disk space depending on connection speed.
3. The final output is a Kaplan–Meier survival plot with an accompanying risk table, annotated with the Cox proportional hazards estimate for *TCN1* expression strata.

## Output

The script produces a stratified Kaplan–Meier survival plot ("TCGA-LUAD Survival Analysis for TCN1") comparing overall survival between patients with high versus low *TCN1* expression, with the hazard ratio, 95% confidence interval, and Cox model p-value displayed on the plot.

## Known Issues / Notes

- Line `query_clinical$results[[1]] <- df` assigns an undefined variable `df` and will raise an error as written; clinical results should instead be assigned from `query_clinical$results[[1]]` itself or the relevant downloaded object before calling `GDCprepare_clinic()`.
- The final `survdiff()` result (`fit2`) is computed but not explicitly printed or reported alongside the Cox model results; consider adding `print(fit2)` or extracting its p-value for direct comparison with the Cox model's Wald p-value.
- The cohort is limited to the first 200 primary tumor cases returned by the GDC query rather than the full TCGA-LUAD cohort; results should be interpreted as exploratory given this reduced sample size.

## Repository Structure

```
.
├── survival_analysis.R
└── README.md
```
