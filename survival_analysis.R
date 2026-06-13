 #Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required packages
BiocManager::install(c("TCGAbiolinks", "EDASeq"), update = FALSE, ask = FALSE)
install.packages(c("survival", "survminer", "tidyverse"), dependencies = TRUE)

# Load libraries
library(TCGAbiolinks)
library(survival)
library(survminer)
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)

# Fetching clinical data for Lung Adenocarcinoma (LUAD)
query_clinical <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "Clinical", 
  data.type = "Clinical Supplement"
)

# 2. Download the XML metadata files
GDCdownload(query_clinical)
df <- query_clinical$results[[1]]
df <- df[!grepl("\\.txt$", df$file_name), ]
query_clinical$results[[1]] <- df

# 3. Prepare and parse the data into a clean dataframe
# This safely parses the XML elements without throwing the row-mismatch error

clinical_data <- GDCprepare_clinic(query_clinical, clinical.info = "patient")

# Take a quick look to verify it worked
head(clinical_data[, c("vital_status", "days_to_death", "days_to_last_followup")])

any(colnames(clinical_data) %in% c("vital_status", "days_to_death", "days_to_last_followup"))

clinical_data$vital_status[is.na(clinical_data$vital_status) & !is.na(clinical_data$days_to_death)] <- "Dead"
clinical_data$vital_status[is.na(clinical_data$vital_status) & !is.na(clinical_data$days_to_last_known_alive)] <- "Alive"
clinical_data$vital_status[is.na(clinical_data$vital_status) & !is.na(clinical_data$days_to_last_followup)] <- "Alive"

clinical_data <- clinical_data %>% drop_na(vital_status)

clinical_data$deceased <- ifelse(clinical_data$vital_status == "Alive", FALSE, TRUE)

clinical_data$overall_survival <- ifelse(clinical_data$vital_status == "Alive",
                                         clinical_data$days_to_last_followup,
                                         clinical_data$days_to_death)


# get gene expression data -----------

# build a query to get gene expression data for entire cohort
query_brca_all = GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")

rnaseq_data <- getResults(query_brca_all)

tumor <- rnaseq_data$cases[1:200]

# get gene expression data from 200 primary tumors 
query_brca <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open",
  barcode = tumor)

# download data
GDCdownload(query_brca)

# get counts
tcga_brca_data <- GDCprepare(query_brca, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, "unstranded")

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_brca_data))

coldata <- as.data.frame(colData(tcga_brca_data))

# vst transform counts to be used in survival analysis ---------------
# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = brca_matrix,
                              colData = coldata,
                              design = ~ 1)
# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# vst 
vsd <- vst(dds, blind=FALSE)
brca_matrix_vst <- assay(vsd)

brca_egfr <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "TCN1")

# get median value
median_value <- median(brca_egfr$counts)

# denote which cases have higher or lower expression than median count
brca_egfr$strata <- ifelse(brca_egfr$counts >= median_value, "HIGH", "LOW")

# Add clinical information to brca_tp53
colnames(clinical_data)[1] <- "case_id"
brca_egfr$case_id <- gsub('-01.*', '', brca_egfr$case_id)
brca_egfr <- merge(brca_egfr, clinical_data, by.x = 'case_id', by.y = 'case_id')

brca_egfr_clean <- brca_egfr[ !is.na(brca_egfr$overall_survival) & 
                                !is.na(brca_egfr$deceased) & 
                                !is.na(brca_egfr$strata), ]

# fitting survival curve -----------
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = brca_egfr_clean)

brca_egfr_clean$strata <- factor(brca_egfr_clean$strata, levels = c("LOW", "HIGH"))
cox_fit <- coxph(Surv(overall_survival, deceased) ~ strata, data = brca_egfr_clean)
cox_summary <- summary(cox_fit)

# 3. Extract HR, Confidence Intervals, and Wald p-value
hr <- round(cox_summary$conf.int[1], 2)
hr_lower <- round(cox_summary$conf.int[3], 2)
hr_upper <- round(cox_summary$conf.int[4], 2)
wald_p <- format.pval(cox_summary$waldtest[3], digits = 3)

# 4. Create a clean string for the plot annotation
hr_string <- paste0("HR = ", hr, " (95% CI: ", hr_lower, "-", hr_upper, ")", "\nCox p = ", wald_p)

p <- ggsurvplot(fit,
                data = brca_egfr_clean,
                pval = FALSE,        # Turned off default log-rank p-value to avoid layout clutter
                risk.table = TRUE,
                conf.int = TRUE,
                title = "TCGA-LUAD Survival Analysis for TCN1",
                xlab = "Time (Days)",
                legend.labs = c("LOW", "HIGH"))

p$plot <- p$plot + 
  theme(plot.title = element_text(hjust = 0.5))

p$plot <- p$plot +
  annotate(
    "text",
    x = 5000, y = 0.9,              # Coordinates on the plot area
    label = hr_string,
    hjust = 0,                      # Left-aligns the text block
    size = 4.5,                     # Font size
    fontface = "bold.italic",
    color = "black"
  )

print(p)


fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = brca_egfr_clean)
print(fit2)