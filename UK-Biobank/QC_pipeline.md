### UKBB Study Criteria
#### Target: 130,950 White-British males, aged 40-69, no prevalent PrCa/T2DM, ≥1yr follow-up
#### Genotyping: v3 (March 2018), Affymetrix UK BiLEVE/UKBB arrays
#### Application: #90981

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(dplyr)
library(data.table)
library(ukbtools)
library(lubridate)

# === UK BIOBANK QC PIPELINE (Application #90981) ===

# Load core datasets - REPLACE WITH YOUR FILE PATHS
pheno <- ukb_df_read([PHENO_FILE_PATH])  # ukbXXXX.enc_ukb (fields 31,34,20001,20002...)
pca <- fread([PCA_FILE_PATH])            # ukb_cal_chr1_v2.fam_pca.txt (40 PCs)
king <- fread([KING_FILE_PATH])          # KING kinship output
cancer_data <- fread([CANCER_FILE_PATH]) # Field 20001 (cancer ascertained), 20002 (ICD10)
t2dm_data <- fread([T2DM_FILE_PATH])     # ICD9:250.*, ICD10:E11,E14 + lab criteria

# Reference dates (Methods: median 12.2y follow-up to 2020)
baseline_date <- as.Date("2006-01-01") + [MEDIAN_FOLLOWUP_DAYS]  # ~4463 days
censor_date <- as.Date("2020-07-01")

# STEP 1: Genotype availability (v3 pipeline, >800k SNPs)
genotype_qc <- fread([UKB_GENOTYPE_QC_PATH])  # genetic_showcase.results
pheno <- pheno %>% 
  filter(eid %in% genotype_qc$eid[genotype_qc$genotyping_array == 1])

# STEP 2: White-British ancestry (PCA cluster 1)
pheno <- pheno %>% 
  left_join(pca %>% select(eid, PC1 = genetic_principal_components.f22009_0_0), 
            by = "eid") %>%
  filter(PC1 >= [PC1_MIN_WB] & PC1 <= [PC1_MAX_WB])  # e.g., -0.01 to 0.02

# STEP 3: Sex mismatch + males only
pheno <- pheno %>% 
  filter(sex_f31_0_0 == 1 &  # Self-reported male
         genetic_sex.f22001_0_0 == 1)  # Genetic sex male

# STEP 4: Remove 2nd-degree relatives (KING kinship > 0.0884)
related_ids <- king %>% 
  filter(kinship > [KING_2ND_THRESH]) %>%  # 0.0884 = 2nd degree
  pull(fid) %>% unique()
pheno <- pheno %>% filter(!eid %in% related_ids)

# STEP 5: EXCLUDE prevalent cancer (Field 20001) & prevalent T2DM
pheno <- pheno %>%
  filter(!eid %in% cancer_data$cancer_eid[cancer_data$cancer_ascertained == 1]) %>%  # No cancer
  filter(!eid %in% t2dm_data$t2dm_eid[t2dm_data$t2dm_baseline == 1])  # No T2DM

# STEP 6: Age 40-69 at baseline + minimum 1yr follow-up
pheno <- pheno %>%
  mutate(
    yob = year_of_birth.f34_0_0,
    age_baseline = as.numeric(baseline_date - as.Date(paste0(yob, "-07-01"))) / 365.25
  ) %>%
  filter(age_baseline >= 40 & age_baseline <= 69) %>%
  filter(as.numeric(date_first_alive_follow_up - baseline_date) >= 365)  # >=1yr

# FINAL OUTPUT: 130,950 White-British males
fwrite(pheno, [UKB_CLEANED_COHORT_PATH], row.names = FALSE)
cat("UKBB Final N:", nrow(pheno), "White-British males\n")
