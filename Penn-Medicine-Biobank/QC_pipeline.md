## Penn Medicine Biobank QC Pipeline

#### PMBB Study Criteria
#### Target: 17,970 males (14,279 Non-Hispanic EUR + 3,691 AA)
#### Genotyping: GSA array (43,623 total samples)
#### EHR: Median 7.0y coverage to July 2020, IRB #813913

```markdown
---
title: "PMBB QC Pipeline - Prostate Cancer Study"
author: "[Your Name]"
date: "`r Sys.Date()`"
output: 
  html_document:
    toc: true
    toc_float: true
---

library(dplyr)
library(data.table)

# === PENN MEDICINE BIOBANK QC PIPELINE ===

# Load PMBB datasets - REPLACE WITH YOUR FILE PATHS
pmbb_pheno <- fread([PMBB_PHENO_PATH])     # EHR phenotypes
pmbb_pca <- fread([PMBB_PCA_PATH])         # Ancestry PCs
king_pmbb <- fread([PMBB_KING_PATH])       # KING output
pmbb_cancer <- fread([PMBB_CANCER_PATH])   # Prevalent cancer
pmbb_t2dm <- fread([PMBB_T2DM_PATH])       # Baseline T2DM (ICD9/10 + labs)

# Baseline: first clinical encounter (median ~2010)
baseline_date_pmbb <- as.Date([PMBB_BASELINE_DATE])  # e.g., "2010-01-01"

# STEP 1: Genotype availability (GSA array)
pmbb_pheno <- pmbb_pheno %>% 
  filter(!is.na(genotype_id) & !is.na(eid) & genotyping_batch != "failed")

# STEP 2: Ancestry - Non-Hispanic EUR OR African American ONLY
pmbb_pheno <- pmbb_pheno %>% 
  filter(ancestry_category %in% c("Non-Hispanic European", "African American"))

# STEP 3: Males only + sex concordance
pmbb_pheno <- pmbb_pheno %>% 
  filter(phenotypic_sex == "Male" & 
         phenotypic_sex == genetic_sex)

# STEP 4: Remove 2nd-degree relatives (KING)
related_pmbb <- king_pmbb %>% 
  filter(kinship > [KING_2ND_THRESH]) %>%  # 0.0884
  pull(IID) %>% unique()
pmbb_pheno <- pmbb_pheno %>% filter(!sample_id %in% related_pmbb)

# STEP 5: EXCLUDE: No medical history + prevalent cancer + prevalent T2DM
pmbb_pheno <- pmbb_pheno %>%
  filter(has_medical_history == TRUE) %>%  # ≥1 EHR record
  anti_join(pmbb_cancer %>% filter(prevalent_cancer == TRUE), by = "sample_id") %>%
  anti_join(pmbb_t2dm %>% filter(t2dm_baseline == TRUE), by = "sample_id")

# STEP 6: Age 40-69 + minimum 1yr post-enrollment follow-up
pmbb_pheno <- pmbb_pheno %>%
  filter(sample_age_years >= 40 & sample_age_years <= 69) %>%
  filter(followup_years_post_enroll >= 1.0)  # Methods: 4.2y median post-enroll

# FINAL OUTPUT: 17,970 males (14,279 EUR + 3,691 AA)
fwrite(pmbb_pheno, [PMBB_CLEANED_COHORT_PATH], row.names = FALSE)

# Ancestry breakdown
cat("PMBB Final N:", nrow(pmbb_pheno), "\n")
cat("EUR:", nrow(pmbb_pheno[ancestry_category == "Non-Hispanic European",]), "\n")
cat("AA: ", nrow(pmbb_pheno[ancestry_category == "African American",]), "\n")
