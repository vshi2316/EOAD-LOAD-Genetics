# ==============================================================================
# A4/HABS/AIBL Independent Cohort Validation Analysis
# ==============================================================================
#
# Analysis Contents:
#   Part 1: Environment Setup
#   Part 2: A4 Cohort - Data Integration and APOE Validation
#   Part 3: HABS Cohort - pTau217 and WMH Mechanism Validation
#   Part 4: AIBL Cohort - Survival Analysis
#   Part 5: Multi-Group SEM Analysis (lavaan, MLR, FIML)
#   Part 6: Clinical Utility Assessment (Firth, AUC, AUPRC, NRI, IDI, DCA)
#   Part 7: AIBL Cox Regression and AFT Models
#   Part 8: Visualization
#   Part 9: Results Summary
#
# Methods Reference:
#   - Multi-group SEM: lavaan with MLR estimator, FIML for missing data
#   - Clinical utility: Firth-corrected logistic regression, DeLong AUC
#   - AIBL survival: Cox regression and AFT with Weibull distribution
#
# ==============================================================================

rm(list = ls())
gc()

# ==============================================================================
# Part 1: Environment Setup
# ==============================================================================

packages_cran <- c("data.table", "dplyr", "tidyr", "ggplot2", "cowplot",
                   "stringr", "survival", "survminer", "pROC", "lavaan",
                   "lme4", "lmerTest", "semTools", "logistf", "PRROC", "readxl")

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, quiet = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

dir.create("results/", showWarnings = FALSE)
dir.create("results/figures/", showWarnings = FALSE)
dir.create("results/tables/", showWarnings = FALSE)

## Utility Functions
find_column <- function(df, patterns, ignore_case = TRUE) {
  for (pattern in patterns) {
    matches <- grep(pattern, colnames(df), value = TRUE, ignore.case = ignore_case)
    if (length(matches) > 0) return(matches[1])
  }
  return(NULL)
}

safe_extract <- function(params_df, label_name, col_name) {
  idx <- which(params_df$label == label_name)
  if (length(idx) > 0) return(params_df[[col_name]][idx[1]])
  return(NA_real_)
}

theme_nc <- theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

baseline_codes <- c("bl", "BL", "sc", "SC", "screening", "baseline", "0", "1")

# ==============================================================================
# Part 2: A4 Cohort - Data Integration and APOE Validation
# ==============================================================================

## Section 2.1: Load A4 Subject Information (APOE from A4_SUBJINFO)
a4_subjinfo <- fread("A4_SUBJINFO_PRV2_08Oct2025.csv")
a4_subjinfo <- a4_subjinfo[!duplicated(BID)]
a4_subjinfo$APOE4_Carrier <- as.integer(grepl("4", a4_subjinfo$APOEGN))

count_apoe4 <- function(geno) {
  if (is.na(geno) || geno == "") return(NA)
  alleles <- unlist(strsplit(as.character(geno), "/"))
  sum(alleles == "4", na.rm = TRUE)
}
a4_subjinfo$APOE4_Count <- sapply(a4_subjinfo$APOEGN, count_apoe4)

## Section 2.2: Load A4 Demographics
a4_demog <- fread("A4_PTDEMOG_PRV2_08Oct2025.csv")
a4_demog_bl <- a4_demog[VISCODE %in% baseline_codes]
if (nrow(a4_demog_bl) == 0) a4_demog_bl <- a4_demog[, .SD[1], by = BID]
a4_demog_bl <- a4_demog_bl[!duplicated(BID)]

## Section 2.3: Load A4 PACC Data
a4_pacc_raw <- fread("A4_SPPACC_PRV2_08Oct2025.csv")
a4_pacc_total <- a4_pacc_raw[PACCQSNUM == "PACCTS"]
a4_pacc_bl <- a4_pacc_total[VISCODE %in% baseline_codes]
if (nrow(a4_pacc_bl) == 0) a4_pacc_bl <- a4_pacc_total[, .SD[1], by = BID]
a4_pacc_bl <- a4_pacc_bl[!duplicated(BID)]
a4_pacc_data <- data.frame(BID = a4_pacc_bl$BID, PACC = as.numeric(a4_pacc_bl$PACCRN))

## Section 2.4: Load A4 Amyloid PET (Centiloid from Composite_Summary)
a4_pet <- fread("A4_PETSUVR_PRV2_08Oct2025.csv")
a4_pet_composite <- a4_pet[brain_region == "Composite_Summary"]
a4_pet_bl <- a4_pet_composite[VISCODE %in% baseline_codes]
if (nrow(a4_pet_bl) == 0) a4_pet_bl <- a4_pet_composite[, .SD[1], by = BID]
a4_pet_bl <- a4_pet_bl[!duplicated(BID)]
a4_pet_data <- data.frame(BID = a4_pet_bl$BID, Centiloid = as.numeric(a4_pet_bl$centiloid))

## Section 2.5: Load A4 MRI Data (WMH as oligodendrocyte proxy)
a4_mri <- fread("A4_VMRI_PRV2_08Jan2026.csv")
a4_mri_bl <- a4_mri[VISCODE %in% c(baseline_codes, 0, 1)]
if (nrow(a4_mri_bl) == 0) a4_mri_bl <- a4_mri[, .SD[1], by = BID]
a4_mri_bl <- a4_mri_bl[!duplicated(BID)]

a4_mri_data <- data.frame(BID = a4_mri_bl$BID)
if (all(c("LeftHippocampus", "RightHippocampus") %in% colnames(a4_mri_bl))) {
  a4_mri_data$Hippocampus <- a4_mri_bl$LeftHippocampus + a4_mri_bl$RightHippocampus
}
if (all(c("LeftWMHypo", "RightWMHypo") %in% colnames(a4_mri_bl))) {
  a4_mri_data$Total_WMH <- a4_mri_bl$LeftWMHypo + a4_mri_bl$RightWMHypo
  a4_mri_data$Log_WMH <- log(a4_mri_data$Total_WMH + 1)
}
if ("IntraCranialVolume" %in% colnames(a4_mri_bl)) {
  a4_mri_data$ICV <- a4_mri_bl$IntraCranialVolume
}

## Section 2.6: Integrate A4 Data
a4_merged <- data.frame(
  BID = a4_subjinfo$BID,
  Age = a4_subjinfo$AGEYR,
  APOE_Genotype = a4_subjinfo$APOEGN,
  APOE4_Carrier = a4_subjinfo$APOE4_Carrier,
  APOE4_Count = a4_subjinfo$APOE4_Count
)

demog_subset <- data.frame(
  BID = a4_demog_bl$BID,
  Gender = a4_demog_bl$PTGENDER,
  Education = a4_demog_bl$PTEDUCAT
)

a4_merged <- merge(a4_merged, demog_subset, by = "BID", all.x = TRUE)
a4_merged <- merge(a4_merged, a4_pacc_data, by = "BID", all.x = TRUE)
a4_merged <- merge(a4_merged, a4_pet_data, by = "BID", all.x = TRUE)
a4_merged <- merge(a4_merged, a4_mri_data, by = "BID", all.x = TRUE)

# Age stratification: <70 vs >=70 per methods
a4_merged$Age_Group <- factor(
  ifelse(a4_merged$Age < 70, "Younger (<70)", "Older (>=70)"),
  levels = c("Younger (<70)", "Older (>=70)")
)
a4_merged$APOE4_Status <- factor(
  ifelse(a4_merged$APOE4_Carrier == 1, "APOE4+", "APOE4-"),
  levels = c("APOE4-", "APOE4+")
)

if ("Log_WMH" %in% colnames(a4_merged)) {
  wmh_median <- median(a4_merged$Log_WMH, na.rm = TRUE)
  a4_merged$WMH_Group <- factor(
    ifelse(a4_merged$Log_WMH > wmh_median, "High WMH", "Low WMH"),
    levels = c("Low WMH", "High WMH")
  )
}

fwrite(a4_merged, "results/A4_Integrated_Data.csv")

# ==============================================================================
# Part 3: HABS Cohort - pTau217 and WMH Mechanism Validation
# ==============================================================================

## Section 3.1: Load HABS Baseline Integrated Data
habs_data <- fread("HABS_Baseline_Integrated_FINAL.csv")

## Section 3.2: Load HABS WMH Data (HD Release 7, QC-passed baseline)
habs_wmh_raw <- read_xlsx("HABS_HD_RP_7_WMH_Combination_Report_Analysis.xlsx", sheet = "x1")
habs_wmh_qc <- habs_wmh_raw %>% filter(tolower(`02_WMH_Fail_QC`) == "pass")
habs_wmh_bl <- habs_wmh_qc %>% distinct(Med_ID, .keep_all = TRUE)

habs_wmh <- habs_wmh_bl %>%
  dplyr::select(Med_ID, Total_WMH = `02_WMH_Volume_Raw`, Log_WMH = `02_WMH_Volume_Log`) %>%
  mutate(
    Med_ID = as.character(Med_ID),
    Total_WMH = as.numeric(as.character(Total_WMH)),
    Log_WMH = as.numeric(as.character(Log_WMH))
  ) %>%
  filter(!is.na(Log_WMH) & is.finite(Log_WMH))

## Section 3.3: Standardize Variable Names
apoe_col <- find_column(habs_data, c("APOE4_Positive", "APOE4_Carrier", "APOE4_Positivity"))
if (!is.null(apoe_col)) habs_data$APOE4_Carrier <- as.integer(habs_data[[apoe_col]])

mmse_col <- find_column(habs_data, c("MMSE_Baseline", "MMSE", "MMSE_Total"))
if (!is.null(mmse_col)) habs_data$MMSE <- as.numeric(habs_data[[mmse_col]])

cdr_col <- find_column(habs_data, c("CDR_Global", "CDR", "CDR_Baseline"))
if (!is.null(cdr_col)) habs_data$CDR <- as.numeric(habs_data[[cdr_col]])

gender_col <- find_column(habs_data, c("Gender", "ID_Gender", "Sex"))
if (!is.null(gender_col)) habs_data$Gender <- habs_data[[gender_col]]

## Section 3.4: Merge WMH Data
merge_key <- find_column(habs_data, c("Med_ID", "MedID", "Subject_ID"))
if (!is.null(merge_key)) {
  habs_data[[merge_key]] <- as.character(habs_data[[merge_key]])
  if (merge_key != "Med_ID") {
    habs_wmh <- habs_wmh %>% dplyr::rename(!!merge_key := Med_ID)
  }
  if ("Log_WMH" %in% colnames(habs_data)) {
    habs_data <- habs_data %>% dplyr::select(-Log_WMH)
  }
  habs_data <- habs_data %>% dplyr::left_join(habs_wmh, by = merge_key)
}

## Section 3.5: Define Cognitive Impairment (CDR >= 0.5 OR MMSE < 24)
habs_data$Cognitive_Impairment <- as.integer(
  (habs_data$CDR >= 0.5 & !is.na(habs_data$CDR)) |
    (habs_data$MMSE < 24 & !is.na(habs_data$MMSE))
)

## Section 3.6: Create Stratification Variables
habs_data$Age_Group <- factor(
  ifelse(habs_data$Age < 70, "Younger (<70)", "Older (>=70)"),
  levels = c("Younger (<70)", "Older (>=70)")
)
habs_data$APOE4_Status <- factor(
  ifelse(habs_data$APOE4_Carrier == 1, "APOE4+", "APOE4-"),
  levels = c("APOE4-", "APOE4+")
)

if ("Log_WMH" %in% colnames(habs_data)) {
  wmh_median <- median(habs_data$Log_WMH, na.rm = TRUE)
  habs_data$WMH_Group <- factor(
    ifelse(habs_data$Log_WMH > wmh_median, "High WMH", "Low WMH"),
    levels = c("Low WMH", "High WMH")
  )
}

habs_data$Gender_num <- as.numeric(as.factor(habs_data$Gender))

fwrite(habs_data, "results/HABS_Integrated_Data.csv")

# ==============================================================================
# Part 4: AIBL Cohort - Data Integration
# ==============================================================================

## Section 4.1: Load AIBL Demographics
baseline_codes_aibl <- c("bl", "BL", "M00", "m00", "0", "sc", "screening")

aibl_demog <- fread("aibl_ptdemog_01-Jun-2018.csv")
aibl_demog_bl <- aibl_demog[VISCODE %in% baseline_codes_aibl]
if (nrow(aibl_demog_bl) == 0) aibl_demog_bl <- aibl_demog[, .SD[1], by = RID]
aibl_demog_bl <- aibl_demog_bl[!duplicated(RID)]

## Section 4.2: Load AIBL APOE Data (from APGEN1 and APGEN2)
aibl_apoe <- fread("aibl_apoeres_01-Jun-2018.csv")
aibl_apoe$APOE4_Carrier <- as.integer(aibl_apoe$APGEN1 == 4 | aibl_apoe$APGEN2 == 4)
aibl_apoe <- aibl_apoe[!duplicated(RID)]

## Section 4.3: Load AIBL Diagnosis Conversion Data (DXCURREN == 3 = AD)
aibl_dx <- fread("aibl_pdxconv_01-Jun-2018.csv")

aibl_dx$Visit_Month <- NA_real_
aibl_dx$Visit_Month[tolower(aibl_dx$VISCODE) == "bl"] <- 0
month_pattern <- grepl("^m[0-9]+$", tolower(aibl_dx$VISCODE))
aibl_dx$Visit_Month[month_pattern] <- as.numeric(gsub("^m", "", tolower(aibl_dx$VISCODE[month_pattern])))

aibl_followup <- aibl_dx[!is.na(Visit_Month), .(
  Max_Followup = max(Visit_Month, na.rm = TRUE),
  Conversion_Month = {
    ad_visits <- Visit_Month[DXCURREN == 3 & !is.na(DXCURREN)]
    if (length(ad_visits) > 0) min(ad_visits) else NA_real_
  },
  Ever_AD = as.integer(any(DXCURREN == 3, na.rm = TRUE))
), by = RID]

aibl_followup$Max_Followup[!is.finite(aibl_followup$Max_Followup)] <- NA

## Section 4.4: Load AIBL MMSE Data
aibl_mmse <- fread("aibl_mmse_01-Jun-2018.csv")
aibl_mmse_bl <- aibl_mmse[VISCODE %in% baseline_codes_aibl]
if (nrow(aibl_mmse_bl) == 0) aibl_mmse_bl <- aibl_mmse[, .SD[1], by = RID]
aibl_mmse_bl <- aibl_mmse_bl[!duplicated(RID)]

## Section 4.5: Integrate AIBL Data
aibl_merged <- data.frame(RID = aibl_demog_bl$RID)

if ("PTGENDER" %in% colnames(aibl_demog_bl)) {
  aibl_merged$Gender <- aibl_demog_bl$PTGENDER
}

if ("PTDOB" %in% colnames(aibl_demog_bl)) {
  raw_dob <- as.character(aibl_demog_bl$PTDOB)
  birth_years <- as.numeric(str_extract(raw_dob, "\\d{4}"))
  if ("EXAMDATE" %in% colnames(aibl_demog_bl)) {
    exam_years <- as.numeric(str_extract(as.character(aibl_demog_bl$EXAMDATE), "\\d{4}"))
    aibl_merged$Age <- exam_years - birth_years
  } else {
    aibl_merged$Age <- 2008 - birth_years
  }
  aibl_merged$Age[aibl_merged$Age < 40 | aibl_merged$Age > 110] <- NA
}

aibl_merged <- merge(aibl_merged, aibl_apoe[, c("RID", "APOE4_Carrier")], by = "RID", all.x = TRUE)
aibl_merged <- merge(aibl_merged, aibl_followup, by = "RID", all.x = TRUE)

# Time to event: conversion month or last follow-up
aibl_merged$Time <- ifelse(
  !is.na(aibl_merged$Conversion_Month) & aibl_merged$Conversion_Month > 0,
  aibl_merged$Conversion_Month,
  aibl_merged$Max_Followup
)
aibl_merged$Event <- aibl_merged$Ever_AD

mmse_subset <- data.frame(RID = aibl_mmse_bl$RID, MMSE = as.numeric(aibl_mmse_bl$MMSCORE))
aibl_merged <- merge(aibl_merged, mmse_subset, by = "RID", all.x = TRUE)

# Age stratification
aibl_merged$Age_Group <- factor(
  ifelse(aibl_merged$Age < 70, "Younger (<70)", "Older (>=70)"),
  levels = c("Younger (<70)", "Older (>=70)")
)
aibl_merged$APOE4_Status <- factor(
  ifelse(aibl_merged$APOE4_Carrier == 1, "APOE4+", "APOE4-"),
  levels = c("APOE4-", "APOE4+")
)

fwrite(aibl_merged, "results/AIBL_Integrated_Data.csv")

# ==============================================================================
# Part 5: Multi-Group SEM Analysis (lavaan, MLR, FIML)
# ==============================================================================

## Section 5.1: A4 Cohort SEM - APOE4 -> Centiloid -> PACC
a4_sem <- subset(a4_merged, !is.na(APOE4_Carrier) & !is.na(Age) & !is.na(Centiloid) & !is.na(PACC))
a4_sem$APOE4 <- as.numeric(a4_sem$APOE4_Carrier)
a4_sem$Age_z <- scale(a4_sem$Age)[,1]
a4_sem$Centiloid_z <- scale(a4_sem$Centiloid)[,1]
a4_sem$PACC_z <- scale(a4_sem$PACC)[,1]
a4_sem$Gender_num <- as.numeric(a4_sem$Gender)

# Overall model
sem_model_a4_overall <- '
  Centiloid_z ~ a * APOE4 + Age_z + Gender_num
  PACC_z ~ b * Centiloid_z + c * APOE4 + Age_z + Gender_num
  indirect := a * b
  direct := c
  total := a * b + c
'

fit_a4_overall <- tryCatch({
  sem(sem_model_a4_overall, data = a4_sem, estimator = "MLR", missing = "FIML")
}, error = function(e) NULL)

# Multi-group model with c(a1, a2) syntax per methods
sem_model_a4_free <- '
  Centiloid_z ~ c(a1, a2) * APOE4 + Age_z + Gender_num
  PACC_z ~ c(b1, b2) * Centiloid_z + c(c1, c2) * APOE4 + Age_z + Gender_num
  indirect1 := a1 * b1
  indirect2 := a2 * b2
  diff_a := a1 - a2
  diff_b := b1 - b2
  diff_indirect := indirect1 - indirect2
'

fit_a4_free <- tryCatch({
  sem(sem_model_a4_free, data = a4_sem, group = "Age_Group", 
      estimator = "MLR", missing = "FIML")
}, error = function(e) NULL)

# Extract model fit indices (CFI, TLI, RMSEA, SRMR)
if (!is.null(fit_a4_free)) {
  fit_a4_indices <- fitMeasures(fit_a4_free, c("cfi", "tli", "rmsea", "srmr"))
  params_a4 <- parameterEstimates(fit_a4_free, standardized = TRUE)
  
  # Bootstrap CI for indirect effects (1000 iterations)
  fit_a4_boot <- tryCatch({
    sem(sem_model_a4_free, data = a4_sem, group = "Age_Group",
        estimator = "MLR", missing = "FIML", se = "bootstrap", bootstrap = 1000)
  }, error = function(e) NULL)
}

## Section 5.2: HABS Cohort SEM - APOE4 -> pTau217 -> MMSE
if ("pTau217" %in% colnames(habs_data) && "MMSE" %in% colnames(habs_data)) {
  
  habs_sem <- subset(habs_data, !is.na(APOE4_Carrier) & !is.na(Age) & 
                       !is.na(pTau217) & !is.na(MMSE))
  habs_sem$APOE4 <- as.numeric(habs_sem$APOE4_Carrier)
  habs_sem$Age_z <- scale(habs_sem$Age)[,1]
  habs_sem$pTau_z <- scale(habs_sem$pTau217)[,1]
  habs_sem$MMSE_z <- scale(habs_sem$MMSE)[,1]
  
  # Overall model
  sem_model_habs_overall <- '
    pTau_z ~ a * APOE4 + Age_z + Gender_num
    MMSE_z ~ b * pTau_z + c * APOE4 + Age_z + Gender_num
    indirect := a * b
  '
  
  fit_habs_overall <- tryCatch({
    sem(sem_model_habs_overall, data = habs_sem, estimator = "MLR", missing = "FIML")
  }, error = function(e) NULL)
  
  # Multi-group model
  sem_model_habs_free <- '
    pTau_z ~ c(a1, a2) * APOE4 + Age_z + Gender_num
    MMSE_z ~ c(b1, b2) * pTau_z + c(c1, c2) * APOE4 + Age_z + Gender_num
    indirect1 := a1 * b1
    indirect2 := a2 * b2
    diff_a := a1 - a2
    diff_indirect := indirect1 - indirect2
  '
  
  fit_habs_free <- tryCatch({
    sem(sem_model_habs_free, data = habs_sem, group = "Age_Group",
        estimator = "MLR", missing = "FIML")
  }, error = function(e) NULL)
  
  if (!is.null(fit_habs_free)) {
    fit_habs_indices <- fitMeasures(fit_habs_free, c("cfi", "tli", "rmsea", "srmr"))
    params_habs <- parameterEstimates(fit_habs_free, standardized = TRUE)
  }
}

## Section 5.3: Save SEM Results
sem_results <- list(
  A4_Overall = if(exists("fit_a4_overall")) fit_a4_overall else NULL,
  A4_Free = if(exists("fit_a4_free")) fit_a4_free else NULL,
  A4_Fit_Indices = if(exists("fit_a4_indices")) fit_a4_indices else NULL,
  HABS_Overall = if(exists("fit_habs_overall")) fit_habs_overall else NULL,
  HABS_Free = if(exists("fit_habs_free")) fit_habs_free else NULL,
  HABS_Fit_Indices = if(exists("fit_habs_indices")) fit_habs_indices else NULL
)

saveRDS(sem_results, "results/SEM_Results.rds")

# Export SEM fit indices table
if (exists("fit_a4_indices") || exists("fit_habs_indices")) {
  fit_table <- data.frame(
    Cohort = c("A4", "HABS"),
    CFI = c(if(exists("fit_a4_indices")) fit_a4_indices["cfi"] else NA,
            if(exists("fit_habs_indices")) fit_habs_indices["cfi"] else NA),
    TLI = c(if(exists("fit_a4_indices")) fit_a4_indices["tli"] else NA,
            if(exists("fit_habs_indices")) fit_habs_indices["tli"] else NA),
    RMSEA = c(if(exists("fit_a4_indices")) fit_a4_indices["rmsea"] else NA,
              if(exists("fit_habs_indices")) fit_habs_indices["rmsea"] else NA),
    SRMR = c(if(exists("fit_a4_indices")) fit_a4_indices["srmr"] else NA,
             if(exists("fit_habs_indices")) fit_habs_indices["srmr"] else NA)
  )
  fwrite(fit_table, "results/tables/Table_SEM_Fit_Indices.csv")
}

# ==============================================================================
# Part 6: Clinical Utility Assessment (Firth, AUC, AUPRC, NRI, IDI, DCA)
# ==============================================================================

## Section 6.1: HABS Clinical Utility Analysis
## Outcome: Cognitive impairment (CDR >= 0.5 OR MMSE < 24)
if ("pTau217" %in% colnames(habs_data) && "Log_WMH" %in% colnames(habs_data) &&
    "Cognitive_Impairment" %in% colnames(habs_data)) {
  
  habs_clinical <- subset(habs_data,
                          !is.na(Age) & !is.na(Gender_num) & !is.na(APOE4_Carrier) &
                            !is.na(pTau217) & !is.na(Log_WMH) & !is.na(Cognitive_Impairment)
  )
  
  if (nrow(habs_clinical) >= 100) {
    
    ## Section 6.2: Train/Test Split (70/30 with fixed seed)
    set.seed(42)
    train_idx <- sample(1:nrow(habs_clinical), size = floor(0.7 * nrow(habs_clinical)))
    habs_train <- habs_clinical[train_idx, ]
    habs_test <- habs_clinical[-train_idx, ]
    
    ## Section 6.3: Firth-Corrected Logistic Regression
    model_base <- logistf(Cognitive_Impairment ~ Age + Gender_num + APOE4_Carrier + pTau217,
                          data = habs_train)
    model_complete <- logistf(Cognitive_Impairment ~ Age + Gender_num + APOE4_Carrier + 
                                pTau217 + Log_WMH, data = habs_train)
    
    pred_base <- predict(model_base, newdata = habs_test, type = "response")
    pred_complete <- predict(model_complete, newdata = habs_test, type = "response")
    
    ## Section 6.4: ROC/AUC Analysis with DeLong CI
    roc_base <- roc(habs_test$Cognitive_Impairment, pred_base, quiet = TRUE)
    roc_complete <- roc(habs_test$Cognitive_Impairment, pred_complete, quiet = TRUE)
    
    auc_base <- auc(roc_base)
    auc_complete <- auc(roc_complete)
    ci_base <- ci.auc(roc_base)
    ci_complete <- ci.auc(roc_complete)
    
    roc_test <- roc.test(roc_base, roc_complete, method = "delong")
    
    ## Section 6.5: AUPRC (Precision-Recall Curve) per methods
    pr_base <- pr.curve(scores.class0 = pred_base[habs_test$Cognitive_Impairment == 1],
                        scores.class1 = pred_base[habs_test$Cognitive_Impairment == 0],
                        curve = TRUE)
    pr_complete <- pr.curve(scores.class0 = pred_complete[habs_test$Cognitive_Impairment == 1],
                            scores.class1 = pred_complete[habs_test$Cognitive_Impairment == 0],
                            curve = TRUE)
    
    auprc_base <- pr_base$auc.integral
    auprc_complete <- pr_complete$auc.integral
    
    ## Section 6.6: NRI Analysis with Three Threshold Sets per methods
    calc_nri <- function(pred_old, pred_new, outcome, thresholds) {
      cat_old <- cut(pred_old, breaks = thresholds, include.lowest = TRUE)
      cat_new <- cut(pred_new, breaks = thresholds, include.lowest = TRUE)
      
      events_idx <- which(outcome == 1)
      nonevents_idx <- which(outcome == 0)
      
      if (length(events_idx) == 0 || length(nonevents_idx) == 0) {
        return(list(events_nri = NA, nonevents_nri = NA, total_nri = NA))
      }
      
      events_table <- table(Old = cat_old[events_idx], New = cat_new[events_idx])
      nonevents_table <- table(Old = cat_old[nonevents_idx], New = cat_new[nonevents_idx])
      
      events_up <- sum(events_table[lower.tri(events_table)])
      events_down <- sum(events_table[upper.tri(events_table)])
      events_nri <- (events_up - events_down) / length(events_idx)
      
      nonevents_down <- sum(nonevents_table[upper.tri(nonevents_table)])
      nonevents_up <- sum(nonevents_table[lower.tri(nonevents_table)])
      nonevents_nri <- (nonevents_down - nonevents_up) / length(nonevents_idx)
      
      return(list(events_nri = events_nri, nonevents_nri = nonevents_nri,
                  total_nri = events_nri + nonevents_nri))
    }
    
    # Three threshold sets per methods
    thresholds_high <- c(0, 0.20, 0.40, 0.60, 1)      # High event rate
    thresholds_standard <- c(0, 0.10, 0.20, 1)        # Standard
    thresholds_fine <- c(0, 0.15, 0.25, 0.35, 0.50, 1) # Fine-grained
    
    nri_high <- calc_nri(pred_base, pred_complete, habs_test$Cognitive_Impairment, thresholds_high)
    nri_standard <- calc_nri(pred_base, pred_complete, habs_test$Cognitive_Impairment, thresholds_standard)
    nri_fine <- calc_nri(pred_base, pred_complete, habs_test$Cognitive_Impairment, thresholds_fine)
    
    ## Section 6.7: IDI Analysis with Bootstrap (1000 iterations)
    events_idx <- habs_test$Cognitive_Impairment == 1
    nonevents_idx <- habs_test$Cognitive_Impairment == 0
    
    events_improvement <- mean(pred_complete[events_idx]) - mean(pred_base[events_idx])
    nonevents_improvement <- mean(pred_complete[nonevents_idx]) - mean(pred_base[nonevents_idx])
    idi_value <- events_improvement - nonevents_improvement
    
    set.seed(42)
    n_bootstrap <- 1000
    idi_bootstrap <- numeric(n_bootstrap)
    
    for (i in 1:n_bootstrap) {
      boot_idx <- sample(1:nrow(habs_test), replace = TRUE)
      boot_data <- habs_test[boot_idx, ]
      boot_events_idx <- boot_data$Cognitive_Impairment == 1
      boot_nonevents_idx <- boot_data$Cognitive_Impairment == 0
      
      if (sum(boot_events_idx) > 0 && sum(boot_nonevents_idx) > 0) {
        boot_events_imp <- mean(pred_complete[boot_idx][boot_events_idx]) -
          mean(pred_base[boot_idx][boot_events_idx])
        boot_nonevents_imp <- mean(pred_complete[boot_idx][boot_nonevents_idx]) -
          mean(pred_base[boot_idx][boot_nonevents_idx])
        idi_bootstrap[i] <- boot_events_imp - boot_nonevents_imp
      } else {
        idi_bootstrap[i] <- NA
      }
    }
    
    idi_bootstrap <- idi_bootstrap[!is.na(idi_bootstrap)]
    idi_ci <- quantile(idi_bootstrap, c(0.025, 0.975))
    idi_p_value <- 2 * min(mean(idi_bootstrap > 0), mean(idi_bootstrap < 0))
    
    ## Section 6.8: DCA Analysis (threshold 0.01 to 0.80)
    calc_net_benefit <- function(pred, outcome, threshold) {
      if (threshold == 0) return(mean(outcome))
      if (threshold >= 1) return(0)
      tp <- sum(pred >= threshold & outcome == 1)
      fp <- sum(pred >= threshold & outcome == 0)
      n <- length(outcome)
      net_benefit <- (tp / n) - (fp / n) * (threshold / (1 - threshold))
      return(net_benefit)
    }
    
    thresholds_dca <- seq(0.01, 0.80, by = 0.01)
    nb_base <- sapply(thresholds_dca, function(t) 
      calc_net_benefit(pred_base, habs_test$Cognitive_Impairment, t))
    nb_complete <- sapply(thresholds_dca, function(t) 
      calc_net_benefit(pred_complete, habs_test$Cognitive_Impairment, t))
    
    max_benefit_idx <- which.max(nb_complete)
    optimal_threshold <- thresholds_dca[max_benefit_idx]
    
    ## Section 6.9: Save Clinical Utility Results
    clinical_results <- list(
      AUC = list(Base = auc_base, Complete = auc_complete, 
                 Delta = auc_complete - auc_base, P_value = roc_test$p.value,
                 CI_Base = ci_base, CI_Complete = ci_complete),
      AUPRC = list(Base = auprc_base, Complete = auprc_complete,
                   Delta = auprc_complete - auprc_base),
      NRI = list(High = nri_high, Standard = nri_standard, Fine = nri_fine),
      IDI = list(Value = idi_value, CI = idi_ci, P_value = idi_p_value),
      DCA = list(Optimal_Threshold = optimal_threshold)
    )
    
    saveRDS(clinical_results, "results/Clinical_Utility_Results.rds")
    
    # Export clinical utility table
    clinical_table <- data.frame(
      Metric = c("AUC (Base)", "AUC (Complete)", "ΔAUC", "DeLong P",
                 "AUPRC (Base)", "AUPRC (Complete)", "ΔAUPRC",
                 "NRI (High)", "NRI (Standard)", "NRI (Fine)",
                 "IDI", "IDI 95% CI"),
      Value = c(
        sprintf("%.3f (%.3f-%.3f)", auc_base, ci_base[1], ci_base[3]),
        sprintf("%.3f (%.3f-%.3f)", auc_complete, ci_complete[1], ci_complete[3]),
        sprintf("%.3f", auc_complete - auc_base),
        format(roc_test$p.value, scientific = TRUE, digits = 3),
        sprintf("%.3f", auprc_base),
        sprintf("%.3f", auprc_complete),
        sprintf("%.3f", auprc_complete - auprc_base),
        sprintf("%.3f", nri_high$total_nri),
        sprintf("%.3f", nri_standard$total_nri),
        sprintf("%.3f", nri_fine$total_nri),
        sprintf("%.4f", idi_value),
        sprintf("%.4f-%.4f", idi_ci[1], idi_ci[2])
      )
    )
    fwrite(clinical_table, "results/tables/Table_Clinical_Utility.csv")
  }
}

# ==============================================================================
# Part 7: AIBL Cox Regression and AFT Models
# ==============================================================================

## Section 7.1: Prepare Survival Data
aibl_surv <- subset(aibl_merged,
                    !is.na(Time) & is.finite(Time) & Time > 0 & 
                      !is.na(APOE4_Carrier) & !is.na(Age) & !is.na(Event)
)
aibl_surv$Gender_num <- as.numeric(as.factor(aibl_surv$Gender))

if (nrow(aibl_surv) >= 50 && sum(aibl_surv$Event, na.rm = TRUE) >= 10) {
  
  ## Section 7.2: Cox Proportional Hazards Regression
  cox_base <- coxph(Surv(Time, Event) ~ APOE4_Carrier + Age + Gender_num, 
                    data = aibl_surv)
  
  # Cox with age interaction
  cox_interaction <- coxph(Surv(Time, Event) ~ APOE4_Carrier * Age + Gender_num, 
                           data = aibl_surv)
  
  # Extract HR and 95% CI
  hr_apoe4 <- exp(coef(cox_base)["APOE4_Carrier"])
  hr_ci <- exp(confint(cox_base)["APOE4_Carrier", ])
  
  # Test proportional hazards assumption
  ph_test <- cox.zph(cox_base)
  
  ## Section 7.3: Accelerated Failure Time Model (Weibull distribution)
  aft_model <- survreg(Surv(Time, Event) ~ APOE4_Carrier + Age + Gender_num,
                       data = aibl_surv, dist = "weibull")
  
  # Acceleration factor (time ratio)
  af_apoe4 <- exp(-coef(aft_model)["APOE4_Carrier"])
  
  ## Section 7.4: Age-Stratified Cox Regression
  cox_younger <- tryCatch({
    coxph(Surv(Time, Event) ~ APOE4_Carrier + Age + Gender_num,
          data = aibl_surv[aibl_surv$Age_Group == "Younger (<70)", ])
  }, error = function(e) NULL)
  
  cox_older <- tryCatch({
    coxph(Surv(Time, Event) ~ APOE4_Carrier + Age + Gender_num,
          data = aibl_surv[aibl_surv$Age_Group == "Older (>=70)", ])
  }, error = function(e) NULL)
  
  ## Section 7.5: Save Survival Results
  survival_results <- list(
    Cox_Base = cox_base,
    Cox_Interaction = cox_interaction,
    HR_APOE4 = hr_apoe4,
    HR_CI = hr_ci,
    PH_Test = ph_test,
    AFT = aft_model,
    AF_APOE4 = af_apoe4,
    Cox_Younger = cox_younger,
    Cox_Older = cox_older
  )
  
  saveRDS(survival_results, "results/Survival_Results.rds")
  
  # Export survival results table
  surv_table <- data.frame(
    Model = c("Cox Overall", "AFT Weibull"),
    N = rep(nrow(aibl_surv), 2),
    N_Events = rep(sum(aibl_surv$Event), 2),
    APOE4_Effect = c(
      sprintf("HR=%.2f (%.2f-%.2f)", hr_apoe4, hr_ci[1], hr_ci[2]),
      sprintf("AF=%.2f", af_apoe4)
    ),
    P_value = c(
      format(summary(cox_base)$coefficients["APOE4_Carrier", "Pr(>|z|)"], 
             scientific = TRUE, digits = 3),
      format(summary(aft_model)$table["APOE4_Carrier", "p"], 
             scientific = TRUE, digits = 3)
    )
  )
  fwrite(surv_table, "results/tables/Table_Survival_Results.csv")
}

# ==============================================================================
# Part 8: Visualization
# ==============================================================================

## Section 8.1: A4 PACC vs Age by APOE4 Status
if (exists("a4_merged") && "PACC" %in% colnames(a4_merged)) {
  plot_data <- a4_merged[!is.na(a4_merged$PACC) & !is.na(a4_merged$Age) & 
                           !is.na(a4_merged$APOE4_Status), ]
  
  if (nrow(plot_data) > 50) {
    p_a4_pacc <- ggplot(plot_data, aes(x = Age, y = PACC, color = APOE4_Status)) +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
      scale_color_manual(values = c("APOE4-" = "#2166AC", "APOE4+" = "#B2182B")) +
      labs(title = "A4: PACC vs Age by APOE4 Status",
           x = "Age (years)", y = "PACC Score", color = "APOE4 Status") +
      theme_nc
    
    ggsave("results/figures/Figure_A4_PACC_Age_APOE4.pdf", p_a4_pacc, 
           width = 8, height = 6, dpi = 300)
  }
}

## Section 8.2: A4 Centiloid vs Age by APOE4 Status
if (exists("a4_merged") && "Centiloid" %in% colnames(a4_merged)) {
  plot_data <- a4_merged[!is.na(a4_merged$Centiloid) & !is.na(a4_merged$Age) & 
                           !is.na(a4_merged$APOE4_Status), ]
  
  if (nrow(plot_data) > 50) {
    p_a4_centiloid <- ggplot(plot_data, aes(x = Age, y = Centiloid, color = APOE4_Status)) +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
      scale_color_manual(values = c("APOE4-" = "#2166AC", "APOE4+" = "#B2182B")) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "gray50") +
      labs(title = "A4: Amyloid Burden vs Age by APOE4 Status",
           x = "Age (years)", y = "Centiloid", color = "APOE4 Status") +
      theme_nc
    
    ggsave("results/figures/Figure_A4_Centiloid_Age_APOE4.pdf", p_a4_centiloid, 
           width = 8, height = 6, dpi = 300)
  }
}

## Section 8.3: HABS pTau217 vs Age by APOE4 Status
if (exists("habs_data") && "pTau217" %in% colnames(habs_data)) {
  plot_data <- habs_data[!is.na(habs_data$pTau217) & !is.na(habs_data$Age) & 
                           !is.na(habs_data$APOE4_Status), ]
  
  if (nrow(plot_data) > 30) {
    p_habs_ptau <- ggplot(plot_data, aes(x = Age, y = pTau217, color = APOE4_Status)) +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
      scale_color_manual(values = c("APOE4-" = "#2166AC", "APOE4+" = "#B2182B")) +
      labs(title = "HABS: Plasma pTau217 vs Age by APOE4 Status",
           x = "Age (years)", y = "Plasma pTau217 (pg/mL)", color = "APOE4 Status") +
      theme_nc
    
    ggsave("results/figures/Figure_HABS_pTau_Age_APOE4.pdf", p_habs_ptau, 
           width = 8, height = 6, dpi = 300)
  }
}

## Section 8.4: AIBL Kaplan-Meier Survival Curves
if (exists("aibl_surv") && nrow(aibl_surv) > 50) {
  km_apoe <- survfit(Surv(Time, Event) ~ APOE4_Status, data = aibl_surv)
  
  p_km <- ggsurvplot(
    km_apoe, data = aibl_surv,
    pval = TRUE, pval.method = TRUE,
    conf.int = TRUE, risk.table = TRUE,
    palette = c("#2166AC", "#B2182B"),
    xlab = "Time (months)", ylab = "Conversion-free Probability",
    title = "AIBL: AD Conversion by APOE4 Status",
    legend.title = "APOE4", legend.labs = c("APOE4-", "APOE4+"),
    ggtheme = theme_nc
  )
  
  ggsave("results/figures/Figure_AIBL_KM_Curve.pdf", p_km$plot, 
         width = 10, height = 8, dpi = 300)
}

## Section 8.5: ROC Curve Comparison (HABS)
if (exists("roc_base") && exists("roc_complete")) {
  p_roc <- ggplot() +
    geom_line(aes(x = 1 - roc_base$specificities, y = roc_base$sensitivities,
                  color = "Base Model (pTau217)"), linewidth = 1.2) +
    geom_line(aes(x = 1 - roc_complete$specificities, y = roc_complete$sensitivities,
                  color = "Complete Model (+WMH)"), linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("Base Model (pTau217)" = "#E69F00",
                                  "Complete Model (+WMH)" = "#009E73")) +
    labs(title = "HABS: ROC Curve Comparison",
         subtitle = sprintf("Base AUC=%.3f, Complete AUC=%.3f, ΔAUC=%.3f (p=%.4e)",
                            auc_base, auc_complete, auc_complete - auc_base, roc_test$p.value),
         x = "1 - Specificity", y = "Sensitivity") +
    theme_nc +
    theme(legend.position = c(0.7, 0.2), legend.title = element_blank())
  
  ggsave("results/figures/Figure_HABS_ROC_Comparison.pdf", p_roc, 
         width = 8, height = 6, dpi = 300)
}

## Section 8.6: DCA Decision Curve
if (exists("thresholds_dca") && exists("nb_base") && exists("nb_complete")) {
  nb_all <- sapply(thresholds_dca, function(t) {
    event_rate <- mean(habs_test$Cognitive_Impairment)
    if (t >= 1) return(0)
    event_rate - (1 - event_rate) * (t / (1 - t))
  })
  nb_none <- rep(0, length(thresholds_dca))
  
  dca_data <- data.frame(
    Threshold = rep(thresholds_dca, 4),
    Net_Benefit = c(nb_base, nb_complete, nb_all, nb_none),
    Strategy = rep(c("Base Model", "Complete Model", "Treat All", "Treat None"),
                   each = length(thresholds_dca))
  )
  
  p_dca <- ggplot(dca_data, aes(x = Threshold, y = Net_Benefit, color = Strategy)) +
    geom_line(linewidth = 1.2) +
    geom_vline(xintercept = optimal_threshold, linetype = "dotted", color = "red") +
    scale_color_manual(values = c("Base Model" = "#E69F00", "Complete Model" = "#009E73",
                                  "Treat All" = "#999999", "Treat None" = "#000000")) +
    labs(title = "Decision Curve Analysis",
         subtitle = sprintf("Optimal threshold: %.3f", optimal_threshold),
         x = "Risk Threshold", y = "Net Benefit") +
    theme_nc +
    theme(legend.position = c(0.7, 0.8))
  
  ggsave("results/figures/Figure_HABS_DCA.pdf", p_dca, width = 10, height = 6, dpi = 300)
}

# ==============================================================================
# Part 9: Results Summary
# ==============================================================================

## Section 9.1: Sample Size Summary
sample_summary <- data.frame(
  Cohort = c("A4", "HABS", "AIBL"),
  N_Total = c(
    if(exists("a4_merged")) nrow(a4_merged) else NA,
    if(exists("habs_data")) nrow(habs_data) else NA,
    if(exists("aibl_merged")) nrow(aibl_merged) else NA
  ),
  N_APOE4_Pos = c(
    if(exists("a4_merged")) sum(a4_merged$APOE4_Carrier == 1, na.rm = TRUE) else NA,
    if(exists("habs_data")) sum(habs_data$APOE4_Carrier == 1, na.rm = TRUE) else NA,
    if(exists("aibl_merged")) sum(aibl_merged$APOE4_Carrier == 1, na.rm = TRUE) else NA
  ),
  N_Younger = c(
    if(exists("a4_merged")) sum(a4_merged$Age_Group == "Younger (<70)", na.rm = TRUE) else NA,
    if(exists("habs_data")) sum(habs_data$Age_Group == "Younger (<70)", na.rm = TRUE) else NA,
    if(exists("aibl_merged")) sum(aibl_merged$Age_Group == "Younger (<70)", na.rm = TRUE) else NA
  ),
  N_Older = c(
    if(exists("a4_merged")) sum(a4_merged$Age_Group == "Older (>=70)", na.rm = TRUE) else NA,
    if(exists("habs_data")) sum(habs_data$Age_Group == "Older (>=70)", na.rm = TRUE) else NA,
    if(exists("aibl_merged")) sum(aibl_merged$Age_Group == "Older (>=70)", na.rm = TRUE) else NA
  ),
  Primary_Analysis = c("Multi-group SEM", "Clinical Utility", "Survival Analysis")
)

fwrite(sample_summary, "results/tables/Table_Sample_Summary.csv")

## Section 9.2: Save Complete Results
all_results <- list(
  a4_data = if(exists("a4_merged")) a4_merged else NULL,
  habs_data = if(exists("habs_data")) habs_data else NULL,
  aibl_data = if(exists("aibl_merged")) aibl_merged else NULL,
  sem_results = if(exists("sem_results")) sem_results else NULL,
  clinical_results = if(exists("clinical_results")) clinical_results else NULL,
  survival_results = if(exists("survival_results")) survival_results else NULL,
  analysis_date = Sys.time(),
  r_version = R.version.string
)

saveRDS(all_results, "results/Complete_Analysis_Results.rds")

## Section 9.3: Print Summary
cat("\n")
cat("================================================================================\n")
cat("  A4/HABS/AIBL Validation Analysis Complete\n")
cat("================================================================================\n")
cat("\nOutput files:\n")
cat("  Data:\n")
cat("    - results/A4_Integrated_Data.csv\n")
cat("    - results/HABS_Integrated_Data.csv\n")
cat("    - results/AIBL_Integrated_Data.csv\n")
cat("  Tables:\n")
cat("    - results/tables/Table_Sample_Summary.csv\n")
cat("    - results/tables/Table_SEM_Fit_Indices.csv\n")
cat("    - results/tables/Table_Clinical_Utility.csv\n")
cat("    - results/tables/Table_Survival_Results.csv\n")
cat("  Figures:\n")
cat("    - results/figures/Figure_A4_PACC_Age_APOE4.pdf\n")
cat("    - results/figures/Figure_A4_Centiloid_Age_APOE4.pdf\n")
cat("    - results/figures/Figure_HABS_pTau_Age_APOE4.pdf\n")
cat("    - results/figures/Figure_AIBL_KM_Curve.pdf\n")
cat("    - results/figures/Figure_HABS_ROC_Comparison.pdf\n")
cat("    - results/figures/Figure_HABS_DCA.pdf\n")
cat("\n")

sessionInfo()


