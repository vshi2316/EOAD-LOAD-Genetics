# ==============================================================================
# Cross-Cohort Mechanistic Validation: A4, HABS, and AIBL
# ==============================================================================
#
# Purpose: Validate WMH-mediated neuroinflammatory mechanisms across three
#          independent cohorts (A4, HABS, AIBL) for Alzheimer's disease
#
# Analyses:
#   1.  Data loading for A4, HABS, and AIBL cohorts
#   2.  A4: WMH-cognition association (WMH -> PACC, HC3 robust SE, BP test, Q-Q)
#   3.  HABS: WMH-pTau217 association (HC3 robust SE, BP test, Q-Q)
#   4.  HABS: Mediation analysis (WMH -> pTau217 -> MMSE, 5000 bootstrap, z-scored)
#   5.  HABS: SEM confirmation (5000 bootstrap, z-scored variables)
#   6.  HABS: Age-stratified mediation AND SEM (<75 vs >=75)
#   7.  HABS: Age interaction tests (Age x WMH, Age x pTau217)
#   8.  HABS: Clinical utility (Firth, train/test, AUC, PR-AUC, NRI/IDI bootstrap, DCA 0.01-0.80)
#   9.  AIBL: Survival analysis (Cox PH, KM curves, AFT models, time in months)
#   10. Missing data and selection bias assessment (numeric + categorical SMD)
#   11. Cross-cohort forest plot visualization
#   12. Mediation path diagram and comprehensive visualizations
#   13. Main execution script
#
# Data requirements:
#   - A4: SUBJINFO, PTDEMOG, SPPACC, PETSUVR, VMRI (5 separate files)
#     Column mapping: APOEGN -> APOE4_Carrier, PTGENDER (1=M,2=F),
#     PTEDUCAT, PACCRN, centiloid, LeftWMHypo + RightWMHypo -> Total_WMH
#   - HABS: Genomics xlsx, Biomarkers csv, WMH xlsx, Clinical csv (4 files)
#     Column mapping: APOE4_Positivity, ID_Gender, ID_Education,
#     r7_QTX_Plasma_pTau217_V2/PLUS/standard (coalesce),
#     02_WMH_Fail_QC (pass), 02_WMH_Volume_Raw, 02_WMH_Volume_Log,
#     MMSE_Total, CDR_Global
#   - AIBL: ptdemog, apoeres, pdxconv (3 csv files)
#     Column mapping: PTGENDER (1=M,2=F), PTDOB (birth year),
#     APGEN1/APGEN2 -> APOE4, VISCODE -> months, DXCURREN==3 -> AD
#
# Statistical notes:
#   - All continuous variables z-scored for mediation/SEM (per manuscript)
#   - WMH: log(raw+1) then z-score; effect = per 1 SD of log(WMH)
#   - pTau217: log(raw) then z-score for mediation/SEM
#   - Bootstrap: 5000 iterations for mediation and SEM
#   - AIBL time unit: months (derived from VISCODE m06, m18, etc.)
#   - Clinical utility: 70/30 train/test split, seed=2026
#
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(sandwich)
  library(lmtest)
  library(mediation)
  library(lavaan)
  library(survival)
  library(survminer)
  library(cowplot)
  library(logistf)
  library(pROC)
  library(PRROC)  # for PR-AUC
})


# ==============================================================================
# 1. DATA LOADING
# ==============================================================================

# --- 1.1 A4 Cohort ---
# Reads 5 separate files: SUBJINFO, PTDEMOG, SPPACC, PETSUVR, VMRI
# Uses VISCODE baseline filtering and specific column names directly
load_a4_data <- function(subjinfo_file, ptdemog_file, sppacc_file,
                         petsuvr_file, vmri_file) {

  baseline_codes <- c("bl", "BL", "sc", "SC", "screening", "baseline", "0", "1")

  # 1) APOE from SUBJINFO: APOEGN column, grepl("4") for carrier status
  subjinfo <- fread(subjinfo_file)
  subjinfo$APOE4_Carrier <- as.integer(grepl("4", subjinfo$APOEGN))
  subjinfo <- subjinfo[!duplicated(BID)]

  # 2) Demographics from PTDEMOG (baseline only)
  demog <- fread(ptdemog_file)
  demog_bl <- demog[VISCODE %in% baseline_codes]
  if (nrow(demog_bl) == 0) demog_bl <- demog[, .SD[1], by = BID]
  demog_bl <- demog_bl[!duplicated(BID)]

  # 3) PACC from SPPACC (PACCTS at baseline)
  pacc <- fread(sppacc_file)
  pacc_total <- pacc[PACCQSNUM == "PACCTS"]
  pacc_bl <- pacc_total[VISCODE %in% baseline_codes]
  if (nrow(pacc_bl) == 0) pacc_bl <- pacc_total[, .SD[1], by = BID]
  pacc_bl <- pacc_bl[!duplicated(BID)]

  # 4) Centiloid from PETSUVR (Composite_Summary at baseline)
  pet <- fread(petsuvr_file)
  pet_composite <- pet[brain_region == "Composite_Summary"]
  pet_bl <- pet_composite[VISCODE %in% baseline_codes]
  if (nrow(pet_bl) == 0) pet_bl <- pet_composite[, .SD[1], by = BID]
  pet_bl <- pet_bl[!duplicated(BID)]

  # 5) WMH from VMRI (baseline): LeftWMHypo + RightWMHypo
  mri <- fread(vmri_file)
  mri_bl <- mri[VISCODE %in% c(baseline_codes, 0, 1)]
  if (nrow(mri_bl) == 0) mri_bl <- mri[, .SD[1], by = BID]
  mri_bl <- mri_bl[!duplicated(BID)]

  # --- Merge all on BID ---
  merged <- data.frame(
    BID = subjinfo$BID,
    Age = subjinfo$AGEYR,
    APOE4_Carrier = subjinfo$APOE4_Carrier
  )

  # Demographics: PTGENDER (1=Male, 2=Female), PTEDUCAT
  demog_sub <- data.frame(
    BID = demog_bl$BID,
    Gender_Raw = demog_bl$PTGENDER,
    Education = demog_bl$PTEDUCAT
  )
  merged <- merge(merged, demog_sub, by = "BID", all.x = TRUE)

  # Explicit gender recoding (A4: 1=Male, 2=Female)
  merged$Gender <- factor(
    ifelse(merged$Gender_Raw == 1, "Male",
           ifelse(merged$Gender_Raw == 2, "Female", NA)),
    levels = c("Male", "Female")
  )

  # PACC: PACCRN column
  pacc_sub <- data.frame(BID = pacc_bl$BID, PACC = as.numeric(pacc_bl$PACCRN))
  merged <- merge(merged, pacc_sub, by = "BID", all.x = TRUE)

  # Centiloid
  pet_sub <- data.frame(BID = pet_bl$BID, Centiloid = as.numeric(pet_bl$centiloid))
  merged <- merge(merged, pet_sub, by = "BID", all.x = TRUE)

  # WMH: LeftWMHypo + RightWMHypo -> Total_WMH -> log(+1) -> z-score
  if (all(c("LeftWMHypo", "RightWMHypo") %in% colnames(mri_bl))) {
    mri_sub <- data.frame(
      BID = mri_bl$BID,
      Total_WMH = mri_bl$LeftWMHypo + mri_bl$RightWMHypo,
      Log_WMH_Raw = log(mri_bl$LeftWMHypo + mri_bl$RightWMHypo + 1)
    )
    merged <- merge(merged, mri_sub, by = "BID", all.x = TRUE)
  } else {
    warning("A4: LeftWMHypo/RightWMHypo columns not found in VMRI file")
    merged$Total_WMH <- NA
    merged$Log_WMH_Raw <- NA
  }

  # WMH standardization (per 1 SD of log-transformed WMH)
  wmh_sd <- sd(merged$Log_WMH_Raw, na.rm = TRUE)
  merged$Log_WMH_z <- as.numeric(scale(merged$Log_WMH_Raw))

  # Save original SD for effect interpretation
  cat(sprintf("  A4: N=%d, WMH SD(log)=%.3f\n", nrow(merged), wmh_sd))
  cat(sprintf("  Gender: %s\n",
              paste(names(table(merged$Gender)), table(merged$Gender),
                    sep = "=", collapse = ", ")))

  attr(merged, "wmh_sd") <- wmh_sd
  return(merged)
}


# --- 1.2 HABS Cohort ---
# Reads 4 files: Genomics xlsx, Biomarkers csv, WMH xlsx, Clinical csv
# Uses multi-source pTau217 coalesce strategy and WMH QC filtering
load_habs_data <- function(genomics_xlsx, biomarkers_csv, wmh_xlsx,
                           clinical_csv, wmh_sheet = "x1") {

  # 1) Demographics + APOE from Genomics xlsx
  genomics <- read_xlsx(genomics_xlsx)
  apoe <- genomics %>%
    filter(APOE4_Positivity >= 0) %>%
    dplyr::select(Med_ID, Age, Gender_Raw = ID_Gender,
                  Education = ID_Education,
                  APOE4_Carrier = APOE4_Positivity) %>%
    distinct(Med_ID, .keep_all = TRUE)

  # Explicit gender recoding for HABS
  if (is.character(apoe$Gender_Raw)) {
    if (all(apoe$Gender_Raw %in% c("Male", "Female", NA))) {
      apoe$Gender <- factor(apoe$Gender_Raw, levels = c("Male", "Female"))
    } else if (all(apoe$Gender_Raw %in% c("M", "F", NA))) {
      apoe$Gender <- factor(
        ifelse(apoe$Gender_Raw == "M", "Male",
               ifelse(apoe$Gender_Raw == "F", "Female", NA)),
        levels = c("Male", "Female"))
    } else {
      apoe$Gender <- factor(apoe$Gender_Raw)
    }
  } else {
    apoe$Gender <- factor(
      ifelse(apoe$Gender_Raw == 1, "Male",
             ifelse(apoe$Gender_Raw == 0, "Female", NA)),
      levels = c("Male", "Female"))
  }

  # 2) pTau217 from Biomarkers csv (multi-source coalesce strategy)
  biomarkers <- fread(biomarkers_csv)

  # Priority order for pTau217 columns (V2 > PLUS > standard > generic)
  priority_cols <- c(
    "r7_QTX_Plasma_pTau217_V2",
    "r7_QTX_Plasma_pTau217_PLUS",
    "r7_QTX_Plasma_pTau217",
    "pTau217",
    "Plasma_pTau217"
  )
  available_ptau_cols <- intersect(priority_cols, colnames(biomarkers))

  if (length(available_ptau_cols) == 0) {
    ptau_cols <- grep("217", grep("pTau|ptau|PTAU", colnames(biomarkers),
                                   value = TRUE, ignore.case = TRUE),
                      value = TRUE)
    available_ptau_cols <- ptau_cols
  }

  if (length(available_ptau_cols) == 0) {
    stop("No pTau217 columns found in HABS Biomarkers file")
  }

  cat(sprintf("  HABS pTau217 columns found: %s\n",
              paste(available_ptau_cols, collapse = ", ")))

  ptau_df <- biomarkers[, c("Med_ID", available_ptau_cols), with = FALSE]

  # Convert to numeric and clean outliers (<=0 or >=1000)
  for (col in available_ptau_cols) {
    ptau_df[[col]] <- as.numeric(ptau_df[[col]])
    ptau_df[[col]][ptau_df[[col]] <= 0 | ptau_df[[col]] >= 1000] <- NA
  }

  # Coalesce: use first non-NA value in priority order
  if (length(available_ptau_cols) == 1) {
    ptau_df$pTau217_Raw <- ptau_df[[available_ptau_cols[1]]]
  } else {
    ptau_df$pTau217_Raw <- do.call(coalesce, lapply(available_ptau_cols,
                                                     function(x) ptau_df[[x]]))
  }

  ptau_df <- ptau_df[!is.na(pTau217_Raw)][!duplicated(Med_ID)]
  ptau_df$Log_pTau217 <- log(ptau_df$pTau217_Raw)

  cat(sprintf("  HABS pTau217: %d samples after coalesce (%.1f%%)\n",
              nrow(ptau_df), 100 * nrow(ptau_df) / nrow(biomarkers)))

  # 3) WMH from xlsx with QC filtering (02_WMH_Fail_QC == "pass")
  wmh_raw <- read_xlsx(wmh_xlsx, sheet = wmh_sheet)
  wmh <- wmh_raw %>%
    filter(tolower(`02_WMH_Fail_QC`) == "pass") %>%
    distinct(Med_ID, .keep_all = TRUE) %>%
    dplyr::select(Med_ID,
                  Total_WMH = `02_WMH_Volume_Raw`,
                  Log_WMH_Raw = `02_WMH_Volume_Log`) %>%
    mutate(Med_ID = as.character(Med_ID),
           Total_WMH = as.numeric(as.character(Total_WMH)),
           Log_WMH_Raw = as.numeric(as.character(Log_WMH_Raw))) %>%
    filter(!is.na(Log_WMH_Raw) & is.finite(Log_WMH_Raw))

  wmh_sd <- sd(wmh$Log_WMH_Raw, na.rm = TRUE)
  wmh$Log_WMH_z <- as.numeric(scale(wmh$Log_WMH_Raw))

  # 4) Clinical data (MMSE_Total, CDR_Global)
  clinical <- fread(clinical_csv)
  clinical_bl <- clinical[!duplicated(Med_ID)]
  clinical_data <- data.frame(
    Med_ID = clinical_bl$Med_ID,
    MMSE = as.numeric(clinical_bl$MMSE_Total),
    CDR = as.numeric(clinical_bl$CDR_Global)
  )
  clinical_data$MMSE[clinical_data$MMSE < 0] <- NA

  # --- Merge all on Med_ID ---
  merged <- as.data.frame(apoe)
  merged <- merge(merged, as.data.frame(ptau_df), by = "Med_ID", all.x = TRUE)
  merged <- merge(merged, as.data.frame(wmh), by = "Med_ID", all.x = TRUE)
  merged <- merge(merged, clinical_data, by = "Med_ID", all.x = TRUE)

  # Construct AD_Conversion per manuscript: CDR >= 0.5 OR MMSE < 24
  merged$AD_Conversion <- as.integer(
    (!is.na(merged$CDR) & merged$CDR >= 0.5) |
    (!is.na(merged$MMSE) & merged$MMSE < 24)
  )
  # Set to NA if both CDR and MMSE are missing
  merged$AD_Conversion[is.na(merged$CDR) & is.na(merged$MMSE)] <- NA

  # Z-scored variables for mediation/SEM (per manuscript: all continuous z-scored)
  merged$Log_pTau217_z <- as.numeric(scale(merged$Log_pTau217))
  merged$MMSE_z <- as.numeric(scale(merged$MMSE))

  cat(sprintf("  HABS: N=%d, WMH SD(log)=%.3f\n", nrow(merged), wmh_sd))
  cat(sprintf("  pTau217 available: %d (%.1f%%)\n",
              sum(!is.na(merged$Log_pTau217)),
              100 * mean(!is.na(merged$Log_pTau217))))
  cat(sprintf("  WMH available: %d (%.1f%%)\n",
              sum(!is.na(merged$Log_WMH_z)),
              100 * mean(!is.na(merged$Log_WMH_z))))
  cat(sprintf("  AD_Conversion events: %d\n",
              sum(merged$AD_Conversion, na.rm = TRUE)))

  attr(merged, "wmh_sd") <- wmh_sd
  return(merged)
}


# --- 1.3 AIBL Cohort ---
# Reads 3 files: ptdemog, apoeres, pdxconv
# Time unit: MONTHS (from VISCODE m06, m18, m36, etc.)
load_aibl_data <- function(ptdemog_file, apoeres_file, pdxconv_file,
                           baseline_year = 2008) {

  baseline_codes <- c("bl", "BL", "M00", "m00", "0", "sc", "screening")

  # 1) Demographics from ptdemog (baseline)
  demog <- fread(ptdemog_file)
  demog_bl <- demog[VISCODE %in% baseline_codes]
  if (nrow(demog_bl) == 0) demog_bl <- demog[, .SD[1], by = RID]
  demog_bl <- demog_bl[!duplicated(RID)]

  # 2) APOE from apoeres (APGEN1/APGEN2 -> APOE4 carrier)
  apoe <- fread(apoeres_file)
  apoe$APOE4_Carrier <- as.integer(apoe$APGEN1 == 4 | apoe$APGEN2 == 4)
  apoe <- apoe[!duplicated(RID)]

  # 3) Conversion tracking from pdxconv (VISCODE -> months)
  dx <- fread(pdxconv_file)
  dx$Visit_Month <- NA_real_
  dx$Visit_Month[tolower(dx$VISCODE) == "bl"] <- 0
  month_pattern <- grepl("^m[0-9]+$", tolower(dx$VISCODE))
  dx$Visit_Month[month_pattern] <- as.numeric(
    gsub("^m", "", tolower(dx$VISCODE[month_pattern])))

  followup <- dx[!is.na(Visit_Month), .(
    Max_Followup = max(Visit_Month, na.rm = TRUE),
    Conversion_Month = {
      ad_visits <- Visit_Month[DXCURREN == 3 & !is.na(DXCURREN)]
      if (length(ad_visits) > 0) min(ad_visits) else NA_real_
    },
    Ever_AD = as.integer(any(DXCURREN == 3, na.rm = TRUE))
  ), by = RID]

  # Time = conversion month if converted, else max follow-up (in MONTHS)
  followup$Time <- ifelse(
    !is.na(followup$Conversion_Month) & followup$Conversion_Month > 0,
    followup$Conversion_Month,
    followup$Max_Followup)
  followup$Event <- followup$Ever_AD

  # --- Merge ---
  merged <- data.frame(RID = demog_bl$RID)

  # Explicit gender recoding (AIBL: 1=Male, 2=Female)
  merged$Gender <- factor(
    ifelse(demog_bl$PTGENDER == 1, "Male",
           ifelse(demog_bl$PTGENDER == 2, "Female", NA)),
    levels = c("Male", "Female"))

  # Age calculation from PTDOB (extract 4-digit year, subtract from baseline)
  birth_years <- as.numeric(str_extract(as.character(demog_bl$PTDOB), "\\d{4}"))
  merged$Age <- baseline_year - birth_years
  merged$Age[merged$Age < 40 | merged$Age > 110 | is.na(merged$Age)] <- NA

  # Merge APOE
  merged <- merge(merged, apoe[, .(RID, APOE4_Carrier)],
                  by = "RID", all.x = TRUE)

  # Merge follow-up
  merged <- merge(merged, followup, by = "RID", all.x = TRUE)

  # Remove Time == 0 samples (QC: no follow-up or baseline-only events)
  n_before <- nrow(merged)
  merged <- merged[is.na(merged$Time) | merged$Time > 0, ]
  n_removed <- n_before - nrow(merged)
  if (n_removed > 0) {
    cat(sprintf("  AIBL: Removed %d samples with Time==0 (QC)\n", n_removed))
  }

  cat(sprintf("  AIBL: N=%d, Events=%d (%.1f%%)\n",
              nrow(merged), sum(merged$Event, na.rm = TRUE),
              100 * mean(merged$Event, na.rm = TRUE)))
  cat(sprintf("  Age: %.1f +/- %.1f (range: %.0f-%.0f)\n",
              mean(merged$Age, na.rm = TRUE), sd(merged$Age, na.rm = TRUE),
              min(merged$Age, na.rm = TRUE), max(merged$Age, na.rm = TRUE)))
  cat(sprintf("  Time unit: MONTHS (median=%.1f, max=%.1f)\n",
              median(merged$Time, na.rm = TRUE),
              max(merged$Time, na.rm = TRUE)))

  return(merged)
}


# ==============================================================================
# 2. A4: WMH-COGNITION ASSOCIATION (HC3 ROBUST SE + BP TEST + Q-Q)
# ==============================================================================

run_a4_wmh_cognition <- function(data, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  df <- data[complete.cases(data[, c("Log_WMH_z", "PACC", "Age", "Gender",
                                      "APOE4_Carrier", "Education",
                                      "Centiloid")]), ]

  # Base model (without WMH)
  fit_base <- lm(PACC ~ Age + Gender + APOE4_Carrier + Education + Centiloid,
                 data = df)

  # Full model (with standardized WMH)
  fit_full <- lm(PACC ~ Log_WMH_z + Age + Gender + APOE4_Carrier +
                   Education + Centiloid, data = df)
  s_ols <- summary(fit_full)$coefficients["Log_WMH_z", ]

  # HC3 robust standard errors
  fit_hc3 <- coeftest(fit_full, vcov = vcovHC(fit_full, type = "HC3"))
  s_hc3 <- fit_hc3["Log_WMH_z", ]

  # R-squared change
  r2_base <- summary(fit_base)$r.squared
  r2_full <- summary(fit_full)$r.squared

  # Breusch-Pagan heteroscedasticity test
  bp_test <- bptest(fit_full)

  # Q-Q diagnostic plot
  pdf(file.path(output_dir, "A4_qq_plot.pdf"), width = 6, height = 6)
  qqnorm(residuals(fit_full), main = "A4: Normal Q-Q Plot of Residuals")
  qqline(residuals(fit_full), col = "red", lwd = 2)
  dev.off()

  list(
    ols = data.frame(
      Beta = s_ols["Estimate"], SE = s_ols["Std. Error"],
      t = s_ols["t value"], P = s_ols["Pr(>|t|)"],
      CI_Lower = s_ols["Estimate"] - 1.96 * s_ols["Std. Error"],
      CI_Upper = s_ols["Estimate"] + 1.96 * s_ols["Std. Error"],
      N = nrow(df), stringsAsFactors = FALSE),
    hc3 = data.frame(
      Beta = s_hc3["Estimate"], SE = s_hc3["Std. Error"],
      t = s_hc3["t value"], P = s_hc3["Pr(>|t|)"],
      CI_Lower = s_hc3["Estimate"] - 1.96 * s_hc3["Std. Error"],
      CI_Upper = s_hc3["Estimate"] + 1.96 * s_hc3["Std. Error"],
      N = nrow(df), stringsAsFactors = FALSE),
    model_fit = data.frame(
      R2_Base = r2_base, R2_Full = r2_full,
      Delta_R2 = r2_full - r2_base, stringsAsFactors = FALSE),
    diagnostics = data.frame(
      BP_Statistic = as.numeric(bp_test$statistic),
      BP_P = bp_test$p.value,
      Heteroscedasticity = bp_test$p.value < 0.05,
      stringsAsFactors = FALSE),
    wmh_sd = attr(data, "wmh_sd")
  )
}


# ==============================================================================
# 3. HABS: WMH-pTau217 ASSOCIATION (HC3 + BP + Q-Q)
# ==============================================================================

run_habs_wmh_ptau217 <- function(data, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  df <- data[complete.cases(data[, c("Log_WMH_z", "Log_pTau217", "Age",
                                      "Gender", "APOE4_Carrier",
                                      "Education")]), ]
  df <- df[is.finite(df$Log_pTau217), ]

  # Base model
  fit_base <- lm(Log_pTau217 ~ Age + Gender + APOE4_Carrier + Education,
                 data = df)

  # Full model
  fit_full <- lm(Log_pTau217 ~ Log_WMH_z + Age + Gender + APOE4_Carrier +
                   Education, data = df)
  s_ols <- summary(fit_full)$coefficients["Log_WMH_z", ]

  # HC3 robust standard errors
  fit_hc3 <- coeftest(fit_full, vcov = vcovHC(fit_full, type = "HC3"))
  s_hc3 <- fit_hc3["Log_WMH_z", ]

  # R-squared change
  r2_base <- summary(fit_base)$r.squared
  r2_full <- summary(fit_full)$r.squared

  # Breusch-Pagan test
  bp_test <- bptest(fit_full)

  # Q-Q diagnostic plot
  pdf(file.path(output_dir, "HABS_qq_plot.pdf"), width = 6, height = 6)
  qqnorm(residuals(fit_full), main = "HABS: Normal Q-Q Plot of Residuals")
  qqline(residuals(fit_full), col = "red", lwd = 2)
  dev.off()

  # APOE4-stratified analysis
  strat_results <- list()
  for (apoe in c(0, 1)) {
    sub <- df[df$APOE4_Carrier == apoe, ]
    if (nrow(sub) < 20) next
    fit <- lm(Log_pTau217 ~ Log_WMH_z + Age + Gender + Education, data = sub)
    s <- summary(fit)$coefficients["Log_WMH_z", ]
    strat_results[[length(strat_results) + 1]] <- data.frame(
      APOE4 = ifelse(apoe == 1, "Carrier", "Non-carrier"),
      Beta = s["Estimate"], SE = s["Std. Error"],
      P = s["Pr(>|t|)"], N = nrow(sub), stringsAsFactors = FALSE)
  }

  list(
    ols = data.frame(
      Beta = s_ols["Estimate"], SE = s_ols["Std. Error"],
      P = s_ols["Pr(>|t|)"], N = nrow(df), stringsAsFactors = FALSE),
    hc3 = data.frame(
      Beta = s_hc3["Estimate"], SE = s_hc3["Std. Error"],
      P = s_hc3["Pr(>|t|)"], N = nrow(df), stringsAsFactors = FALSE),
    model_fit = data.frame(
      R2_Base = r2_base, R2_Full = r2_full,
      Delta_R2 = r2_full - r2_base, stringsAsFactors = FALSE),
    diagnostics = data.frame(
      BP_Statistic = as.numeric(bp_test$statistic),
      BP_P = bp_test$p.value,
      Heteroscedasticity = bp_test$p.value < 0.05,
      stringsAsFactors = FALSE),
    apoe_stratified = do.call(rbind, strat_results),
    wmh_sd = attr(data, "wmh_sd")
  )
}


# ==============================================================================
# 4. HABS: MEDIATION ANALYSIS (Z-SCORED, 5000 BOOTSTRAP)
# ==============================================================================
# Per manuscript: "all continuous variables standardized to z-score"
# Uses Log_WMH_z (already z-scored), Log_pTau217_z, MMSE_z

run_habs_mediation <- function(data, n_boot = 5000) {

  covars <- c("Age", "Gender", "APOE4_Carrier", "Education")
  req_vars <- c("Log_WMH_z", "Log_pTau217", "MMSE", covars)
  df <- data[complete.cases(data[, req_vars]), ]
  df <- df[is.finite(df$Log_pTau217), ]

  # Z-score all continuous variables for mediation
  df$Log_pTau217_z <- as.numeric(scale(df$Log_pTau217))
  df$MMSE_z <- as.numeric(scale(df$MMSE))
  df$Log_WMH_z <- as.numeric(scale(df$Log_WMH_z))  # re-scale within complete cases
  df$Age_z <- as.numeric(scale(df$Age))
  df$Education_z <- as.numeric(scale(df$Education))

  # Path a: WMH -> pTau217 (z-scored)
  fit_m <- lm(Log_pTau217_z ~ Log_WMH_z + Age_z + Gender + APOE4_Carrier +
                Education_z, data = df)

  # Path b + c': pTau217 + WMH -> MMSE (z-scored)
  fit_y <- lm(MMSE_z ~ Log_pTau217_z + Log_WMH_z + Age_z + Gender +
                APOE4_Carrier + Education_z, data = df)

  # Total effect (path c): WMH -> MMSE
  fit_total <- lm(MMSE_z ~ Log_WMH_z + Age_z + Gender + APOE4_Carrier +
                    Education_z, data = df)

  # Bootstrap mediation (Imai et al., 5000 iterations)
  set.seed(2026)
  med_result <- mediate(fit_m, fit_y, treat = "Log_WMH_z",
                        mediator = "Log_pTau217_z",
                        boot = TRUE, sims = n_boot,
                        boot.ci.type = "bca")

  # Extract path coefficients
  path_a <- summary(fit_m)$coefficients["Log_WMH_z", ]
  path_b <- summary(fit_y)$coefficients["Log_pTau217_z", ]
  path_c_prime <- summary(fit_y)$coefficients["Log_WMH_z", ]
  path_c <- summary(fit_total)$coefficients["Log_WMH_z", ]

  list(
    paths = data.frame(
      Path = c("a (WMH->pTau217)", "b (pTau217->MMSE)",
               "c' (WMH->MMSE direct)", "c (WMH->MMSE total)"),
      Beta = c(path_a["Estimate"], path_b["Estimate"],
               path_c_prime["Estimate"], path_c["Estimate"]),
      SE = c(path_a["Std. Error"], path_b["Std. Error"],
             path_c_prime["Std. Error"], path_c["Std. Error"]),
      P = c(path_a["Pr(>|t|)"], path_b["Pr(>|t|)"],
            path_c_prime["Pr(>|t|)"], path_c["Pr(>|t|)"]),
      stringsAsFactors = FALSE),
    mediation = data.frame(
      ACME = med_result$d0,
      ACME_CI_Lower = med_result$d0.ci[1],
      ACME_CI_Upper = med_result$d0.ci[2],
      ACME_P = med_result$d0.p,
      ADE = med_result$z0,
      ADE_CI_Lower = med_result$z0.ci[1],
      ADE_CI_Upper = med_result$z0.ci[2],
      ADE_P = med_result$z0.p,
      Total = med_result$tau.coef,
      Total_P = med_result$tau.p,
      Prop_Mediated = med_result$n0,
      N = nrow(df), N_Boot = n_boot,
      stringsAsFactors = FALSE),
    note = "All continuous variables z-scored per manuscript methods"
  )
}


# ==============================================================================
# 5. HABS: SEM CONFIRMATION (5000 BOOTSTRAP, Z-SCORED)
# ==============================================================================

run_habs_sem <- function(data, n_boot = 5000) {

  covars <- c("Age", "Gender", "APOE4_Carrier", "Education")
  req_vars <- c("Log_WMH_z", "Log_pTau217", "MMSE", covars)
  df <- data[complete.cases(data[, req_vars]), ]
  df <- df[is.finite(df$Log_pTau217), ]

  # Z-score all continuous variables
  df$Log_pTau217_z <- as.numeric(scale(df$Log_pTau217))
  df$MMSE_z <- as.numeric(scale(df$MMSE))
  df$Log_WMH_z <- as.numeric(scale(df$Log_WMH_z))
  df$Age_z <- as.numeric(scale(df$Age))
  df$Education_z <- as.numeric(scale(df$Education))
  df$Gender_Num <- as.numeric(df$Gender) - 1

  sem_model <- '
    # Direct paths (all z-scored continuous variables)
    Log_pTau217_z ~ a * Log_WMH_z + Age_z + Gender_Num + APOE4_Carrier + Education_z
    MMSE_z ~ b * Log_pTau217_z + cp * Log_WMH_z + Age_z + Gender_Num + APOE4_Carrier + Education_z

    # Indirect and total effects
    indirect := a * b
    total := cp + a * b
    prop_mediated := (a * b) / (cp + a * b)
  '

  set.seed(2026)
  fit <- sem(sem_model, data = df, se = "bootstrap", bootstrap = n_boot,
             estimator = "ML")

  params <- parameterEstimates(fit, boot.ci.type = "perc")
  fit_indices <- fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr",
                                     "chisq", "df", "pvalue"))

  list(
    parameters = params,
    fit_indices = as.data.frame(t(fit_indices)),
    N = nrow(df),
    note = "5000 bootstrap, all continuous variables z-scored"
  )
}


# ==============================================================================
# 6. HABS: AGE-STRATIFIED MEDIATION AND SEM (<75 vs >=75)
# ==============================================================================

run_habs_age_stratified <- function(data, n_boot = 5000, age_cutoff = 75) {

  covars <- c("Age", "Gender", "APOE4_Carrier", "Education")
  req_vars <- c("Log_WMH_z", "Log_pTau217", "MMSE", covars)
  df <- data[complete.cases(data[, req_vars]), ]
  df <- df[is.finite(df$Log_pTau217), ]

  med_results <- list()
  sem_results <- list()

  for (grp in c("Younger", "Older")) {
    sub <- if (grp == "Younger") df[df$Age < age_cutoff, ] else df[df$Age >= age_cutoff, ]
    if (nrow(sub) < 30) next

    # Z-score within each stratum
    sub$Log_pTau217_z <- as.numeric(scale(sub$Log_pTau217))
    sub$MMSE_z <- as.numeric(scale(sub$MMSE))
    sub$Log_WMH_z <- as.numeric(scale(sub$Log_WMH_z))
    sub$Age_z <- as.numeric(scale(sub$Age))
    sub$Education_z <- as.numeric(scale(sub$Education))
    sub$Gender_Num <- as.numeric(sub$Gender) - 1

    # --- Mediation ---
    fit_m <- lm(Log_pTau217_z ~ Log_WMH_z + Age_z + Gender + APOE4_Carrier +
                  Education_z, data = sub)
    fit_y <- lm(MMSE_z ~ Log_pTau217_z + Log_WMH_z + Age_z + Gender +
                  APOE4_Carrier + Education_z, data = sub)

    set.seed(2026)
    med <- tryCatch(
      mediate(fit_m, fit_y, treat = "Log_WMH_z",
              mediator = "Log_pTau217_z", boot = TRUE, sims = n_boot,
              boot.ci.type = "bca"),
      error = function(e) NULL)

    if (!is.null(med)) {
      path_a <- summary(fit_m)$coefficients["Log_WMH_z", ]
      path_b <- summary(fit_y)$coefficients["Log_pTau217_z", ]

      med_results[[grp]] <- data.frame(
        Age_Group = grp,
        Age_Range = ifelse(grp == "Younger",
                           paste0("<", age_cutoff), paste0(">=", age_cutoff)),
        N = nrow(sub),
        Path_a_Beta = path_a["Estimate"], Path_a_P = path_a["Pr(>|t|)"],
        Path_b_Beta = path_b["Estimate"], Path_b_P = path_b["Pr(>|t|)"],
        ACME = med$d0,
        ACME_CI_Lower = med$d0.ci[1], ACME_CI_Upper = med$d0.ci[2],
        ACME_P = med$d0.p,
        ADE = med$z0, Total = med$tau.coef,
        Prop_Mediated = med$n0,
        stringsAsFactors = FALSE)
    }

    # --- Age-stratified SEM ---
    sem_model <- '
      Log_pTau217_z ~ a * Log_WMH_z + Age_z + Gender_Num + APOE4_Carrier + Education_z
      MMSE_z ~ b * Log_pTau217_z + cp * Log_WMH_z + Age_z + Gender_Num + APOE4_Carrier + Education_z
      indirect := a * b
      total := cp + a * b
      prop_mediated := (a * b) / (cp + a * b)
    '

    set.seed(2026)
    sem_fit <- tryCatch(
      sem(sem_model, data = sub, se = "bootstrap", bootstrap = n_boot,
          estimator = "ML"),
      error = function(e) NULL)

    if (!is.null(sem_fit)) {
      params <- parameterEstimates(sem_fit, boot.ci.type = "perc")
      indirect_row <- params[params$label == "indirect", ]
      sem_results[[grp]] <- data.frame(
        Age_Group = grp, N = nrow(sub),
        SEM_Indirect = indirect_row$est,
        SEM_Indirect_CI_Lower = indirect_row$ci.lower,
        SEM_Indirect_CI_Upper = indirect_row$ci.upper,
        SEM_Indirect_P = indirect_row$pvalue,
        stringsAsFactors = FALSE)
    }
  }

  list(
    mediation = do.call(rbind, med_results),
    sem = do.call(rbind, sem_results)
  )
}


# ==============================================================================
# 7. HABS: AGE INTERACTION TESTS
# ==============================================================================

run_habs_age_interaction <- function(data) {

  covars <- c("Gender", "APOE4_Carrier", "Education")
  df <- data[complete.cases(data[, c("Log_WMH_z", "Log_pTau217", "MMSE",
                                      "Age", covars)]), ]
  df <- df[is.finite(df$Log_pTau217), ]

  # Age x WMH interaction on pTau217
  fit_wmh_int <- lm(Log_pTau217 ~ Log_WMH_z * Age + Gender +
                      APOE4_Carrier + Education, data = df)
  s_wmh <- summary(fit_wmh_int)$coefficients

  # Age x pTau217 interaction on MMSE
  fit_ptau_int <- lm(MMSE ~ Log_pTau217 * Age + Gender +
                       APOE4_Carrier + Education, data = df)
  s_ptau <- summary(fit_ptau_int)$coefficients

  results <- list()
  if ("Log_WMH_z:Age" %in% rownames(s_wmh)) {
    results$wmh_age <- data.frame(
      Interaction = "WMH x Age -> pTau217",
      Beta = s_wmh["Log_WMH_z:Age", "Estimate"],
      SE = s_wmh["Log_WMH_z:Age", "Std. Error"],
      P = s_wmh["Log_WMH_z:Age", "Pr(>|t|)"],
      N = nrow(df), stringsAsFactors = FALSE)
  }
  if ("Log_pTau217:Age" %in% rownames(s_ptau)) {
    results$ptau_age <- data.frame(
      Interaction = "pTau217 x Age -> MMSE",
      Beta = s_ptau["Log_pTau217:Age", "Estimate"],
      SE = s_ptau["Log_pTau217:Age", "Std. Error"],
      P = s_ptau["Log_pTau217:Age", "Pr(>|t|)"],
      N = nrow(df), stringsAsFactors = FALSE)
  }

  do.call(rbind, results)
}


# ==============================================================================
# 8. HABS: CLINICAL UTILITY (TRAIN/TEST, FIRTH, AUC, PR-AUC, NRI/IDI, DCA)
# ==============================================================================
# Per manuscript: train/test split, PR-AUC, NRI multi-threshold with bootstrap,
# IDI with bootstrap CI, DCA to 0.80

run_habs_clinical_utility <- function(data, seed = 2026, n_boot_nri = 1000) {

  # AD_Conversion constructed in load_habs_data: CDR>=0.5 OR MMSE<24
  if (!"AD_Conversion" %in% colnames(data) || all(is.na(data$AD_Conversion))) {
    return(list(message = "AD_Conversion variable not available"))
  }

  covars <- c("Age", "Gender", "APOE4_Carrier", "Education")
  df <- data[complete.cases(data[, c("Log_WMH_z", "Log_pTau217",
                                      "AD_Conversion", covars)]), ]
  if (nrow(df) < 30 || sum(df$AD_Conversion) < 5) {
    return(list(message = "Insufficient events for clinical utility analysis"))
  }

  # --- Train/test split (70/30, fixed seed) ---
  set.seed(seed)
  train_idx <- sample(seq_len(nrow(df)), size = floor(0.7 * nrow(df)))
  train <- df[train_idx, ]
  test <- df[-train_idx, ]

  cat(sprintf("  Train: N=%d (events=%d), Test: N=%d (events=%d)\n",
              nrow(train), sum(train$AD_Conversion),
              nrow(test), sum(test$AD_Conversion)))

  # --- Firth-corrected logistic regression ---
  fit_base <- logistf(AD_Conversion ~ Age + Gender + APOE4_Carrier + Education,
                      data = train)
  fit_full <- logistf(AD_Conversion ~ Age + Gender + APOE4_Carrier + Education +
                        Log_WMH_z + Log_pTau217, data = train)

  # --- Predictions on TEST set ---
  pred_base <- predict(fit_base, newdata = test, type = "response")
  pred_full <- predict(fit_full, newdata = test, type = "response")

  # --- AUC (ROC) ---
  roc_base <- roc(test$AD_Conversion, pred_base, quiet = TRUE)
  roc_full <- roc(test$AD_Conversion, pred_full, quiet = TRUE)
  delong_test <- roc.test(roc_base, roc_full, method = "delong")

  auc_comparison <- data.frame(
    Model = c("Base", "Base+WMH+pTau217"),
    AUC = c(auc(roc_base), auc(roc_full)),
    stringsAsFactors = FALSE)

  # --- PR-AUC ---
  pr_base <- pr.curve(scores.class0 = pred_base[test$AD_Conversion == 1],
                      scores.class1 = pred_base[test$AD_Conversion == 0])
  pr_full <- pr.curve(scores.class0 = pred_full[test$AD_Conversion == 1],
                      scores.class1 = pred_full[test$AD_Conversion == 0])

  prauc_comparison <- data.frame(
    Model = c("Base", "Base+WMH+pTau217"),
    PR_AUC = c(pr_base$auc.integral, pr_full$auc.integral),
    stringsAsFactors = FALSE)

  # --- NRI with bootstrap CI ---
  thresholds <- c(0.10, 0.20, 0.30)
  nri_results <- compute_nri_bootstrap(test$AD_Conversion, pred_base, pred_full,
                                        thresholds = thresholds,
                                        n_boot = n_boot_nri, seed = seed)

  # --- IDI with bootstrap CI ---
  idi_results <- compute_idi_bootstrap(test$AD_Conversion, pred_base, pred_full,
                                       n_boot = n_boot_nri, seed = seed)

  # --- Decision Curve Analysis (0.01 to 0.80 per manuscript) ---
  dca_results <- compute_dca(test$AD_Conversion, pred_base, pred_full,
                             thresholds = seq(0.01, 0.80, by = 0.01))

  list(
    auc = auc_comparison,
    prauc = prauc_comparison,
    delong_p = delong_test$p.value,
    nri = nri_results,
    idi = idi_results,
    dca = dca_results,
    N_Train = nrow(train),
    N_Test = nrow(test),
    N_Events_Train = sum(train$AD_Conversion),
    N_Events_Test = sum(test$AD_Conversion),
    split_seed = seed
  )
}

# --- NRI with bootstrap CI/P ---
compute_nri_bootstrap <- function(outcome, pred_old, pred_new,
                                   thresholds = c(0.10, 0.20, 0.30),
                                   n_boot = 1000, seed = 2026) {

  compute_nri_single <- function(outcome, pred_old, pred_new, thresh) {
    class_old <- as.integer(pred_old >= thresh)
    class_new <- as.integer(pred_new >= thresh)
    events <- outcome == 1
    nonevents <- outcome == 0
    if (sum(events) == 0 || sum(nonevents) == 0) return(c(NA, NA, NA))

    up_event <- sum(class_new[events] > class_old[events])
    down_event <- sum(class_new[events] < class_old[events])
    nri_event <- (up_event - down_event) / sum(events)

    up_nonevent <- sum(class_new[nonevents] > class_old[nonevents])
    down_nonevent <- sum(class_new[nonevents] < class_old[nonevents])
    nri_nonevent <- (down_nonevent - up_nonevent) / sum(nonevents)

    c(nri_event, nri_nonevent, nri_event + nri_nonevent)
  }

  results <- list()
  for (thresh in thresholds) {
    obs_nri <- compute_nri_single(outcome, pred_old, pred_new, thresh)

    # Bootstrap
    set.seed(seed)
    boot_nri <- replicate(n_boot, {
      idx <- sample(length(outcome), replace = TRUE)
      compute_nri_single(outcome[idx], pred_old[idx], pred_new[idx], thresh)
    })

    boot_total <- boot_nri[3, ]
    boot_total <- boot_total[!is.na(boot_total)]

    results[[length(results) + 1]] <- data.frame(
      Threshold = thresh,
      NRI_Event = obs_nri[1], NRI_NonEvent = obs_nri[2],
      NRI_Total = obs_nri[3],
      NRI_CI_Lower = quantile(boot_total, 0.025, na.rm = TRUE),
      NRI_CI_Upper = quantile(boot_total, 0.975, na.rm = TRUE),
      NRI_P = mean(boot_total <= 0, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }

  do.call(rbind, results)
}

# --- IDI with bootstrap CI/P ---
compute_idi_bootstrap <- function(outcome, pred_old, pred_new,
                                   n_boot = 1000, seed = 2026) {

  compute_idi_single <- function(outcome, pred_old, pred_new) {
    idi_event <- mean(pred_new[outcome == 1]) - mean(pred_old[outcome == 1])
    idi_nonevent <- mean(pred_new[outcome == 0]) - mean(pred_old[outcome == 0])
    c(idi_event - idi_nonevent, idi_event, idi_nonevent)
  }

  obs_idi <- compute_idi_single(outcome, pred_old, pred_new)

  set.seed(seed)
  boot_idi <- replicate(n_boot, {
    idx <- sample(length(outcome), replace = TRUE)
    compute_idi_single(outcome[idx], pred_old[idx], pred_new[idx])
  })

  boot_total <- boot_idi[1, ]

  data.frame(
    IDI = obs_idi[1], IDI_Event = obs_idi[2], IDI_NonEvent = obs_idi[3],
    IDI_CI_Lower = quantile(boot_total, 0.025, na.rm = TRUE),
    IDI_CI_Upper = quantile(boot_total, 0.975, na.rm = TRUE),
    IDI_P = mean(boot_total <= 0, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

# --- DCA (0.01 to 0.80 per manuscript) ---
compute_dca <- function(outcome, pred_base, pred_full,
                        thresholds = seq(0.01, 0.80, by = 0.01)) {

  n <- length(outcome)
  prevalence <- mean(outcome)

  dca_results <- list()
  for (pt in thresholds) {
    odds <- pt / (1 - pt)
    nb_all <- prevalence - (1 - prevalence) * odds

    tp_base <- sum(pred_base >= pt & outcome == 1)
    fp_base <- sum(pred_base >= pt & outcome == 0)
    nb_base <- tp_base / n - fp_base / n * odds

    tp_full <- sum(pred_full >= pt & outcome == 1)
    fp_full <- sum(pred_full >= pt & outcome == 0)
    nb_full <- tp_full / n - fp_full / n * odds

    dca_results[[length(dca_results) + 1]] <- data.frame(
      Threshold = pt,
      NB_TreatAll = nb_all, NB_Base = nb_base, NB_Full = nb_full,
      stringsAsFactors = FALSE)
  }

  do.call(rbind, dca_results)
}


# ==============================================================================
# 9. AIBL: SURVIVAL ANALYSIS (COX PH, KM, AFT; TIME IN MONTHS)
# ==============================================================================

run_aibl_survival <- function(data) {

  df <- data[complete.cases(data[, c("Age", "Gender", "APOE4_Carrier",
                                      "Time", "Event")]), ]

  # Remove Time == 0 (no follow-up, QC per original code)
  df <- df[df$Time > 0, ]

  # Time is in MONTHS (from VISCODE m06, m18, m36, etc.)
  cat(sprintf("  AIBL survival: Time in MONTHS (median=%.1f, max=%.1f)\n",
              median(df$Time), max(df$Time)))

  # Cox proportional hazards
  cox_fit <- coxph(Surv(Time, Event) ~ APOE4_Carrier + Age + Gender, data = df)
  cox_summary <- summary(cox_fit)

  cox_results <- data.frame(
    Variable = rownames(cox_summary$coefficients),
    HR = cox_summary$coefficients[, "exp(coef)"],
    HR_Lower = cox_summary$conf.int[, "lower .95"],
    HR_Upper = cox_summary$conf.int[, "upper .95"],
    P = cox_summary$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE)

  # PH assumption test (Schoenfeld residuals)
  ph_test <- cox.zph(cox_fit)
  ph_results <- data.frame(
    Variable = rownames(ph_test$table),
    Chisq = ph_test$table[, "chisq"],
    P = ph_test$table[, "p"],
    stringsAsFactors = FALSE)

  # Kaplan-Meier by APOE4
  km_fit <- survfit(Surv(Time, Event) ~ APOE4_Carrier, data = df)
  logrank <- survdiff(Surv(Time, Event) ~ APOE4_Carrier, data = df)
  logrank_p <- 1 - pchisq(logrank$chisq, df = 1)

  # Age-stratified Cox
  age_median <- median(df$Age, na.rm = TRUE)
  strat_results <- list()
  for (grp in c("Younger", "Older")) {
    sub <- if (grp == "Younger") df[df$Age < age_median, ] else df[df$Age >= age_median, ]
    if (nrow(sub) < 20 || sum(sub$Event) < 5) next

    fit <- coxph(Surv(Time, Event) ~ APOE4_Carrier + Age + Gender, data = sub)
    s <- summary(fit)
    apoe_row <- which(rownames(s$coefficients) == "APOE4_Carrier")
    if (length(apoe_row) == 0) next

    strat_results[[grp]] <- data.frame(
      Age_Group = grp,
      HR = s$coefficients[apoe_row, "exp(coef)"],
      HR_Lower = s$conf.int[apoe_row, "lower .95"],
      HR_Upper = s$conf.int[apoe_row, "upper .95"],
      P = s$coefficients[apoe_row, "Pr(>|z|)"],
      N = nrow(sub), N_Events = sum(sub$Event),
      stringsAsFactors = FALSE)
  }

  # AFT model (Weibull)
  aft_fit <- survreg(Surv(Time, Event) ~ APOE4_Carrier + Age + Gender,
                     data = df, dist = "weibull")
  aft_summary <- summary(aft_fit)
  aft_results <- data.frame(
    Variable = rownames(aft_summary$table),
    Coefficient = aft_summary$table[, "Value"],
    SE = aft_summary$table[, "Std. Error"],
    P = aft_summary$table[, "p"],
    TR = exp(aft_summary$table[, "Value"]),
    stringsAsFactors = FALSE)

  list(
    cox = cox_results,
    ph_test = ph_results,
    km_fit = km_fit,
    logrank_p = logrank_p,
    age_stratified = do.call(rbind, strat_results),
    aft = aft_results,
    N = nrow(df),
    N_Events = sum(df$Event),
    Median_Followup_Months = median(df$Time, na.rm = TRUE),
    Time_Unit = "months"
  )
}


# ==============================================================================
# 10. MISSING DATA AND SELECTION BIAS ASSESSMENT
# ==============================================================================
# Uses SMD (Standardized Mean Difference) for both numeric and categorical
# variables. SMD < 0.1 = negligible, 0.1-0.2 = small, 0.2-0.5 = moderate

run_missing_data_assessment <- function(a4_data, habs_data, aibl_data) {

  # --- Helper: compute SMD for numeric variable ---
  smd_numeric <- function(x1, x2) {
    m1 <- mean(x1, na.rm = TRUE)
    m2 <- mean(x2, na.rm = TRUE)
    s1 <- sd(x1, na.rm = TRUE)
    s2 <- sd(x2, na.rm = TRUE)
    n1 <- sum(!is.na(x1))
    n2 <- sum(!is.na(x2))
    pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
    if (pooled_sd == 0) return(0)
    (m1 - m2) / pooled_sd
  }

  # --- Helper: compute SMD for binary/categorical variable ---
  smd_binary <- function(x1, x2) {
    p1 <- mean(x1, na.rm = TRUE)
    p2 <- mean(x2, na.rm = TRUE)
    pooled_p <- (sum(x1, na.rm = TRUE) + sum(x2, na.rm = TRUE)) /
                (sum(!is.na(x1)) + sum(!is.na(x2)))
    denom <- sqrt(pooled_p * (1 - pooled_p))
    if (denom == 0) return(0)
    (p1 - p2) / denom
  }

  # --- Helper: assess one cohort ---
  assess_cohort <- function(data, complete_flag, cohort_name,
                            numeric_vars, binary_vars) {
    included <- data[complete_flag, ]
    excluded <- data[!complete_flag, ]

    if (nrow(excluded) == 0) {
      return(list(
        summary = data.frame(
          Cohort = cohort_name, Total_N = nrow(data),
          Complete_N = nrow(included), Excluded_N = 0,
          Completion_Rate = 100.0, Max_SMD = 0, Mean_SMD = 0,
          Bias_Level = "None", stringsAsFactors = FALSE),
        details = data.frame(
          Cohort = cohort_name, Variable = "N/A", SMD = 0,
          Note = "100% complete cases", stringsAsFactors = FALSE)
      ))
    }

    smd_vals <- c()
    detail_rows <- list()

    # Numeric variables
    for (v in numeric_vars) {
      if (v %in% colnames(data) && sum(!is.na(included[[v]])) > 1 &&
          sum(!is.na(excluded[[v]])) > 1) {
        s <- smd_numeric(included[[v]], excluded[[v]])
        smd_vals <- c(smd_vals, abs(s))
        detail_rows[[length(detail_rows) + 1]] <- data.frame(
          Cohort = cohort_name, Variable = v, SMD = round(s, 3),
          Type = "Numeric", stringsAsFactors = FALSE)
      }
    }

    # Binary/categorical variables (Gender as Female proportion, APOE4)
    for (v in binary_vars) {
      if (v %in% colnames(data)) {
        vals_inc <- included[[v]]
        vals_exc <- excluded[[v]]

        # Convert factor Gender to numeric (Female=1)
        if (is.factor(vals_inc)) {
          vals_inc <- as.integer(vals_inc == "Female")
          vals_exc <- as.integer(vals_exc == "Female")
        }

        if (sum(!is.na(vals_inc)) > 0 && sum(!is.na(vals_exc)) > 0) {
          s <- smd_binary(vals_inc, vals_exc)
          smd_vals <- c(smd_vals, abs(s))
          detail_rows[[length(detail_rows) + 1]] <- data.frame(
            Cohort = cohort_name, Variable = v, SMD = round(s, 3),
            Type = "Categorical", stringsAsFactors = FALSE)
        }
      }
    }

    max_smd <- if (length(smd_vals) > 0) max(smd_vals) else NA
    mean_smd <- if (length(smd_vals) > 0) mean(smd_vals) else NA

    bias_level <- if (is.na(max_smd)) "Unknown"
                  else if (max_smd < 0.1) "Negligible"
                  else if (max_smd < 0.2) "Small"
                  else if (max_smd < 0.5) "Moderate"
                  else "Large"

    list(
      summary = data.frame(
        Cohort = cohort_name, Total_N = nrow(data),
        Complete_N = nrow(included), Excluded_N = nrow(excluded),
        Completion_Rate = round(100 * nrow(included) / nrow(data), 1),
        Max_SMD = round(max_smd, 3), Mean_SMD = round(mean_smd, 3),
        Bias_Level = bias_level, stringsAsFactors = FALSE),
      details = do.call(rbind, detail_rows)
    )
  }

  # --- A4: complete case definition ---
  a4_complete <- complete.cases(a4_data[, c("Log_WMH_z", "PACC", "Age",
                                             "Gender", "APOE4_Carrier",
                                             "Education", "Centiloid")])
  a4_res <- assess_cohort(a4_data, a4_complete, "A4",
                          numeric_vars = c("Age", "Education", "Centiloid",
                                           "Log_WMH_Raw", "PACC"),
                          binary_vars = c("Gender", "APOE4_Carrier"))

  # --- HABS: complete case definition ---
  habs_complete <- complete.cases(habs_data[, c("Log_WMH_z", "Log_pTau217",
                                                 "Age", "Gender",
                                                 "APOE4_Carrier", "Education")])
  habs_res <- assess_cohort(habs_data, habs_complete, "HABS",
                            numeric_vars = c("Age", "Education", "Log_WMH_Raw",
                                             "Log_pTau217", "MMSE"),
                            binary_vars = c("Gender", "APOE4_Carrier"))

  # --- AIBL: complete case definition ---
  aibl_complete <- complete.cases(aibl_data[, c("Age", "Gender",
                                                 "APOE4_Carrier",
                                                 "Time", "Event")])
  aibl_res <- assess_cohort(aibl_data, aibl_complete, "AIBL",
                            numeric_vars = c("Age", "Time"),
                            binary_vars = c("Gender", "APOE4_Carrier"))

  list(
    summary = rbind(a4_res$summary, habs_res$summary, aibl_res$summary),
    details = rbind(a4_res$details, habs_res$details, aibl_res$details)
  )
}


# ==============================================================================
# 11. VISUALIZATION FUNCTIONS
# ==============================================================================

# --- 11.1 Cross-cohort forest plot ---
plot_forest <- function(a4_results, habs_results, aibl_results,
                        output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Build forest plot data from results
  forest_data <- data.frame(
    Cohort = c("A4 Study", "HABS"),
    Outcome = c("Cognitive Performance (PACC)",
                "Tau Pathology (Plasma pTau217)"),
    Beta = c(a4_results$hc3$Beta, habs_results$hc3$Beta),
    SE = c(a4_results$hc3$SE, habs_results$hc3$SE),
    P = c(a4_results$hc3$P, habs_results$hc3$P),
    N = c(a4_results$hc3$N, habs_results$hc3$N),
    stringsAsFactors = FALSE
  )
  forest_data$CI_Lower <- forest_data$Beta - 1.96 * forest_data$SE
  forest_data$CI_Upper <- forest_data$Beta + 1.96 * forest_data$SE
  forest_data$Significant <- forest_data$P < 0.05
  forest_data$Label <- sprintf("%s (N=%d)\n%s",
                               forest_data$Cohort, forest_data$N,
                               forest_data$Outcome)

  # Add AIBL APOE4 HR if available
  if (!is.null(aibl_results$cox)) {
    apoe_row <- aibl_results$cox[aibl_results$cox$Variable == "APOE4_Carrier", ]
    if (nrow(apoe_row) > 0) {
      # Convert HR to log scale for forest plot
      aibl_entry <- data.frame(
        Cohort = "AIBL",
        Outcome = sprintf("AD Conversion (HR=%.2f)", apoe_row$HR),
        Beta = log(apoe_row$HR),
        SE = (log(apoe_row$HR_Upper) - log(apoe_row$HR_Lower)) / (2 * 1.96),
        P = apoe_row$P,
        N = aibl_results$N,
        CI_Lower = log(apoe_row$HR_Lower),
        CI_Upper = log(apoe_row$HR_Upper),
        Significant = apoe_row$P < 0.05,
        Label = sprintf("AIBL (N=%d)\nAD Conversion (Cox PH)",
                        aibl_results$N),
        stringsAsFactors = FALSE
      )
      forest_data <- rbind(forest_data, aibl_entry)
    }
  }

  color_sig <- "#0072B2"
  color_nonsig <- "#999999"

  p <- ggplot(forest_data, aes(x = Beta, y = reorder(Label, Beta))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#D55E00",
               linewidth = 0.8, alpha = 0.7) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper,
                       color = Significant),
                   height = 0.25, linewidth = 1.2) +
    geom_point(aes(color = Significant), size = 5, shape = 18) +
    scale_color_manual(values = c("TRUE" = color_sig,
                                  "FALSE" = color_nonsig),
                       guide = "none") +
    labs(title = "WMH Effects Across Disease Stages",
         subtitle = "Per 1 SD increase in log-transformed WMH volume",
         x = "Standardized Beta (95% CI)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          panel.grid.minor = element_blank())

  ggsave(file.path(output_dir, "forest_plot_cross_cohort.pdf"),
         p, width = 10, height = 6)
  return(p)
}


# --- 11.2 Mediation path diagram ---
plot_mediation_path <- function(mediation_results, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  paths <- mediation_results$paths
  med <- mediation_results$mediation

  pdf(file.path(output_dir, "mediation_path_diagram.pdf"),
      width = 10, height = 7)
  par(mar = c(2, 2, 3, 2))
  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 7))

  # Boxes
  rect(0.5, 2.5, 3.5, 3.5, border = "black", lwd = 2)
  text(2, 3, "WMH\n(Log, z-scored)", cex = 1.0)

  rect(3.5, 5.0, 6.5, 6.0, border = "black", lwd = 2)
  text(5, 5.5, "pTau217\n(Log, z-scored)", cex = 1.0)

  rect(6.5, 2.5, 9.5, 3.5, border = "black", lwd = 2)
  text(8, 3, "MMSE\n(z-scored)", cex = 1.0)

  # Path a: WMH -> pTau217
  arrows(3.5, 3.3, 3.8, 5.0, lwd = 2, col = "blue")
  a_beta <- paths$Beta[paths$Path == "a (WMH->pTau217)"]
  a_p <- paths$P[paths$Path == "a (WMH->pTau217)"]
  text(3.0, 4.3, sprintf("a = %.3f\np = %.4f", a_beta, a_p),
       cex = 0.9, col = "blue")

  # Path b: pTau217 -> MMSE
  arrows(6.5, 5.2, 6.8, 3.5, lwd = 2, col = "blue")
  b_beta <- paths$Beta[paths$Path == "b (pTau217->MMSE)"]
  b_p <- paths$P[paths$Path == "b (pTau217->MMSE)"]
  text(7.5, 4.3, sprintf("b = %.3f\np = %.4f", b_beta, b_p),
       cex = 0.9, col = "blue")

  # Path c': WMH -> MMSE (direct)
  arrows(3.5, 2.8, 6.5, 2.8, lwd = 2, col = "red", lty = 2)
  cp_beta <- paths$Beta[paths$Path == "c' (WMH->MMSE direct)"]
  cp_p <- paths$P[paths$Path == "c' (WMH->MMSE direct)"]
  text(5, 2.2, sprintf("c' = %.3f (p = %.4f)", cp_beta, cp_p),
       cex = 0.9, col = "red")

  # Indirect effect
  text(5, 1.0, sprintf("Indirect (ACME) = %.4f [%.4f, %.4f], p = %.4f",
                        med$ACME, med$ACME_CI_Lower, med$ACME_CI_Upper,
                        med$ACME_P),
       cex = 0.9, font = 2)
  text(5, 0.5, sprintf("Proportion Mediated = %.1f%%, N = %d",
                        med$Prop_Mediated * 100, med$N),
       cex = 0.9)

  title("HABS: WMH -> pTau217 -> MMSE Mediation (5000 Bootstrap, z-scored)",
        cex.main = 1.2)
  dev.off()
}


# --- 11.3 Kaplan-Meier curves ---
plot_km_curves <- function(aibl_results, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(aibl_results$km_fit)) return(NULL)

  p <- ggsurvplot(
    aibl_results$km_fit,
    pval = TRUE,
    risk.table = TRUE,
    xlab = "Time (months)",
    ylab = "AD-Free Survival Probability",
    title = sprintf("AIBL: AD Conversion by APOE4 Status (N=%d)",
                    aibl_results$N),
    palette = c("#0072B2", "#D55E00"),
    legend.labs = c("APOE4 Non-carrier", "APOE4 Carrier"),
    ggtheme = theme_minimal(base_size = 12)
  )

  pdf(file.path(output_dir, "AIBL_KM_curves.pdf"), width = 10, height = 8)
  print(p)
  dev.off()

  return(p)
}


# --- 11.4 DCA plot ---
plot_dca <- function(dca_data, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(dca_data)) return(NULL)

  dca_long <- data.frame(
    Threshold = rep(dca_data$Threshold, 3),
    Net_Benefit = c(dca_data$NB_TreatAll, dca_data$NB_Base, dca_data$NB_Full),
    Model = rep(c("Treat All", "Base Model", "Base + WMH + pTau217"),
                each = nrow(dca_data)),
    stringsAsFactors = FALSE
  )

  p <- ggplot(dca_long, aes(x = Threshold, y = Net_Benefit, color = Model)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Treat All" = "#999999",
                                  "Base Model" = "#0072B2",
                                  "Base + WMH + pTau217" = "#D55E00")) +
    labs(title = "Decision Curve Analysis: HABS Clinical Utility",
         x = "Threshold Probability",
         y = "Net Benefit") +
    coord_cartesian(xlim = c(0.01, 0.80)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, "DCA_plot.pdf"), p, width = 8, height = 6)
  return(p)
}


# --- 11.5 Age-stratified comparison plot ---
plot_age_stratified <- function(age_results, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(age_results$mediation)) return(NULL)

  med_df <- age_results$mediation

  p <- ggplot(med_df, aes(x = Age_Group, y = ACME)) +
    geom_point(size = 4, color = "#0072B2") +
    geom_errorbar(aes(ymin = ACME_CI_Lower, ymax = ACME_CI_Upper),
                  width = 0.2, linewidth = 1, color = "#0072B2") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text(aes(label = sprintf("N=%d\np=%.4f", N, ACME_P)),
              vjust = -1.5, size = 3.5) +
    labs(title = "Age-Stratified Mediation: ACME (Indirect Effect)",
         subtitle = "WMH -> pTau217 -> MMSE (z-scored, 5000 bootstrap)",
         x = "Age Group", y = "ACME (95% CI)") +
    theme_minimal(base_size = 12)

  ggsave(file.path(output_dir, "age_stratified_mediation.pdf"),
         p, width = 7, height = 6)
  return(p)
}


# ==============================================================================
# 12. Q-Q DIAGNOSTIC HELPER
# ==============================================================================

plot_qq_diagnostic <- function(model, title_prefix, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  fname <- gsub("[^A-Za-z0-9_]", "_", title_prefix)
  pdf(file.path(output_dir, paste0(fname, "_qq_plot.pdf")),
      width = 6, height = 6)
  qqnorm(residuals(model),
         main = paste0(title_prefix, ": Normal Q-Q Plot"))
  qqline(residuals(model), col = "red", lwd = 2)
  dev.off()

  # Shapiro-Wilk test (subsample if N > 5000)
  resid <- residuals(model)
  if (length(resid) > 5000) {
    set.seed(2026)
    resid <- sample(resid, 5000)
  }
  sw <- shapiro.test(resid)

  data.frame(
    Test = "Shapiro-Wilk",
    Statistic = sw$statistic,
    P = sw$p.value,
    Normal = sw$p.value > 0.05,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 13. MAIN EXECUTION FUNCTION
# ==============================================================================

run_cross_cohort_validation <- function(
    # A4 files
    a4_subjinfo, a4_ptdemog, a4_sppacc, a4_petsuvr, a4_vmri,
    # HABS files
    habs_genomics_xlsx, habs_biomarkers_csv, habs_wmh_xlsx, habs_clinical_csv,
    # AIBL files
    aibl_ptdemog, aibl_apoeres, aibl_pdxconv,
    # Options
    habs_wmh_sheet = "x1",
    aibl_baseline_year = 2008,
    output_dir = "figures",
    n_boot = 5000,
    save_results = TRUE
) {

  results <- list()
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Step 1: Load data ----
  cat("\n====== Step 1: Loading Data ======\n")

  cat("  Loading A4...\n")
  a4_data <- load_a4_data(a4_subjinfo, a4_ptdemog, a4_sppacc,
                          a4_petsuvr, a4_vmri)

  cat("  Loading HABS...\n")
  habs_data <- load_habs_data(habs_genomics_xlsx, habs_biomarkers_csv,
                              habs_wmh_xlsx, habs_clinical_csv,
                              wmh_sheet = habs_wmh_sheet)

  cat("  Loading AIBL...\n")
  aibl_data <- load_aibl_data(aibl_ptdemog, aibl_apoeres, aibl_pdxconv,
                              baseline_year = aibl_baseline_year)

  results$data <- list(a4 = a4_data, habs = habs_data, aibl = aibl_data)

  # ---- Step 2: A4 WMH-Cognition ----
  cat("\n====== Step 2: A4 WMH-Cognition Association ======\n")
  results$a4_wmh <- run_a4_wmh_cognition(a4_data, output_dir)
  cat(sprintf("  Beta(HC3) = %.4f, P = %.4e, N = %d\n",
              results$a4_wmh$hc3$Beta, results$a4_wmh$hc3$P,
              results$a4_wmh$hc3$N))

  # ---- Step 3: HABS WMH-pTau217 ----
  cat("\n====== Step 3: HABS WMH-pTau217 Association ======\n")
  results$habs_wmh_ptau <- run_habs_wmh_ptau217(habs_data, output_dir)
  cat(sprintf("  Beta(HC3) = %.4f, P = %.4e, N = %d\n",
              results$habs_wmh_ptau$hc3$Beta, results$habs_wmh_ptau$hc3$P,
              results$habs_wmh_ptau$hc3$N))

  # ---- Step 4: HABS Mediation ----
  cat("\n====== Step 4: HABS Mediation (5000 bootstrap, z-scored) ======\n")
  results$habs_mediation <- run_habs_mediation(habs_data, n_boot = n_boot)
  cat(sprintf("  ACME = %.4f [%.4f, %.4f], P = %.4f\n",
              results$habs_mediation$mediation$ACME,
              results$habs_mediation$mediation$ACME_CI_Lower,
              results$habs_mediation$mediation$ACME_CI_Upper,
              results$habs_mediation$mediation$ACME_P))
  cat(sprintf("  Prop Mediated = %.1f%%\n",
              results$habs_mediation$mediation$Prop_Mediated * 100))

  # ---- Step 5: HABS SEM ----
  cat("\n====== Step 5: HABS SEM Confirmation (5000 bootstrap) ======\n")
  results$habs_sem <- run_habs_sem(habs_data, n_boot = n_boot)
  indirect <- results$habs_sem$parameters[
    results$habs_sem$parameters$label == "indirect", ]
  if (nrow(indirect) > 0) {
    cat(sprintf("  SEM Indirect = %.4f [%.4f, %.4f], P = %.4f\n",
                indirect$est, indirect$ci.lower, indirect$ci.upper,
                indirect$pvalue))
  }

  # ---- Step 6: HABS Age-Stratified ----
  cat("\n====== Step 6: HABS Age-Stratified Mediation + SEM ======\n")
  results$habs_age_stratified <- run_habs_age_stratified(habs_data,
                                                          n_boot = n_boot)
  if (!is.null(results$habs_age_stratified$mediation)) {
    for (i in seq_len(nrow(results$habs_age_stratified$mediation))) {
      row <- results$habs_age_stratified$mediation[i, ]
      cat(sprintf("  %s (N=%d): ACME=%.4f, P=%.4f\n",
                  row$Age_Group, row$N, row$ACME, row$ACME_P))
    }
  }

  # ---- Step 7: HABS Age Interaction ----
  cat("\n====== Step 7: HABS Age Interaction Tests ======\n")
  results$habs_age_interaction <- run_habs_age_interaction(habs_data)
  if (!is.null(results$habs_age_interaction)) {
    print(results$habs_age_interaction)
  }

  # ---- Step 8: HABS Clinical Utility ----
  cat("\n====== Step 8: HABS Clinical Utility ======\n")
  results$habs_clinical <- run_habs_clinical_utility(habs_data)
  if (!is.null(results$habs_clinical$auc)) {
    cat(sprintf("  AUC Base: %.3f, AUC Full: %.3f, DeLong P: %.4f\n",
                results$habs_clinical$auc$AUC[1],
                results$habs_clinical$auc$AUC[2],
                results$habs_clinical$delong_p))
  }

  # ---- Step 9: AIBL Survival ----
  cat("\n====== Step 9: AIBL Survival Analysis ======\n")
  results$aibl_survival <- run_aibl_survival(aibl_data)
  cat(sprintf("  N=%d, Events=%d, Median follow-up=%.1f months\n",
              results$aibl_survival$N, results$aibl_survival$N_Events,
              results$aibl_survival$Median_Followup_Months))

  # ---- Step 10: Missing Data Assessment ----
  cat("\n====== Step 10: Missing Data Assessment ======\n")
  results$missing_data <- run_missing_data_assessment(a4_data, habs_data,
                                                      aibl_data)
  print(results$missing_data$summary)

  # ---- Step 11: Visualizations ----
  cat("\n====== Step 11: Generating Visualizations ======\n")

  tryCatch({
    plot_forest(results$a4_wmh, results$habs_wmh_ptau,
                results$aibl_survival, output_dir)
    cat("  Forest plot saved\n")
  }, error = function(e) cat(sprintf("  Forest plot error: %s\n", e$message)))

  tryCatch({
    plot_mediation_path(results$habs_mediation, output_dir)
    cat("  Mediation path diagram saved\n")
  }, error = function(e) cat(sprintf("  Mediation path error: %s\n", e$message)))

  tryCatch({
    plot_km_curves(results$aibl_survival, output_dir)
    cat("  KM curves saved\n")
  }, error = function(e) cat(sprintf("  KM curves error: %s\n", e$message)))

  tryCatch({
    if (!is.null(results$habs_clinical$dca)) {
      plot_dca(results$habs_clinical$dca, output_dir)
      cat("  DCA plot saved\n")
    }
  }, error = function(e) cat(sprintf("  DCA plot error: %s\n", e$message)))

  tryCatch({
    plot_age_stratified(results$habs_age_stratified, output_dir)
    cat("  Age-stratified plot saved\n")
  }, error = function(e) cat(sprintf("  Age-stratified error: %s\n", e$message)))

  # ---- Save results ----
  if (save_results) {
    saveRDS(results, file.path(output_dir, "cross_cohort_validation_results.rds"))
    cat(sprintf("\n  Results saved to: %s\n",
                file.path(output_dir, "cross_cohort_validation_results.rds")))
  }

  # ---- Session info ----
  cat("\n====== Session Info ======\n")
  print(sessionInfo())

  cat("\n====== Cross-Cohort Validation Complete ======\n")
  return(results)
}


# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================
# Uncomment and modify paths to run:
#
# results <- run_cross_cohort_validation(
#   # A4 files
#   a4_subjinfo   = "data/A4/SUBJINFO.csv",
#   a4_ptdemog    = "data/A4/PTDEMOG.csv",
#   a4_sppacc     = "data/A4/SPPACC.csv",
#   a4_petsuvr    = "data/A4/PETSUVR.csv",
#   a4_vmri       = "data/A4/VMRI.csv",
#   # HABS files
#   habs_genomics_xlsx  = "data/HABS/Genomics.xlsx",
#   habs_biomarkers_csv = "data/HABS/Biomarkers.csv",
#   habs_wmh_xlsx       = "data/HABS/WMH.xlsx",
#   habs_clinical_csv   = "data/HABS/Clinical.csv",
#   # AIBL files
#   aibl_ptdemog  = "data/AIBL/ptdemog.csv",
#   aibl_apoeres  = "data/AIBL/apoeres.csv",
#   aibl_pdxconv  = "data/AIBL/pdxconv.csv",
#   # Options
#   habs_wmh_sheet     = "x1",
#   aibl_baseline_year = 2008,
#   output_dir         = "figures",
#   n_boot             = 5000,
#   save_results       = TRUE
# )
#
# # Access individual results:
# results$a4_wmh$hc3                    # A4 HC3 robust results
# results$habs_mediation$mediation      # HABS mediation ACME
# results$habs_sem$parameters           # HABS SEM parameters
# results$habs_age_stratified$mediation # Age-stratified mediation
# results$habs_age_stratified$sem       # Age-stratified SEM
# results$habs_clinical$auc             # Clinical utility AUC
# results$habs_clinical$nri             # NRI with bootstrap CI
# results$habs_clinical$idi             # IDI with bootstrap CI
# results$aibl_survival$cox             # AIBL Cox PH results
# results$missing_data$summary          # Missing data SMD summary
