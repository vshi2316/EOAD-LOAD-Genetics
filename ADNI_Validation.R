# ==============================================================================
# ADNI Independent Cohort Validation of Pathway-Specific Genetic Risk
# ==============================================================================
#
# Purpose: Validate age-dependent oligodendrocyte pathway genetic risk in
#          Alzheimer's disease using the ADNI cohort
#
# Analyses:
#   1.  Data loading and cohort preparation (including CSF sTREM2)
#   2.  APOE-stratified PRS-biomarker associations with dose-response
#   3.  Longitudinal MCI-to-AD conversion (logistic regression)
#   4.  PRS-neuroimaging associations with age interactions
#   5.  Sliding-window age-dependency mapping
#   6.  Unsupervised genetic subtyping (k-means, k=3)
#   7.  Regional cortical thickness analysis (68 DK regions)
#   8.  Subcortical structure analysis (thalamus, putamen, etc.)
#   9.  Cross-sectional PRS-cognition analysis
#   10. Pathway enrichment comparison (microglia vs Abeta clearance)
#   11. Exploratory boundary condition analyses
#   12. Demographics table
#   13. Visualization functions
#   14. Main execution script
#
# Data requirements:
#   - ADNIMERGE clinical data (longitudinal, all visits)
#   - Pathway-specific PRS (LDpred2-derived)
#   - FreeSurfer regional measures (UCSFFSX51)
#   - CSF biomarkers (AlzBio3 / sTREM2)
#
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(broom)
  library(cluster)
  library(sandwich)
  library(lmtest)
  library(cowplot)
  library(viridis)
  library(RColorBrewer)
  library(pheatmap)
})


# ==============================================================================
# UTILITY: Detect longitudinal month and diagnosis columns in ADNIMERGE
# ==============================================================================

detect_month_col <- function(df) {
  # ADNIMERGE versions use M, Month, MONTH, or Month_bl
  candidates <- c("M", "Month", "MONTH", "Month_bl")
  found <- intersect(candidates, colnames(df))
  if (length(found) == 0) {
    # Fallback: look for any column containing "month" (case-insensitive)
    found <- grep("month", colnames(df), ignore.case = TRUE, value = TRUE)
  }
  if (length(found) == 0) stop("Cannot find month column in ADNIMERGE longitudinal data")
  return(found[1])
}

detect_dx_col <- function(df) {
  # Longitudinal diagnosis: DX, DXCHANGE, DX.bl, DX_bl
  if ("DX" %in% colnames(df)) return("DX")
  if ("DXCHANGE" %in% colnames(df)) return("DXCHANGE")
  stop("Cannot find longitudinal diagnosis column (DX or DXCHANGE) in ADNIMERGE")
}

detect_dx_bl_col <- function(df) {
  if ("DX.bl" %in% colnames(df)) return("DX.bl")
  if ("DX_bl" %in% colnames(df)) return("DX_bl")
  stop("Cannot find baseline diagnosis column (DX.bl or DX_bl)")
}

# Map DXCHANGE numeric codes to text labels if needed
map_dx_to_dementia <- function(dx_values, col_name) {
  if (col_name == "DXCHANGE") {
    # DXCHANGE: 4 = MCI to Dementia, 5 = CN to Dementia, 6 = Dementia
    return(dx_values %in% c(4, 5, 6))
  }
  # Text-based DX column
  return(dx_values %in% c("Dementia", "AD", "Alzheimer's Disease"))
}


# ==============================================================================
# 1. DATA LOADING AND INTEGRATION
# ==============================================================================

load_adni_data <- function(clinical_file, prs_file, mri_file,
                           csf_file, strem2_file = NULL) {

  # PRS data
  adni_prs <- fread(prs_file)
  if (!"RID" %in% colnames(adni_prs) && "PTID" %in% colnames(adni_prs)) {
    adni_prs$RID <- as.numeric(gsub(".*_S_0*", "", adni_prs$PTID))
  }

  # Full longitudinal ADNIMERGE (needed for conversion tracking)
  adnimerge <- fread(clinical_file)

  # Baseline clinical data
  clin_bl <- adnimerge[VISCODE == "bl"]
  clin_bl <- clin_bl[!duplicated(clin_bl$RID)]

  # Normalize baseline DX column
  dx_bl <- detect_dx_bl_col(clin_bl)
  if (dx_bl != "DX.bl") clin_bl[, DX.bl := get(dx_bl)]

  merged <- merge(adni_prs, clin_bl, by = "RID", all.x = TRUE)

  # FreeSurfer MRI
  mri_st <- fread(mri_file)
  mri_bl <- mri_st[VISCODE %in% c("bl", "sc")]
  mri_bl <- mri_bl[!duplicated(mri_bl$RID)]
  st_cols <- c("RID", grep("^ST", colnames(mri_bl), value = TRUE))
  merged <- merge(merged, mri_bl[, ..st_cols], by = "RID", all.x = TRUE)

  # CSF biomarkers (AlzBio3)
  csf <- fread(csf_file)
  csf_bl <- csf[!duplicated(csf$RID)]
  csf_merge_cols <- intersect(c("RID", "ABETA", "TAU", "PTAU"), colnames(csf_bl))
  merged <- merge(merged, csf_bl[, ..csf_merge_cols], by = "RID", all.x = TRUE)

  # CSF sTREM2 (separate file if available)
  if (!is.null(strem2_file) && file.exists(strem2_file)) {
    strem2 <- fread(strem2_file)
    strem2_bl <- strem2[!duplicated(strem2$RID)]
    strem2_col <- grep("STREM2|sTREM2|strem2", colnames(strem2_bl), value = TRUE)
    if (length(strem2_col) > 0) {
      strem2_bl$STREM2 <- as.numeric(strem2_bl[[strem2_col[1]]])
      merged <- merge(merged, strem2_bl[, .(RID, STREM2)],
                      by = "RID", all.x = TRUE)
    }
  }

  # Attach longitudinal data for conversion tracking
  attr(merged, "adnimerge_long") <- adnimerge

  return(as.data.frame(merged))
}

create_derived_variables <- function(data) {

  data$Age_Group <- factor(
    ifelse(data$AGE < 70, "Younger", "Older"),
    levels = c("Younger", "Older")
  )
  data$Age_Centered <- data$AGE - 70

  data$APOE_Group <- factor(
    ifelse(data$APOE4 == 2, "e4/e4",
           ifelse(data$APOE4 == 1, "e4 heterozygote", "Non-carrier")),
    levels = c("Non-carrier", "e4 heterozygote", "e4/e4")
  )

  icv <- data$ST10CV
  data$WMH_Log        <- log(data$ST128SV + 1)
  data$Hippo_Norm      <- (data$ST29SV + data$ST88SV) / icv * 1000
  data$Ventricle_Norm  <- (data$ST37SV + data$ST96SV) / icv * 1000
  data$Entorhinal_Norm <- (data$ST24CV + data$ST83CV) / icv * 1000

  for (col in c("ABETA", "TAU", "PTAU")) {
    if (col %in% colnames(data)) {
      data[[col]] <- as.numeric(gsub("[^0-9.]", "", data[[col]]))
    }
  }

  # Longitudinal MCI-to-AD conversion from multi-visit data
  adnimerge_long <- attr(data, "adnimerge_long")
  data$Converted <- 0L
  data$Time_to_Conversion <- NA_real_

  if (!is.null(adnimerge_long)) {
    adnimerge_long <- as.data.frame(adnimerge_long)

    # Auto-detect column names
    month_col <- detect_month_col(adnimerge_long)
    dx_long_col <- detect_dx_col(adnimerge_long)
    dx_bl <- detect_dx_bl_col(adnimerge_long)
    if (dx_bl != "DX.bl") adnimerge_long$DX.bl <- adnimerge_long[[dx_bl]]

    mci_rids <- data$RID[data$DX.bl %in% c("LMCI", "EMCI", "MCI")]

    for (rid in mci_rids) {
      visits <- adnimerge_long[adnimerge_long$RID == rid, ]
      visits <- visits[order(visits[[month_col]]), ]

      # Find first visit where DX indicates dementia/AD
      is_dem <- map_dx_to_dementia(visits[[dx_long_col]], dx_long_col)
      dem_visits <- visits[is_dem, ]
      if (nrow(dem_visits) > 0) {
        first_dem <- dem_visits[1, ]
        idx <- which(data$RID == rid)
        data$Converted[idx] <- 1L
        conv_month <- first_dem[[month_col]]
        if (!is.na(conv_month)) {
          data$Time_to_Conversion[idx] <- as.numeric(conv_month) / 12
        }
      }
    }
  }

  # Median follow-up for MCI subjects
  if (!is.null(adnimerge_long)) {
    month_col <- detect_month_col(adnimerge_long)
    mci_data <- data[data$DX.bl %in% c("LMCI", "EMCI", "MCI"), ]
    follow_up_months <- sapply(mci_data$RID, function(rid) {
      visits <- adnimerge_long[adnimerge_long$RID == rid, ]
      if (nrow(visits) > 0) max(as.numeric(visits[[month_col]]), na.rm = TRUE) else NA
    })
    data$Median_Followup_Years <- median(follow_up_months / 12, na.rm = TRUE)
  }

  return(data)
}


# ==============================================================================
# 2. APOE-STRATIFIED PRS-BIOMARKER ASSOCIATIONS WITH DOSE-RESPONSE
# ==============================================================================

run_apoe_stratified_analysis <- function(data) {

  prs_vars <- c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin", "PRS_EOAD_Microglia",
                 "PRS_EOAD_Abeta", "PRS_EOAD_APP", "PRS_EOAD_Global",
                 "PRS_EOAD_Lipid", "PRS_EOAD_TCell")
  outcomes <- c("ABETA", "TAU", "PTAU", "WMH_Log", "Hippo_Norm")
  if ("STREM2" %in% colnames(data)) outcomes <- c(outcomes, "STREM2")

  results <- list()
  for (apoe_grp in levels(data$APOE_Group)) {
    sub <- data[data$APOE_Group == apoe_grp, ]
    if (nrow(sub) < 20) next

    for (prs in prs_vars) {
      if (!prs %in% colnames(data)) next
      for (outcome in outcomes) {
        if (!outcome %in% colnames(data)) next
        covars <- c("AGE", "PTGENDER", "PTEDUCAT")
        df <- sub[complete.cases(sub[, c(prs, outcome, covars)]), ]
        if (nrow(df) < 15) next

        fml <- as.formula(paste(outcome, "~", prs, "+",
                                paste(covars, collapse = " + ")))
        fit <- tryCatch(lm(fml, data = df), error = function(e) NULL)
        if (is.null(fit)) next
        s <- summary(fit)$coefficients
        if (!prs %in% rownames(s)) next

        results[[length(results) + 1]] <- data.frame(
          APOE_Group = apoe_grp, PRS = prs, Outcome = outcome,
          Beta = s[prs, "Estimate"], SE = s[prs, "Std. Error"],
          P = s[prs, "Pr(>|t|)"], N = nrow(df),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  res <- do.call(rbind, results)
  res$P_FDR <- p.adjust(res$P, method = "BH")
  return(res)
}


# ==============================================================================
# 2b. DOSE-RESPONSE: PRS QUARTILE TREND ANALYSIS (sTREM2 focus)
# ==============================================================================

run_dose_response_analysis <- function(data, n_quantiles = 4) {

  prs_vars <- c("PRS_EOAD_Microglia", "PRS_EOAD_Oligo", "PRS_EOAD_Myelin",
                 "PRS_EOAD_Abeta")
  prs_vars <- prs_vars[prs_vars %in% colnames(data)]

  outcomes <- c("ABETA", "TAU", "PTAU", "WMH_Log", "Hippo_Norm")
  if ("STREM2" %in% colnames(data)) outcomes <- c(outcomes, "STREM2")

  covars <- c("AGE", "PTGENDER", "PTEDUCAT")

  trend_results <- list()
  quantile_means <- list()

  for (apoe_grp in levels(data$APOE_Group)) {
    sub <- data[data$APOE_Group == apoe_grp, ]
    if (nrow(sub) < 20) next

    for (prs in prs_vars) {
      if (!prs %in% colnames(sub)) next
      df <- sub[!is.na(sub[[prs]]), ]
      if (nrow(df) < n_quantiles * 5) next

      # Create PRS quantile groups
      df$PRS_Quantile <- ntile(df[[prs]], n_quantiles)

      for (outcome in outcomes) {
        if (!outcome %in% colnames(df)) next
        df_out <- df[complete.cases(df[, c("PRS_Quantile", outcome, covars)]), ]
        if (nrow(df_out) < n_quantiles * 3) next

        # Quantile-level means
        q_means <- df_out %>%
          group_by(PRS_Quantile) %>%
          summarise(
            Mean = mean(.data[[outcome]], na.rm = TRUE),
            SD = sd(.data[[outcome]], na.rm = TRUE),
            N = n(),
            .groups = "drop"
          ) %>%
          mutate(APOE_Group = apoe_grp, PRS = prs, Outcome = outcome)
        quantile_means[[length(quantile_means) + 1]] <- q_means

        # Linear trend test: PRS_Quantile as ordered numeric predictor
        fml <- as.formula(paste(outcome, "~ PRS_Quantile +",
                                paste(covars, collapse = " + ")))
        fit <- tryCatch(lm(fml, data = df_out), error = function(e) NULL)
        if (is.null(fit)) next
        s <- summary(fit)$coefficients
        if (!"PRS_Quantile" %in% rownames(s)) next

        trend_results[[length(trend_results) + 1]] <- data.frame(
          APOE_Group = apoe_grp, PRS = prs, Outcome = outcome,
          Beta_Trend = s["PRS_Quantile", "Estimate"],
          SE = s["PRS_Quantile", "Std. Error"],
          P_Trend = s["PRS_Quantile", "Pr(>|t|)"],
          N = nrow(df_out),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  trend_df <- if (length(trend_results) > 0) do.call(rbind, trend_results) else data.frame()
  means_df <- if (length(quantile_means) > 0) do.call(rbind, quantile_means) else data.frame()

  if (nrow(trend_df) > 0) trend_df$P_FDR <- p.adjust(trend_df$P_Trend, method = "BH")

  list(trend = trend_df, quantile_means = means_df)
}


# ==============================================================================
# 3. LONGITUDINAL MCI-TO-AD CONVERSION (LOGISTIC REGRESSION)
# ==============================================================================

run_conversion_analysis <- function(data) {

  dx_bl_col <- detect_dx_bl_col(data)
  mci <- data[data[[dx_bl_col]] %in% c("LMCI", "EMCI", "MCI"), ]

  prs_vars <- c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin", "PRS_EOAD_Microglia",
                 "PRS_EOAD_Abeta", "PRS_EOAD_Global",
                 "PRS_EOAD_Lipid", "PRS_EOAD_TCell")
  prs_vars <- prs_vars[prs_vars %in% colnames(mci)]
  covars <- c("AGE", "PTGENDER", "APOE4", "PTEDUCAT")

  # --- Full-cohort logistic regression ---
  full_results <- list()
  for (prs in prs_vars) {
    df <- mci[complete.cases(mci[, c(prs, "Converted", covars)]), ]
    if (nrow(df) < 20 || sum(df$Converted) < 5) next
    fml <- as.formula(paste("Converted ~", prs, "+",
                            paste(covars, collapse = " + ")))
    fit <- glm(fml, data = df, family = binomial)
    ci <- confint.default(fit)
    s <- summary(fit)$coefficients[prs, ]

    full_results[[length(full_results) + 1]] <- data.frame(
      Stratum = "Full_MCI", PRS = prs,
      OR = exp(s["Estimate"]),
      OR_Lower = exp(ci[prs, 1]), OR_Upper = exp(ci[prs, 2]),
      P = s["Pr(>|z|)"], N = nrow(df),
      N_Converted = sum(df$Converted, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  # --- Age-stratified logistic regression ---
  strat_results <- list()
  for (age_grp in c("Younger", "Older")) {
    sub <- mci[mci$Age_Group == age_grp, ]
    for (prs in prs_vars) {
      df <- sub[complete.cases(sub[, c(prs, "Converted", covars)]), ]
      if (nrow(df) < 20 || sum(df$Converted) < 5) next

      fml <- as.formula(paste("Converted ~", prs, "+",
                              paste(covars, collapse = " + ")))
      fit <- glm(fml, data = df, family = binomial)
      ci <- confint.default(fit)
      s <- summary(fit)$coefficients[prs, ]

      strat_results[[length(strat_results) + 1]] <- data.frame(
        Stratum = age_grp, PRS = prs,
        OR = exp(s["Estimate"]),
        OR_Lower = exp(ci[prs, 1]), OR_Upper = exp(ci[prs, 2]),
        P = s["Pr(>|z|)"], N = nrow(df),
        N_Converted = sum(df$Converted, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Age x PRS interaction ---
  interaction_results <- list()
  for (prs in prs_vars) {
    df <- mci[complete.cases(mci[, c(prs, "Converted", covars, "Age_Group")]), ]
    if (nrow(df) < 30) next
    fml <- as.formula(paste("Converted ~", prs, "* Age_Group +",
                            paste(covars, collapse = " + ")))
    fit <- glm(fml, data = df, family = binomial)
    int_term <- grep(paste0(prs, ":"), names(coef(fit)), value = TRUE)
    if (length(int_term) > 0) {
      ci <- confint.default(fit)
      s <- summary(fit)$coefficients[int_term, ]
      interaction_results[[length(interaction_results) + 1]] <- data.frame(
        PRS = prs, Interaction_Term = int_term,
        OR_Interaction = exp(s["Estimate"]),
        OR_Lower = exp(ci[int_term, 1]), OR_Upper = exp(ci[int_term, 2]),
        P = s["Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Sensitivity: <65 years ---
  sens_results <- list()
  mci_young <- mci[mci$AGE < 65, ]
  for (prs in c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin")) {
    if (!prs %in% colnames(mci_young)) next
    df <- mci_young[complete.cases(mci_young[, c(prs, "Converted", covars)]), ]
    if (nrow(df) < 15 || sum(df$Converted) < 3) next
    fml <- as.formula(paste("Converted ~", prs, "+",
                            paste(covars, collapse = " + ")))
    fit <- glm(fml, data = df, family = binomial)
    ci <- confint.default(fit)
    s <- summary(fit)$coefficients[prs, ]
    sens_results[[length(sens_results) + 1]] <- data.frame(
      Stratum = "Under65", PRS = prs,
      OR = exp(s["Estimate"]),
      OR_Lower = exp(ci[prs, 1]), OR_Upper = exp(ci[prs, 2]),
      P = s["Pr(>|z|)"], N = nrow(df),
      stringsAsFactors = FALSE
    )
  }

  list(
    full = do.call(rbind, full_results),
    stratified = do.call(rbind, strat_results),
    interaction = do.call(rbind, interaction_results),
    sensitivity_under65 = if (length(sens_results) > 0) do.call(rbind, sens_results) else NULL
  )
}


# ==============================================================================
# 4. PRS-NEUROIMAGING ASSOCIATIONS WITH AGE INTERACTIONS
# ==============================================================================

run_neuroimaging_analysis <- function(data) {

  prs_vars <- c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin", "PRS_EOAD_Microglia",
                 "PRS_EOAD_Abeta", "PRS_EOAD_Global")
  prs_vars <- prs_vars[prs_vars %in% colnames(data)]
  outcomes <- c("WMH_Log", "Hippo_Norm", "Ventricle_Norm", "Entorhinal_Norm")
  covars <- c("AGE", "PTGENDER", "APOE4", "PTEDUCAT")

  # --- Main effects by age stratum ---
  strat_results <- list()
  for (age_grp in c("Younger", "Older")) {
    sub <- data[data$Age_Group == age_grp, ]
    for (prs in prs_vars) {
      for (outcome in outcomes) {
        df <- sub[complete.cases(sub[, c(prs, outcome, covars)]), ]
        if (nrow(df) < 20) next

        fml <- as.formula(paste(outcome, "~", prs, "+",
                                paste(covars, collapse = " + ")))
        fit <- lm(fml, data = df)
        s <- summary(fit)$coefficients[prs, ]

        strat_results[[length(strat_results) + 1]] <- data.frame(
          Age_Group = age_grp, PRS = prs, Outcome = outcome,
          Beta = s["Estimate"], SE = s["Std. Error"],
          P = s["Pr(>|t|)"], N = nrow(df),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  strat_df <- do.call(rbind, strat_results)
  strat_df$P_FDR <- p.adjust(strat_df$P, method = "BH")

  # --- Age x PRS interaction (continuous age) ---
  interaction_results <- list()
  for (prs in prs_vars) {
    for (outcome in outcomes) {
      df <- data[complete.cases(data[, c(prs, outcome, covars)]), ]
      if (nrow(df) < 50) next

      fml <- as.formula(paste(outcome, "~", prs, "* AGE + PTGENDER + APOE4 + PTEDUCAT"))
      fit <- lm(fml, data = df)
      int_term <- paste0(prs, ":AGE")
      if (int_term %in% rownames(summary(fit)$coefficients)) {
        s <- summary(fit)$coefficients[int_term, ]
        interaction_results[[length(interaction_results) + 1]] <- data.frame(
          PRS = prs, Outcome = outcome,
          Beta_Interaction = s["Estimate"], SE = s["Std. Error"],
          P_Interaction = s["Pr(>|t|)"], N = nrow(df),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  int_df <- do.call(rbind, interaction_results)
  if (!is.null(int_df) && nrow(int_df) > 0) {
    int_df$P_FDR <- p.adjust(int_df$P_Interaction, method = "BH")
  }

  list(stratified = strat_df, interaction = int_df)
}


# ==============================================================================
# 5. SLIDING-WINDOW AGE-DEPENDENCY MAPPING
# ==============================================================================

run_sliding_window <- function(data, window_width = 10, step = 2,
                               min_age = 55, max_age = 85, min_n = 30) {

  centre_ages <- seq(min_age + window_width / 2,
                     max_age - window_width / 2, by = step)
  prs_var <- "PRS_EOAD_Oligo"

  dx_bl_col <- detect_dx_bl_col(data)
  mci <- data[data[[dx_bl_col]] %in% c("LMCI", "EMCI", "MCI"), ]

  wmh_results <- list()
  conv_results <- list()

  for (ca in centre_ages) {
    lo <- ca - window_width / 2
    hi <- ca + window_width / 2

    # WMH linear regression (full sample, ICV as covariate)
    w_all <- data[data$AGE >= lo & data$AGE < hi, ]
    df_wmh <- w_all[complete.cases(w_all[, c(prs_var, "WMH_Log", "AGE",
                                              "PTGENDER", "ST10CV")]), ]
    if (nrow(df_wmh) >= min_n) {
      fml <- as.formula(paste("WMH_Log ~", prs_var, "+ AGE + PTGENDER + ST10CV"))
      fit <- lm(fml, data = df_wmh)
      s <- summary(fit)$coefficients[prs_var, ]
      wmh_results[[length(wmh_results) + 1]] <- data.frame(
        Centre_Age = ca, Beta = s["Estimate"], SE = s["Std. Error"],
        P = s["Pr(>|t|)"], N = nrow(df_wmh),
        stringsAsFactors = FALSE
      )
    }

    # Conversion logistic regression (MCI only)
    w_mci <- mci[mci$AGE >= lo & mci$AGE < hi, ]
    conv_covars <- c(prs_var, "Converted", "AGE", "PTGENDER", "APOE4", "PTEDUCAT")
    df_conv <- w_mci[complete.cases(w_mci[, conv_covars]), ]
    if (nrow(df_conv) >= min_n && sum(df_conv$Converted) >= 5) {
      fml <- as.formula(paste("Converted ~", prs_var,
                              "+ AGE + PTGENDER + APOE4 + PTEDUCAT"))
      fit <- glm(fml, data = df_conv, family = binomial)
      ci <- confint.default(fit)
      s <- summary(fit)$coefficients[prs_var, ]
      conv_results[[length(conv_results) + 1]] <- data.frame(
        Centre_Age = ca,
        OR = exp(s["Estimate"]),
        OR_Lower = exp(ci[prs_var, 1]), OR_Upper = exp(ci[prs_var, 2]),
        P = s["Pr(>|z|)"], N = nrow(df_conv),
        N_Converted = sum(df_conv$Converted),
        stringsAsFactors = FALSE
      )
    }
  }

  wmh_df <- if (length(wmh_results) > 0) do.call(rbind, wmh_results) else data.frame()
  conv_df <- if (length(conv_results) > 0) do.call(rbind, conv_results) else data.frame()

  if (nrow(wmh_df) > 0) {
    wmh_df$P_FDR <- p.adjust(wmh_df$P, method = "BH")
    # Peak = window with largest positive beta (strongest risk effect)
    peak_beta <- max(abs(wmh_df$Beta))
    late_beta <- abs(wmh_df$Beta[wmh_df$Centre_Age == max(wmh_df$Centre_Age)])
    wmh_df$Dissipation_Pct <- round((1 - late_beta / peak_beta) * 100, 1)
    wmh_df$Peak_Centre_Age <- wmh_df$Centre_Age[which.max(abs(wmh_df$Beta))]
  }

  if (nrow(conv_df) > 0) {
    conv_df$P_FDR <- p.adjust(conv_df$P, method = "BH")
    # Peak = window with highest OR (strongest risk for conversion)
    peak_idx <- which.max(conv_df$OR)
    peak_or <- conv_df$OR[peak_idx]
    late_or <- conv_df$OR[conv_df$Centre_Age == max(conv_df$Centre_Age)]
    # Dissipation based on log(OR) reduction from peak
    conv_df$Dissipation_Pct <- round(
      (1 - log(late_or) / log(peak_or)) * 100, 1
    )
    conv_df$Peak_Centre_Age <- conv_df$Centre_Age[peak_idx]
  }

  list(wmh = wmh_df, conversion = conv_df)
}


# ==============================================================================
# 6. UNSUPERVISED GENETIC SUBTYPING (K-MEANS, K=3)
# ==============================================================================

run_genetic_subtyping <- function(data) {

  prs_vars <- c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin", "PRS_EOAD_Microglia",
                 "PRS_EOAD_Abeta", "PRS_EOAD_APP", "PRS_EOAD_Global")
  prs_vars <- prs_vars[prs_vars %in% colnames(data)]

  df <- data[complete.cases(data[, prs_vars]), ]
  prs_mat <- scale(as.matrix(df[, prs_vars]))

  set.seed(42)
  km <- kmeans(prs_mat, centers = 3, nstart = 25, iter.max = 100)
  df$Cluster_KM <- km$cluster

  hc <- hclust(dist(prs_mat), method = "ward.D2")
  df$Cluster_HC <- cutree(hc, k = 3)

  tab <- table(df$Cluster_KM, df$Cluster_HC)
  concordance <- sum(apply(tab, 1, max)) / nrow(df)

  # Cohen's kappa for clustering agreement
  kappa_value <- NA
  kappa_p <- NA
  if (requireNamespace("irr", quietly = TRUE)) {
    kappa_result <- irr::kappa2(data.frame(km = df$Cluster_KM, hc = df$Cluster_HC))
    kappa_value <- kappa_result$value
    kappa_p <- kappa_result$p.value
  }

  sil <- silhouette(km$cluster, dist(prs_mat))
  avg_sil <- mean(sil[, "sil_width"])

  # Label clusters by PRS profile
  centres <- as.data.frame(km$centers)
  centres$Cluster <- 1:3
  oligo_myelin_score <- centres$PRS_EOAD_Oligo + centres$PRS_EOAD_Myelin
  abeta_score <- centres$PRS_EOAD_Abeta

  labels <- rep("Background_Risk", 3)
  labels[which.max(oligo_myelin_score)] <- "Oligo_Driven"
  labels[which.max(abeta_score)] <- "High_Abeta"
  if (sum(labels == "Background_Risk") == 0) {
    bg_idx <- which.min(oligo_myelin_score + abeta_score)
    labels[bg_idx] <- "Background_Risk"
  }

  label_map <- setNames(labels, 1:3)
  df$Subtype <- factor(label_map[as.character(df$Cluster_KM)],
                        levels = c("Background_Risk", "Oligo_Driven", "High_Abeta"))

  # Pathway burden as percentage of total variance
  pathway_burden <- df %>%
    group_by(Subtype) %>%
    summarise(
      N = n(),
      across(all_of(prs_vars), ~ var(.x, na.rm = TRUE), .names = "Var_{.col}"),
      .groups = "drop"
    ) %>%
    mutate(
      Total_Var = rowSums(select(., starts_with("Var_"))),
      Myelin_Pct = Var_PRS_EOAD_Myelin / Total_Var * 100,
      Oligo_Pct = Var_PRS_EOAD_Oligo / Total_Var * 100
    )

  # ANOVA across subtypes for key outcomes
  anova_vars <- c("WMH_Log", "Ventricle_Norm", "AGE")
  subtype_anova <- do.call(rbind, lapply(anova_vars, function(v) {
    df_v <- df[!is.na(df[[v]]), ]
    if (nrow(df_v) < 30) return(NULL)
    aov_fit <- aov(as.formula(paste(v, "~ Subtype")), data = df_v)
    data.frame(
      Variable = v,
      F_stat = summary(aov_fit)[[1]]["Subtype", "F value"],
      P = summary(aov_fit)[[1]]["Subtype", "Pr(>F)"],
      N = nrow(df_v)
    )
  }))

  # Descriptive statistics by subtype
  subtype_summary <- df %>%
    group_by(Subtype) %>%
    summarise(
      N = n(),
      Age_Mean = mean(AGE, na.rm = TRUE),
      Age_SD = sd(AGE, na.rm = TRUE),
      WMH_Mean = mean(WMH_Log, na.rm = TRUE),
      WMH_SD = sd(WMH_Log, na.rm = TRUE),
      Ventricle_Mean = mean(Ventricle_Norm, na.rm = TRUE),
      Ventricle_SD = sd(Ventricle_Norm, na.rm = TRUE),
      Hippo_Mean = mean(Hippo_Norm, na.rm = TRUE),
      Hippo_SD = sd(Hippo_Norm, na.rm = TRUE),
      Hippo_Vent_Ratio_Mean = mean(Hippo_Norm / Ventricle_Norm, na.rm = TRUE),
      Hippo_Vent_Ratio_SD = sd(Hippo_Norm / Ventricle_Norm, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    data = df,
    centres = centres,
    label_map = label_map,
    concordance = concordance,
    kappa = kappa_value,
    kappa_p = kappa_p,
    silhouette = avg_sil,
    pathway_burden = pathway_burden,
    subtype_anova = subtype_anova,
    subtype_summary = subtype_summary
  )
}


# ==============================================================================
# 7. REGIONAL CORTICAL THICKNESS ANALYSIS (68 DK REGIONS)
# ==============================================================================
#
# ST codes from UCSFFSX51 data dictionary (FreeSurfer cross-sectional v5.1)
# Cortical thickness uses TA suffix. Left hemisphere: ST12-ST56, Right: ST71-ST115
# Reference: ADNI FreeSurfer Methods, UCSFFSX51_DICT.csv
#

get_dk_regions <- function() {
  # Complete Desikan-Killiany atlas: 34 bilateral cortical regions (68 total)
  # Each entry: c(left_hemisphere_ST_code, right_hemisphere_ST_code)
  list(
    bankssts                 = c("ST14TA", "ST73TA"),
    caudalanteriorcingulate   = c("ST15TA", "ST74TA"),
    caudalmiddlefrontal      = c("ST16TA", "ST75TA"),
    cuneus                   = c("ST18TA", "ST77TA"),
    entorhinal               = c("ST24TA", "ST83TA"),
    fusiform                 = c("ST25TA", "ST84TA"),
    inferiorparietal         = c("ST26TA", "ST85TA"),
    inferiortemporal         = c("ST27TA", "ST86TA"),
    isthmuscingulate         = c("ST28TA", "ST87TA"),
    lateraloccipital         = c("ST30TA", "ST89TA"),
    lateralorbitofrontal     = c("ST31TA", "ST90TA"),
    lingual                  = c("ST32TA", "ST91TA"),
    medialorbitofrontal      = c("ST33TA", "ST92TA"),
    middletemporal           = c("ST34TA", "ST93TA"),
    paracentral              = c("ST35TA", "ST94TA"),
    parsopercularis          = c("ST36TA", "ST95TA"),
    parsorbitalis            = c("ST38TA", "ST97TA"),
    parstriangularis         = c("ST39TA", "ST98TA"),
    pericalcarine            = c("ST40TA", "ST99TA"),
    postcentral              = c("ST41TA", "ST100TA"),
    posteriorcingulate       = c("ST42TA", "ST101TA"),
    precentral               = c("ST43TA", "ST102TA"),
    precuneus                = c("ST44TA", "ST103TA"),
    rostralanteriorcingulate = c("ST45TA", "ST104TA"),
    rostralmiddlefrontal     = c("ST46TA", "ST105TA"),
    superiorfrontal          = c("ST47TA", "ST106TA"),
    superiorparietal         = c("ST48TA", "ST107TA"),
    superiortemporal         = c("ST49TA", "ST108TA"),
    supramarginal            = c("ST50TA", "ST109TA"),
    frontalpole              = c("ST23TA", "ST82TA"),
    temporalpole             = c("ST52TA", "ST111TA"),
    transversetemporal       = c("ST53TA", "ST112TA"),
    insula                   = c("ST56TA", "ST115TA"),
    parahippocampal          = c("ST37TA", "ST96TA")
  )
}

run_brain_region_analysis <- function(data) {

  dk_regions <- get_dk_regions()
  ref_subtype <- "Background_Risk"

  results <- list()
  for (region in names(dk_regions)) {
    codes <- dk_regions[[region]]
    for (i in seq_along(codes)) {
      hemi <- c("lh", "rh")[i]
      st_code <- codes[i]
      if (!st_code %in% colnames(data)) next

      df <- data.frame(
        Thickness = data[[st_code]],
        Subtype = data$Subtype,
        AGE = data$AGE,
        PTGENDER = data$PTGENDER,
        APOE4 = data$APOE4,
        PTEDUCAT = data$PTEDUCAT,
        stringsAsFactors = FALSE
      )
      df <- df[complete.cases(df), ]
      if (nrow(df) < 30) next

      fit <- lm(Thickness ~ Subtype + AGE + PTGENDER + APOE4 + PTEDUCAT, data = df)
      s <- summary(fit)$coefficients

      for (sub in c("Oligo_Driven", "High_Abeta")) {
        term <- paste0("Subtype", sub)
        if (!term %in% rownames(s)) next
        coefs <- s[term, ]
        results[[length(results) + 1]] <- data.frame(
          Region = region, Hemisphere = hemi, ST_Code = st_code,
          Subtype = sub, Reference = ref_subtype,
          Beta = coefs["Estimate"], SE = coefs["Std. Error"],
          t_stat = coefs["t value"], P = coefs["Pr(>|t|)"],
          N = nrow(df),
          Cohens_d = coefs["t value"] / sqrt(nrow(df)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  res <- do.call(rbind, results)
  if (!is.null(res) && nrow(res) > 0) {
    res$P_FDR <- p.adjust(res$P, method = "BH")
  }
  return(res)
}

run_regional_age_prs_interaction <- function(data) {

  dk_regions <- get_dk_regions()
  prs_var <- "PRS_EOAD_Oligo"

  results <- list()
  for (region in names(dk_regions)) {
    codes <- dk_regions[[region]]
    for (i in seq_along(codes)) {
      hemi <- c("lh", "rh")[i]
      st_code <- codes[i]
      if (!st_code %in% colnames(data)) next

      df <- data.frame(
        Thickness = data[[st_code]],
        PRS = data[[prs_var]],
        AGE = data$AGE,
        PTGENDER = data$PTGENDER,
        APOE4 = data$APOE4,
        PTEDUCAT = data$PTEDUCAT,
        stringsAsFactors = FALSE
      )
      df <- df[complete.cases(df), ]
      if (nrow(df) < 50) next

      fit <- lm(Thickness ~ PRS * AGE + PTGENDER + APOE4 + PTEDUCAT, data = df)
      int_term <- "PRS:AGE"
      if (int_term %in% rownames(summary(fit)$coefficients)) {
        s <- summary(fit)$coefficients[int_term, ]
        results[[length(results) + 1]] <- data.frame(
          Region = region, Hemisphere = hemi, ST_Code = st_code,
          Beta_Interaction = s["Estimate"], SE = s["Std. Error"],
          P_Interaction = s["Pr(>|t|)"], N = nrow(df),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  res <- do.call(rbind, results)
  if (!is.null(res) && nrow(res) > 0) {
    res$P_FDR <- p.adjust(res$P_Interaction, method = "BH")
  }
  return(res)
}


# ==============================================================================
# 8. SUBCORTICAL STRUCTURE ANALYSIS
# ==============================================================================
#
# ST codes from UCSFFSX51 data dictionary for subcortical volumes (SV suffix)
# Verified against UCSFFSX51_DICT.csv:
#   ST29SV = Left-Hippocampus, ST88SV = Right-Hippocampus
#   ST13SV = Left-Amygdala,    ST72SV = Right-Amygdala
#   ST19SV = Left-Caudate,     ST78SV = Right-Caudate
#   ST36SV = Left-Putamen,     ST95SV = Right-Putamen
#   ST34SV = Left-Pallidum,    ST93SV = Right-Pallidum
#   ST54SV = Left-Thalamus,    ST113SV = Right-Thalamus
#   ST12SV = Left-Accumbens,   ST71SV = Right-Accumbens
#

get_subcortical_regions <- function() {
  list(
    hippocampus      = c("ST29SV", "ST88SV"),
    amygdala         = c("ST13SV", "ST72SV"),
    caudate          = c("ST19SV", "ST78SV"),
    putamen          = c("ST36SV", "ST95SV"),
    pallidum         = c("ST34SV", "ST93SV"),
    thalamus         = c("ST54SV", "ST113SV"),
    accumbens        = c("ST12SV", "ST71SV")
  )
}

run_subcortical_analysis <- function(data) {

  subcort_regions <- get_subcortical_regions()
  icv <- data$ST10CV
  prs_var <- "PRS_EOAD_Oligo"
  covars <- c("AGE", "PTGENDER", "APOE4", "PTEDUCAT")

  # --- Subtype-based analysis ---
  subtype_results <- list()
  for (region in names(subcort_regions)) {
    codes <- subcort_regions[[region]]
    for (i in seq_along(codes)) {
      hemi <- c("lh", "rh")[i]
      st_code <- codes[i]
      if (!st_code %in% colnames(data)) next

      vol_norm <- data[[st_code]] / icv * 1000
      df <- data.frame(
        Vol = vol_norm,
        Subtype = data$Subtype,
        AGE = data$AGE,
        PTGENDER = data$PTGENDER,
        APOE4 = data$APOE4,
        PTEDUCAT = data$PTEDUCAT,
        stringsAsFactors = FALSE
      )
      df <- df[complete.cases(df), ]
      if (nrow(df) < 30) next

      fit <- lm(Vol ~ Subtype + AGE + PTGENDER + APOE4 + PTEDUCAT, data = df)
      s <- summary(fit)$coefficients

      for (sub in c("Oligo_Driven", "High_Abeta")) {
        term <- paste0("Subtype", sub)
        if (!term %in% rownames(s)) next
        coefs <- s[term, ]
        subtype_results[[length(subtype_results) + 1]] <- data.frame(
          Region = region, Hemisphere = hemi, ST_Code = st_code,
          Subtype = sub, Reference = "Background_Risk",
          Beta = coefs["Estimate"], SE = coefs["Std. Error"],
          t_stat = coefs["t value"], P = coefs["Pr(>|t|)"],
          N = nrow(df),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # --- Age x PRS interaction ---
  interaction_results <- list()
  for (region in names(subcort_regions)) {
    codes <- subcort_regions[[region]]
    for (i in seq_along(codes)) {
      hemi <- c("lh", "rh")[i]
      st_code <- codes[i]
      if (!st_code %in% colnames(data)) next

      vol_norm <- data[[st_code]] / icv * 1000
      df <- data.frame(
        Vol = vol_norm,
        PRS = data[[prs_var]],
        AGE = data$AGE,
        PTGENDER = data$PTGENDER,
        APOE4 = data$APOE4,
        PTEDUCAT = data$PTEDUCAT,
        stringsAsFactors = FALSE
      )
      df <- df[complete.cases(df), ]
      if (nrow(df) < 50) next

      fit <- lm(Vol ~ PRS * AGE + PTGENDER + APOE4 + PTEDUCAT, data = df)
      int_term <- "PRS:AGE"
      if (int_term %in% rownames(summary(fit)$coefficients)) {
        s <- summary(fit)$coefficients[int_term, ]
        interaction_results[[length(interaction_results) + 1]] <- data.frame(
          Region = region, Hemisphere = hemi, ST_Code = st_code,
          Beta_Interaction = s["Estimate"], SE = s["Std. Error"],
          P_Interaction = s["Pr(>|t|)"], N = nrow(df),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  sub_df <- if (length(subtype_results) > 0) do.call(rbind, subtype_results) else data.frame()
  int_df <- if (length(interaction_results) > 0) do.call(rbind, interaction_results) else data.frame()

  if (nrow(sub_df) > 0) sub_df$P_FDR <- p.adjust(sub_df$P, method = "BH")
  if (nrow(int_df) > 0) int_df$P_FDR <- p.adjust(int_df$P_Interaction, method = "BH")

  list(subtype = sub_df, interaction = int_df)
}


# ==============================================================================
# 9. CROSS-SECTIONAL PRS-COGNITION ANALYSIS
# ==============================================================================

run_cross_sectional_cognition <- function(data) {

  prs_vars <- c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin", "PRS_EOAD_Microglia",
                 "PRS_EOAD_Abeta", "PRS_EOAD_APP", "PRS_EOAD_Global",
                 "PRS_EOAD_Lipid", "PRS_EOAD_TCell")
  prs_vars <- prs_vars[prs_vars %in% colnames(data)]
  outcomes <- c("MMSE", "CDRSB", "ADAS13", "RAVLT_immediate")
  outcomes <- outcomes[outcomes %in% colnames(data)]
  covars <- c("AGE", "PTGENDER", "APOE4", "PTEDUCAT")

  results <- list()
  for (prs in prs_vars) {
    for (outcome in outcomes) {
      df <- data[complete.cases(data[, c(prs, outcome, covars)]), ]
      if (nrow(df) < 30) next

      fml <- as.formula(paste(outcome, "~", prs, "+",
                              paste(covars, collapse = " + ")))
      fit <- lm(fml, data = df)
      s <- summary(fit)$coefficients
      if (!prs %in% rownames(s)) next

      results[[length(results) + 1]] <- data.frame(
        PRS = prs, Outcome = outcome,
        Beta = s[prs, "Estimate"], SE = s[prs, "Std. Error"],
        P = s[prs, "Pr(>|t|)"], N = nrow(df),
        stringsAsFactors = FALSE
      )
    }
  }

  # Age-stratified cognition
  strat_results <- list()
  for (age_grp in c("Younger", "Older")) {
    sub <- data[data$Age_Group == age_grp, ]
    for (prs in prs_vars) {
      for (outcome in outcomes) {
        df <- sub[complete.cases(sub[, c(prs, outcome, covars)]), ]
        if (nrow(df) < 20) next

        fml <- as.formula(paste(outcome, "~", prs, "+",
                                paste(covars, collapse = " + ")))
        fit <- lm(fml, data = df)
        s <- summary(fit)$coefficients
        if (!prs %in% rownames(s)) next

        strat_results[[length(strat_results) + 1]] <- data.frame(
          Age_Group = age_grp, PRS = prs, Outcome = outcome,
          Beta = s[prs, "Estimate"], SE = s[prs, "Std. Error"],
          P = s[prs, "Pr(>|t|)"], N = nrow(df),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  full_df <- do.call(rbind, results)
  strat_df <- do.call(rbind, strat_results)
  if (!is.null(full_df) && nrow(full_df) > 0) full_df$P_FDR <- p.adjust(full_df$P, method = "BH")
  if (!is.null(strat_df) && nrow(strat_df) > 0) strat_df$P_FDR <- p.adjust(strat_df$P, method = "BH")

  list(full = full_df, stratified = strat_df)
}


# ==============================================================================
# 10. PATHWAY ENRICHMENT COMPARISON (MICROGLIA vs ABETA CLEARANCE)
# ==============================================================================
#
# Compares standardized PRS effect sizes across pathways to quantify
# relative enrichment of microglia/oligodendrocyte vs Abeta clearance
# pathways in EOAD genetic architecture.
#

run_pathway_enrichment_comparison <- function(data) {

  # Pathway PRS variables grouped by biological function
  microglia_oligo_prs <- c("PRS_EOAD_Microglia", "PRS_EOAD_Oligo", "PRS_EOAD_Myelin")
  abeta_prs <- c("PRS_EOAD_Abeta", "PRS_EOAD_APP")
  all_prs <- c(microglia_oligo_prs, abeta_prs, "PRS_EOAD_Global",
               "PRS_EOAD_Lipid", "PRS_EOAD_TCell")
  all_prs <- all_prs[all_prs %in% colnames(data)]

  outcomes <- c("WMH_Log", "Hippo_Norm", "ABETA", "TAU", "PTAU")
  if ("STREM2" %in% colnames(data)) outcomes <- c(outcomes, "STREM2")
  outcomes <- outcomes[outcomes %in% colnames(data)]
  covars <- c("AGE", "PTGENDER", "APOE4", "PTEDUCAT")

  # --- Standardized effect sizes for each PRS on each outcome ---
  effect_results <- list()
  for (prs in all_prs) {
    for (outcome in outcomes) {
      df <- data[complete.cases(data[, c(prs, outcome, covars)]), ]
      if (nrow(df) < 30) next

      # Standardize PRS and outcome for comparable effect sizes
      df$PRS_z <- scale(df[[prs]])
      df$Outcome_z <- scale(df[[outcome]])

      fml <- as.formula(paste("Outcome_z ~ PRS_z +",
                              paste(covars, collapse = " + ")))
      fit <- tryCatch(lm(fml, data = df), error = function(e) NULL)
      if (is.null(fit)) next
      s <- summary(fit)$coefficients
      if (!"PRS_z" %in% rownames(s)) next

      pathway_group <- ifelse(prs %in% microglia_oligo_prs, "Microglia_Oligo",
                       ifelse(prs %in% abeta_prs, "Abeta_Clearance", "Other"))

      effect_results[[length(effect_results) + 1]] <- data.frame(
        PRS = prs, Pathway_Group = pathway_group, Outcome = outcome,
        Std_Beta = s["PRS_z", "Estimate"],
        SE = s["PRS_z", "Std. Error"],
        P = s["PRS_z", "Pr(>|t|)"],
        N = nrow(df),
        stringsAsFactors = FALSE
      )
    }
  }

  effect_df <- do.call(rbind, effect_results)
  if (nrow(effect_df) == 0) return(list(effects = data.frame(), enrichment = data.frame()))
  effect_df$P_FDR <- p.adjust(effect_df$P, method = "BH")

  # --- Enrichment ratio: mean |Std_Beta| for microglia/oligo vs Abeta ---
  enrichment <- effect_df %>%
    filter(Pathway_Group %in% c("Microglia_Oligo", "Abeta_Clearance")) %>%
    group_by(Pathway_Group, Outcome) %>%
    summarise(
      Mean_Abs_Beta = mean(abs(Std_Beta), na.rm = TRUE),
      N_PRS = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = Pathway_Group,
                values_from = c(Mean_Abs_Beta, N_PRS)) %>%
    mutate(
      Enrichment_Ratio = Mean_Abs_Beta_Microglia_Oligo / Mean_Abs_Beta_Abeta_Clearance
    )

  # --- Age-stratified enrichment ---
  age_enrichment <- list()
  for (age_grp in c("Younger", "Older")) {
    sub <- data[data$Age_Group == age_grp, ]
    for (prs in all_prs) {
      for (outcome in outcomes) {
        df <- sub[complete.cases(sub[, c(prs, outcome, covars)]), ]
        if (nrow(df) < 20) next
        df$PRS_z <- scale(df[[prs]])
        df$Outcome_z <- scale(df[[outcome]])
        fml <- as.formula(paste("Outcome_z ~ PRS_z +",
                                paste(covars, collapse = " + ")))
        fit <- tryCatch(lm(fml, data = df), error = function(e) NULL)
        if (is.null(fit)) next
        s <- summary(fit)$coefficients
        if (!"PRS_z" %in% rownames(s)) next

        pathway_group <- ifelse(prs %in% microglia_oligo_prs, "Microglia_Oligo",
                         ifelse(prs %in% abeta_prs, "Abeta_Clearance", "Other"))

        age_enrichment[[length(age_enrichment) + 1]] <- data.frame(
          Age_Group = age_grp, PRS = prs, Pathway_Group = pathway_group,
          Outcome = outcome,
          Std_Beta = s["PRS_z", "Estimate"],
          P = s["PRS_z", "Pr(>|t|)"],
          N = nrow(df),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  age_df <- if (length(age_enrichment) > 0) do.call(rbind, age_enrichment) else data.frame()

  list(effects = effect_df, enrichment = enrichment, age_stratified = age_df)
}


# ==============================================================================
# 11. EXPLORATORY BOUNDARY CONDITION ANALYSES
# ==============================================================================

run_boundary_analyses <- function(data) {

  covars <- c("AGE", "PTGENDER", "APOE4", "PTEDUCAT")

  # --- Unstratified conversion (pooled logistic) ---
  dx_bl_col <- detect_dx_bl_col(data)
  mci <- data[data[[dx_bl_col]] %in% c("LMCI", "EMCI", "MCI"), ]
  unstrat_results <- list()
  for (prs in c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin")) {
    if (!prs %in% colnames(mci)) next
    df <- mci[complete.cases(mci[, c(prs, "Converted", covars)]), ]
    if (nrow(df) < 20 || sum(df$Converted) < 5) next
    fml <- as.formula(paste("Converted ~", prs, "+",
                            paste(covars, collapse = " + ")))
    fit <- glm(fml, data = df, family = binomial)
    ci <- confint.default(fit)
    s <- summary(fit)$coefficients[prs, ]
    unstrat_results[[prs]] <- data.frame(
      PRS = prs,
      OR = exp(s["Estimate"]),
      OR_Lower = exp(ci[prs, 1]), OR_Upper = exp(ci[prs, 2]),
      P = s["Pr(>|z|)"], N = nrow(df),
      stringsAsFactors = FALSE
    )
  }
  unstrat_df <- do.call(rbind, unstrat_results)

  # --- Pathway specificity (APP, T-cell as null results) ---
  specificity_results <- list()
  alt_prs <- c("PRS_EOAD_APP", "PRS_EOAD_TCell")
  alt_prs <- alt_prs[alt_prs %in% colnames(data)]
  spec_outcomes <- c("ABETA", "TAU", "PTAU", "Hippo_Norm", "WMH_Log")
  for (prs in alt_prs) {
    for (outcome in spec_outcomes) {
      if (!outcome %in% colnames(data)) next
      df <- data[complete.cases(data[, c(prs, outcome, covars)]), ]
      if (nrow(df) < 30) next
      fml <- as.formula(paste(outcome, "~", prs, "+",
                              paste(covars, collapse = " + ")))
      fit <- lm(fml, data = df)
      s <- summary(fit)$coefficients[prs, ]
      specificity_results[[length(specificity_results) + 1]] <- data.frame(
        PRS = prs, Outcome = outcome,
        Beta = s["Estimate"], SE = s["Std. Error"],
        P = s["Pr(>|t|)"], N = nrow(df),
        stringsAsFactors = FALSE
      )
    }
  }
  specificity_df <- do.call(rbind, specificity_results)

  # --- Continuous age x Oligo-PRS interaction on WMH ---
  wmh_interaction <- NULL
  df_wmh <- data[complete.cases(data[, c("PRS_EOAD_Oligo", "WMH_Log",
                                          "AGE", "PTGENDER", "APOE4", "ST10CV")]), ]
  if (nrow(df_wmh) > 50) {
    fit <- lm(WMH_Log ~ PRS_EOAD_Oligo * AGE + PTGENDER + APOE4 + ST10CV,
              data = df_wmh)
    int_s <- summary(fit)$coefficients["PRS_EOAD_Oligo:AGE", ]
    wmh_interaction <- data.frame(
      Beta_Interaction = int_s["Estimate"],
      SE = int_s["Std. Error"],
      P = int_s["Pr(>|t|)"],
      N = nrow(df_wmh)
    )
  }

  list(
    unstratified_conversion = unstrat_df,
    pathway_specificity = specificity_df,
    wmh_age_interaction = wmh_interaction
  )
}


# ==============================================================================
# 12. DEMOGRAPHICS TABLE
# ==============================================================================

generate_demographics <- function(data) {

  dx_bl_col <- detect_dx_bl_col(data)

  demo <- data %>%
    mutate(DX_Group = case_when(
      .data[[dx_bl_col]] %in% c("CN", "NL") ~ "CN",
      .data[[dx_bl_col]] %in% c("LMCI", "EMCI", "MCI") ~ "MCI",
      .data[[dx_bl_col]] %in% c("Dementia", "AD") ~ "Dementia",
      TRUE ~ "Other"
    )) %>%
    group_by(DX_Group) %>%
    summarise(
      N = n(),
      Age_Mean = mean(AGE, na.rm = TRUE),
      Age_SD = sd(AGE, na.rm = TRUE),
      Female_N = sum(PTGENDER == "Female", na.rm = TRUE),
      Female_Pct = mean(PTGENDER == "Female", na.rm = TRUE) * 100,
      Education_Mean = mean(PTEDUCAT, na.rm = TRUE),
      Education_SD = sd(PTEDUCAT, na.rm = TRUE),
      APOE4_0 = sum(APOE4 == 0, na.rm = TRUE),
      APOE4_1 = sum(APOE4 == 1, na.rm = TRUE),
      APOE4_2 = sum(APOE4 == 2, na.rm = TRUE),
      MMSE_Mean = mean(MMSE, na.rm = TRUE),
      MMSE_SD = sd(MMSE, na.rm = TRUE),
      .groups = "drop"
    )

  # MCI conversion summary
  mci_data <- data[data[[dx_bl_col]] %in% c("LMCI", "EMCI", "MCI"), ]
  conversion_summary <- data.frame(
    N_MCI = nrow(mci_data),
    N_Converted = sum(mci_data$Converted, na.rm = TRUE),
    Conversion_Rate = mean(mci_data$Converted, na.rm = TRUE) * 100,
    Median_Followup_Years = unique(data$Median_Followup_Years)[1]
  )

  list(demographics = demo, conversion_summary = conversion_summary)
}


# ==============================================================================
# 13. VISUALIZATION FUNCTIONS
# ==============================================================================

plot_sliding_window <- function(sw_results, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  wmh <- sw_results$wmh
  if (is.null(wmh) || nrow(wmh) == 0) return(NULL)

  wmh$Sig <- ifelse(wmh$P_FDR < 0.05, "FDR < 0.05",
                    ifelse(wmh$P < 0.05, "Nominal", "NS"))

  p_wmh <- ggplot(wmh, aes(x = Centre_Age, y = Beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = Beta - 1.96 * SE, ymax = Beta + 1.96 * SE),
                alpha = 0.2, fill = "#2166AC") +
    geom_line(colour = "#2166AC", linewidth = 1) +
    geom_point(aes(shape = Sig, fill = Sig), size = 3) +
    scale_shape_manual(values = c("FDR < 0.05" = 21, "Nominal" = 22, "NS" = 23)) +
    scale_fill_manual(values = c("FDR < 0.05" = "#2166AC",
                                  "Nominal" = "#92C5DE", "NS" = "white")) +
    labs(x = "Window centre age (years)", y = "Beta (Oligo-PRS on log-WMH)",
         title = "Sliding-window: WMH burden") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank())

  conv <- sw_results$conversion
  p_conv <- NULL
  if (!is.null(conv) && nrow(conv) > 0) {
    conv$Sig <- ifelse(conv$P_FDR < 0.05, "FDR < 0.05",
                       ifelse(conv$P < 0.05, "Nominal", "NS"))

    p_conv <- ggplot(conv, aes(x = Centre_Age, y = OR)) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
      geom_ribbon(aes(ymin = OR_Lower, ymax = OR_Upper),
                  alpha = 0.2, fill = "#B2182B") +
      geom_line(colour = "#B2182B", linewidth = 1) +
      geom_point(aes(shape = Sig, fill = Sig), size = 3) +
      scale_shape_manual(values = c("FDR < 0.05" = 21, "Nominal" = 22, "NS" = 23)) +
      scale_fill_manual(values = c("FDR < 0.05" = "#B2182B",
                                    "Nominal" = "#F4A582", "NS" = "white")) +
      labs(x = "Window centre age (years)", y = "Odds ratio (MCI-to-AD)",
           title = "Sliding-window: MCI conversion risk") +
      theme_classic(base_size = 12) +
      theme(legend.position = "bottom", legend.title = element_blank())
  }

  if (!is.null(p_conv)) {
    combined <- plot_grid(p_wmh, p_conv, ncol = 2, labels = c("A", "B"))
    ggsave(file.path(output_dir, "sliding_window_combined.pdf"),
           combined, width = 12, height = 5)
  } else {
    ggsave(file.path(output_dir, "sliding_window_wmh.pdf"),
           p_wmh, width = 6, height = 5)
  }

  list(wmh = p_wmh, conversion = p_conv)
}

plot_forest <- function(conversion_results, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  strat <- conversion_results$stratified
  if (is.null(strat) || nrow(strat) == 0) return(NULL)

  strat$Label <- paste0(strat$PRS, " (", strat$Stratum, ")")
  strat <- strat[order(strat$PRS, strat$Stratum), ]

  p <- ggplot(strat, aes(x = OR, y = reorder(Label, OR))) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = OR_Lower, xmax = OR_Upper), height = 0.2) +
    geom_point(aes(colour = Stratum), size = 3) +
    scale_colour_manual(values = c("Younger" = "#D73027", "Older" = "#4575B4")) +
    labs(x = "Odds ratio (95% CI)", y = NULL,
         title = "Age-stratified MCI-to-AD conversion") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, "forest_plot_conversion.pdf"), p, width = 8, height = 6)
  return(p)
}

plot_cluster_heatmap <- function(subtype_results, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  centres <- subtype_results$centres
  label_map <- subtype_results$label_map
  rownames(centres) <- label_map[as.character(centres$Cluster)]
  mat <- as.matrix(centres[, grep("^PRS_", colnames(centres))])
  colnames(mat) <- gsub("PRS_EOAD_", "", colnames(mat))

  pdf(file.path(output_dir, "cluster_heatmap.pdf"), width = 8, height = 4)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           display_numbers = TRUE, number_format = "%.2f",
           main = "Genetic subtype PRS profiles")
  dev.off()
}

plot_dose_response <- function(dose_results, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  means_df <- dose_results$quantile_means
  if (is.null(means_df) || nrow(means_df) == 0) return(NULL)

  # Focus on sTREM2 if available, otherwise first outcome
  if ("STREM2" %in% means_df$Outcome) {
    plot_data <- means_df[means_df$Outcome == "STREM2", ]
    y_label <- "CSF sTREM2 (pg/mL)"
  } else {
    first_outcome <- unique(means_df$Outcome)[1]
    plot_data <- means_df[means_df$Outcome == first_outcome, ]
    y_label <- first_outcome
  }

  if (nrow(plot_data) == 0) return(NULL)

  p <- ggplot(plot_data, aes(x = factor(PRS_Quantile), y = Mean,
                              fill = APOE_Group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = Mean - SD / sqrt(N), ymax = Mean + SD / sqrt(N)),
                  position = position_dodge(width = 0.8), width = 0.2) +
    facet_wrap(~ PRS, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "PRS quartile", y = y_label,
         title = "Dose-response: PRS quartile vs biomarker level") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, "dose_response_quartiles.pdf"),
         p, width = 12, height = 8)
  return(p)
}

plot_brain_surface <- function(region_results, output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (!requireNamespace("ggseg", quietly = TRUE)) {
    message("ggseg package required for brain surface maps. Install with:")
    message("  install.packages('ggseg')")
    # Fallback: save a heatmap-style table instead
    if (!is.null(region_results) && nrow(region_results) > 0) {
      fwrite(region_results, file.path(output_dir, "brain_region_results_table.csv"))
    }
    return(NULL)
  }

  library(ggseg)

  for (sub in unique(region_results$Subtype)) {
    sub_data <- region_results[region_results$Subtype == sub, ]

    # Build ggseg-compatible label: "lh_bankssts" or "rh_bankssts"
    sub_data$label <- paste0(sub_data$Hemisphere, "_", sub_data$Region)

    # Merge with dk atlas data
    atlas_data <- as_tibble(dk) %>%
      select(label, hemi, region) %>%
      distinct()

    plot_df <- left_join(atlas_data, sub_data, by = "label")

    p <- ggplot(plot_df) +
      geom_brain(atlas = dk, aes(fill = t_stat),
                 position = position_brain(hemi ~ side)) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "t-statistic",
                           na.value = "grey90") +
      labs(title = paste("Regional cortical thickness:", sub, "vs Background_Risk")) +
      theme_void(base_size = 12)

    ggsave(file.path(output_dir, paste0("brain_surface_", sub, ".pdf")),
           p, width = 10, height = 6)
  }
}


# ==============================================================================
# 14. MAIN EXECUTION
# ==============================================================================

run_adni_validation <- function(clinical_file, prs_file, mri_file, csf_file,
                                 strem2_file = NULL,
                                 output_dir = "results/adni_validation") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  cat("=== ADNI Independent Cohort Validation ===\n\n")

  # --- Step 1: Load and prepare data ---
  cat("[1/12] Loading and integrating data...\n")
  data <- load_adni_data(clinical_file, prs_file, mri_file, csf_file, strem2_file)
  data <- create_derived_variables(data)

  dx_bl_col <- detect_dx_bl_col(data)
  n_mci <- sum(data[[dx_bl_col]] %in% c("LMCI", "EMCI", "MCI"))
  cat(sprintf("  Total N = %d\n", nrow(data)))
  cat(sprintf("  MCI N = %d\n", n_mci))
  cat(sprintf("  Longitudinal conversions = %d\n", sum(data$Converted, na.rm = TRUE)))
  cat(sprintf("  Median follow-up = %.1f years\n", unique(data$Median_Followup_Years)[1]))
  if ("STREM2" %in% colnames(data)) {
    cat(sprintf("  sTREM2 available: N = %d\n", sum(!is.na(data$STREM2))))
  }

  # --- Step 2: APOE-stratified PRS-biomarker ---
  cat("[2/12] APOE-stratified PRS-biomarker analysis...\n")
  apoe_results <- run_apoe_stratified_analysis(data)
  fwrite(apoe_results, file.path(output_dir, "apoe_stratified_results.csv"))

  # --- Step 2b: Dose-response analysis ---
  cat("  Running dose-response (PRS quartile) analysis...\n")
  dose_results <- run_dose_response_analysis(data)
  if (nrow(dose_results$trend) > 0) {
    fwrite(dose_results$trend, file.path(output_dir, "dose_response_trend.csv"))
  }
  if (nrow(dose_results$quantile_means) > 0) {
    fwrite(dose_results$quantile_means, file.path(output_dir, "dose_response_quantile_means.csv"))
  }
  plot_dose_response(dose_results, fig_dir)

  # --- Step 3: Longitudinal MCI conversion ---
  cat("[3/12] Longitudinal MCI-to-AD conversion analysis...\n")
  conv_results <- run_conversion_analysis(data)
  fwrite(conv_results$full, file.path(output_dir, "conversion_full.csv"))
  fwrite(conv_results$stratified, file.path(output_dir, "conversion_stratified.csv"))
  fwrite(conv_results$interaction, file.path(output_dir, "conversion_interaction.csv"))
  if (!is.null(conv_results$sensitivity_under65)) {
    fwrite(conv_results$sensitivity_under65,
           file.path(output_dir, "conversion_sensitivity_under65.csv"))
  }

  # --- Step 4: PRS-neuroimaging with age interactions ---
  cat("[4/12] PRS-neuroimaging associations...\n")
  neuro_results <- run_neuroimaging_analysis(data)
  fwrite(neuro_results$stratified, file.path(output_dir, "neuroimaging_stratified.csv"))
  fwrite(neuro_results$interaction, file.path(output_dir, "neuroimaging_interaction.csv"))

  # --- Step 5: Sliding-window analysis ---
  cat("[5/12] Sliding-window age-dependency mapping...\n")
  sw_results <- run_sliding_window(data)
  if (nrow(sw_results$wmh) > 0) {
    fwrite(sw_results$wmh, file.path(output_dir, "sliding_window_wmh.csv"))
    cat(sprintf("  WMH peak centre age: %.0f years\n", sw_results$wmh$Peak_Centre_Age[1]))
    cat(sprintf("  WMH dissipation: %.1f%%\n", sw_results$wmh$Dissipation_Pct[1]))
  }
  if (nrow(sw_results$conversion) > 0) {
    fwrite(sw_results$conversion, file.path(output_dir, "sliding_window_conversion.csv"))
    cat(sprintf("  Conversion peak centre age: %.0f years\n",
                sw_results$conversion$Peak_Centre_Age[1]))
    cat(sprintf("  Conversion dissipation: %.1f%%\n",
                sw_results$conversion$Dissipation_Pct[1]))
  }
  plot_sliding_window(sw_results, fig_dir)

  # --- Step 6: Genetic subtyping ---
  cat("[6/12] Unsupervised genetic subtyping...\n")
  subtype_results <- run_genetic_subtyping(data)
  data <- subtype_results$data
  cat(sprintf("  Concordance (k-means vs hierarchical): %.1f%%\n",
              subtype_results$concordance * 100))
  cat(sprintf("  Silhouette: %.3f\n", subtype_results$silhouette))
  if (!is.na(subtype_results$kappa)) {
    cat(sprintf("  Cohen's kappa: %.3f (P = %s)\n",
                subtype_results$kappa,
                format.pval(subtype_results$kappa_p, digits = 3)))
  }
  fwrite(subtype_results$pathway_burden,
         file.path(output_dir, "pathway_burden.csv"))
  fwrite(subtype_results$subtype_anova,
         file.path(output_dir, "subtype_anova.csv"))
  fwrite(subtype_results$subtype_summary,
         file.path(output_dir, "subtype_summary.csv"))
  plot_cluster_heatmap(subtype_results, fig_dir)
  plot_forest(conv_results, fig_dir)

  # --- Step 7: Regional cortical thickness ---
  cat("[7/12] Regional cortical thickness analysis (68 DK regions)...\n")
  region_results <- run_brain_region_analysis(data)
  if (!is.null(region_results) && nrow(region_results) > 0) {
    fwrite(region_results, file.path(output_dir, "cortical_thickness_subtype.csv"))
    n_regions <- length(unique(paste(region_results$Region, region_results$Hemisphere)))
    cat(sprintf("  Regions analysed: %d\n", n_regions))
    plot_brain_surface(region_results, fig_dir)
  }

  region_interaction <- run_regional_age_prs_interaction(data)
  if (!is.null(region_interaction) && nrow(region_interaction) > 0) {
    fwrite(region_interaction,
           file.path(output_dir, "cortical_thickness_age_prs_interaction.csv"))
  }

  # --- Step 8: Subcortical structures ---
  cat("[8/12] Subcortical structure analysis...\n")
  subcort_results <- run_subcortical_analysis(data)
  if (nrow(subcort_results$subtype) > 0) {
    fwrite(subcort_results$subtype,
           file.path(output_dir, "subcortical_volume_subtype.csv"))
  }
  if (nrow(subcort_results$interaction) > 0) {
    fwrite(subcort_results$interaction,
           file.path(output_dir, "subcortical_volume_age_prs_interaction.csv"))
  }

  # --- Step 9: Cross-sectional PRS-cognition ---
  cat("[9/12] Cross-sectional PRS-cognition analysis...\n")
  cognition_results <- run_cross_sectional_cognition(data)
  if (!is.null(cognition_results$full) && nrow(cognition_results$full) > 0) {
    fwrite(cognition_results$full,
           file.path(output_dir, "cognition_prs_full.csv"))
  }
  if (!is.null(cognition_results$stratified) && nrow(cognition_results$stratified) > 0) {
    fwrite(cognition_results$stratified,
           file.path(output_dir, "cognition_prs_stratified.csv"))
  }

  # --- Step 10: Pathway enrichment comparison ---
  cat("[10/12] Pathway enrichment comparison (microglia vs Abeta)...\n")
  enrichment_results <- run_pathway_enrichment_comparison(data)
  if (nrow(enrichment_results$effects) > 0) {
    fwrite(enrichment_results$effects,
           file.path(output_dir, "pathway_enrichment_effects.csv"))
  }
  if (nrow(enrichment_results$enrichment) > 0) {
    fwrite(enrichment_results$enrichment,
           file.path(output_dir, "pathway_enrichment_ratio.csv"))
  }
  if (!is.null(enrichment_results$age_stratified) && nrow(enrichment_results$age_stratified) > 0) {
    fwrite(enrichment_results$age_stratified,
           file.path(output_dir, "pathway_enrichment_age_stratified.csv"))
  }

  # --- Step 11: Boundary condition analyses ---
  cat("[11/12] Exploratory boundary condition analyses...\n")
  boundary_results <- run_boundary_analyses(data)
  if (!is.null(boundary_results$unstratified_conversion)) {
    fwrite(boundary_results$unstratified_conversion,
           file.path(output_dir, "boundary_unstratified_conversion.csv"))
  }
  if (!is.null(boundary_results$pathway_specificity)) {
    fwrite(boundary_results$pathway_specificity,
           file.path(output_dir, "boundary_pathway_specificity.csv"))
  }

  # --- Step 12: Demographics ---
  cat("[12/12] Generating demographics table...\n")
  demo_results <- generate_demographics(data)
  fwrite(demo_results$demographics, file.path(output_dir, "demographics.csv"))
  fwrite(demo_results$conversion_summary,
         file.path(output_dir, "conversion_summary.csv"))

  # --- Summary ---
  cat("\n=== Analysis complete ===\n")
  cat(sprintf("Results saved to: %s\n", output_dir))
  cat(sprintf("Figures saved to: %s\n", fig_dir))

  invisible(list(
    data = data,
    apoe = apoe_results,
    dose_response = dose_results,
    conversion = conv_results,
    neuroimaging = neuro_results,
    sliding_window = sw_results,
    subtypes = subtype_results,
    cortical_thickness = region_results,
    cortical_interaction = region_interaction,
    subcortical = subcort_results,
    cognition = cognition_results,
    pathway_enrichment = enrichment_results,
    boundary = boundary_results,
    demographics = demo_results
  ))
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================
#
# results <- run_adni_validation(
#   clinical_file = "data/ADNIMERGE.csv",
#   prs_file      = "data/ADNI_pathway_PRS.csv",
#   mri_file      = "data/UCSFFSX51_ADNI1GO2.csv",
#   csf_file      = "data/UPENNBIOMK_MASTER.csv",
#   strem2_file   = "data/ADNI_CSF_sTREM2.csv",
#   output_dir    = "results/adni_validation"
# )
