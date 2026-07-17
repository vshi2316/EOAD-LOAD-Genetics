# ==============================================================================
# ADNI Scanner-Aware MRI / WMH Sensitivity Analyses
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(broom)
})


# ==============================================================================
# Utility functions
# ==============================================================================

sex_numeric <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x))
  out[grepl("^m", x)] <- 1
  out[grepl("^f", x)] <- 0
  out
}

normalize_manufacturer <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep("Other", length(x))
  out[grepl("siemens", x)] <- "Siemens"
  out[grepl("^ge$|general electric|\\bge\\b", x)] <- "GE"
  out[grepl("philips", x)] <- "Philips"
  out
}

collapse_scanner_model <- function(x, min_n = 10) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- "Unknown"
  counts <- table(x, useNA = "no")
  keep <- names(counts[counts >= min_n])
  ifelse(x %in% keep, x, "Other")
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, format.pval(p, digits = 3, eps = 1e-300))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


# ==============================================================================
# OLS fitting helper
# ==============================================================================

fit_lm_term <- function(data, outcome, term, covariates,
                        categorical_covariates = NULL, min_n = 30) {
  categorical_covariates <- categorical_covariates %||% character(0)
  needed <- unique(c(outcome, term, covariates, categorical_covariates))
  missing <- setdiff(needed, colnames(data))

  if (length(missing) > 0) {
    return(NULL)
  }

  d <- data[, needed, drop = FALSE]
  d <- d[complete.cases(d), , drop = FALSE]

  if (nrow(d) < min_n) {
    return(NULL)
  }

  for (col in setdiff(c(outcome, term, covariates), categorical_covariates)) {
    d[[col]] <- suppressWarnings(as.numeric(d[[col]]))
  }

  for (col in categorical_covariates) {
    d[[col]] <- factor(d[[col]])
  }

  if (sum(is.finite(d[[outcome]])) < min_n ||
      sum(is.finite(d[[term]])) < min_n) {
    return(NULL)
  }

  rhs <- c(term, covariates, categorical_covariates)
  rhs <- rhs[rhs != outcome]
  formula <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))

  fit <- tryCatch(lm(formula, data = d), error = function(e) NULL)
  if (is.null(fit)) {
    return(NULL)
  }

  tidy_fit <- broom::tidy(fit, conf.int = TRUE)
  row <- tidy_fit[tidy_fit$term == term, , drop = FALSE]

  if (nrow(row) == 0) {
    return(NULL)
  }

  data.frame(
    N = nobs(fit),
    DF_residual = df.residual(fit),
    Beta = row$estimate,
    SE = row$std.error,
    CI95_low = row$conf.low,
    CI95_high = row$conf.high,
    t = row$statistic,
    P = row$p.value,
    R2 = summary(fit)$r.squared,
    Covariates = paste(c(covariates, categorical_covariates), collapse = "; "),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Exact matching to imaging metadata
# ==============================================================================

exact_match_freesurfer <- function(final_data, freesurfer_metadata) {
  f <- final_data %>%
    filter(!is.na(ST128SV)) %>%
    select(RID, ST10CV, ST128SV)

  required <- c("RID", "ST10CV", "ST128SV", "PHASE", "FIELD_STRENGTH",
                "VISCODE", "VISCODE2", "STATUS")
  missing <- setdiff(required, colnames(freesurfer_metadata))
  if (length(missing) > 0) {
    stop("FreeSurfer metadata is missing columns: ", paste(missing, collapse = ", "))
  }

  matched <- f %>%
    left_join(
      freesurfer_metadata %>% select(all_of(required)),
      by = "RID",
      suffix = c("_final", "_raw")
    ) %>%
    filter(
      round(as.numeric(ST10CV_final), 3) == round(as.numeric(ST10CV_raw), 3),
      round(as.numeric(ST128SV_final), 3) == round(as.numeric(ST128SV_raw), 3)
    ) %>%
    distinct(RID, .keep_all = TRUE) %>%
    transmute(
      RID,
      FS_PHASE = PHASE,
      FS_FIELD_STRENGTH = FIELD_STRENGTH,
      FS_VISCODE = VISCODE,
      FS_VISCODE2 = VISCODE2,
      FS_STATUS = STATUS
    )

  matched
}

exact_match_ucd_wmh <- function(stratified_data, ucd_wmh_metadata) {
  f <- stratified_data %>%
    filter(!is.na(WMH_Volume)) %>%
    select(RID, WMH_Volume)

  required <- c("RID", "TOTAL_WMH", "PHASE", "VISCODE", "VISCODE2",
                "MANUFACTURER", "MANUFACTURERSMODELNAME",
                "MAGNETICFIELDSTRENGTH")
  missing <- setdiff(required, colnames(ucd_wmh_metadata))
  if (length(missing) > 0) {
    stop("UCD WMH metadata is missing columns: ", paste(missing, collapse = ", "))
  }

  matched <- f %>%
    left_join(ucd_wmh_metadata %>% select(all_of(required)), by = "RID") %>%
    filter(round(as.numeric(WMH_Volume), 4) ==
             round(as.numeric(TOTAL_WMH), 4)) %>%
    distinct(RID, .keep_all = TRUE) %>%
    mutate(
      UCD_MANUFACTURER = normalize_manufacturer(MANUFACTURER),
      UCD_SCANNER_MODEL = collapse_scanner_model(MANUFACTURERSMODELNAME,
                                                 min_n = 10)
    ) %>%
    transmute(
      RID,
      UCD_PHASE = PHASE,
      UCD_VISCODE = VISCODE,
      UCD_VISCODE2 = VISCODE2,
      UCD_MANUFACTURER,
      UCD_SCANNER_MODEL,
      UCD_FIELD_STRENGTH = MAGNETICFIELDSTRENGTH
    )

  matched
}

summarize_batch_counts <- function(data, dataset_label, batch_variables) {
  rows <- list()

  for (var in batch_variables) {
    if (!var %in% colnames(data)) next
    tab <- as.data.frame(table(data[[var]], useNA = "ifany"),
                         stringsAsFactors = FALSE)
    colnames(tab) <- c("Level", "N")
    tab$Dataset <- dataset_label
    tab$Batch_variable <- var
    rows[[length(rows) + 1]] <- tab %>%
      select(Dataset, Batch_variable, Level, N)
  }

  bind_rows(rows)
}


# ==============================================================================
# Sensitivity model runner
# ==============================================================================

append_model_result <- function(results, dataset_label, outcome, prs_term,
                                model_label, data, covariates,
                                categorical_covariates = NULL, notes = "") {
  fit <- fit_lm_term(
    data = data,
    outcome = outcome,
    term = prs_term,
    covariates = covariates,
    categorical_covariates = categorical_covariates
  )

  if (is.null(fit)) {
    return(results)
  }

  row <- data.frame(
    Dataset = dataset_label,
    Outcome = outcome,
    PRS_term = prs_term,
    Model = model_label,
    Notes = notes,
    stringsAsFactors = FALSE
  ) %>%
    bind_cols(fit)

  bind_rows(results, row)
}

residualize_batch_preserving_term <- function(data, outcome, term, covariates,
                                              batch, new_col) {
  needed <- unique(c(outcome, term, covariates, batch))
  missing <- setdiff(needed, colnames(data))
  data[[new_col]] <- NA_real_

  if (length(missing) > 0) {
    return(data)
  }

  d <- data[, needed, drop = FALSE]
  d <- d[complete.cases(d), , drop = FALSE]

  if (nrow(d) < 30 || length(unique(d[[batch]])) < 2) {
    return(data)
  }

  for (col in c(outcome, term, covariates)) {
    d[[col]] <- suppressWarnings(as.numeric(d[[col]]))
  }
  d[[batch]] <- factor(d[[batch]])

  formula <- as.formula(
    paste(outcome, "~", paste(c(term, covariates, batch), collapse = " + "))
  )
  fit <- lm(formula, data = d)

  mm <- model.matrix(fit)
  batch_cols <- grep(paste0("^", batch), colnames(mm), value = TRUE)

  if (length(batch_cols) == 0) {
    return(data)
  }

  batch_effect <- as.numeric(mm[, batch_cols, drop = FALSE] %*%
                               coef(fit)[batch_cols])
  adjusted <- d[[outcome]] - batch_effect
  data[rownames(d), new_col] <- adjusted
  data
}


# ==============================================================================
# Main execution wrapper
# ==============================================================================

run_adni_scanner_aware_sensitivity <- function(
    final_analytic_file,
    stratified_analytic_file,
    freesurfer_metadata_file,
    ucd_wmh_metadata_file,
    output_dir = "results/adni_scanner_sensitivity",
    predictors = c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin"),
    covariates = c("AGE", "Sex_Numeric", "APOE4", "PTEDUCAT")) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  cat("=== ADNI Scanner-Aware Sensitivity Analyses ===\n\n")
  cat(sprintf("Output directory: %s\n", output_dir))

  final <- fread(final_analytic_file) %>% as.data.frame()
  stratified <- fread(stratified_analytic_file) %>% as.data.frame()
  freesurfer <- fread(freesurfer_metadata_file) %>% as.data.frame()
  ucd_wmh <- fread(ucd_wmh_metadata_file) %>% as.data.frame()

  final$Sex_Numeric <- sex_numeric(final$PTGENDER)
  stratified$Sex_Numeric <- sex_numeric(stratified$PTGENDER)

  final$Age_Group_70_clean <- ifelse(final$AGE < 70,
                                     "Younger (<70)", "Older (>=70)")
  stratified$Age_Group_70_clean <- ifelse(stratified$AGE < 70,
                                          "Younger (<70)", "Older (>=70)")

  cat("[1/5] Exact-matching FreeSurfer MRI metadata...\n")
  fs_batch <- exact_match_freesurfer(final, freesurfer)
  final <- final %>% left_join(fs_batch, by = "RID")

  cat("[2/5] Exact-matching UCD FLAIR WMH metadata...\n")
  ucd_batch <- exact_match_ucd_wmh(stratified, ucd_wmh)
  stratified <- stratified %>% left_join(ucd_batch, by = "RID")
  stratified$UCD_Log_WMH <- log1p(suppressWarnings(as.numeric(stratified$WMH_Volume)))

  fwrite(final, file.path(output_dir, "ADNI_FreeSurfer_batch_enriched.csv"))
  fwrite(stratified, file.path(output_dir, "ADNI_UCD_WMH_batch_enriched.csv"))

  fs_nonmissing <- final %>% filter(!is.na(ST128SV))
  ucd_nonmissing <- stratified %>% filter(!is.na(WMH_Volume))

  batch_counts <- bind_rows(
    summarize_batch_counts(
      fs_nonmissing,
      "FreeSurfer_ST128SV",
      c("FS_PHASE", "FS_FIELD_STRENGTH", "FS_VISCODE2", "FS_STATUS")
    ),
    summarize_batch_counts(
      ucd_nonmissing,
      "UCD_WMH",
      c("UCD_PHASE", "UCD_FIELD_STRENGTH", "UCD_MANUFACTURER",
        "UCD_SCANNER_MODEL", "UCD_VISCODE2")
    )
  )
  fwrite(batch_counts, file.path(output_dir, "ADNI_scanner_batch_counts.csv"))

  cat("[3/5] Running FreeSurfer MRI field-strength sensitivity models...\n")
  fs_outcomes <- c("WMH_Log", "WMH_ICV_Ratio", "Ventricle_ICV_Ratio",
                   "Hippo_ICV_Ratio", "Entorhinal_ICV_Ratio")
  results <- data.frame()

  for (outcome in fs_outcomes) {
    if (!outcome %in% colnames(final)) next
    base_df <- final %>% filter(!is.na(.data[[outcome]]))

    for (term in predictors) {
      results <- append_model_result(
        results,
        dataset_label = "ADNI FreeSurfer MRI",
        outcome = outcome,
        prs_term = term,
        model_label = "Original covariate model",
        data = base_df,
        covariates = covariates
      )

      results <- append_model_result(
        results,
        dataset_label = "ADNI FreeSurfer MRI",
        outcome = outcome,
        prs_term = term,
        model_label = "Field-strength adjusted",
        data = base_df,
        covariates = covariates,
        categorical_covariates = "FS_FIELD_STRENGTH",
        notes = "Exact-matched FreeSurfer MRI rows were ADNI1 scans; field strength was 1.5T or 3T."
      )

      results <- append_model_result(
        results,
        dataset_label = "ADNI FreeSurfer MRI",
        outcome = outcome,
        prs_term = term,
        model_label = "1.5T-only restriction",
        data = base_df %>% filter(FS_FIELD_STRENGTH == "1.5T"),
        covariates = covariates,
        notes = "Restriction to the dominant ADNI1 1.5T subgroup."
      )

      for (group in c("Younger (<70)", "Older (>=70)")) {
        gdf <- base_df %>% filter(Age_Group_70_clean == group)

        results <- append_model_result(
          results,
          dataset_label = "ADNI FreeSurfer MRI",
          outcome = outcome,
          prs_term = term,
          model_label = paste("Age-stratified original:", group),
          data = gdf,
          covariates = covariates
        )

        if (length(unique(na.omit(gdf$FS_FIELD_STRENGTH))) > 1) {
          results <- append_model_result(
            results,
            dataset_label = "ADNI FreeSurfer MRI",
            outcome = outcome,
            prs_term = term,
            model_label = paste("Age-stratified field-strength adjusted:", group),
            data = gdf,
            covariates = covariates,
            categorical_covariates = "FS_FIELD_STRENGTH"
          )
        }
      }
    }
  }

  cat("[4/5] Running UCD FLAIR WMH manufacturer/model sensitivity models...\n")
  ucd_df <- stratified %>% filter(!is.na(UCD_Log_WMH))

  for (term in predictors) {
    results <- append_model_result(
      results,
      dataset_label = "ADNI UCD FLAIR WMH",
      outcome = "UCD_Log_WMH",
      prs_term = term,
      model_label = "Original covariate model",
      data = ucd_df,
      covariates = covariates
    )

    results <- append_model_result(
      results,
      dataset_label = "ADNI UCD FLAIR WMH",
      outcome = "UCD_Log_WMH",
      prs_term = term,
      model_label = "Manufacturer adjusted",
      data = ucd_df,
      covariates = covariates,
      categorical_covariates = "UCD_MANUFACTURER",
      notes = "Exact-matched UCD WMH rows were ADNIGO 3T scans."
    )

    results <- append_model_result(
      results,
      dataset_label = "ADNI UCD FLAIR WMH",
      outcome = "UCD_Log_WMH",
      prs_term = term,
      model_label = "Scanner-model adjusted",
      data = ucd_df,
      covariates = covariates,
      categorical_covariates = "UCD_SCANNER_MODEL",
      notes = "Scanner model collapsed to levels with n >= 10; smaller models grouped as Other."
    )

    new_col <- paste0("UCD_Log_WMH_manufacturer_preserved_", term)
    hdf <- residualize_batch_preserving_term(
      data = ucd_df,
      outcome = "UCD_Log_WMH",
      term = term,
      covariates = covariates,
      batch = "UCD_MANUFACTURER",
      new_col = new_col
    )

    results <- append_model_result(
      results,
      dataset_label = "ADNI UCD FLAIR WMH",
      outcome = new_col,
      prs_term = term,
      model_label = "Manufacturer-effect residualized outcome",
      data = hdf,
      covariates = covariates,
      notes = paste(
        "Linear batch-effect removal preserving age, sex, education,",
        "APOE4 and the PRS term of interest."
      )
    )

    for (manufacturer in c("Siemens", "GE")) {
      results <- append_model_result(
        results,
        dataset_label = "ADNI UCD FLAIR WMH",
        outcome = "UCD_Log_WMH",
        prs_term = term,
        model_label = paste0(manufacturer, "-only restriction"),
        data = ucd_df %>% filter(UCD_MANUFACTURER == manufacturer),
        covariates = covariates
      )
    }
  }

  cat("[5/5] Writing results...\n")
  if (nrow(results) > 0) {
    results <- results %>%
      mutate(
        Direction = ifelse(Beta > 0, "positive", "negative"),
        P_formatted = format_p(P)
      ) %>%
      select(Dataset, Outcome, PRS_term, Model, N, Beta, SE, CI95_low,
             CI95_high, P, R2, Covariates, Notes)
  }

  fwrite(results, file.path(output_dir, "ADNI_scanner_aware_sensitivity_results.csv"))

  summary_lines <- c(
    "ADNI scanner-aware sensitivity analysis",
    paste0("FreeSurfer ST128SV exact-matched MRI rows: ", nrow(fs_nonmissing)),
    paste0("FreeSurfer field strength counts: ",
           paste(capture.output(print(table(fs_nonmissing$FS_FIELD_STRENGTH,
                                            useNA = "ifany"))),
                 collapse = " ")),
    paste0("UCD FLAIR WMH exact-matched rows: ", nrow(ucd_nonmissing)),
    paste0("UCD manufacturer counts: ",
           paste(capture.output(print(table(ucd_nonmissing$UCD_MANUFACTURER,
                                            useNA = "ifany"))),
                 collapse = " "))
  )

  writeLines(summary_lines,
             con = file.path(output_dir, "ADNI_scanner_sensitivity_summary.txt"))

  cat("\n=== ADNI Scanner-Aware Sensitivity Complete ===\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  print(sessionInfo())

  invisible(list(
    freesurfer_enriched = final,
    ucd_enriched = stratified,
    batch_counts = batch_counts,
    sensitivity_results = results
  ))
}


# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================
#
# results <- run_adni_scanner_aware_sensitivity(
#   final_analytic_file =
#     "data/ADNI_Age_Stratified_Results/ADNI_Age_Stratified_Final.csv",
#   stratified_analytic_file =
#     "data/ADNI_Validation_Results_Stratified/ADNI_Merged_Data_Stratified.csv",
#   freesurfer_metadata_file =
#     "data/ADNI/data/MRI/MRI_FreeSurfer.csv",
#   ucd_wmh_metadata_file =
#     "data/ADNI/data/MRI/UCD_WMH_03Jan2026.csv",
#   output_dir = "results/adni_scanner_sensitivity"
# )
