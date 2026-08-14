# HABS-HD structural and cognition models
# Primary covariates: age, sex, education, ethnicity, APOE4, and MRI scanner.
# Diagnosis is added only in prespecified sensitivity models.

options(stringsAsFactors = FALSE, scipen = 999)
script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)
if (is.na(script_path)) {
  stop("Run with source('10_3_habs_cognition_models.R').", call. = FALSE)
}

ROOT <- dirname(script_path)
DATA <- Sys.getenv("HABS_DATA_DIR", unset = file.path(ROOT, "data", "HABS_HD"))
RESULTS <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(ROOT, "results"))
INPUT <- file.path(RESULTS, "habs_longitudinal", "habs_oriented_cognition_longitudinal_dataset.tsv")
EDUCATION_FILE <- Sys.getenv("HABS_EDUCATION_FILE", unset = file.path(DATA, "education.csv"))
OUT <- file.path(RESULTS, "habs_cognition")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table.", call. = FALSE)
}
library(data.table)

if (!file.exists(INPUT) || !file.exists(EDUCATION_FILE)) {
  stop("Required HABS-HD inputs are missing.", call. = FALSE)
}

write_tsv <- function(x, name) {
  fwrite(x, file.path(OUT, name), sep = "\t", quote = FALSE, na = "NA")
}
as_numeric_clean <- function(x) {
  x <- trimws(gsub(",", "", as.character(x), fixed = TRUE))
  x[x %in% c("", "NA", "N/A", ".", "-9999", "-8888", "-777777")] <- NA_character_
  suppressWarnings(as.numeric(x))
}
safe_z <- function(x) {
  x <- as_numeric_clean(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
month_date <- function(year, month) {
  y <- suppressWarnings(as.integer(as.character(year)))
  m_chr <- trimws(as.character(month))
  m_num <- suppressWarnings(as.integer(m_chr))
  missing_month <- is.na(m_num)
  m_num[missing_month] <- match(
    toupper(substr(m_chr[missing_month], 1L, 3L)),
    toupper(month.abb)
  )
  out <- rep(as.IDate(NA), length(y))
  keep <- is.finite(y) & is.finite(m_num) & m_num >= 1L & m_num <= 12L
  out[keep] <- as.IDate(sprintf("%04d-%02d-15", y[keep], m_num[keep]))
  out
}

analysis <- fread(INPUT, na.strings = c("", "NA", "NaN"))
analysis[, `:=`(
  baseline_mri_date = as.IDate(baseline_mri_date),
  followup_mri_date = as.IDate(followup_mri_date)
)]

education <- fread(
  EDUCATION_FILE,
  na.strings = c("", "NA", "N/A", "-9999", "-8888", "-777777")
)
education[, `:=`(
  Med_ID = as.integer(MED_ID),
  date = month_date(VISIT_YEAR, VISIT_MONTH),
  education_years = as_numeric_clean(ID_EDUCATION)
)]
# Values outside a plausible completed-education range are treated as missing.
education[education_years < 0 | education_years > 30, education_years := NA_real_]
education <- education[!is.na(Med_ID) & !is.na(date) & is.finite(education_years)]
setorder(education, Med_ID, date)

# Education is matched to the closest record within five years of baseline MRI.
# It is a stable historical covariate, so the wider window avoids discarding
# otherwise complete imaging-pathology records.
education_match <- vector("list", nrow(analysis))
for (i in seq_len(nrow(analysis))) {
  sid <- analysis$Med_ID[i]
  target <- analysis$baseline_mri_date[i]
  candidates <- education[Med_ID == sid]
  if (!nrow(candidates) || is.na(target)) {
    education_match[[i]] <- data.table(
      Med_ID = sid, education_years = NA_real_, education_date = as.IDate(NA),
      education_gap_days = NA_integer_
    )
    next
  }
  gap <- abs(as.integer(candidates$date - target))
  eligible <- which(is.finite(gap) & gap <= 1826L)
  if (!length(eligible)) {
    education_match[[i]] <- data.table(
      Med_ID = sid, education_years = NA_real_, education_date = as.IDate(NA),
      education_gap_days = NA_integer_
    )
    next
  }
  selected <- candidates[eligible[which.min(gap[eligible])]]
  education_match[[i]] <- data.table(
    Med_ID = sid,
    education_years = selected$education_years,
    education_date = selected$date,
    education_gap_days = min(gap[eligible])
  )
}
education_match <- rbindlist(education_match, fill = TRUE)
analysis <- merge(analysis, education_match, by = "Med_ID", all.x = TRUE)
analysis[, education_z := safe_z(education_years)]

# Apply the prespecified coding and exclude special missing-value levels.
analysis[, `:=`(
  sex = factor(ifelse(as.character(sex) %in% c("0", "1"), as.character(sex), NA)),
  ethnicity = factor(ifelse(
    as.character(ethnicity) %in% c("Black", "Hispanic", "White"),
    as.character(ethnicity), NA
  )),
  scanner = factor(ifelse(
    as.character(scanner) %in% c("Skyra", "Vida1", "Vida1a", "Vida2"),
    as.character(scanner), NA
  )),
  diagnosis_factor = factor(ifelse(
    is.finite(as_numeric_clean(diagnosis)) & as_numeric_clean(diagnosis) %in% c(0, 1, 2),
    as.character(as_numeric_clean(diagnosis)), NA
  )),
  apoe4 = as_numeric_clean(apoe4)
)]

hc3_table <- function(fit, model_id, focal_terms, model_role) {
  beta <- coef(fit)
  keep <- is.finite(beta)
  X <- model.matrix(fit)[, keep, drop = FALSE]
  beta <- beta[keep]
  adjusted_residual <- residuals(fit) / pmax(1 - hatvalues(fit), 1e-8)
  bread <- tryCatch(
    solve(crossprod(X)),
    error = function(err) qr.solve(crossprod(X), diag(ncol(X)), tol = 1e-10)
  )
  variance <- bread %*% crossprod(X, X * adjusted_residual^2) %*% bread
  se <- sqrt(pmax(diag(variance), 0))
  statistic <- beta / se
  data.table(
    model = model_id,
    role = model_role,
    term = names(beta),
    estimate = as.numeric(beta),
    standard_error_hc3 = se,
    statistic = statistic,
    p_value = 2 * pnorm(-abs(statistic)),
    conf_low = as.numeric(beta) - qnorm(.975) * se,
    conf_high = as.numeric(beta) + qnorm(.975) * se,
    n = nobs(fit),
    focal = names(beta) %in% focal_terms
  )
}

fit_model <- function(data, outcome, exposure, covariates, model_id, model_role) {
  variables <- unique(c(outcome, exposure, covariates))
  complete_rows <- which(complete.cases(data[, variables, with = FALSE]))
  d <- copy(data[complete_rows])
  covariates <- covariates[vapply(
    covariates,
    function(v) uniqueN(d[[v]]) > 1L,
    logical(1)
  )]
  formula <- reformulate(c(exposure, covariates), response = outcome)
  parameter_count <- length(attr(terms(formula), "term.labels")) + 1L
  if (!nrow(d) || nrow(d) <= parameter_count + 1L) {
    return(list(
      table = data.table(
        model = model_id, role = model_role, term = NA_character_,
        estimate = NA_real_, standard_error_hc3 = NA_real_, statistic = NA_real_,
        p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
        n = nrow(d), focal = NA
      ),
      manifest = data.table(
        model = model_id, role = model_role,
        formula = paste(deparse(formula), collapse = " "),
        n = nrow(d), status = "skipped_insufficient_cases"
      )
    ))
  }
  fit <- tryCatch(lm(formula, data = d), error = function(err) err)
  if (inherits(fit, "error")) {
    return(list(
      table = data.table(
        model = model_id, role = model_role, term = NA_character_,
        estimate = NA_real_, standard_error_hc3 = NA_real_, statistic = NA_real_,
        p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
        n = nrow(d), focal = NA
      ),
      manifest = data.table(
        model = model_id, role = model_role,
        formula = paste(deparse(formula), collapse = " "),
        n = nrow(d), status = paste0("skipped_lm_error: ", conditionMessage(fit))
      )
    ))
  }
  list(
    table = hc3_table(fit, model_id, exposure, model_role),
    manifest = data.table(
      model = model_id, role = model_role,
      formula = paste(deparse(formula), collapse = " "),
      n = nrow(d), status = "completed"
    )
  )
}

primary_covariates <- c(
  "age_z", "sex", "education_z", "ethnicity", "apoe4", "scanner"
)
diagnosis_sensitivity_covariates <- c(primary_covariates, "diagnosis_factor")

specifications <- list(
  # Confirmatory cross-sectional structural-cognition model.
  list(
    y = "baseline_global_z", x = "structure_z",
    cov = c("tau_z", "amyloid_z", primary_covariates),
    id = "primary_cross_sectional_global_cognition",
    role = "primary"
  ),
  # Same model with cognitive diagnosis to quantify possible overadjustment.
  list(
    y = "baseline_global_z", x = "structure_z",
    cov = c("tau_z", "amyloid_z", diagnosis_sensitivity_covariates),
    id = "sensitivity_cross_sectional_with_diagnosis",
    role = "sensitivity"
  ),
  # Longitudinal component models do not require amyloid or tau measurements.
  list(
    y = "followup_global_z", x = "structure_change_z",
    cov = c(
      "baseline_global_z", primary_covariates,
      "mri_interval_years", "cognitive_followup_years"
    ),
    id = "primary_longitudinal_global_cognition",
    role = "primary"
  ),
  list(
    y = "followup_global_z", x = "structure_change_z",
    cov = c(
      "baseline_global_z", diagnosis_sensitivity_covariates,
      "mri_interval_years", "cognitive_followup_years"
    ),
    id = "sensitivity_longitudinal_with_diagnosis",
    role = "sensitivity"
  ),
  list(
    y = "followup_memory_z", x = "structure_change_z",
    cov = c(
      "baseline_memory_z", primary_covariates,
      "mri_interval_years", "cognitive_followup_years"
    ),
    id = "sensitivity_longitudinal_memory",
    role = "sensitivity"
  ),
  list(
    y = "followup_executive_z", x = "structure_change_z",
    cov = c(
      "baseline_executive_z", primary_covariates,
      "mri_interval_years", "cognitive_followup_years"
    ),
    id = "sensitivity_longitudinal_executive",
    role = "sensitivity"
  ),
  list(
    y = "followup_mmse_z", x = "structure_change_z",
    cov = c(
      "baseline_mmse_z", primary_covariates,
      "mri_interval_years"
    ),
    id = "sensitivity_longitudinal_mmse",
    role = "sensitivity"
  )
)

results <- list()
manifests <- list()
for (spec in specifications) {
  fitted <- fit_model(
    analysis, spec$y, spec$x, spec$cov, spec$id, spec$role
  )
  results[[length(results) + 1L]] <- fitted$table
  manifests[[length(manifests) + 1L]] <- fitted$manifest
}
effects <- rbindlist(results, fill = TRUE)
manifest <- rbindlist(manifests, fill = TRUE)

# Multiplicity is controlled separately within primary and sensitivity families.
effects[focal == TRUE & role == "primary",
        p_fdr_primary := p.adjust(p_value, method = "BH")]
effects[focal == TRUE & role == "sensitivity",
        p_fdr_sensitivity := p.adjust(p_value, method = "BH")]

education_audit <- data.table(
  metric = c(
    "MRI participants", "Education matched", "Education missing",
    "Median education years", "Education IQR lower", "Education IQR upper"
  ),
  value = c(
    nrow(analysis),
    sum(is.finite(analysis$education_years)),
    sum(!is.finite(analysis$education_years)),
    median(analysis$education_years, na.rm = TRUE),
    quantile(analysis$education_years, 0.25, na.rm = TRUE),
    quantile(analysis$education_years, 0.75, na.rm = TRUE)
  )
)

write_tsv(analysis, "habs_cognition_analysis_dataset.tsv")
write_tsv(effects, "habs_cognition_effects_hc3.tsv")
write_tsv(manifest, "habs_cognition_model_manifest.tsv")
write_tsv(education_audit, "habs_education_matching_audit.tsv")

capture.output({
  cat("HABS-HD structural and cognition models\n\n")
  cat("Education matching audit\n")
  print(education_audit)
  cat("\nModel manifest\n")
  print(manifest)
  cat("\nFocal effects\n")
  print(effects[focal == TRUE])
}, file = file.path(OUT, "analysis_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))

cat("\nHABS-HD cognition models completed.\nResults: ", OUT, "\n", sep = "")
