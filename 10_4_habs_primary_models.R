# HABS-HD primary multimodal analyses.
# This is the single analysis entry point for the HABS results.
# Primary evidence: amyloid-to-tau, tau-to-temporal structure, and
# temporal structure-to-concurrent cognition. Longitudinal analyses are
# prespecified secondary analyses.

options(stringsAsFactors = FALSE, scipen = 999)
ofile <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
if (is.na(ofile)) stop("Run with source('10_4_habs_primary_models.R').", call. = FALSE)
ROOT <- dirname(ofile)
DATA <- Sys.getenv("HABS_DATA_DIR", unset = file.path(ROOT, "data", "HABS_HD"))
RESULTS <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(ROOT, "results"))
CHAIN_FILE <- file.path(RESULTS, "habs_multimodal", "habs_temporal_chain.tsv")
COG_FILE <- file.path(RESULTS, "habs_cognition", "habs_cognition_analysis_dataset.tsv")
OUT <- file.path(RESULTS, "habs_primary_models")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("readxl", quietly = TRUE)) stop("Install data.table and readxl.", call. = FALSE)
library(data.table)
if (!file.exists(CHAIN_FILE) || !file.exists(COG_FILE)) stop("Required prior HABS outputs are missing.", call. = FALSE)

write_tsv <- function(x, name) fwrite(x, file.path(OUT, name), sep = "\t", quote = FALSE, na = "NA")
as_numeric_clean <- function(x) {
  x <- trimws(gsub(",", "", as.character(x), fixed = TRUE))
  x[x %in% c("", "NA", "N/A", ".", "-9999", "-8888", "-777777")] <- NA_character_
  suppressWarnings(as.numeric(x))
}
safe_z <- function(x) {
  x <- as_numeric_clean(x); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
month_date <- function(year, month) {
  y <- suppressWarnings(as.integer(as.character(year)))
  m <- suppressWarnings(as.integer(as.character(month)))
  out <- rep(as.IDate(NA), length(y)); ok <- is.finite(y) & is.finite(m) & m >= 1 & m <= 12
  out[ok] <- as.IDate(sprintf("%04d-%02d-15", y[ok], m[ok])); out
}
chain <- fread(CHAIN_FILE, na.strings = c("", "NA", "NaN"))
cog <- fread(COG_FILE, na.strings = c("", "NA", "NaN"), check.names = FALSE)

# Remove the obsolete upper-case MMSE columns before merging. The lower-case
# fields are the uniquely named, oriented fields used by the final models.
cog <- cog[, which(!duplicated(names(cog))), with = FALSE]
drop_old <- intersect(c("baseline_MMSE", "followup_MMSE"), names(cog))
if (length(drop_old)) cog[, (drop_old) := NULL]
chain <- chain[, which(!duplicated(names(chain))), with = FALSE]
analysis <- copy(cog)

to_factor <- function(x, allowed = NULL) {
  x <- as.character(x)
  if (!is.null(allowed)) x[!x %in% allowed] <- NA_character_
  factor(x)
}
analysis[, `:=`(
  sex = to_factor(sex),
  ethnicity = to_factor(ethnicity, c("Black", "Hispanic", "White")),
  scanner = to_factor(scanner, c("Skyra", "Vida1", "Vida1a", "Vida2")),
  apoe4 = suppressWarnings(as.numeric(apoe4)),
  age_z = suppressWarnings(as.numeric(age_z)),
  education_z = suppressWarnings(as.numeric(education_z)),
  diagnosis_factor = to_factor(diagnosis, c("0", "1", "2"))
)]

# First eligible follow-up MRI sensitivity, preserving the prespecified eight-region
# score and fold-specific standardization used in the primary analysis.
MRI_FILE <- Sys.getenv("HABS_CORTICAL_THICKNESS_FILE", unset = file.path(DATA, "cortical_thickness.xlsx"))
temporal_fields <- c(
  "L_entorhinal_thickavg", "R_entorhinal_thickavg",
  "L_parahippocampal_thickavg", "R_parahippocampal_thickavg",
  "L_inferiortemporal_thickavg", "R_inferiortemporal_thickavg",
  "L_middletemporal_thickavg", "R_middletemporal_thickavg"
)
mri_sheets <- readxl::excel_sheets(MRI_FILE)
mri_sheets <- mri_sheets[grepl("inclusion", tolower(mri_sheets)) & !grepl("exclusion|pending", tolower(mri_sheets))]
mri <- rbindlist(lapply(mri_sheets, function(sheet) {
  x <- as.data.table(suppressWarnings(readxl::read_excel(MRI_FILE, sheet = sheet)))
  x[, source_sheet := sheet]; x
}), fill = TRUE)
mri[, `:=`(Med_ID = as.integer(Med_ID), date = month_date(MRI_year, MRI_month))]
for (v in temporal_fields) mri[, (v) := as_numeric_clean(get(v))]
mri <- mri[!is.na(Med_ID) & !is.na(date) & complete.cases(mri[, temporal_fields, with = FALSE])]
setorder(mri, Med_ID, date)
analysis[, `:=`(first_followup_structure = NA_real_, first_followup_mri_date = as.IDate(NA))]
for (fold_id in sort(unique(analysis$score_fold))) {
  train_ids <- analysis[score_fold != fold_id, Med_ID]
  training <- mri[Med_ID %in% train_ids]
  mu <- vapply(temporal_fields, function(v) mean(training[[v]], na.rm = TRUE), numeric(1))
  sigma <- vapply(temporal_fields, function(v) sd(training[[v]], na.rm = TRUE), numeric(1))
  sigma[!is.finite(sigma) | sigma <= 0] <- 1
  for (i in which(analysis$score_fold == fold_id)) {
    sid <- analysis$Med_ID[i]; baseline_date <- as.IDate(analysis$baseline_mri_date[i])
    followup <- mri[Med_ID == sid & date >= baseline_date + 365L & date <= baseline_date + 1825L]
    if (!nrow(followup)) next
    followup <- followup[which.min(date)]
    values <- as_numeric_clean(unlist(followup[, temporal_fields, with = FALSE], use.names = FALSE))
    analysis[i, `:=`(
      first_followup_structure = -mean((values - mu) / sigma),
      first_followup_mri_date = followup$date
    )]
  }
}
analysis[, first_mri_interval_years := as.numeric(first_followup_mri_date - as.IDate(baseline_mri_date)) / 365.25]
analysis[, first_structure_change := (first_followup_structure - baseline_structure) / first_mri_interval_years]
analysis[, first_structure_change_z := safe_z(first_structure_change)]

# Match cognition specifically to the first eligible follow-up MRI.
ZS_FILE <- Sys.getenv("HABS_COGNITION_FILE", unset = file.path(DATA, "cognition_zscores.csv"))
positive_tests <- c(
  "DS_ZSCORE", "LM1_AB_ZSCORE", "LM2_AB_ZSCORE", "FAS_ZSCORE",
  "ANIMAL_ZSCORE", "DIGIT_SYMBOL_SUBSTITUTION_ZSCORE",
  "SEVLT_T1235_ZSCORE", "SEVLT_DR_ZSCORE"
)
reverse_tests <- c("TRAILS_A_ZSCORE", "TRAILS_B_ZSCORE")
cognitive_fields <- c(positive_tests, reverse_tests)
zs <- fread(ZS_FILE, na.strings = c("", "NA", "N/A", "-9999", "-9999.0"))
zs[, `:=`(Med_ID = as.integer(MED_ID), date = {
  month_number <- match(toupper(substr(trimws(VISIT_MONTH), 1L, 3L)), toupper(month.abb))
  month_date(VISIT_YEAR, month_number)
})]
for (v in cognitive_fields) zs[, (v) := as_numeric_clean(get(v))]
for (v in reverse_tests) zs[, (v) := -get(v)]
zs[, `:=`(
  oriented_global = rowMeans(.SD, na.rm = TRUE),
  oriented_global_n = rowSums(!is.na(.SD))
), .SDcols = cognitive_fields]
zs[oriented_global_n < 4L, oriented_global := NA_real_]
zs <- zs[is.finite(oriented_global) & !is.na(date)]
analysis[, `:=`(first_followup_global = NA_real_, first_followup_cog_date = as.IDate(NA))]
for (i in which(!is.na(analysis$first_followup_mri_date))) {
  sid <- analysis$Med_ID[i]; anchor <- analysis$first_followup_mri_date[i]
  candidates <- zs[Med_ID == sid & date >= anchor + 90L & date <= anchor + 730L]
  if (!nrow(candidates)) next
  selected <- candidates[which.min(date)]
  analysis[i, `:=`(
    first_followup_global = selected$oriented_global,
    first_followup_cog_date = selected$date
  )]
}
analysis[, first_followup_global_z := safe_z(first_followup_global)]
analysis[, first_cognitive_followup_years := as.numeric(first_followup_cog_date - as.IDate(baseline_cog_date)) / 365.25]

hc3 <- function(fit, model, role, focal) {
  b <- coef(fit); ok <- is.finite(b); b <- b[ok]
  X <- model.matrix(fit)[, ok, drop = FALSE]
  e <- residuals(fit) / pmax(1 - hatvalues(fit), 1e-8)
  bread <- tryCatch(solve(crossprod(X)), error = function(z) qr.solve(crossprod(X), diag(ncol(X))))
  se <- sqrt(pmax(diag(bread %*% crossprod(X, X * e^2) %*% bread), 0))
  st <- b / se
  data.table(model = model, role = role, term = names(b), estimate = as.numeric(b),
             standard_error_hc3 = se, statistic = st, p_value = 2 * pnorm(-abs(st)),
             conf_low = as.numeric(b) - qnorm(.975) * se,
             conf_high = as.numeric(b) + qnorm(.975) * se,
             n = nobs(fit), focal = names(b) %in% focal)
}
run_model <- function(d, y, x, cov, id, role) {
  vars <- unique(c(y, x, cov)); keep <- which(complete.cases(d[, vars, with = FALSE])); d <- d[keep]
  cov <- cov[vapply(cov, function(v) uniqueN(d[[v]]) > 1L, logical(1))]
  f <- reformulate(c(x, cov), response = y)
  if (nrow(d) <= length(attr(terms(f), "term.labels")) + 2L) {
    return(list(table = data.table(model = id, role = role, term = NA_character_,
      estimate = NA_real_, standard_error_hc3 = NA_real_, statistic = NA_real_,
      p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_, n = nrow(d), focal = NA),
      manifest = data.table(model = id, role = role, formula = paste(deparse(f), collapse = " "),
        n = nrow(d), status = "skipped_insufficient_cases")))
  }
  fit <- tryCatch(lm(f, data = d), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(table = data.table(model = id, role = role, term = NA_character_,
      estimate = NA_real_, standard_error_hc3 = NA_real_, statistic = NA_real_,
      p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_, n = nrow(d), focal = NA),
      manifest = data.table(model = id, role = role, formula = paste(deparse(f), collapse = " "),
        n = nrow(d), status = paste0("skipped_lm_error: ", conditionMessage(fit)))))
  }
  list(table = hc3(fit, id, role, x),
       manifest = data.table(model = id, role = role, formula = paste(deparse(f), collapse = " "),
         n = nrow(d), status = "completed"))
}

base_cov <- c("age_z", "sex", "education_z", "ethnicity", "apoe4", "scanner")
primary <- list(
  list(d = analysis, y = "tau_z", x = "amyloid_z", cov = base_cov, id = "primary_amyloid_to_tau", role = "primary"),
  list(d = analysis, y = "structure_z", x = "tau_z", cov = c("amyloid_z", base_cov), id = "primary_tau_to_structure", role = "primary"),
  list(d = analysis, y = "baseline_global_z", x = "structure_z", cov = c("tau_z", "amyloid_z", base_cov), id = "primary_structure_to_concurrent_cognition", role = "primary")
)
secondary <- list(
  list(d = analysis, y = "baseline_global_z", x = "structure_z", cov = base_cov, id = "secondary_structure_total_association", role = "secondary"),
  list(d = analysis, y = "baseline_global_z", x = "structure_z", cov = c("tau_z", "amyloid_z", base_cov, "diagnosis_factor"), id = "secondary_structure_with_diagnosis", role = "secondary"),
  list(d = analysis, y = "followup_global_z", x = "structure_change_z", cov = c("baseline_global_z", base_cov, "mri_interval_years", "cognitive_followup_years"), id = "secondary_longitudinal_global_cognition", role = "secondary"),
  list(d = analysis, y = "followup_memory_z", x = "structure_change_z", cov = c("baseline_memory_z", base_cov, "mri_interval_years", "cognitive_followup_years"), id = "secondary_longitudinal_memory", role = "secondary"),
  list(d = analysis, y = "followup_executive_z", x = "structure_change_z", cov = c("baseline_executive_z", base_cov, "mri_interval_years", "cognitive_followup_years"), id = "secondary_longitudinal_executive", role = "secondary"),
  list(d = analysis, y = "followup_mmse_z", x = "structure_change_z", cov = c("baseline_mmse_z", base_cov, "mri_interval_years"), id = "secondary_longitudinal_mmse", role = "secondary"),
  list(d = analysis, y = "first_followup_global_z", x = "first_structure_change_z", cov = c("baseline_global_z", base_cov, "first_mri_interval_years", "first_cognitive_followup_years"), id = "secondary_first_followup_mri_global_cognition", role = "secondary")
)
all_specs <- c(primary, secondary); tabs <- list(); mans <- list()
for (s in all_specs) { z <- run_model(s$d, s$y, s$x, s$cov, s$id, s$role); tabs[[length(tabs) + 1L]] <- z$table; mans[[length(mans) + 1L]] <- z$manifest }
effects <- rbindlist(tabs, fill = TRUE); manifest <- rbindlist(mans, fill = TRUE)
effects[focal == TRUE & role == "primary", p_fdr_primary := p.adjust(p_value, "BH")]
effects[focal == TRUE & role == "secondary", p_fdr_secondary := p.adjust(p_value, "BH")]

# Amyloid-to-tau time-window sensitivity uses only the two pathology dates.
time_sens <- list()
for (w in c(90L, 180L, 365L)) {
  direct_pet_gap <- abs(as.integer(as.IDate(analysis$baseline_tau_date) - as.IDate(analysis$baseline_amyloid_date)))
  d <- analysis[is.finite(direct_pet_gap) & direct_pet_gap <= w]
  z <- run_model(d, "tau_z", "amyloid_z", base_cov, paste0("secondary_amyloid_to_tau_", w, "d"), "secondary")
  time_sens[[length(time_sens) + 1L]] <- z$table
}
time_sens <- rbindlist(time_sens, fill = TRUE)

# Scanner-stratified models are run only with at least 100 complete cases.
scanner_sens <- list(); scanner_audit <- list()
for (sc in levels(droplevels(analysis$scanner))) {
  d <- analysis[scanner == sc]
  vars <- c("tau_z", "amyloid_z", setdiff(base_cov, "scanner"))
  n_complete <- sum(complete.cases(d[, vars, with = FALSE]))
  scanner_audit[[length(scanner_audit) + 1L]] <- data.table(
    scanner = sc, complete_cases = n_complete,
    status = if (n_complete >= 100L) "eligible" else "excluded_fewer_than_100_complete_cases"
  )
  if (n_complete < 100L) next
  z <- run_model(d, "tau_z", "amyloid_z", setdiff(base_cov, "scanner"),
                 paste0("secondary_amyloid_to_tau_scanner_", sc), "secondary")
  scanner_sens[[length(scanner_sens) + 1L]] <- z$table
}
scanner_sens <- rbindlist(scanner_sens, fill = TRUE)
scanner_audit <- rbindlist(scanner_audit, fill = TRUE)

audit <- data.table(
  metric = c("chain_rows", "analysis_rows", "duplicate_column_names_after_cleaning", "primary_models_completed"),
  value = c(nrow(chain), nrow(analysis), sum(duplicated(names(analysis))), sum(manifest[role == "primary", status == "completed"]))
)
write_tsv(analysis, "habs_primary_analysis_dataset.tsv")
write_tsv(effects, "habs_primary_model_effects_hc3.tsv")
write_tsv(manifest, "habs_primary_model_manifest.tsv")
write_tsv(time_sens, "habs_amyloid_tau_time_sensitivity.tsv")
write_tsv(scanner_sens, "habs_scanner_sensitivity.tsv")
write_tsv(scanner_audit, "habs_scanner_audit.tsv")
write_tsv(audit, "habs_primary_analysis_audit.tsv")
capture.output({
  cat("HABS-HD primary multimodal analysis\n\n"); print(audit)
  cat("\nPrimary focal effects\n"); print(effects[focal == TRUE & role == "primary"])
  cat("\nSecondary focal effects\n"); print(effects[focal == TRUE & role == "secondary"])
  cat("\nAmyloid-to-tau time-window sensitivity\n"); print(time_sens[focal == TRUE])
  cat("\nScanner audit\n"); print(scanner_audit)
  cat("\nEligible scanner focal effects\n"); print(scanner_sens[focal == TRUE])
}, file = file.path(OUT, "analysis_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("\nHABS-HD primary multimodal analysis completed.\nResults: ", OUT, "\n", sep = "")
