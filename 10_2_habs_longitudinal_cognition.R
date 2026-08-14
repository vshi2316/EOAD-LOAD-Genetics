# HABS-HD validation of the temporal structural phenotype against cognition
# The structural score is specified before fitting the cognitive models.

options(stringsAsFactors = FALSE, scipen = 999)
script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)
if (is.na(script_path)) {
  stop("Run with source('10_2_habs_longitudinal_cognition.R').", call. = FALSE)
}

ROOT <- dirname(script_path)
DATA <- Sys.getenv("HABS_DATA_DIR", unset = file.path(ROOT, "data", "HABS_HD"))
RESULTS <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(ROOT, "results"))
INPUT <- file.path(RESULTS, "habs_multimodal", "habs_temporal_chain.tsv")
OUT <- file.path(RESULTS, "habs_longitudinal")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table.", call. = FALSE)
}
library(data.table)

ZS_FILE <- Sys.getenv("HABS_COGNITION_FILE", unset = file.path(DATA, "cognition_zscores.csv"))
MMSE_FILE <- Sys.getenv("HABS_MMSE_FILE", unset = file.path(DATA, "mmse.csv"))
if (!file.exists(INPUT) || !file.exists(ZS_FILE) || !file.exists(MMSE_FILE)) {
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
  m_num[is.na(m_num)] <- match(toupper(substr(m_chr[is.na(m_num)], 1L, 3L)), toupper(month.abb))
  out <- rep(as.IDate(NA), length(y))
  keep <- is.finite(y) & is.finite(m_num) & m_num >= 1L & m_num <= 12L
  out[keep] <- as.IDate(sprintf("%04d-%02d-15", y[keep], m_num[keep]))
  out
}
nearest_record <- function(x, target, before = 180L, after = 180L) {
  if (!nrow(x) || is.na(target)) return(NULL)
  delta <- as.integer(x$date - target)
  eligible <- which(is.finite(delta) & delta >= -before & delta <= after)
  if (!length(eligible)) return(NULL)
  x[eligible[which.min(abs(delta[eligible]))]]
}

chain <- fread(INPUT, na.strings = c("", "NA", "NaN"))
for (field in c("baseline_mri_date", "followup_mri_date")) {
  chain[, (field) := as.IDate(get(field))]
}

# Tests for which higher z scores indicate better performance.
positive_tests <- c(
  "DS_ZSCORE", "LM1_AB_ZSCORE", "LM2_AB_ZSCORE", "FAS_ZSCORE",
  "ANIMAL_ZSCORE", "DIGIT_SYMBOL_SUBSTITUTION_ZSCORE",
  "SEVLT_T1235_ZSCORE", "SEVLT_DR_ZSCORE"
)
# Trails A and B are completion-time z scores and are reversed before averaging.
reverse_tests <- c("TRAILS_A_ZSCORE", "TRAILS_B_ZSCORE")
memory_tests <- c("LM1_AB_ZSCORE", "LM2_AB_ZSCORE", "SEVLT_T1235_ZSCORE", "SEVLT_DR_ZSCORE")
executive_tests <- c(
  "DS_ZSCORE", "FAS_ZSCORE", "ANIMAL_ZSCORE",
  "DIGIT_SYMBOL_SUBSTITUTION_ZSCORE", "TRAILS_A_ZSCORE", "TRAILS_B_ZSCORE"
)

zs <- fread(ZS_FILE, na.strings = c("", "NA", "N/A", "-9999.0", "-9999"))
required_z <- unique(c(positive_tests, reverse_tests))
missing_z <- setdiff(required_z, names(zs))
if (length(missing_z)) stop("Missing cognitive z-score fields: ", paste(missing_z, collapse = ", "), call. = FALSE)
zs[, `:=`(
  Med_ID = as.integer(MED_ID),
  date = month_date(VISIT_YEAR, VISIT_MONTH)
)]
for (field in required_z) zs[, (field) := as_numeric_clean(get(field))]
for (field in reverse_tests) zs[, (field) := -get(field)]

zs[, `:=`(
  cognitive_global = rowMeans(.SD, na.rm = TRUE),
  cognitive_global_n = rowSums(!is.na(.SD))
), .SDcols = required_z]
zs[cognitive_global_n < 4L, cognitive_global := NA_real_]
zs[, `:=`(
  cognitive_memory = rowMeans(.SD, na.rm = TRUE),
  cognitive_memory_n = rowSums(!is.na(.SD))
), .SDcols = memory_tests]
zs[cognitive_memory_n < 2L, cognitive_memory := NA_real_]
zs[, `:=`(
  cognitive_executive = rowMeans(.SD, na.rm = TRUE),
  cognitive_executive_n = rowSums(!is.na(.SD))
), .SDcols = executive_tests]
zs[cognitive_executive_n < 3L, cognitive_executive := NA_real_]
zs <- zs[!is.na(Med_ID) & !is.na(date)]

mmse <- fread(MMSE_FILE, na.strings = c("", "NA", "N/A", "-9999.0", "-9999"))
mmse[, `:=`(
  Med_ID = as.integer(MED_ID),
  date = month_date(VISIT_YEAR, VISIT_MONTH),
  MMSE = as_numeric_clean(MMSE_TOTAL)
)]
mmse <- mmse[is.finite(MMSE) & MMSE >= 0 & MMSE <= 30 & !is.na(date)]

# Match cognition to the prespecified MRI dates. Baseline cognition is within 180
# days of baseline MRI. Follow-up cognition occurs 90 to 730 days after the
# selected follow-up MRI.
matched <- vector("list", nrow(chain))
for (i in seq_len(nrow(chain))) {
  sid <- chain$Med_ID[i]
  baseline_date <- chain$baseline_mri_date[i]
  followup_date <- chain$followup_mri_date[i]
  b <- nearest_record(zs[Med_ID == sid], baseline_date, before = 180L, after = 180L)
  bm <- nearest_record(mmse[Med_ID == sid], baseline_date, before = 180L, after = 180L)
  f <- if (is.na(followup_date)) NULL else nearest_record(zs[Med_ID == sid], followup_date, before = -90L, after = 730L)
  fm <- if (is.na(followup_date)) NULL else nearest_record(mmse[Med_ID == sid], followup_date, before = -90L, after = 730L)
  matched[[i]] <- data.table(
    Med_ID = sid,
    baseline_cog_date = if (is.null(b)) as.IDate(NA) else b$date,
    followup_cog_date = if (is.null(f)) as.IDate(NA) else f$date,
    baseline_global = if (is.null(b)) NA_real_ else b$cognitive_global,
    followup_global = if (is.null(f)) NA_real_ else f$cognitive_global,
    baseline_memory = if (is.null(b)) NA_real_ else b$cognitive_memory,
    followup_memory = if (is.null(f)) NA_real_ else f$cognitive_memory,
    baseline_executive = if (is.null(b)) NA_real_ else b$cognitive_executive,
    followup_executive = if (is.null(f)) NA_real_ else f$cognitive_executive,
    baseline_mmse = if (is.null(bm)) NA_real_ else bm$MMSE,
    followup_mmse = if (is.null(fm)) NA_real_ else fm$MMSE
  )
}
matched <- rbindlist(matched, fill = TRUE)
drop_fields <- intersect(names(chain), setdiff(names(matched), "Med_ID"))
if (length(drop_fields)) chain[, (drop_fields) := NULL]
analysis <- merge(chain, matched, by = "Med_ID", all.x = TRUE)

analysis[, `:=`(
  baseline_global_z = safe_z(baseline_global),
  followup_global_z = safe_z(followup_global),
  baseline_memory_z = safe_z(baseline_memory),
  followup_memory_z = safe_z(followup_memory),
  baseline_executive_z = safe_z(baseline_executive),
  followup_executive_z = safe_z(followup_executive),
  baseline_mmse_z = safe_z(baseline_mmse),
  followup_mmse_z = safe_z(followup_mmse),
  cognitive_followup_years = as.numeric(followup_cog_date - baseline_cog_date) / 365.25,
  diagnosis_factor = factor(ifelse(is.finite(diagnosis) & diagnosis >= 0, diagnosis, NA)),
  scanner = factor(scanner), sex = factor(sex), ethnicity = factor(ethnicity)
)]

hc3_table <- function(fit, model_id, focal_terms) {
  beta <- coef(fit)
  keep <- is.finite(beta)
  X <- model.matrix(fit)[, keep, drop = FALSE]
  beta <- beta[keep]
  e <- residuals(fit) / pmax(1 - hatvalues(fit), 1e-8)
  bread <- tryCatch(
    solve(crossprod(X)),
    error = function(err) qr.solve(crossprod(X), diag(ncol(X)), tol = 1e-10)
  )
  variance <- bread %*% crossprod(X, X * e^2) %*% bread
  se <- sqrt(pmax(diag(variance), 0))
  statistic <- beta / se
  data.table(
    model = model_id, term = names(beta), estimate = as.numeric(beta),
    standard_error_hc3 = se, statistic = statistic,
    p_value = 2 * pnorm(-abs(statistic)),
    conf_low = as.numeric(beta) - qnorm(.975) * se,
    conf_high = as.numeric(beta) + qnorm(.975) * se,
    n = nobs(fit), focal = names(beta) %in% focal_terms
  )
}
fit_model <- function(data, outcome, exposure, covariates, model_id) {
  variables <- unique(c(outcome, exposure, covariates))
  d <- copy(data[which(complete.cases(data[, variables, with = FALSE]))])
  covariates <- covariates[vapply(covariates, function(v) uniqueN(d[[v]]) > 1L, logical(1))]
  formula <- reformulate(c(exposure, covariates), response = outcome)
  if (nrow(d) <= length(attr(terms(formula), "term.labels")) + 2L) {
    return(list(
      table = data.table(model = model_id, term = NA_character_, estimate = NA_real_,
                         standard_error_hc3 = NA_real_, statistic = NA_real_, p_value = NA_real_,
                         conf_low = NA_real_, conf_high = NA_real_, n = nrow(d), focal = NA),
      manifest = data.table(model = model_id, formula = paste(deparse(formula), collapse = " "),
                            n = nrow(d), status = "skipped_insufficient_cases")
    ))
  }
  fit <- lm(formula, data = d)
  list(
    table = hc3_table(fit, model_id, exposure),
    manifest = data.table(model = model_id, formula = paste(deparse(formula), collapse = " "),
                          n = nrow(d), status = "completed")
  )
}

base_cov <- c("age_z", "sex", "ethnicity", "apoe4", "scanner", "diagnosis_factor")
specifications <- list(
  list(y = "baseline_global_z", x = "structure_z",
       cov = c("tau_z", "amyloid_z", base_cov), id = "cross_sectional_global_cognition"),
  list(y = "followup_global_z", x = "structure_change_z",
       cov = c("baseline_global_z", base_cov, "mri_interval_years", "cognitive_followup_years"),
       id = "longitudinal_global_cognition"),
  list(y = "followup_memory_z", x = "structure_change_z",
       cov = c("baseline_memory_z", base_cov, "mri_interval_years", "cognitive_followup_years"),
       id = "longitudinal_memory_sensitivity"),
  list(y = "followup_executive_z", x = "structure_change_z",
       cov = c("baseline_executive_z", base_cov, "mri_interval_years", "cognitive_followup_years"),
       id = "longitudinal_executive_sensitivity"),
  list(y = "followup_mmse_z", x = "structure_change_z",
       cov = c("baseline_mmse_z", base_cov, "mri_interval_years"),
       id = "longitudinal_mmse_sensitivity")
)

results <- list(); manifests <- list()
for (spec in specifications) {
  fitted <- fit_model(analysis, spec$y, spec$x, spec$cov, spec$id)
  results[[length(results) + 1L]] <- fitted$table
  manifests[[length(manifests) + 1L]] <- fitted$manifest
}
effects <- rbindlist(results, fill = TRUE)
effects[focal == TRUE, p_fdr_across_cognitive_models := p.adjust(p_value, "BH")]
manifest <- rbindlist(manifests, fill = TRUE)

# Empirical direction check against MMSE within the same matched visits.
direction_qc <- data.table(
  comparison = c("baseline_global_vs_MMSE", "followup_global_vs_MMSE"),
  n = c(
    sum(complete.cases(analysis[, .(baseline_global, baseline_mmse)])),
    sum(complete.cases(analysis[, .(followup_global, followup_mmse)]))
  ),
  spearman_rho = c(
    suppressWarnings(cor(analysis$baseline_global, analysis$baseline_mmse, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(analysis$followup_global, analysis$followup_mmse, method = "spearman", use = "complete.obs"))
  )
)

write_tsv(analysis, "habs_oriented_cognition_longitudinal_dataset.tsv")
write_tsv(effects, "habs_cognition_model_effects_hc3.tsv")
write_tsv(manifest, "habs_cognition_model_manifest.tsv")
write_tsv(direction_qc, "habs_cognitive_direction_qc.tsv")
capture.output({
  cat("HABS-HD cognition and longitudinal structural validation\n\n")
  cat("Cognitive direction quality control\n")
  print(direction_qc)
  cat("\nModel manifest\n")
  print(manifest)
  cat("\nFocal effects\n")
  print(effects[focal == TRUE])
}, file = file.path(OUT, "analysis_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("\nHABS-HD cognitive validation completed.\nResults: ", OUT, "\n", sep = "")
