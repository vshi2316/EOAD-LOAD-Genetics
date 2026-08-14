# HABS-HD multimodal amyloid-tau-structure-cognition cascade
# HABS is analyzed independently. No ADNI structural weights are used.

options(stringsAsFactors = FALSE, scipen = 999)
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) stop("Run with source('10_1_habs_multimodal_data.R').", call. = FALSE)
ROOT <- dirname(script_path)
DATA <- Sys.getenv("HABS_DATA_DIR", unset = file.path(ROOT, "data", "HABS_HD"))
OUT <- file.path(Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(ROOT, "results")), "habs_multimodal")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("readxl", quietly = TRUE)) {
  stop("Install data.table and readxl.", call. = FALSE)
}
library(data.table)
library(readxl)

input_file <- function(variable, default_name) {
  Sys.getenv(variable, unset = file.path(DATA, default_name))
}
FBB_FILE <- input_file("HABS_FBB_FILE", "FBB_PET.xlsx")
TAU_FILE <- input_file("HABS_TAU_FILE", "PI2620_PET.xlsx")
MRI_FILE <- input_file("HABS_CORTICAL_THICKNESS_FILE", "cortical_thickness.xlsx")
GEN_FILE <- input_file("HABS_GENOMICS_FILE", "genomics.xlsx")
DX_FILE <- input_file("HABS_DIAGNOSIS_FILE", "diagnosis.csv")
DEM_FILE <- input_file("HABS_DEMOGRAPHICS_FILE", "demographics.csv")
ZS_FILE <- input_file("HABS_COGNITION_FILE", "cognition_zscores.csv")
MMSE_FILE <- input_file("HABS_MMSE_FILE", "mmse.csv")
required_files <- c(FBB_FILE, TAU_FILE, MRI_FILE, GEN_FILE, DX_FILE, DEM_FILE, ZS_FILE, MMSE_FILE)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing HABS-HD inputs:\n", paste(missing_files, collapse = "\n"), call. = FALSE)

write_tsv <- function(x, name) fwrite(x, file.path(OUT, name), sep = "\t", quote = FALSE, na = "NA")
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
read_habs_csv <- function(path) fread(path, na.strings = c("", "NA", "N/A", "-9999", "-8888", "-777777"))
read_inclusions <- function(path) {
  sheets <- excel_sheets(path)
  sheets <- sheets[grepl("inclusion", tolower(sheets)) & !grepl("exclusion|pending", tolower(sheets))]
  if (!length(sheets)) stop("No inclusion sheets found in ", path, call. = FALSE)
  rbindlist(lapply(sheets, function(s) {
    x <- as.data.table(suppressWarnings(read_excel(path, sheet = s)))
    x[, source_sheet := s]
    x
  }), fill = TRUE)
}
add_id_date <- function(x, year_field, month_field) {
  x[, Med_ID := as.integer(Med_ID)]
  x[, date := month_date(get(year_field), get(month_field))]
  x[!is.na(Med_ID) & !is.na(date)]
}
nearest_record <- function(x, target, max_days = 180L) {
  if (!nrow(x) || is.na(target)) return(NULL)
  gap <- abs(as.integer(x$date - target))
  keep <- which(is.finite(gap) & gap <= max_days)
  if (!length(keep)) return(NULL)
  x[keep[which.min(gap[keep])]]
}

# Prespecified HABS regions and primary pathology markers.
temporal_fields <- c(
  "L_entorhinal_thickavg", "R_entorhinal_thickavg",
  "L_parahippocampal_thickavg", "R_parahippocampal_thickavg",
  "L_inferiortemporal_thickavg", "R_inferiortemporal_thickavg",
  "L_middletemporal_thickavg", "R_middletemporal_thickavg"
)

fbb <- add_id_date(read_inclusions(FBB_FILE), "FBB_PET_year", "FBB_PET_month")
fbb[, amyloid := as_numeric_clean(Centiloid_Units_GlobalAmyloid)]
fbb <- fbb[is.finite(amyloid)]

tau <- add_id_date(read_inclusions(TAU_FILE), "PI2620_PET_year", "PI2620_PET_month")
tau[, tau_mesial := as_numeric_clean(Median_Medial_Temporal_PI2620_SUVR)]
tau <- tau[is.finite(tau_mesial)]

mri <- add_id_date(read_inclusions(MRI_FILE), "MRI_year", "MRI_month")
if (!all(temporal_fields %in% names(mri))) stop("The temporal cortical fields are incomplete.", call. = FALSE)
for (field in temporal_fields) mri[, (field) := as_numeric_clean(get(field))]
mri[, complete_temporal := complete.cases(.SD), .SDcols = temporal_fields]
mri <- mri[complete_temporal == TRUE]

hip[, hippocampus := rowMeans(cbind(as_numeric_clean(L_hippocampus), as_numeric_clean(R_hippocampus)), na.rm = TRUE)]
hip[!is.finite(hippocampus), hippocampus := NA_real_]

gen <- as.data.table(suppressWarnings(read_excel(GEN_FILE, sheet = "Sheet1")))
gen[, Med_ID := as.integer(Med_ID)]
gen <- gen[, .(Med_ID, apoe4 = as_numeric_clean(APOE4_Positivity))]

dx <- read_habs_csv(DX_FILE)
dx[, `:=`(Med_ID = as.integer(MED_ID), date = month_date(VISIT_YEAR, VISIT_MONTH), diagnosis = as_numeric_clean(CDX_COG))]
dx <- dx[!is.na(Med_ID) & !is.na(date)]

dem <- read_habs_csv(DEM_FILE)
dem[, `:=`(Med_ID = as.integer(MED_ID), date = month_date(VISIT_YEAR, VISIT_MONTH), age = as_numeric_clean(AGE), sex = factor(ID_GENDER), ethnicity = factor(ETHNICITY))]
dem <- dem[!is.na(Med_ID) & !is.na(date)]

zs <- read_habs_csv(ZS_FILE)
zs[, `:=`(Med_ID = as.integer(MED_ID), date = month_date(VISIT_YEAR, VISIT_MONTH))]
z_fields <- grep("_ZSCORE$", names(zs), value = TRUE)
for (field in z_fields) {
  zs[, (field) := as_numeric_clean(get(field))]
  zs[get(field) <= -100, (field) := NA_real_]
}
zs[, `:=`(cognitive_z = rowMeans(.SD, na.rm = TRUE), n_tests = rowSums(!is.na(.SD))), .SDcols = z_fields]
zs[n_tests < 4L, cognitive_z := NA_real_]
zs <- zs[!is.na(Med_ID) & !is.na(date)]

mmse <- read_habs_csv(MMSE_FILE)
mmse[, `:=`(Med_ID = as.integer(MED_ID), date = month_date(VISIT_YEAR, VISIT_MONTH), MMSE = as_numeric_clean(MMSE_TOTAL))]
mmse[MMSE < 0 | MMSE > 30, MMSE := NA_real_]
mmse <- mmse[!is.na(Med_ID) & !is.na(date) & is.finite(MMSE)]

# Select one baseline MRI per participant and match pathology and cognition within
# +/-180 days. Follow-up MRI is at least 365 days later, with later cognition
# collected 90 to 730 days after that MRI.
setorder(mri, Med_ID, date)
baseline_mri <- mri[, .SD[1L], by = Med_ID]
anchor_rows <- vector("list", nrow(baseline_mri))
for (i in seq_len(nrow(baseline_mri))) {
  b <- baseline_mri[i]
  sid <- b$Med_ID
  p <- nearest_record(fbb[Med_ID == sid], b$date)
  t <- nearest_record(tau[Med_ID == sid], b$date)
  c <- nearest_record(zs[Med_ID == sid], b$date)
  m <- nearest_record(mmse[Med_ID == sid], b$date)
  d <- nearest_record(dx[Med_ID == sid], b$date)
  e <- nearest_record(dem[Med_ID == sid], b$date)
  anchor_rows[[i]] <- data.table(
    Med_ID = sid,
    baseline_mri_date = b$date,
    baseline_amyloid_date = if (is.null(p)) as.IDate(NA) else p$date,
    baseline_amyloid = if (is.null(p)) NA_real_ else p$amyloid,
    baseline_tau_date = if (is.null(t)) as.IDate(NA) else t$date,
    baseline_tau = if (is.null(t)) NA_real_ else t$tau_mesial,
    baseline_cognition_date = if (is.null(c)) as.IDate(NA) else c$date,
    baseline_cognition = if (is.null(c)) NA_real_ else c$cognitive_z,
    baseline_MMSE = if (is.null(m)) NA_real_ else m$MMSE,
    diagnosis = if (is.null(d)) NA_real_ else d$diagnosis,
    age = if (is.null(e)) b$Age else e$age,
    sex = if (is.null(e)) NA_character_ else as.character(e$sex),
    ethnicity = if (is.null(e)) NA_character_ else as.character(e$ethnicity),
    scanner = as.character(b$Scanner_MRI)
  )
}
chain <- rbindlist(anchor_rows, fill = TRUE)
chain[, `:=`(
  amyloid_gap = abs(as.integer(baseline_amyloid_date - baseline_mri_date)),
  tau_gap = abs(as.integer(baseline_tau_date - baseline_mri_date)),
  cognition_gap = abs(as.integer(baseline_cognition_date - baseline_mri_date))
)]
chain <- merge(chain, gen, by = "Med_ID", all.x = TRUE)

# Cross-fitted score. Higher values denote thinner temporal cortex. The score is
# uses eight prespecified bilateral temporal regions
# against HABS outcomes.
set.seed(813L)
chain[, score_fold := sample(rep(1:5, length.out = .N))]
chain[, `:=`(
  baseline_structure = NA_real_, followup_structure = NA_real_,
  followup_mri_date = as.IDate(NA), followup_cognition = NA_real_,
  followup_MMSE = NA_real_, mri_interval_years = NA_real_
)]

for (fold_id in 1:5) {
  test_ids <- chain[score_fold == fold_id, Med_ID]
  train_ids <- chain[score_fold != fold_id, Med_ID]
  training <- mri[Med_ID %in% train_ids]
  mu <- vapply(temporal_fields, function(v) mean(training[[v]], na.rm = TRUE), numeric(1))
  sigma <- vapply(temporal_fields, function(v) sd(training[[v]], na.rm = TRUE), numeric(1))
  sigma[!is.finite(sigma) | sigma <= 0] <- 1

  baseline_rows <- mri[Med_ID %in% test_ids, .SD[1L], by = Med_ID]
  for (j in seq_along(temporal_fields)) {
    v <- temporal_fields[j]
    baseline_rows[, (v) := (as_numeric_clean(get(v)) - mu[j]) / sigma[j]]
  }
  baseline_rows[, score := -rowMeans(.SD), .SDcols = temporal_fields]
  chain[score_fold == fold_id, baseline_structure :=
          baseline_rows$score[match(Med_ID, baseline_rows$Med_ID)]]

  for (i in which(chain$score_fold == fold_id)) {
    sid <- chain$Med_ID[i]
    bdate <- chain$baseline_mri_date[i]
    followup <- mri[Med_ID == sid & date >= bdate + 365L & date <= bdate + 1825L]
    if (!nrow(followup)) next
    followup <- followup[which.max(date)]
    values <- as_numeric_clean(unlist(followup[, temporal_fields, with = FALSE], use.names = FALSE))
    values <- (values - mu) / sigma
    chain[i, `:=`(
      followup_mri_date = followup$date,
      followup_structure = -mean(values),
      mri_interval_years = as.numeric(followup$date - bdate) / 365.25
    )]
    later_cog <- zs[Med_ID == sid & date >= followup$date + 90L & date <= followup$date + 730L]
    if (nrow(later_cog)) chain[i, followup_cognition := later_cog[which.min(date), cognitive_z]]
    later_mmse <- mmse[Med_ID == sid & date >= followup$date + 90L & date <= followup$date + 730L]
    if (nrow(later_mmse)) chain[i, followup_MMSE := later_mmse[which.min(date), MMSE]]
  }
}

chain[, structure_change := (followup_structure - baseline_structure) / mri_interval_years]
chain[, `:=`(
  amyloid_z = safe_z(baseline_amyloid),
  tau_z = safe_z(baseline_tau),
  structure_z = safe_z(baseline_structure),
  cognition_z = safe_z(baseline_cognition),
  structure_change_z = safe_z(structure_change),
  later_cognition_z = safe_z(followup_cognition),
  age_z = safe_z(age),
  apoe4 = as_numeric_clean(apoe4),
  scanner = factor(scanner), sex = factor(sex), ethnicity = factor(ethnicity)
)]
write_tsv(chain, "habs_temporal_chain.tsv")

hc3_table <- function(fit, model_id, focal_terms) {
  beta <- coef(fit)
  keep <- is.finite(beta)
  X <- model.matrix(fit)[, keep, drop = FALSE]
  beta <- beta[keep]
  e <- residuals(fit) / pmax(1 - hatvalues(fit), 1e-8)
  bread <- tryCatch(solve(crossprod(X)), error = function(err) qr.solve(crossprod(X), diag(ncol(X)), tol = 1e-10))
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
fit_model <- function(data, outcome, exposures, covariates, model_id) {
  variables <- unique(c(outcome, exposures, covariates))
  complete_rows <- which(complete.cases(data[, variables, with = FALSE]))
  d <- copy(data[complete_rows])
  covariates <- covariates[vapply(covariates, function(v) length(unique(d[[v]])) > 1L, logical(1))]
  if (!nrow(d)) {
    empty <- data.table(
      model = model_id, term = character(), estimate = numeric(),
      standard_error_hc3 = numeric(), statistic = numeric(), p_value = numeric(),
      conf_low = numeric(), conf_high = numeric(), n = integer(), focal = logical(),
      status = "skipped_no_complete_cases"
    )
    return(list(data = d, fit = NULL, formula = NA_character_, table = empty,
                status = "skipped_no_complete_cases"))
  }
  formula <- reformulate(c(exposures, covariates), response = outcome)
  parameter_count <- length(attr(terms(formula), "term.labels")) + 1L
  if (nrow(d) <= parameter_count) {
    empty <- data.table(
      model = model_id, term = character(), estimate = numeric(),
      standard_error_hc3 = numeric(), statistic = numeric(), p_value = numeric(),
      conf_low = numeric(), conf_high = numeric(), n = integer(), focal = logical(),
      status = "skipped_insufficient_cases"
    )
    return(list(data = d, fit = NULL, formula = paste(deparse(formula), collapse = " "),
                table = empty, status = "skipped_insufficient_cases"))
  }
  fit <- tryCatch(
    lm(formula, data = d),
    error = function(err) err
  )
  if (inherits(fit, "error")) {
    empty <- data.table(
      model = model_id, term = character(), estimate = numeric(),
      standard_error_hc3 = numeric(), statistic = numeric(), p_value = numeric(),
      conf_low = numeric(), conf_high = numeric(), n = nrow(d), focal = logical(),
      status = paste0("skipped_lm_error: ", conditionMessage(fit))
    )
    return(list(data = d, fit = NULL, formula = paste(deparse(formula), collapse = " "),
                table = empty, status = paste0("skipped_lm_error: ", conditionMessage(fit))))
  }
  list(data = d, fit = fit, formula = paste(deparse(formula), collapse = " "),
       table = hc3_table(fit, model_id, exposures), status = "completed")
}

covariates <- c("age_z", "sex", "ethnicity", "apoe4", "scanner")
models <- list(
  list(y = "tau_z", x = "amyloid_z", cov = covariates, id = "HABS_amyloid_to_tau"),
  list(y = "structure_z", x = "tau_z", cov = covariates, id = "HABS_tau_to_temporal_structure"),
  list(y = "cognition_z", x = "structure_z", cov = c("tau_z", "amyloid_z", covariates), id = "HABS_structure_to_cognition"),
  list(y = "later_cognition_z", x = "structure_change_z", cov = c("tau_z", "amyloid_z", "cognition_z", covariates, "mri_interval_years"), id = "HABS_structure_change_to_later_cognition")
)
model_results <- list(); model_manifest <- list()
for (spec in models) {
  fitted <- fit_model(chain, spec$y, spec$x, spec$cov, spec$id)
  model_results[[length(model_results) + 1L]] <- fitted$table
  model_manifest[[length(model_manifest) + 1L]] <- data.table(
    model = spec$id, formula = fitted$formula, n = nrow(fitted$data), status = fitted$status
  )
}
effects <- rbindlist(model_results, fill = TRUE)
effects[focal == TRUE, p_fdr_within_model := p.adjust(p_value, "BH"), by = model]
write_tsv(effects, "habs_component_model_effects_hc3.tsv")
write_tsv(rbindlist(model_manifest, fill = TRUE), "habs_component_model_manifest.tsv")

# Fixed-window and scanner sensitivity analyses for the amyloid-to-tau relation.
sensitivity <- list()
for (window in c(90L, 180L, 365L)) {
  x <- chain[amyloid_gap <= window & tau_gap <= window & cognition_gap <= window]
  fitted <- fit_model(x, "tau_z", "amyloid_z", covariates, paste0("HABS_amyloid_to_tau_window_", window, "d"))
  if (nrow(fitted$table)) sensitivity[[length(sensitivity) + 1L]] <- fitted$table
  else sensitivity[[length(sensitivity) + 1L]] <- data.table(
    model = paste0("HABS_amyloid_to_tau_window_", window, "d"),
    term = NA_character_, estimate = NA_real_, standard_error_hc3 = NA_real_,
    statistic = NA_real_, p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
    n = nrow(fitted$data), focal = NA, status = fitted$status
  )
}
for (scanner_id in levels(droplevels(chain$scanner))) {
  x <- chain[scanner == scanner_id]
  if (nrow(x) < 100L) next
  fitted <- fit_model(x, "tau_z", "amyloid_z", setdiff(covariates, "scanner"), paste0("HABS_amyloid_to_tau_scanner_", scanner_id))
  if (nrow(fitted$data) >= 100L && nrow(fitted$table)) sensitivity[[length(sensitivity) + 1L]] <- fitted$table
  else sensitivity[[length(sensitivity) + 1L]] <- data.table(
    model = paste0("HABS_amyloid_to_tau_scanner_", scanner_id),
    term = NA_character_, estimate = NA_real_, standard_error_hc3 = NA_real_,
    statistic = NA_real_, p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
    n = nrow(fitted$data), focal = NA,
    status = if (nrow(fitted$data) < 100L) "skipped_fewer_than_100_complete_cases" else fitted$status
  )
}
write_tsv(rbindlist(sensitivity, fill = TRUE), "habs_time_scanner_sensitivity.tsv")

# Complete four-stage chain is explicitly exploratory.
complete_chain <- chain[is.finite(amyloid_z) & is.finite(tau_z) & is.finite(structure_change_z) & is.finite(later_cognition_z)]
exploratory_gate <- data.table(
  analysis = "FBB_to_tau_to_structure_change_to_later_cognition",
  subjects = nrow(complete_chain),
  interpretation = "Exploratory only; not a confirmatory serial mediation result"
)
write_tsv(exploratory_gate, "habs_complete_chain_exploratory_gate.tsv")

cohort_summary <- data.table(
  cohort = c(
    "Complete baseline temporal MRI",
    "MRI with FBB and tau within 180 days",
    "Repeated MRI with later cognitive composite",
    "Strict complete longitudinal chain"
  ),
  subjects = c(
    uniqueN(mri$Med_ID),
    uniqueN(chain[is.finite(amyloid_z) & is.finite(tau_z), Med_ID]),
    uniqueN(chain[is.finite(structure_change_z) & is.finite(later_cognition_z), Med_ID]),
    nrow(complete_chain)
  )
)
write_tsv(cohort_summary, "habs_cohort_summary.tsv")
capture.output({
  cat("HABS-HD multimodal cascade\n")
  cat("Prespecified regions: bilateral entorhinal, parahippocampal, inferior temporal, middle temporal\n")
  cat("Primary amyloid: Centiloid_Units_GlobalAmyloid\n")
  cat("Primary tau: Median_Medial_Temporal_PI2620_SUVR\n")
  cat("Primary cognition: mean of at least four available z-score tests\n\n")
  print(cohort_summary)
  cat("\nPrimary component effects\n")
  print(effects[focal == TRUE])
  cat("\nExploratory complete chain\n")
  print(exploratory_gate)
}, file = file.path(OUT, "analysis_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("\nHABS-HD multimodal cascade completed.\nResults: ", OUT, "\n", sep = "")
