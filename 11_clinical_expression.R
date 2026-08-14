# Clinical expression of Alzheimer pathology across A4, ADNI, and HABS-HD

options(stringsAsFactors = FALSE, scipen = 999)
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) stop("Run with source('11_clinical_expression.R').", call. = FALSE)
ROOT <- dirname(script_path)
RESULTS <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(ROOT, "results"))
OUT <- file.path(RESULTS, "clinical_expression")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("sandwich", quietly = TRUE)) {
  stop("Install data.table and sandwich.", call. = FALSE)
}
library(data.table)

A4_FILE <- Sys.getenv("A4_CLINICAL_ANALYSIS_FILE", unset = file.path(ROOT, "data", "A4", "A4_Centiloid_PACC.tsv.gz"))
ADNI_FILE <- Sys.getenv("ADNI_CLINICAL_ANALYSIS_FILE", unset = file.path(ROOT, "data", "ADNI", "ADNI_Centiloid_ADAS13.tsv.gz"))
HABS_EFFECTS <- file.path(RESULTS, "habs_primary_models", "habs_primary_model_effects_hc3.tsv")
required <- c(A4_FILE, ADNI_FILE, HABS_EFFECTS)
if (any(!file.exists(required))) stop("Clinical expression inputs are incomplete.", call. = FALSE)

safe_z <- function(x) {
  x <- as.numeric(x); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
hc3 <- function(fit, term, cohort, model) {
  v <- sandwich::vcovHC(fit, type = "HC3")
  b <- coef(fit)[term]; se <- sqrt(v[term, term]); p <- 2 * pnorm(abs(b / se), lower.tail = FALSE)
  data.table(cohort, model, term, estimate = b, std_error = se,
             ci_low = b - 1.96 * se, ci_high = b + 1.96 * se,
             p_value = p, n = nobs(fit))
}
prepare <- function(x, cohort) {
  required <- c("subject_id", "pathology", "cognition", "age", "sex", "education", "apoe4")
  if (!all(required %in% names(x))) stop(cohort, " input fields are incomplete.", call. = FALSE)
  x <- copy(x)
  x[, `:=`(pathology_z = safe_z(pathology), cognition_z = safe_z(cognition),
           age_z = safe_z(age), education_z = safe_z(education), apoe4_z = safe_z(apoe4))]
  x[, sex := factor(sex)]
  x
}

a4 <- prepare(fread(A4_FILE), "A4")
adni <- prepare(fread(ADNI_FILE), "ADNI")
if ("diagnosis" %in% names(adni)) adni[, diagnosis := factor(diagnosis)]

fit_a4 <- lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4_z, data = a4)
fit_adni <- if ("diagnosis" %in% names(adni)) {
  lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4_z + diagnosis, data = adni)
} else {
  lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4_z, data = adni)
}
effects <- rbindlist(list(
  hc3(fit_a4, "pathology_z", "A4", "Centiloid_to_PACC_impairment"),
  hc3(fit_adni, "pathology_z", "ADNI", "Centiloid_to_ADAS13")
))

stage <- list()
if ("diagnosis" %in% names(adni)) {
  for (level in levels(adni$diagnosis)) {
    d <- adni[diagnosis == level]
    if (nrow(d) >= 50L) {
      fit <- lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4_z, data = d)
      stage[[level]] <- hc3(fit, "pathology_z", "ADNI", paste0("diagnosis_", level))
    }
  }
}
stage_effects <- if (length(stage)) {
  rbindlist(stage, fill = TRUE)
} else {
  data.table(cohort = character(), model = character(), term = character(),
             estimate = numeric(), std_error = numeric(), ci_low = numeric(),
             ci_high = numeric(), p_value = numeric(), n = integer())
}
habs <- fread(HABS_EFFECTS)

fwrite(effects, file.path(OUT, "clinical_expression_primary_effects.tsv"), sep = "\t")
fwrite(stage_effects, file.path(OUT, "adni_diagnostic_stage_effects.tsv"), sep = "\t")
fwrite(habs, file.path(OUT, "habs_structural_expression_effects.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
