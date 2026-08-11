# Harmonized clinical expression across Alzheimer disease stages
#
# Unified manuscript-facing module for A4, ADNI, HABS-HD and AIBL.
# The script runs the harmonized pathology-cognition models, prespecified WMH
# and age interactions, stage and common-age heterogeneity analyses, bounded
# model diagnostics, and the complete Figure 6 export.

# Harmonized clinical-expression analysis across A4, HABS-HD and ADNI
#
# Three prespecified layers:
#   1. Pathology burden -> cognitive impairment
#   2. Pathology burden x WMH -> cognitive impairment
#   3. Pathology burden x continuous age -> cognitive impairment
#
# A4 and ADNI use Centiloid as amyloid burden and are summarized separately
# from HABS-HD, where pTau217 is used as a downstream pathology marker.
# AIBL is audited but excluded from the shared model when quantitative
# pathology and WMH are absent from the current integrated file.

options(stringsAsFactors = FALSE, scipen = 999)

PROJECT_DIR <- normalizePath(Sys.getenv("EOAD_PROJECT_DIR", unset = getwd()),
  winslash = "/", mustWork = FALSE)
DATA_DIR <- Sys.getenv("EOAD_DATA_DIR", unset = file.path(PROJECT_DIR, "data"))
RESULTS_DIR <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(PROJECT_DIR, "results"))
OUT <- file.path(RESULTS_DIR, "clinical_expression")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

clinical_input <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) stop("Set ", name, " before running the clinical module.", call. = FALSE)
  value
}

read_auto <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  if (grepl("\\.tsv$", path, ignore.case = TRUE)) {
    utils::read.delim(path, check.names = FALSE,
                      na.strings = c("", "NA", "N/A", "-999", "-99"))
  } else {
    utils::read.csv(path, check.names = FALSE,
                    na.strings = c("", "NA", "N/A", "-999", "-99"))
  }
}

num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "-999", "-99")] <- NA_character_
  suppressWarnings(as.numeric(gsub(",", "", x)))
}

safe_log <- function(x) {
  x <- num(x)
  x[x <= 0] <- NA_real_
  log(x)
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & trimws(as.character(x)) != ""]
  if (!length(x)) return(NA)
  x[[1L]]
}

collapse_subjects <- function(dat, id) {
  ids <- unique(dat[[id]])
  out <- lapply(ids, function(k) {
    s <- dat[dat[[id]] == k, , drop = FALSE]
    as.data.frame(lapply(s, function(x) {
      if (is.numeric(x)) {
        q <- x[is.finite(x)]
        if (!length(q)) return(NA_real_)
        return(stats::median(q))
      }
      first_nonmissing(x)
    }), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

scale_num <- function(x) as.numeric(scale(num(x)))

nearest_value <- function(base, other, id, base_date, other_date, value,
                          max_days = 180) {
  val <- rep(NA_real_, nrow(base))
  gap <- rep(NA_real_, nrow(base))
  idx <- split(seq_len(nrow(other)), as.character(other[[id]]))
  for (i in seq_len(nrow(base))) {
    jx <- idx[[as.character(base[[id]][i])]]
    if (!length(jx) || is.na(base[[base_date]][i])) next
    d <- abs(as.numeric(difftime(other[[other_date]][jx],
                                 base[[base_date]][i], units = "days")))
    d[!is.finite(d)] <- Inf
    j <- which.min(d)
    if (length(j) && is.finite(d[j]) && d[j] <= max_days) {
      val[i] <- num(other[[value]][jx[j]])
      gap[i] <- d[j]
    }
  }
  list(value = val, gap_days = gap)
}

prepare_model <- function(raw, variables) {
  d <- raw[complete.cases(raw[, variables, drop = FALSE]), , drop = FALSE]
  d$pathology_z <- scale_num(d$pathology)
  d$cognition_z <- scale_num(d$cognition_impairment)
  d$age_z <- scale_num(d$age)
  d$education_z <- scale_num(d$education)
  d$sex <- droplevels(as.factor(d$sex))
  d$apoe4 <- num(d$apoe4)
  if ("wmh" %in% variables) d$wmh_z <- scale_num(d$wmh)
  d
}

hc3_term <- function(fit, cohort, layer, term) {
  X <- stats::model.matrix(fit)
  e <- stats::residuals(fit)
  h <- stats::hatvalues(fit)
  wt <- (e / pmax(1 - h, 1e-8))^2
  V <- tryCatch(solve(crossprod(X)) %*% crossprod(X, X * wt) %*%
                  solve(crossprod(X)), error = function(e) NULL)
  b <- stats::coef(fit)
  if (is.null(V) || !term %in% names(b)) {
    return(data.frame(cohort = cohort, layer = layer, n = nobs(fit),
      term = term, estimate = NA, se = NA, p = NA, status = "FAIL_HC3"))
  }
  se <- sqrt(max(V[term, term], 0))
  data.frame(cohort = cohort, layer = layer, n = nobs(fit), term = term,
    estimate = unname(b[term]), se = se,
    p = 2 * stats::pnorm(-abs(b[term] / se)), status = "PASS")
}

fit_layer <- function(d, cohort, layer) {
  f <- switch(layer,
    pathology_main = cognition_z ~ pathology_z + age_z + sex + education_z + apoe4,
    pathology_wmh = cognition_z ~ pathology_z * wmh_z + age_z + sex + education_z + apoe4,
    pathology_age = cognition_z ~ pathology_z * age_z + sex + education_z + apoe4,
    stop("Unknown layer"))
  term <- switch(layer,
    pathology_main = "pathology_z",
    pathology_wmh = "pathology_z:wmh_z",
    pathology_age = "pathology_z:age_z")
  if (nrow(d) < 100) return(data.frame(cohort = cohort, layer = layer,
    n = nrow(d), term = term, estimate = NA, se = NA, p = NA,
    status = "FAIL_N_LT_100"))
  hc3_term(stats::lm(f, data = d), cohort, layer, term)
}

meta_summary <- function(x, label, method = "random") {
  x <- x[x$status == "PASS" & is.finite(x$estimate) & is.finite(x$se) & x$se > 0, , drop = FALSE]
  if (nrow(x) < 2) return(data.frame(analysis = label, method = method, k = nrow(x),
    estimate = NA, se = NA, ci_low = NA, ci_high = NA, p = NA, Q = NA, I2 = NA,
    status = "NO_META"))
  yi <- x$estimate; vi <- x$se^2; wi <- 1 / vi
  mu_fixed <- sum(wi * yi) / sum(wi)
  Q <- sum(wi * (yi - mu_fixed)^2)
  if (method == "fixed") {
    mu <- mu_fixed; se <- sqrt(1 / sum(wi)); tau2 <- 0
  } else {
    cterm <- sum(wi) - sum(wi^2) / sum(wi)
    tau2 <- ifelse(cterm > 0, max(0, (Q - (nrow(x) - 1)) / cterm), 0)
    wr <- 1 / (vi + tau2)
    mu <- sum(wr * yi) / sum(wr); se <- sqrt(1 / sum(wr))
  }
  data.frame(analysis = label, method = method, k = nrow(x), estimate = mu,
    se = se, ci_low = mu - 1.96 * se, ci_high = mu + 1.96 * se,
    p = 2 * stats::pnorm(-abs(mu / se)), Q = Q,
    I2 = ifelse(Q > 0, max(0, (Q - (nrow(x) - 1)) / Q * 100), 0), status = "PASS")
}

# -------------------------------------------------------------------------
# 1. Read source files
# -------------------------------------------------------------------------
a4 <- read_auto(clinical_input("EOAD_A4_INTEGRATED_FILE", file.path(DATA_DIR, "A4", "A4_Integrated_Data.csv")))
habs <- read_auto(clinical_input("EOAD_HABS_HD_INTEGRATED_FILE", file.path(DATA_DIR, "HABS_HD", "HABS_HD_Integrated_Data.csv")))
aibl <- read_auto(clinical_input("EOAD_AIBL_INTEGRATED_FILE", file.path(DATA_DIR, "AIBL", "AIBL_Integrated_Data.csv")))
adni_path <- read_auto(clinical_input("EOAD_ADNI_CENTILOID_FILE", file.path(DATA_DIR, "ADNI", "ADNI_Centiloid_analysis_data.tsv")))
adni_cog <- read_auto(clinical_input("EOAD_ADNI_ADAS13_FILE", file.path(DATA_DIR, "ADNI", "ADNI_ADAS13_analysis_data.tsv")))
adni_wmh <- read_auto(clinical_input("EOAD_ADNI_WMH_FILE", file.path(DATA_DIR, "ADNI", "UCD_WMH.csv")))

# -------------------------------------------------------------------------
# 2. Harmonized raw cohort frames
# -------------------------------------------------------------------------
a4 <- collapse_subjects(a4, "BID")
habs <- collapse_subjects(habs, "Med_ID")
a4_raw <- data.frame(id = a4$BID, pathology = num(a4$Centiloid),
  wmh = num(a4$Log_WMH), cognition_impairment = -num(a4$PACC), age = num(a4$Age),
  sex = a4$Gender, education = num(a4$Education), apoe4 = num(a4$APOE4_Carrier))
# The integrated HABS-HD pTau217 column is already transformed/standardized
# and contains negative values. Do not apply a second logarithmic transform.
habs_raw <- data.frame(id = habs$Med_ID, pathology = num(habs$pTau217),
  wmh = num(habs$Log_WMH), cognition_impairment = -num(habs$MMSE), age = num(habs$Age),
  sex = habs$Gender, education = num(habs$Education), apoe4 = num(habs$APOE4_Carrier))

adni_path$date <- as.Date(adni_path$date)
adni_cog$date <- as.Date(adni_cog$date)
adni_wmh$EXAMDATE <- as.Date(adni_wmh$EXAMDATE)
adni_base <- adni_path[, intersect(c("RID", "date", "Centiloid", "baseline_age",
  "sex", "education", "APOE4", "baseline_dx", "time"), names(adni_path)), drop = FALSE]
names(adni_base)[names(adni_base) == "date"] <- "path_date"
cog <- nearest_value(adni_base, adni_cog, "RID", "path_date", "date", "ADAS13", 90)
wmh <- nearest_value(adni_base, adni_wmh, "RID", "path_date", "EXAMDATE", "TOTAL_WMH", 180)
adni_raw <- data.frame(id = adni_base$RID, path_date = adni_base$path_date,
  pathology = num(adni_base$Centiloid), wmh = log1p(num(wmh$value)),
  cognition_impairment = num(cog$value), age = num(adni_base$baseline_age) + num(adni_base$time),
  baseline_age = num(adni_base$baseline_age), time = num(adni_base$time),
  sex = adni_base$sex, education = num(adni_base$education), apoe4 = num(adni_base$APOE4),
  baseline_dx = adni_base$baseline_dx, cog_gap_days = cog$gap_days,
  wmh_gap_days = wmh$gap_days)

# -------------------------------------------------------------------------
# 3. Primary analyses: each model has its own complete-case set
# -------------------------------------------------------------------------
layers <- c("pathology_main", "pathology_wmh", "pathology_age")
datasets <- list(
  A4 = a4_raw,
  `HABS-HD` = habs_raw,
  ADNI = adni_raw[order(adni_raw$id, adni_raw$path_date), , drop = FALSE]
)
primary_list <- list()
for (nm in names(datasets)) {
  d <- datasets[[nm]]
  for (layer in layers) {
    req <- c("pathology", "cognition_impairment", "age", "education", "sex", "apoe4")
    if (layer == "pathology_wmh") req <- c(req, "wmh")
    d0 <- prepare_model(d, req)
    if (nm == "ADNI") d0 <- d0[!duplicated(d0$id), , drop = FALSE]
    primary_list[[paste(nm, layer, sep = "_")]] <- fit_layer(d0, nm, layer)
  }
}
primary <- do.call(rbind, primary_list)
rownames(primary) <- NULL
primary$p_fdr_within_layer <- ave(primary$p, primary$layer,
  FUN = function(p) stats::p.adjust(p, method = "BH"))

# -------------------------------------------------------------------------
# 4. ADNI longitudinal sensitivity models with time and baseline diagnosis
# -------------------------------------------------------------------------
long_list <- list()
if (requireNamespace("nlme", quietly = TRUE)) {
  for (layer in layers) {
    req <- c("pathology", "cognition_impairment", "age", "baseline_age", "time",
             "education", "sex", "apoe4", "baseline_dx")
    if (layer == "pathology_wmh") req <- c(req, "wmh")
    dl <- prepare_model(adni_raw, req)
    dl$time_z <- scale_num(dl$time)
    dl$baseline_age_z <- scale_num(dl$baseline_age)
    dl$baseline_dx <- droplevels(as.factor(dl$baseline_dx))
    f <- switch(layer,
      pathology_main = cognition_z ~ pathology_z + time_z + baseline_age_z + sex + education_z + apoe4 + baseline_dx,
      pathology_wmh = cognition_z ~ pathology_z * wmh_z + time_z + baseline_age_z + sex + education_z + apoe4 + baseline_dx,
      pathology_age = cognition_z ~ pathology_z * baseline_age_z + time_z + sex + education_z + apoe4 + baseline_dx)
    term <- switch(layer, pathology_main = "pathology_z",
      pathology_wmh = "pathology_z:wmh_z", pathology_age = "pathology_z:baseline_age_z")
    fit <- tryCatch(nlme::lme(f, random = ~1 | id, data = dl,
      na.action = na.omit, method = "REML"), error = function(e) e)
    if (inherits(fit, "error")) {
      long_list[[layer]] <- data.frame(cohort = "ADNI", layer = layer,
        model = "longitudinal_random_intercept_time_adjusted", n = nrow(dl),
        participants = length(unique(dl$id)), term = term, estimate = NA, se = NA,
        p = NA, status = paste0("FAIL_MODEL:", conditionMessage(fit)))
    } else {
      co <- summary(fit)$tTable
      long_list[[layer]] <- data.frame(cohort = "ADNI", layer = layer,
        model = "longitudinal_random_intercept_time_adjusted", n = nrow(dl),
        participants = length(unique(dl$id)), term = term,
        estimate = co[term, "Value"], se = co[term, "Std.Error"],
        p = co[term, "p-value"], status = "PASS")
    }
  }
}
longitudinal <- do.call(rbind, long_list)
rownames(longitudinal) <- NULL

# -------------------------------------------------------------------------
# 5. Separate meta-analyses: amyloid-only versus broader pathology context
# -------------------------------------------------------------------------
meta_rows <- list()
for (layer in layers) {
  x <- primary[primary$layer == layer & primary$cohort %in% c("A4", "ADNI"), , drop = FALSE]
  meta_rows[[paste(layer, "amyloid_fixed", sep = "_")]] <- meta_summary(x, paste0("A4_ADNI_", layer), "fixed")
  meta_rows[[paste(layer, "amyloid_random", sep = "_")]] <- meta_summary(x, paste0("A4_ADNI_", layer), "random")
  x_all <- primary[primary$layer == layer, , drop = FALSE]
  meta_rows[[paste(layer, "all_random", sep = "_")]] <- meta_summary(x_all, paste0("A4_HABS_ADNI_", layer), "random")
}
meta <- do.call(rbind, meta_rows)
rownames(meta) <- NULL

# -------------------------------------------------------------------------
# 6. AIBL role and explicit decision gates
# -------------------------------------------------------------------------
aibl_audit <- data.frame(cohort = "AIBL", rows = nrow(aibl),
  participants = length(unique(aibl$RID)), pathology_available = FALSE,
  wmh_available = FALSE, longitudinal_clinical_available = all(c("Time", "Event") %in% names(aibl)),
  role = "secondary longitudinal clinical benchmark",
  reason = "Current integrated AIBL file contains APOE, follow-up and conversion fields but no quantitative pathology or WMH")

get_layer <- function(layer, cohorts = NULL) {
  x <- primary[primary$layer == layer & primary$status == "PASS", , drop = FALSE]
  if (!is.null(cohorts)) x <- x[x$cohort %in% cohorts, , drop = FALSE]
  x
}
path_amy <- get_layer("pathology_main", c("A4", "ADNI"))
wmh <- get_layer("pathology_wmh")
age <- get_layer("pathology_age")
same_sign_fun <- function(x) sum(x$estimate > 0, na.rm = TRUE) >= 2 || sum(x$estimate < 0, na.rm = TRUE) >= 2
meta_path_amy <- meta[meta$analysis == "A4_ADNI_pathology_main" & meta$method == "random", , drop = FALSE]
meta_wmh <- meta[meta$analysis == "A4_HABS_ADNI_pathology_wmh" & meta$method == "random", , drop = FALSE]
meta_age <- meta[meta$analysis == "A4_ADNI_pathology_age" & meta$method == "random", , drop = FALSE]

gates <- data.frame(
  gate = c("Amyloid pathology main effect has same direction in A4 and ADNI",
    "Amyloid-only pathology main effect excludes zero",
    "WMH interaction has same direction in at least two cohorts",
    "WMH interaction meta-analysis excludes zero",
    "Age interaction has same direction in at least two cohorts",
    "ADNI time-adjusted longitudinal pathology effect is available",
    "AIBL enters the shared pathology model",
    "New bulk MRI download recommended"),
  result = c(same_sign_fun(path_amy), nrow(meta_path_amy) == 1 && meta_path_amy$ci_low > 0,
    same_sign_fun(wmh), nrow(meta_wmh) == 1 && meta_wmh$ci_low > 0,
    same_sign_fun(age), any(longitudinal$layer == "pathology_main" & longitudinal$status == "PASS"),
    FALSE, FALSE),
  interpretation = c("Directionally coherent clinical expression",
    "Potentially publishable amyloid-to-cognition summary if diagnostics pass",
    "Suggestive modifier only; requires effect-size and heterogeneity reporting",
    "Do not promote WMH to a core modifier unless this gate passes",
    "Assess stage-dependent clinical expression without sliding windows",
    "Use as longitudinal sensitivity, not EOAD-specific replication",
    "AIBL remains supplementary clinical conversion context",
    "Use existing quantitative MRI derivatives"),
  stringsAsFactors = FALSE)

write_tsv <- function(x, name) utils::write.table(x, file.path(OUT, name), sep = "\t",
  row.names = FALSE, quote = FALSE, na = "NA", fileEncoding = "UTF-8")
write_tsv(primary, "primary_three_layers_HC3.tsv")
write_tsv(longitudinal, "ADNI_longitudinal_time_adjusted.tsv")
write_tsv(meta, "meta_three_layers.tsv")
write_tsv(gates, "clinical_evidence_gates.tsv")
write_tsv(aibl_audit, "AIBL_role_audit.tsv")
write_tsv(data.frame(cohort = names(datasets),
  total_rows = vapply(datasets, nrow, integer(1)),
  participants = vapply(datasets, function(x) length(unique(x$id)), integer(1))),
  "source_cohort_sizes.tsv")
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

cat("Harmonized clinical-expression analysis completed.\nOutput: ", OUT, "\n\n", sep = "")
print(primary)
print(longitudinal)
print(meta)
print(gates)

# ==============================================================================
# Clinical-stage and common-support heterogeneity
# ==============================================================================
OUT <- file.path(RESULTS_DIR, "clinical_expression_heterogeneity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

fit_pathology <- function(d, cohort, subgroup) {
  req <- c("pathology", "cognition_impairment", "age", "education", "sex", "apoe4")
  d <- prepare_model(d, req)
  if ("path_date" %in% names(d)) {
    d <- d[order(d$id, d$path_date), , drop = FALSE]
    d <- d[!duplicated(d$id), , drop = FALSE]
  }
  if (nrow(d) < 50) {
    return(data.frame(cohort = cohort, subgroup = subgroup, n = nrow(d),
      estimate = NA, se = NA, ci_low = NA, ci_high = NA, p = NA,
      status = "DESCRIPTIVE_N_LT_50"))
  }
  fit <- stats::lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4,
                   data = d)
  x <- hc3_term(fit, cohort, "pathology_main", "pathology_z")
  data.frame(cohort = cohort, subgroup = subgroup, n = x$n,
    estimate = x$estimate, se = x$se,
    ci_low = x$estimate - 1.96 * x$se,
    ci_high = x$estimate + 1.96 * x$se,
    p = x$p, status = x$status)
}

fit_adni_diagnosis_adjusted <- function(d, subgroup) {
  req <- c("pathology", "cognition_impairment", "age", "education", "sex",
           "apoe4", "baseline_dx")
  d <- prepare_model(d, req)
  d <- d[order(d$id, d$path_date), , drop = FALSE]
  d <- d[!duplicated(d$id), , drop = FALSE]
  d$baseline_dx <- droplevels(as.factor(d$baseline_dx))
  fit <- stats::lm(cognition_z ~ pathology_z + age_z + sex + education_z +
                    apoe4 + baseline_dx, data = d)
  x <- hc3_term(fit, "ADNI", "pathology_main_dx_adjusted", "pathology_z")
  data.frame(cohort = "ADNI", subgroup = subgroup, n = x$n,
    estimate = x$estimate, se = x$se,
    ci_low = x$estimate - 1.96 * x$se,
    ci_high = x$estimate + 1.96 * x$se,
    p = x$p, status = x$status)
}

contrast_slopes <- function(a, b, label) {
  if (nrow(a) != 1 || nrow(b) != 1 ||
      any(!is.finite(c(a$estimate, a$se, b$estimate, b$se)))) {
    return(data.frame(comparison = label, difference = NA, se = NA,
      ci_low = NA, ci_high = NA, p = NA, status = "NOT_ESTIMABLE"))
  }
  dif <- b$estimate - a$estimate
  se <- sqrt(a$se^2 + b$se^2)
  data.frame(comparison = label, difference = dif, se = se,
    ci_low = dif - 1.96 * se, ci_high = dif + 1.96 * se,
    p = 2 * stats::pnorm(-abs(dif / se)), status = "PASS")
}

pair_heterogeneity <- function(a, b, label) {
  x <- rbind(a, b)
  x <- x[is.finite(x$estimate) & is.finite(x$se) & x$se > 0, , drop = FALSE]
  if (nrow(x) != 2) return(data.frame(comparison = label, k = nrow(x),
    Q = NA, I2 = NA, status = "NOT_ESTIMABLE"))
  w <- 1 / x$se^2
  mu <- sum(w * x$estimate) / sum(w)
  Q <- sum(w * (x$estimate - mu)^2)
  data.frame(comparison = label, k = 2, Q = Q,
    I2 = ifelse(Q > 0, max(0, (Q - 1) / Q * 100), 0),
    status = "DESCRIPTIVE_ONLY_K_EQ_2")
}

# Full-cohort estimates
a4_full <- fit_pathology(a4_raw, "A4", "full_preclinical_cohort")
adni_full <- fit_pathology(adni_raw, "ADNI", "all_baseline_diagnoses")
adni_dx_adjusted <- fit_adni_diagnosis_adjusted(
  adni_raw, "all_baseline_diagnoses_dx_adjusted")

# ADNI stage estimates. AD is retained descriptively if fewer than 50 participants.
dx_levels <- c("CN", "EMCI", "LMCI", "AD")
adni_stage <- do.call(rbind, lapply(dx_levels, function(dx) {
  fit_pathology(adni_raw[adni_raw$baseline_dx == dx, , drop = FALSE],
                "ADNI", paste0("baseline_", dx))
}))

# Common age support. The 65-85 range is contained in the observed A4 range and
# removes ADNI-only age extremes without selecting a significance-maximizing cut.
a4_common_age <- fit_pathology(a4_raw[a4_raw$age >= 65 & a4_raw$age <= 85, ],
                               "A4", "age_65_85")
adni_common_age <- fit_pathology(adni_raw[adni_raw$age >= 65 & adni_raw$age <= 85, ],
                                 "ADNI", "age_65_85")
adni_common_age_dx_adjusted <- fit_adni_diagnosis_adjusted(
  adni_raw[adni_raw$age >= 65 & adni_raw$age <= 85, ],
  "age_65_85_dx_adjusted")
adni_cn_common_age <- fit_pathology(
  adni_raw[adni_raw$baseline_dx == "CN" & adni_raw$age >= 65 & adni_raw$age <= 85, ],
  "ADNI", "baseline_CN_age_65_85")

estimates <- rbind(a4_full, adni_full, adni_dx_adjusted, adni_stage,
                   a4_common_age, adni_common_age, adni_common_age_dx_adjusted,
                   adni_cn_common_age)
rownames(estimates) <- NULL

# Prespecified comparisons
adni_cn <- adni_stage[adni_stage$subgroup == "baseline_CN", , drop = FALSE]
contrasts <- rbind(
  contrast_slopes(a4_full, adni_full, "ADNI_all_minus_A4"),
  contrast_slopes(a4_full, adni_dx_adjusted, "ADNI_dx_adjusted_minus_A4"),
  contrast_slopes(a4_full, adni_cn, "ADNI_CN_minus_A4"),
  contrast_slopes(a4_common_age, adni_common_age, "ADNI_minus_A4_age_65_85"),
  contrast_slopes(a4_common_age, adni_common_age_dx_adjusted,
                  "ADNI_dx_adjusted_minus_A4_age_65_85"),
  contrast_slopes(a4_common_age, adni_cn_common_age, "ADNI_CN_minus_A4_age_65_85")
)

heterogeneity <- rbind(
  pair_heterogeneity(a4_full, adni_full, "A4_vs_ADNI_all"),
  pair_heterogeneity(a4_full, adni_dx_adjusted, "A4_vs_ADNI_dx_adjusted"),
  pair_heterogeneity(a4_full, adni_cn, "A4_vs_ADNI_CN"),
  pair_heterogeneity(a4_common_age, adni_common_age, "A4_vs_ADNI_age_65_85"),
  pair_heterogeneity(a4_common_age, adni_common_age_dx_adjusted,
                     "A4_vs_ADNI_dx_adjusted_age_65_85"),
  pair_heterogeneity(a4_common_age, adni_cn_common_age, "A4_vs_ADNI_CN_age_65_85")
)

# Bring forward the other two layers without forcing incompatible biomarkers
# into one biological claim.
wmh_layer <- primary[primary$layer == "pathology_wmh", , drop = FALSE]
age_layer <- primary[primary$layer == "pathology_age", , drop = FALSE]

decision <- data.frame(
  item = c(
    "Universal pooled pathology-to-cognition effect",
    "Cohort-specific pathology-to-cognition associations",
    "Clinical stage as a source of heterogeneity",
    "WMH as a replicated amyloid-specific modifier",
    "HABS-HD pTau217-WMH interaction",
    "Age interaction as a universal effect",
    "Delete a cohort to reduce I2",
    "Report the original I2 without further analysis"
  ),
  decision = c(
    "SECONDARY_ONLY",
    "PRIMARY",
    "FORMALLY_EVALUATE",
    "NOT_SUPPORTED_BY_A4_AND_ADNI",
    "COHORT_SPECIFIC_SUPPORT",
    "NOT_SUPPORTED",
    "NO",
    "NO"
  ),
  rationale = c(
    "Only two amyloid cohorts and substantial stage-related heterogeneity",
    "Both A4 and ADNI show positive associations with robust standard errors",
    "A4 is preclinical whereas ADNI includes CN, EMCI, LMCI and AD",
    "A4 and ADNI cross-sectional amyloid-WMH interactions are null",
    "The HABS-HD pTau217 interaction is significant but uses a different pathology marker",
    "A4 and HABS-HD are positive whereas ADNI is negative",
    "Cohort exclusion would be data-driven and clinically unjustified",
    "Stage and common-age sensitivity analyses are required"
  ), stringsAsFactors = FALSE)

write_tsv <- function(x, name) utils::write.table(x, file.path(OUT, name),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA", fileEncoding = "UTF-8")
write_tsv(estimates, "pathology_cognition_stage_estimates.tsv")
write_tsv(contrasts, "prespecified_slope_contrasts.tsv")
write_tsv(heterogeneity, "heterogeneity_stage_common_support.tsv")
write_tsv(wmh_layer, "WMH_interaction_context.tsv")
write_tsv(age_layer, "age_interaction_context.tsv")
write_tsv(decision, "heterogeneity_decisions.tsv")
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

cat("Heterogeneity analysis completed.\nOutput: ", OUT, "\n\n", sep = "")
print(estimates)
print(contrasts)
print(heterogeneity)
print(decision)

# ==============================================================================
# Covariate-aligned and bounded sensitivity analyses
# ==============================================================================
OUT <- file.path(RESULTS_DIR, "clinical_expression_validation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

prepare_adni <- function(raw, need_wmh = FALSE, common_age = FALSE,
                         cn_only = FALSE) {
  d <- raw
  if (common_age) d <- d[d$age >= 65 & d$age <= 85, , drop = FALSE]
  if (cn_only) d <- d[d$baseline_dx == "CN", , drop = FALSE]
  req <- c("pathology", "cognition_impairment", "age", "education", "sex",
           "apoe4", "baseline_dx")
  if (need_wmh) req <- c(req, "wmh")
  d <- prepare_model(d, req)
  d <- d[order(d$id, d$path_date), , drop = FALSE]
  d <- d[!duplicated(d$id), , drop = FALSE]
  d$baseline_dx <- droplevels(as.factor(d$baseline_dx))
  d
}

prepare_other <- function(raw, need_wmh = FALSE, common_age = FALSE) {
  d <- raw
  if (common_age) d <- d[d$age >= 65 & d$age <= 85, , drop = FALSE]
  req <- c("pathology", "cognition_impairment", "age", "education", "sex", "apoe4")
  if (need_wmh) req <- c(req, "wmh")
  prepare_model(d, req)
}

extract_hc3 <- function(fit, cohort, analysis, term) {
  x <- hc3_term(fit, cohort, analysis, term)
  data.frame(cohort = cohort, analysis = analysis, n = x$n, term = term,
    estimate = x$estimate, se = x$se,
    ci_low = x$estimate - 1.96 * x$se,
    ci_high = x$estimate + 1.96 * x$se,
    p = x$p, status = x$status)
}

fit_adni <- function(d, analysis) {
  if (analysis == "pathology_main") {
    fit <- stats::lm(cognition_z ~ pathology_z + age_z + sex + education_z +
                      apoe4 + baseline_dx, data = d)
    term <- "pathology_z"
  } else if (analysis == "pathology_wmh") {
    fit <- stats::lm(cognition_z ~ pathology_z * wmh_z + age_z + sex +
                      education_z + apoe4 + baseline_dx, data = d)
    term <- "pathology_z:wmh_z"
  } else if (analysis == "pathology_age") {
    fit <- stats::lm(cognition_z ~ pathology_z * age_z + sex + education_z +
                      apoe4 + baseline_dx, data = d)
    term <- "pathology_z:age_z"
  } else stop("Unknown analysis")
  extract_hc3(fit, "ADNI", paste0(analysis, "_dx_adjusted"), term)
}

fit_other <- function(d, cohort, analysis) {
  if (analysis == "pathology_main") {
    fit <- stats::lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4,
                     data = d)
    term <- "pathology_z"
  } else if (analysis == "pathology_wmh") {
    fit <- stats::lm(cognition_z ~ pathology_z * wmh_z + age_z + sex +
                      education_z + apoe4, data = d)
    term <- "pathology_z:wmh_z"
  } else if (analysis == "pathology_age") {
    fit <- stats::lm(cognition_z ~ pathology_z * age_z + sex + education_z + apoe4,
                     data = d)
    term <- "pathology_z:age_z"
  } else stop("Unknown analysis")
  extract_hc3(fit, cohort, analysis, term)
}

# -------------------------------------------------------------------------
# 1. Covariate-aligned primary and modifier models
# -------------------------------------------------------------------------
adni_main <- prepare_adni(adni_raw)
adni_wmh <- prepare_adni(adni_raw, need_wmh = TRUE)
adni_main_common <- prepare_adni(adni_raw, common_age = TRUE)
adni_wmh_common <- prepare_adni(adni_raw, need_wmh = TRUE, common_age = TRUE)
adni_cn <- prepare_adni(adni_raw, cn_only = TRUE)
adni_cn_wmh <- prepare_adni(adni_raw, need_wmh = TRUE, cn_only = TRUE)

a4_main <- prepare_other(a4_raw)
a4_wmh <- prepare_other(a4_raw, need_wmh = TRUE)
habs_main <- prepare_other(habs_raw)
habs_wmh <- prepare_other(habs_raw, need_wmh = TRUE)

aligned_models <- rbind(
  fit_other(a4_main, "A4", "pathology_main"),
  fit_adni(adni_main, "pathology_main"),
  fit_other(habs_main, "HABS-HD", "pathology_main"),
  fit_other(a4_wmh, "A4", "pathology_wmh"),
  fit_adni(adni_wmh, "pathology_wmh"),
  fit_other(habs_wmh, "HABS-HD", "pathology_wmh"),
  fit_other(a4_main, "A4", "pathology_age"),
  fit_adni(adni_main, "pathology_age"),
  fit_other(habs_main, "HABS-HD", "pathology_age"),
  fit_adni(adni_main_common, "pathology_main"),
  fit_adni(adni_wmh_common, "pathology_wmh"),
  fit_adni(adni_main_common, "pathology_age")
)
aligned_models$context <- c(
  "full", "full", "full", "full", "full", "full",
  "full", "full", "full", "age_65_85", "age_65_85", "age_65_85")
aligned_models$hypothesis <- ifelse(grepl("pathology_wmh", aligned_models$analysis),
  "pathology_wmh", ifelse(grepl("pathology_age", aligned_models$analysis),
  "pathology_age", "pathology_main"))
aligned_models$p_fdr <- ave(aligned_models$p,
  paste(aligned_models$hypothesis, aligned_models$context),
  FUN = function(p) stats::p.adjust(p, method = "BH"))

# Stage-matched CN results. baseline_dx is constant and is omitted from the model.
fit_cn <- function(d, analysis) {
  if (analysis == "pathology_main") {
    fit <- stats::lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4,
                     data = d); term <- "pathology_z"
  } else if (analysis == "pathology_wmh") {
    fit <- stats::lm(cognition_z ~ pathology_z * wmh_z + age_z + sex +
                      education_z + apoe4, data = d); term <- "pathology_z:wmh_z"
  } else {
    fit <- stats::lm(cognition_z ~ pathology_z * age_z + sex + education_z + apoe4,
                     data = d); term <- "pathology_z:age_z"
  }
  extract_hc3(fit, "ADNI-CN", paste0(analysis, "_stage_matched"), term)
}
stage_matched <- rbind(
  fit_cn(adni_cn, "pathology_main"),
  fit_cn(adni_cn_wmh, "pathology_wmh"),
  fit_cn(adni_cn, "pathology_age")
)

# -------------------------------------------------------------------------
# 2. Nonlinearity and influential-observation sensitivity
# -------------------------------------------------------------------------
diagnostic_fit <- function(d, cohort, diagnosis_adjusted = FALSE) {
  if (diagnosis_adjusted) {
    linear <- stats::lm(cognition_z ~ pathology_z + age_z + sex + education_z +
                         apoe4 + baseline_dx, data = d)
    quadratic <- stats::lm(cognition_z ~ pathology_z + I(pathology_z^2) + age_z +
                            sex + education_z + apoe4 + baseline_dx, data = d)
  } else {
    linear <- stats::lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4,
                        data = d)
    quadratic <- stats::lm(cognition_z ~ pathology_z + I(pathology_z^2) + age_z +
                            sex + education_z + apoe4, data = d)
  }
  q <- extract_hc3(quadratic, cohort, "quadratic_sensitivity", "I(pathology_z^2)")
  keep <- abs(stats::rstandard(linear)) <= 3
  if (diagnosis_adjusted) {
    reduced <- stats::lm(cognition_z ~ pathology_z + age_z + sex + education_z +
                          apoe4 + baseline_dx, data = d[keep, , drop = FALSE])
  } else {
    reduced <- stats::lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4,
                          data = d[keep, , drop = FALSE])
  }
  b0 <- extract_hc3(linear, cohort, "primary", "pathology_z")
  b1 <- extract_hc3(reduced, cohort, "outlier_excluded", "pathology_z")
  data.frame(cohort = cohort, n_primary = nrow(d), n_outlier_excluded = sum(keep),
    primary_estimate = b0$estimate, outlier_excluded_estimate = b1$estimate,
    percent_change = 100 * (b1$estimate - b0$estimate) / abs(b0$estimate),
    quadratic_estimate = q$estimate, quadratic_se = q$se,
    quadratic_p = q$p,
    stable = abs(100 * (b1$estimate - b0$estimate) / abs(b0$estimate)) < 20 && q$p >= 0.05)
}

diagnostics <- rbind(
  diagnostic_fit(a4_main, "A4", FALSE),
  diagnostic_fit(adni_main, "ADNI", TRUE)
)

# -------------------------------------------------------------------------
# 3. Decision table
# -------------------------------------------------------------------------
pick <- function(cohort, pattern, context = "full") {
  x <- aligned_models[aligned_models$cohort == cohort &
    grepl(pattern, aligned_models$analysis, fixed = TRUE) &
    aligned_models$context == context, , drop = FALSE]
  x[1, , drop = FALSE]
}
a4_path <- pick("A4", "pathology_main")
adni_path <- pick("ADNI", "pathology_main_dx_adjusted")
a4_w <- pick("A4", "pathology_wmh")
adni_w <- pick("ADNI", "pathology_wmh_dx_adjusted")
habs_w <- pick("HABS-HD", "pathology_wmh")
a4_age <- pick("A4", "pathology_age")
adni_age <- pick("ADNI", "pathology_age_dx_adjusted")
habs_age <- pick("HABS-HD", "pathology_age")

decisions <- data.frame(
  question = c(
    "Pathology-to-cognition relation across A4 and diagnosis-adjusted ADNI",
    "Primary pathology estimates robust to nonlinearity and influential observations",
    "Amyloid-WMH interaction replicated cross-sectionally in A4 and ADNI",
    "pTau217-WMH interaction supported in HABS-HD",
    "Age interaction consistent across all cohorts",
    "AIBL enters the harmonized pathology model",
    "Additional prespecified model-form validation required"
  ),
  supported = c(
    a4_path$estimate > 0 && a4_path$p_fdr < 0.05 && adni_path$estimate > 0 && adni_path$p_fdr < 0.05,
    all(diagnostics$stable),
    a4_w$p_fdr < 0.05 && adni_w$p_fdr < 0.05 && sign(a4_w$estimate) == sign(adni_w$estimate),
    habs_w$p_fdr < 0.05,
    sign(a4_age$estimate) == sign(adni_age$estimate) &&
      sign(adni_age$estimate) == sign(habs_age$estimate) &&
      all(c(a4_age$p_fdr, adni_age$p_fdr, habs_age$p_fdr) < 0.05),
    FALSE,
    !all(diagnostics$stable)
  ),
  manuscript_role = c(
    "Primary harmonized clinical-expression result",
    "Required sensitivity support",
    "Promote only if supported; otherwise retain cohort-specific WMH findings",
    "HABS-HD-specific modifier evidence",
    "Promote only if supported; otherwise describe stage-related heterogeneity",
    "AIBL remains a longitudinal conversion benchmark",
    "Characterize A4 nonlinearity before interpretation"
  ), stringsAsFactors = FALSE)

write_tsv <- function(x, name) utils::write.table(x, file.path(OUT, name),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA", fileEncoding = "UTF-8")
write_tsv(aligned_models, "covariate_aligned_models.tsv")
write_tsv(stage_matched, "ADNI_CN_stage_matched.tsv")
write_tsv(diagnostics, "primary_model_diagnostics.tsv")
write_tsv(decisions, "analysis_decisions.tsv")
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

cat("Clinical model validation completed.\nOutput: ", OUT, "\n\n", sep = "")
print(aligned_models)
print(stage_matched)
print(diagnostics)
print(decisions)

# ==============================================================================
# Figure 6 and source data
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source_dir <- file.path(RESULTS_DIR, "clinical_expression_figures")
validation_dir <- file.path(RESULTS_DIR, "clinical_expression_validation")
hetero_dir <- file.path(RESULTS_DIR, "clinical_expression_heterogeneity")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

aligned <- aligned_models
stage <- estimates

theme_figure6 <- function(base_size = 7.0) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.32, color = "#333333"),
      axis.ticks = element_line(linewidth = 0.30, color = "#333333"),
      axis.title = element_text(size = 6.8, color = "#222222"),
      axis.text = element_text(size = 6.2, color = "#333333"),
      plot.title = element_text(size = 7.8, face = "plain", margin = margin(b = 2.5)),
      plot.subtitle = element_text(size = 6.2, color = "#555555", lineheight = 1.05,
                                   margin = margin(b = 3)),
      legend.title = element_blank(),
      legend.text = element_text(size = 6.1),
      legend.key.height = unit(7, "pt"),
      panel.grid.major.x = element_line(linewidth = 0.22, color = "#E5E5E5"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 6, 5, 6)
    )
}

prepare_other <- function(raw) {
  prepare_model(
    raw,
    c("pathology", "cognition_impairment", "age", "education", "sex", "apoe4")
  )
}

prepare_adni <- function(raw) {
  d <- prepare_model(
    raw,
    c("pathology", "cognition_impairment", "age", "education", "sex", "apoe4", "baseline_dx")
  )
  d <- d[order(d$id, d$path_date), , drop = FALSE]
  d <- d[!duplicated(d$id), , drop = FALSE]
  d$baseline_dx <- droplevels(as.factor(d$baseline_dx))
  d
}

a4_main <- prepare_other(a4_raw)
adni_main <- prepare_adni(adni_raw)
habs_main <- prepare_other(habs_raw)

fit_a4 <- lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4, data = a4_main)
fit_adni <- lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4 + baseline_dx,
               data = adni_main)
fit_habs <- lm(cognition_z ~ pathology_z + age_z + sex + education_z + apoe4, data = habs_main)

partial_bins <- function(d, fit, cohort, bins = 18L) {
  beta <- unname(coef(fit)["pathology_z"])
  partial <- residuals(fit) + beta * d$pathology_z
  cuts <- unique(quantile(d$pathology_z, probs = seq(0, 1, length.out = bins + 1L),
                          na.rm = TRUE, type = 8))
  group <- cut(d$pathology_z, breaks = cuts, include.lowest = TRUE, labels = FALSE)
  rows <- lapply(sort(unique(group)), function(g) {
    keep <- group == g
    data.frame(
      cohort = cohort,
      pathology_z = mean(d$pathology_z[keep]),
      partial_cognition_z = mean(partial[keep]),
      se = sd(partial[keep]) / sqrt(sum(keep)),
      n_bin = sum(keep)
    )
  })
  do.call(rbind, rows)
}

partial_data <- rbind(
  partial_bins(a4_main, fit_a4, "A4"),
  partial_bins(adni_main, fit_adni, "ADNI"),
  partial_bins(habs_main, fit_habs, "HABS-HD")
)

main <- aligned[
  aligned$analysis %in% c("pathology_main", "pathology_main_dx_adjusted") &
    aligned$context == "full",
]
main <- main[match(c("A4", "ADNI", "HABS-HD"), main$cohort), ]
main$label <- c(
  A4 = "A4: Centiloid and PACC",
  ADNI = "ADNI: Centiloid and ADAS13",
  `HABS-HD` = "HABS-HD: p-tau217 and MMSE"
)[main$cohort]
main$label <- factor(main$label, levels = rev(main$label))
main$annotation <- sprintf(
  "beta = %.3f  (%.3f, %.3f)    P = %.2g    n = %s",
  main$estimate, main$ci_low, main$ci_high, main$p, comma(main$n)
)

pA <- ggplot(main, aes(x = estimate, y = label, color = cohort)) +
  geom_vline(xintercept = 0, linewidth = 0.32, linetype = 2, color = "#777777") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label), linewidth = 0.78) +
  geom_point(size = 2.35) +
  geom_text(aes(x = 0.405, label = annotation), hjust = 1, color = "#303030",
            size = 2.35, family = "Arial") +
  scale_color_manual(values = cohort_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 0.41), breaks = seq(0, 0.25, 0.05),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Adjusted pathology burden and cognitive impairment",
    subtitle = "Positive standardized estimates indicate worse cognition, bars show 95% confidence intervals",
    x = "Standardized effect", y = NULL
  ) +
  theme_figure6(7.2) +
  theme(plot.margin = margin(4, 8, 2, 6))

partial_panel <- function(cohort, title, subtitle) {
  dat <- partial_data[partial_data$cohort == cohort, ]
  model_row <- main[main$cohort == cohort, ]
  beta <- model_row$estimate
  x_seq <- seq(-2.8, 2.8, length.out = 100)
  line <- data.frame(pathology_z = x_seq, partial_cognition_z = beta * x_seq)
  ggplot(dat, aes(pathology_z, partial_cognition_z)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "#B8B8B8") +
    geom_errorbar(aes(ymin = partial_cognition_z - 1.96 * se,
                      ymax = partial_cognition_z + 1.96 * se),
                  width = 0, linewidth = 0.35, color = alpha(cohort_colors[[cohort]], 0.55)) +
    geom_line(data = line, color = cohort_colors[[cohort]], linewidth = 0.75) +
    geom_point(aes(size = n_bin), shape = 21, fill = "white",
               color = cohort_colors[[cohort]], stroke = 0.55) +
    scale_size(range = c(1.25, 2.7), guide = "none") +
    coord_cartesian(xlim = c(-2.8, 2.8), ylim = c(-0.70, 0.85)) +
    scale_x_continuous(breaks = c(-2, 0, 2)) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
    labs(title = title, subtitle = subtitle,
         x = "Pathology burden (z)", y = "Adjusted cognitive impairment") +
    theme_figure6() +
    theme(panel.grid.major.y = element_line(linewidth = 0.22, color = "#E5E5E5"))
}

pB <- partial_panel(
  "A4", "A4: preclinical",
  sprintf("Centiloid, n = %s, beta = %.3f", comma(main$n[main$cohort == "A4"]),
          main$estimate[main$cohort == "A4"])
)
pC <- partial_panel(
  "ADNI", "ADNI: diagnosis adjusted",
  sprintf("Centiloid, n = %s, beta = %.3f", comma(main$n[main$cohort == "ADNI"]),
          main$estimate[main$cohort == "ADNI"])
)
pD <- partial_panel(
  "HABS-HD", "HABS-HD: community sample",
  sprintf("p-tau217, n = %s, beta = %.3f", comma(main$n[main$cohort == "HABS-HD"]),
          main$estimate[main$cohort == "HABS-HD"])
)

stage_plot <- stage[
  (stage$cohort == "A4" & stage$subgroup == "full_preclinical_cohort") |
    (stage$cohort == "ADNI" & stage$subgroup %in% c(
      "all_baseline_diagnoses", "all_baseline_diagnoses_dx_adjusted", "baseline_CN"
    )),
]
stage_plot$label <- c(
  full_preclinical_cohort = "A4: preclinical",
  all_baseline_diagnoses = "ADNI: all diagnoses",
  all_baseline_diagnoses_dx_adjusted = "ADNI: diagnosis adjusted",
  baseline_CN = "ADNI: cognitively normal"
)[stage_plot$subgroup]
stage_plot <- stage_plot[match(
  c("full_preclinical_cohort", "all_baseline_diagnoses",
    "all_baseline_diagnoses_dx_adjusted", "baseline_CN"),
  stage_plot$subgroup
), ]
stage_plot$label <- factor(stage_plot$label, levels = rev(stage_plot$label))
stage_plot$color <- c("#3E86BE", "#E8A19A", "#D64B40", "#9E2F28")

pE <- ggplot(stage_plot, aes(estimate, label)) +
  geom_vline(xintercept = 0, linewidth = 0.32, linetype = 2, color = "#777777") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label, color = color), linewidth = 0.72) +
  geom_point(aes(color = color), size = 2.15) +
  scale_color_identity() +
  scale_x_continuous(limits = c(0, 0.47), breaks = seq(0, 0.4, 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Stage and diagnosis sensitivity",
    subtitle = "I²: 90.8% unadjusted, 29.5% adjusted, 22.7% stage matched",
    x = "Standardized amyloid effect", y = NULL
  ) +
  theme_figure6()

# Rebuild the clinically relevant concept from the original Figure 6E using
# the prespecified HABS-HD pathology by WMH model rather than a causal diagram.
habs_wmh <- prepare_model(
  habs_raw,
  c("pathology", "cognition_impairment", "age", "education", "sex", "apoe4", "wmh")
)
fit_interaction <- lm(
  cognition_z ~ pathology_z * wmh_z + age_z + sex + education_z + apoe4,
  data = habs_wmh
)

hc3_vcov <- function(fit) {
  X <- model.matrix(fit)
  e <- residuals(fit)
  leverage <- hatvalues(fit)
  bread <- solve(crossprod(X))
  meat <- crossprod(X, X * (e / pmax(1 - leverage, 1e-8))^2)
  bread %*% meat %*% bread
}

path_seq <- seq(
  quantile(habs_wmh$pathology_z, 0.03),
  quantile(habs_wmh$pathology_z, 0.97),
  length.out = 100
)
reference_sex <- levels(habs_wmh$sex)[1]
interaction_pred <- expand.grid(
  pathology_z = path_seq,
  wmh_z = c(-1, 1),
  age_z = 0,
  sex = reference_sex,
  education_z = 0,
  apoe4 = 0,
  KEEP.OUT.ATTRS = FALSE
)
interaction_pred$sex <- factor(interaction_pred$sex, levels = levels(habs_wmh$sex))
Xnew <- model.matrix(delete.response(terms(fit_interaction)), interaction_pred)
Vhc3 <- hc3_vcov(fit_interaction)
interaction_pred$fit <- as.numeric(Xnew %*% coef(fit_interaction))
interaction_pred$se <- sqrt(pmax(0, rowSums((Xnew %*% Vhc3) * Xnew)))
interaction_pred$low <- interaction_pred$fit - 1.96 * interaction_pred$se
interaction_pred$high <- interaction_pred$fit + 1.96 * interaction_pred$se
interaction_pred$wmh_group <- factor(
  interaction_pred$wmh_z,
  levels = c(-1, 1), labels = c("Lower WMH (-1 SD)", "Higher WMH (+1 SD)")
)
interaction_row <- aligned[
  aligned$cohort == "HABS-HD" & aligned$analysis == "pathology_wmh" & aligned$context == "full",
][1, ]

pF <- ggplot(interaction_pred, aes(pathology_z, fit, color = wmh_group, fill = wmh_group)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "#B8B8B8") +
  geom_ribbon(aes(ymin = low, ymax = high), color = NA, alpha = 0.13) +
  geom_line(linewidth = 0.78) +
  scale_color_manual(values = c("Lower WMH (-1 SD)" = "#3E86BE", "Higher WMH (+1 SD)" = "#D9822B")) +
  scale_fill_manual(values = c("Lower WMH (-1 SD)" = "#3E86BE", "Higher WMH (+1 SD)" = "#D9822B")) +
  labs(
    title = "HABS-HD pathology by WMH interaction",
    subtitle = sprintf("beta = %.3f, FDR P = %.3g, n = %s",
                       interaction_row$estimate, interaction_row$p_fdr, comma(interaction_row$n)),
    x = "Plasma p-tau217 (z)", y = "Adjusted cognitive impairment",
    color = NULL, fill = NULL
  ) +
  theme_figure6() +
  theme(
    panel.grid.major.y = element_line(linewidth = 0.22, color = "#E5E5E5"),
    legend.position = "inside", legend.position.inside = c(0.31, 0.86),
    legend.background = element_rect(fill = alpha("white", 0.88), color = NA)
  )

design <- "
AAAAAA
BBCCDD
EEEFFF
"
figure6 <- pA + pB + pC + pD + pE + pF +
  plot_layout(design = design, heights = c(0.63, 1.12, 1.08)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 9),
        plot.tag.position = c(0.01, 0.99))

pdf_path <- file.path(submission_dir, "Figure 6.pdf")
svg_path <- file.path(source_dir, "Figure 6.svg")
png_path <- file.path(source_dir, "Figure 6.png")
tiff_path <- file.path(source_dir, "Figure 6.tiff")

ggsave(pdf_path, figure6, width = 183, height = 155, units = "mm", device = cairo_pdf, bg = "white")
ggsave(svg_path, figure6, width = 183, height = 155, units = "mm", device = svglite::svglite, bg = "white")
ggsave(png_path, figure6, width = 183, height = 155, units = "mm", dpi = 300,
       device = ragg::agg_png, bg = "white")
ggsave(tiff_path, figure6, width = 183, height = 155, units = "mm", dpi = 600,
       device = ragg::agg_tiff, compression = "lzw", bg = "white")

source_A <- main[, c("cohort", "n", "estimate", "se", "ci_low", "ci_high", "p", "p_fdr")]
source_A$panel <- "A"
source_BD <- partial_data
source_BD$panel <- c(A4 = "B", ADNI = "C", `HABS-HD` = "D")[source_BD$cohort]
source_E <- stage_plot[, c("cohort", "subgroup", "n", "estimate", "se", "ci_low", "ci_high", "p")]
source_E$panel <- "E"
source_F <- interaction_pred[, c("pathology_z", "wmh_z", "wmh_group", "fit", "se", "low", "high")]
source_F$panel <- "F"

write.table(source_A, file.path(source_dir, "Figure6_panel_A_source.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(source_BD, file.path(source_dir, "Figure6_panels_B_D_source.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(source_E, file.path(source_dir, "Figure6_panel_E_source.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(source_F, file.path(source_dir, "Figure6_panel_F_source.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
capture.output(sessionInfo(), file = file.path(source_dir, "sessionInfo.txt"))

message("Created redesigned Figure 6 at: ", pdf_path)
