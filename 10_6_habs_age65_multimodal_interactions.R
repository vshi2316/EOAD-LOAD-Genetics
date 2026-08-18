# HABS-HD age modification of prespecified multimodal associations

options(stringsAsFactors = FALSE, scipen = 999)

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)
ROOT <- if (!is.na(script_path)) dirname(script_path) else {
  Sys.getenv("EOAD_CODE_DIR", unset = getwd())
}
RESULTS <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(ROOT, "results"))
OUT <- file.path(RESULTS, "habs_age65_multimodal_interactions")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table before running this script.", call. = FALSE)
}
library(data.table)

input_candidates <- c(
  Sys.getenv("HABS_ANALYSIS_FILE", unset = NA_character_),
  file.path(RESULTS, "HABS_Submission_Freeze", "habs_submission_analysis_dataset.tsv"),
  file.path(RESULTS, "habs_primary_models", "habs_primary_analysis_dataset.tsv"),
  file.path(ROOT, "data", "HABS_HD", "habs_submission_analysis_dataset.tsv")
)
input_candidates <- unique(input_candidates[!is.na(input_candidates) & nzchar(input_candidates)])
INPUT <- input_candidates[file.exists(input_candidates)][1L]
if (!length(INPUT) || is.na(INPUT)) {
  stop("The frozen HABS-HD analysis dataset was not found.", call. = FALSE)
}

write_tsv <- function(x, name) {
  fwrite(x, file.path(OUT, name), sep = "\t", quote = FALSE, na = "NA")
}

as_numeric_clean <- function(x) {
  x <- trimws(gsub(",", "", as.character(x), fixed = TRUE))
  x[x %in% c("", "NA", "N/A", ".", "-99", "-999", "-9999")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

hc3_table <- function(fit, model_id, role, focal_term) {
  beta <- coef(fit)
  keep <- is.finite(beta)
  beta <- beta[keep]
  x <- model.matrix(fit)[, keep, drop = FALSE]
  leverage <- pmin(hatvalues(fit), 1 - 1e-8)
  residual_hc3 <- residuals(fit) / pmax(1 - leverage, 1e-8)
  cross_product <- crossprod(x)
  bread <- tryCatch(
    solve(cross_product),
    error = function(e) qr.solve(cross_product, diag(ncol(x)), tol = 1e-10)
  )
  variance <- bread %*% crossprod(x, x * residual_hc3^2) %*% bread
  standard_error <- sqrt(pmax(diag(variance), 0))
  statistic <- beta / standard_error
  data.table(
    model = model_id,
    role = role,
    term = names(beta),
    estimate = as.numeric(beta),
    standard_error_hc3 = standard_error,
    statistic = statistic,
    p_value = 2 * pnorm(-abs(statistic)),
    conf_low = as.numeric(beta) - qnorm(0.975) * standard_error,
    conf_high = as.numeric(beta) + qnorm(0.975) * standard_error,
    n = nobs(fit),
    focal = names(beta) == focal_term,
    covariance = "HC3"
  )
}

fit_interaction <- function(data, outcome, exposure, modifier, covariates,
                            model_id, role) {
  variables <- unique(c(outcome, exposure, modifier, covariates))
  model_data <- copy(data[complete.cases(data[, ..variables])])
  if (!nrow(model_data)) stop("No complete cases for ", model_id, call. = FALSE)
  covariates <- covariates[vapply(
    covariates,
    function(variable) uniqueN(model_data[[variable]]) > 1L,
    logical(1)
  )]
  modifier_expression <- paste0(exposure, " * ", modifier)
  formula <- as.formula(paste(
    outcome, "~", modifier_expression,
    if (length(covariates)) paste("+", paste(covariates, collapse = " + ")) else ""
  ))
  fit <- lm(formula, data = model_data)
  interaction_candidates <- grep(
    paste0("^", exposure, ":", modifier, "|^", modifier, ".*:", exposure, "$"),
    names(coef(fit)), value = TRUE
  )
  if (length(interaction_candidates) != 1L) {
    stop("The focal interaction term could not be identified for ", model_id, call. = FALSE)
  }
  list(
    table = hc3_table(fit, model_id, role, interaction_candidates),
    data = model_data,
    formula = paste(deparse(formula), collapse = " "),
    focal_term = interaction_candidates
  )
}

fit_stratified <- function(data, outcome, exposure, covariates, model_id, group_name) {
  variables <- unique(c(outcome, exposure, covariates))
  model_data <- copy(data[complete.cases(data[, ..variables])])
  covariates <- covariates[vapply(
    covariates,
    function(variable) uniqueN(model_data[[variable]]) > 1L,
    logical(1)
  )]
  formula <- reformulate(c(exposure, covariates), response = outcome)
  fit <- lm(formula, data = model_data)
  result <- hc3_table(fit, model_id, "stratified_description", exposure)[focal == TRUE]
  result[, age_group := group_name]
  result
}

analysis <- fread(INPUT, na.strings = c("", "NA", "N/A", "NaN", "-9999"))
required_fields <- c(
  "age", "age_z", "amyloid_z", "tau_z", "structure_z", "baseline_global_z",
  "sex", "education_z", "ethnicity", "apoe4", "scanner", "diagnosis_factor"
)
missing_fields <- setdiff(required_fields, names(analysis))
if (length(missing_fields)) {
  stop("Required HABS-HD fields are missing: ", paste(missing_fields, collapse = ", "), call. = FALSE)
}

analysis[, `:=`(
  age = as_numeric_clean(age),
  age_z = as_numeric_clean(age_z),
  amyloid_z = as_numeric_clean(amyloid_z),
  tau_z = as_numeric_clean(tau_z),
  structure_z = as_numeric_clean(structure_z),
  baseline_global_z = as_numeric_clean(baseline_global_z),
  education_z = as_numeric_clean(education_z),
  apoe4 = as_numeric_clean(apoe4),
  sex = factor(sex),
  ethnicity = factor(ethnicity),
  scanner = factor(scanner),
  diagnosis_factor = factor(diagnosis_factor)
)]
analysis[, age_group := factor(
  fifelse(age < 65, "Younger_under_65", "Older_65_plus"),
  levels = c("Older_65_plus", "Younger_under_65")
)]

common_covariates <- c(
  "age_z", "sex", "education_z", "ethnicity", "apoe4", "scanner"
)
specifications <- list(
  list(
    id = "amyloid_to_tau", outcome = "tau_z", exposure = "amyloid_z",
    covariates = common_covariates
  ),
  list(
    id = "tau_to_temporal_structure", outcome = "structure_z", exposure = "tau_z",
    covariates = c("amyloid_z", common_covariates)
  ),
  list(
    id = "temporal_structure_to_cognition", outcome = "baseline_global_z",
    exposure = "structure_z", covariates = c("tau_z", "amyloid_z", common_covariates)
  )
)

group_tables <- list()
continuous_tables <- list()
stratified_tables <- list()
sensitivity_tables <- list()
manifest <- list()
model_sample_audit <- list()

for (specification in specifications) {
  grouped <- fit_interaction(
    analysis, specification$outcome, specification$exposure, "age_group",
    specification$covariates,
    paste0("age65_group_interaction_", specification$id),
    "primary_age_group_interaction"
  )
  group_tables[[length(group_tables) + 1L]] <- grouped$table
  manifest[[length(manifest) + 1L]] <- data.table(
    model = paste0("age65_group_interaction_", specification$id),
    role = "primary_age_group_interaction",
    formula = grouped$formula,
    focal_term = grouped$focal_term,
    n = nrow(grouped$data)
  )
  model_sample_audit[[length(model_sample_audit) + 1L]] <- grouped$data[, .N, by = age_group][
    , model := paste0("age65_group_interaction_", specification$id)
  ][, .(model, age_group, n = N)]

  continuous <- fit_interaction(
    analysis, specification$outcome, specification$exposure, "age_z",
    setdiff(specification$covariates, "age_z"),
    paste0("continuous_age_interaction_", specification$id),
    "continuous_age_sensitivity"
  )
  continuous_tables[[length(continuous_tables) + 1L]] <- continuous$table
  manifest[[length(manifest) + 1L]] <- data.table(
    model = paste0("continuous_age_interaction_", specification$id),
    role = "continuous_age_sensitivity",
    formula = continuous$formula,
    focal_term = continuous$focal_term,
    n = nrow(continuous$data)
  )

  diagnosis_adjusted <- fit_interaction(
    analysis, specification$outcome, specification$exposure, "age_group",
    c(specification$covariates, "diagnosis_factor"),
    paste0("diagnosis_adjusted_age65_interaction_", specification$id),
    "diagnosis_adjusted_sensitivity"
  )
  sensitivity_tables[[length(sensitivity_tables) + 1L]] <- diagnosis_adjusted$table
  manifest[[length(manifest) + 1L]] <- data.table(
    model = paste0("diagnosis_adjusted_age65_interaction_", specification$id),
    role = "diagnosis_adjusted_sensitivity",
    formula = diagnosis_adjusted$formula,
    focal_term = diagnosis_adjusted$focal_term,
    n = nrow(diagnosis_adjusted$data)
  )

  vida1 <- fit_interaction(
    analysis[scanner == "Vida1"], specification$outcome, specification$exposure,
    "age_group", setdiff(specification$covariates, "scanner"),
    paste0("vida1_age65_interaction_", specification$id),
    "vida1_sensitivity"
  )
  sensitivity_tables[[length(sensitivity_tables) + 1L]] <- vida1$table
  manifest[[length(manifest) + 1L]] <- data.table(
    model = paste0("vida1_age65_interaction_", specification$id),
    role = "vida1_sensitivity",
    formula = vida1$formula,
    focal_term = vida1$focal_term,
    n = nrow(vida1$data)
  )

  for (group_name in levels(analysis$age_group)) {
    stratified_tables[[length(stratified_tables) + 1L]] <- fit_stratified(
      analysis[age_group == group_name], specification$outcome,
      specification$exposure, specification$covariates,
      paste0("stratified_", specification$id), group_name
    )
  }
}

group_effects <- rbindlist(group_tables, fill = TRUE)
continuous_effects <- rbindlist(continuous_tables, fill = TRUE)
stratified_effects <- rbindlist(stratified_tables, fill = TRUE)
sensitivity_effects <- rbindlist(sensitivity_tables, fill = TRUE)
group_effects[focal == TRUE, p_fdr_primary_interactions := p.adjust(p_value, "BH")]
continuous_effects[focal == TRUE, p_fdr_continuous_interactions := p.adjust(p_value, "BH")]
sensitivity_effects[focal == TRUE, p_fdr_sensitivity := p.adjust(p_value, "BH"), by = role]

primary_focal <- group_effects[focal == TRUE]
continuous_focal <- continuous_effects[focal == TRUE]
sensitivity_focal <- sensitivity_effects[focal == TRUE]
primary_focal[, model_key := sub("^age65_group_interaction_", "", model)]
continuous_focal[, model_key := sub("^continuous_age_interaction_", "", model)]
decision <- merge(
  primary_focal[, .(
    model, model_key,
    group_interaction_estimate = estimate,
    group_interaction_p = p_value,
    group_interaction_fdr = p_fdr_primary_interactions
  )],
  continuous_focal[, .(
    model_continuous = model,
    model_key,
    continuous_interaction_estimate = estimate,
    continuous_interaction_p = p_value,
    continuous_interaction_fdr = p_fdr_continuous_interactions
  )],
  by = "model_key", all.x = TRUE, sort = FALSE
)
decision[, continuous_direction_consistent := fifelse(
  is.finite(group_interaction_estimate) & is.finite(continuous_interaction_estimate),
  sign(group_interaction_estimate) == -sign(continuous_interaction_estimate),
  FALSE
)]
diagnosis_focal <- sensitivity_focal[role == "diagnosis_adjusted_sensitivity"]
diagnosis_focal[, model_key := sub("^diagnosis_adjusted_age65_interaction_", "", model)]
vida1_focal <- sensitivity_focal[role == "vida1_sensitivity"]
vida1_focal[, model_key := sub("^vida1_age65_interaction_", "", model)]
decision <- merge(
  decision,
  diagnosis_focal[, .(
    model_key,
    diagnosis_adjusted_estimate = estimate,
    diagnosis_adjusted_p = p_value,
    diagnosis_adjusted_fdr = p_fdr_sensitivity
  )],
  by = "model_key", all.x = TRUE, sort = FALSE
)
decision <- merge(
  decision,
  vida1_focal[, .(
    model_key,
    vida1_estimate = estimate,
    vida1_p = p_value,
    vida1_fdr = p_fdr_sensitivity
  )],
  by = "model_key", all.x = TRUE, sort = FALSE
)
decision[, sensitivity_direction_consistent :=
  sign(group_interaction_estimate) == sign(diagnosis_adjusted_estimate) &
  sign(group_interaction_estimate) == sign(vida1_estimate)]
decision[, eligible_for_main_text :=
  group_interaction_fdr < 0.05 & continuous_direction_consistent &
    sensitivity_direction_consistent]
decision[, interpretation := fifelse(
  eligible_for_main_text,
  "Age modification is supported by the prespecified group interaction and continuous-age sensitivity.",
  "The result does not meet the prespecified criterion for a main-text age-modification claim."
)]

sample_audit <- analysis[, .(
  participants = .N,
  age_min = min(age, na.rm = TRUE),
  age_median = median(age, na.rm = TRUE),
  age_max = max(age, na.rm = TRUE)
), by = age_group]
ethnicity_audit <- analysis[, .N, by = .(age_group, ethnicity)]
diagnosis_audit <- analysis[, .N, by = .(age_group, diagnosis_factor)]
scanner_audit <- analysis[, .N, by = .(age_group, scanner)]

write_tsv(group_effects, "age65_group_interaction_effects.tsv")
write_tsv(continuous_effects, "continuous_age_interaction_effects.tsv")
write_tsv(stratified_effects, "age65_stratified_effects.tsv")
write_tsv(sensitivity_effects, "age65_interaction_sensitivity_effects.tsv")
write_tsv(rbindlist(manifest, fill = TRUE), "model_manifest.tsv")
write_tsv(rbindlist(model_sample_audit, fill = TRUE), "model_sample_audit.tsv")
write_tsv(decision, "decision.tsv")
write_tsv(sample_audit, "sample_audit.tsv")
write_tsv(ethnicity_audit, "ethnicity_audit.tsv")
write_tsv(diagnosis_audit, "diagnosis_audit.tsv")
write_tsv(scanner_audit, "scanner_audit.tsv")

capture.output({
  cat("HABS-HD age modification of prespecified multimodal associations\n\n")
  cat("Participants younger than 65 years are an age-defined community subgroup and are not classified as EOAD cases.\n\n")
  cat("Sample audit\n")
  print(sample_audit)
  cat("\nPrimary age-group interactions\n")
  print(primary_focal)
  cat("\nContinuous-age interaction sensitivity\n")
  print(continuous_focal)
  cat("\nAge-stratified descriptive effects\n")
  print(stratified_effects)
  cat("\nDiagnosis-adjusted and Vida1 sensitivity interactions\n")
  print(sensitivity_focal)
  cat("\nDecision\n")
  print(decision)
}, file = file.path(OUT, "analysis_report.txt"))

writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("\nHABS-HD age interaction analysis completed.\nResults: ", OUT, "\n", sep = "")
