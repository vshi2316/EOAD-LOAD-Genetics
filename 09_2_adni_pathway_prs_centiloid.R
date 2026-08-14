# Robustness analysis of replicated pathway PRS and amyloid burden in ADNI

options(stringsAsFactors = FALSE, scipen = 999)

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)
if (is.na(script_path)) {
  stop(
    "Run this analysis with source('09_2_adni_pathway_prs_centiloid.R').",
    call. = FALSE
  )
}

ROOT <- dirname(script_path)
RESULTS <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(ROOT, "results"))
INPUT <- Sys.getenv(
  "ADNI_PATHWAY_ANALYSIS_FILE",
  unset = file.path(ROOT, "data", "ADNI", "ADNI_pathway_PRS_Centiloid.tsv.gz")
)
OUT <- file.path(RESULTS, "adni_pathway_centiloid")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT)) {
  stop("The ADNI pathway PRS and Centiloid analysis file was not found.", call. = FALSE)
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table before running this script.", call. = FALSE)
}
library(data.table)

write_tsv <- function(x, name) {
  fwrite(x, file.path(OUT, name), sep = "\t", quote = FALSE, na = "NA")
}

safe_z <- function(x) {
  x <- as.numeric(x)
  sigma <- sd(x, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / sigma
}

winsorize <- function(x, probabilities = c(0.01, 0.99)) {
  limits <- quantile(x, probabilities, na.rm = TRUE, names = FALSE, type = 7)
  pmin(pmax(x, limits[1]), limits[2])
}

cluster_vcov <- function(fit, cluster) {
  X_full <- model.matrix(fit)
  beta <- coef(fit)
  keep <- is.finite(beta)
  X <- X_full[, keep, drop = FALSE]
  cluster <- droplevels(factor(cluster))
  n <- nrow(X)
  k <- ncol(X)
  g <- nlevels(cluster)
  if (g < 20L) {
    stop("Fewer than 20 research centers remain in a fitted model.", call. = FALSE)
  }

  if (inherits(fit, "glm")) {
    mu <- fitted(fit)
    score_residual <- model.response(model.frame(fit)) - mu
    weights <- as.numeric(fit$weights)
    bread <- tryCatch(
      solve(crossprod(X, X * weights)),
      error = function(e) qr.solve(
        crossprod(X, X * weights), diag(k), tol = 1e-10
      )
    )
  } else {
    score_residual <- residuals(fit)
    bread <- tryCatch(
      solve(crossprod(X)),
      error = function(e) qr.solve(crossprod(X), diag(k), tol = 1e-10)
    )
  }

  score_rows <- X * score_residual
  cluster_scores <- rowsum(score_rows, cluster, reorder = FALSE)
  meat <- crossprod(cluster_scores)
  correction <- (g / (g - 1)) * ((n - 1) / (n - k))
  variance <- correction * bread %*% meat %*% bread
  dimnames(variance) <- list(colnames(X), colnames(X))
  variance
}

cluster_table <- function(fit, cluster, model_id, focal_terms) {
  variance <- cluster_vcov(fit, cluster)
  beta <- coef(fit)[rownames(variance)]
  standard_error <- sqrt(pmax(diag(variance), 0))
  statistic <- beta / standard_error
  degrees_freedom <- nlevels(droplevels(factor(cluster))) - 1L
  p_value <- 2 * pt(-abs(statistic), df = degrees_freedom)
  critical <- qt(0.975, df = degrees_freedom)
  answer <- data.table(
    model = model_id,
    term = names(beta),
    estimate = as.numeric(beta),
    standard_error_cluster = standard_error,
    statistic = statistic,
    degrees_freedom = degrees_freedom,
    p_value = p_value,
    conf_low = as.numeric(beta) - critical * standard_error,
    conf_high = as.numeric(beta) + critical * standard_error,
    n = nobs(fit),
    centers = nlevels(droplevels(factor(cluster)))
  )
  answer[, focal := term %in% focal_terms]
  answer
}

fit_cluster_model <- function(data, formula, model_id, focal_terms,
                              family = NULL) {
  variables <- all.vars(formula)
  d <- copy(data[complete.cases(data[, ..variables]) & !is.na(site)])
  d[, site := droplevels(factor(site))]
  fit <- if (is.null(family)) {
    lm(formula, data = d)
  } else {
    glm(formula, data = d, family = family)
  }
  list(
    fit = fit,
    data = d,
    table = cluster_table(fit, d$site, model_id, focal_terms)
  )
}

if (!file.exists(INPUT)) stop("Missing input: ", INPUT, call. = FALSE)
d <- fread(INPUT, na.strings = c("", "NA", "N/A", "-99", "-999"))

prs_clearance <- "PRS_AmyloidBetaClearance_z"
prs_catabolism <- "PRS_NegativeRegulationAPPCatabolism_z"
focal_terms <- c(prs_clearance, prs_catabolism)
pc_covariates <- paste0("PC", 1:10, "_z")
covariates <- c(
  "age_z", "sex", "education_z", "apoe4_z", "diagnosis", pc_covariates
)
required <- unique(c(
  "RID", "Centiloid", "Centiloid_z", "site", focal_terms, covariates
))
absent <- setdiff(required, names(d))
if (length(absent)) {
  stop("Missing analysis fields: ", paste(absent, collapse = ", "),
       call. = FALSE)
}

d[, `:=`(
  RID = as.integer(RID),
  site = factor(site),
  sex = factor(sex),
  diagnosis = factor(diagnosis)
)]
for (field in c("Centiloid", "Centiloid_z", focal_terms, "age_z",
                "education_z", "apoe4_z", pc_covariates)) {
  d[, (field) := as.numeric(get(field))]
}
if (anyDuplicated(d$RID)) {
  stop("The robustness analysis expects one baseline record per RID.",
       call. = FALSE)
}

joint_predictors <- c(focal_terms, covariates)
joint_formula <- reformulate(joint_predictors, response = "Centiloid_z")

# Primary joint model: both pathway scores compete for the same amyloid outcome.
primary <- fit_cluster_model(
  d, joint_formula, "Primary_joint_PRS_Centiloid",
  focal_terms
)

# Sensitivity model without diagnosis avoids conditioning on clinical stage.
formula_without_diagnosis <- reformulate(
  c(focal_terms, setdiff(covariates, "diagnosis")),
  response = "Centiloid_z"
)
without_diagnosis <- fit_cluster_model(
  d, formula_without_diagnosis,
  "Sensitivity_joint_PRS_Centiloid_without_diagnosis",
  focal_terms
)

# Center fixed effects complement clustered inference.
site_fixed_formula <- reformulate(
  c(focal_terms, covariates, "site"), response = "Centiloid_z"
)
site_fixed_data <- copy(d[complete.cases(d[, c(
  "Centiloid_z", focal_terms, covariates, "site"
), with = FALSE])])
site_fixed_fit <- lm(site_fixed_formula, data = site_fixed_data)
site_fixed_table <- cluster_table(
  site_fixed_fit, site_fixed_data$site,
  "Sensitivity_joint_PRS_Centiloid_center_fixed_effects",
  focal_terms
)

# Amyloid positivity thresholds complement the continuous outcome.
d[, amyloid_positive_20 := as.integer(Centiloid > 20)]
d[, amyloid_positive_25 := as.integer(Centiloid > 25)]
binary_20 <- fit_cluster_model(
  d,
  reformulate(joint_predictors, response = "amyloid_positive_20"),
  "Sensitivity_amyloid_positive_CL20", focal_terms,
  family = binomial(link = "logit")
)
binary_25 <- fit_cluster_model(
  d,
  reformulate(joint_predictors, response = "amyloid_positive_25"),
  "Sensitivity_amyloid_positive_CL25", focal_terms,
  family = binomial(link = "logit")
)

# Winsorized and trimmed outcomes assess influence from extreme Centiloid values.
limits <- quantile(d$Centiloid, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
d[, Centiloid_winsorized_z := safe_z(winsorize(Centiloid))]
winsorized <- fit_cluster_model(
  d,
  reformulate(joint_predictors, response = "Centiloid_winsorized_z"),
  "Sensitivity_winsorized_Centiloid", focal_terms
)
trimmed <- fit_cluster_model(
  d[Centiloid >= limits[1] & Centiloid <= limits[2]],
  joint_formula,
  "Sensitivity_trimmed_Centiloid", focal_terms
)

all_tables <- rbindlist(list(
  primary$table,
  without_diagnosis$table,
  site_fixed_table,
  binary_20$table,
  binary_25$table,
  winsorized$table,
  trimmed$table
), fill = TRUE)
all_tables[focal == TRUE, p_fdr_within_model := p.adjust(p_value, method = "BH"),
           by = model]
all_tables[grepl("amyloid_positive", model), `:=`(
  odds_ratio = exp(estimate),
  odds_ratio_low = exp(conf_low),
  odds_ratio_high = exp(conf_high)
)]
write_tsv(all_tables, "pathway_prs_amyloid_robustness_models.tsv")

# Leave-one-center-out estimates show whether one center determines the result.
primary_data <- primary$data
center_results <- lapply(levels(primary_data$site), function(center) {
  z <- droplevels(primary_data[site != center])
  fit <- lm(joint_formula, data = z)
  coefficient <- coef(summary(fit))
  rbindlist(lapply(focal_terms, function(term) {
    data.table(
      omitted_center = center,
      pathway_term = term,
      estimate = coefficient[term, "Estimate"],
      standard_error = coefficient[term, "Std. Error"],
      p_value = coefficient[term, "Pr(>|t|)"],
      n = nrow(z),
      centers_remaining = nlevels(z$site)
    )
  }))
})
center_results <- rbindlist(center_results)
write_tsv(center_results, "leave_one_center_out_estimates.tsv")

primary_focal <- primary$table[focal == TRUE]
loo_summary <- center_results[, .(
  omissions = .N,
  minimum_estimate = min(estimate),
  maximum_estimate = max(estimate),
  same_direction_fraction = mean(sign(estimate) == sign(
    primary_focal$estimate[match(pathway_term[1], primary_focal$term)]
  )),
  nominal_p_below_0_05_fraction = mean(p_value < 0.05)
), by = pathway_term]
write_tsv(loo_summary, "leave_one_center_out_summary.tsv")

focal_summary <- all_tables[focal == TRUE, .(
  model, term, estimate, standard_error_cluster, p_value,
  p_fdr_within_model, conf_low, conf_high, odds_ratio,
  odds_ratio_low, odds_ratio_high, n, centers
)]
write_tsv(focal_summary, "focal_robustness_summary.tsv")

decision <- data.table(
  criterion = c(
    "Catabolism pathway is independently associated with continuous Centiloid",
    "Catabolism association survives center-clustered inference",
    "Catabolism direction is unchanged after every center omission",
    "Clearance pathway is independently associated with continuous Centiloid"
  ),
  pass = c(
    primary_focal[term == prs_catabolism, p_value] < 0.05,
    primary_focal[term == prs_catabolism, conf_low] > 0,
    loo_summary[pathway_term == prs_catabolism, same_direction_fraction] == 1,
    primary_focal[term == prs_clearance, p_value] < 0.05
  )
)
write_tsv(decision, "robustness_evidence_gates.tsv")

capture.output(
  {
    cat("Pathway PRS and amyloid robustness analysis\n\n")
    cat("Primary mutually adjusted model\n")
    print(primary_focal[, .(
      term, estimate, standard_error_cluster, conf_low, conf_high,
      p_value, n, centers
    )])
    cat("\nFocal sensitivity results\n")
    print(focal_summary)
    cat("\nLeave-one-center-out summary\n")
    print(loo_summary)
    cat("\nEvidence gates\n")
    print(decision)
  },
  file = file.path(OUT, "analysis_report.txt")
)

saveRDS(
  list(
    primary = primary$fit,
    without_diagnosis = without_diagnosis$fit,
    center_fixed_effects = site_fixed_fit,
    amyloid_positive_20 = binary_20$fit,
    amyloid_positive_25 = binary_25$fit,
    winsorized = winsorized$fit,
    trimmed = trimmed$fit
  ),
  file.path(OUT, "robustness_model_objects.rds")
)
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))

cat("\nPathway PRS and amyloid robustness analysis completed.\n")
cat("Results:", OUT, "\n")
