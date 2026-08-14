# HABS-HD temporal vulnerability map and primary multimodal evidence figure
# Regional models decompose the prespecified eight-region cortical score.

options(stringsAsFactors = FALSE, scipen = 999)
ofile <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
if (is.na(ofile)) stop("Run with source('10_5_habs_regional_analysis.R').", call. = FALSE)
ROOT <- dirname(ofile)
DATA <- Sys.getenv("HABS_DATA_DIR", unset = file.path(ROOT, "data", "HABS_HD"))
RESULTS <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(ROOT, "results"))
PRIMARY <- file.path(RESULTS, "habs_primary_models")
INPUT <- file.path(PRIMARY, "habs_primary_analysis_dataset.tsv")
PRIMARY_EFFECTS <- file.path(PRIMARY, "habs_primary_model_effects_hc3.tsv")
MRI_FILE <- Sys.getenv("HABS_CORTICAL_THICKNESS_FILE", unset = file.path(DATA, "cortical_thickness.xlsx"))
OUT <- file.path(RESULTS, "habs_regional_analysis")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("data.table", "readxl", "ggplot2", "ggseg", "patchwork", "scales", "svglite", "ragg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
library(data.table)
library(ggplot2)
library(ggseg)
library(patchwork)

write_tsv <- function(x, name) fwrite(x, file.path(OUT, name), sep = "\t", quote = FALSE, na = "NA")
as_numeric_clean <- function(x) {
  x <- trimws(gsub(",", "", as.character(x), fixed = TRUE))
  x[x %in% c("", "NA", "N/A", ".", "-9999", "-8888", "-777777")] <- NA_character_
  suppressWarnings(as.numeric(x))
}
month_date <- function(year, month) {
  y <- suppressWarnings(as.integer(as.character(year)))
  m <- suppressWarnings(as.integer(as.character(month)))
  out <- rep(as.IDate(NA), length(y))
  ok <- is.finite(y) & is.finite(m) & m >= 1L & m <= 12L
  out[ok] <- as.IDate(sprintf("%04d-%02d-15", y[ok], m[ok]))
  out
}

analysis <- fread(INPUT, na.strings = c("", "NA", "NaN"))
primary_effects <- fread(PRIMARY_EFFECTS, na.strings = c("", "NA", "NaN"))
analysis[, baseline_mri_date := as.IDate(baseline_mri_date)]
analysis[, `:=`(
  sex = factor(sex), ethnicity = factor(ethnicity), scanner = factor(scanner),
  apoe4 = as_numeric_clean(apoe4), age_z = as_numeric_clean(age_z),
  education_z = as_numeric_clean(education_z), tau_z = as_numeric_clean(tau_z),
  amyloid_z = as_numeric_clean(amyloid_z)
)]

# The eight regions were specified before regional effect estimation.
region_dictionary <- data.table(
  field = c(
    "L_entorhinal_thickavg", "R_entorhinal_thickavg",
    "L_parahippocampal_thickavg", "R_parahippocampal_thickavg",
    "L_inferiortemporal_thickavg", "R_inferiortemporal_thickavg",
    "L_middletemporal_thickavg", "R_middletemporal_thickavg"
  ),
  hemisphere = rep(c("left", "right"), 4),
  region_key = rep(c("entorhinal", "parahippocampal", "inferiortemporal", "middletemporal"), each = 2),
  label = c(
    "Left entorhinal", "Right entorhinal",
    "Left parahippocampal", "Right parahippocampal",
    "Left inferior temporal", "Right inferior temporal",
    "Left middle temporal", "Right middle temporal"
  ),
  order = 1:8
)

# Read the processed cortical-thickness inclusion sheets and recover the exact
# baseline MRI record used in the primary analysis.
mri_sheets <- readxl::excel_sheets(MRI_FILE)
mri_sheets <- mri_sheets[grepl("inclusion", tolower(mri_sheets)) & !grepl("exclusion|pending", tolower(mri_sheets))]
mri <- rbindlist(lapply(mri_sheets, function(sheet) {
  x <- as.data.table(suppressWarnings(readxl::read_excel(MRI_FILE, sheet = sheet)))
  x[, source_sheet := sheet]
  x
}), fill = TRUE)
if (!all(region_dictionary$field %in% names(mri))) stop("One or more prespecified cortical fields are missing.", call. = FALSE)
mri[, `:=`(Med_ID = as.integer(Med_ID), baseline_mri_date = month_date(MRI_year, MRI_month))]
for (field in region_dictionary$field) mri[, (field) := as_numeric_clean(get(field))]
setorder(mri, Med_ID, baseline_mri_date, source_sheet)
mri <- mri[, .SD[1L], by = .(Med_ID, baseline_mri_date)]

regional <- merge(
  analysis[, .(Med_ID, baseline_mri_date, tau_z, amyloid_z, age_z, sex, education_z, ethnicity, apoe4, scanner)],
  mri[, c("Med_ID", "baseline_mri_date", region_dictionary$field), with = FALSE],
  by = c("Med_ID", "baseline_mri_date"), all.x = TRUE
)

# Standardize thickness across the full eligible baseline MRI cohort, then
# reverse its sign so higher values consistently indicate greater vulnerability.
for (field in region_dictionary$field) {
  values <- regional[[field]]
  regional[, (paste0(field, "_vulnerability_z")) := -(values - mean(values, na.rm = TRUE)) / sd(values, na.rm = TRUE)]
}

hc3_term <- function(fit, term) {
  beta <- coef(fit); keep <- is.finite(beta); beta <- beta[keep]
  X <- model.matrix(fit)[, keep, drop = FALSE]
  adjusted_residual <- residuals(fit) / pmax(1 - hatvalues(fit), 1e-8)
  bread <- tryCatch(solve(crossprod(X)), error = function(e) qr.solve(crossprod(X), diag(ncol(X)), tol = 1e-10))
  variance <- bread %*% crossprod(X, X * adjusted_residual^2) %*% bread
  se <- sqrt(pmax(diag(variance), 0))
  index <- match(term, names(beta))
  data.table(
    estimate = as.numeric(beta[index]), standard_error_hc3 = se[index],
    conf_low = as.numeric(beta[index]) - qnorm(.975) * se[index],
    conf_high = as.numeric(beta[index]) + qnorm(.975) * se[index],
    p_value = 2 * pnorm(-abs(beta[index] / se[index])), n = nobs(fit)
  )
}

regional_effects <- vector("list", nrow(region_dictionary))
regional_covariates <- c("amyloid_z", "age_z", "sex", "education_z", "ethnicity", "apoe4", "scanner")
for (i in seq_len(nrow(region_dictionary))) {
  outcome <- paste0(region_dictionary$field[i], "_vulnerability_z")
  variables <- c(outcome, "tau_z", regional_covariates)
  model_data <- regional[which(complete.cases(regional[, variables, with = FALSE]))]
  varying_covariates <- regional_covariates[vapply(regional_covariates, function(v) uniqueN(model_data[[v]]) > 1L, logical(1))]
  fit <- lm(reformulate(c("tau_z", varying_covariates), response = outcome), data = model_data)
  regional_effects[[i]] <- cbind(region_dictionary[i], hc3_term(fit, "tau_z"))
}
regional_effects <- rbindlist(regional_effects)
regional_effects[, p_fdr_eight_regions := p.adjust(p_value, method = "BH")]
regional_effects[, significant_fdr := p_fdr_eight_regions < 0.05]
write_tsv(regional_effects, "source_data_regional_tau_effects.tsv")

# Resolve region names against the installed Desikan-Killiany atlas instead of
# relying on package-version-specific spacing in labels.
normalize_region <- function(x) gsub("[^a-z]", "", tolower(as.character(x)))
atlas_data <- as.data.table(as.data.frame(ggseg::dk))
atlas_regions <- unique(na.omit(atlas_data$region))
atlas_lookup <- data.table(region = atlas_regions, region_key = normalize_region(atlas_regions))
regional_effects[, region_key := normalize_region(region_key)]
regional_effects <- merge(regional_effects, atlas_lookup, by = "region_key", all.x = TRUE)
if (anyNA(regional_effects$region)) {
  stop("DK atlas mapping failed for: ", paste(regional_effects[is.na(region), unique(region_key)], collapse = ", "), call. = FALSE)
}
brain_source <- regional_effects[, .(region, hemi = hemisphere, estimate)]
write_tsv(brain_source, "source_data_brain_surface.tsv")

primary_plot_data <- primary_effects[focal == TRUE & role == "primary", .(
  model, estimate, conf_low, conf_high, p_value, p_fdr_primary, n
)]
primary_plot_data[, path := factor(model,
  levels = c("primary_structure_to_concurrent_cognition", "primary_tau_to_structure", "primary_amyloid_to_tau"),
  labels = c("Temporal structure to cognition", "Tau to temporal structure", "Amyloid to tau")
)]
write_tsv(primary_plot_data, "source_data_primary_paths.tsv")

palette <- c(negative = "#3B6FB6", zero = "#F2F2F2", positive = "#B44A4A", ink = "#252525", muted = "#777777")
theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = palette["ink"]),
      axis.ticks = element_line(linewidth = 0.35, colour = palette["ink"]),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      axis.title = element_text(size = base_size), axis.text = element_text(size = base_size - 0.3),
      legend.title = element_text(size = base_size - 0.2), legend.text = element_text(size = base_size - 0.5),
      panel.grid = element_blank(), plot.margin = margin(4, 5, 4, 5)
    )
}

map_limit <- max(regional_effects$estimate, na.rm = TRUE)
p_brain <- ggplot(brain_source, aes(fill = estimate)) +
  ggseg::geom_brain(
    atlas = ggseg::dk, position = ggseg::position_brain(hemi ~ side),
    colour = "white", size = 0.18, show.legend = TRUE
  ) +
  scale_fill_gradient(
    low = "#F3E6E3", high = palette["positive"],
    limits = c(0, map_limit), oob = scales::squish,
    na.value = "#E6E6E6", name = expression(beta)
  ) +
  ggplot2::coord_sf(clip = "off", datum = NA) +
  labs(title = "Prespecified temporal vulnerability map") +
  theme_void(base_family = "Arial", base_size = 7) +
  theme(
    plot.title = element_text(size = 7.5, face = "bold", hjust = 0),
    legend.position = "bottom", legend.key.width = grid::unit(15, "mm"),
    legend.key.height = grid::unit(2.5, "mm"), plot.margin = margin(4, 5, 2, 5)
  )

regional_effects[, label_factor := factor(label, levels = rev(region_dictionary$label))]
p_forest <- ggplot(regional_effects, aes(x = estimate, y = label_factor)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#9A9A9A") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0, linewidth = 0.55, colour = palette["ink"]) +
  geom_point(aes(fill = significant_fdr), shape = 21, size = 2.1, stroke = 0.35, colour = palette["ink"]) +
  scale_fill_manual(values = c(`TRUE` = palette["positive"], `FALSE` = "white"), guide = "none") +
  labs(title = "Regional tau associations", x = expression("Standardized " * beta * " (tau to vulnerability)"), y = NULL) +
  theme_pub() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

primary_plot_data[, path_colour := factor(path, levels = levels(path))]
primary_plot_data[, fdr_label := ifelse(
  p_fdr_primary < 1e-4,
  sprintf("n=%d; FDR P<0.0001", n),
  sprintf("n=%d; FDR P=%.4f", n, p_fdr_primary)
)]
p_paths <- ggplot(primary_plot_data, aes(x = estimate, y = path)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#9A9A9A") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0, linewidth = 0.65, colour = palette["ink"]) +
  geom_point(aes(fill = path_colour), shape = 21, size = 2.6, stroke = 0.4, colour = palette["ink"]) +
  scale_fill_manual(values = c("#3E8C8C", "#D08A36", "#B44A4A"), guide = "none") +
  geom_text(
    aes(x = conf_high + 0.025, label = fdr_label),
    hjust = 0, size = 2.15, family = "Arial"
  ) +
  scale_x_continuous(limits = c(-0.26, 0.72), breaks = c(-0.2, 0, 0.2, 0.4, 0.6), expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "Primary multimodal associations", x = "Standardized coefficient (95% CI)", y = NULL) +
  theme_pub() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

figure <- p_brain + p_forest + p_paths +
  patchwork::plot_layout(design = "AAB\nAAB\nCCC", heights = c(1, 1, 0.85)) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8.5, face = "bold", family = "Arial"))

base_name <- file.path(OUT, "Figure_HABS_temporal_vulnerability")
width_in <- 183 / 25.4
height_in <- 132 / 25.4
svglite::svglite(paste0(base_name, ".svg"), width = width_in, height = height_in)
print(figure); dev.off()
grDevices::cairo_pdf(paste0(base_name, ".pdf"), width = width_in, height = height_in, family = "Arial")
print(figure); dev.off()
ragg::agg_tiff(paste0(base_name, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
print(figure); dev.off()
ragg::agg_png(paste0(base_name, ".png"), width = width_in, height = height_in, units = "in", res = 300)
print(figure); dev.off()

writeLines(c(
  "Figure contract",
  "Core conclusion: Amyloid and tau converge on a prespecified temporal structural vulnerability phenotype associated with concurrent cognition.",
  "Panel a: Regional standardized tau coefficients from eight prespecified bilateral cortical regions.",
  "Panel b: HC3 95% confidence intervals for the same regional coefficients; filled points survive BH correction across eight regions.",
  "Panel c: Three prespecified primary standardized associations with HC3 confidence intervals and BH correction across the primary family.",
  "Regional models adjust for amyloid, age, sex, education, ethnicity, APOE4, and MRI scanner.",
  "The brain map is a regional decomposition of the prespecified score and does not represent a whole-cortex discovery scan."
), file.path(OUT, "figure_contract_and_legend_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("\nHABS-HD temporal vulnerability figure completed.\nResults: ", OUT, "\n", sep = "")
