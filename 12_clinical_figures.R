# Figures for the ADNI genetic bridge and cross-cohort clinical expression

options(stringsAsFactors = FALSE, scipen = 999)
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) stop("Run with source('12_clinical_figures.R').", call. = FALSE)
ROOT <- dirname(script_path)
RESULTS <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(ROOT, "results"))
OUT <- file.path(RESULTS, "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

packages <- c("data.table", "ggplot2", "patchwork")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "), call. = FALSE)
library(data.table); library(ggplot2); library(patchwork)

theme_manuscript <- function() {
  theme_classic(base_family = "Arial", base_size = 8) +
    theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_text(size = 9, face = "plain"), legend.title = element_blank())
}
save_plot <- function(plot, name, width, height) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".tiff")), plot, width = width, height = height, units = "in", dpi = 600, compression = "lzw")
}

adni_dir <- file.path(RESULTS, "adni_pathway_centiloid")
adni_effects <- fread(file.path(adni_dir, "pathway_prs_amyloid_robustness_models.tsv"))
adni_loo <- fread(file.path(adni_dir, "leave_one_center_out_estimates.tsv"))
focal <- adni_effects[term %chin% c("PRS_AmyloidBetaClearance_z", "PRS_NegativeRegulationAPPCatabolism_z")]
focal[, label := fifelse(term == "PRS_AmyloidBetaClearance_z", "Amyloid beta clearance", "Negative regulation of APP catabolism")]
p4a <- ggplot(focal, aes(estimate, reorder(label, estimate))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18, linewidth = 0.45) +
  geom_point(size = 2.1, colour = "#B4473D") + labs(x = "Standardized effect on Centiloid", y = NULL, title = "Pathway PRS associations") + theme_manuscript()
site_field <- intersect(c("omitted_center", "excluded_site"), names(adni_loo))[1]
if (is.na(site_field)) stop("The leave-one-center-out identifier was not found.", call. = FALSE)
p4b <- ggplot(adni_loo, aes(estimate, reorder(as.factor(get(site_field)), estimate))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_point(size = 1.1, colour = "#3E6F8F") + labs(x = "Leave-one-center-out effect", y = "Research center", title = "Center sensitivity") + theme_manuscript()
save_plot(p4a + p4b + plot_annotation(tag_levels = "A"), "Figure_4", 7.2, 4.6)

clinical_dir <- file.path(RESULTS, "clinical_expression")
clinical <- fread(file.path(clinical_dir, "clinical_expression_primary_effects.tsv"))
stage <- fread(file.path(clinical_dir, "adni_diagnostic_stage_effects.tsv"))
habs <- fread(file.path(clinical_dir, "habs_structural_expression_effects.tsv"))
p6a <- ggplot(clinical, aes(estimate, reorder(cohort, estimate), colour = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.45) +
  geom_point(size = 2.1) + labs(x = "Standardized pathology effect", y = NULL, title = "Clinical expression by cohort") + theme_manuscript() + guides(colour = "none")
p6b <- ggplot(stage, aes(estimate, reorder(model, estimate))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.45) +
  geom_point(size = 2.1, colour = "#B4473D") + labs(x = "Standardized Centiloid effect", y = NULL, title = "ADNI diagnostic settings") + theme_manuscript()
habs_plot <- habs[focal == TRUE]
p6c <- ggplot(habs_plot, aes(estimate, reorder(model, estimate))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.45) +
  geom_point(size = 2.1, colour = "#3A8A63") + labs(x = "Standardized effect", y = NULL, title = "HABS-HD structural expression") + theme_manuscript()
save_plot((p6a | p6b) / p6c + plot_annotation(tag_levels = "A"), "Figure_6", 7.2, 6.0)

writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
