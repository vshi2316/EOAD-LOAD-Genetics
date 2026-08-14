# Execute the manuscript analyses in evidence order.

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (file.exists(file.path(root, "config.R"))) source(file.path(root, "config.R"), encoding = "UTF-8")

modules <- c(
  "00_preflight.R",
  "01_genetic_architecture.R",
  "02_polygenic_and_transcriptomic_analysis.R",
  "03_independent_pathway_replication.R",
  "04_variant_replication.R",
  "05_meta_enhanced_analysis.R",
  "06_genetic_robustness.R",
  "07_white_matter_analysis.R",
  "08_single_nucleus_analysis.R",
  "09_adni_pathway_prs_analysis.R",
  "10_habs_multimodal_analysis.R",
  "11_clinical_expression.R",
  "12_clinical_figures.R"
)

for (module in modules) {
  source(file.path(root, module), encoding = "UTF-8", chdir = TRUE)
}
