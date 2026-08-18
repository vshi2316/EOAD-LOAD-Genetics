# HABS-HD multimodal analysis entry point

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)
if (is.na(script_path)) {
  stop("Run with source('10_habs_multimodal_analysis.R').", call. = FALSE)
}

root <- dirname(script_path)
modules <- c(
  "10_1_habs_multimodal_data.R",
  "10_2_habs_longitudinal_cognition.R",
  "10_3_habs_cognition_models.R",
  "10_4_habs_primary_models.R",
  "10_5_habs_regional_analysis.R",
  "10_6_habs_age65_multimodal_interactions.R"
)

for (module in modules) {
  source(file.path(root, module), encoding = "UTF-8", chdir = TRUE)
}
