# ADNI pathway PRS analysis entry point

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) stop("Run with source('09_adni_pathway_prs_analysis.R').", call. = FALSE)
root <- dirname(script_path)
source(file.path(root, "09_1_adni_pathway_prs_data.R"), encoding = "UTF-8", chdir = TRUE)
source(file.path(root, "09_2_adni_pathway_prs_centiloid.R"), encoding = "UTF-8", chdir = TRUE)
