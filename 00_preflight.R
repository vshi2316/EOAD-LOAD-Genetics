required_packages <- c(
  "data.table", "dplyr", "tidyr", "readr", "stringr", "tibble", "purrr",
  "ggplot2", "ggrepel", "scales", "forcats", "broom", "survival",
  "lmtest", "sandwich", "pROC", "cowplot", "RColorBrewer", "bigsnpr",
  "bigstatsr", "nlme", "patchwork", "svglite", "ragg", "GenomicRanges", "IRanges", "rtracklayer", "BSgenome",
  "SNPlocs.Hsapiens.dbSNP155.GRCh37", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "ComplexHeatmap", "edgeR"
)

required_environment <- c(
  "EOAD_PROJECT_DIR", "EOAD_DATA_DIR", "EOAD_RESULTS_DIR",
  "FINNGEN_EOAD_SUMSTATS", "EADB_ADRD_SUMSTATS", "BRADLEY_RAW_DIR",
  "MAGMA_EXE", "MAGMA_BFILE", "MAGMA_GENE_LOC", "MAGMA_SET_ANNOT"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

environment_values <- Sys.getenv(required_environment, unset = "")
missing_environment <- names(environment_values)[!nzchar(environment_values)]

path_like_variables <- setdiff(required_environment, "MAGMA_BFILE")
missing_paths <- path_like_variables[
  nzchar(environment_values[path_like_variables]) &
    !file.exists(environment_values[path_like_variables]) &
    !dir.exists(environment_values[path_like_variables])
]

magma_prefix <- environment_values[["MAGMA_BFILE"]]
missing_magma_prefix <- nzchar(magma_prefix) &&
  !all(file.exists(paste0(magma_prefix, c(".bed", ".bim", ".fam"))))

cat("EOAD-LOAD Genetics preflight\n")
cat("=============================\n")
cat("R version:", R.version.string, "\n")
cat("Missing packages:",
    if (length(missing_packages)) paste(missing_packages, collapse = ", ") else "none",
    "\n")
cat("Unset environment variables:",
    if (length(missing_environment)) paste(missing_environment, collapse = ", ") else "none",
    "\n")
cat("Configured paths not found:",
    if (length(missing_paths)) paste(missing_paths, collapse = ", ") else "none",
    "\n")
cat("MAGMA reference prefix complete:", !missing_magma_prefix, "\n")

if (length(missing_packages) || length(missing_environment) ||
    length(missing_paths) || missing_magma_prefix) {
  stop("Preflight failed. Resolve the items above before running analyses.")
}

cat("Preflight passed.\n")
