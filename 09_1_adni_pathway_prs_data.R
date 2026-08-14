# Construct ADNI pathway polygenic scores and merge continuous Centiloid

options(stringsAsFactors = FALSE, scipen = 999)
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) stop("Run with source('09_1_adni_pathway_prs_data.R').", call. = FALSE)
ROOT <- dirname(script_path)
DATA <- Sys.getenv("EOAD_DATA_DIR", unset = file.path(ROOT, "data"))
RESULTS <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(ROOT, "results"))
OUT <- file.path(RESULTS, "adni_pathway_prs_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table.", call. = FALSE)
library(data.table)

PLINK <- Sys.getenv("PLINK_EXE", unset = "plink")
BFILE <- Sys.getenv("ADNI_PLINK_PREFIX", unset = file.path(DATA, "ADNI", "genotype", "adni_qc"))
WEIGHTS <- Sys.getenv("ADNI_PATHWAY_WEIGHT_FILE", unset = file.path(DATA, "ADNI", "pathway_score_weights.tsv"))
PHENOTYPE <- Sys.getenv("ADNI_CENTILOID_COVARIATE_FILE", unset = file.path(DATA, "ADNI", "ADNI_Centiloid_covariates.tsv.gz"))
OUTPUT <- Sys.getenv("ADNI_PATHWAY_ANALYSIS_FILE", unset = file.path(DATA, "ADNI", "ADNI_pathway_PRS_Centiloid.tsv.gz"))

required <- c(paste0(BFILE, c(".bed", ".bim", ".fam")), WEIGHTS, PHENOTYPE)
if (any(!file.exists(required))) stop("ADNI genotype, weight, or phenotype inputs are incomplete.", call. = FALSE)

weights <- fread(WEIGHTS)
needed <- c("pathway", "variant_id", "effect_allele", "weight")
if (!all(needed %in% names(weights))) stop("Pathway weight fields are incomplete.", call. = FALSE)
pathways <- c("AmyloidBetaClearance", "NegativeRegulationAPPCatabolism")
if (!all(pathways %in% weights$pathway)) stop("Both replicated pathway weights are required.", call. = FALSE)

read_profile <- function(path) {
  x <- fread(path)
  id <- intersect(c("IID", "#IID"), names(x))[1]
  score <- intersect(c("SCORE", "SCORESUM", "SCORE1_SUM"), names(x))[1]
  if (is.na(id) || is.na(score)) stop("PLINK score fields were not recognized.", call. = FALSE)
  x[, .(IID = as.character(get(id)), score = as.numeric(get(score)))]
}

scores <- list()
for (pathway_name in pathways) {
  w <- weights[pathway == pathway_name, .(variant_id, effect_allele, weight)]
  score_file <- file.path(OUT, paste0("weights_", pathway_name, ".txt"))
  fwrite(w, score_file, sep = "\t", col.names = FALSE, quote = FALSE)
  prefix <- file.path(OUT, paste0("plink_", pathway_name))
  status <- system2(PLINK, c("--bfile", BFILE, "--score", score_file, "1", "2", "3", "sum", "--out", prefix))
  if (!identical(status, 0L)) stop("PLINK scoring failed for ", pathway_name, call. = FALSE)
  profile <- paste0(prefix, ".profile")
  if (!file.exists(profile)) stop("PLINK profile was not produced for ", pathway_name, call. = FALSE)
  x <- read_profile(profile)
  setnames(x, "score", pathway_name)
  scores[[pathway_name]] <- x
}

prs <- Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), scores)
phenotype <- fread(PHENOTYPE)
if (!"IID" %in% names(phenotype)) {
  candidate <- intersect(c("RID", "subject_id"), names(phenotype))[1]
  if (is.na(candidate)) stop("The ADNI phenotype identifier was not found.", call. = FALSE)
  setnames(phenotype, candidate, "IID")
}
phenotype[, IID := as.character(IID)]
analysis <- merge(phenotype, prs, by = "IID", all.x = TRUE)
for (pathway_name in pathways) {
  sigma <- sd(analysis[[pathway_name]], na.rm = TRUE)
  analysis[, paste0("PRS_", pathway_name, "_z") := (get(pathway_name) - mean(get(pathway_name), na.rm = TRUE)) / sigma]
}

dir.create(dirname(OUTPUT), recursive = TRUE, showWarnings = FALSE)
fwrite(analysis, OUTPUT, sep = "\t", quote = FALSE, na = "NA")
counts <- data.table(pathway = pathways, variants = vapply(pathways, function(x) weights[pathway == x, uniqueN(variant_id)], integer(1)))
fwrite(counts, file.path(OUT, "pathway_score_variant_counts.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
