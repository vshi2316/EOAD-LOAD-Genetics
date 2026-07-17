# ==============================================================================
# 07. Bradley 2025 EOAD GWAS Replication
# ==============================================================================
# Bradley/Pottier Age65 is the primary external EOAD dataset. Age70 is retained
# as a sensitivity definition. The five tested pathways are frozen from FinnGen
# EOAD discovery before evaluation in Bradley/Pottier.

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260613)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------------------------
# 1. User configuration
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv("EOAD_PROJECT_DIR", unset = getwd())
DATA_DIR <- Sys.getenv("EOAD_DATA_DIR", unset = file.path(PROJECT_DIR, "data"))
RESULTS_DIR <- Sys.getenv("EOAD_RESULTS_DIR", unset = file.path(PROJECT_DIR, "results"))

OUT_DIR <- Sys.getenv(
  "BRADLEY_REPLICATION_OUT_DIR",
  unset = file.path(RESULTS_DIR, "07_Bradley_Pottier_2025_EOAD_GWAS_replication")
)
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
MAGMA_DIR <- file.path(OUT_DIR, "magma")
INPUT_DIR <- file.path(OUT_DIR, "magma_input")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAGMA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)

BRADLEY_RAW_DIR <- Sys.getenv(
  "BRADLEY_RAW_DIR",
  unset = file.path(DATA_DIR, "Bradley_Pottier_2025")
)
BRADLEY_BUILD <- Sys.getenv("BRADLEY_BUILD", unset = "hg38")
TARGET_BUILD <- Sys.getenv("TARGET_MAGMA_BUILD", unset = "hg19")
CHAIN_FILE <- Sys.getenv(
  "HG38_TO_HG19_CHAIN",
  unset = file.path(DATA_DIR, "reference", "hg38ToHg19.over.chain.gz")
)

BRADLEY_AGE65_FILE <- Sys.getenv(
  "BRADLEY_AGE65_FILE",
  unset = file.path(BRADLEY_RAW_DIR, "NHW_Age65_sumstats_base_model_SNP-SEX-10PCs-Array.txt")
)
BRADLEY_AGE70_FILE <- Sys.getenv(
  "BRADLEY_AGE70_FILE",
  unset = file.path(BRADLEY_RAW_DIR, "NHW_Age70_sumstats_base_model_SNP-SEX-10PCs-Array.txt")
)

BRADLEY_AGE65_PRECOMPUTED_INPUT <- Sys.getenv("BRADLEY_AGE65_MAGMA_INPUT", unset = "")
BRADLEY_AGE70_PRECOMPUTED_INPUT <- Sys.getenv("BRADLEY_AGE70_MAGMA_INPUT", unset = "")

CURRENT_EOAD_GENES <- Sys.getenv(
  "FINNGEN_EOAD_MAGMA_GENES",
  unset = file.path(RESULTS_DIR, "01_Pathway_Discovery", "MAGMA_EOAD", "EOAD.genes.out.txt")
)
CURRENT_EADB_GENES <- Sys.getenv(
  "EADB_MAGMA_GENES",
  unset = file.path(RESULTS_DIR, "01_Pathway_Discovery", "MAGMA_EADB", "EADB.genes.out.txt")
)
CURRENT_EOAD_GSA <- Sys.getenv(
  "FINNGEN_EOAD_MAGMA_GSA",
  unset = file.path(RESULTS_DIR, "01_Pathway_Discovery", "MAGMA_EOAD", "EOAD.gsa.out.txt")
)
CURRENT_EADB_GSA <- Sys.getenv(
  "EADB_MAGMA_GSA",
  unset = file.path(RESULTS_DIR, "01_Pathway_Discovery", "MAGMA_EADB", "EADB.gsa.out.txt")
)
CANDIDATE_GENESET_FILE <- Sys.getenv(
  "CANDIDATE_GENESET_FILE",
  unset = file.path(DATA_DIR, "candidate_gene_sets", "eoad_eadb_candidate_gene_sets.tsv")
)

MAGMA_EXE <- Sys.getenv("MAGMA_EXE", unset = "magma")
MAGMA_BFILE <- Sys.getenv(
  "MAGMA_BFILE",
  unset = file.path(DATA_DIR, "reference", "magma", "g1000_eur", "g1000_eur")
)
MAGMA_GENE_LOC <- Sys.getenv(
  "MAGMA_GENE_LOC",
  unset = file.path(DATA_DIR, "reference", "magma", "ENSGv110.coding.genes.txt")
)
MAGMA_SET_ANNOT <- Sys.getenv(
  "MAGMA_SET_ANNOT",
  unset = file.path(DATA_DIR, "reference", "magma", "MSigDB_20231Hs_MAGMA.txt")
)
RUN_MAGMA <- tolower(Sys.getenv("RUN_MAGMA", unset = "true")) %in% c("true", "1", "yes", "y")
MIN_GENE_SET_OVERLAP <- as.integer(Sys.getenv("MIN_GENE_SET_OVERLAP", unset = "5"))

APOE_CHR <- 19L
APOE_START <- 44000000L
APOE_END <- 46000000L

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

message_step <- function(...) {
  message("\n", paste0(...))
}

safe_file <- function(path) {
  nzchar(path) && file.exists(path)
}

not_apoe <- function(x) {
  is.na(x) | !x
}

first_existing <- function(...) {
  x <- unlist(list(...), use.names = FALSE)
  x <- x[nzchar(x)]
  hit <- x[file.exists(x)]
  if (length(hit) == 0) return(x[1])
  hit[1]
}

detect_col <- function(dt, candidates, required = TRUE, label = "column") {
  hit <- intersect(candidates, names(dt))
  if (length(hit) > 0) return(hit[1])
  lower_map <- setNames(names(dt), tolower(names(dt)))
  hit_lower <- intersect(tolower(candidates), names(lower_map))
  if (length(hit_lower) > 0) return(lower_map[[hit_lower[1]]])
  if (required) {
    stop("Could not identify ", label, ". Tried: ", paste(candidates, collapse = ", "))
  }
  NA_character_
}

clean_chr <- function(x) {
  x <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
  x[x %in% c("X", "x")] <- "23"
  suppressWarnings(as.integer(x))
}

clean_allele <- function(x) {
  toupper(gsub("[^ACGT]", "", as.character(x)))
}

clean_rsid <- function(x) {
  x <- as.character(x)
  x <- ifelse(grepl("^rs[0-9]+$", x, ignore.case = TRUE), paste0("rs", sub("^rs", "", x, ignore.case = TRUE)), NA_character_)
  x
}

standardise_sumstats <- function(path, tag) {
  message_step("Reading raw summary statistics: ", tag)
  dt <- fread(path, nThread = max(1, parallel::detectCores() - 1), showProgress = TRUE)
  chr_col <- detect_col(dt, c("CHR", "CHROM", "#CHROM", "chromosome", "chrom"), label = "chromosome column")
  bp_col <- detect_col(dt, c("BP", "POS", "position", "base_pair_location", "bp"), label = "position column")
  snp_col <- detect_col(dt, c("SNP", "rsid", "rsids", "ID", "variant_id", "MarkerName"), required = FALSE)
  beta_col <- detect_col(dt, c("BETA", "beta", "Effect", "effect"), label = "beta column")
  se_col <- detect_col(dt, c("SE", "se", "standard_error", "SEBETA", "sebeta"), label = "SE column")
  p_col <- detect_col(dt, c("P", "p", "PVAL", "pval", "p_value"), required = FALSE)
  log10p_col <- detect_col(dt, c("LOG10_P", "LOG10P", "neglog10_p", "minus_log10_p"), required = FALSE)
  a1_col <- detect_col(dt, c("A1", "EA", "effect_allele", "ALT", "alt"), label = "effect allele column")
  a2_col <- detect_col(dt, c("A2", "NEA", "other_allele", "REF", "ref"), required = FALSE)
  ref_col <- detect_col(dt, c("REF", "ref"), required = FALSE)
  alt_col <- detect_col(dt, c("ALT", "alt"), required = FALSE)
  n_col <- detect_col(dt, c("N", "OBS_CT", "N_total", "N_TOTAL", "sample_size"), required = FALSE)

  out <- data.table(
    row_id = seq_len(nrow(dt)),
    SNP_raw = if (!is.na(snp_col)) as.character(dt[[snp_col]]) else NA_character_,
    SNP = if (!is.na(snp_col)) clean_rsid(dt[[snp_col]]) else NA_character_,
    CHR = clean_chr(dt[[chr_col]]),
    BP = suppressWarnings(as.integer(dt[[bp_col]])),
    A1 = clean_allele(dt[[a1_col]]),
    A2 = if (!is.na(a2_col)) clean_allele(dt[[a2_col]]) else NA_character_,
    BETA = suppressWarnings(as.numeric(dt[[beta_col]])),
    SE = suppressWarnings(as.numeric(dt[[se_col]])),
    P = if (!is.na(p_col)) suppressWarnings(as.numeric(dt[[p_col]])) else NA_real_,
    N = if (!is.na(n_col)) suppressWarnings(as.numeric(dt[[n_col]])) else NA_real_
  )

  if (all(is.na(out$P)) && !is.na(log10p_col)) {
    out[, P := 10^(-suppressWarnings(as.numeric(dt[[log10p_col]])))]
  }
  if (all(is.na(out$A2)) && !is.na(ref_col) && !is.na(alt_col)) {
    ref <- clean_allele(dt[[ref_col]])
    alt <- clean_allele(dt[[alt_col]])
    out[, A2 := ifelse(A1 == alt, ref, alt)]
  }
  if (all(is.na(out$N))) {
    out[, N := max(nrow(out), 1)]
  }

  out <- out[
    !is.na(CHR) & CHR >= 1 & CHR <= 22 &
      !is.na(BP) & !is.na(BETA) & !is.na(SE) & SE > 0 &
      !is.na(P) & P > 0 & P <= 1 &
      nchar(A1) == 1 & nchar(A2) == 1 & A1 != A2
  ]
  out[, P := pmax(P, .Machine$double.xmin)]
  out[]
}

liftover_to_hg19 <- function(dt, tag) {
  if (tolower(BRADLEY_BUILD) %in% c("hg19", "grch37") || tolower(TARGET_BUILD) %in% c("hg38", "grch38")) {
    dt[, `:=`(CHR_ORIGINAL = CHR, BP_ORIGINAL = BP)]
    return(dt)
  }
  if (!file.exists(CHAIN_FILE)) {
    stop("Chain file is required for hg38-to-hg19 liftover: ", CHAIN_FILE)
  }
  for (pkg in c("GenomicRanges", "IRanges", "rtracklayer")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package required for liftover is missing: ", pkg)
    }
  }

  message_step("Lifting coordinates to hg19: ", tag)
  chain <- rtracklayer::import.chain(CHAIN_FILE)
  dt[, `:=`(CHR_ORIGINAL = CHR, BP_ORIGINAL = BP)]
  lifted_list <- lapply(split(dt, dt$CHR), function(x) {
    gr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", x$CHR),
      ranges = IRanges::IRanges(start = x$BP, width = 1),
      row_id = x$row_id
    )
    lo <- rtracklayer::liftOver(gr, chain)
    keep <- lengths(lo) == 1
    if (!any(keep)) return(data.table(row_id = integer(), CHR_lifted = integer(), BP_lifted = integer()))
    flat <- unlist(lo[keep], use.names = FALSE)
    data.table(
      row_id = S4Vectors::mcols(gr)$row_id[keep],
      CHR_lifted = clean_chr(as.character(GenomicRanges::seqnames(flat))),
      BP_lifted = as.integer(GenomicRanges::start(flat))
    )
  })
  lifted <- rbindlist(lifted_list, use.names = TRUE, fill = TRUE)
  dt <- merge(dt, lifted, by = "row_id", all.x = FALSE, all.y = FALSE)
  dt[, `:=`(CHR = CHR_lifted, BP = BP_lifted, CHR_lifted = NULL, BP_lifted = NULL)]
  dt[!is.na(CHR) & !is.na(BP)]
}

map_chrpos_to_rsid <- function(dt, tag) {
  if (sum(!is.na(dt$SNP)) > 0.9 * nrow(dt)) return(dt)
  for (pkg in c("GenomicRanges", "IRanges", "BSgenome", "SNPlocs.Hsapiens.dbSNP155.GRCh37")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package required for chr:pos to rsID mapping is missing: ", pkg,
        "\nInstall it with BiocManager before rerunning this script."
      )
    }
  }

  message_step("Mapping hg19 chr:pos to rsID with dbSNP155 GRCh37: ", tag)
  suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP155.GRCh37))
  snplocs <- get("SNPlocs.Hsapiens.dbSNP155.GRCh37")
  snps_by_overlaps <- getExportedValue("BSgenome", "snpsByOverlaps")

  todo <- dt[is.na(SNP), .(row_id, CHR, BP)]
  mapped <- rbindlist(lapply(split(todo, todo$CHR), function(x) {
    if (nrow(x) == 0) return(data.table(row_id = integer(), SNP_mapped = character()))
    try_prefix <- function(prefix) {
      gr <- GenomicRanges::GRanges(
        seqnames = paste0(prefix, x$CHR),
        ranges = IRanges::IRanges(start = x$BP, width = 1),
        row_id = x$row_id
      )
      hits <- try(snps_by_overlaps(snplocs, gr, columns = c("RefSNP_id")), silent = TRUE)
      if (inherits(hits, "try-error") || length(hits) == 0) return(NULL)
      h <- as.data.table(hits)
      if (!"RefSNP_id" %in% names(h)) {
        rid <- grep("RefSNP|rsid|snp", names(h), ignore.case = TRUE, value = TRUE)
        if (length(rid) == 0) return(NULL)
        setnames(h, rid[1], "RefSNP_id")
      }
      h[, query_pos := as.integer(start)]
      merge(
        x,
        h[, .(CHR = clean_chr(as.character(seqnames)), BP = query_pos, RefSNP_id)],
        by = c("CHR", "BP"),
        allow.cartesian = TRUE
      )[, .(SNP_mapped = paste0("rs", sub("^rs", "", as.character(RefSNP_id), ignore.case = TRUE))[1]), by = row_id]
    }
    out <- try_prefix("")
    if (is.null(out) || nrow(out) == 0) out <- try_prefix("chr")
    if (is.null(out)) data.table(row_id = integer(), SNP_mapped = character()) else out
  }), use.names = TRUE, fill = TRUE)

  dt <- merge(dt, mapped, by = "row_id", all.x = TRUE)
  dt[is.na(SNP) & !is.na(SNP_mapped), SNP := SNP_mapped]
  dt[, SNP_mapped := NULL]
  dt[!is.na(SNP)]
}

write_standard_magma_input <- function(dt, tag) {
  out <- copy(dt)
  out <- out[!duplicated(SNP)]
  out <- out[order(P)]
  out_file <- file.path(INPUT_DIR, paste0(tag, "_MAGMA_input_rsID_hg19.tsv"))
  fwrite(out[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, N)], out_file, sep = "\t", quote = FALSE, na = "NA")
  qc <- data.table(
    dataset = tag,
    n_variants = nrow(out),
    n_unique_rsids = uniqueN(out$SNP),
    min_p = min(out$P, na.rm = TRUE),
    n_p_lt_5e_8 = sum(out$P < 5e-8, na.rm = TRUE),
    n_p_lt_1e_5 = sum(out$P < 1e-5, na.rm = TRUE),
    mean_n = mean(out$N, na.rm = TRUE),
    input_file = out_file
  )
  fwrite(qc, file.path(TAB_DIR, paste0(tag, "_MAGMA_input_QC.tsv")), sep = "\t")
  list(data = out, file = out_file, qc = qc)
}

prepare_bradley_input <- function(raw_file, precomputed_file, tag) {
  if (safe_file(precomputed_file)) {
    message_step("Using precomputed MAGMA input: ", tag)
    dt <- fread(precomputed_file, nThread = max(1, parallel::detectCores() - 1))
    return(write_standard_magma_input(dt, tag))
  }
  if (!file.exists(raw_file)) {
    stop("Raw Bradley summary-statistics file not found for ", tag, ": ", raw_file)
  }
  dt <- standardise_sumstats(raw_file, tag)
  qc0 <- data.table(dataset = tag, n_after_basic_qc = nrow(dt), n_rsid_before_mapping = sum(!is.na(dt$SNP)))
  dt <- liftover_to_hg19(dt, tag)
  dt <- map_chrpos_to_rsid(dt, tag)
  qc1 <- data.table(dataset = tag, n_after_liftover_rsid_mapping = nrow(dt), n_rsid_after_mapping = sum(!is.na(dt$SNP)))
  fwrite(cbind(qc0, qc1[, -1]), file.path(TAB_DIR, paste0(tag, "_coordinate_rsid_QC.tsv")), sep = "\t")
  write_standard_magma_input(dt, tag)
}

run_magma_for_input <- function(input_file, tag) {
  out_prefix <- file.path(MAGMA_DIR, tag, tag)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

  genes_out <- first_existing(paste0(out_prefix, ".genes.out"), paste0(out_prefix, ".genes.out.txt"))
  genes_raw <- first_existing(paste0(out_prefix, ".genes.raw"))
  gsa_out <- first_existing(paste0(out_prefix, ".gsa.out"), paste0(out_prefix, ".gsa.out.txt"))

  if (!RUN_MAGMA) {
    return(data.table(dataset = tag, genes_out = genes_out, genes_raw = genes_raw, gsa_out = gsa_out, status = "not_run"))
  }
  if (!file.exists(MAGMA_GENE_LOC)) stop("MAGMA gene-location file not found: ", MAGMA_GENE_LOC)
  if (!file.exists(paste0(MAGMA_BFILE, ".bed"))) stop("MAGMA bfile prefix not found: ", MAGMA_BFILE)

  dt <- fread(input_file)
  snp_loc <- file.path(INPUT_DIR, paste0(tag, "_MAGMA_snp_loc.tsv"))
  pval_file <- file.path(INPUT_DIR, paste0(tag, "_MAGMA_pval.tsv"))
  fwrite(dt[, .(SNP, CHR, BP)], snp_loc, sep = "\t", col.names = FALSE, quote = FALSE)
  fwrite(dt[, .(SNP, P, N)], pval_file, sep = "\t", quote = FALSE)

  if (!file.exists(genes_raw) || !file.exists(genes_out)) {
    message_step("Running MAGMA gene analysis: ", tag)
    system2(MAGMA_EXE, c("--annotate", "--snp-loc", snp_loc, "--gene-loc", MAGMA_GENE_LOC, "--out", out_prefix))
    annot_file <- first_existing(paste0(out_prefix, ".genes.annot"), paste0(out_prefix, ".genes.annot.txt"))
    system2(MAGMA_EXE, c(
      "--bfile", MAGMA_BFILE,
      "--pval", pval_file, "use=SNP,P", "ncol=N",
      "--gene-annot", annot_file,
      "--out", out_prefix
    ))
  }
  if (file.exists(MAGMA_SET_ANNOT) && !file.exists(gsa_out)) {
    message_step("Running MAGMA gene-set analysis: ", tag)
    genes_raw <- first_existing(paste0(out_prefix, ".genes.raw"))
    system2(MAGMA_EXE, c("--gene-results", genes_raw, "--set-annot", MAGMA_SET_ANNOT, "--out", paste0(out_prefix, ".gsa")))
    gsa_out <- first_existing(paste0(out_prefix, ".gsa.gsa.out"), paste0(out_prefix, ".gsa.out"), paste0(out_prefix, ".gsa.out.txt"))
  }
  data.table(dataset = tag, genes_out = first_existing(genes_out), genes_raw = first_existing(genes_raw), gsa_out = first_existing(gsa_out), status = "completed")
}

read_magma_genes <- function(path, dataset) {
  if (!file.exists(path)) stop("MAGMA gene output not found: ", path)
  dt <- fread(path)
  gene_col <- detect_col(dt, c("GENE", "gene", "GENE_ID", "gene_id"), label = "MAGMA gene column")
  p_col <- detect_col(dt, c("P", "p", "PVALUE", "p_value"), label = "MAGMA P column")
  z_col <- detect_col(dt, c("ZSTAT", "Z", "z"), required = FALSE)
  chr_col <- detect_col(dt, c("CHR", "chromosome"), required = FALSE)
  start_col <- detect_col(dt, c("START", "start"), required = FALSE)
  stop_col <- detect_col(dt, c("STOP", "END", "stop", "end"), required = FALSE)
  symbol_col <- detect_col(dt, c("SYMBOL", "symbol", "GENE_NAME", "gene_name"), required = FALSE)

  out <- data.table(
    dataset = dataset,
    GENE = as.character(dt[[gene_col]]),
    SYMBOL = if (!is.na(symbol_col)) toupper(as.character(dt[[symbol_col]])) else toupper(as.character(dt[[gene_col]])),
    P = as.numeric(dt[[p_col]]),
    ZSTAT = if (!is.na(z_col)) as.numeric(dt[[z_col]]) else qnorm(pmax(.Machine$double.xmin, 1 - as.numeric(dt[[p_col]]) / 2)),
    CHR = if (!is.na(chr_col)) clean_chr(dt[[chr_col]]) else NA_integer_,
    START = if (!is.na(start_col)) as.integer(dt[[start_col]]) else NA_integer_,
    STOP = if (!is.na(stop_col)) as.integer(dt[[stop_col]]) else NA_integer_
  )
  out[, FDR := p.adjust(P, "BH")]
  out[, APOE_REGION := CHR == APOE_CHR & !is.na(START) & !is.na(STOP) & START <= APOE_END & STOP >= APOE_START]
  out[]
}

read_magma_gsa <- function(path, dataset) {
  if (!file.exists(path)) return(NULL)
  dt <- tryCatch(fread(path, fill = TRUE), error = function(e) fread(path, skip = "VARIABLE", fill = TRUE))
  setnames(dt, names(dt), gsub("^#", "", names(dt)))
  set_col <- detect_col(dt, c("VARIABLE", "FULL_NAME", "GENESET", "SET", "NAME"), label = "MAGMA gene-set column")
  p_col <- detect_col(dt, c("P", "P-value", "PVALUE"), label = "MAGMA gene-set P column")
  beta_col <- detect_col(dt, c("BETA", "BETA_STD", "ZSTAT", "Z"), required = FALSE)
  data.table(
    dataset = dataset,
    gene_set = as.character(dt[[set_col]]),
    statistic = if (!is.na(beta_col)) suppressWarnings(as.numeric(dt[[beta_col]])) else NA_real_,
    P = suppressWarnings(as.numeric(dt[[p_col]]))
  )[
    !is.na(P)
  ][
    , FDR := p.adjust(P, "BH")
  ][]
}

read_candidate_gene_sets <- function(path) {
  if (!file.exists(path)) {
    warning("Candidate gene-set file not found; candidate pathway tests will be skipped: ", path)
    return(NULL)
  }
  gs <- fread(path)
  set_col <- detect_col(gs, c("set_name", "pathway", "gene_set", "module"), label = "candidate set column")
  gene_col <- detect_col(gs, c("GENE", "gene", "gene_id", "magma_gene_id", "symbol", "gene_symbol"), label = "candidate gene column")
  out <- data.table(
    set_name = as.character(gs[[set_col]]),
    gene_token = toupper(as.character(gs[[gene_col]]))
  )
  unique(out[!is.na(set_name) & !is.na(gene_token) & nzchar(gene_token)])
}

candidate_rank_test <- function(magma_dt, candidate_sets, dataset, analysis) {
  if (is.null(candidate_sets)) return(NULL)
  dt <- copy(magma_dt)
  if (analysis == "APOE_excluded") dt <- dt[not_apoe(APOE_REGION)]
  rbindlist(lapply(split(candidate_sets, candidate_sets$set_name), function(gs) {
    in_set <- dt[GENE %in% gs$gene_token | SYMBOL %in% gs$gene_token]
    bg <- dt[!(GENE %in% gs$gene_token | SYMBOL %in% gs$gene_token)]
    if (nrow(in_set) < MIN_GENE_SET_OVERLAP || nrow(bg) < 100) {
      return(data.table(
        dataset = dataset, analysis = analysis, set_name = gs$set_name[1],
        n_genes = nrow(in_set), mean_z = NA_real_, median_z = NA_real_,
        wilcox_p = NA_real_, status = "insufficient_overlap"
      ))
    }
    wt <- wilcox.test(in_set$ZSTAT, bg$ZSTAT, alternative = "greater")
    data.table(
      dataset = dataset,
      analysis = analysis,
      set_name = gs$set_name[1],
      n_genes = nrow(in_set),
      mean_z = mean(in_set$ZSTAT, na.rm = TRUE),
      median_z = median(in_set$ZSTAT, na.rm = TRUE),
      wilcox_p = wt$p.value,
      status = "tested"
    )
  }), use.names = TRUE, fill = TRUE)[
    , wilcox_fdr := p.adjust(wilcox_p, "BH"), by = .(dataset, analysis)
  ][]
}

gene_concordance <- function(a, b, label_a, label_b, analysis) {
  x <- merge(a, b, by = "GENE", suffixes = c("_a", "_b"))
  if (analysis == "APOE_excluded") x <- x[not_apoe(APOE_REGION_a) & not_apoe(APOE_REGION_b)]
  data.table(
    dataset_a = label_a,
    dataset_b = label_b,
    analysis = analysis,
    n_genes = nrow(x),
    pearson_r = cor(x$ZSTAT_a, x$ZSTAT_b, method = "pearson", use = "pairwise.complete.obs"),
    pearson_p = cor.test(x$ZSTAT_a, x$ZSTAT_b, method = "pearson")$p.value,
    spearman_r = cor(x$ZSTAT_a, x$ZSTAT_b, method = "spearman", use = "pairwise.complete.obs"),
    spearman_p = suppressWarnings(cor.test(x$ZSTAT_a, x$ZSTAT_b, method = "spearman")$p.value)
  )
}

gsa_concordance <- function(a, b, label_a, label_b) {
  if (is.null(a) || is.null(b)) return(NULL)
  x <- merge(a, b, by = "gene_set", suffixes = c("_a", "_b"))
  data.table(
    dataset_a = label_a,
    dataset_b = label_b,
    n_gene_sets = nrow(x),
    pearson_r = cor(-log10(x$P_a), -log10(x$P_b), method = "pearson", use = "pairwise.complete.obs"),
    pearson_p = cor.test(-log10(x$P_a), -log10(x$P_b), method = "pearson")$p.value,
    spearman_r = cor(-log10(x$P_a), -log10(x$P_b), method = "spearman", use = "pairwise.complete.obs"),
    spearman_p = suppressWarnings(cor.test(-log10(x$P_a), -log10(x$P_b), method = "spearman")$p.value)
  )
}

plot_gene_concordance <- function(a, b, label_b, analysis, out_name) {
  x <- merge(a, b, by = "GENE", suffixes = c("_FinnGen", "_Bradley"))
  if (analysis == "APOE_excluded") x <- x[not_apoe(APOE_REGION_FinnGen) & not_apoe(APOE_REGION_Bradley)]
  x[, highlight := fifelse(P_FinnGen < 0.05 / nrow(a) | P_Bradley < 0.05 / nrow(b), "Bonferroni gene", "Other genes")]
  p <- ggplot(x, aes(ZSTAT_FinnGen, ZSTAT_Bradley)) +
    geom_hline(yintercept = 0, color = "grey85") +
    geom_vline(xintercept = 0, color = "grey85") +
    geom_point(aes(color = highlight), alpha = 0.55, size = 0.9) +
    scale_color_manual(values = c("Bonferroni gene" = "#b24745", "Other genes" = "#496f8a")) +
    labs(
      title = paste0("FinnGen EOAD vs ", label_b, " MAGMA gene-level concordance"),
      subtitle = analysis,
      x = "FinnGen EOAD MAGMA Z",
      y = paste0(label_b, " MAGMA Z"),
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(file.path(FIG_DIR, paste0(out_name, ".pdf")), p, width = 6.2, height = 5.2)
  ggsave(file.path(FIG_DIR, paste0(out_name, ".png")), p, width = 6.2, height = 5.2, dpi = 300)
}

plot_candidate_heatmap <- function(candidate_results) {
  x <- candidate_results[status == "tested" & analysis == "APOE_excluded"]
  if (nrow(x) == 0) return(invisible(NULL))
  x[, label := ifelse(!is.na(wilcox_fdr) & wilcox_fdr < 0.05, "*", "")]
  p <- ggplot(x, aes(dataset, set_name, fill = mean_z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4.2) +
    scale_fill_gradient2(low = "#2c5c88", mid = "white", high = "#b24745", midpoint = 0, name = "Mean MAGMA Z") +
    labs(
      title = "Candidate pathway replication in independent EOAD GWAS",
      subtitle = "APOE-excluded; asterisk indicates FDR < 0.05 in rank-based gene-set test",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(FIG_DIR, "Figure_07_Bradley_candidate_pathway_replication_APOE_excluded.pdf"), p, width = 8.5, height = 5.6)
  ggsave(file.path(FIG_DIR, "Figure_07_Bradley_candidate_pathway_replication_APOE_excluded.png"), p, width = 8.5, height = 5.6, dpi = 300)
}

# ------------------------------------------------------------------------------
# 3. Prepare Bradley inputs and run MAGMA
# ------------------------------------------------------------------------------

message_step("07. Bradley/Pottier 2025 EOAD GWAS replication")

BRADLEY_AGE65_FILE <- first_existing(
  BRADLEY_AGE65_FILE,
  file.path(BRADLEY_RAW_DIR, "NHW_Age65_sumstats_base_model_SNP-SEX-10PCs-Array (1).txt")
)
BRADLEY_AGE70_FILE <- first_existing(BRADLEY_AGE70_FILE)

age65_input <- prepare_bradley_input(BRADLEY_AGE65_FILE, BRADLEY_AGE65_PRECOMPUTED_INPUT, "Bradley_NHW_Age65")
age70_input <- prepare_bradley_input(BRADLEY_AGE70_FILE, BRADLEY_AGE70_PRECOMPUTED_INPUT, "Bradley_NHW_Age70")
fwrite(rbind(age65_input$qc, age70_input$qc, fill = TRUE), file.path(TAB_DIR, "Bradley_MAGMA_input_QC_all.tsv"), sep = "\t")

magma_outputs <- rbindlist(list(
  run_magma_for_input(age65_input$file, "Bradley_NHW_Age65"),
  run_magma_for_input(age70_input$file, "Bradley_NHW_Age70")
), fill = TRUE)
fwrite(magma_outputs, file.path(TAB_DIR, "Bradley_MAGMA_output_files.tsv"), sep = "\t")

# ------------------------------------------------------------------------------
# 4. Read MAGMA outputs and summarize replication
# ------------------------------------------------------------------------------

message_step("Summarising MAGMA gene-level, gene-set, and candidate pathway replication")

bradley65_genes <- read_magma_genes(magma_outputs[dataset == "Bradley_NHW_Age65", genes_out][1], "Bradley_NHW_Age65")
bradley70_genes <- read_magma_genes(magma_outputs[dataset == "Bradley_NHW_Age70", genes_out][1], "Bradley_NHW_Age70")
current_eoad_genes <- read_magma_genes(CURRENT_EOAD_GENES, "FinnGen_EOAD")
current_eadb_genes <- read_magma_genes(CURRENT_EADB_GENES, "EADB")

gene_overview <- rbindlist(lapply(list(current_eoad_genes, current_eadb_genes, bradley65_genes, bradley70_genes), function(x) {
  data.table(
    dataset = x$dataset[1],
    n_genes = nrow(x),
    n_fdr_005 = sum(x$FDR < 0.05, na.rm = TRUE),
    n_bonferroni_005 = sum(x$P < 0.05 / nrow(x), na.rm = TRUE),
    min_p = min(x$P, na.rm = TRUE)
  )
}), fill = TRUE)
fwrite(gene_overview, file.path(TAB_DIR, "Table_07_MAGMA_gene_level_overview.tsv"), sep = "\t")

gene_conc <- rbindlist(list(
  gene_concordance(current_eoad_genes, bradley65_genes, "FinnGen_EOAD", "Bradley_NHW_Age65", "APOE_included"),
  gene_concordance(current_eoad_genes, bradley65_genes, "FinnGen_EOAD", "Bradley_NHW_Age65", "APOE_excluded"),
  gene_concordance(current_eoad_genes, bradley70_genes, "FinnGen_EOAD", "Bradley_NHW_Age70", "APOE_included"),
  gene_concordance(current_eoad_genes, bradley70_genes, "FinnGen_EOAD", "Bradley_NHW_Age70", "APOE_excluded")
), fill = TRUE)
fwrite(gene_conc, file.path(TAB_DIR, "Table_07_Bradley_vs_FinnGen_EOAD_gene_level_concordance.tsv"), sep = "\t")

bradley65_gsa <- read_magma_gsa(magma_outputs[dataset == "Bradley_NHW_Age65", gsa_out][1], "Bradley_NHW_Age65")
bradley70_gsa <- read_magma_gsa(magma_outputs[dataset == "Bradley_NHW_Age70", gsa_out][1], "Bradley_NHW_Age70")
current_eoad_gsa <- read_magma_gsa(CURRENT_EOAD_GSA, "FinnGen_EOAD")
current_eadb_gsa <- read_magma_gsa(CURRENT_EADB_GSA, "EADB")
gsa_conc <- rbindlist(list(
  gsa_concordance(current_eoad_gsa, bradley65_gsa, "FinnGen_EOAD", "Bradley_NHW_Age65"),
  gsa_concordance(current_eoad_gsa, bradley70_gsa, "FinnGen_EOAD", "Bradley_NHW_Age70")
), fill = TRUE)
if (nrow(gsa_conc) > 0) {
  fwrite(gsa_conc, file.path(TAB_DIR, "Table_07_Bradley_vs_FinnGen_EOAD_MAGMA_gene_set_concordance.tsv"), sep = "\t")
}

candidate_sets <- read_candidate_gene_sets(CANDIDATE_GENESET_FILE)
candidate_results <- rbindlist(list(
  candidate_rank_test(current_eoad_genes, candidate_sets, "FinnGen_EOAD", "APOE_included"),
  candidate_rank_test(current_eoad_genes, candidate_sets, "FinnGen_EOAD", "APOE_excluded"),
  candidate_rank_test(bradley65_genes, candidate_sets, "Bradley_NHW_Age65", "APOE_included"),
  candidate_rank_test(bradley65_genes, candidate_sets, "Bradley_NHW_Age65", "APOE_excluded"),
  candidate_rank_test(bradley70_genes, candidate_sets, "Bradley_NHW_Age70", "APOE_included"),
  candidate_rank_test(bradley70_genes, candidate_sets, "Bradley_NHW_Age70", "APOE_excluded"),
  candidate_rank_test(current_eadb_genes, candidate_sets, "EADB", "APOE_included"),
  candidate_rank_test(current_eadb_genes, candidate_sets, "EADB", "APOE_excluded")
), fill = TRUE)
if (nrow(candidate_results) > 0) {
  fwrite(candidate_results, file.path(TAB_DIR, "Table_07_candidate_pathway_replication_by_MAGMA_rank.tsv"), sep = "\t")
  candidate_wide <- dcast(
    candidate_results[status == "tested"],
    set_name + analysis ~ dataset,
    value.var = c("mean_z", "wilcox_p", "wilcox_fdr")
  )
  fwrite(candidate_wide, file.path(TAB_DIR, "Table_07_candidate_pathway_replication_wide.tsv"), sep = "\t")
}

# ------------------------------------------------------------------------------
# 5. Figures
# ------------------------------------------------------------------------------

plot_gene_concordance(current_eoad_genes, bradley65_genes, "Bradley Age65", "APOE_excluded", "Figure_07_Bradley_Age65_gene_level_concordance_APOE_excluded")
plot_gene_concordance(current_eoad_genes, bradley70_genes, "Bradley Age70", "APOE_excluded", "Figure_07_Bradley_Age70_gene_level_concordance_APOE_excluded")
if (nrow(candidate_results) > 0) plot_candidate_heatmap(candidate_results)

# ------------------------------------------------------------------------------
# 6. Run summary
# ------------------------------------------------------------------------------

summary_lines <- c(
  "07. Bradley/Pottier 2025 EOAD GWAS replication completed.",
  paste0("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE)),
  "",
  "Primary manuscript-facing outputs:",
  file.path(TAB_DIR, "Table_07_Bradley_vs_FinnGen_EOAD_gene_level_concordance.tsv"),
  file.path(TAB_DIR, "Table_07_Bradley_vs_FinnGen_EOAD_MAGMA_gene_set_concordance.tsv"),
  file.path(TAB_DIR, "Table_07_candidate_pathway_replication_by_MAGMA_rank.tsv"),
  file.path(FIG_DIR, "Figure_07_Bradley_Age65_gene_level_concordance_APOE_excluded.png"),
  file.path(FIG_DIR, "Figure_07_Bradley_candidate_pathway_replication_APOE_excluded.png"),
  "",
  "Interpretation note:",
  "Bradley Age65 is used as primary independent EOAD genetic replication.",
  "Bradley Age70 is used as sensitivity analysis."
)
writeLines(summary_lines, file.path(OUT_DIR, "README_07_Bradley_EOAD_GWAS_replication.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
