# ==============================================================================
# 08. Bradley 2025 SNP_Locus-Level Replication
# ==============================================================================
# Directional tests use LD-independent non-APOE loci; the APOE LD block is
# represented by one index variant rather than multiple correlated SNPs.
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
  "SNP_LOCUS_REPLICATION_OUT_DIR",
  unset = file.path(RESULTS_DIR, "08_Bradley_Pottier_2025_SNP_locus_replication")
)
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

FINNGEN_EOAD_SUMSTATS <- Sys.getenv(
  "FINNGEN_EOAD_SUMSTATS",
  unset = file.path(DATA_DIR, "GWAS", "finngen_R11_AD_EO_EXMORE.gz")
)
BRADLEY_AGE65_INPUT <- Sys.getenv(
  "BRADLEY_AGE65_MAGMA_INPUT",
  unset = file.path(RESULTS_DIR, "07_Bradley_Pottier_2025_EOAD_GWAS_replication", "magma_input", "Bradley_NHW_Age65_MAGMA_input_rsID_hg19.tsv")
)
BRADLEY_AGE70_INPUT <- Sys.getenv(
  "BRADLEY_AGE70_MAGMA_INPUT",
  unset = file.path(RESULTS_DIR, "07_Bradley_Pottier_2025_EOAD_GWAS_replication", "magma_input", "Bradley_NHW_Age70_MAGMA_input_rsID_hg19.tsv")
)
CANDIDATE_GENESET_FILE <- Sys.getenv(
  "CANDIDATE_GENESET_FILE",
  unset = file.path(DATA_DIR, "candidate_gene_sets", "eoad_eadb_candidate_gene_sets.tsv")
)
GENE_LOC_FILE <- Sys.getenv(
  "MAGMA_GENE_LOC",
  unset = file.path(DATA_DIR, "reference", "magma", "ENSGv110.coding.genes.txt")
)

TOP_N_FOR_HARMONIZED_EXPORT <- as.integer(Sys.getenv("TOP_N_FOR_HARMONIZED_EXPORT", unset = "200000"))
LEAD_P_THRESHOLD <- as.numeric(Sys.getenv("LEAD_P_THRESHOLD", unset = "1e-5"))
LEAD_WINDOW_BP <- as.integer(Sys.getenv("LEAD_WINDOW_BP", unset = "500000"))
LOCUS_SUPPORT_WINDOW_BP <- as.integer(Sys.getenv("LOCUS_SUPPORT_WINDOW_BP", unset = "500000"))
GENE_WINDOW_BP <- as.integer(Sys.getenv("GENE_WINDOW_BP", unset = "100000"))

APOE_CHR <- 19L
APOE_START <- 44000000L
APOE_END <- 46000000L

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

message_step <- function(...) {
  message("\n", paste0(...))
}

not_apoe <- function(x) {
  is.na(x) | !x
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

is_palindromic <- function(a1, a2) {
  pair <- paste0(a1, a2)
  pair %in% c("AT", "TA", "CG", "GC")
}

standardise_for_harmonization <- function(path, dataset) {
  message_step("Reading summary statistics for harmonization: ", dataset)
  dt <- fread(path, nThread = max(1, parallel::detectCores() - 1), showProgress = TRUE)
  snp_col <- detect_col(dt, c("SNP", "rsid", "rsids", "ID", "variant_id", "MarkerName"), label = "SNP column")
  chr_col <- detect_col(dt, c("CHR", "CHROM", "#CHROM", "chromosome", "chrom"), required = FALSE)
  bp_col <- detect_col(dt, c("BP", "POS", "position", "base_pair_location", "bp"), required = FALSE)
  a1_col <- detect_col(dt, c("A1", "EA", "effect_allele", "ALT", "alt"), label = "effect allele column")
  a2_col <- detect_col(dt, c("A2", "NEA", "other_allele", "REF", "ref"), label = "other allele column")
  beta_col <- detect_col(dt, c("BETA", "beta", "Effect", "effect"), label = "beta column")
  se_col <- detect_col(dt, c("SE", "se", "standard_error", "SEBETA", "sebeta"), label = "SE column")
  p_col <- detect_col(dt, c("P", "p", "PVAL", "pval", "p_value"), required = FALSE)
  log10p_col <- detect_col(dt, c("LOG10_P", "LOG10P", "neglog10_p", "minus_log10_p"), required = FALSE)
  n_col <- detect_col(dt, c("N", "OBS_CT", "N_total", "N_TOTAL", "sample_size"), required = FALSE)

  out <- data.table(
    dataset = dataset,
    SNP = clean_rsid(dt[[snp_col]]),
    CHR = if (!is.na(chr_col)) clean_chr(dt[[chr_col]]) else NA_integer_,
    BP = if (!is.na(bp_col)) suppressWarnings(as.integer(dt[[bp_col]])) else NA_integer_,
    A1 = clean_allele(dt[[a1_col]]),
    A2 = clean_allele(dt[[a2_col]]),
    BETA = suppressWarnings(as.numeric(dt[[beta_col]])),
    SE = suppressWarnings(as.numeric(dt[[se_col]])),
    P = if (!is.na(p_col)) suppressWarnings(as.numeric(dt[[p_col]])) else NA_real_,
    N = if (!is.na(n_col)) suppressWarnings(as.numeric(dt[[n_col]])) else NA_real_
  )
  if (all(is.na(out$P)) && !is.na(log10p_col)) {
    out[, P := 10^(-suppressWarnings(as.numeric(dt[[log10p_col]])))]
  }
  out <- out[
    !is.na(SNP) &
      !is.na(BETA) & !is.na(SE) & SE > 0 &
      !is.na(P) & P > 0 & P <= 1 &
      nchar(A1) == 1 & nchar(A2) == 1 & A1 != A2
  ]
  out[, P := pmax(P, .Machine$double.xmin)]
  out <- out[order(P)]
  out <- out[!duplicated(SNP)]
  out[]
}

harmonize_pair <- function(finngen, bradley, label) {
  message_step("Harmonizing FinnGen EOAD with ", label)
  x <- merge(
    finngen,
    bradley,
    by = "SNP",
    suffixes = c("_FINNGEN", "_BRADLEY")
  )
  x[, allele_match := A1_FINNGEN == A1_BRADLEY & A2_FINNGEN == A2_BRADLEY]
  x[, allele_flip := A1_FINNGEN == A2_BRADLEY & A2_FINNGEN == A1_BRADLEY]
  x <- x[allele_match | allele_flip]
  x[, BETA_BRADLEY_ALIGNED := fifelse(allele_match, BETA_BRADLEY, -BETA_BRADLEY)]
  x[, A1 := A1_FINNGEN]
  x[, A2 := A2_FINNGEN]
  x[, CHR := fifelse(!is.na(CHR_FINNGEN), CHR_FINNGEN, CHR_BRADLEY)]
  x[, BP := fifelse(!is.na(BP_FINNGEN), BP_FINNGEN, BP_BRADLEY)]
  x[, PALINDROMIC := is_palindromic(A1, A2)]
  x[, SAME_DIRECTION := sign(BETA_FINNGEN) == sign(BETA_BRADLEY_ALIGNED)]
  x[, APOE_REGION := CHR == APOE_CHR & BP >= APOE_START & BP <= APOE_END]
  x[, BRADLEY_DATASET := label]
  x[]
}

directional_summary <- function(h, label) {
  thresholds <- c(5e-8, 1e-6, 1e-5, 1e-4, 1e-3)
  rbindlist(lapply(thresholds, function(th) {
    rbindlist(lapply(c("APOE_included", "APOE_excluded"), function(analysis) {
      x <- h[P_FINNGEN < th]
      if (analysis == "APOE_excluded") x <- x[not_apoe(APOE_REGION)]
      if (nrow(x) == 0) {
        return(data.table(
          bradley_dataset = label, analysis = analysis, finngen_p_threshold = th,
          n_snps = 0L, same_direction_n = 0L, same_direction_fraction = NA_real_,
          binomial_p = NA_real_, pearson_r = NA_real_, pearson_p = NA_real_
        ))
      }
      ct <- binom.test(sum(x$SAME_DIRECTION, na.rm = TRUE), nrow(x), p = 0.5, alternative = "greater")
      cor_test <- if (nrow(x) >= 3) cor.test(x$BETA_FINNGEN, x$BETA_BRADLEY_ALIGNED, method = "pearson") else NULL
      data.table(
        bradley_dataset = label,
        analysis = analysis,
        finngen_p_threshold = th,
        n_snps = nrow(x),
        same_direction_n = sum(x$SAME_DIRECTION, na.rm = TRUE),
        same_direction_fraction = mean(x$SAME_DIRECTION, na.rm = TRUE),
        binomial_p = ct$p.value,
        pearson_r = if (!is.null(cor_test)) unname(cor_test$estimate) else NA_real_,
        pearson_p = if (!is.null(cor_test)) cor_test$p.value else NA_real_
      )
    }), fill = TRUE)
  }), fill = TRUE)[
    , directional_fdr := p.adjust(binomial_p, "BH"), by = .(bradley_dataset, analysis)
  ][]
}

define_lead_loci <- function(finngen) {
  x <- finngen[!is.na(CHR) & !is.na(BP) & P < LEAD_P_THRESHOLD]
  if (nrow(x) == 0) return(data.table())
  setorder(x, CHR, P)
  leads <- list()
  used <- rep(FALSE, nrow(x))
  for (i in seq_len(nrow(x))) {
    if (used[i]) next
    lead <- x[i]
    in_window <- x$CHR == lead$CHR & abs(x$BP - lead$BP) <= LEAD_WINDOW_BP
    used <- used | in_window
    leads[[length(leads) + 1]] <- data.table(
      locus_id = paste0("chr", lead$CHR, ":", lead$BP - LEAD_WINDOW_BP, "-", lead$BP + LEAD_WINDOW_BP),
      lead_snp = lead$SNP,
      CHR = lead$CHR,
      lead_bp = lead$BP,
      finngen_beta = lead$BETA,
      finngen_p = lead$P,
      apoe_region = lead$CHR == APOE_CHR & lead$BP >= APOE_START & lead$BP <= APOE_END
    )
  }
  rbindlist(leads, fill = TRUE)
}

locus_replication <- function(leads, harmonized, label) {
  if (nrow(leads) == 0) return(data.table())
  rbindlist(lapply(seq_len(nrow(leads)), function(i) {
    lead <- leads[i]
    x <- harmonized[
      CHR == lead$CHR &
        BP >= lead$lead_bp - LOCUS_SUPPORT_WINDOW_BP &
        BP <= lead$lead_bp + LOCUS_SUPPORT_WINDOW_BP
    ]
    lead_hit <- harmonized[SNP == lead$lead_snp]
    best <- if (nrow(x) > 0) x[which.min(P_BRADLEY)] else data.table()
    data.table(
      bradley_dataset = label,
      locus_id = lead$locus_id,
      lead_snp = lead$lead_snp,
      CHR = lead$CHR,
      lead_bp = lead$lead_bp,
      finngen_p = lead$finngen_p,
      finngen_beta = lead$finngen_beta,
      apoe_region = lead$apoe_region,
      n_harmonized_snps_in_window = nrow(x),
      lead_snp_present_in_bradley = nrow(lead_hit) > 0,
      lead_snp_bradley_p = if (nrow(lead_hit) > 0) lead_hit$P_BRADLEY[1] else NA_real_,
      lead_snp_same_direction = if (nrow(lead_hit) > 0) lead_hit$SAME_DIRECTION[1] else NA,
      best_bradley_snp = if (nrow(best) > 0) best$SNP[1] else NA_character_,
      best_bradley_p = if (nrow(best) > 0) best$P_BRADLEY[1] else NA_real_,
      best_bradley_same_direction = if (nrow(best) > 0) sign(lead$finngen_beta) == sign(best$BETA_BRADLEY_ALIGNED[1]) else NA
    )
  }), fill = TRUE)
}

read_gene_loc <- function(path) {
  if (!file.exists(path)) {
    warning("Gene-location file not found; candidate gene locus support will be skipped: ", path)
    return(NULL)
  }
  gm <- fread(path, header = FALSE, fill = TRUE)
  if (ncol(gm) < 4) return(NULL)
  setnames(gm, seq_len(ncol(gm)), paste0("V", seq_len(ncol(gm))))
  data.table(
    GENE = toupper(as.character(gm$V1)),
    CHR = clean_chr(gm$V2),
    START = as.integer(gm$V3),
    STOP = as.integer(gm$V4),
    STRAND = if (ncol(gm) >= 5) as.character(gm$V5) else NA_character_,
    SYMBOL = if (ncol(gm) >= 6) toupper(as.character(gm$V6)) else toupper(as.character(gm$V1))
  )[!is.na(CHR) & !is.na(START) & !is.na(STOP)]
}

read_candidate_gene_sets <- function(path) {
  if (!file.exists(path)) {
    warning("Candidate gene-set file not found; candidate locus support will be skipped: ", path)
    return(NULL)
  }
  gs <- fread(path)
  set_col <- detect_col(gs, c("set_name", "pathway", "gene_set", "module"), label = "candidate set column")
  gene_col <- detect_col(gs, c("GENE", "gene", "gene_id", "magma_gene_id", "symbol", "gene_symbol"), label = "candidate gene column")
  unique(data.table(set_name = as.character(gs[[set_col]]), gene_token = toupper(as.character(gs[[gene_col]]))))
}

candidate_gene_locus_support <- function(h, label, candidate_sets, gene_loc) {
  if (is.null(candidate_sets) || is.null(gene_loc)) return(NULL)
  gs_loc <- merge(
    candidate_sets,
    gene_loc,
    by.x = "gene_token",
    by.y = "GENE",
    all.x = FALSE,
    allow.cartesian = TRUE
  )
  if (nrow(gs_loc) == 0) {
    gs_loc <- merge(
      candidate_sets,
      gene_loc,
      by.x = "gene_token",
      by.y = "SYMBOL",
      all.x = FALSE,
      allow.cartesian = TRUE
    )
  }
  if (nrow(gs_loc) == 0) return(NULL)
  gs_loc[, `:=`(WINDOW_START = pmax(1L, START - GENE_WINDOW_BP), WINDOW_STOP = STOP + GENE_WINDOW_BP)]

  out <- rbindlist(lapply(seq_len(nrow(gs_loc)), function(i) {
    g <- gs_loc[i]
    x <- h[
      not_apoe(APOE_REGION) &
        CHR == g$CHR &
        BP >= g$WINDOW_START &
        BP <= g$WINDOW_STOP
    ]
    if (nrow(x) == 0) {
      return(data.table(
        bradley_dataset = label,
        set_name = g$set_name,
        gene = g$gene_token,
        CHR = g$CHR,
        START = g$START,
        STOP = g$STOP,
        n_snps = 0L,
        min_finngen_p = NA_real_,
        min_bradley_p = NA_real_,
        same_direction_fraction = NA_real_,
        status = "no_harmonized_snps"
      ))
    }
    data.table(
      bradley_dataset = label,
      set_name = g$set_name,
      gene = g$gene_token,
      CHR = g$CHR,
      START = g$START,
      STOP = g$STOP,
      n_snps = nrow(x),
      min_finngen_p = min(x$P_FINNGEN, na.rm = TRUE),
      min_bradley_p = min(x$P_BRADLEY, na.rm = TRUE),
      same_direction_fraction = mean(x$SAME_DIRECTION, na.rm = TRUE),
      status = "tested"
    )
  }), fill = TRUE)
  out[, bradley_min_p_fdr := p.adjust(min_bradley_p, "BH"), by = .(bradley_dataset, set_name)]
  out[]
}

plot_snp_concordance <- function(h, label, out_name) {
  x <- h[P_FINNGEN < 1e-5 & not_apoe(APOE_REGION)]
  if (nrow(x) == 0) return(invisible(NULL))
  p <- ggplot(x, aes(BETA_FINNGEN, BETA_BRADLEY_ALIGNED)) +
    geom_hline(yintercept = 0, color = "grey85") +
    geom_vline(xintercept = 0, color = "grey85") +
    geom_point(aes(color = P_BRADLEY < 0.05), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("TRUE" = "#b24745", "FALSE" = "#496f8a"), name = "Bradley P < 0.05") +
    labs(
      title = paste0("SNP-level effect concordance: FinnGen EOAD vs ", label),
      subtitle = "FinnGen P < 1e-5; APOE region excluded",
      x = "FinnGen EOAD beta",
      y = paste0(label, " aligned beta")
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(file.path(FIG_DIR, paste0(out_name, ".pdf")), p, width = 6, height = 5)
  ggsave(file.path(FIG_DIR, paste0(out_name, ".png")), p, width = 6, height = 5, dpi = 300)
}

plot_candidate_locus_summary <- function(candidate_support) {
  x <- candidate_support[status == "tested"]
  if (nrow(x) == 0) return(invisible(NULL))
  x_sum <- x[, .(
    n_genes_tested = uniqueN(gene),
    min_bradley_p = min(min_bradley_p, na.rm = TRUE),
    mean_same_direction_fraction = mean(same_direction_fraction, na.rm = TRUE),
    n_genes_bradley_p_lt_005 = sum(min_bradley_p < 0.05, na.rm = TRUE)
  ), by = .(bradley_dataset, set_name)]
  x_sum[, neglog10_min_p := -log10(pmax(min_bradley_p, .Machine$double.xmin))]
  p <- ggplot(x_sum, aes(bradley_dataset, set_name, fill = neglog10_min_p)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = n_genes_bradley_p_lt_005), size = 3.5) +
    scale_fill_gradient(low = "#e9eef2", high = "#b24745", name = "-log10 minimum Bradley P") +
    labs(
      title = "Candidate pathway gene-locus SNP support",
      subtitle = "APOE-excluded; tile label counts genes with minimum Bradley SNP P < 0.05",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(FIG_DIR, "Figure_08_candidate_gene_locus_SNP_support_APOE_excluded.pdf"), p, width = 8, height = 5.5)
  ggsave(file.path(FIG_DIR, "Figure_08_candidate_gene_locus_SNP_support_APOE_excluded.png"), p, width = 8, height = 5.5, dpi = 300)
}

# ------------------------------------------------------------------------------
# 3. Read, harmonize, and summarize
# ------------------------------------------------------------------------------

message_step("08. Bradley/Pottier 2025 SNP/locus-level replication")

if (!file.exists(FINNGEN_EOAD_SUMSTATS)) stop("FinnGen EOAD summary statistics not found: ", FINNGEN_EOAD_SUMSTATS)
if (!file.exists(BRADLEY_AGE65_INPUT)) stop("Bradley Age65 input not found. Run script 07 first or set BRADLEY_AGE65_MAGMA_INPUT.")
if (!file.exists(BRADLEY_AGE70_INPUT)) stop("Bradley Age70 input not found. Run script 07 first or set BRADLEY_AGE70_MAGMA_INPUT.")

finngen <- standardise_for_harmonization(FINNGEN_EOAD_SUMSTATS, "FinnGen_EOAD")
bradley65 <- standardise_for_harmonization(BRADLEY_AGE65_INPUT, "Bradley_NHW_Age65")
bradley70 <- standardise_for_harmonization(BRADLEY_AGE70_INPUT, "Bradley_NHW_Age70")

h65 <- harmonize_pair(finngen, bradley65, "Bradley_NHW_Age65")
h70 <- harmonize_pair(finngen, bradley70, "Bradley_NHW_Age70")

harm_qc <- data.table(
  bradley_dataset = c("Bradley_NHW_Age65", "Bradley_NHW_Age70"),
  n_finngen_variants = nrow(finngen),
  n_bradley_variants = c(nrow(bradley65), nrow(bradley70)),
  n_harmonized_variants = c(nrow(h65), nrow(h70)),
  n_harmonized_non_palindromic = c(sum(!h65$PALINDROMIC), sum(!h70$PALINDROMIC)),
  n_harmonized_apoe_excluded = c(sum(!h65$APOE_REGION), sum(!h70$APOE_REGION))
)
fwrite(harm_qc, file.path(TAB_DIR, "Table_08_SNP_harmonization_QC.tsv"), sep = "\t")

fwrite(h65[order(P_FINNGEN)][seq_len(min(nrow(h65), TOP_N_FOR_HARMONIZED_EXPORT))], file.path(TAB_DIR, "Table_08_harmonized_FinnGen_Bradley_Age65_top_by_FinnGen_p.tsv"), sep = "\t")
fwrite(h70[order(P_FINNGEN)][seq_len(min(nrow(h70), TOP_N_FOR_HARMONIZED_EXPORT))], file.path(TAB_DIR, "Table_08_harmonized_FinnGen_Bradley_Age70_top_by_FinnGen_p.tsv"), sep = "\t")

direction_summary <- rbindlist(list(
  directional_summary(h65, "Bradley_NHW_Age65"),
  directional_summary(h70, "Bradley_NHW_Age70")
), fill = TRUE)
fwrite(direction_summary, file.path(TAB_DIR, "Table_08_SNP_directional_replication_summary.tsv"), sep = "\t")

leads <- define_lead_loci(finngen)
fwrite(leads, file.path(TAB_DIR, "Table_08_FinnGen_EOAD_independent_lead_SNPs.tsv"), sep = "\t")
locus_results <- rbindlist(list(
  locus_replication(leads, h65, "Bradley_NHW_Age65"),
  locus_replication(leads, h70, "Bradley_NHW_Age70")
), fill = TRUE)
if (nrow(locus_results) > 0) {
  fwrite(locus_results, file.path(TAB_DIR, "Table_08_Bradley_locus_level_replication_of_FinnGen_leads.tsv"), sep = "\t")
  locus_summary <- locus_results[, .(
    n_loci = .N,
    n_non_apoe_loci = sum(!apoe_region),
    n_lead_snp_present = sum(lead_snp_present_in_bradley, na.rm = TRUE),
    n_lead_snp_same_direction = sum(lead_snp_same_direction, na.rm = TRUE),
    n_loci_best_bradley_p_lt_005 = sum(best_bradley_p < 0.05, na.rm = TRUE),
    n_loci_best_same_direction = sum(best_bradley_same_direction, na.rm = TRUE)
  ), by = bradley_dataset]
  fwrite(locus_summary, file.path(TAB_DIR, "Table_08_locus_replication_summary.tsv"), sep = "\t")
}

candidate_sets <- read_candidate_gene_sets(CANDIDATE_GENESET_FILE)
gene_loc <- read_gene_loc(GENE_LOC_FILE)
candidate_support <- rbindlist(list(
  candidate_gene_locus_support(h65, "Bradley_NHW_Age65", candidate_sets, gene_loc),
  candidate_gene_locus_support(h70, "Bradley_NHW_Age70", candidate_sets, gene_loc)
), fill = TRUE)
if (nrow(candidate_support) > 0) {
  fwrite(candidate_support, file.path(TAB_DIR, "Table_08_candidate_gene_locus_SNP_support_APOE_excluded.tsv"), sep = "\t")
  candidate_summary <- candidate_support[status == "tested", .(
    n_genes_tested = uniqueN(gene),
    n_genes_with_bradley_p_lt_005 = sum(min_bradley_p < 0.05, na.rm = TRUE),
    min_bradley_p = min(min_bradley_p, na.rm = TRUE),
    mean_same_direction_fraction = mean(same_direction_fraction, na.rm = TRUE)
  ), by = .(bradley_dataset, set_name)]
  candidate_summary[, min_bradley_p_fdr := p.adjust(min_bradley_p, "BH"), by = bradley_dataset]
  fwrite(candidate_summary, file.path(TAB_DIR, "Table_08_candidate_gene_locus_SNP_support_summary.tsv"), sep = "\t")
}

# ------------------------------------------------------------------------------
# 4. Figures
# ------------------------------------------------------------------------------

plot_snp_concordance(h65, "Bradley Age65", "Figure_08_SNP_effect_concordance_FinnGen_vs_Bradley_Age65_p1e5")
plot_snp_concordance(h70, "Bradley Age70", "Figure_08_SNP_effect_concordance_FinnGen_vs_Bradley_Age70_p1e5")
if (nrow(candidate_support) > 0) plot_candidate_locus_summary(candidate_support)

# ------------------------------------------------------------------------------
# 5. Run summary
# ------------------------------------------------------------------------------

summary_lines <- c(
  "08. Bradley/Pottier 2025 SNP/locus-level replication completed.",
  paste0("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE)),
  "",
  "Primary outputs:",
  file.path(TAB_DIR, "Table_08_SNP_directional_replication_summary.tsv"),
  file.path(TAB_DIR, "Table_08_Bradley_locus_level_replication_of_FinnGen_leads.tsv"),
  file.path(TAB_DIR, "Table_08_candidate_gene_locus_SNP_support_summary.tsv"),
  file.path(FIG_DIR, "Figure_08_SNP_effect_concordance_FinnGen_vs_Bradley_Age65_p1e5.png"),
  file.path(FIG_DIR, "Figure_08_candidate_gene_locus_SNP_support_APOE_excluded.png"),
  "",
  "Interpretation note:",
  "SNP/locus-level concordance is reported as directional and locus-window support.",
  "The module complements, but does not replace, the MAGMA gene/pathway replication in script 07."
)
writeLines(summary_lines, file.path(OUT_DIR, "README_08_SNP_locus_replication.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
