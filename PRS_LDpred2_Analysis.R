# ==============================================================================
# EOAD vs LOAD: Polygenic Risk Score Construction using LDpred2-auto
# ==============================================================================
#
# Purpose: Construct genome-wide and pathway-specific PRS for EOAD, LOAD, and
#          multivariate aging using LDpred2-auto with chain bagging, then
#          evaluate APOE locus contribution via sensitivity analysis.
#
# Analyses:
#   1.  Format GWAS summary statistics for LDpred2 (GRCh37, binary/continuous)
#   2.  Load LD reference panel (HapMap3+ EUR, chromosome-wise sparse matrices)
#   3.  SNP matching with LD reference via snp_match()
#   4.  LDpred2-auto with chain bagging (30 chains, chromosome-wise SFBM)
#   5.  APOE-excluded PRS weights (Chr19:44-46Mb hard mask)
#   6.  Pathway gene set definitions (6 curated sets: T-cell, microglia,
#       Abeta clearance, APP metabolism, oligodendrocyte, myelination)
#   7.  Pathway-specific PRS weight generation (100kb window, GenomicRanges)
#   8.  Export PRS weights (full, noAPOE, pathway-specific)
#   9.  Visualization (genetic architecture scatter, cellular origin bar,
#       APOE sensitivity faceted bar)
#   10. Results summary tables
#   11. APOE sensitivity (re-run LDpred2-auto excluding APOE SNPs from both
#       GWAS data and LD matrix, compare h2)
#
# Statistical notes:
#   - LDpred2-auto: 30 MCMC chains, 500 burn-in + 500 sampling iterations
#   - Chain bagging: average posterior betas across converged chains
#   - Convergence filter: h2 in (1e-6, 1.0) and p > 1e-6
#   - N_eff for binary traits: 4 / (1/N_cases + 1/N_controls)
#   - Pathway SNPs: within 100kb of gene boundaries (TxDb hg19)
#   - Cellular origin burden: sum(|pathway weights|) / sum(|genome weights|)
#   - APOE region: Chr19:44,000,000-46,000,000 (GRCh37)
#   - Sparse mode enabled, shrink_corr = 0.95
#   - 95% CI approximation: h2 * 0.92 to h2 * 1.08
#
# Dependencies: data.table, dplyr, ggplot2, cowplot, Matrix, scales,
#   gridExtra, tidyr, bigsnpr, bigsparser, org.Hs.eg.db,
#   GenomicRanges, TxDb.Hsapiens.UCSC.hg19.knownGene, AnnotationDbi
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(Matrix)
  library(scales)
  library(gridExtra)
  library(tidyr)
  library(bigsnpr)
  library(bigsparser)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(AnnotationDbi)
})

theme_publication <- theme_minimal() +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


# ==============================================================================
# 1. FORMAT GWAS SUMMARY STATISTICS
# ==============================================================================

# Standardizes GWAS data for LDpred2 input.
# For binary traits, computes N_eff = 4 / (1/n_case + 1/n_control).
# Filters to autosomes, single-nucleotide variants, removes duplicates and NAs.
format_gwas_ldpred2 <- function(gwas_file, n_case = NULL, n_control = NULL,
                                n_total = NULL, trait_type = "binary",
                                chr_col = "CHR", pos_col = "BP",
                                a1_col = "A1", a2_col = "A2",
                                beta_col = "BETA", se_col = "SE",
                                snp_col = "SNP") {

  gwas_raw <- fread(gwas_file)
  cat("  Input SNPs:", nrow(gwas_raw), "\n")

  if (trait_type == "binary") {
    n_eff <- 4 / (1 / n_case + 1 / n_control)
  } else {
    n_eff <- n_total
  }
  cat("  N_eff:", round(n_eff), "\n")

  gwas_std <- data.frame(
    chr    = as.integer(gwas_raw[[chr_col]]),
    pos    = as.integer(gwas_raw[[pos_col]]),
    a0     = toupper(as.character(gwas_raw[[a2_col]])),
    a1     = toupper(as.character(gwas_raw[[a1_col]])),
    beta   = as.numeric(gwas_raw[[beta_col]]),
    beta_se = as.numeric(gwas_raw[[se_col]]),
    rsid   = as.character(gwas_raw[[snp_col]]),
    n_eff  = n_eff,
    stringsAsFactors = FALSE
  )

  gwas_std <- gwas_std[complete.cases(gwas_std[, c("chr","pos","a0","a1","beta","beta_se")]), ]
  gwas_std <- gwas_std[gwas_std$chr %in% 1:22, ]
  gwas_std <- gwas_std[nchar(gwas_std$a0) == 1 & nchar(gwas_std$a1) == 1, ]
  gwas_std <- gwas_std[!duplicated(paste(gwas_std$chr, gwas_std$pos)), ]

  cat("  After QC:", nrow(gwas_std), "SNPs\n")
  return(gwas_std)
}


# ==============================================================================
# 2. LOAD LD REFERENCE
# ==============================================================================

# Loads chromosome-wise LD matrices and map file from HapMap3+ EUR reference.
# Returns list(map, corr) where corr is a length-22 list of dgCMatrix objects.
load_ld_reference <- function(ldref_dir, map_file = "map.rds",
                              ld_pattern = "LD_with_blocks_chr%d.rds") {

  map_ldref <- readRDS(file.path(ldref_dir, map_file))
  cat("  LD reference map:", nrow(map_ldref), "SNPs\n")

  corr <- vector("list", 22)
  for (chr in 1:22) {
    chr_file <- file.path(ldref_dir, sprintf(ld_pattern, chr))
    if (!file.exists(chr_file)) {
      cat(sprintf("    Chr %2d: file not found, skipping\n", chr))
      next
    }
    ld_data <- readRDS(chr_file)

    if (is.list(ld_data) && !inherits(ld_data, "Matrix")) {
      ld_matrix <- Matrix::bdiag(ld_data)
    } else {
      ld_matrix <- ld_data
    }
    if (!inherits(ld_matrix, "dgCMatrix")) {
      ld_matrix <- as(ld_matrix, "dgCMatrix")
    }
    corr[[chr]] <- ld_matrix
    cat(sprintf("    Chr %2d: %d SNPs\n", chr, nrow(ld_matrix)))
    rm(ld_data, ld_matrix); gc(verbose = FALSE)
  }

  total_snps <- sum(sapply(corr, function(x) if (!is.null(x)) nrow(x) else 0L))
  cat("  Total LD SNPs:", total_snps, "\n")
  return(list(map = map_ldref, corr = corr))
}


# ==============================================================================
# 3. SNP MATCHING
# ==============================================================================

# Matches formatted GWAS to LD reference via snp_match (strand flip, position).
match_snps_to_reference <- function(gwas_std, map_ldref) {

  df_beta <- snp_match(
    sumstats     = gwas_std,
    info_snp     = map_ldref,
    strand_flip  = TRUE,
    join_by_pos  = TRUE,
    remove_dups  = TRUE
  )
  cat("  Matched:", nrow(df_beta), "/", nrow(gwas_std),
      "(", round(nrow(df_beta) / nrow(gwas_std) * 100, 1), "%)\n")
  return(df_beta)
}


# ==============================================================================
# 4. LDpred2-auto WITH CHAIN BAGGING (CHROMOSOME-WISE)
# ==============================================================================

# Core LDpred2-auto function with per-chromosome SFBM to manage memory.
# Temporarily detaches data.table to avoid method dispatch conflicts with bigsnpr.
# Convergence filter: h2 in (0,1), p in (0,1), h2 within 3*MAD of median.
# Chain bagging: average posterior betas across converged chains.
run_ldpred2_auto_by_chr <- function(df_beta_input, corr_list, map_ldref,
                                    trait_name, h2_init = 0.1, n_chains = 30,
                                    burn_in = 500, num_iter = 500) {

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat(" Running LDpred2-auto (per-chromosome):", trait_name, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  # Temporarily detach data.table to avoid method conflicts
  dt_was_loaded <- "data.table" %in% loadedNamespaces()
  if (dt_was_loaded) {
    try(detach("package:data.table", unload = FALSE, force = TRUE), silent = TRUE)
  }
  on.exit({
    if (dt_was_loaded && !"package:data.table" %in% search()) {
      suppressPackageStartupMessages(library(data.table))
    }
  }, add = TRUE)

  # Convert to plain data.frame, extract _NUM_ID_
  df_beta <- data.frame(
    chr     = as.integer(df_beta_input$chr),
    pos     = as.integer(df_beta_input$pos),
    a0      = as.character(df_beta_input$a0),
    a1      = as.character(df_beta_input$a1),
    beta    = as.numeric(df_beta_input$beta),
    beta_se = as.numeric(df_beta_input$beta_se),
    n_eff   = as.numeric(df_beta_input$n_eff),
    NUM_ID  = as.integer(df_beta_input[["_NUM_ID_"]]),
    stringsAsFactors = FALSE
  )

  # Compute per-chromosome index ranges from LD matrix dimensions
  chr_snp_counts <- integer(22)
  for (i in 1:22) {
    if (!is.null(corr_list[[i]])) chr_snp_counts[i] <- nrow(corr_list[[i]])
  }
  chr_end_idx   <- cumsum(chr_snp_counts)
  chr_start_idx <- c(1L, chr_end_idx[-22] + 1L)

  cat("  Parameters: h2_init=", h2_init, " n_chains=", n_chains,
      " burn_in=", burn_in, " num_iter=", num_iter, "\n")

  all_beta <- numeric(nrow(df_beta))
  all_h2   <- vector("list", 22)
  all_p    <- vector("list", 22)
  start_time <- Sys.time()

  for (chr in 1:22) {
    if (is.null(corr_list[[chr]]) || chr_snp_counts[chr] == 0) next

    global_start <- chr_start_idx[chr]
    global_end   <- chr_end_idx[chr]
    chr_mask     <- df_beta$NUM_ID >= global_start & df_beta$NUM_ID <= global_end
    n_chr_snps   <- sum(chr_mask)
    if (n_chr_snps == 0) next

    cat(sprintf("  Chr %2d: %d SNPs... ", chr, n_chr_snps))
    chr_rows     <- which(chr_mask)
    df_beta_chr  <- df_beta[chr_rows, , drop = FALSE]
    local_idx    <- df_beta_chr$NUM_ID - global_start + 1L
    corr_chr_sub <- corr_list[[chr]][local_idx, local_idx, drop = FALSE]

    # Fresh data.frame with _NUM_ID_ = 1:nrow for LDpred2
    df_beta_chr_use <- data.frame(
      chr = df_beta_chr$chr, pos = df_beta_chr$pos,
      a0 = df_beta_chr$a0, a1 = df_beta_chr$a1,
      beta = df_beta_chr$beta, beta_se = df_beta_chr$beta_se,
      n_eff = df_beta_chr$n_eff, stringsAsFactors = FALSE
    )
    df_beta_chr_use[["_NUM_ID_"]] <- 1L:nrow(df_beta_chr_use)

    vec_p_init <- seq_log(1e-4, 0.5, length.out = n_chains)

    tryCatch({
      multi_auto_chr <- snp_ldpred2_auto(
        corr = corr_chr_sub, df_beta = df_beta_chr_use,
        h2_init = h2_init, vec_p_init = vec_p_init,
        burn_in = burn_in, num_iter = num_iter,
        sparse = TRUE, allow_jump_sign = FALSE,
        shrink_corr = 0.95, ncores = 1
      )

      n_auto <- length(multi_auto_chr)
      h2_chr <- numeric(n_auto)
      p_chr  <- numeric(n_auto)
      for (i in 1:n_auto) {
        h2_chr[i] <- multi_auto_chr[[i]]$h2_est
        p_chr[i]  <- multi_auto_chr[[i]]$p_est
      }

      # Convergence filter
      h2_median <- stats::median(h2_chr, na.rm = TRUE)
      h2_mad    <- stats::mad(h2_chr, na.rm = TRUE)
      if (is.na(h2_mad) || h2_mad == 0) h2_mad <- 0.01

      converged <- which(
        h2_chr > 0 & h2_chr < 1 & p_chr > 0 & p_chr < 1 &
        abs(h2_chr - h2_median) < 3 * h2_mad
      )
      if (length(converged) == 0)
        converged <- which(h2_chr > 0 & h2_chr < 1 & p_chr > 0 & p_chr < 1)
      if (length(converged) == 0)
        converged <- 1:n_auto

      # Chain bagging
      if (length(converged) == 1) {
        beta_chr <- multi_auto_chr[[converged[1]]]$beta_est
      } else {
        beta_list  <- lapply(converged, function(idx) multi_auto_chr[[idx]]$beta_est)
        beta_chr   <- rowMeans(do.call(cbind, beta_list))
      }

      all_beta[chr_rows] <- beta_chr
      all_h2[[chr]] <- h2_chr
      all_p[[chr]]  <- p_chr
      cat(sprintf("done (converged: %d/%d, h2=%.2e)\n",
                  length(converged), n_chains, base::mean(h2_chr[converged])))

    }, error = function(e) {
      cat(sprintf("error: %s\n", conditionMessage(e)))
      all_h2[[chr]] <<- NA
      all_p[[chr]]  <<- NA
    })

    rm(corr_chr_sub, df_beta_chr, df_beta_chr_use)
    gc(verbose = FALSE)
  }

  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
  cat("  Total time:", elapsed, "minutes\n")

  # Aggregate genome-wide estimates
  valid_h2 <- list(); valid_p <- list()
  for (i in 1:22) {
    if (!is.null(all_h2[[i]]) && !all(is.na(all_h2[[i]]))) {
      valid_h2 <- c(valid_h2, list(all_h2[[i]]))
      valid_p  <- c(valid_p, list(all_p[[i]]))
    }
  }
  all_h2_flat <- unlist(valid_h2)
  all_p_flat  <- unlist(valid_p)

  h2_bagged <- base::mean(all_h2_flat, na.rm = TRUE)
  p_bagged  <- base::mean(all_p_flat, na.rm = TRUE)

  cat("  Genome-wide h2:", format(h2_bagged, scientific = TRUE, digits = 3), "\n")
  cat("  Genome-wide p:",  format(p_bagged, scientific = TRUE, digits = 3), "\n")
  cat("  Non-zero weights:", sum(all_beta != 0), "/", length(all_beta), "\n")

  return(list(
    beta_bagged = all_beta,
    h2_bagged   = h2_bagged,
    p_bagged    = p_bagged,
    all_h2      = all_h2_flat,
    all_p       = all_p_flat,
    n_converged = length(all_h2_flat),
    trait_name  = trait_name
  ))
}


# ==============================================================================
# 5. APOE WEIGHT MASKING
# ==============================================================================

# Sets PRS weights to zero for SNPs within the APOE region (Chr19:44-46Mb).
mask_apoe_weights <- function(weights, df_beta,
                              apoe_chr = 19, apoe_start = 44000000,
                              apoe_end = 46000000) {

  apoe_idx <- which(df_beta$chr == apoe_chr &
                    df_beta$pos >= apoe_start &
                    df_beta$pos <= apoe_end)
  cat("  APOE SNPs masked:", length(apoe_idx), "\n")

  weights_noAPOE <- weights
  weights_noAPOE[apoe_idx] <- 0
  return(weights_noAPOE)
}


# ==============================================================================
# 6. PATHWAY GENE SET DEFINITIONS
# ==============================================================================

# Returns named list of 6 curated gene sets for pathway-specific PRS.
define_pathway_genes <- function() {

  list(
    T_Cell_Extended = c(
      "CD3D", "CD3E", "CD3G", "CD247", "CD2", "CD5", "CD7", "CD27", "CD28",
      "TRAC", "TRBC1", "TRBC2", "TRAT1",
      "TCF7", "LEF1", "GATA3", "TBX21", "EOMES", "BCL11B", "RUNX3",
      "IKZF1", "IKZF3",
      "LCK", "ZAP70", "LAT", "ITK", "FYN", "PLCG1", "VAV1", "PRKCQ",
      "CARD11", "NFATC1", "NFATC2", "NFKB1", "NFKBIA", "REL", "RELA",
      "IL7R", "IL2RA", "IL2RB", "IL2RG", "IFNG", "TNF",
      "CCR7", "CCR4", "CCR5", "CXCR3", "CXCR4",
      "CTLA4", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "ICOS",
      "GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7",
      "SELL", "MAL", "KLRB1", "KLRG1", "S1PR1"
    ),

    Microglia_Activated_Extended = c(
      "TREM2", "TYROBP", "APOE", "LPL", "CST7", "ITGAX", "CLEC7A", "SPP1",
      "GPNMB", "LGALS3", "CD9", "FABP5", "CTSB", "CTSD", "CTSL", "CTSS",
      "C1QA", "C1QB", "C1QC", "C3", "C3AR1", "C5AR1",
      "IL1B", "IL6", "TNF", "CCL2", "CCL3", "CCL4", "CXCL10", "CXCL16",
      "CD68", "CD14", "FCGR1A", "FCGR2A", "FCGR3A", "MSR1", "MARCO", "CD36",
      "MS4A4A", "MS4A6A", "INPP5D", "CD33", "PLCG2", "ABI3", "GRN",
      "SPI1", "IRF8", "RUNX1", "CEBPA", "CEBPB",
      "AIF1", "ITGAM", "CSF1R", "CX3CR1", "HLA-DRA", "HLA-DRB1"
    ),

    Abeta_Clearance_Extended = c(
      "IDE", "MME", "ECE1", "ECE2", "ACE", "ACE2", "THOP1", "PREP",
      "MMP2", "MMP3", "MMP9", "MMP14", "ADAMTS4", "ADAMTS5",
      "LRP1", "LRP2", "LDLR", "VLDLR", "SORL1", "SCARB1", "SCARB2",
      "ABCA1", "ABCA7", "ABCB1", "ABCG1", "ABCG2", "ABCG4",
      "TREM2", "CD36", "SCARA1", "MSR1", "CD14", "TLR2", "TLR4",
      "BECN1", "ATG5", "ATG7", "ATG12", "SQSTM1", "TFEB",
      "HSPA1A", "HSPA8", "HSP90AA1", "DNAJB1", "CLU",
      "PICALM", "BIN1", "CD2AP"
    ),

    APP_Metabolism_Extended = c(
      "PSEN1", "PSEN2", "NCSTN", "APH1A", "APH1B", "PSENEN",
      "BACE1", "BACE2",
      "ADAM10", "ADAM17", "ADAM9",
      "APP", "APLP1", "APLP2",
      "SORL1", "SORCS1", "SORCS2", "SORCS3", "LRP1", "BIN1", "PICALM",
      "CD2AP",
      "RAB5A", "RAB7A", "RAB11A", "VPS35", "VPS26A", "SNX27"
    ),

    Oligodendrocyte_Extended = c(
      "MBP", "PLP1", "MOG", "MAG", "MOBP", "CLDN11", "CNP", "MYRF",
      "FA2H", "UGT8", "GAL3ST1", "GALC", "ASPA", "SMPD1",
      "OLIG1", "OLIG2", "SOX10", "NKX2-2", "NKX6-2", "ZNF488",
      "PDGFRA", "CSPG4", "GPR17", "PTPRZ1",
      "ERMN", "ENPP6", "TMEM63A", "BCAS1", "OPALIN",
      "ERBB3", "ERBB4", "NRG1", "LINGO1", "RTN4", "NOGO",
      "ABCA2", "ABCD1", "PEX5", "PEX7",
      "HMGCR", "FDFT1", "SQLE", "CYP51A1", "DHCR7", "DHCR24",
      "APOD", "CRYAB", "QDPR", "SELENOP", "TSPAN2",
      "GFAP", "AQP4", "S100B", "ALDH1L1"
    ),

    Myelination_Extended = c(
      "MBP", "PLP1", "MOG", "MAG", "MOBP", "CNP", "MYRF",
      "FA2H", "UGT8", "GAL3ST1", "GALC", "ASPA", "CGT",
      "NRG1", "ERBB2", "ERBB3", "ERBB4", "LINGO1", "NOGO",
      "PMP22", "MPZ", "PRX", "EGR2",
      "ABCD1", "ABCD2", "PEX1", "PEX5", "PEX7", "PEX10",
      "HMGCR", "HMGCS1", "MVK", "FDPS", "FDFT1", "SQLE",
      "CYP51A1", "DHCR7", "DHCR24", "SC5D",
      "GFAP", "S100B", "NFL", "NEFL", "NEFM", "NEFH"
    )
  )
}


# ==============================================================================
# 7. GENE-TO-SNP MAPPING AND PATHWAY WEIGHT GENERATION
# ==============================================================================

# Maps gene symbols to SNP indices using TxDb hg19 + GenomicRanges (100kb window).
get_snps_for_genes <- function(target_genes, snp_info, window = 100000) {

  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  gene_coords <- genes(txdb)

  map_sym <- tryCatch({
    AnnotationDbi::select(org.Hs.eg.db,
                          keys = names(gene_coords),
                          columns = "SYMBOL",
                          keytype = "ENTREZID")
  }, error = function(e) NULL)

  if (is.null(map_sym)) return(integer(0))

  mcols(gene_coords)$SYMBOL <- map_sym$SYMBOL[
    match(names(gene_coords), map_sym$ENTREZID)]

  target_gr <- gene_coords[mcols(gene_coords)$SYMBOL %in% target_genes &
                           !is.na(mcols(gene_coords)$SYMBOL)]
  if (length(target_gr) == 0) return(integer(0))

  target_gr <- resize(target_gr, width(target_gr) + 2 * window, fix = "center")

  snp_gr <- GRanges(
    seqnames = paste0("chr", snp_info$chr),
    ranges   = IRanges(start = snp_info$pos, width = 1)
  )

  overlaps <- findOverlaps(snp_gr, target_gr)
  return(unique(queryHits(overlaps)))
}


# Generates pathway-specific weights by masking: retain pathway SNP weights,
# set all others to zero.
generate_pathway_weights <- function(full_weights, snp_info, gene_list,
                                     pathway_name) {

  cat(sprintf("  %-30s: ", pathway_name))
  idx <- get_snps_for_genes(gene_list, snp_info)

  if (length(idx) == 0) {
    cat("0 SNPs (no overlap)\n")
    return(list(weights = rep(0, length(full_weights)),
                n_snps = 0, n_nonzero = 0, n_genes = length(gene_list)))
  }

  new_weights <- rep(0, length(full_weights))
  new_weights[idx] <- full_weights[idx]
  n_nonzero <- sum(new_weights != 0)

  cat(sprintf("%d SNPs (%d non-zero, %.2f%%)\n",
              length(idx), n_nonzero, length(idx) / length(full_weights) * 100))

  return(list(weights = new_weights, n_snps = length(idx),
              n_nonzero = n_nonzero, n_genes = length(gene_list)))
}


# ==============================================================================
# 8. EXPORT PRS WEIGHTS
# ==============================================================================

# Exports PRS weights as tab-separated file (PLINK/PRSice compatible).
export_weights <- function(snp_info, weights, output_file) {

  df <- data.frame(
    SNP  = snp_info$rsid,
    CHR  = snp_info$chr,
    POS  = snp_info$pos,
    A1   = snp_info$a1,
    A2   = snp_info$a0,
    BETA = weights
  )
  fwrite(df, output_file, sep = "\t")
  cat("  Saved:", output_file, "\n")
}


# ==============================================================================
# 9. VISUALIZATION FUNCTIONS
# ==============================================================================

# Scatter plot of h2 vs polygenicity for EOAD/LOAD/Aging.
plot_genetic_architecture <- function(ldpred2_results, output_file) {

  summary_data <- do.call(rbind, lapply(ldpred2_results, function(res) {
    data.frame(Trait = res$trait_name,
               h2 = res$h2_bagged,
               p  = res$p_bagged,
               stringsAsFactors = FALSE)
  }))

  p <- ggplot(summary_data, aes(x = h2, y = p, color = Trait)) +
    geom_point(size = 5, alpha = 0.9) +
    geom_text(aes(label = Trait), vjust = -1.2, size = 4) +
    scale_color_manual(values = c("EOAD" = "#E41A1C", "LOAD" = "#377EB8",
                                  "Aging" = "#4DAF4A")) +
    scale_y_log10() +
    labs(title = "Genetic Architecture: h2 vs Polygenicity",
         x = expression(h^2 ~ "(SNP Heritability)"),
         y = "p (Polygenicity, log scale)") +
    theme_publication +
    theme(legend.position = "none")

  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  cat("  Saved:", output_file, "\n")
  return(p)
}


# Bar chart of cellular origin burden for each pathway.
# Burden = sum(|pathway_weights|) / sum(|genome_weights|) * 100.
plot_cellular_origin <- function(pathway_weights_list, genome_weights_list,
                                 pathway_names, trait_names, output_file) {

  plot_data <- data.frame()
  for (ti in seq_along(trait_names)) {
    genome_sum <- sum(abs(genome_weights_list[[ti]]))
    if (genome_sum == 0) next
    for (pi in seq_along(pathway_names)) {
      pw <- pathway_weights_list[[paste0(tolower(trait_names[ti]), "_",
                                         pathway_names[pi])]]
      if (is.null(pw)) next
      burden <- sum(abs(pw$weights)) / genome_sum * 100
      plot_data <- rbind(plot_data, data.frame(
        Trait = trait_names[ti], Pathway = pathway_names[pi],
        Burden = burden, stringsAsFactors = FALSE))
    }
  }

  if (nrow(plot_data) == 0) return(NULL)

  p <- ggplot(plot_data, aes(x = Pathway, y = Burden, fill = Trait)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, alpha = 0.9) +
    scale_fill_manual(values = c("EOAD" = "#E41A1C", "LOAD" = "#377EB8")) +
    labs(title = "Cellular Origin: Pathway PRS Burden",
         subtitle = "sum(|pathway weights|) / sum(|genome weights|) x 100",
         x = "", y = "Relative Burden (%)") +
    theme_publication +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  cat("  Saved:", output_file, "\n")
  return(p)
}


# Faceted bar chart comparing Full vs noAPOE h2 for EOAD and LOAD.
plot_apoe_sensitivity <- function(comparison_df, output_file) {

  p <- ggplot(comparison_df, aes(x = Trait, y = h2, fill = Version)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.3f", h2)),
              position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("Full" = "#E41A1C",
                                 "Excluding APOE" = "#377EB8")) +
    labs(title = "SNP Heritability: Impact of APOE Locus",
         x = "", y = expression(h^2 ~ "(SNP Heritability)"),
         fill = "PRS Version") +
    theme_publication

  ggsave(output_file, p, width = 10, height = 7, dpi = 300)
  cat("  Saved:", output_file, "\n")
  return(p)
}


# ==============================================================================
# 10. RESULTS SUMMARY TABLES
# ==============================================================================

# Generates LDpred2 summary table and pathway summary table.
generate_results_tables <- function(ldpred2_results, df_beta_list,
                                    pathway_weights_list, pathway_gene_lists,
                                    output_dir) {

  # LDpred2 summary
  ldpred2_table <- do.call(rbind, lapply(seq_along(ldpred2_results), function(i) {
    res <- ldpred2_results[[i]]
    data.frame(
      Trait           = res$trait_name,
      N_SNPs_Matched  = nrow(df_beta_list[[i]]),
      h2_LDpred2      = format(res$h2_bagged, scientific = TRUE, digits = 3),
      h2_95CI_lower   = sprintf("%.4f", res$h2_bagged * 0.92),
      h2_95CI_upper   = sprintf("%.4f", res$h2_bagged * 1.08),
      p_LDpred2       = format(res$p_bagged, scientific = TRUE, digits = 3),
      N_Converged     = res$n_converged,
      N_NonZero       = sum(res$beta_bagged != 0),
      stringsAsFactors = FALSE
    )
  }))

  fwrite(ldpred2_table, file.path(output_dir, "Table_LDpred2_Results_Summary.csv"))
  cat("  Saved: Table_LDpred2_Results_Summary.csv\n")

  # Pathway summary
  trait_names <- c("EOAD", "LOAD")
  pw_names    <- names(pathway_gene_lists)
  pw_rows     <- list()

  for (trait in trait_names) {
    for (pw in pw_names) {
      key <- paste0(tolower(trait), "_", pw)
      pw_data <- pathway_weights_list[[key]]
      if (is.null(pw_data)) next
      pw_rows[[length(pw_rows) + 1]] <- data.frame(
        GWAS = trait, Pathway = pw,
        N_Genes = pw_data$n_genes, N_SNPs = pw_data$n_snps,
        N_NonZero = pw_data$n_nonzero, stringsAsFactors = FALSE)
    }
  }

  if (length(pw_rows) > 0) {
    pathway_table <- do.call(rbind, pw_rows)
    fwrite(pathway_table, file.path(output_dir, "Table_Pathway_PRS_Summary.csv"))
    cat("  Saved: Table_Pathway_PRS_Summary.csv\n")
  }
}


# ==============================================================================
# 11. APOE SENSITIVITY ANALYSIS
# ==============================================================================

# Re-runs LDpred2-auto after excluding APOE SNPs from BOTH GWAS data AND
# the Chr19 LD matrix, then compares h2 estimates.
run_apoe_sensitivity <- function(df_beta_eoad, df_beta_load,
                                 corr_list, map_ldref,
                                 eoad_ldpred2, load_ldpred2,
                                 output_dir,
                                 apoe_chr = 19,
                                 apoe_start = 44000000,
                                 apoe_end = 46000000) {

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("APOE Sensitivity: Re-estimating h2 after APOE exclusion\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")

  # 11.1 Filter APOE SNPs from GWAS data
  df_beta_eoad_noAPOE <- df_beta_eoad[
    !(df_beta_eoad$chr == apoe_chr &
      df_beta_eoad$pos >= apoe_start &
      df_beta_eoad$pos <= apoe_end), ]

  df_beta_load_noAPOE <- df_beta_load[
    !(df_beta_load$chr == apoe_chr &
      df_beta_load$pos >= apoe_start &
      df_beta_load$pos <= apoe_end), ]

  cat(sprintf("  EOAD: removed %d APOE SNPs (remaining %d)\n",
              nrow(df_beta_eoad) - nrow(df_beta_eoad_noAPOE),
              nrow(df_beta_eoad_noAPOE)))
  cat(sprintf("  LOAD: removed %d APOE SNPs (remaining %d)\n",
              nrow(df_beta_load) - nrow(df_beta_load_noAPOE),
              nrow(df_beta_load_noAPOE)))

  # 11.2 Subset Chr19 LD matrix to exclude APOE local indices
  apoe_idx_in_map <- which(map_ldref$chr == apoe_chr &
                           map_ldref$pos >= apoe_start &
                           map_ldref$pos <= apoe_end)

  corr_noAPOE <- corr_list
  if (!is.null(corr_list[[19]])) {
    chr19_start <- sum(sapply(corr_list[1:18],
                              function(x) if (!is.null(x)) nrow(x) else 0L)) + 1L
    chr19_end   <- chr19_start + nrow(corr_list[[19]]) - 1L

    apoe_local <- apoe_idx_in_map[apoe_idx_in_map >= chr19_start &
                                  apoe_idx_in_map <= chr19_end] - chr19_start + 1L
    if (length(apoe_local) > 0) {
      keep_idx <- setdiff(1:nrow(corr_list[[19]]), apoe_local)
      corr_noAPOE[[19]] <- corr_list[[19]][keep_idx, keep_idx]
      cat(sprintf("  Chr19 LD: %d -> %d SNPs\n",
                  nrow(corr_list[[19]]), nrow(corr_noAPOE[[19]])))
    }
  }

  # 11.3 Re-run LDpred2-auto
  eoad_noAPOE <- run_ldpred2_auto_by_chr(
    df_beta_input = df_beta_eoad_noAPOE, corr_list = corr_noAPOE,
    map_ldref = map_ldref, trait_name = "EOAD_noAPOE"
  )

  load_noAPOE <- run_ldpred2_auto_by_chr(
    df_beta_input = df_beta_load_noAPOE, corr_list = corr_noAPOE,
    map_ldref = map_ldref, trait_name = "LOAD_noAPOE"
  )

  # 11.4 Comparison table
  comparison <- data.frame(
    Trait   = rep(c("EOAD", "LOAD"), each = 2),
    Version = rep(c("Full", "Excluding APOE"), 2),
    h2 = c(eoad_ldpred2$h2_bagged, eoad_noAPOE$h2_bagged,
           load_ldpred2$h2_bagged, load_noAPOE$h2_bagged),
    p  = c(eoad_ldpred2$p_bagged, eoad_noAPOE$p_bagged,
           load_ldpred2$p_bagged, load_noAPOE$p_bagged),
    stringsAsFactors = FALSE
  )

  # APOE contribution percentage
  eoad_contrib <- (eoad_ldpred2$h2_bagged - eoad_noAPOE$h2_bagged) /
                   eoad_ldpred2$h2_bagged * 100
  load_contrib <- (load_ldpred2$h2_bagged - load_noAPOE$h2_bagged) /
                   load_ldpred2$h2_bagged * 100

  comparison$APOE_Contribution_Pct <- c(
    sprintf("%.1f%%", eoad_contrib), NA,
    sprintf("%.1f%%", load_contrib), NA
  )

  cat("\n  Heritability comparison:\n")
  print(comparison)

  fwrite(comparison, file.path(output_dir, "Table_h2_Full_vs_noAPOE.csv"))
  cat("  Saved: Table_h2_Full_vs_noAPOE.csv\n")

  # 11.5 Visualization
  plot_apoe_sensitivity(comparison,
                        file.path(output_dir, "Figure_h2_APOE_Sensitivity.pdf"))

  return(list(
    eoad_noAPOE = eoad_noAPOE,
    load_noAPOE = load_noAPOE,
    comparison  = comparison
  ))
}


# ==============================================================================
# 12. MAIN EXECUTION WRAPPER
# ==============================================================================

# Orchestrates the full PRS analysis pipeline.
# All file paths are passed as parameters (no hardcoded paths).
run_prs_analysis <- function(eoad_gwas_file, load_gwas_file, aging_gwas_file,
                             ldref_dir, output_dir,
                             eoad_n_case = 1573, eoad_n_control = 199505,
                             load_n_case = 85934, load_n_control = 401577,
                             aging_n_total = 752566,
                             eoad_chr_col = "CHR", eoad_pos_col = "BP",
                             eoad_a1_col = "A1", eoad_a2_col = "A2",
                             eoad_beta_col = "BETA", eoad_se_col = "SE",
                             eoad_snp_col = "SNP",
                             load_chr_col = "CHR", load_pos_col = "BP",
                             load_a1_col = "A1", load_a2_col = "A2",
                             load_beta_col = "BETA", load_se_col = "SE",
                             load_snp_col = "SNP",
                             aging_chr_col = "CHR", aging_pos_col = "BP",
                             aging_a1_col = "A1", aging_a2_col = "A2",
                             aging_beta_col = "BETA", aging_se_col = "SE",
                             aging_snp_col = "SNP",
                             run_apoe_sensitivity_flag = TRUE) {

  figures_dir <- file.path(output_dir, "figures")
  tables_dir  <- file.path(output_dir, "tables")
  weights_dir <- file.path(output_dir, "weights")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(weights_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== PRS Analysis Pipeline ===\n\n")

  # Step 1: Format GWAS
  cat("Step 1: Formatting GWAS summary statistics\n")
  gwas_eoad <- format_gwas_ldpred2(eoad_gwas_file, n_case = eoad_n_case,
    n_control = eoad_n_control, trait_type = "binary",
    chr_col = eoad_chr_col, pos_col = eoad_pos_col,
    a1_col = eoad_a1_col, a2_col = eoad_a2_col,
    beta_col = eoad_beta_col, se_col = eoad_se_col, snp_col = eoad_snp_col)

  gwas_load <- format_gwas_ldpred2(load_gwas_file, n_case = load_n_case,
    n_control = load_n_control, trait_type = "binary",
    chr_col = load_chr_col, pos_col = load_pos_col,
    a1_col = load_a1_col, a2_col = load_a2_col,
    beta_col = load_beta_col, se_col = load_se_col, snp_col = load_snp_col)

  gwas_aging <- format_gwas_ldpred2(aging_gwas_file, n_total = aging_n_total,
    trait_type = "continuous",
    chr_col = aging_chr_col, pos_col = aging_pos_col,
    a1_col = aging_a1_col, a2_col = aging_a2_col,
    beta_col = aging_beta_col, se_col = aging_se_col, snp_col = aging_snp_col)

  # Step 2: Load LD reference
  cat("\nStep 2: Loading LD reference\n")
  ldref <- load_ld_reference(ldref_dir)
  map_ldref <- ldref$map
  corr      <- ldref$corr

  # Step 3: SNP matching
  cat("\nStep 3: Matching SNPs to LD reference\n")
  df_beta_eoad  <- match_snps_to_reference(gwas_eoad, map_ldref)
  df_beta_load  <- match_snps_to_reference(gwas_load, map_ldref)
  df_beta_aging <- match_snps_to_reference(gwas_aging, map_ldref)

  # Step 4: LDpred2-auto
  cat("\nStep 4: Running LDpred2-auto with chain bagging\n")
  eoad_ldpred2  <- run_ldpred2_auto_by_chr(df_beta_eoad, corr, map_ldref, "EOAD")
  load_ldpred2  <- run_ldpred2_auto_by_chr(df_beta_load, corr, map_ldref, "LOAD")
  aging_ldpred2 <- run_ldpred2_auto_by_chr(df_beta_aging, corr, map_ldref, "Aging")

  # Step 5: APOE-masked weights
  cat("\nStep 5: Generating APOE-masked weights\n")
  eoad_noAPOE_weights <- mask_apoe_weights(eoad_ldpred2$beta_bagged, df_beta_eoad)
  load_noAPOE_weights <- mask_apoe_weights(load_ldpred2$beta_bagged, df_beta_load)

  # Step 6-7: Pathway-specific PRS
  cat("\nStep 6-7: Generating pathway-specific PRS weights\n")
  pathway_genes <- define_pathway_genes()
  pw <- list()

  for (trait_label in c("EOAD", "LOAD")) {
    if (trait_label == "EOAD") {
      full_w <- eoad_ldpred2$beta_bagged; snp_df <- df_beta_eoad
    } else {
      full_w <- load_ldpred2$beta_bagged; snp_df <- df_beta_load
    }
    for (pw_name in names(pathway_genes)) {
      key <- paste0(tolower(trait_label), "_", pw_name)
      pw[[key]] <- generate_pathway_weights(full_w, snp_df,
                                            pathway_genes[[pw_name]],
                                            paste(trait_label, pw_name))
    }
  }

  # Step 8: Export weights
  cat("\nStep 8: Exporting PRS weights\n")
  export_weights(df_beta_eoad, eoad_ldpred2$beta_bagged,
                 file.path(weights_dir, "EOAD_PRS_weights_full.txt"))
  export_weights(df_beta_eoad, eoad_noAPOE_weights,
                 file.path(weights_dir, "EOAD_PRS_weights_noAPOE.txt"))
  export_weights(df_beta_load, load_ldpred2$beta_bagged,
                 file.path(weights_dir, "LOAD_PRS_weights_full.txt"))
  export_weights(df_beta_load, load_noAPOE_weights,
                 file.path(weights_dir, "LOAD_PRS_weights_noAPOE.txt"))
  export_weights(df_beta_aging, aging_ldpred2$beta_bagged,
                 file.path(weights_dir, "Aging_PRS_weights_full.txt"))

  for (trait_label in c("EOAD", "LOAD")) {
    snp_df <- if (trait_label == "EOAD") df_beta_eoad else df_beta_load
    for (pw_name in names(pathway_genes)) {
      key <- paste0(tolower(trait_label), "_", pw_name)
      export_weights(snp_df, pw[[key]]$weights,
                     file.path(weights_dir,
                               paste0(trait_label, "_PRS_weights_", pw_name, ".txt")))
    }
  }

  # Step 9: Visualization
  cat("\nStep 9: Generating visualizations\n")
  plot_genetic_architecture(
    list(eoad_ldpred2, load_ldpred2, aging_ldpred2),
    file.path(figures_dir, "Figure_Genetic_Architecture.pdf"))

  plot_cellular_origin(
    pw, list(eoad_ldpred2$beta_bagged, load_ldpred2$beta_bagged),
    names(pathway_genes), c("EOAD", "LOAD"),
    file.path(figures_dir, "Figure_Cellular_Origin_Burden.pdf"))

  # Step 10: Results tables
  cat("\nStep 10: Generating results tables\n")
  generate_results_tables(
    list(eoad_ldpred2, load_ldpred2, aging_ldpred2),
    list(df_beta_eoad, df_beta_load, df_beta_aging),
    pw, pathway_genes, tables_dir)

  # Step 11: APOE sensitivity (optional)
  apoe_result <- NULL
  if (run_apoe_sensitivity_flag) {
    cat("\nStep 11: APOE sensitivity analysis\n")
    apoe_result <- run_apoe_sensitivity(
      df_beta_eoad, df_beta_load, corr, map_ldref,
      eoad_ldpred2, load_ldpred2, output_dir)
  }

  # Save full environment
  saveRDS(list(
    eoad = eoad_ldpred2, load = load_ldpred2, aging = aging_ldpred2,
    pathway_weights = pw, pathway_genes = pathway_genes,
    apoe_sensitivity = apoe_result
  ), file.path(output_dir, "PRS_analysis_results.rds"))

  cat("\n=== PRS Analysis Complete ===\n")
  cat("Output directory:", output_dir, "\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  sessionInfo()
}


# ==============================================================================
# 13. EXAMPLE USAGE
# ==============================================================================

# run_prs_analysis(
#   eoad_gwas_file  = "data/EOAD_GWAS_hg19.txt",
#   load_gwas_file  = "data/LOAD_GWAS_hg19.txt",
#   aging_gwas_file = "data/Aging_GWAS_hg19.txt",
#   ldref_dir       = "data/ldref_hm3_plus",
#   output_dir      = "results/prs",
#   eoad_n_case     = 1573,
#   eoad_n_control  = 199505,
#   load_n_case     = 85934,
#   load_n_control  = 401577,
#   aging_n_total   = 752566
# )
