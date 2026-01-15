==============================================================================
  # EOAD vs LOAD Polygenic Risk Score Construction using LDpred2
  # ==============================================================================
#
# Analysis Contents:
#   Part 1: Environment Setup
#   Part 2: Load and Format GWAS Summary Statistics
#   Part 3: LD Reference Panel and SNP Matching
#   Part 4: LDpred2-auto with Chain Bagging
#   Part 5: APOE-Excluded PRS Weights
#   Part 6: Pathway Gene Set Definitions (Methods-Matched)
#   Part 7: Pathway-Specific PRS Weight Generation
#   Part 8: Export PRS Weights
#   Part 9: Visualization
#   Part 10: Results Summary
#
# ==============================================================================
# ==============================================================================
# Part 1: Environment Setup
# ==============================================================================

packages_cran <- c("data.table", "dplyr", "ggplot2", "cowplot", "Matrix",
                   "RColorBrewer", "scales", "gridExtra")
packages_bioc <- c("bigsnpr", "org.Hs.eg.db", "GenomicRanges",
                   "TxDb.Hsapiens.UCSC.hg19.knownGene", "AnnotationDbi")

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

dir.create("results/", showWarnings = FALSE)
dir.create("results/weights/", showWarnings = FALSE)
dir.create("results/figures/", showWarnings = FALSE)
dir.create("results/tables/", showWarnings = FALSE)

# ==============================================================================
# Part 2: Load and Format GWAS Summary Statistics
# ==============================================================================

## Section 2.1: Load GWAS Data
eoad_gwas <- fread("EOAD_sumstats_hg19.txt")
load_gwas <- fread("LOAD_sumstats_hg19.txt")
aging_gwas <- fread("Aging_sumstats_hg19.txt")

eoad_gwas <- eoad_gwas[!is.na(BP), ]
load_gwas <- load_gwas[!is.na(BP), ]
aging_gwas <- aging_gwas[!is.na(BP), ]

## Section 2.2: Format GWAS for LDpred2
format_gwas_ldpred2 <- function(gwas, n_case = NULL, n_control = NULL,
                                n_total = NULL, trait_type = "binary") {
  # Calculate effective sample size
  # N_eff = 4/(1/N_cases + 1/N_controls) for binary traits
  if (trait_type == "binary") {
    n_eff <- 4 / (1/n_case + 1/n_control)
  } else {
    n_eff <- n_total
  }
  
  gwas_std <- data.frame(
    chr = as.integer(gwas$CHR),
    pos = as.integer(gwas$BP),
    a0 = toupper(gwas$A2),
    a1 = toupper(gwas$A1),
    beta = as.numeric(gwas$BETA),
    beta_se = as.numeric(gwas$SE),
    rsid = gwas$SNP,
    n_eff = n_eff,
    stringsAsFactors = FALSE
  )
  
  # Quality control
  gwas_std <- gwas_std[complete.cases(gwas_std), ]
  gwas_std <- gwas_std[gwas_std$chr %in% 1:22, ]
  gwas_std <- gwas_std[nchar(gwas_std$a0) == 1 & nchar(gwas_std$a1) == 1, ]
  gwas_std <- gwas_std[!duplicated(paste(gwas_std$chr, gwas_std$pos)), ]
  
  return(gwas_std)
}

## Section 2.3: Apply Formatting
gwas_eoad <- format_gwas_ldpred2(eoad_gwas, n_case = 1573, n_control = 199505,
                                 trait_type = "binary")
gwas_load <- format_gwas_ldpred2(load_gwas, n_case = 85934, n_control = 401577,
                                 trait_type = "binary")
gwas_aging <- format_gwas_ldpred2(aging_gwas, n_total = 752566,
                                  trait_type = "continuous")

cat("GWAS formatted:\n")
cat("  EOAD:", nrow(gwas_eoad), "SNPs\n")
cat("  LOAD:", nrow(gwas_load), "SNPs\n")
cat("  Aging:", nrow(gwas_aging), "SNPs\n")

# ==============================================================================
# Part 3: LD Reference Panel and SNP Matching
# ==============================================================================

## Section 3.1: Load LD Reference
map_ldref <- readRDS("ldref/map.rds")
corr <- readRDS("ldref/LD_with_blocks.rds")

## Section 3.2: SNP Matching with HapMap3+ Reference
df_beta_eoad <- snp_match(gwas_eoad, map_ldref, strand_flip = TRUE,
                          join_by_pos = TRUE, remove_dups = TRUE)
df_beta_load <- snp_match(gwas_load, map_ldref, strand_flip = TRUE,
                          join_by_pos = TRUE, remove_dups = TRUE)
df_beta_aging <- snp_match(gwas_aging, map_ldref, strand_flip = TRUE,
                           join_by_pos = TRUE, remove_dups = TRUE)

cat("\nSNPs matched to LD reference:\n")
cat("  EOAD:", nrow(df_beta_eoad), "SNPs\n")
cat("  LOAD:", nrow(df_beta_load), "SNPs\n")
cat("  Aging:", nrow(df_beta_aging), "SNPs\n")

# ==============================================================================
# Part 4: LDpred2-auto with Chain Bagging
# ==============================================================================

## Section 4.1: LDpred2-auto Function
run_ldpred2_auto <- function(df_beta, corr, trait_name,
                             h2_init = 0.1, n_chains = 30,
                             burn_in = 500, num_iter = 500) {
  
  # Initialize polygenicity across range 10^-4 to 0.5
  vec_p_init <- seq_log(1e-4, 0.5, length.out = n_chains)
  
  # Run LDpred2-auto with sparse option
  multi_auto <- snp_ldpred2_auto(
    corr = corr,
    df_beta = df_beta,
    h2_init = h2_init,
    vec_p_init = vec_p_init,
    burn_in = burn_in,
    num_iter = num_iter,
    sparse = TRUE,
    allow_jump_sign = TRUE,
    shrink_corr = 0.95,
    ncores = 4
  )
  
  # Extract estimates from converged chains
  h2_est <- sapply(multi_auto, function(x) x$h2_est)
  p_est <- sapply(multi_auto, function(x) x$p_est)
  
  # Filter converged chains
  keep <- (h2_est > 1e-6) & (h2_est < 1) & (p_est > 1e-6)
  if (sum(keep) == 0) keep <- rep(TRUE, length(h2_est))
  
  # Chain bagging: average posterior estimates across converged chains
  beta_list <- lapply(multi_auto[keep], function(x) x$beta_est)
  beta_bagged <- Reduce("+", beta_list) / length(beta_list)
  
  return(list(
    beta_bagged = beta_bagged,
    h2_bagged = mean(h2_est[keep]),
    p_bagged = mean(p_est[keep]),
    n_converged = sum(keep),
    trait_name = trait_name
  ))
}

## Section 4.2: Run LDpred2-auto for All Traits
cat("\nRunning LDpred2-auto (30 chains, 500 burn-in, 500 iterations)...\n")

eoad_ldpred2 <- run_ldpred2_auto(df_beta_eoad, corr, "EOAD")
cat("  EOAD: h2 =", round(eoad_ldpred2$h2_bagged, 4), 
    ", p =", format(eoad_ldpred2$p_bagged, scientific = TRUE),
    ", converged =", eoad_ldpred2$n_converged, "/30\n")

load_ldpred2 <- run_ldpred2_auto(df_beta_load, corr, "LOAD")
cat("  LOAD: h2 =", round(load_ldpred2$h2_bagged, 4),
    ", p =", format(load_ldpred2$p_bagged, scientific = TRUE),
    ", converged =", load_ldpred2$n_converged, "/30\n")

aging_ldpred2 <- run_ldpred2_auto(df_beta_aging, corr, "Aging")
cat("  Aging: h2 =", round(aging_ldpred2$h2_bagged, 4),
    ", p =", format(aging_ldpred2$p_bagged, scientific = TRUE),
    ", converged =", aging_ldpred2$n_converged, "/30\n")

## Section 4.3: Save LDpred2 Results
saveRDS(eoad_ldpred2, "results/EOAD_LDpred2_weights.rds")
saveRDS(load_ldpred2, "results/LOAD_LDpred2_weights.rds")
saveRDS(aging_ldpred2, "results/Aging_LDpred2_weights.rds")

# ==============================================================================
# Part 5: APOE-Excluded PRS Weights
# ==============================================================================

## Section 5.1: Define APOE Region (chr19: 44-46 Mb)
apoe_snps_eoad <- which(df_beta_eoad$chr == 19 &
                          df_beta_eoad$pos >= 44000000 &
                          df_beta_eoad$pos <= 46000000)

apoe_snps_load <- which(df_beta_load$chr == 19 &
                          df_beta_load$pos >= 44000000 &
                          df_beta_load$pos <= 46000000)

apoe_snps_aging <- which(df_beta_aging$chr == 19 &
                           df_beta_aging$pos >= 44000000 &
                           df_beta_aging$pos <= 46000000)

## Section 5.2: Generate APOE-Excluded Weights
eoad_weights_noAPOE <- eoad_ldpred2$beta_bagged
eoad_weights_noAPOE[apoe_snps_eoad] <- 0

load_weights_noAPOE <- load_ldpred2$beta_bagged
load_weights_noAPOE[apoe_snps_load] <- 0

aging_weights_noAPOE <- aging_ldpred2$beta_bagged
aging_weights_noAPOE[apoe_snps_aging] <- 0

cat("\nAPOE region SNPs masked:\n")
cat("  EOAD:", length(apoe_snps_eoad), "SNPs\n")
cat("  LOAD:", length(apoe_snps_load), "SNPs\n")
cat("  Aging:", length(apoe_snps_aging), "SNPs\n")

# ==============================================================================
# Part 6: Pathway Gene Set Definitions (Methods-Matched)
# ==============================================================================

## Section 6.1: Define Six Pathway Gene Sets
## Gene counts matched to supplementary methods:
## T-cell (56), Microglia (58), Abeta (52), APP (38), Oligo (62), Myelin (52)

define_pathway_genes <- function() {
  list(
    # T-cell activation (56 genes)
    TCell = c(
      "CD3D", "CD3E", "CD3G", "CD247", "CD2", "CD5", "CD7", "CD27", "CD28",
      "TRAC", "TRBC1", "TRBC2", "TCF7", "LEF1", "GATA3", "TBX21", "EOMES",
      "LCK", "ZAP70", "LAT", "ITK", "FYN", "PLCG1", "VAV1", "PRKCQ", "CARD11",
      "IL7R", "IL2RA", "IL2RB", "IL2RG", "IFNG", "TNF", "IL17A",
      "CCR7", "CCR4", "CCR5", "CXCR3", "CXCR4", "CXCR5",
      "CTLA4", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "ICOS", "CD40LG",
      "GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7", "FASLG",
      "NFATC1", "NFATC2", "NFKB1"
    ),
    
    # Activated microglia (58 genes)
    Microglia = c(
      "TREM2", "TYROBP", "APOE", "LPL", "CST7", "ITGAX", "CLEC7A", "SPP1",
      "GPNMB", "LGALS3", "CD9", "FABP5", "CTSB", "CTSD", "CTSL", "CTSS",
      "C1QA", "C1QB", "C1QC", "C3", "C3AR1", "C5AR1", "ITGAM",
      "IL1B", "IL6", "TNF", "CCL2", "CCL3", "CCL4", "CXCL10", "CXCL16",
      "CD68", "CD14", "FCGR1A", "FCGR2A", "FCGR3A", "MSR1", "CD36", "MARCO",
      "MS4A4A", "MS4A6A", "INPP5D", "CD33", "PLCG2", "ABI3", "GRN",
      "SPI1", "IRF8", "RUNX1", "AIF1", "CSF1R", "CX3CR1",
      "P2RY12", "TMEM119", "HEXB", "SALL1", "MERTK", "GAS6"
    ),
    
    # Amyloid-beta clearance (52 genes)
    Abeta = c(
      "IDE", "MME", "ECE1", "ECE2", "ACE", "THOP1", "MMP2", "MMP9", "MMP14",
      "LRP1", "LRP2", "LDLR", "VLDLR", "SORL1", "SCARB1", "SCARB2",
      "ABCA1", "ABCA7", "ABCB1", "ABCG1", "ABCG2", "ABCG4",
      "TREM2", "CD36", "SCARA1", "MSR1", "CD14", "TLR2", "TLR4", "TLR6",
      "BECN1", "ATG5", "ATG7", "ATG12", "SQSTM1", "TFEB", "LAMP1", "LAMP2",
      "HSPA1A", "HSPA8", "HSP90AA1", "DNAJB1", "CLU", "PICALM", "BIN1",
      "CD2AP", "EPHA1", "PTK2B", "CASS4", "FERMT2", "SLC24A4", "ZCWPW1"
    ),
    
    # APP metabolism (38 genes)
    APP = c(
      "PSEN1", "PSEN2", "NCSTN", "APH1A", "APH1B", "PSENEN", "PEN2",
      "BACE1", "BACE2", "ADAM10", "ADAM17", "ADAM9", "ADAM12",
      "APP", "APLP1", "APLP2", "ITM2B", "ITM2C",
      "SORL1", "SORCS1", "SORCS2", "SORCS3", "LRP1", "LRP8",
      "BIN1", "PICALM", "CD2AP", "RIN3", "SH3KBP1",
      "RAB5A", "RAB7A", "RAB11A", "RAB11B", "VPS35", "VPS26A", "VPS29",
      "SNX27", "APBB1"
    ),
    
    # Oligodendrocyte function (62 genes)
    Oligo = c(
      "MBP", "PLP1", "MOG", "MAG", "MOBP", "CLDN11", "CNP", "MYRF", "OPALIN",
      "FA2H", "UGT8", "GAL3ST1", "GALC", "ASPA", "ENPP6", "BCAS1",
      "OLIG1", "OLIG2", "SOX10", "NKX2-2", "NKX6-2", "ZNF488",
      "PDGFRA", "CSPG4", "GPR17", "PTPRZ1", "GPR37", "GPR37L1",
      "ERMN", "NFASC", "CNTN2", "LPAR1", "S1PR5",
      "ERBB3", "ERBB4", "NRG1", "NRG2", "LINGO1", "RTN4R",
      "ABCA2", "ABCD1", "PEX1", "PEX5", "PEX7", "PEX10", "PEX13",
      "HMGCR", "HMGCS1", "MVK", "FDPS", "FDFT1", "SQLE", "CYP51A1",
      "DHCR7", "DHCR24", "SC5D", "EBP", "NSDHL", "HSD17B7",
      "QKI", "SIRT2"
    ),
    
    # Myelination (52 genes)
    Myelin = c(
      "MBP", "PLP1", "MOG", "MAG", "MOBP", "CNP", "MYRF", "OPALIN",
      "FA2H", "UGT8", "GAL3ST1", "GALC", "ASPA", "GJC2", "GJB1",
      "NRG1", "NRG2", "ERBB2", "ERBB3", "ERBB4", "LINGO1", "RTN4", "RTN4R",
      "PMP22", "MPZ", "PRX", "EGR2", "SOX10", "KROX20",
      "ABCD1", "ABCD2", "PEX1", "PEX5", "PEX7", "PEX10",
      "HMGCR", "HMGCS1", "MVK", "FDPS", "FDFT1", "SQLE", "CYP51A1",
      "DHCR7", "DHCR24", "SC5D", "EBP", "NSDHL",
      "CERS2", "SPTLC1", "SPTLC2", "SGMS1", "UGCG"
    )
  )
}

pathway_gene_lists <- define_pathway_genes()

cat("\nPathway gene sets defined:\n")
for (pw in names(pathway_gene_lists)) {
  cat("  ", pw, ":", length(pathway_gene_lists[[pw]]), "genes\n")
}

# ==============================================================================
# Part 7: Pathway-Specific PRS Weight Generation
# ==============================================================================

## Section 7.1: Function to Get SNPs Within Gene Window
get_snps_for_genes <- function(target_genes, snp_info, window = 100000) {
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  gene_coords <- genes(txdb)
  
  # Map Entrez IDs to gene symbols
  map_sym <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = names(gene_coords),
                                   columns = "SYMBOL", 
                                   keytype = "ENTREZID")
  mcols(gene_coords)$SYMBOL <- map_sym$SYMBOL[match(names(gene_coords),
                                                    map_sym$ENTREZID)]
  
  # Filter to target genes
  target_gr <- gene_coords[mcols(gene_coords)$SYMBOL %in% target_genes &
                             !is.na(mcols(gene_coords)$SYMBOL)]
  
  if (length(target_gr) == 0) return(integer(0))
  
  # Extend gene boundaries by window (100kb upstream and downstream)
  target_gr <- resize(target_gr, width(target_gr) + 2 * window, fix = "center")
  
  # Create SNP GRanges
  snp_gr <- GRanges(
    seqnames = paste0("chr", snp_info$chr),
    ranges = IRanges(start = snp_info$pos, width = 1)
  )
  
  # Find overlaps
  overlaps <- findOverlaps(snp_gr, target_gr)
  return(unique(queryHits(overlaps)))
}

## Section 7.2: Generate Pathway-Specific Weights
generate_pathway_weights <- function(full_weights, snp_info, gene_list) {
  idx <- get_snps_for_genes(gene_list, snp_info)
  
  if (length(idx) == 0) {
    return(list(
      weights = rep(0, length(full_weights)),
      n_snps = 0,
      n_nonzero = 0
    ))
  }
  
  # Mask non-pathway SNPs to zero, retain LDpred2 weights for pathway SNPs
  new_weights <- rep(0, length(full_weights))
  new_weights[idx] <- full_weights[idx]
  
  return(list(
    weights = new_weights,
    n_snps = length(idx),
    n_nonzero = sum(new_weights != 0)
  ))
}

## Section 7.3: Generate All Pathway Weights
pathway_weights <- list()

for (pathway in names(pathway_gene_lists)) {
  # EOAD pathway weights
  pathway_weights[[paste0("eoad_", pathway)]] <- 
    generate_pathway_weights(eoad_ldpred2$beta_bagged, df_beta_eoad,
                             pathway_gene_lists[[pathway]])
  
  # LOAD pathway weights
  pathway_weights[[paste0("load_", pathway)]] <- 
    generate_pathway_weights(load_ldpred2$beta_bagged, df_beta_load,
                             pathway_gene_lists[[pathway]])
}

cat("\nPathway-specific PRS weights generated:\n")
for (pw in names(pathway_gene_lists)) {
  cat("  ", pw, ": EOAD=", pathway_weights[[paste0("eoad_", pw)]]$n_snps,
      " SNPs, LOAD=", pathway_weights[[paste0("load_", pw)]]$n_snps, " SNPs\n", sep = "")
}

# ==============================================================================
# Part 8: Export PRS Weights
# ==============================================================================

## Section 8.1: Export Function
export_weights <- function(snp_info, weights, filename) {
  df <- data.frame(
    SNP = snp_info$rsid,
    CHR = snp_info$chr,
    POS = snp_info$pos,
    A1 = snp_info$a1,
    A2 = snp_info$a0,
    BETA = weights
  )
  fwrite(df, filename, sep = "\t")
}

## Section 8.2: Export Full PRS Weights
export_weights(df_beta_eoad, eoad_ldpred2$beta_bagged,
               "results/weights/EOAD_PRS_weights_full.txt")
export_weights(df_beta_eoad, eoad_weights_noAPOE,
               "results/weights/EOAD_PRS_weights_noAPOE.txt")

export_weights(df_beta_load, load_ldpred2$beta_bagged,
               "results/weights/LOAD_PRS_weights_full.txt")
export_weights(df_beta_load, load_weights_noAPOE,
               "results/weights/LOAD_PRS_weights_noAPOE.txt")

export_weights(df_beta_aging, aging_ldpred2$beta_bagged,
               "results/weights/Aging_PRS_weights_full.txt")
export_weights(df_beta_aging, aging_weights_noAPOE,
               "results/weights/Aging_PRS_weights_noAPOE.txt")

## Section 8.3: Export Pathway-Specific Weights
for (name in names(pathway_weights)) {
  if (grepl("^eoad_", name)) {
    export_weights(df_beta_eoad, pathway_weights[[name]]$weights,
                   paste0("results/weights/", toupper(name), "_PRS_weights.txt"))
  } else {
    export_weights(df_beta_load, pathway_weights[[name]]$weights,
                   paste0("results/weights/", toupper(name), "_PRS_weights.txt"))
  }
}

cat("\nPRS weights exported to results/weights/\n")

# ==============================================================================
# Part 9: Visualization
# ==============================================================================

theme_nc <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

## Section 9.1: Genetic Architecture Plot
architecture_data <- data.frame(
  Trait = c("EOAD", "LOAD", "Aging"),
  h2 = c(eoad_ldpred2$h2_bagged, load_ldpred2$h2_bagged, aging_ldpred2$h2_bagged),
  p = c(eoad_ldpred2$p_bagged, load_ldpred2$p_bagged, aging_ldpred2$p_bagged),
  N_NonZero = c(sum(eoad_ldpred2$beta_bagged != 0),
                sum(load_ldpred2$beta_bagged != 0),
                sum(aging_ldpred2$beta_bagged != 0))
)

p_architecture <- ggplot(architecture_data, 
                         aes(x = p, y = h2, color = Trait, size = N_NonZero)) +
  geom_point(alpha = 0.9) +
  geom_text(aes(label = Trait), vjust = -1.5, size = 5, fontface = "bold",
            show.legend = FALSE) +
  scale_x_log10(labels = scales::scientific) +
  scale_color_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5", 
                                "Aging" = "#00A087")) +
  scale_size_continuous(range = c(8, 15), name = "Non-zero SNPs",
                        labels = scales::comma) +
  labs(title = "Genetic Architecture: EOAD vs LOAD vs Aging",
       subtitle = "LDpred2-auto estimates",
       x = "Polygenicity (p, log scale)",
       y = "SNP Heritability (h²)") +
  theme_nc

ggsave("results/figures/Figure_Genetic_Architecture.pdf", p_architecture, 
       width = 10, height = 8, dpi = 300)

## Section 9.2: Cellular Origin Map
calculate_pathway_burden <- function(pw, fw) {
  ps <- sum(abs(pw), na.rm = TRUE)
  fs <- sum(abs(fw), na.rm = TRUE)
  if (fs == 0) return(0)
  return(ps / fs * 100)
}

cellular_origin_data <- data.frame(
  GWAS = rep(c("EOAD", "LOAD"), each = 6),
  Pathway = factor(rep(c("T-Cell", "Microglia", "Aβ Clearance", 
                         "APP Metabolism", "Oligodendrocyte", "Myelination"), 2),
                   levels = c("T-Cell", "Microglia", "Aβ Clearance", 
                              "APP Metabolism", "Oligodendrocyte", "Myelination")),
  Burden = c(
    calculate_pathway_burden(pathway_weights$eoad_TCell$weights, 
                             eoad_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$eoad_Microglia$weights, 
                             eoad_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$eoad_Abeta$weights, 
                             eoad_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$eoad_APP$weights, 
                             eoad_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$eoad_Oligo$weights, 
                             eoad_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$eoad_Myelin$weights, 
                             eoad_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$load_TCell$weights, 
                             load_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$load_Microglia$weights, 
                             load_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$load_Abeta$weights, 
                             load_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$load_APP$weights, 
                             load_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$load_Oligo$weights, 
                             load_ldpred2$beta_bagged),
    calculate_pathway_burden(pathway_weights$load_Myelin$weights, 
                             load_ldpred2$beta_bagged)
  )
)

p_cellular_origin <- ggplot(cellular_origin_data, 
                            aes(x = Pathway, y = Burden, fill = GWAS)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", Burden)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5")) +
  labs(title = "Cellular Origin Map: EOAD vs LOAD",
       subtitle = "Pathway-specific genetic burden",
       x = "", y = "Pathway Genetic Burden (%)") +
  theme_nc +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

ggsave("results/figures/Figure_Cellular_Origin_Map.pdf", p_cellular_origin, 
       width = 12, height = 8, dpi = 300)

## Section 9.3: Pathway SNP Distribution
pathway_snp_data <- data.frame(
  GWAS = rep(c("EOAD", "LOAD"), each = 6),
  Pathway = factor(rep(names(pathway_gene_lists), 2),
                   levels = names(pathway_gene_lists)),
  N_SNPs = c(
    sapply(names(pathway_gene_lists), function(p) 
      pathway_weights[[paste0("eoad_", p)]]$n_snps),
    sapply(names(pathway_gene_lists), function(p) 
      pathway_weights[[paste0("load_", p)]]$n_snps)
  ),
  N_NonZero = c(
    sapply(names(pathway_gene_lists), function(p) 
      pathway_weights[[paste0("eoad_", p)]]$n_nonzero),
    sapply(names(pathway_gene_lists), function(p) 
      pathway_weights[[paste0("load_", p)]]$n_nonzero)
  )
)

p_snp_dist <- ggplot(pathway_snp_data, aes(x = Pathway, y = N_SNPs, fill = GWAS)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5")) +
  labs(title = "Pathway SNP Coverage",
       subtitle = "SNPs within 100kb of pathway genes",
       x = "", y = "Number of SNPs") +
  theme_nc +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

ggsave("results/figures/Figure_Pathway_SNP_Distribution.pdf", p_snp_dist, 
       width = 10, height = 6, dpi = 300)

# ==============================================================================
# Part 10: Results Summary
# ==============================================================================

## Section 10.1: LDpred2 Results Summary
results_summary <- data.frame(
  Trait = c("EOAD", "LOAD", "Aging"),
  GWAS_Source = c("FinnGen R11", "Bellenguez 2022", "Timmers 2020"),
  N_Cases = c(1573, 85934, NA),
  N_Controls = c(199505, 401577, NA),
  N_Total = c(201078, 487511, 752566),
  N_SNPs_Matched = c(nrow(df_beta_eoad), nrow(df_beta_load), nrow(df_beta_aging)),
  h2_LDpred2 = c(eoad_ldpred2$h2_bagged, load_ldpred2$h2_bagged, 
                 aging_ldpred2$h2_bagged),
  p_LDpred2 = c(eoad_ldpred2$p_bagged, load_ldpred2$p_bagged, 
                aging_ldpred2$p_bagged),
  N_Chains_Converged = c(eoad_ldpred2$n_converged, load_ldpred2$n_converged, 
                         aging_ldpred2$n_converged)
)

fwrite(results_summary, "results/tables/Table_LDpred2_Results_Summary.csv")

## Section 10.2: Pathway Summary
pathway_summary <- data.frame(
  GWAS = rep(c("EOAD", "LOAD"), each = length(pathway_gene_lists)),
  Pathway = rep(names(pathway_gene_lists), 2),
  N_Genes = rep(sapply(pathway_gene_lists, length), 2),
  N_SNPs = c(
    sapply(names(pathway_gene_lists), function(p) 
      pathway_weights[[paste0("eoad_", p)]]$n_snps),
    sapply(names(pathway_gene_lists), function(p) 
      pathway_weights[[paste0("load_", p)]]$n_snps)
  ),
  N_NonZero = c(
    sapply(names(pathway_gene_lists), function(p) 
      pathway_weights[[paste0("eoad_", p)]]$n_nonzero),
    sapply(names(pathway_gene_lists), function(p) 
      pathway_weights[[paste0("load_", p)]]$n_nonzero)
  )
)

fwrite(pathway_summary, "results/tables/Table_Pathway_Specific_PRS_Summary.csv")
fwrite(cellular_origin_data, "results/tables/Table_Cellular_Origin_Burden.csv")

## Section 10.3: Print Summary
cat("\n")
cat("================================================================================\n")
cat("  PRS Construction using LDpred2-auto - COMPLETE\n")
cat("================================================================================\n")
cat("\nLDpred2 Results:\n")
print(results_summary[, c("Trait", "N_SNPs_Matched", "h2_LDpred2", "p_LDpred2", 
                          "N_Chains_Converged")])
cat("\nPathway Gene Counts:\n")
for (pw in names(pathway_gene_lists)) {
  cat("  ", pw, ":", length(pathway_gene_lists[[pw]]), "genes\n")
}
cat("\nOutput Files:\n")
cat("  Weights:\n")
cat("    - results/weights/EOAD_PRS_weights_full.txt\n")
cat("    - results/weights/EOAD_PRS_weights_noAPOE.txt\n")
cat("    - results/weights/LOAD_PRS_weights_full.txt\n")
cat("    - results/weights/LOAD_PRS_weights_noAPOE.txt\n")
cat("    - results/weights/Aging_PRS_weights_full.txt\n")
cat("    - results/weights/[PATHWAY]_PRS_weights.txt (12 files)\n")
cat("  Tables:\n")
cat("    - results/tables/Table_LDpred2_Results_Summary.csv\n")
cat("    - results/tables/Table_Pathway_Specific_PRS_Summary.csv\n")
cat("    - results/tables/Table_Cellular_Origin_Burden.csv\n")
cat("  Figures:\n")
cat("    - results/figures/Figure_Genetic_Architecture.pdf\n")
cat("    - results/figures/Figure_Cellular_Origin_Map.pdf\n")
cat("    - results/figures/Figure_Pathway_SNP_Distribution.pdf\n")
cat("\n")

sessionInfo()



