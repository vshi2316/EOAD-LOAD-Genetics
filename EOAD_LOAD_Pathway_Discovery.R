# ==============================================================================
# EOAD vs LOAD: Data-Driven Pathway Discovery and Cell Type Specificity
# ==============================================================================
#
# Purpose: Identify divergent biological pathways between early-onset (EOAD)
#          and late-onset (LOAD) Alzheimer's disease using MAGMA gene-based
#          association, curated module enrichment, GSEA, gene-property
#          analysis, conditional regression, and downsampling validation.
#
# Analyses:
#   1. MAGMA gene-based association (read pre-computed .genes.out)
#   2. Graham et al. module enrichment (5 curated modules, Fisher's exact)
#   3. Gene set enrichment analysis (GSEA GO-BP, ranked by MAGMA Z-stats)
#   4. MAGMA gene-property analysis (9 cell type marker sets, linear model)
#   5. Conditional regression (T cell independence from microglia)
#   6. LOAD downsampling validation (deterministic SE scaling)
#   7. Publication-quality visualization
#
# Statistical notes:
#   - ZSTAT = qnorm(1 - P/2) with sign correction for P > 0.5
#   - Gene-property model: Z_MAGMA ~ beta_0 + beta_1 * CellType
#   - Conditional model: Z_MAGMA ~ Target + Microglia_Covariate
#   - Downsampling: SE_new = SE_old / sqrt(N_eff_EOAD / N_eff_LOAD)
#   - GSEA parameters: minGSSize=10, maxGSSize=500, pvalueCutoff=0.25
#   - Multiple testing: Bonferroni for MAGMA genes, BH-FDR for modules/GSEA
#
# Dependencies: data.table, dplyr, tidyr, ggplot2, clusterProfiler,
#   org.Hs.eg.db, enrichplot, ComplexHeatmap, circlize, RColorBrewer,
#   ggrepel, cowplot, scales
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(ggrepel)
  library(cowplot)
  library(scales)
})


# ==============================================================================
# Publication theme
# ==============================================================================

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
# 1. MAGMA Gene-Based Association
# ==============================================================================
# Reads pre-computed MAGMA .genes.out files.
# Expects columns: GENE (or SYMBOL), P (or PVALUE), NSNPS, NPARAM, ZSTAT.
# Computes FDR, Bonferroni, and signed Z-statistics.
# ==============================================================================

read_magma_results <- function(magma_file) {
  magma <- fread(magma_file)
  magma$SYMBOL    <- toupper(magma$GENE)
  magma$PVALUE    <- magma$P
  magma$FDR       <- p.adjust(magma$PVALUE, method = "BH")
  magma$Bonferroni <- p.adjust(magma$PVALUE, method = "bonferroni")
  magma$ZSTAT     <- qnorm(1 - magma$PVALUE / 2)
  magma$ZSTAT[magma$PVALUE > 0.5] <- -magma$ZSTAT[magma$PVALUE > 0.5]
  return(magma)
}


# ==============================================================================
# 2. Graham et al. Module Enrichment (Fisher's Exact Test)
# ==============================================================================
# Five curated modules from Graham et al.:
#   Human Microglial AD-Significant (26 genes)
#   Human Oligodendrocyte AD-Significant (29 genes)
#   Mouse Activated Response Microglia (17 genes)
#   Mouse Phagolysosomal (45 genes)
#   Mouse HM2 Longevity-Significant (38 genes)
# ==============================================================================

get_graham_modules <- function() {
  list(
    Human_Microglial_AD = c(
      "MS4A4A", "MS4A6A", "SPI1", "APOC2", "TREM2", "HLA-DRB1", "CD33",
      "HLA-DRB5", "INPP5D", "ARHGAP45", "HLA-DQB1", "LAPTM5", "HLA-DQA1",
      "HLA-DRA", "NOP2", "ATP8B4", "GPSM3", "GAL3ST4", "CMTM7", "PCED1B",
      "ITGAM", "DOK3", "TMC8", "MARCO", "COX7A1", "LILRB4"
    ),
    Human_Oligodendrocyte_AD = c(
      "CLASRP", "HNRNPA2B1", "TEX22", "ZKSCAN1", "SHC4", "BBIP1", "FOLH1",
      "SLC45A3", "CDR2L", "FAM76B", "ZFPL1", "PDCD4", "TMEM37", "SLC20A2",
      "PVRIG", "CENPC", "TMED7", "SLC26A9", "STAG3", "ACP7", "SNX1",
      "TMEM42", "TARDBP", "UNC5CL", "SLF2", "METTL14", "FAM222A", "SLC4A9",
      "UBA6"
    ),
    Mouse_ARM = c(
      "APOE", "RELB", "MS4A6A", "PVR", "PTK2B", "HLA-DQB1", "PLEKHA1",
      "YDJC", "LILRA5", "TRIM37", "NRP1", "STBD1", "ALKBH2", "ANKH",
      "TNIP2", "WDR55", "LILRB4"
    ),
    Mouse_Phagolysosomal = c(
      "TOMM40", "CLPTM1", "MS4A6A", "SPI1", "PSMC3", "TREM2", "FBXO46",
      "HBEGF", "GEMIN7", "ISYNA1", "PFDN1", "SSBP4", "MTCH2", "CNPY4",
      "FCF1", "NDUFAF6", "LILRA5", "DLST", "ACTB", "CSNK2B", "ARPC1A",
      "ZYX", "PRKRA", "TMEM37", "ZNF655", "PKP4", "PSMC6", "FAM220A",
      "SNX1", "TMEM42", "BNIP3L", "DDAH2", "MRPL43", "EIF1B", "KIAA1671",
      "ZMAT2", "POLR2E", "ECH1", "IK", "LILRB4", "OSER1", "ETF1", "SKP1",
      "MRPL58", "GRN"
    ),
    Mouse_HM2_Longevity = c(
      "FES", "ZKSCAN5", "PTCD1", "JAM3", "CASP8", "MAFK",
      "MCRS1", "GOLPH3L", "RBM4", "PPP2R3C", "BCAP29", "RIT1", "SRP54",
      "C2ORF69", "CLN8", "PLPP3", "MAIP1", "BUD31", "STAT3", "SCAMP5",
      "KLC2", "PTPN1", "ZNF394", "ZNF664", "EXOC3", "DUSP6", "SLC44A2",
      "XPC", "TAP2", "CTSF", "CSNK1A1", "GSTZ1", "FADS1", "USP38", "F11R",
      "BSN", "ABHD16A", "ZSCAN29"
    )
  )
}

test_module_enrichment <- function(magma_results, module_genes, module_name,
                                   p_threshold = 0.05) {
  module_genes <- toupper(module_genes)
  all_genes    <- unique(magma_results$SYMBOL)
  sig_genes    <- unique(magma_results$SYMBOL[magma_results$PVALUE < p_threshold])
  module_in_bg <- intersect(module_genes, all_genes)

  if (length(module_in_bg) < 5) {
    return(data.frame(Module = module_name, Module_Size = length(module_genes),
                      Overlap = NA, OR = NA, P = NA, stringsAsFactors = FALSE))
  }

  overlap <- intersect(module_in_bg, sig_genes)
  a  <- length(overlap)
  b  <- length(sig_genes) - a
  cc <- length(module_in_bg) - a
  d  <- length(all_genes) - a - b - cc

  ft <- fisher.test(matrix(c(a, cc, b, d), nrow = 2), alternative = "greater")

  data.frame(
    Module           = module_name,
    Module_Size      = length(module_genes),
    Genes_In_Background = length(module_in_bg),
    Overlap          = a,
    OR               = as.numeric(ft$estimate),
    CI_Lower         = ft$conf.int[1],
    CI_Upper         = ft$conf.int[2],
    P                = ft$p.value,
    Overlap_Genes    = paste(overlap, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

run_graham_enrichment <- function(eoad_magma, load_magma) {
  graham_modules <- get_graham_modules()
  results <- data.frame()

  for (mod_name in names(graham_modules)) {
    for (gwas_name in c("EOAD", "LOAD")) {
      magma_df <- if (gwas_name == "EOAD") eoad_magma else load_magma
      res <- test_module_enrichment(magma_df, graham_modules[[mod_name]], mod_name)
      res$GWAS <- gwas_name
      results <- rbind(results, res)
    }
  }

  results$FDR <- p.adjust(results$P, method = "BH")
  results$Significant <- results$FDR < 0.05

  cat("\nGraham Module Enrichment:\n")
  for (g in c("EOAD", "LOAD")) {
    cat(sprintf("  %s:\n", g))
    sub <- results[results$GWAS == g, ]
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("    %-30s OR=%.2f FDR=%.2e %s\n",
                  sub$Module[i], sub$OR[i], sub$FDR[i],
                  ifelse(sub$FDR[i] < 0.05, "*", "")))
    }
  }

  return(results)
}


# ==============================================================================
# 3. Gene Set Enrichment Analysis (GSEA GO-BP)
# ==============================================================================
# Genes ranked by MAGMA Z-statistics; GO Biological Process only.
# Parameters: minGSSize=10, maxGSSize=500, pvalueCutoff=0.25, FDR < 0.05
# ==============================================================================

run_gsea_gobp <- function(magma_df, gwas_label) {
  gene_list <- magma_df$ZSTAT
  names(gene_list) <- magma_df$SYMBOL
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  id_map <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
  gene_list_entrez <- gene_list[names(gene_list) %in% id_map$SYMBOL]
  names(gene_list_entrez) <- id_map$ENTREZID[match(names(gene_list_entrez),
                                                    id_map$SYMBOL)]
  gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
  gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

  prop_mapped <- length(gene_list_entrez) / length(gene_list)
  cat(sprintf("  %s gene mapping: %.1f%% (%d / %d)\n",
              gwas_label, 100 * prop_mapped,
              length(gene_list_entrez), length(gene_list)))

  gsea_res <- gseGO(geneList     = gene_list_entrez,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "BP",
                    minGSSize    = 10,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.25,
                    pAdjustMethod = "BH",
                    verbose      = FALSE)

  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    res_df <- gsea_res@result
    res_df$GWAS <- gwas_label
    cat(sprintf("  %s GSEA GO-BP: %d pathways (FDR < 0.05: %d)\n",
                gwas_label, nrow(res_df), sum(res_df$p.adjust < 0.05)))
    return(list(result = gsea_res, table = res_df))
  }
  return(NULL)
}


# ==============================================================================
# 4. MAGMA Gene-Property Analysis (9 Cell Type Marker Sets)
# ==============================================================================
# Linear model: Z_MAGMA = beta_0 + beta_1 * CellType_indicator
# Positive beta_1 indicates genetic association enrichment for that cell type.
# ==============================================================================

get_cell_type_markers <- function() {
  list(
    T_Cell_Pan = c(
      "CD3D", "CD3E", "CD3G", "CD247", "CD2", "CD5", "CD7", "CD27", "CD28",
      "TRAC", "TRBC1", "TRBC2", "TRAT1",
      "TCF7", "LEF1", "GATA3", "TBX21", "EOMES", "BCL11B", "RUNX3",
      "IKZF1", "IKZF3",
      "LCK", "ZAP70", "LAT", "ITK", "FYN", "PLCG1", "VAV1", "PRKCQ",
      "CARD11", "NFATC1", "NFATC2", "NFKB1", "NFKBIA", "REL",
      "IL7R", "IL2RA", "IL2RB", "IL2RG", "IL15RA", "IL21R", "IL23R",
      "IFNGR1",
      "CCR7", "CCR4", "CCR5", "CCR6", "CXCR3", "CXCR4", "CXCR5", "CXCR6",
      "CTLA4", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "ICOS", "CD40LG",
      "GZMA", "GZMB", "GZMK", "PRF1", "GNLY"
    ),
    T_Cell_CD4 = c(
      "CD4", "IL7R", "CCR7", "SELL", "LEF1", "TCF7", "MAL",
      "TBX21", "IFNG", "TNF", "IL2", "CXCR3", "CCR5", "STAT4",
      "IL12RB1", "IL12RB2",
      "GATA3", "IL4", "IL5", "IL13", "CCR4", "STAT6", "IL4R",
      "RORC", "IL17A", "IL17F", "IL22", "IL23R", "CCR6", "STAT3",
      "FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2", "TIGIT",
      "CXCR5", "BCL6", "ICOS", "PDCD1", "IL21", "SH2D1A"
    ),
    T_Cell_CD8 = c(
      "CD8A", "CD8B", "EOMES", "TBX21", "RUNX3", "PRDM1", "ID2",
      "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "NKG7",
      "FASLG", "IFNG", "TNF", "LTA", "CSF2", "IL2",
      "CX3CR1", "CXCR3", "CCR5", "CXCR6", "CCR7",
      "KLRG1", "KLRD1", "KLRC1", "KLRK1", "CD69",
      "PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4",
      "ZEB2", "ZNF683", "TOX", "TCF7"
    ),
    TCR_Signaling = c(
      "CD3D", "CD3E", "CD3G", "CD247", "TRAC", "TRBC1", "TRBC2",
      "LCK", "FYN", "ZAP70", "LAT", "SLP76", "ITK",
      "GRAP2", "GRB2", "VAV1", "NCK1",
      "PLCG1", "PRKCQ", "RASGRP1",
      "MAP3K7", "MAP2K1", "MAPK1", "MAPK3",
      "CARD11", "BCL10", "MALT1", "NFKB1", "RELA", "NFKBIA",
      "NFATC1", "NFATC2", "PPP3CA", "PPP3CB", "PPP3CC",
      "FOS", "JUN"
    ),
    NK_Cell = c(
      "NCAM1", "NCR1", "NCR2", "NCR3", "FCGR3A", "KLRD1", "KLRF1",
      "KIR2DL1", "KIR2DL3", "KIR3DL1", "KIR3DL2", "KIR2DS4",
      "KLRC1", "KLRC2", "KLRB1", "KLRK1", "KLRG1",
      "GZMB", "GZMA", "PRF1", "GNLY", "NKG7", "FASLG",
      "IFNG", "TNF", "XCL1", "CCL3", "CCL4", "CCL5",
      "EOMES", "TBX21", "ID2"
    ),
    Microglia_Activated = c(
      "TREM2", "TYROBP", "APOE", "LPL", "CST7", "ITGAX", "CLEC7A",
      "SPP1", "GPNMB", "LGALS3", "CD9", "FABP5", "CTSB", "CTSD",
      "CTSL", "CTSS",
      "C1QA", "C1QB", "C1QC", "C3", "C3AR1",
      "IL1B", "IL6", "TNF", "CCL2", "CCL3", "CCL4", "CXCL10",
      "CD68", "CD14", "FCGR1A", "FCGR2A", "MSR1", "CD36",
      "MS4A4A", "MS4A6A", "INPP5D", "CD33", "PLCG2", "ABI3", "GRN",
      "SPI1", "IRF8", "RUNX1",
      "AIF1", "ITGAM", "CSF1R"
    ),
    Microglia_Homeostatic = c(
      "P2RY12", "P2RY13", "TMEM119", "SALL1", "HEXB", "CST3", "SPARC",
      "CX3CR1", "SELPLG", "SIGLECH", "OLFML3", "GPR34", "TGFBR1",
      "BDNF", "IGF1", "TGFB1", "IL10",
      "MERTK", "GAS6", "PROS1",
      "SLC2A5", "GLUL", "GPX1", "SOD1", "SOD2",
      "FCRLS", "MEF2C"
    ),
    Abeta_Clearance = c(
      "IDE", "MME", "ECE1", "ECE2", "ACE", "THOP1", "PREP",
      "MMP2", "MMP9", "ADAMTS4",
      "LRP1", "LRP2", "LDLR", "VLDLR", "SORL1", "SCARB1",
      "ABCA1", "ABCA7", "ABCB1", "ABCG1", "ABCG2",
      "TREM2", "CD36", "SCARA1", "MSR1", "CD14", "TLR2", "TLR4",
      "BECN1", "ATG5", "ATG7", "SQSTM1", "TFEB",
      "PSMB5", "PSMB6", "PSMA1",
      "HSPA1A", "HSPA8", "CLU"
    ),
    APP_Metabolism = c(
      "PSEN1", "PSEN2", "NCSTN", "APH1A", "APH1B", "PSENEN",
      "BACE1", "BACE2",
      "ADAM10", "ADAM17",
      "APP", "APLP1", "APLP2",
      "SORL1", "SORCS1", "LRP1", "BIN1", "PICALM", "CD2AP",
      "RAB5A", "VPS35", "SNX27",
      "APBB1"
    )
  )
}

run_gene_property <- function(magma_df, markers, ct_name) {
  df <- magma_df[!is.na(magma_df$SYMBOL) & !is.na(magma_df$ZSTAT), ]
  df$CellType <- as.numeric(toupper(df$SYMBOL) %in% toupper(markers))
  if (sum(df$CellType) < 10) return(NULL)

  mod <- summary(lm(ZSTAT ~ CellType, data = df))
  ct  <- mod$coefficients

  data.frame(
    CellType = ct_name,
    N_Genes  = sum(df$CellType),
    Beta     = ct["CellType", 1],
    SE       = ct["CellType", 2],
    T_Stat   = ct["CellType", 3],
    P        = ct["CellType", 4]
  )
}

run_gene_property_analysis <- function(eoad_magma, load_magma) {
  cell_type_markers <- get_cell_type_markers()
  results <- data.frame()

  for (gwas_name in c("EOAD", "LOAD")) {
    magma_df <- if (gwas_name == "EOAD") eoad_magma else load_magma
    for (ct in names(cell_type_markers)) {
      res <- run_gene_property(magma_df, cell_type_markers[[ct]], ct)
      if (!is.null(res)) {
        res$GWAS <- gwas_name
        results <- rbind(results, res)
      }
    }
  }

  results$Bonferroni_Sig <- results$P < 0.005
  results$FDR <- p.adjust(results$P, method = "BH")

  cat("\nGene-Property Analysis:\n")
  for (g in c("EOAD", "LOAD")) {
    cat(sprintf("  %s:\n", g))
    sub <- results[results$GWAS == g, ]
    sub <- sub[order(sub$P), ]
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("    %-25s beta=%.3f P=%.2e %s\n",
                  sub$CellType[i], sub$Beta[i], sub$P[i],
                  ifelse(sub$P[i] < 0.005, "**",
                         ifelse(sub$P[i] < 0.05, "*", ""))))
    }
  }

  return(results)
}


# ==============================================================================
# 5. Conditional Regression (T Cell Independence Testing)
# ==============================================================================
# Model: Z_MAGMA = beta_0 + beta_1 * Immune_Cell + beta_2 * Microglia
# Tests whether T cell / NK cell associations persist after controlling for
# microglial gene expression.
# ==============================================================================

run_conditional_regression <- function(magma_df, target_genes, target_name,
                                       covariate_genes) {
  df <- magma_df[!is.na(magma_df$SYMBOL) & !is.na(magma_df$ZSTAT), ]
  df$SYMBOL    <- toupper(df$SYMBOL)
  df$Target    <- as.numeric(df$SYMBOL %in% toupper(target_genes))
  df$Covariate <- as.numeric(df$SYMBOL %in% toupper(covariate_genes))

  if (sum(df$Target) < 10) return(NULL)

  mod_u <- summary(lm(ZSTAT ~ Target, data = df))$coefficients
  mod_c <- summary(lm(ZSTAT ~ Target + Covariate, data = df))$coefficients

  if (nrow(mod_c) < 3) return(NULL)

  data.frame(
    CellType      = target_name,
    N_Target      = sum(df$Target),
    N_Overlap     = sum(df$Target == 1 & df$Covariate == 1),
    Beta_Uncond   = mod_u["Target", 1],
    P_Uncond      = mod_u["Target", 4],
    Beta_Cond     = mod_c["Target", 1],
    SE_Cond       = mod_c["Target", 2],
    P_Cond        = mod_c["Target", 4],
    Beta_Microglia = mod_c["Covariate", 1],
    P_Microglia   = mod_c["Covariate", 4]
  )
}

run_conditional_regression_all <- function(eoad_magma, load_magma) {
  cell_type_markers <- get_cell_type_markers()

  microglia_covariate <- unique(toupper(c(
    cell_type_markers$Microglia_Activated,
    cell_type_markers$Microglia_Homeostatic
  )))

  immune_cell_types <- c("T_Cell_Pan", "T_Cell_CD4", "T_Cell_CD8",
                         "TCR_Signaling", "NK_Cell")

  results <- data.frame()
  for (gwas_name in c("EOAD", "LOAD")) {
    magma_df <- if (gwas_name == "EOAD") eoad_magma else load_magma
    for (ct in immune_cell_types) {
      res <- run_conditional_regression(magma_df, cell_type_markers[[ct]],
                                        ct, microglia_covariate)
      if (!is.null(res)) {
        res$GWAS <- gwas_name
        results <- rbind(results, res)
      }
    }
  }

  cat("\nConditional Regression (controlling for microglia):\n")
  for (g in c("EOAD", "LOAD")) {
    cat(sprintf("  %s:\n", g))
    sub <- results[results$GWAS == g, ]
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("    %-20s Uncond P=%.4f -> Cond P=%.4f\n",
                  sub$CellType[i], sub$P_Uncond[i], sub$P_Cond[i]))
    }
  }

  return(results)
}


# ==============================================================================
# 6. LOAD Downsampling Validation (Deterministic SE Scaling)
# ==============================================================================
# Downsample LOAD GWAS to match EOAD effective sample size.
# SE_new = SE_old / sqrt(N_eff_EOAD / N_eff_LOAD).
# This is a closed-form transformation with no stochastic component,
# producing a single canonical downsampled dataset.
# ==============================================================================

downsample_gwas <- function(gwas_df, scale_factor) {
  ds <- copy(gwas_df)
  ds$beta_ds <- ds$beta
  ds$se_ds   <- ds$se / scale_factor
  ds$Z_ds    <- ds$beta_ds / ds$se_ds
  ds$pval_ds <- 2 * pnorm(-abs(ds$Z_ds))
  return(ds)
}

run_downsampling_validation <- function(load_gwas_file,
                                        eoad_n_case = 1573,
                                        eoad_n_ctrl = 199505,
                                        ds_magma_file = NULL,
                                        output_dir = "results") {

  load_gwas <- fread(load_gwas_file)

  N_eff_EOAD <- 4 / (1 / eoad_n_case + 1 / eoad_n_ctrl)
  N_case_LOAD <- load_gwas$n_cases[1]
  N_ctrl_LOAD <- load_gwas$n_controls[1]
  N_eff_LOAD  <- 4 / (1 / N_case_LOAD + 1 / N_ctrl_LOAD)

  scaling_factor <- sqrt(N_eff_EOAD / N_eff_LOAD)

  cat(sprintf("Downsampling: LOAD N_eff=%.0f -> EOAD N_eff=%.0f (scale=%.4f)\n",
              N_eff_LOAD, N_eff_EOAD, scaling_factor))

  ds <- downsample_gwas(load_gwas, scaling_factor)

  ds_dir <- file.path(output_dir, "downsampling")
  dir.create(ds_dir, showWarnings = FALSE, recursive = TRUE)
  ds_snp <- ds[, .(SNP = variant_id, P = pval_ds)]
  fwrite(ds_snp, file.path(ds_dir, "LOAD_downsampled_snps.txt"), sep = "\t")
  cat("  Downsampled SNP-level data written.\n")

  # If pre-computed MAGMA results on downsampled data are available
  downsample_summary <- data.frame()
  if (!is.null(ds_magma_file) && file.exists(ds_magma_file)) {
    ds_magma <- read_magma_results(ds_magma_file)
    n_bonf   <- sum(ds_magma$Bonferroni < 0.05)

    graham_modules <- get_graham_modules()
    ds_graham <- data.frame()
    for (mod_name in names(graham_modules)) {
      res <- test_module_enrichment(ds_magma, graham_modules[[mod_name]], mod_name)
      ds_graham <- rbind(ds_graham, res)
    }

    ds_gsea <- tryCatch({
      gene_list <- ds_magma$ZSTAT
      names(gene_list) <- ds_magma$SYMBOL
      gene_list <- sort(gene_list[!duplicated(names(gene_list))],
                        decreasing = TRUE)
      id_map <- bitr(names(gene_list), fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      gl <- gene_list[names(gene_list) %in% id_map$SYMBOL]
      names(gl) <- id_map$ENTREZID[match(names(gl), id_map$SYMBOL)]
      gl <- sort(gl[!duplicated(names(gl))], decreasing = TRUE)
      res <- gseGO(geneList = gl, OrgDb = org.Hs.eg.db, ont = "BP",
                   minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.25,
                   pAdjustMethod = "BH", verbose = FALSE)
      nrow(res@result[res@result$p.adjust < 0.05, ])
    }, error = function(e) NA)

    downsample_summary <- data.frame(
      N_Bonferroni_Genes = n_bonf,
      N_GSEA_Sig         = ds_gsea
    )

    cat(sprintf("  Bonferroni genes: %d | GSEA pathways: %s\n",
                n_bonf, ifelse(is.na(ds_gsea), "N/A", as.character(ds_gsea))))

    tbl_dir <- file.path(output_dir, "tables")
    dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
    fwrite(downsample_summary, file.path(tbl_dir, "Downsampling_Summary.csv"))
    fwrite(ds_graham, file.path(tbl_dir, "Downsampled_Graham_Enrichment.csv"))
  } else {
    cat("  MAGMA results not yet available for downsampled data.\n")
    cat("  Run MAGMA externally on the downsampled SNP file.\n")
  }

  return(list(ds_snp = ds_snp, summary = downsample_summary))
}


# ==============================================================================
# 7. Publication-Quality Visualization
# ==============================================================================

# --- 7A: Graham Module Enrichment Bar Plot ---
plot_graham_enrichment <- function(graham_results, output_dir = "results") {
  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  plot_data <- graham_results %>%
    mutate(log10_FDR    = -log10(FDR + 1e-20),
           Module_Label = gsub("_", " ", Module))

  p <- ggplot(plot_data, aes(x = reorder(Module_Label, log10_FDR),
                              y = log10_FDR, fill = GWAS)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "gray50") +
    coord_flip() +
    scale_fill_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5")) +
    labs(x = NULL, y = expression(-log[10](FDR)),
         title = "Graham Module Enrichment") +
    theme_publication +
    theme(legend.title = element_blank())

  ggsave(file.path(fig_dir, "Graham_Module_Enrichment.pdf"), p,
         width = 8, height = 5)
  return(p)
}

# --- 7B: GSEA NES Comparison Dot Plot ---
plot_gsea_nes <- function(eoad_gsea, load_gsea, output_dir = "results") {
  if (is.null(eoad_gsea) || is.null(load_gsea)) return(NULL)

  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  eoad_top <- eoad_gsea$table %>%
    filter(p.adjust < 0.05) %>% arrange(p.adjust) %>% head(20) %>%
    mutate(GWAS = "EOAD")
  load_top <- load_gsea$table %>%
    filter(p.adjust < 0.05) %>% arrange(p.adjust) %>% head(20) %>%
    mutate(GWAS = "LOAD")

  plot_data <- bind_rows(eoad_top, load_top) %>%
    mutate(Description = ifelse(nchar(Description) > 50,
                                paste0(substr(Description, 1, 47), "..."),
                                Description))

  p <- ggplot(plot_data, aes(x = NES, y = reorder(Description, NES),
                              color = GWAS, size = -log10(p.adjust))) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5")) +
    scale_size_continuous(range = c(2, 6),
                          name = expression(-log[10](FDR))) +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL,
         title = "GSEA GO-BP Pathways") +
    theme_publication +
    theme(legend.title = element_text(size = 10))

  ggsave(file.path(fig_dir, "GSEA_NES_Comparison.pdf"), p,
         width = 10, height = 8)
  return(p)
}

# --- 7C: Gene-Property Heatmap ---
plot_gene_property_heatmap <- function(gene_property_results,
                                       output_dir = "results") {
  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  gp_matrix <- gene_property_results %>%
    mutate(neg_log10_P = -log10(P + 1e-20)) %>%
    select(CellType, GWAS, neg_log10_P) %>%
    pivot_wider(names_from = GWAS, values_from = neg_log10_P) %>%
    tibble::column_to_rownames("CellType")

  gp_beta <- gene_property_results %>%
    select(CellType, GWAS, Beta) %>%
    pivot_wider(names_from = GWAS, values_from = Beta) %>%
    tibble::column_to_rownames("CellType")

  col_fun <- colorRamp2(
    c(0, -log10(0.005), max(as.matrix(gp_matrix), na.rm = TRUE)),
    c("white", "lightyellow", "red3")
  )

  pdf(file.path(fig_dir, "GeneProperty_Heatmap.pdf"), width = 6, height = 7)
  ht <- Heatmap(
    as.matrix(gp_matrix),
    name = expression(-log[10](P)),
    col  = col_fun,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    row_names_gp    = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 12),
    cell_fun = function(j, i, x, y, width, height, fill) {
      beta_val <- as.matrix(gp_beta)[i, j]
      p_val    <- 10^(-as.matrix(gp_matrix)[i, j])
      sig <- ifelse(p_val < 0.001, "***",
                    ifelse(p_val < 0.005, "**",
                           ifelse(p_val < 0.05, "*", "")))
      grid.text(sprintf("%.2f%s", beta_val, sig), x, y,
                gp = gpar(fontsize = 8))
    },
    column_title = "MAGMA Gene-Property Analysis"
  )
  draw(ht)
  dev.off()
}

# --- 7D: Conditional Regression Forest Plot ---
plot_conditional_forest <- function(conditional_results,
                                    output_dir = "results") {
  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  plot_data <- conditional_results %>%
    mutate(CellType_Label = gsub("_", " ", CellType),
           CI_Lower = Beta_Cond - 1.96 * SE_Cond,
           CI_Upper = Beta_Cond + 1.96 * SE_Cond)

  p <- ggplot(plot_data, aes(x = Beta_Cond,
                              y = reorder(CellType_Label, Beta_Cond),
                              color = GWAS)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper),
                   height = 0.2,
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5")) +
    labs(x = expression(beta ~ "(conditional on microglia)"), y = NULL,
         title = "Conditional Regression: T Cell Independence") +
    theme_publication +
    theme(legend.title = element_blank())

  ggsave(file.path(fig_dir, "Conditional_Forest.pdf"), p,
         width = 8, height = 5)
  return(p)
}

# --- 7E: Downsampling Validation Bar Plot ---
plot_downsampling <- function(eoad_gsea, load_gsea, ds_n_gsea,
                              output_dir = "results") {
  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  n_eoad <- if (!is.null(eoad_gsea)) sum(eoad_gsea$table$p.adjust < 0.05) else 0
  n_load <- if (!is.null(load_gsea)) sum(load_gsea$table$p.adjust < 0.05) else 0

  bar_data <- data.frame(
    Group      = c("EOAD\n(full)", "LOAD\n(full)", "LOAD\n(downsampled)"),
    N_Pathways = c(n_eoad, n_load, ds_n_gsea),
    SD         = c(0, 0, 0)
  )
  bar_data$Group <- factor(bar_data$Group, levels = bar_data$Group)

  p <- ggplot(bar_data, aes(x = Group, y = N_Pathways, fill = Group)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_errorbar(aes(ymin = N_Pathways - SD, ymax = N_Pathways + SD),
                  width = 0.2) +
    scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#91D1C2")) +
    labs(x = NULL, y = "Significant GO-BP Pathways (FDR < 0.05)",
         title = "LOAD Downsampling Validation") +
    theme_publication +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "Downsampling_Validation.pdf"), p,
         width = 6, height = 5)
  return(p)
}

# --- Combined multi-panel figure ---
plot_combined_figure <- function(panels, output_dir = "results") {
  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  valid_panels <- Filter(Negate(is.null), panels)
  if (length(valid_panels) == 0) return(NULL)

  labels <- LETTERS[seq_along(valid_panels)]
  combined <- plot_grid(plotlist = valid_panels, labels = labels,
                        ncol = 2, rel_heights = c(1, 1, 1))
  ggsave(file.path(fig_dir, "Combined_Pathway_Discovery.pdf"), combined,
         width = 16, height = 18)
  return(combined)
}


# ==============================================================================
# 8. Main Execution
# ==============================================================================

run_pathway_discovery <- function(eoad_magma_file,
                                   load_magma_file,
                                   load_gwas_file    = NULL,
                                   ds_magma_file     = NULL,
                                   eoad_n_case       = 1573,
                                   eoad_n_ctrl       = 199505,
                                   output_dir        = "results") {

  dir.create(file.path(output_dir, "tables"),  showWarnings = FALSE,
             recursive = TRUE)
  dir.create(file.path(output_dir, "figures"), showWarnings = FALSE,
             recursive = TRUE)

  # --- Part 1: Read MAGMA results ---
  cat("=== Part 1: MAGMA Gene-Based Association ===\n")
  eoad_magma <- read_magma_results(eoad_magma_file)
  load_magma <- read_magma_results(load_magma_file)

  cat(sprintf("EOAD: %d genes, %d Bonferroni-significant\n",
              nrow(eoad_magma), sum(eoad_magma$Bonferroni < 0.05)))
  cat(sprintf("LOAD: %d genes, %d Bonferroni-significant\n",
              nrow(load_magma), sum(load_magma$Bonferroni < 0.05)))

  fwrite(eoad_magma[eoad_magma$Bonferroni < 0.05, ],
         file.path(output_dir, "tables", "EOAD_Significant_Genes.csv"))
  fwrite(load_magma[load_magma$Bonferroni < 0.05, ],
         file.path(output_dir, "tables", "LOAD_Significant_Genes.csv"))

  # --- Part 2: Graham module enrichment ---
  cat("\n=== Part 2: Graham Module Enrichment ===\n")
  graham_results <- run_graham_enrichment(eoad_magma, load_magma)
  fwrite(graham_results,
         file.path(output_dir, "tables", "Graham_Module_Enrichment.csv"))

  # --- Part 3: GSEA GO-BP ---
  cat("\n=== Part 3: GSEA GO-BP ===\n")
  eoad_gsea <- run_gsea_gobp(eoad_magma, "EOAD")
  load_gsea <- run_gsea_gobp(load_magma, "LOAD")

  if (!is.null(eoad_gsea))
    fwrite(eoad_gsea$table,
           file.path(output_dir, "tables", "EOAD_GSEA_GOBP.csv"))
  if (!is.null(load_gsea))
    fwrite(load_gsea$table,
           file.path(output_dir, "tables", "LOAD_GSEA_GOBP.csv"))

  if (!is.null(eoad_gsea) && !is.null(load_gsea)) {
    eoad_sig <- eoad_gsea$table$ID[eoad_gsea$table$p.adjust < 0.05]
    load_sig <- load_gsea$table$ID[load_gsea$table$p.adjust < 0.05]
    cat(sprintf("  Overlap: %d pathways shared (EOAD %d, LOAD %d)\n",
                length(intersect(eoad_sig, load_sig)),
                length(eoad_sig), length(load_sig)))
  }

  # --- Part 4: Gene-property analysis ---
  cat("\n=== Part 4: Gene-Property Analysis ===\n")
  gp_results <- run_gene_property_analysis(eoad_magma, load_magma)
  fwrite(gp_results,
         file.path(output_dir, "tables", "Gene_Property_CellType.csv"))

  # --- Part 5: Conditional regression ---
  cat("\n=== Part 5: Conditional Regression ===\n")
  cond_results <- run_conditional_regression_all(eoad_magma, load_magma)
  fwrite(cond_results,
         file.path(output_dir, "tables", "Conditional_Regression.csv"))

  # --- Part 6: Downsampling validation ---
  ds_result <- NULL
  if (!is.null(load_gwas_file)) {
    cat("\n=== Part 6: Downsampling Validation ===\n")
    ds_result <- run_downsampling_validation(
      load_gwas_file, eoad_n_case, eoad_n_ctrl, ds_magma_file, output_dir
    )
  }

  # --- Part 7: Visualization ---
  cat("\n=== Part 7: Visualization ===\n")
  p_graham <- plot_graham_enrichment(graham_results, output_dir)
  p_gsea   <- plot_gsea_nes(eoad_gsea, load_gsea, output_dir)
  plot_gene_property_heatmap(gp_results, output_dir)
  p_cond   <- plot_conditional_forest(cond_results, output_dir)

  p_ds <- NULL
  if (!is.null(ds_result) && nrow(ds_result$summary) > 0) {
    p_ds <- plot_downsampling(eoad_gsea, load_gsea,
                              ds_result$summary$N_GSEA_Sig, output_dir)
  }

  # Heatmap placeholder for combined figure (generated as separate PDF)
  panels <- list(p_graham, p_gsea, ggplot() + theme_void(), p_cond, p_ds)
  plot_combined_figure(panels, output_dir)

  cat("\n=== All analyses complete ===\n")
  cat(sprintf("  Tables:  %s/tables/\n", output_dir))
  cat(sprintf("  Figures: %s/figures/\n", output_dir))

  sessionInfo()

  invisible(list(
    eoad_magma       = eoad_magma,
    load_magma       = load_magma,
    graham_results   = graham_results,
    eoad_gsea        = eoad_gsea,
    load_gsea        = load_gsea,
    gene_property    = gp_results,
    conditional      = cond_results,
    downsampling     = ds_result
  ))
}


# ==============================================================================
# Example Usage
# ==============================================================================
# results <- run_pathway_discovery(
#   eoad_magma_file = "data/EOAD_MAGMA.genes.out",
#   load_magma_file = "data/LOAD_MAGMA.genes.out",
#   load_gwas_file  = "data/LOAD_GWAS_summary_stats.tsv.gz",
#   ds_magma_file   = "results/downsampling/LOAD_downsampled.genes.out",
#   eoad_n_case     = 1573,
#   eoad_n_ctrl     = 199505,
#   output_dir      = "results"
# )
