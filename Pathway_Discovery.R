# ==============================================================================
# EOAD vs LOAD Data-Driven Pathway Discovery Analysis
# ==============================================================================
#
# Analysis Contents:
#   Part 1: Environment Setup
#   Part 2: WGCNA Co-expression Network Analysis (GSE137313)
#   Part 3: MAGMA Gene-Based Association Analysis
#   Part 4: Graham et al. 2025 Module Enrichment (Generic)
#   Part 5: GSEA and Data-Driven Pathway Discovery
#   Part 6: Cell Type MAGMA Gene-Property Analysis (Basic)
#   Part 7: Visualization (Basic)
#   Part 8: Results Summary (Initial)
#   Part 9: Extended Cell Type and Conditional Regression Analysis
#   Part 10: Graham et al. 2025 Correct Module Definitions (5 Modules)
#   Part 11: EOAD-Specific vs LOAD-Specific ORA Analysis
#   Part 12: Extended Visualization
#   Part 13: Final Results Summary
#
# ==============================================================================

rm(list = ls())
gc()

# ==============================================================================
# Part 1: Environment Setup
# ==============================================================================

packages_bioc <- c("WGCNA", "limma", "ComplexHeatmap", "biomaRt",
                   "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE")
packages_cran <- c("readxl", "dplyr", "tidyr", "ggplot2", "data.table",
                   "pheatmap", "corrplot", "igraph", "circlize", "tibble",
                   "RColorBrewer", "ggrepel", "cowplot", "gridExtra",
                   "VennDiagram", "UpSetR", "scales")

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

allowWGCNAThreads()

dir.create("results/", showWarnings = FALSE)
dir.create("results/figures/", showWarnings = FALSE)
dir.create("results/tables/", showWarnings = FALSE)

# ==============================================================================
# Part 2: WGCNA Co-expression Network Analysis
# ==============================================================================

## Section 2.1: Load GSE137313 Mouse Hippocampus Data
expr_data <- read_excel("GSE137313_geneQuantification_TPM_sampleinfo_MOD_DS_11May18.xlsx")
expr_data <- as.data.frame(expr_data)
expr_data <- expr_data[!duplicated(expr_data[,1]),]
rownames(expr_data) <- NULL
expr_data <- column_to_rownames(expr_data, var = "SampleName")
colnames(expr_data) <- expr_data[1,]
colnames(expr_data) <- gsub("RNA_", "", colnames(expr_data))

pathology_data <- read.csv("pathology updated for RNAseq Apr17.csv")
pathology_data$Sample.ID <- gsub("_133_2", "", pathology_data$Sample.ID)
expr_data <- expr_data[, colnames(expr_data) %in% pathology_data$Sample.ID]

rownames(pathology_data) <- pathology_data[,5]
allTraits <- pathology_data[,-c(1,2,4,6,7,10)]
allTraits <- allTraits[, -2]
rownames(allTraits) <- rownames(pathology_data)

sample_names <- colnames(expr_data)
traitRows <- match(sample_names, rownames(allTraits))
sample_info <- as.data.frame(allTraits[traitRows, ])

## Section 2.2: Batch Assignment
batch1_samples <- c("1145_HIP", "1462_HIP", "1463_HIP", "1472_HIP", "1816_HIP",
                    "1818_HIP", "1825_HIP", "2152_HIP", "2393A_HIP", "2417_HIP",
                    "454_HIP", "835_HIP", "1150_HIP", "1474_HIP", "1838_HIP",
                    "1839_HIP", "2098_HIP", "2099_HIP", "2337_HIP", "2392_HIP",
                    "451_HIP", "809_HIP", "811_HIP", "812_HIP", "113_HIP",
                    "1151_HIP", "1161_HIP", "154_HIP", "1819_HIP", "2153_HIP",
                    "810_HIP", "97_9_1_HIP", "1105_HIP", "1461_HIP", "1494_HIP",
                    "1808_HIP", "1817_HIP", "2330_HIP", "2416_HIP", "452_HIP", "453_HIP")

batch2_samples <- c("2187_HIP", "2297_HIP", "2314_HIP", "2350_HIP", "2351_HIP",
                    "2431_HIP", "968_HIP", "1160_HIP", "2125_HIP", "2306_HIP",
                    "2344_HIP", "2374_HIP", "2375_HIP", "2434_HIP", "1106_HIP")

sample_info$batch <- rep(NA)
sample_info[rownames(sample_info) %in% batch1_samples, ]$batch <- 1
sample_info[rownames(sample_info) %in% batch2_samples, ]$batch <- 2
sample_info[is.na(sample_info$batch), ]$batch <- 3

sample_info$AGE_cat <- as.numeric(sample_info$AGE_cat)
sample_info$MODELDIS <- as.factor(sample_info$MODELDIS)

sample_info$Age_Group <- factor(
  ifelse(sample_info$AGE_cat == 8, "Young_2m",
         ifelse(sample_info$AGE_cat == 16, "Middle_4m",
                ifelse(sample_info$AGE_cat == 32, "Old_8m", "Aged_18m"))),
  levels = c("Young_2m", "Middle_4m", "Old_8m", "Aged_18m")
)

sample_info$Disease_Model <- factor(
  ifelse(sample_info$MODELDIS == 1, "WT",
         ifelse(sample_info$MODELDIS %in% c(2,3), "AD_Model", "Other")),
  levels = c("WT", "AD_Model", "Other")
)

expr_matrix <- as.matrix(expr_data[-1, ])
expr_matrix <- apply(expr_matrix, 2, as.numeric)
rownames(expr_matrix) <- rownames(expr_data)[-1]

## Section 2.3: Filter Genes and Run WGCNA
expr_var <- apply(expr_matrix, 1, var, na.rm = TRUE)
expr_mean <- apply(expr_matrix, 1, mean, na.rm = TRUE)
keep_genes <- expr_var > quantile(expr_var, 0.25, na.rm = TRUE) &
  expr_mean > quantile(expr_mean, 0.2, na.rm = TRUE) &
  !is.na(expr_var) & !is.na(expr_mean)

expr_filtered <- expr_matrix[keep_genes, ]
datExpr <- t(expr_filtered)

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

## Section 2.4: Soft Threshold Selection
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
soft_power <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)

## Section 2.5: Build Network
net <- blockwiseModules(
  datExpr,
  power = soft_power,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.15,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 0,
  maxBlockSize = 20000,
  networkType = "signed",
  corType = "bicor",
  deepSplit = 2
)

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

## Section 2.6: Module-Trait Correlation
trait_samples <- rownames(sample_info)[rownames(sample_info) %in% rownames(datExpr)]
trait_data <- sample_info[trait_samples, ]
trait_data <- trait_data[match(rownames(datExpr), rownames(trait_data)), ]

trait_matrix <- data.frame(
  Age = as.numeric(trait_data$AGE_cat),
  Age_Young = as.numeric(trait_data$AGE_cat == 8),
  Age_Middle = as.numeric(trait_data$AGE_cat == 16),
  Age_Old = as.numeric(trait_data$AGE_cat == 32),
  Age_Aged = as.numeric(trait_data$AGE_cat == 72),
  AD_Model = as.numeric(trait_data$MODELDIS != 1),
  WT = as.numeric(trait_data$MODELDIS == 1)
)
rownames(trait_matrix) <- rownames(trait_data)

MEs_ordered <- orderMEs(MEs)
moduleTraitCor <- cor(MEs_ordered, trait_matrix, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

age_positive_modules <- rownames(moduleTraitCor)[
  moduleTraitCor[, "Age"] > 0.2 & moduleTraitPvalue[, "Age"] < 0.1
]
age_negative_modules <- rownames(moduleTraitCor)[
  moduleTraitCor[, "Age"] < -0.2 & moduleTraitPvalue[, "Age"] < 0.1
]

## Section 2.7: Save WGCNA Results
module_info <- data.frame(
  Gene = colnames(datExpr),
  Module_Color = moduleColors
)
fwrite(module_info, "results/tables/Table_WGCNA_Module_Assignment.csv")

wgcna_results <- list(
  datExpr = datExpr,
  net = net,
  moduleColors = moduleColors,
  MEs = MEs_ordered,
  moduleTraitCor = moduleTraitCor,
  moduleTraitPvalue = moduleTraitPvalue,
  soft_power = soft_power,
  age_positive_modules = age_positive_modules,
  age_negative_modules = age_negative_modules
)
saveRDS(wgcna_results, "results/WGCNA_Results.rds")


# ==============================================================================
# Part 3: MAGMA Gene-Based Association Analysis
# ==============================================================================

## Section 3.1: Mouse-to-Human Gene Conversion
convert_mouse_to_human <- function(mouse_genes) {
  tryCatch({
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    conversion <- getLDS(
      attributes = c("mgi_symbol"),
      filters = "mgi_symbol",
      values = mouse_genes,
      mart = mouse,
      attributesL = c("hgnc_symbol"),
      martL = human,
      uniqueRows = TRUE
    )
    
    colnames(conversion) <- c("Mouse_Gene", "Human_Gene")
    return(conversion)
  }, error = function(e) {
    message("biomaRt conversion failed, using manual mapping")
    return(NULL)
  })
}

## Section 3.2: Read MAGMA Results
read_magma_results <- function(magma_file) {
  if (!file.exists(magma_file)) {
    warning(paste("MAGMA file not found:", magma_file))
    return(NULL)
  }
  
  magma <- fread(magma_file)
  
  if ("P" %in% colnames(magma)) {
    magma$PVALUE <- magma$P
  }
  
  if ("GENE" %in% colnames(magma) && !("SYMBOL" %in% colnames(magma))) {
    magma$SYMBOL <- magma$GENE
  }
  
  magma$FDR <- p.adjust(magma$PVALUE, method = "BH")
  magma$Bonferroni <- p.adjust(magma$PVALUE, method = "bonferroni")
  
  return(magma)
}

## Section 3.3: Load EOAD and LOAD MAGMA Results
eoad_magma <- read_magma_results("EOAD_MAGMA.genes.out")
load_magma <- read_magma_results("LOAD_MAGMA.genes.out")

## Section 3.4: Module Enrichment Test
test_module_enrichment <- function(module_genes, magma_results, background_genes = NULL) {
  if (is.null(magma_results)) return(NULL)
  
  if (is.null(background_genes)) {
    background_genes <- magma_results$SYMBOL
  }
  
  module_genes <- intersect(module_genes, background_genes)
  if (length(module_genes) < 5) return(NULL)
  
  magma_genes <- magma_results$SYMBOL[magma_results$PVALUE < 0.05]
  
  overlap <- intersect(module_genes, magma_genes)
  
  contingency <- matrix(c(
    length(overlap),
    length(setdiff(module_genes, magma_genes)),
    length(setdiff(magma_genes, module_genes)),
    length(setdiff(background_genes, union(module_genes, magma_genes)))
  ), nrow = 2)
  
  fisher_result <- fisher.test(contingency, alternative = "greater")
  
  return(list(
    n_module = length(module_genes),
    n_overlap = length(overlap),
    overlap_genes = overlap,
    odds_ratio = fisher_result$estimate,
    p_value = fisher_result$p.value
  ))
}

## Section 3.5: Test All Modules
if (exists("wgcna_results")) {
  module_colors_unique <- unique(wgcna_results$moduleColors)
  module_colors_unique <- module_colors_unique[module_colors_unique != "grey"]
  
  module_genes_list <- lapply(module_colors_unique, function(color) {
    colnames(wgcna_results$datExpr)[wgcna_results$moduleColors == color]
  })
  names(module_genes_list) <- module_colors_unique
  
  mouse_genes_all <- unique(unlist(module_genes_list))
  gene_conversion <- convert_mouse_to_human(mouse_genes_all)
  
  if (!is.null(gene_conversion)) {
    module_genes_human <- lapply(module_genes_list, function(genes) {
      human_genes <- gene_conversion$Human_Gene[gene_conversion$Mouse_Gene %in% genes]
      unique(human_genes[human_genes != ""])
    })
    
    eoad_enrichment <- lapply(module_genes_human, function(genes) {
      test_module_enrichment(genes, eoad_magma)
    })
    
    load_enrichment <- lapply(module_genes_human, function(genes) {
      test_module_enrichment(genes, load_magma)
    })
    
    enrichment_summary <- data.frame(
      Module = names(module_genes_human),
      N_Genes = sapply(module_genes_human, length),
      EOAD_Overlap = sapply(eoad_enrichment, function(x) if(!is.null(x)) x$n_overlap else NA),
      EOAD_OR = sapply(eoad_enrichment, function(x) if(!is.null(x)) x$odds_ratio else NA),
      EOAD_P = sapply(eoad_enrichment, function(x) if(!is.null(x)) x$p_value else NA),
      LOAD_Overlap = sapply(load_enrichment, function(x) if(!is.null(x)) x$n_overlap else NA),
      LOAD_OR = sapply(load_enrichment, function(x) if(!is.null(x)) x$odds_ratio else NA),
      LOAD_P = sapply(load_enrichment, function(x) if(!is.null(x)) x$p_value else NA)
    )
    
    enrichment_summary$EOAD_FDR <- p.adjust(enrichment_summary$EOAD_P, method = "BH")
    enrichment_summary$LOAD_FDR <- p.adjust(enrichment_summary$LOAD_P, method = "BH")
    
    fwrite(enrichment_summary, "results/tables/Table_Module_MAGMA_Enrichment.csv")
  }
}


# ==============================================================================
# Part 4: Graham et al. 2025 Module Enrichment
# ==============================================================================

## Section 4.1: Define Graham Modules
define_graham_modules <- function() {
  graham_modules <- list(
    M1_Synaptic = c("SYN1", "SYN2", "SYP", "SNAP25", "VAMP2", "STX1A", "SYT1",
                    "NRXN1", "NLGN1", "SHANK3", "DLG4", "HOMER1", "GRIN1", "GRIN2A",
                    "GRIN2B", "GRIA1", "GRIA2", "CAMK2A", "CAMK2B"),
    
    M2_Myelination = c("MBP", "PLP1", "MOG", "MAG", "MOBP", "CNP", "OLIG1", "OLIG2",
                       "SOX10", "NKX2-2", "MYRF", "CLDN11", "GJC2", "ERMN", "FA2H",
                       "UGT8", "GALC", "ASPA", "ENPP6"),
    
    M3_Immune = c("CD68", "AIF1", "ITGAM", "CX3CR1", "TREM2", "TYROBP", "CSF1R",
                  "C1QA", "C1QB", "C1QC", "C3", "ITGAX", "CD14", "TLR2", "TLR4",
                  "IL1B", "TNF", "IL6", "CCL2", "CXCL10"),
    
    M4_Astrocyte = c("GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A2", "SLC1A3", "GJA1",
                     "KCNJ10", "GLUL", "SOX9", "NFIA", "NFIB", "FABP7", "VIM",
                     "CLU", "APOE", "SERPINA3"),
    
    M5_Neuronal = c("RBFOX3", "MAP2", "TUBB3", "NEFL", "NEFM", "NEFH", "ENO2",
                    "SYN1", "GAD1", "GAD2", "SLC17A7", "SLC17A6", "TH", "CHAT",
                    "SLC6A3", "SLC6A4", "HTR2A", "DRD1", "DRD2"),
    
    M6_Mitochondria = c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", "MT-ND6",
                        "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-CYB",
                        "NDUFA1", "NDUFB1", "SDHA", "SDHB", "UQCRC1", "COX4I1", "ATP5A1"),
    
    M7_Proteostasis = c("HSPA1A", "HSPA5", "HSP90AA1", "DNAJB1", "HSPB1", "BAG3",
                        "SQSTM1", "NBR1", "OPTN", "CALCOCO2", "ATG5", "ATG7", "BECN1",
                        "LC3B", "LAMP1", "LAMP2", "CTSD", "CTSB", "PSMA1", "PSMB1"),
    
    M8_Lipid = c("APOE", "CLU", "ABCA1", "ABCA7", "LDLR", "LRP1", "SREBF1", "SREBF2",
                 "HMGCR", "FDFT1", "SQLE", "CYP46A1", "CH25H", "SOAT1", "NPC1", "NPC2",
                 "PLIN2", "PLIN3", "FASN", "SCD")
  )
  
  return(graham_modules)
}

## Section 4.2: Test Graham Module Enrichment
test_graham_enrichment <- function(magma_results, graham_modules) {
  if (is.null(magma_results)) return(NULL)
  
  results <- lapply(names(graham_modules), function(module_name) {
    module_genes <- graham_modules[[module_name]]
    enrichment <- test_module_enrichment(module_genes, magma_results)
    
    if (!is.null(enrichment)) {
      data.frame(
        Module = module_name,
        N_Module = enrichment$n_module,
        N_Overlap = enrichment$n_overlap,
        OR = enrichment$odds_ratio,
        P_value = enrichment$p_value,
        Overlap_Genes = paste(enrichment$overlap_genes, collapse = ";")
      )
    } else {
      NULL
    }
  })
  
  results_df <- do.call(rbind, results[!sapply(results, is.null)])
  if (!is.null(results_df)) {
    results_df$FDR <- p.adjust(results_df$P_value, method = "BH")
  }
  
  return(results_df)
}

## Section 4.3: Run Graham Enrichment Analysis
graham_modules <- define_graham_modules()

eoad_graham <- test_graham_enrichment(eoad_magma, graham_modules)
load_graham <- test_graham_enrichment(load_magma, graham_modules)

if (!is.null(eoad_graham) && !is.null(load_graham)) {
  graham_comparison <- merge(
    eoad_graham[, c("Module", "N_Overlap", "OR", "P_value", "FDR")],
    load_graham[, c("Module", "N_Overlap", "OR", "P_value", "FDR")],
    by = "Module", suffixes = c("_EOAD", "_LOAD")
  )
  
  graham_comparison$Differential <- log2(graham_comparison$OR_EOAD / graham_comparison$OR_LOAD)
  graham_comparison$Differential[!is.finite(graham_comparison$Differential)] <- 0
  
  fwrite(graham_comparison, "results/tables/Table_Graham_Module_Comparison.csv")
}

graham_results <- list(
  modules = graham_modules,
  eoad = eoad_graham,
  load = load_graham,
  comparison = if(exists("graham_comparison")) graham_comparison else NULL
)
saveRDS(graham_results, "results/Graham_Module_Results.rds")


# ==============================================================================
# Part 5: GSEA and Data-Driven Pathway Discovery
# ==============================================================================

## Section 5.1: Convert Gene Symbols to Entrez IDs
convert_to_entrez <- function(gene_symbols) {
  tryCatch({
    entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, 
                         column = "ENTREZID", keytype = "SYMBOL",
                         multiVals = "first")
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    return(entrez_ids)
  }, error = function(e) {
    warning("Entrez ID conversion failed")
    return(NULL)
  })
}

## Section 5.2: Run GSEA Analysis
run_gsea_analysis <- function(magma_results, ont = "BP", pvalueCutoff = 0.05) {
  if (is.null(magma_results)) return(NULL)
  
  gene_list <- -log10(magma_results$PVALUE)
  names(gene_list) <- magma_results$SYMBOL
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!is.na(gene_list) & is.finite(gene_list)]
  
  entrez_map <- convert_to_entrez(names(gene_list))
  if (is.null(entrez_map) || length(entrez_map) < 100) return(NULL)
  
  gene_list_entrez <- gene_list[names(gene_list) %in% names(entrez_map)]
  names(gene_list_entrez) <- entrez_map[names(gene_list_entrez)]
  
  gsea_result <- tryCatch({
    gseGO(geneList = gene_list_entrez,
          OrgDb = org.Hs.eg.db,
          ont = ont,
          minGSSize = 10,
          maxGSSize = 500,
          pvalueCutoff = pvalueCutoff,
          verbose = FALSE)
  }, error = function(e) NULL)
  
  return(gsea_result)
}

## Section 5.3: Run ORA Enrichment
run_ora_enrichment <- function(gene_list, background = NULL, ont = "BP") {
  if (length(gene_list) < 10) return(NULL)
  
  entrez_genes <- convert_to_entrez(gene_list)
  if (is.null(entrez_genes) || length(entrez_genes) < 5) return(NULL)
  
  entrez_background <- NULL
  if (!is.null(background)) {
    entrez_background <- convert_to_entrez(background)
  }
  
  ora_result <- tryCatch({
    enrichGO(gene = entrez_genes,
             universe = entrez_background,
             OrgDb = org.Hs.eg.db,
             ont = ont,
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2,
             readable = TRUE)
  }, error = function(e) NULL)
  
  return(ora_result)
}

## Section 5.4: Run KEGG Enrichment
run_kegg_enrichment <- function(gene_list, background = NULL) {
  if (length(gene_list) < 10) return(NULL)
  
  entrez_genes <- convert_to_entrez(gene_list)
  if (is.null(entrez_genes) || length(entrez_genes) < 5) return(NULL)
  
  entrez_background <- NULL
  if (!is.null(background)) {
    entrez_background <- convert_to_entrez(background)
  }
  
  kegg_result <- tryCatch({
    enrichKEGG(gene = entrez_genes,
               universe = entrez_background,
               organism = "hsa",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)
  }, error = function(e) NULL)
  
  return(kegg_result)
}

## Section 5.5: Run GSEA for EOAD and LOAD
eoad_gsea_bp <- run_gsea_analysis(eoad_magma, ont = "BP")
load_gsea_bp <- run_gsea_analysis(load_magma, ont = "BP")

eoad_gsea_mf <- run_gsea_analysis(eoad_magma, ont = "MF")
load_gsea_mf <- run_gsea_analysis(load_magma, ont = "MF")

## Section 5.6: Extract Significant Genes for ORA
if (!is.null(eoad_magma)) {
  eoad_sig_genes <- eoad_magma$SYMBOL[eoad_magma$FDR < 0.05]
  eoad_nominal_genes <- eoad_magma$SYMBOL[eoad_magma$PVALUE < 0.05]
  
  eoad_ora_bp <- run_ora_enrichment(eoad_nominal_genes, eoad_magma$SYMBOL, "BP")
  eoad_kegg <- run_kegg_enrichment(eoad_nominal_genes, eoad_magma$SYMBOL)
}

if (!is.null(load_magma)) {
  load_sig_genes <- load_magma$SYMBOL[load_magma$FDR < 0.05]
  load_nominal_genes <- load_magma$SYMBOL[load_magma$PVALUE < 0.05]
  
  load_ora_bp <- run_ora_enrichment(load_nominal_genes, load_magma$SYMBOL, "BP")
  load_kegg <- run_kegg_enrichment(load_nominal_genes, load_magma$SYMBOL)
}

## Section 5.7: Compare EOAD vs LOAD Pathways
compare_pathways <- function(eoad_result, load_result) {
  if (is.null(eoad_result) || is.null(load_result)) return(NULL)
  
  eoad_df <- as.data.frame(eoad_result)
  load_df <- as.data.frame(load_result)
  
  if (nrow(eoad_df) == 0 || nrow(load_df) == 0) return(NULL)
  
  all_pathways <- union(eoad_df$Description, load_df$Description)
  
  comparison <- data.frame(
    Pathway = all_pathways,
    EOAD_NES = sapply(all_pathways, function(p) {
      idx <- which(eoad_df$Description == p)
      if (length(idx) > 0) eoad_df$NES[idx[1]] else NA
    }),
    EOAD_padj = sapply(all_pathways, function(p) {
      idx <- which(eoad_df$Description == p)
      if (length(idx) > 0) eoad_df$p.adjust[idx[1]] else NA
    }),
    LOAD_NES = sapply(all_pathways, function(p) {
      idx <- which(load_df$Description == p)
      if (length(idx) > 0) load_df$NES[idx[1]] else NA
    }),
    LOAD_padj = sapply(all_pathways, function(p) {
      idx <- which(load_df$Description == p)
      if (length(idx) > 0) load_df$p.adjust[idx[1]] else NA
    })
  )
  
  comparison$EOAD_Sig <- comparison$EOAD_padj < 0.05
  comparison$LOAD_Sig <- comparison$LOAD_padj < 0.05
  comparison$Differential <- comparison$EOAD_NES - comparison$LOAD_NES
  
  comparison$Category <- ifelse(
    comparison$EOAD_Sig & !comparison$LOAD_Sig, "EOAD-specific",
    ifelse(!comparison$EOAD_Sig & comparison$LOAD_Sig, "LOAD-specific",
           ifelse(comparison$EOAD_Sig & comparison$LOAD_Sig, "Shared", "Neither"))
  )
  
  return(comparison)
}

pathway_comparison <- compare_pathways(eoad_gsea_bp, load_gsea_bp)
if (!is.null(pathway_comparison)) {
  fwrite(pathway_comparison, "results/tables/Table_GSEA_Pathway_Comparison.csv")
}

## Section 5.8: Save GSEA Results
gsea_results <- list(
  EOAD_BP = eoad_gsea_bp,
  LOAD_BP = load_gsea_bp,
  EOAD_MF = eoad_gsea_mf,
  LOAD_MF = load_gsea_mf,
  EOAD_ORA = if(exists("eoad_ora_bp")) eoad_ora_bp else NULL,
  LOAD_ORA = if(exists("load_ora_bp")) load_ora_bp else NULL,
  EOAD_KEGG = if(exists("eoad_kegg")) eoad_kegg else NULL,
  LOAD_KEGG = if(exists("load_kegg")) load_kegg else NULL,
  comparison = pathway_comparison
)
saveRDS(gsea_results, "results/GSEA_Results.rds")


# ==============================================================================
# Part 6: Cell Type MAGMA Gene-Property Analysis
# ==============================================================================

## Section 6.1: Define Cell Type Marker Genes
define_cell_type_markers <- function() {
  markers <- list(
    Excitatory_Neurons = c("SLC17A7", "CAMK2A", "GRIN1", "GRIN2A", "GRIN2B",
                           "NRGN", "SATB2", "TBR1", "NEUROD2", "NEUROD6"),
    
    Inhibitory_Neurons = c("GAD1", "GAD2", "SLC32A1", "PVALB", "SST", "VIP",
                           "CALB1", "CALB2", "NPY", "CCK"),
    
    Astrocytes = c("GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A2", "SLC1A3",
                   "GJA1", "KCNJ10", "GLUL", "SOX9"),
    
    Oligodendrocytes = c("MBP", "PLP1", "MOG", "MAG", "MOBP", "CNP", "OLIG1",
                         "OLIG2", "SOX10", "MYRF"),
    
    OPCs = c("PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10", "NKX2-2",
             "GPR17", "PCDH15", "LHFPL3", "VCAN"),
    
    Microglia = c("AIF1", "CD68", "CX3CR1", "ITGAM", "TREM2", "TYROBP",
                  "CSF1R", "P2RY12", "TMEM119", "HEXB"),
    
    Endothelial = c("CLDN5", "PECAM1", "VWF", "CDH5", "FLT1", "KDR",
                    "TIE1", "TEK", "ESAM", "ICAM2"),
    
    Pericytes = c("PDGFRB", "RGS5", "KCNJ8", "ABCC9", "DES", "ACTA2",
                  "TAGLN", "MYH11", "NOTCH3", "CSPG4")
  )
  
  return(markers)
}

## Section 6.2: Gene Property Regression
run_gene_property_regression <- function(magma_results, cell_markers) {
  if (is.null(magma_results)) return(NULL)
  
  results <- lapply(names(cell_markers), function(cell_type) {
    marker_genes <- cell_markers[[cell_type]]
    
    magma_results$is_marker <- as.integer(magma_results$SYMBOL %in% marker_genes)
    
    if (sum(magma_results$is_marker) < 5) {
      return(data.frame(
        Cell_Type = cell_type,
        N_Markers = sum(magma_results$is_marker),
        Beta = NA, SE = NA, P_value = NA
      ))
    }
    
    magma_results$ZSTAT <- qnorm(1 - magma_results$PVALUE/2)
    magma_results$ZSTAT[magma_results$PVALUE > 0.5] <- -magma_results$ZSTAT[magma_results$PVALUE > 0.5]
    
    model <- tryCatch({
      lm(ZSTAT ~ is_marker, data = magma_results)
    }, error = function(e) NULL)
    
    if (!is.null(model)) {
      coef_summary <- summary(model)$coefficients
      if ("is_marker" %in% rownames(coef_summary)) {
        return(data.frame(
          Cell_Type = cell_type,
          N_Markers = sum(magma_results$is_marker),
          Beta = coef_summary["is_marker", "Estimate"],
          SE = coef_summary["is_marker", "Std. Error"],
          P_value = coef_summary["is_marker", "Pr(>|t|)"]
        ))
      }
    }
    
    return(data.frame(
      Cell_Type = cell_type,
      N_Markers = sum(magma_results$is_marker),
      Beta = NA, SE = NA, P_value = NA
    ))
  })
  
  results_df <- do.call(rbind, results)
  results_df$FDR <- p.adjust(results_df$P_value, method = "BH")
  
  return(results_df)
}

## Section 6.3: Run Cell Type Analysis
cell_markers <- define_cell_type_markers()

eoad_celltype <- run_gene_property_regression(eoad_magma, cell_markers)
load_celltype <- run_gene_property_regression(load_magma, cell_markers)

## Section 6.4: Compare Cell Type Enrichment
if (!is.null(eoad_celltype) && !is.null(load_celltype)) {
  celltype_comparison <- merge(
    eoad_celltype[, c("Cell_Type", "Beta", "SE", "P_value", "FDR")],
    load_celltype[, c("Cell_Type", "Beta", "SE", "P_value", "FDR")],
    by = "Cell_Type", suffixes = c("_EOAD", "_LOAD")
  )
  
  celltype_comparison$Beta_Diff <- celltype_comparison$Beta_EOAD - celltype_comparison$Beta_LOAD
  
  z_diff <- (celltype_comparison$Beta_EOAD - celltype_comparison$Beta_LOAD) /
    sqrt(celltype_comparison$SE_EOAD^2 + celltype_comparison$SE_LOAD^2)
  celltype_comparison$P_Diff <- 2 * pnorm(-abs(z_diff))
  
  fwrite(celltype_comparison, "results/tables/Table_CellType_MAGMA_Comparison.csv")
}

## Section 6.5: Save Cell Type Results
celltype_results <- list(
  markers = cell_markers,
  eoad = eoad_celltype,
  load = load_celltype,
  comparison = if(exists("celltype_comparison")) celltype_comparison else NULL
)
saveRDS(celltype_results, "results/CellType_Results.rds")


# ==============================================================================
# Part 7: Visualization
# ==============================================================================

theme_nc <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

## Section 7.1: WGCNA Module-Trait Heatmap
plot_wgcna_heatmap <- function(wgcna_results) {
  if (is.null(wgcna_results)) return(NULL)
  
  cor_matrix <- wgcna_results$moduleTraitCor
  pval_matrix <- wgcna_results$moduleTraitPvalue
  
  sig_text <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  sig_text[pval_matrix < 0.001] <- "***"
  sig_text[pval_matrix >= 0.001 & pval_matrix < 0.01] <- "**"
  sig_text[pval_matrix >= 0.01 & pval_matrix < 0.05] <- "*"
  
  display_matrix <- matrix(
    paste0(round(cor_matrix, 2), sig_text),
    nrow = nrow(cor_matrix), ncol = ncol(cor_matrix)
  )
  
  col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
  
  ht <- Heatmap(
    cor_matrix,
    name = "Correlation",
    col = col_fun,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(display_matrix[i, j], x, y, gp = gpar(fontsize = 8))
    },
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_title = "Module-Trait Correlation",
    heatmap_legend_param = list(title = "r")
  )
  
  return(ht)
}

if (exists("wgcna_results")) {
  pdf("results/figures/Figure_WGCNA_ModuleTrait_Heatmap.pdf", width = 10, height = 12)
  ht <- plot_wgcna_heatmap(wgcna_results)
  if (!is.null(ht)) draw(ht)
  dev.off()
}

## Section 7.2: Graham Module Enrichment Comparison
plot_graham_enrichment <- function(graham_results) {
  if (is.null(graham_results$comparison)) return(NULL)
  
  plot_data <- graham_results$comparison
  plot_data$Module <- gsub("_", "\n", plot_data$Module)
  
  plot_data_long <- data.frame(
    Module = rep(plot_data$Module, 2),
    GWAS = rep(c("EOAD", "LOAD"), each = nrow(plot_data)),
    OR = c(plot_data$OR_EOAD, plot_data$OR_LOAD),
    P_value = c(plot_data$P_value_EOAD, plot_data$P_value_LOAD)
  )
  
  plot_data_long$Sig <- ifelse(plot_data_long$P_value < 0.05, "Significant", "NS")
  plot_data_long$OR[is.na(plot_data_long$OR)] <- 1
  plot_data_long$log2OR <- log2(plot_data_long$OR)
  
  p <- ggplot(plot_data_long, aes(x = Module, y = log2OR, fill = GWAS)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(shape = Sig), position = position_dodge(width = 0.8), 
               size = 2, show.legend = FALSE) +
    scale_fill_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5")) +
    scale_shape_manual(values = c("Significant" = 8, "NS" = NA)) +
    labs(title = "Graham Module Enrichment: EOAD vs LOAD",
         subtitle = "MAGMA Gene-Based Association",
         x = "Module", y = "log2(Odds Ratio)", fill = "GWAS") +
    theme_nc +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  
  return(p)
}

if (exists("graham_results")) {
  p_graham <- plot_graham_enrichment(graham_results)
  if (!is.null(p_graham)) {
    ggsave("results/figures/Figure_Graham_Module_Enrichment.pdf", p_graham, 
           width = 12, height = 6, dpi = 300)
  }
}

## Section 7.3: GSEA Comparison Plot
plot_gsea_comparison <- function(gsea_results, top_n = 15) {
  if (is.null(gsea_results$comparison)) return(NULL)
  
  plot_data <- gsea_results$comparison
  plot_data <- plot_data[!is.na(plot_data$EOAD_NES) | !is.na(plot_data$LOAD_NES), ]
  
  plot_data$max_abs_NES <- pmax(abs(plot_data$EOAD_NES), abs(plot_data$LOAD_NES), na.rm = TRUE)
  plot_data <- plot_data[order(-plot_data$max_abs_NES), ]
  plot_data <- head(plot_data, top_n)
  
  plot_data$Pathway <- factor(plot_data$Pathway, levels = rev(plot_data$Pathway))
  
  plot_data_long <- data.frame(
    Pathway = rep(plot_data$Pathway, 2),
    GWAS = rep(c("EOAD", "LOAD"), each = nrow(plot_data)),
    NES = c(plot_data$EOAD_NES, plot_data$LOAD_NES),
    padj = c(plot_data$EOAD_padj, plot_data$LOAD_padj)
  )
  
  plot_data_long$Sig <- ifelse(plot_data_long$padj < 0.05, "FDR<0.05", "NS")
  
  p <- ggplot(plot_data_long, aes(x = NES, y = Pathway, fill = GWAS)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5")) +
    labs(title = "GSEA Pathway Comparison: EOAD vs LOAD",
         subtitle = "Top pathways by normalized enrichment score",
         x = "Normalized Enrichment Score (NES)", y = "", fill = "GWAS") +
    theme_nc +
    theme(axis.text.y = element_text(size = 9))
  
  return(p)
}

if (exists("gsea_results")) {
  p_gsea <- plot_gsea_comparison(gsea_results)
  if (!is.null(p_gsea)) {
    ggsave("results/figures/Figure_GSEA_Pathway_Comparison.pdf", p_gsea, 
           width = 12, height = 8, dpi = 300)
  }
}

## Section 7.4: Cell Type Enrichment Plot
plot_celltype_enrichment <- function(celltype_results) {
  if (is.null(celltype_results$comparison)) return(NULL)
  
  plot_data <- celltype_results$comparison
  
  plot_data_long <- data.frame(
    Cell_Type = rep(plot_data$Cell_Type, 2),
    GWAS = rep(c("EOAD", "LOAD"), each = nrow(plot_data)),
    Beta = c(plot_data$Beta_EOAD, plot_data$Beta_LOAD),
    SE = c(plot_data$SE_EOAD, plot_data$SE_LOAD),
    FDR = c(plot_data$FDR_EOAD, plot_data$FDR_LOAD)
  )
  
  plot_data_long$Sig <- ifelse(plot_data_long$FDR < 0.05, "FDR<0.05", "NS")
  plot_data_long$Lower <- plot_data_long$Beta - 1.96 * plot_data_long$SE
  plot_data_long$Upper <- plot_data_long$Beta + 1.96 * plot_data_long$SE
  
  p <- ggplot(plot_data_long, aes(x = Cell_Type, y = Beta, fill = GWAS)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                  position = position_dodge(width = 0.8), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5")) +
    labs(title = "Cell Type Enrichment: EOAD vs LOAD",
         subtitle = "MAGMA Gene-Property Analysis",
         x = "Cell Type", y = "Enrichment (Beta)", fill = "GWAS") +
    theme_nc +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

if (exists("celltype_results")) {
  p_celltype <- plot_celltype_enrichment(celltype_results)
  if (!is.null(p_celltype)) {
    ggsave("results/figures/Figure_CellType_Enrichment.pdf", p_celltype, 
           width = 10, height = 6, dpi = 300)
  }
}

## Section 7.5: WGCNA Module Size Distribution
if (exists("wgcna_results")) {
  module_sizes <- table(wgcna_results$moduleColors)
  module_sizes <- module_sizes[names(module_sizes) != "grey"]
  
  size_df <- data.frame(
    Module = names(module_sizes),
    Size = as.numeric(module_sizes)
  )
  size_df <- size_df[order(-size_df$Size), ]
  size_df$Module <- factor(size_df$Module, levels = size_df$Module)
  
  p_size <- ggplot(size_df, aes(x = Module, y = Size, fill = Module)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = setNames(size_df$Module, size_df$Module)) +
    labs(title = "WGCNA Module Size Distribution",
         x = "Module", y = "Number of Genes") +
    theme_nc +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave("results/figures/Figure_WGCNA_Module_Sizes.pdf", p_size, 
         width = 10, height = 6, dpi = 300)
}

## Section 7.6: Differential Pathway Volcano Plot
if (exists("pathway_comparison") && !is.null(pathway_comparison)) {
  volcano_data <- pathway_comparison[!is.na(pathway_comparison$Differential), ]
  volcano_data$min_padj <- pmin(volcano_data$EOAD_padj, volcano_data$LOAD_padj, na.rm = TRUE)
  volcano_data$neg_log_p <- -log10(volcano_data$min_padj)
  
  volcano_data$Category <- ifelse(
    abs(volcano_data$Differential) > 1 & volcano_data$min_padj < 0.05, "Differential",
    ifelse(volcano_data$min_padj < 0.05, "Significant", "NS")
  )
  
  p_volcano <- ggplot(volcano_data, aes(x = Differential, y = neg_log_p, color = Category)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Differential" = "#E64B35", "Significant" = "#4DBBD5", "NS" = "gray70")) +
    labs(title = "Differential Pathway Analysis: EOAD vs LOAD",
         x = "NES Difference (EOAD - LOAD)", y = "-log10(min p-adjusted)") +
    theme_nc
  
  ggsave("results/figures/Figure_Pathway_Volcano.pdf", p_volcano, 
         width = 10, height = 8, dpi = 300)
}


# ==============================================================================
# Part 8: Results Summary
# ==============================================================================

## Section 8.1: Generate Summary Tables
summary_stats <- data.frame(
  Analysis = c("WGCNA Modules", "MAGMA Genes (EOAD)", "MAGMA Genes (LOAD)",
               "Graham Modules", "Cell Types", "GSEA Pathways"),
  N_Total = c(
    if(exists("wgcna_results")) length(unique(wgcna_results$moduleColors)) - 1 else NA,
    if(exists("eoad_magma")) nrow(eoad_magma) else NA,
    if(exists("load_magma")) nrow(load_magma) else NA,
    length(graham_modules),
    length(cell_markers),
    if(exists("pathway_comparison")) nrow(pathway_comparison) else NA
  ),
  N_Significant = c(
    if(exists("enrichment_summary")) sum(enrichment_summary$EOAD_FDR < 0.05 | enrichment_summary$LOAD_FDR < 0.05, na.rm = TRUE) else NA,
    if(exists("eoad_magma")) sum(eoad_magma$FDR < 0.05, na.rm = TRUE) else NA,
    if(exists("load_magma")) sum(load_magma$FDR < 0.05, na.rm = TRUE) else NA,
    if(exists("graham_comparison")) sum(graham_comparison$FDR_EOAD < 0.05 | graham_comparison$FDR_LOAD < 0.05, na.rm = TRUE) else NA,
    if(exists("celltype_comparison")) sum(celltype_comparison$FDR_EOAD < 0.05 | celltype_comparison$FDR_LOAD < 0.05, na.rm = TRUE) else NA,
    if(exists("pathway_comparison")) sum(pathway_comparison$EOAD_Sig | pathway_comparison$LOAD_Sig, na.rm = TRUE) else NA
  )
)

fwrite(summary_stats, "results/tables/Table_Analysis_Summary.csv")

## Section 8.2: Key Findings Summary
key_findings <- list()

if (exists("graham_comparison") && !is.null(graham_comparison)) {
  eoad_enriched <- graham_comparison$Module[graham_comparison$FDR_EOAD < 0.05]
  load_enriched <- graham_comparison$Module[graham_comparison$FDR_LOAD < 0.05]
  key_findings$graham <- list(
    eoad_specific = setdiff(eoad_enriched, load_enriched),
    load_specific = setdiff(load_enriched, eoad_enriched),
    shared = intersect(eoad_enriched, load_enriched)
  )
}

if (exists("celltype_comparison") && !is.null(celltype_comparison)) {
  eoad_celltypes <- celltype_comparison$Cell_Type[celltype_comparison$FDR_EOAD < 0.05]
  load_celltypes <- celltype_comparison$Cell_Type[celltype_comparison$FDR_LOAD < 0.05]
  key_findings$celltype <- list(
    eoad_specific = setdiff(eoad_celltypes, load_celltypes),
    load_specific = setdiff(load_celltypes, eoad_celltypes),
    shared = intersect(eoad_celltypes, load_celltypes)
  )
}

if (exists("pathway_comparison") && !is.null(pathway_comparison)) {
  key_findings$pathways <- list(
    eoad_specific = sum(pathway_comparison$Category == "EOAD-specific", na.rm = TRUE),
    load_specific = sum(pathway_comparison$Category == "LOAD-specific", na.rm = TRUE),
    shared = sum(pathway_comparison$Category == "Shared", na.rm = TRUE)
  )
}

saveRDS(key_findings, "results/Key_Findings.rds")

## Section 8.3: Save Complete Results
all_results <- list(
  wgcna = if(exists("wgcna_results")) wgcna_results else NULL,
  magma = list(eoad = eoad_magma, load = load_magma),
  graham = if(exists("graham_results")) graham_results else NULL,
  gsea = if(exists("gsea_results")) gsea_results else NULL,
  celltype = if(exists("celltype_results")) celltype_results else NULL,
  key_findings = key_findings,
  analysis_date = Sys.time(),
  r_version = R.version.string
)

saveRDS(all_results, "results/Complete_Pathway_Analysis_Results.rds")

## Section 8.4: Print Summary
cat("\n")
cat("================================================================================\n")
cat("  EOAD vs LOAD Pathway Discovery Analysis Complete\n")
cat("================================================================================\n")
cat("\nOutput files:\n")
cat("  Tables:\n")
cat("    - results/tables/Table_WGCNA_Module_Assignment.csv\n")
cat("    - results/tables/Table_Module_MAGMA_Enrichment.csv\n")
cat("    - results/tables/Table_Graham_Module_Comparison.csv\n")
cat("    - results/tables/Table_GSEA_Pathway_Comparison.csv\n")
cat("    - results/tables/Table_CellType_MAGMA_Comparison.csv\n")
cat("    - results/tables/Table_Analysis_Summary.csv\n")
cat("  Figures:\n")
cat("    - results/figures/Figure_WGCNA_ModuleTrait_Heatmap.pdf\n")
cat("    - results/figures/Figure_WGCNA_Module_Sizes.pdf\n")
cat("    - results/figures/Figure_Graham_Module_Enrichment.pdf\n")
cat("    - results/figures/Figure_GSEA_Pathway_Comparison.pdf\n")
cat("    - results/figures/Figure_CellType_Enrichment.pdf\n")
cat("    - results/figures/Figure_Pathway_Volcano.pdf\n")
cat("  RDS Objects:\n")
cat("    - results/WGCNA_Results.rds\n")
cat("    - results/Graham_Module_Results.rds\n")
cat("    - results/GSEA_Results.rds\n")
cat("    - results/CellType_Results.rds\n")
cat("    - results/Key_Findings.rds\n")
cat("    - results/Complete_Pathway_Analysis_Results.rds\n")
cat("\n")

sessionInfo()





# ==============================================================================
# Part 9: Extended Cell Type and Conditional Regression Analysis
# ==============================================================================

## Section 9.1: Define Extended Cell Type Markers (Methods-Matched)
define_extended_markers <- function() {
  list(
    T_Cell_Extended = c(
      "CD3D", "CD3E", "CD3G", "CD247", "CD2", "CD5", "CD7", "CD27", "CD28",
      "TRAC", "TRBC1", "TRBC2", "TCF7", "LEF1", "GATA3", "TBX21", "EOMES",
      "LCK", "ZAP70", "LAT", "ITK", "FYN", "PLCG1", "VAV1", "PRKCQ", "CARD11",
      "NFATC1", "NFATC2", "NFKB1", "IL7R", "IL2RA", "IL2RB", "IFNG", "TNF",
      "CCR7", "CCR4", "CCR5", "CXCR3", "CXCR4", "CTLA4", "PDCD1", "LAG3",
      "GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7", "FASLG", "SELL"
    ),
    T_Cell_CD4 = c(
      "CD4", "IL7R", "CCR7", "SELL", "LEF1", "TCF7", "TBX21", "IFNG", "TNF",
      "CXCR3", "CCR5", "GATA3", "IL4", "IL5", "IL13", "CCR4", "RORC", "IL17A",
      "IL17F", "IL22", "IL23R", "CCR6", "FOXP3", "IL2RA", "CTLA4", "CXCR5",
      "BCL6", "ICOS", "PDCD1", "IL21", "CD28", "CD40LG", "LCK", "ZAP70"
    ),
    T_Cell_CD8 = c(
      "CD8A", "CD8B", "EOMES", "TBX21", "RUNX3", "GZMA", "GZMB", "GZMH",
      "GZMK", "PRF1", "GNLY", "NKG7", "FASLG", "IFNG", "TNF", "CX3CR1",
      "CXCR3", "CCR5", "KLRG1", "KLRD1", "KLRC1", "PDCD1", "LAG3", "HAVCR2",
      "TIGIT", "ZEB2", "TOX", "TCF7", "LCK", "ZAP70", "FCGR3A", "FGFBP2"
    ),
    TCR_Signaling = c(
      "CD3D", "CD3E", "CD3G", "CD247", "TRAC", "TRBC1", "LCK", "FYN", "ZAP70",
      "LAT", "ITK", "GRAP2", "GRB2", "SOS1", "VAV1", "PLCG1", "PRKCQ",
      "RASGRP1", "MAP3K7", "MAP2K1", "MAPK1", "MAPK3", "CARD11", "BCL10",
      "MALT1", "NFKB1", "NFATC1", "NFATC2", "PPP3CA", "FOS", "JUN"
    ),
    NK_Cell = c(
      "NCAM1", "NCR1", "NCR2", "NCR3", "FCGR3A", "KLRD1", "KLRF1", "KIR2DL1",
      "KIR2DL3", "KIR3DL1", "KLRC1", "KLRC2", "KLRB1", "KLRK1", "GZMB",
      "GZMA", "PRF1", "GNLY", "NKG7", "FASLG", "IFNG", "TNF", "CCL3", "CCL4",
      "EOMES", "TBX21", "ID2", "TYROBP", "FCER1G", "SYK"
    ),
    Microglia_Activated = c(
      "TREM2", "TYROBP", "APOE", "LPL", "CST7", "ITGAX", "CLEC7A", "SPP1",
      "GPNMB", "LGALS3", "CD9", "FABP5", "CTSB", "CTSD", "CTSL", "CTSS",
      "C1QA", "C1QB", "C1QC", "C3", "IL1B", "IL6", "TNF", "CCL2", "CCL3",
      "CD68", "CD14", "FCGR1A", "MSR1", "CD36", "MS4A4A", "MS4A6A", "INPP5D",
      "CD33", "PLCG2", "ABI3", "GRN", "SPI1", "IRF8", "AIF1", "CSF1R"
    ),
    Microglia_Homeostatic = c(
      "P2RY12", "P2RY13", "TMEM119", "SALL1", "HEXB", "CST3", "SPARC",
      "CX3CR1", "SELPLG", "SIGLECH", "OLFML3", "GPR34", "TGFBR1", "BDNF",
      "IGF1", "TGFB1", "C1QA", "C1QB", "MERTK", "GAS6", "MEF2C"
    ),
    Abeta_Clearance = c(
      "IDE", "MME", "ECE1", "ACE", "THOP1", "MMP2", "MMP9", "LRP1", "LRP2",
      "LDLR", "SORL1", "SCARB1", "ABCA1", "ABCA7", "ABCB1", "ABCG1", "TREM2",
      "CD36", "MSR1", "CD14", "TLR2", "TLR4", "BECN1", "ATG5", "ATG7",
      "SQSTM1", "TFEB", "PSMB5", "HSPA1A", "HSPA8", "HSP90AA1", "CLU"
    ),
    APP_Metabolism = c(
      "PSEN1", "PSEN2", "NCSTN", "APH1A", "APH1B", "PSENEN", "BACE1", "BACE2",
      "ADAM10", "ADAM17", "APP", "APLP1", "APLP2", "SORL1", "SORCS1",
      "SORCS2", "LRP1", "BIN1", "PICALM", "CD2AP", "RAB5A", "RAB7A",
      "VPS35", "SNX27", "APBB1", "APBB2"
    )
  )
}

## Section 9.2: Conditional Regression Function
run_conditional_regression <- function(magma_df, target_genes, target_name,
                                       covariate_genes_list = NULL,
                                       covariate_names = NULL) {
  if (is.null(magma_df) || !"PVALUE" %in% colnames(magma_df)) return(NULL)
  
  df <- magma_df[!is.na(magma_df$SYMBOL) & !is.na(magma_df$PVALUE), ]
  df$SYMBOL <- toupper(df$SYMBOL)
  target_genes <- toupper(target_genes)
  
  df$Target <- as.numeric(df$SYMBOL %in% target_genes)
  n_target <- sum(df$Target)
  if (n_target < 10) return(NULL)
  
  df$ZSTAT <- qnorm(1 - df$PVALUE/2)
  df$ZSTAT[df$PVALUE > 0.5] <- -df$ZSTAT[df$PVALUE > 0.5]
  
  if (!is.null(covariate_genes_list) && length(covariate_genes_list) > 0) {
    for (i in seq_along(covariate_genes_list)) {
      cov_genes <- toupper(covariate_genes_list[[i]])
      cov_name <- ifelse(!is.null(covariate_names), covariate_names[i], paste0("Cov", i))
      df[[cov_name]] <- as.numeric(df$SYMBOL %in% cov_genes)
    }
    formula_str <- paste("ZSTAT ~ Target +", paste(covariate_names, collapse = " + "))
    model <- tryCatch(lm(as.formula(formula_str), data = df), error = function(e) NULL)
  } else {
    model <- tryCatch(lm(ZSTAT ~ Target, data = df), error = function(e) NULL)
  }
  
  if (is.null(model)) return(NULL)
  
  coef_table <- summary(model)$coefficients
  if (!"Target" %in% rownames(coef_table)) return(NULL)
  
  list(
    CellType = target_name,
    N_Target = n_target,
    N_Total = nrow(df),
    Beta = coef_table["Target", 1],
    SE = coef_table["Target", 2],
    T_Stat = coef_table["Target", 3],
    P_Value = coef_table["Target", 4],
    Mean_Z_Target = mean(df$ZSTAT[df$Target == 1], na.rm = TRUE),
    Mean_Z_Other = mean(df$ZSTAT[df$Target == 0], na.rm = TRUE),
    Conditional = !is.null(covariate_genes_list)
  )
}

## Section 9.3: Create Exclusive Gene Sets
create_exclusive_sets <- function(extended_markers) {
  microglia_all <- unique(c(extended_markers$Microglia_Activated,
                            extended_markers$Microglia_Homeostatic))
  abeta_all <- extended_markers$Abeta_Clearance
  confounding <- unique(c(microglia_all, abeta_all))
  
  t_cell_types <- c("T_Cell_Extended", "T_Cell_CD4", "T_Cell_CD8", 
                    "TCR_Signaling", "NK_Cell")
  
  exclusive <- list()
  for (ct in names(extended_markers)) {
    original <- toupper(extended_markers[[ct]])
    if (ct %in% t_cell_types) {
      exclusive[[ct]] <- setdiff(original, confounding)
    } else {
      exclusive[[ct]] <- original
    }
  }
  exclusive
}

## Section 9.4: Run Extended Analysis
extended_markers <- define_extended_markers()
exclusive_markers <- create_exclusive_sets(extended_markers)

microglia_covariate <- list(
  Microglia = unique(c(extended_markers$Microglia_Activated,
                       extended_markers$Microglia_Homeostatic))
)

t_cell_types <- c("T_Cell_Extended", "T_Cell_CD4", "T_Cell_CD8", 
                  "TCR_Signaling", "NK_Cell")

results_unconditional <- data.frame()
results_conditional <- data.frame()
results_exclusive <- data.frame()

for (gwas_name in c("EOAD", "LOAD")) {
  magma_df <- if (gwas_name == "EOAD") eoad_magma else load_magma
  if (is.null(magma_df)) next
  
  for (ct in names(extended_markers)) {
    result <- run_conditional_regression(magma_df, extended_markers[[ct]], ct)
    if (!is.null(result)) {
      result$GWAS <- gwas_name
      result$Analysis <- "Unconditional"
      results_unconditional <- rbind(results_unconditional, as.data.frame(result))
    }
  }
  
  for (ct in t_cell_types) {
    result <- run_conditional_regression(
      magma_df, extended_markers[[ct]], ct,
      covariate_genes_list = microglia_covariate,
      covariate_names = "Microglia"
    )
    if (!is.null(result)) {
      result$GWAS <- gwas_name
      result$Analysis <- "Conditional_Microglia"
      results_conditional <- rbind(results_conditional, as.data.frame(result))
    }
    
    if (length(exclusive_markers[[ct]]) >= 10) {
      result <- run_conditional_regression(magma_df, exclusive_markers[[ct]], ct)
      if (!is.null(result)) {
        result$GWAS <- gwas_name
        result$Analysis <- "Exclusive"
        results_exclusive <- rbind(results_exclusive, as.data.frame(result))
      }
    }
  }
}

extended_results <- rbind(results_unconditional, results_conditional, results_exclusive)
if (nrow(extended_results) > 0) {
  fwrite(extended_results, "results/tables/Table_Extended_CellType_Analysis.csv")
}


# ==============================================================================
# Part 10: Graham et al. 2025 Correct Module Definitions
# ==============================================================================

## Section 10.1: Define Correct Graham Modules (From Methods)
define_graham_modules_correct <- function() {
  list(
    Human_Microglial_AD = c(
      "TREM2", "MS4A4A", "MS4A6A", "SPI1", "CD33", "INPP5D", "HLA-DRA",
      "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
      "TYROBP", "PLCG2", "ABI3", "CD68", "CSF1R", "C1QA", "C1QB", "C1QC",
      "ITGAM", "CX3CR1", "P2RY12", "HEXB", "GRN"
    ),
    Human_Oligodendrocyte_AD = c(
      "HNRNPA2B1", "TARDBP", "SLC45A3", "SLC20A2", "MBP", "PLP1", "MOG",
      "MAG", "MOBP", "CNP", "OLIG1", "OLIG2", "SOX10", "MYRF", "CLDN11",
      "GJC2", "ERMN", "FA2H", "UGT8", "GALC", "ASPA", "ENPP6", "QKI",
      "NKX2-2", "BCAS1", "NFASC", "CNTN2", "LPAR1", "GPR37"
    ),
    Mouse_ARM = c(
      "APOE", "RELB", "PTK2B", "TREM2", "TYROBP", "LPL", "CST7", "ITGAX",
      "CLEC7A", "SPP1", "GPNMB", "LGALS3", "CD9", "FABP5", "CTSB", "CTSD",
      "CTSL"
    ),
    Mouse_Phagolysosomal = c(
      "TOMM40", "TREM2", "GRN", "CTSD", "CTSB", "CTSL", "CTSS", "CTSZ",
      "LAMP1", "LAMP2", "ATP6V0A1", "ATP6V0D1", "ATP6V1A", "ATP6V1B2",
      "ATP6V1C1", "ATP6V1E1", "ATP6V1G1", "ATP6V1H", "PSAP", "NPC1", "NPC2",
      "LIPA", "HEXA", "HEXB", "GBA", "GALC", "ARSA", "ARSB", "NAGLU",
      "GUSB", "HGSNAT", "SGSH", "IDS", "IDUA", "GLB1", "NEU1", "FUCA1",
      "MAN2B1", "MANBA", "AGA", "ASAH1", "SMPD1", "PPT1", "TPP1", "CLN3"
    ),
    Mouse_HM2_Longevity = c(
      "STAT3", "CASP8", "FADS1", "FADS2", "ELOVL2", "ELOVL5", "SCD",
      "ACSL1", "ACSL3", "ACSL4", "ACSL5", "ACSL6", "FASN", "SREBF1",
      "SREBF2", "PPARG", "PPARA", "RXRA", "NR1H3", "NR1H2", "ABCA1",
      "ABCG1", "APOE", "CLU", "LDLR", "LRP1", "SCARB1", "CD36", "MSR1",
      "SOAT1", "HMGCR", "HMGCS1", "MVK", "FDFT1", "SQLE", "CYP51A1",
      "DHCR7", "DHCR24"
    )
  )
}

## Section 10.2: Test Correct Graham Module Enrichment
graham_modules_correct <- define_graham_modules_correct()

test_graham_correct <- function(magma_df, modules) {
  if (is.null(magma_df)) return(NULL)
  
  results <- lapply(names(modules), function(mod_name) {
    mod_genes <- toupper(modules[[mod_name]])
    magma_genes <- toupper(magma_df$SYMBOL)
    sig_genes <- magma_genes[magma_df$PVALUE < 0.05]
    
    overlap <- intersect(mod_genes, sig_genes)
    mod_in_magma <- intersect(mod_genes, magma_genes)
    
    if (length(mod_in_magma) < 5) return(NULL)
    
    contingency <- matrix(c(
      length(overlap),
      length(setdiff(mod_in_magma, sig_genes)),
      length(setdiff(sig_genes, mod_in_magma)),
      length(setdiff(magma_genes, union(mod_in_magma, sig_genes)))
    ), nrow = 2)
    
    fisher_res <- fisher.test(contingency, alternative = "greater")
    
    data.frame(
      Module = mod_name,
      N_Module = length(mod_in_magma),
      N_Overlap = length(overlap),
      OR = fisher_res$estimate,
      P_value = fisher_res$p.value,
      Overlap_Genes = paste(overlap, collapse = ";")
    )
  })
  
  results_df <- do.call(rbind, results[!sapply(results, is.null)])
  if (!is.null(results_df) && nrow(results_df) > 0) {
    results_df$FDR <- p.adjust(results_df$P_value, method = "BH")
  }
  results_df
}

eoad_graham_correct <- test_graham_correct(eoad_magma, graham_modules_correct)
load_graham_correct <- test_graham_correct(load_magma, graham_modules_correct)

if (!is.null(eoad_graham_correct) && !is.null(load_graham_correct)) {
  graham_correct_comparison <- merge(
    eoad_graham_correct[, c("Module", "N_Overlap", "OR", "P_value", "FDR")],
    load_graham_correct[, c("Module", "N_Overlap", "OR", "P_value", "FDR")],
    by = "Module", suffixes = c("_EOAD", "_LOAD")
  )
  graham_correct_comparison$log2OR_Diff <- log2(graham_correct_comparison$OR_EOAD / 
                                                  graham_correct_comparison$OR_LOAD)
  graham_correct_comparison$log2OR_Diff[!is.finite(graham_correct_comparison$log2OR_Diff)] <- 0
  
  fwrite(graham_correct_comparison, "results/tables/Table_Graham_Correct_Module_Comparison.csv")
}

## Section 10.3: Visualize Correct Graham Modules
if (exists("graham_correct_comparison")) {
  plot_data <- graham_correct_comparison
  plot_data$Module <- gsub("_", "\n", plot_data$Module)
  
  plot_long <- data.frame(
    Module = rep(plot_data$Module, 2),
    GWAS = rep(c("EOAD", "LOAD"), each = nrow(plot_data)),
    log2OR = c(log2(plot_data$OR_EOAD), log2(plot_data$OR_LOAD)),
    FDR = c(plot_data$FDR_EOAD, plot_data$FDR_LOAD)
  )
  plot_long$Sig <- ifelse(plot_long$FDR < 0.05, "*", "")
  plot_long$log2OR[!is.finite(plot_long$log2OR)] <- 0
  
  p_graham_correct <- ggplot(plot_long, aes(x = Module, y = log2OR, fill = GWAS)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = Sig, y = log2OR + sign(log2OR) * 0.1),
              position = position_dodge(0.8), size = 5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5")) +
    labs(title = "Graham et al. Module Enrichment",
         subtitle = "Human and Mouse AD-Associated Modules",
         x = "", y = "log2(Odds Ratio)") +
    theme_nc +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  
  ggsave("results/figures/Figure_Graham_Correct_Modules.pdf", p_graham_correct,
         width = 10, height = 6, dpi = 300)
}


# ==============================================================================
# Part 11: EOAD-Specific vs LOAD-Specific ORA Analysis
# ==============================================================================

## Section 11.1: Identify EOAD-Specific and LOAD-Specific Genes
identify_specific_genes <- function(eoad_magma, load_magma, p_threshold = 0.05) {
  if (is.null(eoad_magma) || is.null(load_magma)) return(NULL)
  
  eoad_df <- data.frame(SYMBOL = toupper(eoad_magma$SYMBOL), P_EOAD = eoad_magma$PVALUE)
  load_df <- data.frame(SYMBOL = toupper(load_magma$SYMBOL), P_LOAD = load_magma$PVALUE)
  
  merged <- merge(eoad_df, load_df, by = "SYMBOL", all = TRUE)
  
  list(
    eoad_specific = merged$SYMBOL[merged$P_EOAD < p_threshold & 
                                    (is.na(merged$P_LOAD) | merged$P_LOAD >= p_threshold)],
    load_specific = merged$SYMBOL[merged$P_LOAD < p_threshold & 
                                    (is.na(merged$P_EOAD) | merged$P_EOAD >= p_threshold)],
    shared = merged$SYMBOL[merged$P_EOAD < p_threshold & merged$P_LOAD < p_threshold],
    background = merged$SYMBOL
  )
}

specific_genes <- identify_specific_genes(eoad_magma, load_magma)

## Section 11.2: Run ORA for Specific Gene Sets
if (!is.null(specific_genes)) {
  cat("\nGene set sizes:\n")
  cat("  EOAD-specific:", length(specific_genes$eoad_specific), "\n")
  cat("  LOAD-specific:", length(specific_genes$load_specific), "\n")
  cat("  Shared:", length(specific_genes$shared), "\n")
  
  if (length(specific_genes$eoad_specific) >= 20) {
    eoad_specific_ora <- run_ora_enrichment(
      specific_genes$eoad_specific, 
      specific_genes$background, 
      "BP"
    )
    if (!is.null(eoad_specific_ora) && nrow(as.data.frame(eoad_specific_ora)) > 0) {
      fwrite(as.data.frame(eoad_specific_ora), 
             "results/tables/Table_EOAD_Specific_ORA.csv")
    }
  }
  
  if (length(specific_genes$load_specific) >= 20) {
    load_specific_ora <- run_ora_enrichment(
      specific_genes$load_specific, 
      specific_genes$background, 
      "BP"
    )
    if (!is.null(load_specific_ora) && nrow(as.data.frame(load_specific_ora)) > 0) {
      fwrite(as.data.frame(load_specific_ora), 
             "results/tables/Table_LOAD_Specific_ORA.csv")
    }
  }
  
  if (length(specific_genes$shared) >= 20) {
    shared_ora <- run_ora_enrichment(
      specific_genes$shared, 
      specific_genes$background, 
      "BP"
    )
    if (!is.null(shared_ora) && nrow(as.data.frame(shared_ora)) > 0) {
      fwrite(as.data.frame(shared_ora), 
             "results/tables/Table_Shared_Genes_ORA.csv")
    }
  }
}

## Section 11.3: Compare Specific Pathways
compare_specific_ora <- function(eoad_ora, load_ora) {
  if (is.null(eoad_ora) || is.null(load_ora)) return(NULL)
  
  eoad_df <- as.data.frame(eoad_ora)
  load_df <- as.data.frame(load_ora)
  
  if (nrow(eoad_df) == 0 || nrow(load_df) == 0) return(NULL)
  
  all_terms <- union(eoad_df$Description, load_df$Description)
  
  comparison <- data.frame(
    Term = all_terms,
    EOAD_Count = sapply(all_terms, function(t) {
      idx <- which(eoad_df$Description == t)
      if (length(idx) > 0) eoad_df$Count[idx[1]] else 0
    }),
    EOAD_padj = sapply(all_terms, function(t) {
      idx <- which(eoad_df$Description == t)
      if (length(idx) > 0) eoad_df$p.adjust[idx[1]] else 1
    }),
    LOAD_Count = sapply(all_terms, function(t) {
      idx <- which(load_df$Description == t)
      if (length(idx) > 0) load_df$Count[idx[1]] else 0
    }),
    LOAD_padj = sapply(all_terms, function(t) {
      idx <- which(load_df$Description == t)
      if (length(idx) > 0) load_df$p.adjust[idx[1]] else 1
    })
  )
  
  comparison$Specificity <- ifelse(
    comparison$EOAD_padj < 0.05 & comparison$LOAD_padj >= 0.05, "EOAD-enriched",
    ifelse(comparison$LOAD_padj < 0.05 & comparison$EOAD_padj >= 0.05, "LOAD-enriched",
           ifelse(comparison$EOAD_padj < 0.05 & comparison$LOAD_padj < 0.05, "Both", "Neither"))
  )
  
  comparison
}

if (exists("eoad_specific_ora") && exists("load_specific_ora")) {
  ora_comparison <- compare_specific_ora(eoad_specific_ora, load_specific_ora)
  if (!is.null(ora_comparison)) {
    fwrite(ora_comparison, "results/tables/Table_Specific_ORA_Comparison.csv")
  }
}

## Section 11.4: Visualize Specific Gene Analysis
if (!is.null(specific_genes)) {
  venn_data <- list(
    EOAD = specific_genes$eoad_specific,
    LOAD = specific_genes$load_specific,
    Shared = specific_genes$shared
  )
  
  gene_counts <- data.frame(
    Category = c("EOAD-specific", "LOAD-specific", "Shared"),
    Count = c(length(specific_genes$eoad_specific),
              length(specific_genes$load_specific),
              length(specific_genes$shared))
  )
  gene_counts$Category <- factor(gene_counts$Category, 
                                  levels = c("EOAD-specific", "Shared", "LOAD-specific"))
  
  p_specific <- ggplot(gene_counts, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = Count), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("EOAD-specific" = "#E64B35", 
                                  "Shared" = "#7E6148",
                                  "LOAD-specific" = "#4DBBD5")) +
    labs(title = "MAGMA Significant Genes (P < 0.05)",
         subtitle = "EOAD-specific vs LOAD-specific",
         x = "", y = "Number of Genes") +
    theme_nc +
    theme(legend.position = "none")
  
  ggsave("results/figures/Figure_Specific_Genes_Barplot.pdf", p_specific,
         width = 6, height = 5, dpi = 300)
}


# ==============================================================================
# Part 12: Extended Visualization
# ==============================================================================

## Section 12.1: Conditional Regression Forest Plot
if (exists("extended_results") && nrow(extended_results) > 0) {
  cond_data <- extended_results[extended_results$Analysis == "Conditional_Microglia", ]
  
  if (nrow(cond_data) > 0) {
    cond_data$Lower <- cond_data$Beta - 1.96 * cond_data$SE
    cond_data$Upper <- cond_data$Beta + 1.96 * cond_data$SE
    cond_data$Sig <- ifelse(cond_data$P_Value < 0.005, "Significant", "NS")
    
    p_forest <- ggplot(cond_data, aes(x = Beta, y = CellType, color = GWAS)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2,
                     position = position_dodge(0.4)) +
      geom_point(aes(shape = Sig), size = 3, position = position_dodge(0.4)) +
      scale_color_manual(values = c("EOAD" = "#E64B35", "LOAD" = "#4DBBD5")) +
      scale_shape_manual(values = c("Significant" = 16, "NS" = 1)) +
      labs(title = "T Cell Enrichment Conditional on Microglia",
           subtitle = "Z_MAGMA = β₀ + β₁×T_Cell + β₂×Microglia + ε",
           x = "Beta (Enrichment)", y = "") +
      theme_nc
    
    ggsave("results/figures/Figure_Conditional_Regression_Forest.pdf", p_forest,
           width = 9, height = 5, dpi = 300)
  }
}

## Section 12.2: Heatmap of All Cell Type Results
if (exists("extended_results") && nrow(extended_results) > 0) {
  uncond_data <- extended_results[extended_results$Analysis == "Unconditional", ]
  
  if (nrow(uncond_data) > 0) {
    heatmap_matrix <- reshape2::dcast(uncond_data, CellType ~ GWAS, value.var = "Beta")
    rownames(heatmap_matrix) <- heatmap_matrix$CellType
    heatmap_matrix <- as.matrix(heatmap_matrix[, -1])
    
    pval_matrix <- reshape2::dcast(uncond_data, CellType ~ GWAS, value.var = "P_Value")
    rownames(pval_matrix) <- pval_matrix$CellType
    pval_matrix <- as.matrix(pval_matrix[, -1])
    
    sig_matrix <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
    sig_matrix[pval_matrix < 0.001] <- "***"
    sig_matrix[pval_matrix >= 0.001 & pval_matrix < 0.01] <- "**"
    sig_matrix[pval_matrix >= 0.01 & pval_matrix < 0.05] <- "*"
    
    col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#B2182B"))
    
    pdf("results/figures/Figure_Extended_CellType_Heatmap.pdf", width = 6, height = 8)
    ht <- Heatmap(
      heatmap_matrix,
      name = "Beta",
      col = col_fun,
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sig_matrix[i, j], x, y, gp = gpar(fontsize = 10))
      },
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      row_names_gp = gpar(fontsize = 10),
      column_names_gp = gpar(fontsize = 11),
      column_title = "Extended Cell Type Enrichment",
      heatmap_legend_param = list(title = "Beta")
    )
    draw(ht)
    dev.off()
  }
}


# ==============================================================================
# Part 13: Final Results Summary
# ==============================================================================

## Section 13.1: Comprehensive Summary
final_summary <- list(
  analysis_date = Sys.time(),
  
  wgcna = list(
    n_modules = if(exists("wgcna_results")) length(unique(wgcna_results$moduleColors)) - 1 else NA,
    n_genes = if(exists("wgcna_results")) ncol(wgcna_results$datExpr) else NA,
    age_positive = if(exists("wgcna_results")) length(wgcna_results$age_positive_modules) else NA,
    age_negative = if(exists("wgcna_results")) length(wgcna_results$age_negative_modules) else NA
  ),
  
  magma = list(
    eoad_total = if(exists("eoad_magma")) nrow(eoad_magma) else NA,
    eoad_sig = if(exists("eoad_magma")) sum(eoad_magma$FDR < 0.05, na.rm = TRUE) else NA,
    load_total = if(exists("load_magma")) nrow(load_magma) else NA,
    load_sig = if(exists("load_magma")) sum(load_magma$FDR < 0.05, na.rm = TRUE) else NA
  ),
  
  specific_genes = if(exists("specific_genes")) list(
    eoad_specific = length(specific_genes$eoad_specific),
    load_specific = length(specific_genes$load_specific),
    shared = length(specific_genes$shared)
  ) else NULL,
  
  graham_modules = if(exists("graham_correct_comparison")) list(
    eoad_enriched = sum(graham_correct_comparison$FDR_EOAD < 0.05, na.rm = TRUE),
    load_enriched = sum(graham_correct_comparison$FDR_LOAD < 0.05, na.rm = TRUE)
  ) else NULL,
  
  extended_celltype = if(exists("extended_results")) list(
    n_celltypes = length(unique(extended_results$CellType)),
    eoad_sig = sum(extended_results$P_Value < 0.005 & extended_results$GWAS == "EOAD", na.rm = TRUE),
    load_sig = sum(extended_results$P_Value < 0.005 & extended_results$GWAS == "LOAD", na.rm = TRUE)
  ) else NULL
)

saveRDS(final_summary, "results/Final_Analysis_Summary.rds")

## Section 13.2: Print Final Summary
cat("\n")
cat("================================================================================\n")
cat("  EOAD vs LOAD Pathway Discovery Analysis - COMPLETE\n")
cat("================================================================================\n")
cat("\nAnalysis Summary:\n")
cat("  WGCNA: ", final_summary$wgcna$n_modules, " modules from ", 
    final_summary$wgcna$n_genes, " genes\n", sep = "")
cat("  MAGMA EOAD: ", final_summary$magma$eoad_sig, "/", final_summary$magma$eoad_total, 
    " FDR-significant genes\n", sep = "")
cat("  MAGMA LOAD: ", final_summary$magma$load_sig, "/", final_summary$magma$load_total, 
    " FDR-significant genes\n", sep = "")
if (!is.null(final_summary$specific_genes)) {
  cat("  Specific genes: EOAD=", final_summary$specific_genes$eoad_specific,
      ", LOAD=", final_summary$specific_genes$load_specific,
      ", Shared=", final_summary$specific_genes$shared, "\n", sep = "")
}
cat("\nOutput Files:\n")
cat("  Tables:\n")
cat("    - Table_WGCNA_Module_Assignment.csv\n")
cat("    - Table_Module_MAGMA_Enrichment.csv\n")
cat("    - Table_Graham_Correct_Module_Comparison.csv\n")
cat("    - Table_GSEA_Pathway_Comparison.csv\n")
cat("    - Table_CellType_MAGMA_Comparison.csv\n")
cat("    - Table_Extended_CellType_Analysis.csv\n")
cat("    - Table_EOAD_Specific_ORA.csv\n")
cat("    - Table_LOAD_Specific_ORA.csv\n")
cat("  Figures:\n")
cat("    - Figure_WGCNA_ModuleTrait_Heatmap.pdf\n")
cat("    - Figure_Graham_Correct_Modules.pdf\n")
cat("    - Figure_GSEA_Pathway_Comparison.pdf\n")
cat("    - Figure_CellType_Enrichment.pdf\n")
cat("    - Figure_Conditional_Regression_Forest.pdf\n")
cat("    - Figure_Extended_CellType_Heatmap.pdf\n")
cat("    - Figure_Specific_Genes_Barplot.pdf\n")
cat("\n")

sessionInfo()
