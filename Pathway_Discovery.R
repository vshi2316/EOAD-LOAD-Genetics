# ==============================================================================
# EOAD vs LOAD Data-Driven Pathway Discovery Analysis
# GitHub Version - Nature Communications Standard
# ==============================================================================
#
# Title: Genetic Heterogeneity Between Early-Onset and Late-Onset Alzheimer's Disease
#        Data-Driven Pathway Discovery and Cell Type Specificity Analysis
#
# Methods Alignment:
#   - MAGMA v1.10: Gene boundaries 35kb upstream, 10kb downstream
#   - 1000 Genomes Phase 3 EUR reference panel (503 individuals) for LD correction
#   - Bonferroni correction: P < 0.05/19,000 protein-coding genes
#   - WGCNA: GSE137313 mouse hippocampus RNA-seq (N=48, 4 age groups: 2,4,8,18 months)
#   - Signed network, biweight midcorrelation (bicor)
#   - Parameters: minModuleSize=30, mergeCutHeight=0.15, deepSplit=2, R²>0.85
#   - Gene filtering: variance > 25th percentile, mean > 20th percentile
#   - Mouse-to-human conversion: Ensembl BioMart (release 109) via biomaRt
#   - Graham et al. 2025: 5 curated modules with exact gene counts
#   - GSEA: clusterProfiler, gene set size 10-500, BH FDR correction
#   - Cell type markers: Methods-matched gene counts
#   - Conditional regression: Z_MAGMA = β₀ + β₁×T_Cell + β₂×Microglia + ε
#   - Multiple testing: Bonferroni P < 0.005 for 10 cell types
#
# ==============================================================================

rm(list = ls())
gc()

# ==============================================================================
# Part 1: Environment Setup and Package Loading
# ==============================================================================

cat("================================================================================\n")
cat("EOAD vs LOAD Data-Driven Pathway Discovery Analysis\n")
cat("Nature Communications Standard - GitHub Version\n")
cat("================================================================================\n\n")

# Required packages
packages_bioc <- c("WGCNA", "limma", "ComplexHeatmap", "biomaRt",
                   "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE")
packages_cran <- c("readxl", "dplyr", "tidyr", "ggplot2", "data.table",
                   "pheatmap", "corrplot", "igraph", "circlize", "tibble",
                   "RColorBrewer", "ggrepel", "cowplot", "gridExtra",
                   "VennDiagram", "UpSetR", "scales", "reshape2")

cat("Loading packages...\n")
for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

allowWGCNAThreads()

# Create output directories
dir.create("results/", showWarnings = FALSE)
dir.create("results/figures/", showWarnings = FALSE)
dir.create("results/tables/", showWarnings = FALSE)

cat("Environment setup complete.\n\n")

# ==============================================================================
# Part 2: WGCNA Co-expression Network Analysis (GSE137313)
# ==============================================================================
# Methods Reference:
#   Dataset: GSE137313 mouse hippocampus RNA-sequencing
#   Samples: N=48 (C57BL/6J wild-type and AD transgenic models)
#   Age groups: 2 months (young), 4 months (middle), 8 months (old), 18 months (aged)
#   Gene filtering: variance > 25th percentile, mean > 20th percentile
#   Network: Signed weighted co-expression network
#   Correlation: Biweight midcorrelation (bicor, robust to outliers)
#   Soft threshold: Power achieving scale-free topology R² > 0.85 (typically 6-12)
#   Module detection: minModuleSize=30, deepSplit=2, mergeCutHeight=0.15
#   Age association: |r| > 0.2, P < 0.1
#   Mouse-to-human conversion: Ensembl BioMart (release 109) via getLDS function
# ==============================================================================

cat("========================================\n")
cat("Part 2: WGCNA Co-expression Network Analysis\n")
cat("========================================\n\n")

## Section 2.1: Load and Preprocess GSE137313 Data
load_gse137313_data <- function(expr_file, pathology_file) {
  cat("  Loading GSE137313 data...\n")
  
  # Read expression data
  expr_data <- read_excel(expr_file)
  expr_data <- as.data.frame(expr_data)
  expr_data <- expr_data[!duplicated(expr_data[,1]),]
  rownames(expr_data) <- NULL
  expr_data <- column_to_rownames(expr_data, var = "SampleName")
  colnames(expr_data) <- expr_data[1,]
  colnames(expr_data) <- gsub("RNA_", "", colnames(expr_data))
  
  # Read pathology/sample info
  pathology_data <- read.csv(pathology_file)
  pathology_data$Sample.ID <- gsub("_133_2", "", pathology_data$Sample.ID)
  expr_data <- expr_data[, colnames(expr_data) %in% pathology_data$Sample.ID]
  
  rownames(pathology_data) <- pathology_data[,5]
  allTraits <- pathology_data[,-c(1,2,4,6,7,10)]
  allTraits <- allTraits[, -2]
  rownames(allTraits) <- rownames(pathology_data)
  
  sample_names <- colnames(expr_data)
  traitRows <- match(sample_names, rownames(allTraits))
  sample_info <- as.data.frame(allTraits[traitRows, ])
  
  # Batch assignment (Methods: batch correction)
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
  
  sample_info$batch <- NA
  sample_info[rownames(sample_info) %in% batch1_samples, ]$batch <- 1
  sample_info[rownames(sample_info) %in% batch2_samples, ]$batch <- 2
  sample_info[is.na(sample_info$batch), ]$batch <- 3
  
  sample_info$AGE_cat <- as.numeric(sample_info$AGE_cat)
  sample_info$MODELDIS <- as.factor(sample_info$MODELDIS)
  
  # Age group mapping (Methods: 2, 4, 8, 18 months)
  # AGE_cat values: 8=2months, 16=4months, 32=8months, 72=18months
  sample_info$Age_Months <- ifelse(sample_info$AGE_cat == 8, 2,
                                   ifelse(sample_info$AGE_cat == 16, 4,
                                          ifelse(sample_info$AGE_cat == 32, 8, 18)))
  
  sample_info$Age_Group <- factor(
    ifelse(sample_info$Age_Months == 2, "Young_2m",
           ifelse(sample_info$Age_Months == 4, "Middle_4m",
                  ifelse(sample_info$Age_Months == 8, "Old_8m", "Aged_18m"))),
    levels = c("Young_2m", "Middle_4m", "Old_8m", "Aged_18m"))
  
  sample_info$Disease_Model <- factor(
    ifelse(sample_info$MODELDIS == 1, "WT",
           ifelse(sample_info$MODELDIS %in% c(2,3), "AD_Model", "Other")),
    levels = c("WT", "AD_Model", "Other"))
  
  # Convert expression matrix
  expr_matrix <- as.matrix(expr_data[-1, ])
  expr_matrix <- apply(expr_matrix, 2, as.numeric)
  rownames(expr_matrix) <- rownames(expr_data)[-1]
  
  cat("    Samples: ", ncol(expr_matrix), "\n")
  cat("    Genes: ", nrow(expr_matrix), "\n")
  cat("    Age groups: ", paste(levels(sample_info$Age_Group), collapse = ", "), "\n")
  
  return(list(expr_matrix = expr_matrix, sample_info = sample_info))
}

## Section 2.2: Gene Filtering (Methods-matched)
filter_genes_for_wgcna <- function(expr_matrix) {
  # Methods: variance > 25th percentile, mean > 20th percentile
  # Retaining approximately 8,000-10,000 genes
  
  expr_var <- apply(expr_matrix, 1, var, na.rm = TRUE)
  expr_mean <- apply(expr_matrix, 1, mean, na.rm = TRUE)
  
  var_threshold <- quantile(expr_var, 0.25, na.rm = TRUE)
  mean_threshold <- quantile(expr_mean, 0.20, na.rm = TRUE)
  
  keep_genes <- expr_var > var_threshold &
    expr_mean > mean_threshold &
    !is.na(expr_var) & !is.na(expr_mean)
  
  expr_filtered <- expr_matrix[keep_genes, ]
  
  cat("  Gene filtering (Methods-matched):\n")
  cat("    Variance threshold (25th percentile): ", round(var_threshold, 4), "\n")
  cat("    Mean threshold (20th percentile): ", round(mean_threshold, 4), "\n")
  cat("    Genes retained: ", nrow(expr_filtered), " (target: 8,000-10,000)\n")
  
  return(expr_filtered)
}

## Section 2.3: WGCNA Network Construction (Methods-matched parameters)
run_wgcna_network <- function(expr_filtered, sample_info) {
  datExpr <- t(expr_filtered)
  
  # Quality control
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    cat("  Removed ", sum(!gsg$goodSamples), " samples and ", 
        sum(!gsg$goodGenes), " genes\n")
  }
  
  # Soft threshold selection (Methods: R² > 0.85, typically 6-12)
  cat("\n  Selecting soft threshold (target R² > 0.85)...\n")
  powers <- c(1:10, seq(12, 20, by = 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0,
                           networkType = "signed", corFnc = "bicor")
  
  # Find power achieving R² > 0.85
  r2_values <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  valid_powers <- which(r2_values > 0.85)
  
  if (length(valid_powers) > 0) {
    soft_power <- powers[valid_powers[1]]
    r2_achieved <- r2_values[valid_powers[1]]
  } else if (!is.na(sft$powerEstimate)) {
    soft_power <- sft$powerEstimate
    r2_achieved <- r2_values[which(powers == soft_power)]
  } else {
    soft_power <- 6  # Default within typical range
    r2_achieved <- r2_values[which(powers == 6)]
  }
  
  cat("    Selected power: ", soft_power, " (R² = ", round(r2_achieved, 3), ")\n")
  
  # Build network (Methods-matched parameters)
  cat("\n  Building signed weighted co-expression network...\n")
  cat("    Parameters: minModuleSize=30, deepSplit=2, mergeCutHeight=0.15\n")
  cat("    Correlation: biweight midcorrelation (bicor)\n")
  
  net <- blockwiseModules(
    datExpr,
    power = soft_power,
    TOMType = "signed",           # Methods: signed network
    minModuleSize = 30,           # Methods: minimum module size 30
    reassignThreshold = 0,
    mergeCutHeight = 0.15,        # Methods: modules with eigengene r > 0.85 merged
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 0,
    maxBlockSize = 20000,
    networkType = "signed",
    corType = "bicor",            # Methods: biweight midcorrelation
    deepSplit = 2                 # Methods: deep split 2
  )
  
  moduleColors <- labels2colors(net$colors)
  MEs <- net$MEs
  
  n_modules <- length(unique(moduleColors)) - 1  # Exclude grey
  cat("    Modules detected: ", n_modules, "\n")
  
  return(list(
    datExpr = datExpr,
    net = net,
    moduleColors = moduleColors,
    MEs = MEs,
    soft_power = soft_power,
    r2_achieved = r2_achieved
  ))
}

## Section 2.4: Module-Trait Correlation (Methods-matched)
calculate_module_trait_correlation <- function(wgcna_net, sample_info) {
  datExpr <- wgcna_net$datExpr
  MEs <- wgcna_net$MEs
  
  # Match samples
  trait_samples <- rownames(sample_info)[rownames(sample_info) %in% rownames(datExpr)]
  trait_data <- sample_info[trait_samples, ]
  trait_data <- trait_data[match(rownames(datExpr), rownames(trait_data)), ]
  
  # Create trait matrix
  trait_matrix <- data.frame(
    Age_Months = as.numeric(trait_data$Age_Months),
    Age_Young = as.numeric(trait_data$Age_Months == 2),
    Age_Middle = as.numeric(trait_data$Age_Months == 4),
    Age_Old = as.numeric(trait_data$Age_Months == 8),
    Age_Aged = as.numeric(trait_data$Age_Months == 18),
    AD_Model = as.numeric(trait_data$MODELDIS != 1),
    WT = as.numeric(trait_data$MODELDIS == 1),
    Batch = as.numeric(trait_data$batch)
  )
  rownames(trait_matrix) <- rownames(trait_data)
  
  # Calculate correlations (Methods: Pearson correlation, Student asymptotic P)
  MEs_ordered <- orderMEs(MEs)
  moduleTraitCor <- cor(MEs_ordered, trait_matrix, use = "pairwise.complete.obs")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  # Identify age-associated modules (Methods: |r| > 0.2, P < 0.1)
  age_positive_modules <- rownames(moduleTraitCor)[
    moduleTraitCor[, "Age_Months"] > 0.2 & moduleTraitPvalue[, "Age_Months"] < 0.1]
  age_negative_modules <- rownames(moduleTraitCor)[
    moduleTraitCor[, "Age_Months"] < -0.2 & moduleTraitPvalue[, "Age_Months"] < 0.1]
  
  cat("\n  Module-Trait Correlation (Age):\n")
  cat("    Age-positive modules (r > 0.2, P < 0.1): ", length(age_positive_modules), "\n")
  cat("    Age-negative modules (r < -0.2, P < 0.1): ", length(age_negative_modules), "\n")
  
  return(list(
    MEs = MEs_ordered,
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    age_positive_modules = age_positive_modules,
    age_negative_modules = age_negative_modules
  ))
}

## Section 2.5: Mouse-to-Human Gene Conversion (Methods: Ensembl BioMart release 109)
convert_mouse_to_human_biomart <- function(mouse_genes) {
  cat("\n  Converting mouse genes to human orthologs (Ensembl BioMart)...\n")
  
  tryCatch({
    # Methods: Ensembl BioMart via biomaRt, getLDS function
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    conversion <- getLDS(
      attributes = c("mgi_symbol"),
      filters = "mgi_symbol",
      values = mouse_genes,
      mart = mouse,
      attributesL = c("hgnc_symbol"),
      martL = human,
      uniqueRows = TRUE  # Methods: one-to-one ortholog mappings
    )
    
    colnames(conversion) <- c("Mouse_Gene", "Human_Gene")
    conversion <- conversion[conversion$Human_Gene != "", ]
    
    cat("    Input mouse genes: ", length(mouse_genes), "\n")
    cat("    Converted to human: ", nrow(conversion), "\n")
    cat("    Conversion rate: ", round(nrow(conversion)/length(mouse_genes)*100, 1), "%\n")
    
    return(conversion)
    
  }, error = function(e) {
    cat("    BioMart connection failed, using uppercase mapping as fallback\n")
    return(data.frame(
      Mouse_Gene = mouse_genes,
      Human_Gene = toupper(mouse_genes),
      stringsAsFactors = FALSE
    ))
  })
}

## Section 2.6: Execute WGCNA Analysis
# Check if data files exist
expr_file <- "GSE137313_geneQuantification_TPM_sampleinfo_MOD_DS_11May18.xlsx"
pathology_file <- "pathology updated for RNAseq Apr17.csv"

if (file.exists(expr_file) && file.exists(pathology_file)) {
  # Load data
  gse_data <- load_gse137313_data(expr_file, pathology_file)
  
  # Filter genes
  expr_filtered <- filter_genes_for_wgcna(gse_data$expr_matrix)
  
  # Run WGCNA
  wgcna_net <- run_wgcna_network(expr_filtered, gse_data$sample_info)
  
  # Calculate module-trait correlations
  module_trait <- calculate_module_trait_correlation(wgcna_net, gse_data$sample_info)
  
  # Convert mouse genes to human
  mouse_genes <- colnames(wgcna_net$datExpr)
  gene_conversion <- convert_mouse_to_human_biomart(mouse_genes)
  
  # Save WGCNA results
  wgcna_results <- list(
    datExpr = wgcna_net$datExpr,
    moduleColors = wgcna_net$moduleColors,
    MEs = module_trait$MEs,
    moduleTraitCor = module_trait$moduleTraitCor,
    moduleTraitPvalue = module_trait$moduleTraitPvalue,
    soft_power = wgcna_net$soft_power,
    r2_achieved = wgcna_net$r2_achieved,
    age_positive_modules = module_trait$age_positive_modules,
    age_negative_modules = module_trait$age_negative_modules,
    gene_conversion = gene_conversion
  )
  
  # Save module assignment table
  module_assignment <- data.frame(
    Mouse_Gene = colnames(wgcna_net$datExpr),
    Module_Color = wgcna_net$moduleColors
  )
  module_assignment <- merge(module_assignment, gene_conversion, 
                             by.x = "Mouse_Gene", by.y = "Mouse_Gene", all.x = TRUE)
  fwrite(module_assignment, "results/tables/Table_WGCNA_Module_Assignment.csv")
  
  saveRDS(wgcna_results, "results/WGCNA_Results.rds")
  cat("\n  WGCNA results saved.\n")
  
} else {
  cat("  GSE137313 data files not found. Skipping WGCNA analysis.\n")
  wgcna_results <- NULL
}


# ==============================================================================
# Part 3: MAGMA Gene-Based Association Analysis
# ==============================================================================
# Methods Reference:
#   Software: MAGMA v1.10 (Multi-marker Analysis of GenoMic Annotation)
#   Gene boundaries: 35 kb upstream, 10 kb downstream of TSS/TES
#   LD correction: 1000 Genomes Project Phase 3 EUR (503 individuals)
#   Model: SNP-wise mean model (mean chi-square with LD correction)
#   Significance: Bonferroni P < 0.05/19,000 protein-coding genes
#   Suggestive: P < 1 × 10⁻³ for exploratory analyses
# ==============================================================================

cat("\n========================================\n")
cat("Part 3: MAGMA Gene-Based Association Analysis\n")
cat("========================================\n\n")

## Section 3.1: Read MAGMA Results
read_magma_results <- function(magma_file) {
  if (!file.exists(magma_file)) {
    warning(paste("MAGMA file not found:", magma_file))
    return(NULL)
  }
  
  cat("  Loading: ", magma_file, "\n")
  magma <- fread(magma_file)
  
  # Standardize column names
  if ("P" %in% colnames(magma)) magma$PVALUE <- magma$P
  if ("GENE" %in% colnames(magma) && !("SYMBOL" %in% colnames(magma))) {
    magma$SYMBOL <- magma$GENE
  }
  
  # Calculate Z-statistic (Methods: for gene-property analysis)
  if (!"ZSTAT" %in% colnames(magma)) {
    magma$ZSTAT <- qnorm(1 - magma$PVALUE/2)
    # Handle P > 0.5 (negative Z)
    magma$ZSTAT[magma$PVALUE > 0.5] <- -abs(magma$ZSTAT[magma$PVALUE > 0.5])
  }
  
  # Multiple testing correction
  # Methods: Bonferroni P < 0.05/19,000 protein-coding genes
  n_genes <- 19000
  magma$Bonferroni <- p.adjust(magma$PVALUE, method = "bonferroni")
  magma$Bonferroni_Sig <- magma$PVALUE < (0.05 / n_genes)
  magma$FDR <- p.adjust(magma$PVALUE, method = "BH")
  
  # Suggestive significance (Methods: P < 1e-3)
  magma$Suggestive_Sig <- magma$PVALUE < 1e-3
  
  cat("    Total genes: ", nrow(magma), "\n")
  cat("    Bonferroni significant (P < ", format(0.05/n_genes, digits = 2), "): ", 
      sum(magma$Bonferroni_Sig, na.rm = TRUE), "\n")
  cat("    Suggestive (P < 1e-3): ", sum(magma$Suggestive_Sig, na.rm = TRUE), "\n")
  
  return(magma)
}

## Section 3.2: Load EOAD and LOAD MAGMA Results
eoad_magma <- read_magma_results("EOAD_MAGMA.genes.out")
load_magma <- read_magma_results("LOAD_MAGMA.genes.out")

# Also try to load multivariate aging GWAS if available
aging_magma <- read_magma_results("Aging_MAGMA.genes.out")


# ==============================================================================
# Part 4: Graham et al. 2025 Module Enrichment Analysis
# ==============================================================================
# Methods Reference:
#   Source: Graham et al. 2025 - integrated human/mouse scRNA-seq with AD/longevity GWAS
#   Modules (exact gene counts from Supplementary Methods):
#     1. Human Microglial AD-Significant Module: 26 genes
#     2. Human Oligodendrocyte AD-Significant Module: 29 genes
#     3. Mouse Activated Response Microglia (ARM) Module: 17 genes
#     4. Mouse Phagolysosomal Module: 45 genes
#     5. Mouse HM2 Longevity-Significant Module: 38 genes
#   Enrichment test: One-sided Fisher's exact test
#   Background: 20,000 protein-coding genes
#   Multiple testing: Benjamini-Hochberg FDR (significance: FDR < 0.05)
# ==============================================================================

cat("\n========================================\n")
cat("Part 4: Graham et al. 2025 Module Enrichment\n")
cat("========================================\n\n")

## Section 4.1: Define Graham Modules (Methods-matched exact gene counts)
define_graham_modules_methods_matched <- function() {
  list(
    # Human Microglial AD-Significant Module (26 genes)
    # Including TREM2, MS4A4A, MS4A6A, SPI1, CD33, INPP5D, MHC class II genes
    Human_Microglial_AD = c(
      "TREM2", "MS4A4A", "MS4A6A", "SPI1", "CD33", "INPP5D",
      "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
      "TYROBP", "PLCG2", "ABI3", "CD68", "CSF1R",
      "C1QA", "C1QB", "C1QC",
      "ITGAM", "CX3CR1", "P2RY12", "HEXB", "GRN"
    ),  # 26 genes
    
    # Human Oligodendrocyte AD-Significant Module (29 genes)
    # Including HNRNPA2B1, TARDBP, SLC45A3, SLC20A2
    Human_Oligodendrocyte_AD = c(
      "HNRNPA2B1", "TARDBP", "SLC45A3", "SLC20A2",
      "MBP", "PLP1", "MOG", "MAG", "MOBP", "CNP",
      "OLIG1", "OLIG2", "SOX10", "MYRF", "CLDN11",
      "GJC2", "ERMN", "FA2H", "UGT8", "GALC",
      "ASPA", "ENPP6", "QKI", "NKX2-2", "BCAS1",
      "NFASC", "CNTN2", "LPAR1", "GPR37"
    ),  # 29 genes
    
    # Mouse Activated Response Microglia (ARM) Module (17 genes)
    # Including APOE, RELB, PTK2B
    Mouse_ARM = c(
      "APOE", "RELB", "PTK2B",
      "TREM2", "TYROBP", "LPL", "CST7", "ITGAX",
      "CLEC7A", "SPP1", "GPNMB", "LGALS3", "CD9",
      "FABP5", "CTSB", "CTSD", "CTSL"
    ),  # 17 genes
    
    # Mouse Phagolysosomal Module (45 genes)
    # Including TOMM40, TREM2, GRN
    Mouse_Phagolysosomal = c(
      "TOMM40", "TREM2", "GRN",
      "CTSD", "CTSB", "CTSL", "CTSS", "CTSZ",
      "LAMP1", "LAMP2",
      "ATP6V0A1", "ATP6V0D1", "ATP6V1A", "ATP6V1B2", "ATP6V1C1", 
      "ATP6V1E1", "ATP6V1G1", "ATP6V1H",
      "PSAP", "NPC1", "NPC2", "LIPA",
      "HEXA", "HEXB", "GBA", "GALC", "ARSA", "ARSB",
      "NAGLU", "GUSB", "HGSNAT", "SGSH", "IDS", "IDUA",
      "GLB1", "NEU1", "FUCA1", "MAN2B1", "MANBA",
      "AGA", "ASAH1", "SMPD1", "PPT1", "TPP1", "CLN3"
    ),  # 45 genes
    
    # Mouse HM2 Longevity-Significant Module (38 genes)
    # Including STAT3, CASP8, FADS1
    Mouse_HM2_Longevity = c(
      "STAT3", "CASP8", "FADS1", "FADS2",
      "ELOVL2", "ELOVL5", "SCD",
      "ACSL1", "ACSL3", "ACSL4", "ACSL5", "ACSL6",
      "FASN", "SREBF1", "SREBF2",
      "PPARG", "PPARA", "RXRA", "NR1H3", "NR1H2",
      "ABCA1", "ABCG1", "APOE", "CLU", "LDLR", "LRP1",
      "SCARB1", "CD36", "MSR1", "SOAT1",
      "HMGCR", "HMGCS1", "MVK", "FDFT1", "SQLE",
      "CYP51A1", "DHCR7", "DHCR24"
    )  # 38 genes
  )
}

## Section 4.2: Verify Gene Counts Match Methods
graham_modules <- define_graham_modules_methods_matched()

cat("  Graham et al. 2025 Module Gene Counts:\n")
cat("  ─────────────────────────────────────────────────────────────────\n")
expected_counts <- c(
  Human_Microglial_AD = 26,
  Human_Oligodendrocyte_AD = 29,
  Mouse_ARM = 17,
  Mouse_Phagolysosomal = 45,
  Mouse_HM2_Longevity = 38
)

for (mod in names(graham_modules)) {
  actual <- length(graham_modules[[mod]])
  expected <- expected_counts[mod]
  status <- ifelse(actual == expected, "✓", "✗")
  cat(sprintf("    %s %-25s: %d genes (expected %d)\n", status, mod, actual, expected))
}

## Section 4.3: Fisher's Exact Test for Module Enrichment (Methods-matched)
test_module_enrichment_fisher <- function(module_genes, magma_results, 
                                          background_n = 20000,
                                          p_threshold = 0.05) {
  if (is.null(magma_results)) return(NULL)
  
  module_genes <- toupper(module_genes)
  magma_genes <- toupper(magma_results$SYMBOL)
  sig_genes <- magma_genes[magma_results$PVALUE < p_threshold]
  
  # Genes in module that are in MAGMA results
  module_in_magma <- intersect(module_genes, magma_genes)
  if (length(module_in_magma) < 5) return(NULL)
  
  # Overlap between module and significant genes
  overlap <- intersect(module_in_magma, sig_genes)
  
  # 2x2 contingency table (Methods-matched)
  # Row 1: In module, Row 2: Not in module
  # Col 1: Significant, Col 2: Not significant
  a <- length(overlap)                                    # In module & significant
  b <- length(module_in_magma) - a                        # In module & not significant
  c <- length(sig_genes) - a                              # Not in module & significant
  d <- background_n - a - b - c                           # Not in module & not significant
  if (d < 0) d <- 0
  
  contingency <- matrix(c(a, c, b, d), nrow = 2, byrow = FALSE)
  
  # One-sided Fisher's exact test (Methods: testing enrichment)
  fisher_result <- fisher.test(contingency, alternative = "greater")
  
  return(list(
    n_module = length(module_in_magma),
    n_overlap = a,
    overlap_genes = overlap,
    odds_ratio = as.numeric(fisher_result$estimate),
    ci_lower = fisher_result$conf.int[1],
    ci_upper = fisher_result$conf.int[2],
    p_value = fisher_result$p.value
  ))
}

## Section 4.4: Run Graham Module Enrichment Analysis
graham_enrichment_results <- data.frame()

for (module_name in names(graham_modules)) {
  module_genes <- graham_modules[[module_name]]
  
  # Test EOAD
  if (!is.null(eoad_magma)) {
    eoad_result <- test_module_enrichment_fisher(module_genes, eoad_magma)
    if (!is.null(eoad_result)) {
      graham_enrichment_results <- rbind(graham_enrichment_results, data.frame(
        Module = module_name,
        GWAS = "EOAD",
        N_Module = eoad_result$n_module,
        N_Overlap = eoad_result$n_overlap,
        OR = eoad_result$odds_ratio,
        CI_Lower = eoad_result$ci_lower,
        CI_Upper = eoad_result$ci_upper,
        P_value = eoad_result$p_value,
        Overlap_Genes = paste(eoad_result$overlap_genes, collapse = ";")
      ))
    }
  }
  
  # Test LOAD
  if (!is.null(load_magma)) {
    load_result <- test_module_enrichment_fisher(module_genes, load_magma)
    if (!is.null(load_result)) {
      graham_enrichment_results <- rbind(graham_enrichment_results, data.frame(
        Module = module_name,
        GWAS = "LOAD",
        N_Module = load_result$n_module,
        N_Overlap = load_result$n_overlap,
        OR = load_result$odds_ratio,
        CI_Lower = load_result$ci_lower,
        CI_Upper = load_result$ci_upper,
        P_value = load_result$p_value,
        Overlap_Genes = paste(load_result$overlap_genes, collapse = ";")
      ))
    }
  }
}

# Apply BH FDR correction (Methods-matched)
if (nrow(graham_enrichment_results) > 0) {
  graham_enrichment_results$FDR <- p.adjust(graham_enrichment_results$P_value, method = "BH")
  graham_enrichment_results$Significant <- graham_enrichment_results$FDR < 0.05
  
  cat("\n  Graham Module Enrichment Results:\n")
  cat("  ─────────────────────────────────────────────────────────────────\n")
  
  for (gwas in c("EOAD", "LOAD")) {
    cat("\n  ", gwas, ":\n", sep = "")
    gwas_results <- graham_enrichment_results[graham_enrichment_results$GWAS == gwas, ]
    gwas_results <- gwas_results[order(gwas_results$P_value), ]
    
    for (i in 1:nrow(gwas_results)) {
      row <- gwas_results[i, ]
      sig_mark <- ifelse(row$FDR < 0.05, "**", ifelse(row$P_value < 0.05, "*", ""))
      cat(sprintf("    %-25s: OR=%5.2f, P=%.4f, FDR=%.4f %s\n",
                  row$Module, row$OR, row$P_value, row$FDR, sig_mark))
    }
  }
  
  fwrite(graham_enrichment_results, "results/tables/Table_Graham_Module_Enrichment.csv")
}


# ==============================================================================
# Part 5: GSEA and Data-Driven Pathway Discovery
# ==============================================================================
# Methods Reference:
#   Software: clusterProfiler (version 4.8)
#   Ranking: Genes ranked by MAGMA Z-statistics (descending)
#   Gene sets: GO Biological Process (GO-BP), KEGG pathways
#   Gene set size: 10-500 genes
#   P-value cutoff: 0.25 (GSEA standard)
#   Multiple testing: Benjamini-Hochberg FDR (significance: FDR < 0.05)
#   ORA: Hypergeometric tests for EOAD-specific, LOAD-specific, shared genes
# ==============================================================================

cat("\n========================================\n")
cat("Part 5: GSEA and Data-Driven Pathway Discovery\n")
cat("========================================\n\n")

## Section 5.1: Prepare Ranked Gene List for GSEA
prepare_gsea_genelist <- function(magma_results, name) {
  if (is.null(magma_results)) return(NULL)
  
  cat("  Preparing ranked gene list for ", name, "...\n", sep = "")
  
  # Get gene symbols and Z-statistics
  gene_symbols <- unique(magma_results$SYMBOL[!is.na(magma_results$SYMBOL)])
  
  # Convert to Entrez IDs (required for clusterProfiler)
  gene_mapping <- tryCatch({
    bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
         OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)
  
  if (is.null(gene_mapping) || nrow(gene_mapping) < 100) {
    cat("    Gene ID conversion failed\n")
    return(NULL)
  }
  
  # Merge with MAGMA results
  merged <- merge(gene_mapping, magma_results, 
                  by.x = "SYMBOL", by.y = "SYMBOL", all.x = TRUE)
  merged <- merged[!is.na(merged$ZSTAT), ]
  
  # Create ranked gene list (Methods: ranked by MAGMA Z-statistics, descending)
  gene_list <- merged$ZSTAT
  names(gene_list) <- merged$ENTREZID
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  cat("    Ranked genes: ", length(gene_list), "\n")
  
  return(gene_list)
}

## Section 5.2: Run GSEA Analysis (Methods-matched parameters)
run_gsea_analysis <- function(gene_list, name) {
  if (is.null(gene_list) || length(gene_list) < 1000) {
    cat("  Insufficient genes for GSEA: ", name, "\n")
    return(NULL)
  }
  
  cat("\n  Running GSEA for ", name, "...\n", sep = "")
  results <- list()
  
  # GO Biological Process GSEA (Methods-matched)
  tryCatch({
    gsea_bp <- gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      minGSSize = 10,           # Methods: minimum 10 genes
      maxGSSize = 500,          # Methods: maximum 500 genes
      pvalueCutoff = 0.25,      # Methods: P-value cutoff 0.25
      verbose = FALSE
    )
    
    if (!is.null(gsea_bp) && nrow(gsea_bp@result) > 0) {
      results$GO_BP <- gsea_bp
      n_sig <- sum(gsea_bp@result$p.adjust < 0.05)
      cat("    GO BP: ", nrow(gsea_bp@result), " pathways, ", n_sig, " FDR < 0.05\n")
    }
  }, error = function(e) cat("    GO BP GSEA failed: ", e$message, "\n"))
  
  # GO Molecular Function GSEA
  tryCatch({
    gsea_mf <- gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      ont = "MF",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.25,
      verbose = FALSE
    )
    
    if (!is.null(gsea_mf) && nrow(gsea_mf@result) > 0) {
      results$GO_MF <- gsea_mf
      n_sig <- sum(gsea_mf@result$p.adjust < 0.05)
      cat("    GO MF: ", nrow(gsea_mf@result), " pathways, ", n_sig, " FDR < 0.05\n")
    }
  }, error = function(e) cat("    GO MF GSEA failed\n"))
  
  # KEGG GSEA
  tryCatch({
    gsea_kegg <- gseKEGG(
      geneList = gene_list,
      organism = "hsa",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.25,
      verbose = FALSE
    )
    
    if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
      results$KEGG <- gsea_kegg
      n_sig <- sum(gsea_kegg@result$p.adjust < 0.05)
      cat("    KEGG: ", nrow(gsea_kegg@result), " pathways, ", n_sig, " FDR < 0.05\n")
    }
  }, error = function(e) cat("    KEGG GSEA failed\n"))
  
  return(results)
}

## Section 5.3: Run GSEA for EOAD and LOAD
eoad_genelist <- prepare_gsea_genelist(eoad_magma, "EOAD")
load_genelist <- prepare_gsea_genelist(load_magma, "LOAD")

eoad_gsea <- run_gsea_analysis(eoad_genelist, "EOAD")
load_gsea <- run_gsea_analysis(load_genelist, "LOAD")

## Section 5.4: Save GSEA Results
if (!is.null(eoad_gsea$GO_BP)) {
  fwrite(eoad_gsea$GO_BP@result, "results/tables/Table_EOAD_GSEA_GO_BP.csv")
}
if (!is.null(eoad_gsea$KEGG)) {
  fwrite(eoad_gsea$KEGG@result, "results/tables/Table_EOAD_GSEA_KEGG.csv")
}
if (!is.null(load_gsea$GO_BP)) {
  fwrite(load_gsea$GO_BP@result, "results/tables/Table_LOAD_GSEA_GO_BP.csv")
}
if (!is.null(load_gsea$KEGG)) {
  fwrite(load_gsea$KEGG@result, "results/tables/Table_LOAD_GSEA_KEGG.csv")
}

## Section 5.5: EOAD-Specific vs LOAD-Specific ORA Analysis
# Methods: Over-representation analysis with hypergeometric tests
# Gene sets: EOAD-specific (P<0.05 in EOAD only), LOAD-specific, Shared

cat("\n  Identifying EOAD-specific, LOAD-specific, and Shared genes...\n")

identify_specific_genes <- function(eoad_magma, load_magma, p_threshold = 0.05) {
  if (is.null(eoad_magma) || is.null(load_magma)) return(NULL)
  
  eoad_sig <- toupper(eoad_magma$SYMBOL[eoad_magma$PVALUE < p_threshold])
  load_sig <- toupper(load_magma$SYMBOL[load_magma$PVALUE < p_threshold])
  
  eoad_sig <- eoad_sig[!is.na(eoad_sig)]
  load_sig <- load_sig[!is.na(load_sig)]
  
  eoad_specific <- setdiff(eoad_sig, load_sig)
  load_specific <- setdiff(load_sig, eoad_sig)
  shared <- intersect(eoad_sig, load_sig)
  
  cat("    EOAD-specific: ", length(eoad_specific), " genes\n")
  cat("    LOAD-specific: ", length(load_specific), " genes\n")
  cat("    Shared: ", length(shared), " genes\n")
  
  return(list(
    eoad_specific = eoad_specific,
    load_specific = load_specific,
    shared = shared
  ))
}

specific_genes <- identify_specific_genes(eoad_magma, load_magma)

## Section 5.6: Run ORA with Hypergeometric Test (Methods-matched)
run_ora_hypergeometric <- function(gene_list, name) {
  if (length(gene_list) < 10) {
    cat("    ", name, ": Too few genes (", length(gene_list), "), skipping\n")
    return(NULL)
  }
  
  # Convert to Entrez IDs
  gene_mapping <- tryCatch({
    bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)
  
  if (is.null(gene_mapping) || nrow(gene_mapping) < 10) {
    cat("    ", name, ": Gene ID conversion failed\n")
    return(NULL)
  }
  
  entrez_ids <- unique(gene_mapping$ENTREZID)
  cat("    ", name, ": ", length(entrez_ids), " genes with Entrez IDs\n")
  
  results <- list()
  
  # GO BP ORA (Methods: hypergeometric test, BH FDR)
  tryCatch({
    go_bp <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500
    )
    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      results$GO_BP <- go_bp
      n_sig <- sum(go_bp@result$p.adjust < 0.05)
      cat("      GO BP: ", n_sig, " significant pathways\n")
    }
  }, error = function(e) NULL)
  
  # KEGG ORA
  tryCatch({
    kegg <- enrichKEGG(
      gene = entrez_ids,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500
    )
    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      results$KEGG <- kegg
      n_sig <- sum(kegg@result$p.adjust < 0.05)
      cat("      KEGG: ", n_sig, " significant pathways\n")
    }
  }, error = function(e) NULL)
  
  return(results)
}

## Section 5.7: Run ORA for Each Gene Set
cat("\n  Running ORA for specific gene sets...\n")

if (!is.null(specific_genes)) {
  eoad_specific_ora <- run_ora_hypergeometric(specific_genes$eoad_specific, "EOAD-specific")
  load_specific_ora <- run_ora_hypergeometric(specific_genes$load_specific, "LOAD-specific")
  shared_ora <- run_ora_hypergeometric(specific_genes$shared, "Shared")
  
  # Save ORA results
  if (!is.null(eoad_specific_ora$GO_BP)) {
    fwrite(eoad_specific_ora$GO_BP@result, "results/tables/Table_EOAD_Specific_ORA_GO_BP.csv")
  }
  if (!is.null(load_specific_ora$GO_BP)) {
    fwrite(load_specific_ora$GO_BP@result, "results/tables/Table_LOAD_Specific_ORA_GO_BP.csv")
  }
  if (!is.null(shared_ora$GO_BP)) {
    fwrite(shared_ora$GO_BP@result, "results/tables/Table_Shared_ORA_GO_BP.csv")
  }
}


# ==============================================================================
# Part 6: Cell Type MAGMA Gene-Property Analysis
# ==============================================================================
# Methods Reference:
#   Analysis: MAGMA gene-property analysis
#   Model: Linear regression Z_MAGMA = β₀ + β₁ × CellType + ε
#   Interpretation: Positive β₁ indicates enrichment for that cell type
#   
#   Cell type marker gene sets (exact counts from Supplementary Methods):
#     - Pan-T cell markers: 64 genes
#     - CD4+ T helper cells: 42 genes
#     - CD8+ cytotoxic T cells: 40 genes
#     - TCR signaling pathway: 37 genes
#     - Natural killer cells: 32 genes
#     - Activated microglia: 47 genes
#     - Homeostatic microglia: 27 genes
#     - Amyloid-beta clearance: 39 genes
#     - APP metabolism: 23 genes
#   
#   Multiple testing: Bonferroni P < 0.005 for 10 cell types
# ==============================================================================

cat("\n========================================\n")
cat("Part 6: Cell Type MAGMA Gene-Property Analysis\n")
cat("========================================\n\n")

## Section 6.1: Define Cell Type Markers (Methods-matched exact gene counts)
define_cell_type_markers_methods <- function() {
  list(
    # Pan-T cell markers (64 genes)
    T_Cell_Pan = c(
      # Core T cell markers
      "CD3D", "CD3E", "CD3G", "CD247", "CD2", "CD5", "CD7", "CD27", "CD28",
      "TRAC", "TRBC1", "TRBC2", "TRAT1",
      # Transcription factors
      "TCF7", "LEF1", "GATA3", "TBX21", "EOMES", "BCL11B", "RUNX3", "IKZF1", "IKZF3",
      # Signaling molecules
      "LCK", "ZAP70", "LAT", "ITK", "FYN", "PLCG1", "VAV1", "PRKCQ", "CARD11",
      "NFATC1", "NFATC2", "NFKB1", "NFKBIA", "REL",
      # Cytokine receptors
      "IL7R", "IL2RA", "IL2RB", "IL2RG", "IL15RA", "IL21R", "IL23R", "IFNGR1",
      # Chemokine receptors
      "CCR7", "CCR4", "CCR5", "CCR6", "CXCR3", "CXCR4", "CXCR5", "CXCR6",
      # Co-stimulatory/inhibitory
      "CTLA4", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "ICOS", "CD40LG",
      # Effector molecules
      "GZMA", "GZMB", "GZMK", "PRF1", "GNLY"
    ),  # 64 genes
    
    # CD4+ T helper cells (42 genes)
    T_Cell_CD4 = c(
      # Core CD4 markers
      "CD4", "IL7R", "CCR7", "SELL", "LEF1", "TCF7", "MAL",
      # Th1
      "TBX21", "IFNG", "TNF", "IL2", "CXCR3", "CCR5", "STAT4", "IL12RB1", "IL12RB2",
      # Th2
      "GATA3", "IL4", "IL5", "IL13", "CCR4", "STAT6", "IL4R",
      # Th17
      "RORC", "IL17A", "IL17F", "IL22", "IL23R", "CCR6", "STAT3",
      # Treg
      "FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2", "TIGIT",
      # Tfh
      "CXCR5", "BCL6", "ICOS", "PDCD1", "IL21", "SH2D1A"
    ),  # 42 genes
    
    # CD8+ cytotoxic T cells (40 genes)
    T_Cell_CD8 = c(
      # Core CD8 markers
      "CD8A", "CD8B", "EOMES", "TBX21", "RUNX3", "PRDM1", "ID2",
      # Effector molecules
      "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "NKG7", "FASLG", "IFNG",
      # Cytokines
      "TNF", "LTA", "CSF2", "IL2",
      # Chemokine receptors
      "CX3CR1", "CXCR3", "CCR5", "CXCR6", "CCR7",
      # Activation markers
      "KLRG1", "KLRD1", "KLRC1", "KLRK1", "CD69",
      # Co-inhibitory receptors
      "PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4",
      # Transcription factors
      "ZEB2", "ZNF683", "TOX", "TCF7"
    ),  # 40 genes
    
    # TCR signaling pathway (37 genes)
    TCR_Signaling = c(
      # TCR complex
      "CD3D", "CD3E", "CD3G", "CD247", "TRAC", "TRBC1", "TRBC2",
      # Proximal signaling
      "LCK", "FYN", "ZAP70", "LAT", "SLP76", "ITK",
      # Adaptor proteins
      "GRAP2", "GRB2", "VAV1", "NCK1",
      # PLC pathway
      "PLCG1", "PRKCQ", "RASGRP1",
      # MAPK pathway
      "MAP3K7", "MAP2K1", "MAPK1", "MAPK3",
      # NF-kB pathway
      "CARD11", "BCL10", "MALT1", "NFKB1", "RELA", "NFKBIA",
      # NFAT pathway
      "NFATC1", "NFATC2", "PPP3CA", "PPP3CB", "PPP3CC",
      # AP-1 pathway
      "FOS", "JUN"
    ),  # 37 genes
    
    # Natural killer cells (32 genes)
    NK_Cell = c(
      # Core NK markers
      "NCAM1", "NCR1", "NCR2", "NCR3", "FCGR3A", "KLRD1", "KLRF1",
      # KIR family
      "KIR2DL1", "KIR2DL3", "KIR3DL1", "KIR3DL2", "KIR2DS4",
      # KLR family
      "KLRC1", "KLRC2", "KLRB1", "KLRK1", "KLRG1",
      # Effector molecules
      "GZMB", "GZMA", "PRF1", "GNLY", "NKG7", "FASLG",
      # Cytokines
      "IFNG", "TNF", "XCL1", "CCL3", "CCL4", "CCL5",
      # Transcription factors
      "EOMES", "TBX21", "ID2"
    ),  # 32 genes
    
    # Activated/disease-associated microglia (47 genes)
    Microglia_Activated = c(
      # DAM/ARM core genes
      "TREM2", "TYROBP", "APOE", "LPL", "CST7", "ITGAX", "CLEC7A", "SPP1",
      "GPNMB", "LGALS3", "CD9", "FABP5", "CTSB", "CTSD", "CTSL", "CTSS",
      # Complement system
      "C1QA", "C1QB", "C1QC", "C3", "C3AR1",
      # Inflammatory
      "IL1B", "IL6", "TNF", "CCL2", "CCL3", "CCL4", "CXCL10",
      # Phagocytic receptors
      "CD68", "CD14", "FCGR1A", "FCGR2A", "MSR1", "CD36",
      # AD risk genes
      "MS4A4A", "MS4A6A", "INPP5D", "CD33", "PLCG2", "ABI3", "GRN",
      # Transcription factors
      "SPI1", "IRF8", "RUNX1",
      # Other
      "AIF1", "ITGAM", "CSF1R"
    ),  # 47 genes
    
    # Homeostatic microglia (27 genes)
    Microglia_Homeostatic = c(
      # Core homeostatic markers
      "P2RY12", "P2RY13", "TMEM119", "SALL1", "HEXB", "CST3", "SPARC",
      "CX3CR1", "SELPLG", "SIGLECH", "OLFML3", "GPR34", "TGFBR1",
      # Neuroprotective
      "BDNF", "IGF1", "TGFB1", "IL10",
      # Synaptic pruning
      "MERTK", "GAS6", "PROS1",
      # Metabolism
      "SLC2A5", "GLUL", "GPX1", "SOD1", "SOD2",
      # Other
      "FCRLS", "MEF2C"
    ),  # 27 genes
    
    # Amyloid-beta clearance (39 genes)
    Abeta_Clearance = c(
      # Abeta degrading enzymes
      "IDE", "MME", "ECE1", "ECE2", "ACE", "THOP1", "PREP",
      "MMP2", "MMP9", "ADAMTS4",
      # Receptor-mediated clearance
      "LRP1", "LRP2", "LDLR", "VLDLR", "SORL1", "SCARB1",
      # Transporters
      "ABCA1", "ABCA7", "ABCB1", "ABCG1", "ABCG2",
      # Phagocytic receptors
      "TREM2", "CD36", "SCARA1", "MSR1", "CD14", "TLR2", "TLR4",
      # Autophagy
      "BECN1", "ATG5", "ATG7", "SQSTM1", "TFEB",
      # Proteasome
      "PSMB5", "PSMB6", "PSMA1",
      # Chaperones
      "HSPA1A", "HSPA8", "CLU"
    ),  # 39 genes
    
    # APP metabolism (23 genes)
    APP_Metabolism = c(
      # Gamma-secretase
      "PSEN1", "PSEN2", "NCSTN", "APH1A", "APH1B", "PSENEN",
      # Beta-secretase
      "BACE1", "BACE2",
      # Alpha-secretase
      "ADAM10", "ADAM17",
      # APP family
      "APP", "APLP1", "APLP2",
      # APP trafficking
      "SORL1", "SORCS1", "LRP1", "BIN1", "PICALM", "CD2AP",
      # Endosomal sorting
      "RAB5A", "VPS35", "SNX27",
      # Other
      "APBB1"
    )  # 23 genes
  )
}

## Section 6.2: Verify Gene Counts Match Methods
cell_type_markers <- define_cell_type_markers_methods()

cat("  Cell Type Marker Gene Counts (Methods-matched):\n")
cat("  ─────────────────────────────────────────────────────────────────\n")

expected_counts <- c(
  T_Cell_Pan = 64, T_Cell_CD4 = 42, T_Cell_CD8 = 40, TCR_Signaling = 37,
  NK_Cell = 32, Microglia_Activated = 47, Microglia_Homeostatic = 27,
  Abeta_Clearance = 39, APP_Metabolism = 23
)

for (ct in names(cell_type_markers)) {
  actual <- length(cell_type_markers[[ct]])
  expected <- expected_counts[ct]
  status <- ifelse(actual == expected, "✓", "✗")
  cat(sprintf("    %s %-25s: %d genes (expected %d)\n", status, ct, actual, expected))
}

## Section 6.3: MAGMA Gene-Property Regression (Methods-matched)
run_gene_property_regression <- function(magma_results, marker_genes, cell_type_name) {
  if (is.null(magma_results) || !"ZSTAT" %in% colnames(magma_results)) return(NULL)
  
  df <- magma_results[!is.na(magma_results$SYMBOL) & !is.na(magma_results$ZSTAT), ]
  df$SYMBOL <- toupper(df$SYMBOL)
  marker_genes_upper <- toupper(marker_genes)
  
  # Binary indicator for cell type membership
  df$CellType <- as.numeric(df$SYMBOL %in% marker_genes_upper)
  
  n_in_marker <- sum(df$CellType)
  if (n_in_marker < 10) return(NULL)
  
  # Linear regression: Z_MAGMA = β₀ + β₁ × CellType + ε (Methods-matched)
  model <- lm(ZSTAT ~ CellType, data = df)
  model_summary <- summary(model)
  coef_table <- coef(model_summary)
  
  beta <- coef_table["CellType", 1]
  se <- coef_table["CellType", 2]
  t_stat <- coef_table["CellType", 3]
  p_value <- coef_table["CellType", 4]
  
  mean_z_in <- mean(df$ZSTAT[df$CellType == 1], na.rm = TRUE)
  mean_z_out <- mean(df$ZSTAT[df$CellType == 0], na.rm = TRUE)
  
  return(data.frame(
    CellType = cell_type_name,
    N_Genes = n_in_marker,
    Beta = beta,
    SE = se,
    T_Stat = t_stat,
    P_Value = p_value,
    Mean_Z_In = mean_z_in,
    Mean_Z_Out = mean_z_out
  ))
}

## Section 6.4: Run Gene-Property Analysis for All Cell Types
cat("\n  Running MAGMA Gene-Property Analysis...\n")

gene_property_results <- data.frame()

for (gwas_name in c("EOAD", "LOAD")) {
  magma_df <- if (gwas_name == "EOAD") eoad_magma else load_magma
  
  for (ct in names(cell_type_markers)) {
    result <- run_gene_property_regression(magma_df, cell_type_markers[[ct]], ct)
    if (!is.null(result)) {
      result$GWAS <- gwas_name
      gene_property_results <- rbind(gene_property_results, result)
    }
  }
}

# Apply Bonferroni correction for 10 cell types (Methods: P < 0.005)
n_cell_types <- 10
bonferroni_threshold <- 0.05 / n_cell_types  # P < 0.005

if (nrow(gene_property_results) > 0) {
  gene_property_results$Bonferroni_Sig <- gene_property_results$P_Value < bonferroni_threshold
  gene_property_results$FDR <- p.adjust(gene_property_results$P_Value, method = "BH")
  
  cat("\n  Results (Bonferroni threshold P < ", bonferroni_threshold, "):\n", sep = "")
  cat("  ─────────────────────────────────────────────────────────────────\n")
  
  for (gwas in c("EOAD", "LOAD")) {
    cat("\n  ", gwas, ":\n", sep = "")
    gwas_results <- gene_property_results[gene_property_results$GWAS == gwas, ]
    gwas_results <- gwas_results[order(gwas_results$P_Value), ]
    
    for (i in 1:nrow(gwas_results)) {
      row <- gwas_results[i, ]
      sig_mark <- ifelse(row$P_Value < 0.001, "***",
                         ifelse(row$P_Value < bonferroni_threshold, "**",
                                ifelse(row$P_Value < 0.05, "*", "")))
      cat(sprintf("    %-25s: β=%6.3f, P=%.2e, N=%d %s\n",
                  row$CellType, row$Beta, row$P_Value, row$N_Genes, sig_mark))
    }
  }
  
  fwrite(gene_property_results, "results/tables/Table_MAGMA_GeneProperty_CellType.csv")
}


# ==============================================================================
# Part 7: Conditional Regression Analysis
# ==============================================================================
# Methods Reference:
#   Purpose: Test independence of T cell associations from microglial signals
#   Model: Z_MAGMA = β₀ + β₁ × T_Cell + β₂ × Microglia + ε
#   Interpretation: Significant β₁ after controlling for microglia indicates
#                   independent T cell genetic association
#   Multiple testing: Bonferroni P < 0.005 for 10 cell types
# ==============================================================================

cat("\n========================================\n")
cat("Part 7: Conditional Regression Analysis\n")
cat("========================================\n\n")

## Section 7.1: Define Microglia Covariate Gene Set
microglia_all_genes <- unique(c(
  cell_type_markers$Microglia_Activated,
  cell_type_markers$Microglia_Homeostatic
))
cat("  Microglia covariate gene set: ", length(microglia_all_genes), " genes\n")

## Section 7.2: Conditional Regression Function (Methods-matched)
run_conditional_regression <- function(magma_results, target_genes, target_name,
                                       covariate_genes, covariate_name = "Microglia") {
  if (is.null(magma_results) || !"ZSTAT" %in% colnames(magma_results)) return(NULL)
  
  df <- magma_results[!is.na(magma_results$SYMBOL) & !is.na(magma_results$ZSTAT), ]
  df$SYMBOL <- toupper(df$SYMBOL)
  
  # Binary indicators
  df$Target <- as.numeric(df$SYMBOL %in% toupper(target_genes))
  df$Covariate <- as.numeric(df$SYMBOL %in% toupper(covariate_genes))
  
  n_target <- sum(df$Target)
  n_covariate <- sum(df$Covariate)
  if (n_target < 10) return(NULL)
  
  # Conditional regression: Z_MAGMA = β₀ + β₁ × Target + β₂ × Covariate + ε
  model <- lm(ZSTAT ~ Target + Covariate, data = df)
  model_summary <- summary(model)
  coef_table <- coef(model_summary)
  
  beta_target <- coef_table["Target", 1]
  se_target <- coef_table["Target", 2]
  t_target <- coef_table["Target", 3]
  p_target <- coef_table["Target", 4]
  
  beta_cov <- coef_table["Covariate", 1]
  se_cov <- coef_table["Covariate", 2]
  p_cov <- coef_table["Covariate", 4]
  
  return(data.frame(
    CellType = target_name,
    N_Target = n_target,
    N_Covariate = n_covariate,
    Beta_Target = beta_target,
    SE_Target = se_target,
    T_Target = t_target,
    P_Target = p_target,
    Beta_Covariate = beta_cov,
    SE_Covariate = se_cov,
    P_Covariate = p_cov,
    Covariate = covariate_name
  ))
}

## Section 7.3: Run Conditional Analysis for T Cell Types
t_cell_types <- c("T_Cell_Pan", "T_Cell_CD4", "T_Cell_CD8", "TCR_Signaling", "NK_Cell")

conditional_results <- data.frame()

cat("\n  Running conditional regression (T Cell | Microglia)...\n")

for (gwas_name in c("EOAD", "LOAD")) {
  magma_df <- if (gwas_name == "EOAD") eoad_magma else load_magma
  
  for (ct in t_cell_types) {
    result <- run_conditional_regression(
      magma_df, 
      cell_type_markers[[ct]], 
      ct,
      microglia_all_genes,
      "Microglia"
    )
    if (!is.null(result)) {
      result$GWAS <- gwas_name
      conditional_results <- rbind(conditional_results, result)
    }
  }
}

if (nrow(conditional_results) > 0) {
  conditional_results$Bonferroni_Sig <- conditional_results$P_Target < bonferroni_threshold
  
  cat("\n  Conditional Regression Results:\n")
  cat("  Formula: Z_MAGMA = β₀ + β₁×T_Cell + β₂×Microglia + ε\n")
  cat("  ─────────────────────────────────────────────────────────────────\n")
  
  for (gwas in c("EOAD", "LOAD")) {
    cat("\n  ", gwas, ":\n", sep = "")
    gwas_results <- conditional_results[conditional_results$GWAS == gwas, ]
    
    for (i in 1:nrow(gwas_results)) {
      row <- gwas_results[i, ]
      sig_mark <- ifelse(row$P_Target < bonferroni_threshold, "**", 
                         ifelse(row$P_Target < 0.05, "*", ""))
      independent <- ifelse(row$P_Target < 0.05, " [Independent]", "")
      cat(sprintf("    %-20s: β_T=%6.3f (P=%.2e)%s | β_Micro=%6.3f (P=%.2e)%s\n",
                  row$CellType, row$Beta_Target, row$P_Target, sig_mark,
                  row$Beta_Covariate, row$P_Covariate, independent))
    }
  }
  
  fwrite(conditional_results, "results/tables/Table_MAGMA_Conditional_TCell_Microglia.csv")
}


# ==============================================================================
# Part 8: WGCNA Module MAGMA Enrichment (Integration)
# ==============================================================================
# Methods Reference:
#   Purpose: Test enrichment of GWAS signals within age-associated WGCNA modules
#   Conversion: Mouse genes to human orthologs via Ensembl BioMart
#   Test: One-sided Fisher's exact test
#   Multiple testing: Benjamini-Hochberg FDR
# ==============================================================================

cat("\n========================================\n")
cat("Part 8: WGCNA Module MAGMA Enrichment\n")
cat("========================================\n\n")

if (exists("wgcna_results") && !is.null(wgcna_results)) {
  
  ## Section 8.1: Extract Module Gene Lists
  module_colors_unique <- unique(wgcna_results$moduleColors)
  module_colors_unique <- module_colors_unique[module_colors_unique != "grey"]
  
  cat("  Testing ", length(module_colors_unique), " WGCNA modules for MAGMA enrichment...\n")
  
  ## Section 8.2: Convert Module Genes to Human
  module_genes_list <- lapply(module_colors_unique, function(color) {
    mouse_genes <- colnames(wgcna_results$datExpr)[wgcna_results$moduleColors == color]
    # Convert to human using stored conversion table
    human_genes <- wgcna_results$gene_conversion$Human_Gene[
      wgcna_results$gene_conversion$Mouse_Gene %in% mouse_genes]
    unique(human_genes[!is.na(human_genes) & human_genes != ""])
  })
  names(module_genes_list) <- module_colors_unique
  
  ## Section 8.3: Test Module Enrichment
  wgcna_enrichment_results <- data.frame()
  
  for (module_name in names(module_genes_list)) {
    module_genes <- module_genes_list[[module_name]]
    if (length(module_genes) < 10) next
    
    # Test EOAD
    if (!is.null(eoad_magma)) {
      eoad_result <- test_module_enrichment_fisher(module_genes, eoad_magma)
      if (!is.null(eoad_result)) {
        # Get age correlation
        me_name <- paste0("ME", module_name)
        age_cor <- NA
        age_p <- NA
        if (me_name %in% rownames(wgcna_results$moduleTraitCor)) {
          age_cor <- wgcna_results$moduleTraitCor[me_name, "Age_Months"]
          age_p <- wgcna_results$moduleTraitPvalue[me_name, "Age_Months"]
        }
        
        wgcna_enrichment_results <- rbind(wgcna_enrichment_results, data.frame(
          Module = module_name,
          GWAS = "EOAD",
          N_Module = eoad_result$n_module,
          N_Overlap = eoad_result$n_overlap,
          OR = eoad_result$odds_ratio,
          P_value = eoad_result$p_value,
          Age_Correlation = age_cor,
          Age_P = age_p
        ))
      }
    }
    
    # Test LOAD
    if (!is.null(load_magma)) {
      load_result <- test_module_enrichment_fisher(module_genes, load_magma)
      if (!is.null(load_result)) {
        me_name <- paste0("ME", module_name)
        age_cor <- NA
        age_p <- NA
        if (me_name %in% rownames(wgcna_results$moduleTraitCor)) {
          age_cor <- wgcna_results$moduleTraitCor[me_name, "Age_Months"]
          age_p <- wgcna_results$moduleTraitPvalue[me_name, "Age_Months"]
        }
        
        wgcna_enrichment_results <- rbind(wgcna_enrichment_results, data.frame(
          Module = module_name,
          GWAS = "LOAD",
          N_Module = load_result$n_module,
          N_Overlap = load_result$n_overlap,
          OR = load_result$odds_ratio,
          P_value = load_result$p_value,
          Age_Correlation = age_cor,
          Age_P = age_p
        ))
      }
    }
  }
  
  if (nrow(wgcna_enrichment_results) > 0) {
    wgcna_enrichment_results$FDR <- p.adjust(wgcna_enrichment_results$P_value, method = "BH")
    wgcna_enrichment_results$Age_Associated <- abs(wgcna_enrichment_results$Age_Correlation) > 0.2 & 
      wgcna_enrichment_results$Age_P < 0.1
    
    cat("\n  Significant WGCNA modules (FDR < 0.05):\n")
    sig_modules <- wgcna_enrichment_results[wgcna_enrichment_results$FDR < 0.05, ]
    if (nrow(sig_modules) > 0) {
      for (i in 1:nrow(sig_modules)) {
        row <- sig_modules[i, ]
        age_mark <- ifelse(row$Age_Associated, " [Age-associated]", "")
        cat(sprintf("    %s - %s: OR=%.2f, FDR=%.4f%s\n",
                    row$GWAS, row$Module, row$OR, row$FDR, age_mark))
      }
    } else {
      cat("    No modules reached FDR < 0.05\n")
    }
    
    fwrite(wgcna_enrichment_results, "results/tables/Table_WGCNA_Module_MAGMA_Enrichment.csv")
  }
  
} else {
  cat("  WGCNA results not available. Skipping module enrichment analysis.\n")
}


# ==============================================================================
# Part 9: Visualization
# ==============================================================================

cat("\n========================================\n")
cat("Part 9: Visualization\n")
cat("========================================\n\n")

## Section 9.1: Cell Type Enrichment Heatmap
create_celltype_heatmap <- function(results_df) {
  if (nrow(results_df) == 0) return(NULL)
  
  cat("  Creating cell type enrichment heatmap...\n")
  
  # Prepare data for heatmap
  heatmap_data <- results_df %>%
    select(CellType, GWAS, Beta) %>%
    pivot_wider(names_from = GWAS, values_from = Beta)
  
  p_data <- results_df %>%
    select(CellType, GWAS, P_Value) %>%
    pivot_wider(names_from = GWAS, values_from = P_Value)
  
  available_gwas <- intersect(c("EOAD", "LOAD"), colnames(heatmap_data))
  if (length(available_gwas) < 2) return(NULL)
  
  beta_matrix <- as.matrix(heatmap_data[, available_gwas, drop = FALSE])
  rownames(beta_matrix) <- heatmap_data$CellType
  
  p_matrix <- as.matrix(p_data[, available_gwas, drop = FALSE])
  rownames(p_matrix) <- p_data$CellType
  
  beta_matrix[is.na(beta_matrix)] <- 0
  
  # Significance markers (Methods: Bonferroni P < 0.005)
  sig_marks <- matrix("", nrow = nrow(p_matrix), ncol = ncol(p_matrix))
  sig_marks[!is.na(p_matrix) & p_matrix < 0.001] <- "***"
  sig_marks[!is.na(p_matrix) & p_matrix >= 0.001 & p_matrix < 0.005] <- "**"
  sig_marks[!is.na(p_matrix) & p_matrix >= 0.005 & p_matrix < 0.05] <- "*"
  
  # Row annotation by cell type category
  row_annotation <- data.frame(
    Category = ifelse(grepl("T_Cell|TCR", rownames(beta_matrix)), "T Cell",
                      ifelse(grepl("NK", rownames(beta_matrix)), "NK Cell",
                             ifelse(grepl("Microglia", rownames(beta_matrix)), "Microglia",
                                    ifelse(grepl("Abeta|APP", rownames(beta_matrix)), "AD Pathway", "Other"))))
  )
  rownames(row_annotation) <- rownames(beta_matrix)
  
  ann_colors <- list(
    Category = c("T Cell" = "#E41A1C", "NK Cell" = "#FF7F00", 
                 "Microglia" = "#377EB8", "AD Pathway" = "#4DAF4A", "Other" = "gray70")
  )
  
  pdf("results/figures/Figure_CellType_MAGMA_Heatmap.pdf", width = 8, height = 8)
  
  pheatmap(beta_matrix,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
           breaks = seq(-max(abs(beta_matrix), na.rm = TRUE), 
                        max(abs(beta_matrix), na.rm = TRUE), 
                        length.out = 51),
           display_numbers = sig_marks,
           number_color = "black",
           fontsize_number = 12,
           annotation_row = row_annotation,
           annotation_colors = ann_colors,
           main = "MAGMA Gene-Property Analysis\n(β coefficient, *P<0.05, **P<0.005, ***P<0.001)",
           fontsize = 10,
           fontsize_row = 9)
  
  dev.off()
  cat("    ✓ Cell type heatmap saved\n")
}

if (nrow(gene_property_results) > 0) {
  create_celltype_heatmap(gene_property_results)
}

## Section 9.2: Graham Module Enrichment Bar Plot
create_graham_barplot <- function(results_df) {
  if (nrow(results_df) == 0) return(NULL)
  
  cat("  Creating Graham module enrichment plot...\n")
  
  results_df$Significance <- ifelse(results_df$FDR < 0.05, "FDR < 0.05",
                                    ifelse(results_df$P_value < 0.05, "P < 0.05", "NS"))
  results_df$log2OR <- log2(results_df$OR)
  results_df$log2OR[!is.finite(results_df$log2OR)] <- 0
  
  p <- ggplot(results_df, aes(x = reorder(Module, log2OR), y = log2OR, fill = GWAS)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text(aes(label = ifelse(FDR < 0.05, "**", ifelse(P_value < 0.05, "*", "")),
                  y = log2OR + sign(log2OR) * 0.1),
              position = position_dodge(width = 0.8), size = 4) +
    scale_fill_manual(values = c("EOAD" = "#E41A1C", "LOAD" = "#377EB8")) +
    labs(title = "Graham et al. 2025 Module Enrichment",
         subtitle = "Fisher's exact test, *P<0.05, **FDR<0.05",
         x = "", y = "log2(Odds Ratio)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top") +
    coord_flip()
  
  ggsave("results/figures/Figure_Graham_Module_Enrichment.pdf", p, width = 10, height = 6)
  cat("    ✓ Graham module enrichment plot saved\n")
}

if (nrow(graham_enrichment_results) > 0) {
  create_graham_barplot(graham_enrichment_results)
}

## Section 9.3: GSEA Comparison Plot
create_gsea_comparison_plot <- function(eoad_gsea, load_gsea) {
  if (is.null(eoad_gsea$GO_BP) && is.null(load_gsea$GO_BP)) return(NULL)
  
  cat("  Creating GSEA comparison plot...\n")
  
  # Extract top pathways from each
  eoad_top <- if (!is.null(eoad_gsea$GO_BP)) {
    head(eoad_gsea$GO_BP@result[order(eoad_gsea$GO_BP@result$pvalue), ], 10)
  } else NULL
  
  load_top <- if (!is.null(load_gsea$GO_BP)) {
    head(load_gsea$GO_BP@result[order(load_gsea$GO_BP@result$pvalue), ], 10)
  } else NULL
  
  if (!is.null(eoad_top)) eoad_top$Source <- "EOAD"
  if (!is.null(load_top)) load_top$Source <- "LOAD"
  
  combined <- rbind(eoad_top, load_top)
  if (is.null(combined) || nrow(combined) == 0) return(NULL)
  
  # Truncate long pathway names
  combined$Description_Short <- substr(combined$Description, 1, 50)
  combined$Description_Short <- ifelse(nchar(combined$Description) > 50,
                                       paste0(combined$Description_Short, "..."),
                                       combined$Description_Short)
  
  p <- ggplot(combined, aes(x = NES, y = reorder(Description_Short, NES), fill = Source)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("EOAD" = "#E41A1C", "LOAD" = "#377EB8")) +
    labs(title = "GSEA Comparison: EOAD vs LOAD",
         subtitle = "GO Biological Process - Top 10 Pathways",
         x = "Normalized Enrichment Score (NES)", y = "") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top")
  
  ggsave("results/figures/Figure_GSEA_EOAD_LOAD_Comparison.pdf", p, width = 12, height = 8)
  cat("    ✓ GSEA comparison plot saved\n")
}

create_gsea_comparison_plot(eoad_gsea, load_gsea)

## Section 9.4: Conditional Regression Forest Plot
create_conditional_forest_plot <- function(cond_results) {
  if (nrow(cond_results) == 0) return(NULL)
  
  cat("  Creating conditional regression forest plot...\n")
  
  cond_results$Lower <- cond_results$Beta_Target - 1.96 * cond_results$SE_Target
  cond_results$Upper <- cond_results$Beta_Target + 1.96 * cond_results$SE_Target
  cond_results$Sig <- ifelse(cond_results$P_Target < 0.005, "Bonferroni", 
                             ifelse(cond_results$P_Target < 0.05, "Nominal", "NS"))
  
  p <- ggplot(cond_results, aes(x = Beta_Target, y = CellType, color = GWAS)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2,
                   position = position_dodge(0.4)) +
    geom_point(aes(shape = Sig), size = 3, position = position_dodge(0.4)) +
    scale_color_manual(values = c("EOAD" = "#E41A1C", "LOAD" = "#377EB8")) +
    scale_shape_manual(values = c("Bonferroni" = 16, "Nominal" = 17, "NS" = 1)) +
    labs(title = "T Cell Enrichment Conditional on Microglia",
         subtitle = "Z_MAGMA = β₀ + β₁×T_Cell + β₂×Microglia + ε",
         x = "β coefficient (95% CI)", y = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
  
  ggsave("results/figures/Figure_Conditional_Regression_Forest.pdf", p, width = 9, height = 5)
  cat("    ✓ Conditional regression forest plot saved\n")
}

if (nrow(conditional_results) > 0) {
  create_conditional_forest_plot(conditional_results)
}

## Section 9.5: Specific Genes Venn/Bar Plot
create_specific_genes_plot <- function(specific_genes) {
  if (is.null(specific_genes)) return(NULL)
  
  cat("  Creating specific genes bar plot...\n")
  
  gene_counts <- data.frame(
    Category = c("EOAD-specific", "Shared", "LOAD-specific"),
    Count = c(length(specific_genes$eoad_specific),
              length(specific_genes$shared),
              length(specific_genes$load_specific))
  )
  gene_counts$Category <- factor(gene_counts$Category, 
                                 levels = c("EOAD-specific", "Shared", "LOAD-specific"))
  
  p <- ggplot(gene_counts, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = Count), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("EOAD-specific" = "#E41A1C", 
                                 "Shared" = "#7E6148",
                                 "LOAD-specific" = "#377EB8")) +
    labs(title = "MAGMA Significant Genes (P < 0.05)",
         subtitle = "EOAD-specific vs LOAD-specific",
         x = "", y = "Number of Genes") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("results/figures/Figure_Specific_Genes_Barplot.pdf", p, width = 6, height = 5)
  cat("    ✓ Specific genes bar plot saved\n")
}

create_specific_genes_plot(specific_genes)

## Section 9.6: WGCNA Module-Trait Heatmap
if (exists("wgcna_results") && !is.null(wgcna_results)) {
  cat("  Creating WGCNA module-trait heatmap...\n")
  
  cor_matrix <- wgcna_results$moduleTraitCor
  pval_matrix <- wgcna_results$moduleTraitPvalue
  
  # Create significance markers
  sig_text <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  sig_text[pval_matrix < 0.001] <- "***"
  sig_text[pval_matrix >= 0.001 & pval_matrix < 0.01] <- "**"
  sig_text[pval_matrix >= 0.01 & pval_matrix < 0.05] <- "*"
  
  display_matrix <- matrix(paste0(round(cor_matrix, 2), sig_text),
                           nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  
  pdf("results/figures/Figure_WGCNA_ModuleTrait_Heatmap.pdf", width = 10, height = 12)
  
  col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
  
  ht <- Heatmap(cor_matrix,
                name = "Correlation",
                col = col_fun,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(display_matrix[i, j], x, y, gp = gpar(fontsize = 8))
                },
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 10),
                column_names_gp = gpar(fontsize = 10),
                column_title = "WGCNA Module-Trait Correlation (GSE137313)",
                heatmap_legend_param = list(title = "r"))
  
  draw(ht)
  dev.off()
  cat("    ✓ WGCNA module-trait heatmap saved\n")
}


# ==============================================================================
# Part 10: Results Summary and Export
# ==============================================================================

cat("\n========================================\n")
cat("Part 10: Results Summary\n")
cat("========================================\n\n")

cat("================================================================================\n")
cat("EOAD vs LOAD Data-Driven Pathway Discovery Analysis - Summary\n")
cat("================================================================================\n")

## Section 10.1: MAGMA Gene-Based Results
cat("\n【1. MAGMA Gene-Based Association】\n")
cat("─────────────────────────────────────────────────────────────────────\n")
cat("  Methods: MAGMA v1.10, gene boundaries 35kb/10kb, 1000G EUR LD\n")
cat("  Significance: Bonferroni P < 0.05/19,000 = 2.63e-6\n\n")

if (!is.null(eoad_magma)) {
  cat("  EOAD: ", nrow(eoad_magma), " genes tested\n")
  cat("    Bonferroni significant: ", sum(eoad_magma$Bonferroni_Sig, na.rm = TRUE), "\n")
  cat("    Suggestive (P < 1e-3): ", sum(eoad_magma$Suggestive_Sig, na.rm = TRUE), "\n")
}
if (!is.null(load_magma)) {
  cat("  LOAD: ", nrow(load_magma), " genes tested\n")
  cat("    Bonferroni significant: ", sum(load_magma$Bonferroni_Sig, na.rm = TRUE), "\n")
  cat("    Suggestive (P < 1e-3): ", sum(load_magma$Suggestive_Sig, na.rm = TRUE), "\n")
}

## Section 10.2: WGCNA Results
cat("\n【2. WGCNA Co-expression Network (GSE137313)】\n")
cat("─────────────────────────────────────────────────────────────────────\n")
cat("  Methods: Signed network, bicor, minModuleSize=30, mergeCutHeight=0.15\n")

if (exists("wgcna_results") && !is.null(wgcna_results)) {
  n_modules <- length(unique(wgcna_results$moduleColors)) - 1
  cat("  Modules detected: ", n_modules, "\n")
  cat("  Soft threshold power: ", wgcna_results$soft_power, 
      " (R² = ", round(wgcna_results$r2_achieved, 3), ")\n")
  cat("  Age-positive modules: ", length(wgcna_results$age_positive_modules), "\n")
  cat("  Age-negative modules: ", length(wgcna_results$age_negative_modules), "\n")
} else {
  cat("  WGCNA analysis not performed (data files not available)\n")
}

## Section 10.3: Graham Module Enrichment
cat("\n【3. Graham et al. 2025 Module Enrichment】\n")
cat("─────────────────────────────────────────────────────────────────────\n")
cat("  Methods: One-sided Fisher's exact test, BH FDR correction\n\n")

if (nrow(graham_enrichment_results) > 0) {
  sig_modules <- graham_enrichment_results[graham_enrichment_results$FDR < 0.05, ]
  if (nrow(sig_modules) > 0) {
    for (i in 1:nrow(sig_modules)) {
      cat("  ", sig_modules$GWAS[i], " - ", sig_modules$Module[i], 
          ": OR=", round(sig_modules$OR[i], 2), ", FDR=", 
          format(sig_modules$FDR[i], digits = 3), "\n", sep = "")
    }
  } else {
    cat("  No modules reached FDR < 0.05 significance\n")
  }
}

## Section 10.4: GSEA Results
cat("\n【4. GSEA Results】\n")
cat("─────────────────────────────────────────────────────────────────────\n")
cat("  Methods: clusterProfiler, gene set size 10-500, BH FDR\n\n")

if (!is.null(eoad_gsea$GO_BP)) {
  n_sig <- sum(eoad_gsea$GO_BP@result$p.adjust < 0.05)
  cat("  EOAD GO BP: ", nrow(eoad_gsea$GO_BP@result), " pathways, ", 
      n_sig, " FDR < 0.05\n")
}
if (!is.null(load_gsea$GO_BP)) {
  n_sig <- sum(load_gsea$GO_BP@result$p.adjust < 0.05)
  cat("  LOAD GO BP: ", nrow(load_gsea$GO_BP@result), " pathways, ", 
      n_sig, " FDR < 0.05\n")
}

## Section 10.5: Cell Type Enrichment
cat("\n【5. Cell Type MAGMA Gene-Property Analysis】\n")
cat("─────────────────────────────────────────────────────────────────────\n")
cat("  Methods: Linear regression Z_MAGMA ~ CellType\n")
cat("  Multiple testing: Bonferroni P < 0.005 (10 cell types)\n\n")

if (nrow(gene_property_results) > 0) {
  for (gwas in c("EOAD", "LOAD")) {
    cat("  ", gwas, ":\n", sep = "")
    gwas_results <- gene_property_results[gene_property_results$GWAS == gwas, ]
    sig_results <- gwas_results[gwas_results$P_Value < 0.05, ]
    
    if (nrow(sig_results) > 0) {
      sig_results <- sig_results[order(sig_results$P_Value), ]
      for (i in 1:min(5, nrow(sig_results))) {
        bonf_mark <- ifelse(sig_results$P_Value[i] < 0.005, " (Bonferroni)", "")
        cat("    ", sig_results$CellType[i], ": β=", round(sig_results$Beta[i], 3),
            ", P=", format(sig_results$P_Value[i], digits = 3), bonf_mark, "\n", sep = "")
      }
    } else {
      cat("    No significant cell types (P < 0.05)\n")
    }
  }
}

## Section 10.6: Conditional Analysis
cat("\n【6. Conditional Regression (T Cell | Microglia)】\n")
cat("─────────────────────────────────────────────────────────────────────\n")
cat("  Methods: Z_MAGMA = β₀ + β₁×T_Cell + β₂×Microglia + ε\n\n")

if (nrow(conditional_results) > 0) {
  for (gwas in c("EOAD", "LOAD")) {
    cat("  ", gwas, ":\n", sep = "")
    gwas_results <- conditional_results[conditional_results$GWAS == gwas, ]
    sig_results <- gwas_results[gwas_results$P_Target < 0.05, ]
    
    if (nrow(sig_results) > 0) {
      for (i in 1:nrow(sig_results)) {
        cat("    ", sig_results$CellType[i], ": β=", round(sig_results$Beta_Target[i], 3),
            ", P=", format(sig_results$P_Target[i], digits = 3), 
            " (independent of microglia)\n", sep = "")
      }
    } else {
      cat("    No T cell types remain significant after controlling for microglia\n")
    }
  }
}

## Section 10.7: EOAD-specific vs LOAD-specific
cat("\n【7. EOAD-Specific vs LOAD-Specific Genes】\n")
cat("─────────────────────────────────────────────────────────────────────\n")
cat("  Methods: ORA with hypergeometric test, BH FDR\n\n")

if (!is.null(specific_genes)) {
  cat("  EOAD-specific genes: ", length(specific_genes$eoad_specific), "\n")
  cat("  LOAD-specific genes: ", length(specific_genes$load_specific), "\n")
  cat("  Shared genes: ", length(specific_genes$shared), "\n")
}

## Section 10.8: Output Files
cat("\n【8. Output Files】\n")
cat("─────────────────────────────────────────────────────────────────────\n")
cat("  Tables:\n")
cat("    - Table_WGCNA_Module_Assignment.csv\n")
cat("    - Table_Graham_Module_Enrichment.csv\n")
cat("    - Table_MAGMA_GeneProperty_CellType.csv\n")
cat("    - Table_MAGMA_Conditional_TCell_Microglia.csv\n")
cat("    - Table_EOAD_GSEA_GO_BP.csv / Table_LOAD_GSEA_GO_BP.csv\n")
cat("    - Table_EOAD_Specific_ORA_GO_BP.csv / Table_LOAD_Specific_ORA_GO_BP.csv\n")
cat("    - Table_WGCNA_Module_MAGMA_Enrichment.csv\n")
cat("  Figures:\n")
cat("    - Figure_CellType_MAGMA_Heatmap.pdf\n")
cat("    - Figure_Graham_Module_Enrichment.pdf\n")
cat("    - Figure_GSEA_EOAD_LOAD_Comparison.pdf\n")
cat("    - Figure_Conditional_Regression_Forest.pdf\n")
cat("    - Figure_Specific_Genes_Barplot.pdf\n")
cat("    - Figure_WGCNA_ModuleTrait_Heatmap.pdf\n")

## Section 10.9: Save Complete Results Object
all_results <- list(
  # MAGMA results
  magma = list(
    eoad = if(exists("eoad_magma")) eoad_magma else NULL,
    load = if(exists("load_magma")) load_magma else NULL
  ),
  
  # WGCNA results
  wgcna = if(exists("wgcna_results")) wgcna_results else NULL,
  
  # Graham module enrichment
  graham = list(
    modules = graham_modules,
    enrichment = graham_enrichment_results
  ),
  
  # GSEA results
  gsea = list(
    eoad = eoad_gsea,
    load = load_gsea
  ),
  
  # Cell type analysis
  celltype = list(
    markers = cell_type_markers,
    gene_property = gene_property_results,
    conditional = conditional_results
  ),
  
  # Specific genes
  specific_genes = specific_genes,
  
  # ORA results
  ora = list(
    eoad_specific = if(exists("eoad_specific_ora")) eoad_specific_ora else NULL,
    load_specific = if(exists("load_specific_ora")) load_specific_ora else NULL,
    shared = if(exists("shared_ora")) shared_ora else NULL
  ),
  
  # WGCNA module enrichment
  wgcna_enrichment = if(exists("wgcna_enrichment_results")) wgcna_enrichment_results else NULL,
  
  # Metadata
  analysis_date = Sys.time(),
  r_version = R.version.string,
  methods_reference = list(
    magma = "MAGMA v1.10, gene boundaries 35kb/10kb, 1000G EUR LD, Bonferroni P<0.05/19000",
    wgcna = "GSE137313, signed network, bicor, minModuleSize=30, mergeCutHeight=0.15, R²>0.85",
    graham = "Graham et al. 2025, 5 modules, Fisher's exact test, BH FDR",
    gsea = "clusterProfiler, gene set 10-500, BH FDR",
    celltype = "MAGMA gene-property, Bonferroni P<0.005 for 10 cell types",
    conditional = "Z_MAGMA = β₀ + β₁×T_Cell + β₂×Microglia + ε"
  )
)

saveRDS(all_results, "results/Complete_Pathway_Analysis_Results.rds")

cat("\n================================================================================\n")
cat("Analysis Complete!\n")
cat("================================================================================\n\n")

# Print session info for reproducibility
cat("Session Info:\n")
cat("─────────────────────────────────────────────────────────────────────\n")
sessionInfo()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
