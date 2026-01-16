# ==============================================================================
# ADNI Independent Cohort Validation Analysis
# ==============================================================================
#
# Analysis Contents:
#   Part 1: Environment Setup
#   Part 2: PRS Calculation from LDpred2 Weights
#   Part 3: FreeSurfer UCSFFSX51 Neuroimaging Data Integration
#   Part 4: CSF Biomarker Integration (AlzBio3, Phases 1/GO/2)
#   Part 5: PRS-Biomarker Association Analysis
#   Part 6: Age-Stratified Analysis (<70 vs ≥70)
#   Part 7: MCI-to-Dementia Conversion (Cox Regression)
#   Part 8: Unsupervised K-means Clustering (k=3)
#   Part 9: Longitudinal Cognitive Trajectory (LMM)
#   Part 10: Visualization
#   Part 11: Results Summary
#
# FreeSurfer Variable Reference:
#   https://adni.bitbucket.io/reference/ucsffsx51.html
#
# ==============================================================================

# ==============================================================================
# Part 1: Environment Setup
# ==============================================================================

packages_cran <- c("data.table", "dplyr", "tidyr", "ggplot2", "cowplot",
                   "survival", "survminer", "lme4", "lmerTest",
                   "cluster", "factoextra", "pheatmap", "RColorBrewer",
                   "corrplot", "viridis")
packages_bioc <- c("bigsnpr")

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

dir.create("results/", showWarnings = FALSE)
dir.create("results/figures/", showWarnings = FALSE)
dir.create("results/tables/", showWarnings = FALSE)

# ==============================================================================
# Part 2: PRS Calculation from LDpred2 Weights
# ==============================================================================

## Section 2.1: Load Genotype Data
adni_rds_file <- "results/ADNI_genotypes.rds"
adni_bk_file <- "results/ADNI_genotypes.bk"

if (file.exists(adni_rds_file) && file.exists(adni_bk_file)) {
  adni_geno <- snp_attach(adni_rds_file)
} else {
  rds_path <- snp_readBed(
    bedfile = "WGS_Omni25_BIN_wo_ConsentsIssues.bed",
    backingfile = "results/ADNI_genotypes"
  )
  adni_geno <- snp_attach(rds_path)
}

## Section 2.2: Load PRS Weights
prs_weights <- readRDS("PRS_weights_all.rds")

## Section 2.3: SNP Matching (GRCh37, strand flip enabled)
adni_map <- adni_geno$map
colnames(adni_map) <- c("chr", "rsid", "genetic_dist", "pos", "a0", "a1")

snp_info_eoad <- prs_weights$snp_info_eoad
snp_info_load <- prs_weights$snp_info_load

match_eoad <- snp_match(
  sumstats = data.frame(
    chr = snp_info_eoad$chr, pos = snp_info_eoad$pos,
    a0 = snp_info_eoad$a0, a1 = snp_info_eoad$a1,
    beta = prs_weights$eoad_full
  ),
  info_snp = adni_map, strand_flip = TRUE,
  join_by_pos = TRUE, remove_dups = TRUE
)

match_load <- snp_match(
  sumstats = data.frame(
    chr = snp_info_load$chr, pos = snp_info_load$pos,
    a0 = snp_info_load$a0, a1 = snp_info_load$a1,
    beta = prs_weights$load_full
  ),
  info_snp = adni_map, strand_flip = TRUE,
  join_by_pos = TRUE, remove_dups = TRUE
)

## Section 2.4: Calculate PRS Function
G <- adni_geno$genotypes
G_imputed <- snp_fastImputeSimple(G, method = "mean0")

calculate_prs <- function(G, match_df, weights) {
  snp_idx <- match_df[["_NUM_ID_"]]
  weight_idx <- match_df[["_NUM_ID_.ss"]]
  matched_weights <- weights[weight_idx]
  matched_weights[is.na(matched_weights)] <- 0
  prs <- big_prodVec(G, matched_weights, ind.col = snp_idx, ncores = 1)
  return(as.vector(prs))
}


## Section 2.5: Calculate All PRS Types and Standardize to Z-scores
adni_prs <- data.frame(
  IID = adni_geno$fam$sample.ID,
  FID = adni_geno$fam$family.ID,
  PRS_EOAD_Global = scale(calculate_prs(G_imputed, match_eoad, prs_weights$eoad_full))[,1],
  PRS_EOAD_noAPOE = scale(calculate_prs(G_imputed, match_eoad, prs_weights$eoad_noAPOE))[,1],
  PRS_LOAD_Global = scale(calculate_prs(G_imputed, match_load, prs_weights$load_full))[,1],
  PRS_LOAD_noAPOE = scale(calculate_prs(G_imputed, match_load, prs_weights$load_noAPOE))[,1],
  PRS_EOAD_TCell = scale(calculate_prs(G_imputed, match_eoad, prs_weights$eoad_T_Cell))[,1],
  PRS_EOAD_Microglia = scale(calculate_prs(G_imputed, match_eoad, prs_weights$eoad_Microglia))[,1],
  PRS_EOAD_Abeta = scale(calculate_prs(G_imputed, match_eoad, prs_weights$eoad_Abeta))[,1],
  PRS_EOAD_APP = scale(calculate_prs(G_imputed, match_eoad, prs_weights$eoad_APP))[,1],
  PRS_EOAD_Oligo = scale(calculate_prs(G_imputed, match_eoad, prs_weights$eoad_Oligo))[,1],
  PRS_EOAD_Myelin = scale(calculate_prs(G_imputed, match_eoad, prs_weights$eoad_Myelin))[,1],
  PRS_LOAD_TCell = scale(calculate_prs(G_imputed, match_load, prs_weights$load_T_Cell))[,1],
  PRS_LOAD_Microglia = scale(calculate_prs(G_imputed, match_load, prs_weights$load_Microglia))[,1],
  PRS_LOAD_Abeta = scale(calculate_prs(G_imputed, match_load, prs_weights$load_Abeta))[,1],
  PRS_LOAD_APP = scale(calculate_prs(G_imputed, match_load, prs_weights$load_APP))[,1],
  PRS_LOAD_Oligo = scale(calculate_prs(G_imputed, match_load, prs_weights$load_Oligo))[,1],
  PRS_LOAD_Myelin = scale(calculate_prs(G_imputed, match_load, prs_weights$load_Myelin))[,1]
)

adni_prs$PTID <- adni_prs$IID
fwrite(adni_prs, "results/ADNI_PRS_All.csv")

# ==============================================================================
# Part 3: FreeSurfer UCSFFSX51 Neuroimaging Data Integration
# ==============================================================================

## Section 3.1: Load ADNIMERGE for Clinical Data
adnimerge <- fread("ADNIMERGE.csv")

## Section 3.2: Load FreeSurfer UCSFFSX51 Data
## Variable nomenclature: https://adni.bitbucket.io/reference/ucsffsx51.html
ucsffsx51_file <- "UCSFFSX51_08_01_22.csv"

if (file.exists(ucsffsx51_file)) {
  freesurfer <- fread(ucsffsx51_file)
  
  # Extract baseline FreeSurfer data
  fs_bl <- freesurfer[freesurfer$VISCODE2 == "bl", ]
  
  # FreeSurfer variable definitions from supplementary methods:
  # ST10CV: Intracranial Volume (ICV) - for normalization
  # ST128SV: Total WMH volume (WMHypoIntensities, mm³)
  # ST29SV: Left hippocampus subcortical volume
  # ST88SV: Right hippocampus subcortical volume
  # ST24CV: Left entorhinal cortical volume
  # ST83CV: Right entorhinal cortical volume
  # ST37SV: Left lateral ventricle
  # ST96SV: Right lateral ventricle
  # ST44CV: Left parahippocampal cortical volume
  # ST103CV: Right parahippocampal cortical volume
  
  fs_vars <- c("RID", "ST10CV", "ST128SV", "ST29SV", "ST88SV", "ST24CV", "ST83CV", 
               "ST37SV", "ST96SV", "ST44CV", "ST103CV")
  fs_vars_available <- fs_vars[fs_vars %in% colnames(fs_bl)]
  
  if (length(fs_vars_available) > 1) {
    fs_subset <- fs_bl[, ..fs_vars_available]
    fs_subset <- fs_subset[!duplicated(fs_subset$RID), ]
    
    # Convert to numeric
    for (col in fs_vars_available[-1]) {
      fs_subset[[col]] <- as.numeric(fs_subset[[col]])
    }
    
    # Calculate derived variables per methods
    # ICV for normalization
    if ("ST10CV" %in% colnames(fs_subset)) {
      fs_subset$ICV <- fs_subset$ST10CV
    }
    
    # Total WMH volume
    if ("ST128SV" %in% colnames(fs_subset)) {
      fs_subset$WMH_Total <- fs_subset$ST128SV
      if ("ICV" %in% colnames(fs_subset)) {
        fs_subset$WMH_Normalized <- fs_subset$ST128SV / fs_subset$ICV * 100
      }
    }
    
    # Bilateral hippocampal volume = ST29SV + ST88SV
    if (all(c("ST29SV", "ST88SV") %in% colnames(fs_subset))) {
      fs_subset$Hippocampus_Bilateral <- fs_subset$ST29SV + fs_subset$ST88SV
      if ("ICV" %in% colnames(fs_subset)) {
        fs_subset$Hippo_Normalized <- fs_subset$Hippocampus_Bilateral / fs_subset$ICV * 100
      }
    }
    
    # Bilateral entorhinal cortical volume = ST24CV + ST83CV
    if (all(c("ST24CV", "ST83CV") %in% colnames(fs_subset))) {
      fs_subset$Entorhinal_Bilateral <- fs_subset$ST24CV + fs_subset$ST83CV
      if ("ICV" %in% colnames(fs_subset)) {
        fs_subset$Entorhinal_Normalized <- fs_subset$Entorhinal_Bilateral / fs_subset$ICV * 100
      }
    }
    
    # Bilateral lateral ventricular volume = ST37SV + ST96SV
    if (all(c("ST37SV", "ST96SV") %in% colnames(fs_subset))) {
      fs_subset$Ventricle_Bilateral <- fs_subset$ST37SV + fs_subset$ST96SV
      if ("ICV" %in% colnames(fs_subset)) {
        fs_subset$Ventricle_Normalized <- fs_subset$Ventricle_Bilateral / fs_subset$ICV * 100
      }
    }
    
    # Bilateral parahippocampal = ST44CV + ST103CV
    if (all(c("ST44CV", "ST103CV") %in% colnames(fs_subset))) {
      fs_subset$Parahippo_Bilateral <- fs_subset$ST44CV + fs_subset$ST103CV
      if ("ICV" %in% colnames(fs_subset)) {
        fs_subset$Parahippo_Normalized <- fs_subset$Parahippo_Bilateral / fs_subset$ICV * 100
      }
    }
  }
}

## Section 3.3: Load Hippocampal Subfields
## ST131HS: Left CA1, ST137HS: Left subiculum
## ST139HS: Right CA1, ST145HS: Right subiculum
hippo_subfield_vars <- c("RID", "ST131HS", "ST137HS", "ST139HS", "ST145HS",
                         "ST132HS", "ST133HS", "ST140HS", "ST141HS")

if (exists("freesurfer")) {
  sf_vars_available <- hippo_subfield_vars[hippo_subfield_vars %in% colnames(fs_bl)]
  if (length(sf_vars_available) > 1) {
    sf_subset <- fs_bl[, ..sf_vars_available]
    sf_subset <- sf_subset[!duplicated(sf_subset$RID), ]
    for (col in sf_vars_available[-1]) {
      sf_subset[[col]] <- as.numeric(sf_subset[[col]])
    }
    
    # Bilateral CA1 = ST131HS + ST139HS
    if (all(c("ST131HS", "ST139HS") %in% colnames(sf_subset))) {
      sf_subset$CA1_Bilateral <- sf_subset$ST131HS + sf_subset$ST139HS
    }
    
    # Bilateral Subiculum = ST137HS + ST145HS
    if (all(c("ST137HS", "ST145HS") %in% colnames(sf_subset))) {
      sf_subset$Subiculum_Bilateral <- sf_subset$ST137HS + sf_subset$ST145HS
    }
    
    if (exists("fs_subset")) {
      fs_subset <- merge(fs_subset, sf_subset, by = "RID", all.x = TRUE)
    }
  }
}

## Section 3.4: Merge Clinical Data
clin_bl <- adnimerge[adnimerge$VISCODE == "bl", ]
clin_cols <- c("RID", "PTID", "AGE", "PTGENDER", "APOE4", "PTEDUCAT",
               "DX_bl", "DX", "MMSE", "CDRSB", "ADAS11", "ADAS13", "PHASE")
clin_cols_available <- clin_cols[clin_cols %in% colnames(clin_bl)]
clin_subset <- clin_bl[!duplicated(clin_bl$RID), ..clin_cols_available]

## Section 3.5: Merge All Data
merged_data <- as.data.frame(adni_prs)

# Merge clinical
if ("PTID" %in% colnames(clin_subset)) {
  merged_data <- merge(merged_data, clin_subset, by = "PTID", all.x = TRUE)
}

# Merge FreeSurfer
if (exists("fs_subset") && "RID" %in% colnames(merged_data)) {
  merged_data <- merge(merged_data, fs_subset, by = "RID", all.x = TRUE)
}

## Section 3.6: Create Derived Variables
# Age stratification (<70 vs ≥70) per methods
if ("AGE" %in% colnames(merged_data)) {
  merged_data$Age_Group <- factor(
    ifelse(merged_data$AGE < 70, "Younger (<70)", "Older (≥70)"),
    levels = c("Younger (<70)", "Older (≥70)")
  )
  merged_data$Age_Centered <- merged_data$AGE - 70
}

# APOE4 status and count
if ("APOE4" %in% colnames(merged_data)) {
  merged_data$APOE4_Status <- factor(
    ifelse(merged_data$APOE4 > 0, "APOE4+", "APOE4-"),
    levels = c("APOE4-", "APOE4+")
  )
  merged_data$APOE4_Count <- merged_data$APOE4
}

# Sex coding
if ("PTGENDER" %in% colnames(merged_data)) {
  merged_data$Sex_Numeric <- ifelse(merged_data$PTGENDER == "Male", 1, 0)
}

# ==============================================================================
# Part 4: CSF Biomarker Integration (AlzBio3, Phases 1/GO/2)
# ==============================================================================

## Section 4.1: Load CSF Biomarker Data
csf_file <- "UPENNBIOMK_MASTER.csv"

if (file.exists(csf_file)) {
  csf_data <- fread(csf_file)
  
  # Extract baseline CSF data
  csf_bl <- csf_data[csf_data$VISCODE2 == "bl", ]
  
  # CSF biomarker variables:
  # ABETA: Amyloid-beta 1-42 (pg/mL)
  # TAU: Total tau (pg/mL)
  # PTAU: Phosphorylated tau (pg/mL)
  
  csf_vars <- c("RID", "ABETA", "TAU", "PTAU")
  csf_vars_available <- csf_vars[csf_vars %in% colnames(csf_bl)]
  
  if (length(csf_vars_available) > 1) {
    csf_subset <- csf_bl[, ..csf_vars_available]
    csf_subset <- csf_subset[!duplicated(csf_subset$RID), ]
    
    # Convert to numeric
    for (col in csf_vars_available[-1]) {
      csf_subset[[col]] <- as.numeric(as.character(csf_subset[[col]]))
    }
    
    # Calculate derived ratios
    if (all(c("ABETA", "TAU") %in% colnames(csf_subset))) {
      csf_subset$TAU_ABETA_Ratio <- csf_subset$TAU / csf_subset$ABETA
    }
    if (all(c("ABETA", "PTAU") %in% colnames(csf_subset))) {
      csf_subset$PTAU_ABETA_Ratio <- csf_subset$PTAU / csf_subset$ABETA
    }
    
    # Merge with main data
    if ("RID" %in% colnames(merged_data)) {
      merged_data <- merge(merged_data, csf_subset, by = "RID", all.x = TRUE)
    }
  }
}

## Section 4.2: Create Amyloid Positivity Status
if ("ABETA" %in% colnames(merged_data)) {
  # Standard cutoff: Abeta < 192 pg/mL = amyloid positive
  merged_data$Amyloid_Positive <- factor(
    ifelse(merged_data$ABETA < 192, "A+", "A-"),
    levels = c("A-", "A+")
  )
}

# ==============================================================================
# Part 5: PRS-Biomarker Association Analysis
# ==============================================================================

## Section 5.1: Define Analysis Variables
prs_vars <- c("PRS_EOAD_Global", "PRS_EOAD_noAPOE", "PRS_LOAD_Global", "PRS_LOAD_noAPOE",
              "PRS_EOAD_TCell", "PRS_EOAD_Microglia", "PRS_EOAD_Abeta", "PRS_EOAD_APP",
              "PRS_EOAD_Oligo", "PRS_EOAD_Myelin",
              "PRS_LOAD_TCell", "PRS_LOAD_Microglia", "PRS_LOAD_Abeta", "PRS_LOAD_APP",
              "PRS_LOAD_Oligo", "PRS_LOAD_Myelin")

biomarker_vars <- c("ABETA", "TAU", "PTAU", "TAU_ABETA_Ratio", "PTAU_ABETA_Ratio",
                    "Hippo_Normalized", "Entorhinal_Normalized", "WMH_Normalized",
                    "Ventricle_Normalized", "MMSE", "CDRSB", "ADAS13")

## Section 5.2: Run PRS-Biomarker Associations
prs_biomarker_results <- list()

for (prs_var in prs_vars) {
  if (!prs_var %in% colnames(merged_data)) next
  
  for (bio_var in biomarker_vars) {
    if (!bio_var %in% colnames(merged_data)) next
    
    # Prepare analysis data
    analysis_cols <- c(prs_var, bio_var, "AGE", "PTGENDER", "APOE4")
    analysis_cols <- analysis_cols[analysis_cols %in% colnames(merged_data)]
    analysis_data <- merged_data[, analysis_cols, with = FALSE]
    analysis_data <- na.omit(analysis_data)
    
    if (nrow(analysis_data) < 50) next
    
    tryCatch({
      # Model: Biomarker ~ PRS + Age + Sex + APOE4
      formula_str <- paste(bio_var, "~", prs_var, "+ AGE + PTGENDER")
      if ("APOE4" %in% colnames(analysis_data)) {
        formula_str <- paste(formula_str, "+ APOE4")
      }
      
      model <- lm(as.formula(formula_str), data = analysis_data)
      model_summary <- summary(model)
      
      if (prs_var %in% rownames(coef(model_summary))) {
        coef_row <- coef(model_summary)[prs_var, ]
        
        prs_biomarker_results[[paste(prs_var, bio_var, sep = "_")]] <- data.frame(
          PRS = prs_var,
          Biomarker = bio_var,
          N = nrow(analysis_data),
          Beta = coef_row["Estimate"],
          SE = coef_row["Std. Error"],
          t_value = coef_row["t value"],
          P_value = coef_row["Pr(>|t|)"],
          R2 = model_summary$r.squared,
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) NULL)
  }
}

if (length(prs_biomarker_results) > 0) {
  prs_biomarker_df <- do.call(rbind, prs_biomarker_results)
  rownames(prs_biomarker_df) <- NULL
  prs_biomarker_df$P_FDR <- p.adjust(prs_biomarker_df$P_value, method = "fdr")
  fwrite(prs_biomarker_df, "results/tables/PRS_Biomarker_Associations.csv")
}

# ==============================================================================
# Part 6: Age-Stratified Analysis (<70 vs ≥70)
# ==============================================================================

## Section 6.1: Age-Stratified PRS-Biomarker Analysis
age_stratified_results <- list()

for (prs_var in c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin", "PRS_LOAD_Oligo", "PRS_LOAD_Myelin")) {
  if (!prs_var %in% colnames(merged_data)) next
  
  for (bio_var in c("WMH_Normalized", "Hippo_Normalized", "MMSE", "CDRSB")) {
    if (!bio_var %in% colnames(merged_data)) next
    
    for (age_grp in c("Younger (<70)", "Older (≥70)")) {
      
      # Subset by age group
      grp_data <- merged_data[merged_data$Age_Group == age_grp, ]
      
      analysis_cols <- c(prs_var, bio_var, "AGE", "PTGENDER", "APOE4")
      analysis_cols <- analysis_cols[analysis_cols %in% colnames(grp_data)]
      analysis_data <- grp_data[, analysis_cols, with = FALSE]
      analysis_data <- na.omit(analysis_data)
      
      if (nrow(analysis_data) < 30) next
      
      tryCatch({
        formula_str <- paste(bio_var, "~", prs_var, "+ AGE + PTGENDER")
        if ("APOE4" %in% colnames(analysis_data)) {
          formula_str <- paste(formula_str, "+ APOE4")
        }
        
        model <- lm(as.formula(formula_str), data = analysis_data)
        model_summary <- summary(model)
        
        if (prs_var %in% rownames(coef(model_summary))) {
          coef_row <- coef(model_summary)[prs_var, ]
          
          age_stratified_results[[paste(prs_var, bio_var, age_grp, sep = "_")]] <- data.frame(
            PRS = prs_var,
            Biomarker = bio_var,
            Age_Group = age_grp,
            N = nrow(analysis_data),
            Beta = coef_row["Estimate"],
            SE = coef_row["Std. Error"],
            t_value = coef_row["t value"],
            P_value = coef_row["Pr(>|t|)"],
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) NULL)
    }
  }
}

if (length(age_stratified_results) > 0) {
  age_stratified_df <- do.call(rbind, age_stratified_results)
  rownames(age_stratified_df) <- NULL
  age_stratified_df$P_FDR <- p.adjust(age_stratified_df$P_value, method = "fdr")
  fwrite(age_stratified_df, "results/tables/Age_Stratified_PRS_Associations.csv")
}

## Section 6.2: Interaction Analysis (PRS x Age)
interaction_results <- list()

for (prs_var in c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin", "PRS_LOAD_Oligo")) {
  if (!prs_var %in% colnames(merged_data)) next
  
  for (bio_var in c("WMH_Normalized", "Hippo_Normalized", "MMSE")) {
    if (!bio_var %in% colnames(merged_data)) next
    
    analysis_cols <- c(prs_var, bio_var, "AGE", "Age_Centered", "PTGENDER", "APOE4")
    analysis_cols <- analysis_cols[analysis_cols %in% colnames(merged_data)]
    analysis_data <- merged_data[, analysis_cols, with = FALSE]
    analysis_data <- na.omit(analysis_data)
    
    if (nrow(analysis_data) < 100) next
    
    tryCatch({
      # Model with interaction: Biomarker ~ PRS * Age_Centered + Sex + APOE4
      formula_str <- paste(bio_var, "~", prs_var, "* Age_Centered + PTGENDER")
      if ("APOE4" %in% colnames(analysis_data)) {
        formula_str <- paste(formula_str, "+ APOE4")
      }
      
      model <- lm(as.formula(formula_str), data = analysis_data)
      model_summary <- summary(model)
      
      int_term <- paste0(prs_var, ":Age_Centered")
      if (int_term %in% rownames(coef(model_summary))) {
        coef_row <- coef(model_summary)[int_term, ]
        
        interaction_results[[paste(prs_var, bio_var, sep = "_")]] <- data.frame(
          PRS = prs_var,
          Biomarker = bio_var,
          N = nrow(analysis_data),
          Interaction_Beta = coef_row["Estimate"],
          Interaction_SE = coef_row["Std. Error"],
          Interaction_t = coef_row["t value"],
          Interaction_P = coef_row["Pr(>|t|)"],
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) NULL)
  }
}

if (length(interaction_results) > 0) {
  interaction_df <- do.call(rbind, interaction_results)
  rownames(interaction_df) <- NULL
  fwrite(interaction_df, "results/tables/PRS_Age_Interaction_Results.csv")
}

# ==============================================================================
# Part 7: MCI-to-Dementia Conversion (Cox Regression)
# ==============================================================================

## Section 7.1: Prepare Longitudinal Data for Survival Analysis
# Load longitudinal diagnosis data
adnimerge_long <- fread("ADNIMERGE.csv")

# Identify MCI subjects at baseline
mci_baseline <- adnimerge_long[adnimerge_long$VISCODE == "bl" & 
                                 adnimerge_long$DX_bl %in% c("LMCI", "EMCI", "MCI"), ]

# Track conversion to dementia
conversion_data <- list()

for (rid in unique(mci_baseline$RID)) {
  subj_data <- adnimerge_long[adnimerge_long$RID == rid, ]
  subj_data <- subj_data[order(subj_data$M), ]
  
  # Find first dementia diagnosis
  dementia_visits <- subj_data[subj_data$DX %in% c("Dementia", "AD"), ]
  
  if (nrow(dementia_visits) > 0) {
    # Converted to dementia
    first_dementia <- dementia_visits[1, ]
    conversion_data[[as.character(rid)]] <- data.frame(
      RID = rid,
      Event = 1,
      Time_Months = first_dementia$M,
      stringsAsFactors = FALSE
    )
  } else {
    # Censored (did not convert)
    last_visit <- subj_data[nrow(subj_data), ]
    conversion_data[[as.character(rid)]] <- data.frame(
      RID = rid,
      Event = 0,
      Time_Months = last_visit$M,
      stringsAsFactors = FALSE
    )
  }
}

if (length(conversion_data) > 0) {
  conversion_df <- do.call(rbind, conversion_data)
  rownames(conversion_df) <- NULL
  
  # Merge with PRS data
  if ("RID" %in% colnames(merged_data)) {
    survival_data <- merge(conversion_df, merged_data, by = "RID", all.x = TRUE)
    survival_data <- survival_data[!is.na(survival_data$PRS_EOAD_Oligo), ]
  }
}

## Section 7.2: Cox Proportional Hazards Models
cox_results <- list()

if (exists("survival_data") && nrow(survival_data) > 50) {
  
  # Ensure time is positive
  survival_data <- survival_data[survival_data$Time_Months > 0, ]
  
  for (prs_var in c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin", "PRS_LOAD_Oligo", 
                    "PRS_EOAD_Global", "PRS_LOAD_Global")) {
    if (!prs_var %in% colnames(survival_data)) next
    
    tryCatch({
      # Standardize PRS
      survival_data$PRS_z <- scale(survival_data[[prs_var]])[, 1]
      
      # Cox model: Time ~ PRS + Age + Sex + APOE4
      formula_str <- "Surv(Time_Months, Event) ~ PRS_z + AGE + PTGENDER"
      if ("APOE4" %in% colnames(survival_data)) {
        formula_str <- paste(formula_str, "+ APOE4")
      }
      
      cox_model <- coxph(as.formula(formula_str), data = survival_data)
      cox_summary <- summary(cox_model)
      
      if ("PRS_z" %in% rownames(cox_summary$coefficients)) {
        coef_row <- cox_summary$coefficients["PRS_z", ]
        
        cox_results[[prs_var]] <- data.frame(
          PRS = prs_var,
          N = cox_summary$n,
          N_Events = cox_summary$nevent,
          HR = exp(coef_row["coef"]),
          HR_Lower = exp(coef_row["coef"] - 1.96 * coef_row["se(coef)"]),
          HR_Upper = exp(coef_row["coef"] + 1.96 * coef_row["se(coef)"]),
          P_value = coef_row["Pr(>|z|)"],
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) NULL)
  }
  
  if (length(cox_results) > 0) {
    cox_df <- do.call(rbind, cox_results)
    rownames(cox_df) <- NULL
    fwrite(cox_df, "results/tables/Cox_MCI_Conversion_Results.csv")
  }
}

## Section 7.3: Age-Stratified Cox Analysis
cox_age_results <- list()

if (exists("survival_data") && "Age_Group" %in% colnames(survival_data)) {
  
  for (age_grp in c("Younger (<70)", "Older (≥70)")) {
    grp_data <- survival_data[survival_data$Age_Group == age_grp, ]
    grp_data <- grp_data[grp_data$Time_Months > 0, ]
    
    if (nrow(grp_data) < 30) next
    
    for (prs_var in c("PRS_EOAD_Oligo", "PRS_LOAD_Oligo")) {
      if (!prs_var %in% colnames(grp_data)) next
      
      tryCatch({
        grp_data$PRS_z <- scale(grp_data[[prs_var]])[, 1]
        
        formula_str <- "Surv(Time_Months, Event) ~ PRS_z + AGE + PTGENDER"
        if ("APOE4" %in% colnames(grp_data)) {
          formula_str <- paste(formula_str, "+ APOE4")
        }
        
        cox_model <- coxph(as.formula(formula_str), data = grp_data)
        cox_summary <- summary(cox_model)
        
        if ("PRS_z" %in% rownames(cox_summary$coefficients)) {
          coef_row <- cox_summary$coefficients["PRS_z", ]
          
          cox_age_results[[paste(prs_var, age_grp, sep = "_")]] <- data.frame(
            PRS = prs_var,
            Age_Group = age_grp,
            N = cox_summary$n,
            N_Events = cox_summary$nevent,
            HR = exp(coef_row["coef"]),
            HR_Lower = exp(coef_row["coef"] - 1.96 * coef_row["se(coef)"]),
            HR_Upper = exp(coef_row["coef"] + 1.96 * coef_row["se(coef)"]),
            P_value = coef_row["Pr(>|z|)"],
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) NULL)
    }
  }
  
  if (length(cox_age_results) > 0) {
    cox_age_df <- do.call(rbind, cox_age_results)
    rownames(cox_age_df) <- NULL
    fwrite(cox_age_df, "results/tables/Cox_Age_Stratified_Results.csv")
  }
}

# ==============================================================================
# Part 8: Unsupervised K-means Clustering (k=3)
# ==============================================================================

## Section 8.1: Prepare Clustering Data
cluster_vars <- c("PRS_EOAD_Oligo", "PRS_EOAD_Myelin", "PRS_EOAD_TCell", 
                  "PRS_EOAD_Microglia", "PRS_EOAD_Abeta")
cluster_vars <- cluster_vars[cluster_vars %in% colnames(merged_data)]

if (length(cluster_vars) >= 3) {
  
  cluster_data <- merged_data[, c("PTID", cluster_vars), with = FALSE]
  cluster_data <- na.omit(cluster_data)
  
  # Standardize for clustering
  cluster_matrix <- scale(cluster_data[, cluster_vars, with = FALSE])
  
  ## Section 8.2: Determine Optimal k
  # Silhouette analysis
  sil_scores <- sapply(2:6, function(k) {
    km <- kmeans(cluster_matrix, centers = k, nstart = 25, iter.max = 100)
    mean(silhouette(km$cluster, dist(cluster_matrix))[, 3])
  })
  
  optimal_k <- which.max(sil_scores) + 1
  
  ## Section 8.3: Final Clustering with k=3 (per methods)
  set.seed(42)
  km_result <- kmeans(cluster_matrix, centers = 3, nstart = 50, iter.max = 100)
  
  cluster_data$Cluster <- factor(km_result$cluster)
  
  # Characterize clusters by PRS profiles
  cluster_profiles <- aggregate(
    cluster_data[, cluster_vars, with = FALSE],
    by = list(Cluster = cluster_data$Cluster),
    FUN = mean
  )
  
  # Assign biological labels based on dominant PRS
  cluster_labels <- sapply(1:3, function(i) {
    profile <- as.numeric(cluster_profiles[i, -1])
    names(profile) <- cluster_vars
    dominant <- names(which.max(profile))
    
    if (grepl("Oligo|Myelin", dominant)) return("Oligo-Driven")
    if (grepl("TCell|Microglia", dominant)) return("Immune-Driven")
    if (grepl("Abeta", dominant)) return("Abeta-Driven")
    return(paste0("Cluster_", i))
  })
  
  cluster_data$Subtype <- factor(cluster_labels[as.numeric(cluster_data$Cluster)])
  
  # Merge back to main data
  merged_data <- merge(merged_data, 
                       cluster_data[, c("PTID", "Cluster", "Subtype")], 
                       by = "PTID", all.x = TRUE)
  
  ## Section 8.4: Save Clustering Results
  fwrite(cluster_profiles, "results/tables/Cluster_PRS_Profiles.csv")
  
  cluster_summary <- data.frame(
    Cluster = 1:3,
    Label = cluster_labels,
    N = as.vector(table(cluster_data$Cluster)),
    stringsAsFactors = FALSE
  )
  fwrite(cluster_summary, "results/tables/Cluster_Summary.csv")
}

# ==============================================================================
# Part 9: Longitudinal Cognitive Trajectory (LMM)
# ==============================================================================

## Section 9.1: Prepare Longitudinal Cognitive Data
cognitive_long <- adnimerge_long[, c("RID", "VISCODE", "M", "MMSE", "CDRSB", "ADAS13")]
cognitive_long <- cognitive_long[!is.na(cognitive_long$MMSE) | 
                                   !is.na(cognitive_long$CDRSB) | 
                                   !is.na(cognitive_long$ADAS13), ]

# Merge with baseline PRS and covariates
baseline_vars <- c("RID", "PRS_EOAD_Oligo", "PRS_EOAD_Myelin", "PRS_LOAD_Oligo",
                   "AGE", "PTGENDER", "APOE4", "Age_Group", "Subtype")
baseline_vars <- baseline_vars[baseline_vars %in% colnames(merged_data)]

if ("RID" %in% colnames(merged_data)) {
  baseline_data <- merged_data[, baseline_vars, with = FALSE]
  baseline_data <- baseline_data[!duplicated(baseline_data$RID), ]
  
  lmm_data <- merge(cognitive_long, baseline_data, by = "RID", all.x = TRUE)
  lmm_data <- lmm_data[!is.na(lmm_data$PRS_EOAD_Oligo), ]
  
  # Convert time to years
  lmm_data$Time_Years <- lmm_data$M / 12
}

## Section 9.2: Linear Mixed Models for Cognitive Decline
lmm_results <- list()

if (exists("lmm_data") && nrow(lmm_data) > 100) {
  
  for (cog_var in c("MMSE", "CDRSB", "ADAS13")) {
    if (!cog_var %in% colnames(lmm_data)) next
    
    cog_data <- lmm_data[!is.na(lmm_data[[cog_var]]), ]
    if (nrow(cog_data) < 100) next
    
    for (prs_var in c("PRS_EOAD_Oligo", "PRS_LOAD_Oligo")) {
      if (!prs_var %in% colnames(cog_data)) next
      
      tryCatch({
        # Standardize PRS
        cog_data$PRS_z <- scale(cog_data[[prs_var]])[, 1]
        
        # LMM: Cognition ~ Time * PRS + Age + Sex + APOE4 + (1 + Time | RID)
        formula_str <- paste(cog_var, "~ Time_Years * PRS_z + AGE + PTGENDER")
        if ("APOE4" %in% colnames(cog_data)) {
          formula_str <- paste(formula_str, "+ APOE4")
        }
        formula_str <- paste(formula_str, "+ (1 + Time_Years | RID)")
        
        lmm_model <- lmer(as.formula(formula_str), data = cog_data,
                          control = lmerControl(optimizer = "bobyqa"))
        lmm_summary <- summary(lmm_model)
        
        # Extract fixed effects
        fixed_effects <- coef(lmm_summary)
        
        # Main effect of PRS
        if ("PRS_z" %in% rownames(fixed_effects)) {
          prs_effect <- fixed_effects["PRS_z", ]
          
          lmm_results[[paste(cog_var, prs_var, "main", sep = "_")]] <- data.frame(
            Cognitive_Measure = cog_var,
            PRS = prs_var,
            Effect = "Main",
            Estimate = prs_effect["Estimate"],
            SE = prs_effect["Std. Error"],
            t_value = prs_effect["t value"],
            stringsAsFactors = FALSE
          )
        }
        
        # Interaction effect (PRS x Time)
        int_term <- "Time_Years:PRS_z"
        if (int_term %in% rownames(fixed_effects)) {
          int_effect <- fixed_effects[int_term, ]
          
          lmm_results[[paste(cog_var, prs_var, "interaction", sep = "_")]] <- data.frame(
            Cognitive_Measure = cog_var,
            PRS = prs_var,
            Effect = "PRS x Time",
            Estimate = int_effect["Estimate"],
            SE = int_effect["Std. Error"],
            t_value = int_effect["t value"],
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) NULL)
    }
  }
  
  if (length(lmm_results) > 0) {
    lmm_df <- do.call(rbind, lmm_results)
    rownames(lmm_df) <- NULL
    fwrite(lmm_df, "results/tables/LMM_Cognitive_Trajectory_Results.csv")
  }
}

## Section 9.3: Age-Stratified LMM Analysis
lmm_age_results <- list()

if (exists("lmm_data") && "Age_Group" %in% colnames(lmm_data)) {
  
  for (age_grp in c("Younger (<70)", "Older (≥70)")) {
    grp_data <- lmm_data[lmm_data$Age_Group == age_grp, ]
    
    for (cog_var in c("MMSE", "CDRSB")) {
      if (!cog_var %in% colnames(grp_data)) next
      
      cog_data <- grp_data[!is.na(grp_data[[cog_var]]), ]
      if (nrow(cog_data) < 50) next
      
      tryCatch({
        cog_data$PRS_z <- scale(cog_data$PRS_EOAD_Oligo)[, 1]
        
        formula_str <- paste(cog_var, "~ Time_Years * PRS_z + AGE + PTGENDER + (1 | RID)")
        
        lmm_model <- lmer(as.formula(formula_str), data = cog_data,
                          control = lmerControl(optimizer = "bobyqa"))
        lmm_summary <- summary(lmm_model)
        fixed_effects <- coef(lmm_summary)
        
        int_term <- "Time_Years:PRS_z"
        if (int_term %in% rownames(fixed_effects)) {
          int_effect <- fixed_effects[int_term, ]
          
          lmm_age_results[[paste(cog_var, age_grp, sep = "_")]] <- data.frame(
            Cognitive_Measure = cog_var,
            Age_Group = age_grp,
            PRS = "PRS_EOAD_Oligo",
            Interaction_Estimate = int_effect["Estimate"],
            Interaction_SE = int_effect["Std. Error"],
            Interaction_t = int_effect["t value"],
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) NULL)
    }
  }
  
  if (length(lmm_age_results) > 0) {
    lmm_age_df <- do.call(rbind, lmm_age_results)
    rownames(lmm_age_df) <- NULL
    fwrite(lmm_age_df, "results/tables/LMM_Age_Stratified_Results.csv")
  }
}

# ==============================================================================
# Part 10: Visualization
# ==============================================================================

## Section 10.1: PRS Distribution by Diagnosis
if ("DX_bl" %in% colnames(merged_data) && "PRS_EOAD_Oligo" %in% colnames(merged_data)) {
  
  p_prs_dx <- ggplot(merged_data[!is.na(merged_data$DX_bl), ], 
                     aes(x = DX_bl, y = PRS_EOAD_Oligo, fill = DX_bl)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "PRS_EOAD_Oligo Distribution by Baseline Diagnosis",
      x = "Diagnosis",
      y = "PRS_EOAD_Oligo (z-score)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  
  ggsave("results/figures/PRS_by_Diagnosis.pdf", p_prs_dx, width = 8, height = 6)
  ggsave("results/figures/PRS_by_Diagnosis.png", p_prs_dx, width = 8, height = 6, dpi = 300)
}

## Section 10.2: Age-Stratified PRS-WMH Association
if ("WMH_Normalized" %in% colnames(merged_data) && 
    "PRS_EOAD_Oligo" %in% colnames(merged_data) &&
    "Age_Group" %in% colnames(merged_data)) {
  
  plot_data <- merged_data[!is.na(merged_data$WMH_Normalized) & 
                             !is.na(merged_data$PRS_EOAD_Oligo) &
                             !is.na(merged_data$Age_Group), ]
  
  p_wmh_age <- ggplot(plot_data, aes(x = PRS_EOAD_Oligo, y = WMH_Normalized, 
                                      color = Age_Group)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
    scale_color_manual(values = c("Younger (<70)" = "#E41A1C", "Older (≥70)" = "#377EB8")) +
    facet_wrap(~Age_Group) +
    labs(
      title = "PRS_EOAD_Oligo vs WMH by Age Group",
      subtitle = "Testing oligodendrocyte hypothesis: Higher PRS -> Greater WMH burden",
      x = "PRS_EOAD_Oligo (z-score)",
      y = "Normalized WMH Volume"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  
  ggsave("results/figures/PRS_WMH_Age_Stratified.pdf", p_wmh_age, width = 10, height = 5)
  ggsave("results/figures/PRS_WMH_Age_Stratified.png", p_wmh_age, width = 10, height = 5, dpi = 300)
}

## Section 10.3: Cluster Visualization
if ("Subtype" %in% colnames(merged_data) && length(cluster_vars) >= 2) {
  
  # PCA for visualization
  pca_data <- merged_data[!is.na(merged_data$Subtype), ]
  pca_matrix <- scale(pca_data[, cluster_vars, with = FALSE])
  pca_result <- prcomp(pca_matrix, center = FALSE, scale. = FALSE)
  
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Subtype = pca_data$Subtype
  )
  
  p_cluster <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype)) +
    geom_point(alpha = 0.6, size = 2) +
    stat_ellipse(level = 0.95, linewidth = 1) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Biological Subtypes Based on PRS Profiles",
      subtitle = "K-means clustering (k=3) on pathway-specific PRS",
      x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")
    ) +
    theme_bw(base_size = 12)
  
  ggsave("results/figures/Biological_Subtypes_PCA.pdf", p_cluster, width = 8, height = 6)
  ggsave("results/figures/Biological_Subtypes_PCA.png", p_cluster, width = 8, height = 6, dpi = 300)
}

## Section 10.4: Cox Survival Curves
if (exists("survival_data") && "PRS_EOAD_Oligo" %in% colnames(survival_data)) {
  
  # Tertile-based groups
  survival_data$PRS_Tertile <- cut(survival_data$PRS_EOAD_Oligo,
                                    breaks = quantile(survival_data$PRS_EOAD_Oligo, 
                                                      probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                    labels = c("Low", "Medium", "High"),
                                    include.lowest = TRUE)
  
  surv_fit <- survfit(Surv(Time_Months, Event) ~ PRS_Tertile, data = survival_data)
  
  p_surv <- ggsurvplot(
    surv_fit,
    data = survival_data,
    pval = TRUE,
    risk.table = TRUE,
    palette = c("#2166AC", "#F7F7F7", "#B2182B"),
    title = "MCI-to-Dementia Conversion by PRS_EOAD_Oligo Tertile",
    xlab = "Time (months)",
    ylab = "Conversion-free Probability",
    legend.title = "PRS Tertile",
    legend.labs = c("Low", "Medium", "High")
  )
  
  pdf("results/figures/Survival_PRS_Tertile.pdf", width = 10, height = 8)
  print(p_surv)
  dev.off()
  
  png("results/figures/Survival_PRS_Tertile.png", width = 10, height = 8, units = "in", res = 300)
  print(p_surv)
  dev.off()
}

## Section 10.5: Heatmap of PRS-Biomarker Associations
if (exists("prs_biomarker_df") && nrow(prs_biomarker_df) > 10) {
  
  # Create matrix for heatmap
  heatmap_data <- prs_biomarker_df[, c("PRS", "Biomarker", "t_value")]
  heatmap_wide <- reshape(heatmap_data, idvar = "PRS", timevar = "Biomarker", 
                          direction = "wide")
  rownames(heatmap_wide) <- heatmap_wide$PRS
  heatmap_matrix <- as.matrix(heatmap_wide[, -1])
  colnames(heatmap_matrix) <- gsub("t_value.", "", colnames(heatmap_matrix))
  
  # Remove rows/cols with all NA
  heatmap_matrix <- heatmap_matrix[rowSums(!is.na(heatmap_matrix)) > 0, 
                                    colSums(!is.na(heatmap_matrix)) > 0]
  
  if (nrow(heatmap_matrix) > 2 && ncol(heatmap_matrix) > 2) {
    pdf("results/figures/PRS_Biomarker_Heatmap.pdf", width = 12, height = 10)
    pheatmap(
      heatmap_matrix,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      breaks = seq(-5, 5, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      main = "PRS-Biomarker Association t-statistics",
      fontsize = 10,
      na_col = "grey90"
    )
    dev.off()
  }
}

# ==============================================================================
# Part 11: ggseg Brain Region Visualization (Optional - requires ggseg package)
# ==============================================================================
# Note: This section provides brain surface visualization using ggseg atlas
# If ggseg is not installed, this section will be skipped
# ==============================================================================

cat("\n========================================\n")
cat("Part 11: ggseg Brain Region Visualization\n")
cat("========================================\n")

# Check if ggseg is available
ggseg_available <- requireNamespace("ggseg", quietly = TRUE)

if (ggseg_available) {
  library(ggseg)
  
  # ST-to-DK atlas mapping for cortical regions
  st_dk_mapping <- data.frame(
    ST_Code = c("ST24CV", "ST83CV", "ST44CV", "ST103CV", "ST40CV", "ST99CV",
                "ST31CV", "ST90CV", "ST32CV", "ST91CV", "ST26CV", "ST85CV"),
    dk_region = c("entorhinal", "entorhinal", "parahippocampal", "parahippocampal",
                  "middle temporal", "middle temporal", "inferior parietal", "inferior parietal",
                  "inferior temporal", "inferior temporal", "fusiform", "fusiform"),
    hemi = c("left", "right", "left", "right", "left", "right",
             "left", "right", "left", "right", "left", "right"),
    stringsAsFactors = FALSE
  )
  
  # Prepare PRS-brain region data for ggseg visualization
  if (exists("prs_biomarker_df") && "Hippo_Normalized" %in% prs_biomarker_df$Biomarker) {
    
    cat("  Preparing ggseg visualization data...\n")
    
    # This is a simplified example - full implementation would analyze all cortical regions
    # For now, we'll create a placeholder message
    cat("  Note: Full ggseg visualization requires cortical region analysis.\n")
    cat("  See ADNI_Brain_Region_Visualization_NC_FINAL.R for complete implementation.\n")
  }
  
} else {
  cat("  ggseg package not installed. Skipping brain surface visualization.\n")
  cat("  To install: install.packages('ggseg')\n")
}

# ==============================================================================
# Part 12: Results Summary
# ==============================================================================

cat("\n================================================================================\n")
cat("ADNI Independent Cohort Validation - Results Summary\n")
cat("================================================================================\n")

## Section 11.1: Sample Characteristics
cat("\n--- Sample Characteristics ---\n")
cat(sprintf("Total N: %d\n", nrow(merged_data)))

if ("DX_bl" %in% colnames(merged_data)) {
  cat("\nBaseline Diagnosis:\n")
  print(table(merged_data$DX_bl, useNA = "ifany"))
}

if ("Age_Group" %in% colnames(merged_data)) {
  cat("\nAge Groups:\n")
  print(table(merged_data$Age_Group, useNA = "ifany"))
}

if ("APOE4_Status" %in% colnames(merged_data)) {
  cat("\nAPOE4 Status:\n")
  print(table(merged_data$APOE4_Status, useNA = "ifany"))
}

if ("Subtype" %in% colnames(merged_data)) {
  cat("\nBiological Subtypes:\n")
  print(table(merged_data$Subtype, useNA = "ifany"))
}

## Section 11.2: Key Findings Summary
cat("\n--- Key Findings ---\n")

# PRS-Biomarker associations
if (exists("prs_biomarker_df")) {
  sig_assoc <- prs_biomarker_df[prs_biomarker_df$P_FDR < 0.05, ]
  cat(sprintf("\nSignificant PRS-Biomarker associations (FDR < 0.05): %d\n", nrow(sig_assoc)))
  
  if (nrow(sig_assoc) > 0) {
    cat("\nTop associations:\n")
    sig_assoc <- sig_assoc[order(sig_assoc$P_FDR), ]
    for (i in 1:min(10, nrow(sig_assoc))) {
      cat(sprintf("  %s ~ %s: Beta=%.3f, P_FDR=%.4f\n",
                  sig_assoc$Biomarker[i], sig_assoc$PRS[i],
                  sig_assoc$Beta[i], sig_assoc$P_FDR[i]))
    }
  }
}

# Age-stratified findings
if (exists("age_stratified_df")) {
  cat("\n--- Age-Stratified Analysis ---\n")
  
  for (age_grp in c("Younger (<70)", "Older (≥70)")) {
    grp_results <- age_stratified_df[age_stratified_df$Age_Group == age_grp, ]
    sig_grp <- grp_results[grp_results$P_value < 0.05, ]
    cat(sprintf("\n%s - Significant associations (P < 0.05): %d\n", age_grp, nrow(sig_grp)))
  }
}

# Cox regression findings
if (exists("cox_df")) {
  cat("\n--- MCI-to-Dementia Conversion (Cox Regression) ---\n")
  for (i in 1:nrow(cox_df)) {
    sig <- ifelse(cox_df$P_value[i] < 0.05, "*", "")
    cat(sprintf("  %s: HR=%.2f (%.2f-%.2f), P=%.4f %s\n",
                cox_df$PRS[i], cox_df$HR[i], cox_df$HR_Lower[i], cox_df$HR_Upper[i],
                cox_df$P_value[i], sig))
  }
}

# LMM findings
if (exists("lmm_df")) {
  cat("\n--- Longitudinal Cognitive Trajectory (LMM) ---\n")
  int_results <- lmm_df[lmm_df$Effect == "PRS x Time", ]
  for (i in 1:nrow(int_results)) {
    cat(sprintf("  %s ~ %s x Time: Estimate=%.4f, t=%.2f\n",
                int_results$Cognitive_Measure[i], int_results$PRS[i],
                int_results$Estimate[i], int_results$t_value[i]))
  }
}

## Section 11.3: Save Final Dataset
fwrite(merged_data, "results/ADNI_Complete_Validation_Data.csv")
cat("\n\nFinal dataset saved: results/ADNI_Complete_Validation_Data.csv\n")

## Section 11.4: Output Files Summary
cat("\n--- Output Files ---\n")
cat("Tables:\n")
table_files <- list.files("results/tables/", pattern = "\\.csv$", full.names = FALSE)
for (f in table_files) cat(sprintf("  - %s\n", f))

cat("\nFigures:\n")
figure_files <- list.files("results/figures/", pattern = "\\.(pdf|png)$", full.names = FALSE)
for (f in unique(gsub("\\.(pdf|png)$", "", figure_files))) {
  cat(sprintf("  - %s\n", f))
}

cat("\n================================================================================\n")
cat("Part 1-11 Complete!\n")
cat("================================================================================\n")

# ==============================================================================
# Part 13: Mechanistic Validation - WMH & Hippocampal Subfield Analysis
# ==============================================================================
# Purpose: Validate the "Oligodendrocyte Hypothesis" core mechanism
# - 13.1: PRS_Oligo vs WMH (White Matter Hyperintensities) - Direct myelin hypothesis validation
# - 13.2: Hippocampal and Subfield-Specific Atrophy - Echo gsMap results
# - 13.3: Subcortical Structure Analysis - Supplementary validation
# - 13.4: Subtype-Specific WMH Analysis
# - 13.5: Combined Figure 7 Generation
# - 13.6: Summary Report
# ==============================================================================

cat("\n================================================================================\n")
cat("Part 13: Mechanistic Validation - WMH & Hippocampal Subfield Analysis\n")
cat("================================================================================\n")

# Color configuration
col_young <- "#E41A1C"
col_old <- "#377EB8"

# Output directory for Part 13
output_dir_p13 <- "results/figures/Part13_Mechanistic"
dir.create(output_dir_p13, showWarnings = FALSE, recursive = TRUE)

# Determine ICV column
icv_col <- ifelse("ST10CV" %in% colnames(merged_data), "ST10CV", 
                  ifelse("ICV" %in% colnames(merged_data), "ICV", NULL))

# Check key variables
cat("\nKey variables check:\n")
cat(sprintf("  ST128SV (WMH): %s\n", ifelse("ST128SV" %in% colnames(merged_data), "Available", "Missing")))
cat(sprintf("  PRS_EOAD_Oligo: %s\n", ifelse("PRS_EOAD_Oligo" %in% colnames(merged_data), "Available", "Missing")))
cat(sprintf("  Age_Group: %s\n", ifelse("Age_Group" %in% colnames(merged_data), "Available", "Missing")))
cat(sprintf("  Subtype: %s\n", ifelse("Subtype" %in% colnames(merged_data), "Available", "Missing")))

# ==============================================================================
# Part 13.1: PRS_Oligo vs WMH (White Matter Hyperintensities) Association
# ==============================================================================
# Core hypothesis: Higher PRS_Oligo -> Greater myelin damage -> Larger WMH volume
# This is the clinical gold-standard evidence for validating the "oligodendrocyte hypothesis"
# ==============================================================================

cat("\n========================================\n")
cat("13.1 Analyzing PRS_Oligo vs WMH (ST128SV)\n")
cat("========================================\n")

p_wmh <- NULL
wmh_res_df <- NULL

has_wmh <- "ST128SV" %in% colnames(merged_data)
has_prs <- "PRS_EOAD_Oligo" %in% colnames(merged_data)
has_age_group <- "Age_Group" %in% colnames(merged_data)

if (has_wmh && has_prs && has_age_group && !is.null(icv_col)) {
  
  # Prepare data
  wmh_cols <- c("ST128SV", "PRS_EOAD_Oligo", "AGE", "PTGENDER", "Age_Group", icv_col)
  wmh_data <- merged_data[, wmh_cols, with = FALSE]
  wmh_data <- na.omit(wmh_data)
  
  cat(sprintf("  WMH analysis sample size: N=%d\n", nrow(wmh_data)))
  
  # Log transform WMH (because WMH is typically right-skewed)
  wmh_data$WMH_Log <- log(wmh_data$ST128SV + 1)
  
  # Standardize PRS and ICV
  wmh_data$PRS_z <- as.numeric(scale(wmh_data$PRS_EOAD_Oligo))
  wmh_data$ICV_z <- as.numeric(scale(wmh_data[[icv_col]]))
  
  # Ensure factor levels are correct
  wmh_data$Age_Group <- factor(wmh_data$Age_Group, 
                                levels = c("Younger (<70)", "Older (≥70)"))
  
  # Age-stratified analysis
  wmh_stats <- list()
  age_groups <- c("Younger (<70)", "Older (≥70)")
  
  cat("\n  Age-stratified regression results:\n")
  cat("  ", paste(rep("-", 60), collapse = ""), "\n")
  
  for (grp in age_groups) {
    sub_dat <- wmh_data[wmh_data$Age_Group == grp, ]
    
    if (nrow(sub_dat) < 30) {
      cat(sprintf("  %s: N=%d (skipped - insufficient sample)\n", grp, nrow(sub_dat)))
      next
    }
    
    tryCatch({
      mod <- lm(WMH_Log ~ PRS_z + AGE + PTGENDER + ICV_z, data = sub_dat)
      res <- summary(mod)
      
      if ("PRS_z" %in% rownames(coef(res))) {
        coef_row <- coef(res)["PRS_z", ]
        
        wmh_stats[[grp]] <- data.frame(
          Age_Group = grp,
          N = nrow(sub_dat),
          Beta = coef_row["Estimate"],
          SE = coef_row["Std. Error"],
          t_value = coef_row["t value"],
          P_value = coef_row["Pr(>|t|)"],
          R2 = res$r.squared,
          stringsAsFactors = FALSE
        )
        
        sig <- ifelse(coef_row["Pr(>|t|)"] < 0.001, "***",
                      ifelse(coef_row["Pr(>|t|)"] < 0.01, "**",
                             ifelse(coef_row["Pr(>|t|)"] < 0.05, "*", "")))
        cat(sprintf("  %s: N=%d, Beta=%.4f, SE=%.4f, t=%.3f, P=%.4f %s\n",
                    grp, nrow(sub_dat), coef_row["Estimate"], 
                    coef_row["Std. Error"], coef_row["t value"], 
                    coef_row["Pr(>|t|)"], sig))
      }
    }, error = function(e) {
      cat(sprintf("  %s: Error - %s\n", grp, e$message))
    })
  }
  
  # Full sample analysis
  cat("\n  Full sample analysis:\n")
  tryCatch({
    mod_all <- lm(WMH_Log ~ PRS_z + AGE + PTGENDER + ICV_z, data = wmh_data)
    res_all <- summary(mod_all)
    
    if ("PRS_z" %in% rownames(coef(res_all))) {
      coef_all <- coef(res_all)["PRS_z", ]
      
      wmh_stats[["All"]] <- data.frame(
        Age_Group = "All",
        N = nrow(wmh_data),
        Beta = coef_all["Estimate"],
        SE = coef_all["Std. Error"],
        t_value = coef_all["t value"],
        P_value = coef_all["Pr(>|t|)"],
        R2 = res_all$r.squared,
        stringsAsFactors = FALSE
      )
      
      sig <- ifelse(coef_all["Pr(>|t|)"] < 0.001, "***",
                    ifelse(coef_all["Pr(>|t|)"] < 0.01, "**",
                           ifelse(coef_all["Pr(>|t|)"] < 0.05, "*", "")))
      cat(sprintf("  All: N=%d, Beta=%.4f, SE=%.4f, t=%.3f, P=%.4f %s\n",
                  nrow(wmh_data), coef_all["Estimate"], 
                  coef_all["Std. Error"], coef_all["t value"], 
                  coef_all["Pr(>|t|)"], sig))
    }
  }, error = function(e) {
    cat(sprintf("  All: Error - %s\n", e$message))
  })
  
  if (length(wmh_stats) > 0) {
    wmh_res_df <- do.call(rbind, wmh_stats)
    rownames(wmh_res_df) <- NULL
    
    fwrite(wmh_res_df, file.path(output_dir_p13, "WMH_PRS_Association_Results.csv"))
    cat("\n  Saved: WMH_PRS_Association_Results.csv\n")
    
    # Visualization: Scatter plot + regression line
    p_wmh <- ggplot(wmh_data, aes(x = PRS_z, y = WMH_Log, color = Age_Group)) +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
      scale_color_manual(values = c("Younger (<70)" = col_young, 
                                    "Older (≥70)" = col_old)) +
      facet_wrap(~Age_Group, scales = "fixed") +
      labs(
        title = "A. Oligodendrocyte Genetic Risk and White Matter Damage",
        subtitle = "Hypothesis: Higher Oligo-PRS -> Greater WMH burden",
        x = "PRS_EOAD_Oligo (z-score)",
        y = "log(WMH Volume + 1)",
        color = "Age Group"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "grey95")
      )
    
    # Add statistical annotations
    wmh_res_stratified <- wmh_res_df[wmh_res_df$Age_Group != "All", ]
    if (nrow(wmh_res_stratified) > 0) {
      ann_data <- data.frame(
        Age_Group = factor(wmh_res_stratified$Age_Group, 
                           levels = c("Younger (<70)", "Older (≥70)")),
        PRS_z = min(wmh_data$PRS_z) + 0.1 * diff(range(wmh_data$PRS_z)),
        WMH_Log = max(wmh_data$WMH_Log) * 0.95,
        label = sprintf("Beta = %.3f\nP = %.4f\nN = %d", 
                        wmh_res_stratified$Beta,
                        wmh_res_stratified$P_value,
                        wmh_res_stratified$N),
        stringsAsFactors = FALSE
      )
      
      p_wmh <- p_wmh +
        geom_text(data = ann_data, aes(label = label), 
                  color = "black", hjust = 0, vjust = 1, size = 3.5)
    }
    
    ggsave(file.path(output_dir_p13, "Figure7A_WMH_PRS_Validation.pdf"), 
           p_wmh, width = 10, height = 5, dpi = 300)
    ggsave(file.path(output_dir_p13, "Figure7A_WMH_PRS_Validation.png"), 
           p_wmh, width = 10, height = 5, dpi = 300)
    cat("  Saved: Figure7A_WMH_PRS_Validation.pdf/png\n")
    
    # Interaction test: PRS x Age
    cat("\n  Testing PRS x Age interaction:\n")
    wmh_data$Age_Binary <- ifelse(wmh_data$Age_Group == "Younger (<70)", 0, 1)
    
    tryCatch({
      mod_int <- lm(WMH_Log ~ PRS_z * Age_Binary + AGE + PTGENDER + ICV_z, 
                    data = wmh_data)
      res_int <- summary(mod_int)
      
      if ("PRS_z:Age_Binary" %in% rownames(coef(res_int))) {
        int_coef <- coef(res_int)["PRS_z:Age_Binary", ]
        cat(sprintf("  Interaction (PRS x Age): Beta=%.4f, P=%.4f\n",
                    int_coef["Estimate"], int_coef["Pr(>|t|)"]))
        
        if (int_coef["Pr(>|t|)"] < 0.05) {
          cat("  *** Significant interaction: PRS effect differs by age group!\n")
        } else {
          cat("  Interaction not significant: PRS effect similar across age groups.\n")
        }
      }
    }, error = function(e) {
      cat(sprintf("  Interaction test error: %s\n", e$message))
    })
  }
  
} else {
  cat("  Required variables missing. Skipping WMH analysis.\n")
  cat(sprintf("    ST128SV: %s, PRS_EOAD_Oligo: %s, Age_Group: %s\n",
              has_wmh, has_prs, has_age_group))
}

# ==============================================================================
# Part 13.2: Hippocampal and Subfield Analysis
# ==============================================================================
# Purpose: Validate gsMap results - CA1/Subiculum should show specific atrophy
# Variables: ST29SV/ST88SV (Left/Right hippocampus), ST131HS-ST146HS (subfields if available)
# ==============================================================================

cat("\n========================================\n")
cat("13.2 Hippocampal and Subfield Analysis\n")
cat("========================================\n")

p_hipp <- NULL
hipp_res_df <- NULL

# Check hippocampus variables
hipp_vars <- c("ST29SV", "ST88SV")
hipp_available <- hipp_vars[hipp_vars %in% colnames(merged_data)]

# Check hippocampal subfield variables (FreeSurfer 7.x)
subfield_vars <- c("ST131HS", "ST132HS", "ST133HS", "ST134HS", "ST135HS", 
                   "ST136HS", "ST137HS", "ST138HS", "ST139HS", "ST140HS",
                   "ST141HS", "ST142HS", "ST143HS", "ST144HS", "ST145HS", "ST146HS")
subfield_available <- subfield_vars[subfield_vars %in% colnames(merged_data)]

cat(sprintf("  Hippocampus variables available: %s\n", 
            ifelse(length(hipp_available) > 0, paste(hipp_available, collapse = ", "), "None")))
cat(sprintf("  Subfield variables available: %d of %d\n", 
            length(subfield_available), length(subfield_vars)))

# Hippocampal subfield mapping table
subfield_mapping <- data.frame(
  ST_Code = c("ST131HS", "ST132HS", "ST133HS", "ST134HS", "ST135HS", 
              "ST136HS", "ST137HS", "ST138HS", "ST139HS", "ST140HS",
              "ST141HS", "ST142HS", "ST143HS", "ST144HS", "ST145HS", "ST146HS"),
  Region_Name = c("Left_CA1", "Left_CA2/3", "Left_CA4", "Left_DG", "Left_Presubiculum",
                  "Left_Subiculum", "Left_Fimbria", "Left_HATA", 
                  "Right_CA1", "Right_CA2/3", "Right_CA4", "Right_DG", "Right_Presubiculum",
                  "Right_Subiculum", "Right_Fimbria", "Right_HATA"),
  Hemi = c(rep("Left", 8), rep("Right", 8)),
  stringsAsFactors = FALSE
)

# Analyze hippocampal volume
if (length(hipp_available) > 0 && has_prs && has_age_group && !is.null(icv_col)) {
  
  cat("\n  Analyzing hippocampal volumes...\n")
  
  hipp_stats <- list()
  
  for (hipp_var in hipp_available) {
    
    hipp_name <- ifelse(hipp_var == "ST29SV", "Left Hippocampus", "Right Hippocampus")
    
    # Prepare data
    hipp_cols <- c(hipp_var, "PRS_EOAD_Oligo", "AGE", "PTGENDER", "Age_Group", icv_col)
    hipp_data <- merged_data[, hipp_cols, with = FALSE]
    hipp_data <- na.omit(hipp_data)
    
    if (nrow(hipp_data) < 50) {
      cat(sprintf("  %s: N=%d (skipped)\n", hipp_name, nrow(hipp_data)))
      next
    }
    
    # ICV normalization
    hipp_data$Vol_Norm <- hipp_data[[hipp_var]] / hipp_data[[icv_col]] * 1000
    hipp_data$PRS_z <- as.numeric(scale(hipp_data$PRS_EOAD_Oligo))
    hipp_data$Age_Group <- factor(hipp_data$Age_Group, 
                                   levels = c("Younger (<70)", "Older (≥70)"))
    
    # Age-stratified analysis
    for (grp in c("Younger (<70)", "Older (≥70)", "All")) {
      
      if (grp == "All") {
        sub_dat <- hipp_data
      } else {
        sub_dat <- hipp_data[hipp_data$Age_Group == grp, ]
      }
      
      if (nrow(sub_dat) < 30) next
      
      tryCatch({
        mod <- lm(Vol_Norm ~ PRS_z + AGE + PTGENDER, data = sub_dat)
        res <- summary(mod)
        
        if ("PRS_z" %in% rownames(coef(res))) {
          coef_row <- coef(res)["PRS_z", ]
          
          hipp_stats[[paste(hipp_var, grp, sep = "_")]] <- data.frame(
            Region = hipp_name,
            ST_Code = hipp_var,
            Age_Group = grp,
            N = nrow(sub_dat),
            Beta = coef_row["Estimate"],
            SE = coef_row["Std. Error"],
            t_value = coef_row["t value"],
            P_value = coef_row["Pr(>|t|)"],
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) NULL)
    }
  }
  
  if (length(hipp_stats) > 0) {
    hipp_res_df <- do.call(rbind, hipp_stats)
    rownames(hipp_res_df) <- NULL
    
    cat("\n  Hippocampus PRS association results:\n")
    cat("  ", paste(rep("-", 70), collapse = ""), "\n")
    
    for (i in 1:nrow(hipp_res_df)) {
      row <- hipp_res_df[i, ]
      sig <- ifelse(row$P_value < 0.001, "***",
                    ifelse(row$P_value < 0.01, "**",
                           ifelse(row$P_value < 0.05, "*", "")))
      cat(sprintf("  %s (%s): Beta=%.4f, t=%.3f, P=%.4f %s\n",
                  row$Region, row$Age_Group, row$Beta, row$t_value, row$P_value, sig))
    }
    
    fwrite(hipp_res_df, file.path(output_dir_p13, "Hippocampus_PRS_Association_Results.csv"))
    cat("\n  Saved: Hippocampus_PRS_Association_Results.csv\n")
  }
  
  # Visualization: Hippocampal volume vs PRS
  if (length(hipp_available) >= 2) {
    
    # Calculate mean hippocampal volume
    hipp_plot_data <- merged_data[, c("ST29SV", "ST88SV", "PRS_EOAD_Oligo", 
                                       "AGE", "PTGENDER", "Age_Group", icv_col), with = FALSE]
    hipp_plot_data <- na.omit(hipp_plot_data)
    hipp_plot_data$Hipp_Mean <- (hipp_plot_data$ST29SV + hipp_plot_data$ST88SV) / 2
    hipp_plot_data$Hipp_Norm <- hipp_plot_data$Hipp_Mean / hipp_plot_data[[icv_col]] * 1000
    hipp_plot_data$PRS_z <- as.numeric(scale(hipp_plot_data$PRS_EOAD_Oligo))
    hipp_plot_data$Age_Group <- factor(hipp_plot_data$Age_Group, 
                                        levels = c("Younger (<70)", "Older (≥70)"))
    
    p_hipp <- ggplot(hipp_plot_data, aes(x = PRS_z, y = Hipp_Norm, color = Age_Group)) +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
      scale_color_manual(values = c("Younger (<70)" = col_young, 
                                    "Older (≥70)" = col_old)) +
      facet_wrap(~Age_Group, scales = "fixed") +
      labs(
        title = "B. Oligodendrocyte Genetic Risk and Hippocampal Volume",
        subtitle = "Hypothesis: Higher Oligo-PRS -> Smaller hippocampal volume",
        x = "PRS_EOAD_Oligo (z-score)",
        y = "Normalized Hippocampal Volume",
        color = "Age Group"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "grey95")
      )
    
    ggsave(file.path(output_dir_p13, "Figure7B_Hippocampus_PRS_Validation.pdf"), 
           p_hipp, width = 10, height = 5, dpi = 300)
    ggsave(file.path(output_dir_p13, "Figure7B_Hippocampus_PRS_Validation.png"), 
           p_hipp, width = 10, height = 5, dpi = 300)
    cat("  Saved: Figure7B_Hippocampus_PRS_Validation.pdf/png\n")
  }
  
} else {
  cat("  Hippocampus variables or PRS not available. Skipping hippocampal analysis.\n")
}

# Analyze hippocampal subfields (if available)
subfield_res_df <- NULL
p_subfield <- NULL

if (length(subfield_available) > 0 && has_prs && has_age_group && !is.null(icv_col)) {
  
  cat("\n  Analyzing hippocampal subfields...\n")
  
  subfield_stats <- list()
  
  for (sf_var in subfield_available) {
    
    sf_name <- subfield_mapping$Region_Name[subfield_mapping$ST_Code == sf_var]
    if (length(sf_name) == 0) sf_name <- sf_var
    
    # Prepare data
    sf_cols <- c(sf_var, "PRS_EOAD_Oligo", "AGE", "PTGENDER", "Age_Group", icv_col)
    sf_data <- merged_data[, sf_cols, with = FALSE]
    sf_data <- na.omit(sf_data)
    
    if (nrow(sf_data) < 50) next
    
    # ICV normalization
    sf_data$Vol_Norm <- sf_data[[sf_var]] / sf_data[[icv_col]] * 1000
    sf_data$PRS_z <- as.numeric(scale(sf_data$PRS_EOAD_Oligo))
    
    # Full sample analysis
    tryCatch({
      mod <- lm(Vol_Norm ~ PRS_z + AGE + PTGENDER, data = sf_data)
      res <- summary(mod)
      
      if ("PRS_z" %in% rownames(coef(res))) {
        coef_row <- coef(res)["PRS_z", ]
        
        subfield_stats[[sf_var]] <- data.frame(
          Region = sf_name,
          ST_Code = sf_var,
          N = nrow(sf_data),
          Beta = coef_row["Estimate"],
          SE = coef_row["Std. Error"],
          t_value = coef_row["t value"],
          P_value = coef_row["Pr(>|t|)"],
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) NULL)
  }
  
  if (length(subfield_stats) > 0) {
    subfield_res_df <- do.call(rbind, subfield_stats)
    rownames(subfield_res_df) <- NULL
    subfield_res_df$P_FDR <- p.adjust(subfield_res_df$P_value, method = "fdr")
    
    cat("\n  Hippocampal subfield PRS association results:\n")
    cat("  ", paste(rep("-", 70), collapse = ""), "\n")
    
    # Sort by t value
    subfield_res_df <- subfield_res_df[order(subfield_res_df$t_value), ]
    
    for (i in 1:nrow(subfield_res_df)) {
      row <- subfield_res_df[i, ]
      sig <- ifelse(row$P_FDR < 0.05, "*FDR", 
                    ifelse(row$P_value < 0.05, "*", ""))
      cat(sprintf("  %s: Beta=%.4f, t=%.3f, P=%.4f, P_FDR=%.4f %s\n",
                  row$Region, row$Beta, row$t_value, row$P_value, row$P_FDR, sig))
    }
    
    fwrite(subfield_res_df, file.path(output_dir_p13, "Hippocampal_Subfield_PRS_Results.csv"))
    cat("\n  Saved: Hippocampal_Subfield_PRS_Results.csv\n")
    
    # Visualization: Subfield effect bar plot
    p_subfield <- ggplot(subfield_res_df, aes(x = reorder(Region, t_value), y = t_value)) +
      geom_bar(stat = "identity", aes(fill = t_value < 0), width = 0.7) +
      geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "grey50") +
      scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"), guide = "none") +
      coord_flip() +
      labs(
        title = "C. Hippocampal Subfield Specificity",
        subtitle = "PRS_Oligo effect on individual subfields (dashed = P<0.05)",
        x = "",
        y = "t-statistic (PRS effect)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        axis.text.y = element_text(size = 10)
      )
    
    ggsave(file.path(output_dir_p13, "Figure7C_Subfield_Specificity.pdf"), 
           p_subfield, width = 8, height = 6, dpi = 300)
    ggsave(file.path(output_dir_p13, "Figure7C_Subfield_Specificity.png"), 
           p_subfield, width = 8, height = 6, dpi = 300)
    cat("  Saved: Figure7C_Subfield_Specificity.pdf/png\n")
  }
  
} else {
  cat("  Hippocampal subfield data not available.\n")
}

# ==============================================================================
# Part 13.3: Subcortical Structure Analysis (Supplementary Validation)
# ==============================================================================
# Variables: ST12SV/ST71SV (Amygdala), ST16SV/ST75SV (Caudate), 
#            ST53SV/ST112SV (Putamen), ST42SV/ST101SV (Pallidum), 
#            ST61SV/ST120SV (Thalamus), ST11SV/ST70SV (Accumbens)
# ==============================================================================

cat("\n========================================\n")
cat("13.3 Subcortical Structure Analysis\n")
cat("========================================\n")

p_subcort <- NULL
p_subcort_age <- NULL
subcort_res_df <- NULL

# Subcortical structure mapping
subcort_mapping <- data.frame(
  ST_Code = c("ST29SV", "ST88SV", "ST12SV", "ST71SV", "ST16SV", "ST75SV",
              "ST53SV", "ST112SV", "ST42SV", "ST101SV", "ST61SV", "ST120SV",
              "ST11SV", "ST70SV"),
  Region_Name = c("Left Hippocampus", "Right Hippocampus",
                  "Left Amygdala", "Right Amygdala",
                  "Left Caudate", "Right Caudate",
                  "Left Putamen", "Right Putamen",
                  "Left Pallidum", "Right Pallidum",
                  "Left Thalamus", "Right Thalamus",
                  "Left Accumbens", "Right Accumbens"),
  Structure = c("Hippocampus", "Hippocampus", "Amygdala", "Amygdala",
                "Caudate", "Caudate", "Putamen", "Putamen",
                "Pallidum", "Pallidum", "Thalamus", "Thalamus",
                "Accumbens", "Accumbens"),
  Hemi = c("Left", "Right", "Left", "Right", "Left", "Right",
           "Left", "Right", "Left", "Right", "Left", "Right",
           "Left", "Right"),
  stringsAsFactors = FALSE
)

# Check available variables
subcort_available <- subcort_mapping$ST_Code[subcort_mapping$ST_Code %in% colnames(merged_data)]
cat(sprintf("  Subcortical variables available: %d of %d\n", 
            length(subcort_available), nrow(subcort_mapping)))

if (length(subcort_available) > 0 && has_prs && has_age_group && !is.null(icv_col)) {
  
  cat("\n  Analyzing subcortical structures...\n")
  
  subcort_stats <- list()
  
  for (sc_var in subcort_available) {
    
    sc_info <- subcort_mapping[subcort_mapping$ST_Code == sc_var, ]
    sc_name <- sc_info$Region_Name
    
    # Prepare data
    sc_cols <- c(sc_var, "PRS_EOAD_Oligo", "AGE", "PTGENDER", "Age_Group", icv_col)
    sc_data <- merged_data[, sc_cols, with = FALSE]
    sc_data <- na.omit(sc_data)
    
    if (nrow(sc_data) < 50) next
    
    # ICV normalization
    sc_data$Vol_Norm <- sc_data[[sc_var]] / sc_data[[icv_col]] * 1000
    sc_data$PRS_z <- as.numeric(scale(sc_data$PRS_EOAD_Oligo))
    sc_data$Age_Group <- factor(sc_data$Age_Group, 
                                 levels = c("Younger (<70)", "Older (≥70)"))
    
    # Age-stratified analysis
    for (grp in c("Younger (<70)", "Older (≥70)", "All")) {
      
      if (grp == "All") {
        sub_dat <- sc_data
      } else {
        sub_dat <- sc_data[sc_data$Age_Group == grp, ]
      }
      
      if (nrow(sub_dat) < 30) next
      
      tryCatch({
        mod <- lm(Vol_Norm ~ PRS_z + AGE + PTGENDER, data = sub_dat)
        res <- summary(mod)
        
        if ("PRS_z" %in% rownames(coef(res))) {
          coef_row <- coef(res)["PRS_z", ]
          
          subcort_stats[[paste(sc_var, grp, sep = "_")]] <- data.frame(
            Region = sc_name,
            Structure = sc_info$Structure,
            Hemi = sc_info$Hemi,
            ST_Code = sc_var,
            Age_Group = grp,
            N = nrow(sub_dat),
            Beta = coef_row["Estimate"],
            SE = coef_row["Std. Error"],
            t_value = coef_row["t value"],
            P_value = coef_row["Pr(>|t|)"],
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) NULL)
    }
  }
  
  if (length(subcort_stats) > 0) {
    subcort_res_df <- do.call(rbind, subcort_stats)
    rownames(subcort_res_df) <- NULL
    
    # Calculate FDR by structure and age group
    structures <- unique(subcort_res_df$Structure)
    subcort_res_df$P_FDR <- NA
    for (st in structures) {
      for (ag in c("Younger (<70)", "Older (≥70)", "All")) {
        idx <- subcort_res_df$Structure == st & subcort_res_df$Age_Group == ag
        if (sum(idx) > 0) {
          subcort_res_df$P_FDR[idx] <- p.adjust(subcort_res_df$P_value[idx], method = "fdr")
        }
      }
    }
    
    cat("\n  Subcortical PRS association results (All samples):\n")
    cat("  ", paste(rep("-", 70), collapse = ""), "\n")
    
    all_results <- subcort_res_df[subcort_res_df$Age_Group == "All", ]
    all_results <- all_results[order(all_results$t_value), ]
    
    for (i in 1:nrow(all_results)) {
      row <- all_results[i, ]
      sig <- ifelse(row$P_value < 0.001, "***",
                    ifelse(row$P_value < 0.01, "**",
                           ifelse(row$P_value < 0.05, "*", "")))
      cat(sprintf("  %s: Beta=%.4f, t=%.3f, P=%.4f %s\n",
                  row$Region, row$Beta, row$t_value, row$P_value, sig))
    }
    
    fwrite(subcort_res_df, file.path(output_dir_p13, "Subcortical_PRS_Association_Results.csv"))
    cat("\n  Saved: Subcortical_PRS_Association_Results.csv\n")
    
    # Visualization: Subcortical structure effects (summarized by structure)
    struct_summary <- aggregate(
      cbind(Beta, t_value, P_value) ~ Structure + Age_Group,
      data = subcort_res_df,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    
    # Plot only full sample results
    struct_all <- struct_summary[struct_summary$Age_Group == "All", ]
    
    p_subcort <- ggplot(struct_all, aes(x = reorder(Structure, t_value), y = t_value)) +
      geom_bar(stat = "identity", aes(fill = t_value < 0), width = 0.7) +
      geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "grey50") +
      scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"), guide = "none") +
      coord_flip() +
      labs(
        title = "D. Subcortical Structure Effects",
        subtitle = "PRS_Oligo effect on subcortical volumes (dashed = P<0.05)",
        x = "",
        y = "t-statistic (PRS effect)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        axis.text.y = element_text(size = 11)
      )
    
    ggsave(file.path(output_dir_p13, "Figure7D_Subcortical_Effects.pdf"), 
           p_subcort, width = 8, height = 5, dpi = 300)
    ggsave(file.path(output_dir_p13, "Figure7D_Subcortical_Effects.png"), 
           p_subcort, width = 8, height = 5, dpi = 300)
    cat("  Saved: Figure7D_Subcortical_Effects.pdf/png\n")
    
    # Age-stratified comparison plot
    struct_stratified <- struct_summary[struct_summary$Age_Group != "All", ]
    struct_stratified$Age_Group <- factor(struct_stratified$Age_Group, 
                                           levels = c("Younger (<70)", "Older (≥70)"))
    
    p_subcort_age <- ggplot(struct_stratified, 
                            aes(x = reorder(Structure, t_value), y = t_value, fill = Age_Group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "grey50") +
      scale_fill_manual(values = c("Younger (<70)" = col_young, "Older (≥70)" = col_old)) +
      coord_flip() +
      labs(
        title = "E. Age-Stratified Subcortical Effects",
        subtitle = "PRS_Oligo effect by age group",
        x = "",
        y = "t-statistic (PRS effect)",
        fill = "Age Group"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        axis.text.y = element_text(size = 11),
        legend.position = "bottom"
      )
    
    ggsave(file.path(output_dir_p13, "Figure7E_Subcortical_Age_Stratified.pdf"), 
           p_subcort_age, width = 9, height = 6, dpi = 300)
    ggsave(file.path(output_dir_p13, "Figure7E_Subcortical_Age_Stratified.png"), 
           p_subcort_age, width = 9, height = 6, dpi = 300)
    cat("  Saved: Figure7E_Subcortical_Age_Stratified.pdf/png\n")
  }
  
} else {
  cat("  Subcortical variables or PRS not available. Skipping subcortical analysis.\n")
}

# ==============================================================================
# Part 13.4: Subtype-Specific WMH Analysis
# ==============================================================================
# Purpose: Validate whether Oligo-Driven subtype has higher WMH burden
# ==============================================================================

cat("\n========================================\n")
cat("13.4 Subtype-Specific WMH Analysis\n")
cat("========================================\n")

p_wmh_subtype <- NULL
wmh_subtype_df <- NULL

has_subtype <- "Subtype" %in% colnames(merged_data)

if (has_wmh && has_subtype && !is.null(icv_col)) {
  
  # Prepare data
  wmh_sub_cols <- c("ST128SV", "Subtype", "AGE", "PTGENDER", icv_col)
  if (has_age_group) wmh_sub_cols <- c(wmh_sub_cols, "Age_Group")
  
  wmh_sub_data <- merged_data[, wmh_sub_cols, with = FALSE]
  wmh_sub_data <- na.omit(wmh_sub_data)
  
  cat(sprintf("  Sample size: N=%d\n", nrow(wmh_sub_data)))
  cat("  Subtype distribution:\n")
  print(table(wmh_sub_data$Subtype))
  
  # Log transform WMH
  wmh_sub_data$WMH_Log <- log(wmh_sub_data$ST128SV + 1)
  wmh_sub_data$ICV_z <- as.numeric(scale(wmh_sub_data[[icv_col]]))
  
  # Subtype comparison
  cat("\n  Subtype WMH comparison:\n")
  cat("  ", paste(rep("-", 60), collapse = ""), "\n")
  
  subtypes <- unique(wmh_sub_data$Subtype)
  wmh_subtype_stats <- list()
  
  for (st in subtypes) {
    st_data <- wmh_sub_data[wmh_sub_data$Subtype == st, ]
    wmh_subtype_stats[[st]] <- data.frame(
      Subtype = st,
      N = nrow(st_data),
      WMH_Mean = mean(st_data$ST128SV, na.rm = TRUE),
      WMH_SD = sd(st_data$ST128SV, na.rm = TRUE),
      WMH_Log_Mean = mean(st_data$WMH_Log, na.rm = TRUE),
      WMH_Log_SD = sd(st_data$WMH_Log, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %s: N=%d, WMH=%.1f +/- %.1f, log(WMH)=%.2f +/- %.2f\n",
                st, nrow(st_data), 
                mean(st_data$ST128SV), sd(st_data$ST128SV),
                mean(st_data$WMH_Log), sd(st_data$WMH_Log)))
  }
  
  wmh_subtype_df <- do.call(rbind, wmh_subtype_stats)
  rownames(wmh_subtype_df) <- NULL
  
  # ANOVA test
  cat("\n  ANOVA test:\n")
  tryCatch({
    mod_anova <- lm(WMH_Log ~ Subtype + AGE + PTGENDER + ICV_z, data = wmh_sub_data)
    anova_res <- anova(mod_anova)
    cat(sprintf("  Subtype effect: F=%.2f, P=%.4f\n", 
                anova_res["Subtype", "F value"], anova_res["Subtype", "Pr(>F)"]))
    
    # Post-hoc comparison (if Oligo-Driven exists)
    if ("Oligo-Driven" %in% subtypes) {
      cat("\n  Pairwise comparisons (vs Oligo-Driven):\n")
      
      # Set Oligo-Driven as reference
      wmh_sub_data$Subtype <- relevel(factor(wmh_sub_data$Subtype), ref = "Oligo-Driven")
      mod_post <- lm(WMH_Log ~ Subtype + AGE + PTGENDER + ICV_z, data = wmh_sub_data)
      res_post <- summary(mod_post)
      
      subtype_terms <- grep("^Subtype", rownames(coef(res_post)), value = TRUE)
      for (term in subtype_terms) {
        coef_row <- coef(res_post)[term, ]
        st_name <- gsub("Subtype", "", term)
        sig <- ifelse(coef_row["Pr(>|t|)"] < 0.001, "***",
                      ifelse(coef_row["Pr(>|t|)"] < 0.01, "**",
                             ifelse(coef_row["Pr(>|t|)"] < 0.05, "*", "")))
        cat(sprintf("  %s vs Oligo-Driven: Beta=%.3f, t=%.2f, P=%.4f %s\n",
                    st_name, coef_row["Estimate"], coef_row["t value"], 
                    coef_row["Pr(>|t|)"], sig))
      }
    }
  }, error = function(e) {
    cat(sprintf("  ANOVA error: %s\n", e$message))
  })
  
  fwrite(wmh_subtype_df, file.path(output_dir_p13, "WMH_Subtype_Summary.csv"))
  cat("\n  Saved: WMH_Subtype_Summary.csv\n")
  
  # Visualization: Box plot
  p_wmh_subtype <- ggplot(wmh_sub_data, aes(x = Subtype, y = WMH_Log, fill = Subtype)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "F. White Matter Hyperintensities by Biological Subtype",
      subtitle = "Hypothesis: Oligo-Driven subtype shows highest WMH burden",
      x = "Biological Subtype",
      y = "log(WMH Volume + 1)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey40"),
      axis.text.x = element_text(angle = 15, hjust = 1, size = 10)
    )
  
  ggsave(file.path(output_dir_p13, "Figure7F_WMH_by_Subtype.pdf"), 
         p_wmh_subtype, width = 8, height = 6, dpi = 300)
  ggsave(file.path(output_dir_p13, "Figure7F_WMH_by_Subtype.png"), 
         p_wmh_subtype, width = 8, height = 6, dpi = 300)
  cat("  Saved: Figure7F_WMH_by_Subtype.pdf/png\n")
  
} else {
  cat("  WMH or Subtype variable not available. Skipping subtype WMH analysis.\n")
}

# ==============================================================================
# Part 13.5: Combined Figure 7 Generation
# ==============================================================================

cat("\n========================================\n")
cat("13.5 Generating Combined Figure 7\n")
cat("========================================\n")

# Load patchwork for combining plots
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

# Collect all available plots
fig7_plots <- list()

if (!is.null(p_wmh)) {
  fig7_plots[["A"]] <- p_wmh
  cat("  Panel A (WMH-PRS): Available\n")
}

if (!is.null(p_hipp)) {
  fig7_plots[["B"]] <- p_hipp
  cat("  Panel B (Hippocampus-PRS): Available\n")
}

if (!is.null(p_subfield)) {
  fig7_plots[["C"]] <- p_subfield
  cat("  Panel C (Subfield specificity): Available\n")
}

if (!is.null(p_subcort)) {
  fig7_plots[["D"]] <- p_subcort
  cat("  Panel D (Subcortical effects): Available\n")
}

if (!is.null(p_subcort_age)) {
  fig7_plots[["E"]] <- p_subcort_age
  cat("  Panel E (Subcortical age-stratified): Available\n")
}

if (!is.null(p_wmh_subtype)) {
  fig7_plots[["F"]] <- p_wmh_subtype
  cat("  Panel F (WMH by subtype): Available\n")
}

# Combine plots
if (length(fig7_plots) >= 2) {
  
  cat(sprintf("\n  Combining %d panels into Figure 7...\n", length(fig7_plots)))
  
  # Adjust layout based on number of available panels
  if (length(fig7_plots) == 2) {
    fig7_combined <- fig7_plots[[1]] / fig7_plots[[2]]
    fig_width <- 12
    fig_height <- 10
  } else if (length(fig7_plots) == 3) {
    fig7_combined <- fig7_plots[[1]] / (fig7_plots[[2]] | fig7_plots[[3]])
    fig_width <- 14
    fig_height <- 12
  } else if (length(fig7_plots) == 4) {
    fig7_combined <- (fig7_plots[[1]] | fig7_plots[[2]]) / (fig7_plots[[3]] | fig7_plots[[4]])
    fig_width <- 16
    fig_height <- 12
  } else if (length(fig7_plots) >= 5) {
    # 3-row layout
    row1 <- fig7_plots[[1]] | fig7_plots[[2]]
    row2 <- fig7_plots[[3]] | fig7_plots[[4]]
    if (length(fig7_plots) >= 6) {
      row3 <- fig7_plots[[5]] | fig7_plots[[6]]
    } else {
      row3 <- fig7_plots[[5]] | plot_spacer()
    }
    fig7_combined <- row1 / row2 / row3
    fig_width <- 16
    fig_height <- 18
  }
  
  fig7_combined <- fig7_combined +
    plot_annotation(
      title = "Figure 7. Mechanistic Validation: Oligodendrocyte Genetic Risk and Brain Structure",
      subtitle = "Validating the oligodendrocyte hypothesis through imaging-genetics associations",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30")
      )
    )
  
  ggsave(file.path(output_dir_p13, "Figure7_Mechanistic_Validation_Combined.pdf"), 
         fig7_combined, width = fig_width, height = fig_height, dpi = 300)
  ggsave(file.path(output_dir_p13, "Figure7_Mechanistic_Validation_Combined.png"), 
         fig7_combined, width = fig_width, height = fig_height, dpi = 300)
  cat("  Saved: Figure7_Mechanistic_Validation_Combined.pdf/png\n")
  
} else {
  cat("  Not enough panels available for combined figure.\n")
}

# ==============================================================================
# Part 13.6: Summary Report
# ==============================================================================

cat("\n========================================\n")
cat("Part 13 Summary Report\n")
cat("========================================\n")

cat("\n  Output files generated:\n")
cat("  ", paste(rep("-", 50), collapse = ""), "\n")

output_files <- list.files(output_dir_p13, pattern = "^(Figure7|WMH|Hippocampus|Subcortical)", full.names = FALSE)
for (f in output_files) {
  cat(sprintf("  - %s\n", f))
}

cat("\n  Key findings:\n")
cat("  ", paste(rep("-", 50), collapse = ""), "\n")

# WMH results summary
if (!is.null(wmh_res_df)) {
  wmh_all <- wmh_res_df[wmh_res_df$Age_Group == "All", ]
  if (nrow(wmh_all) > 0) {
    cat(sprintf("  WMH-PRS (All): Beta=%.4f, P=%.4f\n", wmh_all$Beta, wmh_all$P_value))
    if (wmh_all$P_value < 0.05 && wmh_all$Beta > 0) {
      cat("    -> SUPPORTS oligodendrocyte hypothesis (higher PRS = more WMH)\n")
    }
  }
}

# Hippocampus results summary
if (!is.null(hipp_res_df)) {
  hipp_all <- hipp_res_df[hipp_res_df$Age_Group == "All", ]
  if (nrow(hipp_all) > 0) {
    cat(sprintf("  Hippocampus-PRS (All): Mean Beta=%.4f\n", mean(hipp_all$Beta)))
    if (mean(hipp_all$Beta) < 0) {
      cat("    -> Higher PRS associated with smaller hippocampal volume\n")
    }
  }
}

# Subcortical results summary
if (!is.null(subcort_res_df)) {
  subcort_all <- subcort_res_df[subcort_res_df$Age_Group == "All", ]
  sig_regions <- subcort_all[subcort_all$P_value < 0.05, ]
  if (nrow(sig_regions) > 0) {
    cat(sprintf("  Significant subcortical regions: %d\n", nrow(sig_regions)))
    for (i in 1:min(5, nrow(sig_regions))) {
      cat(sprintf("    - %s: t=%.2f, P=%.4f\n", 
                  sig_regions$Region[i], sig_regions$t_value[i], sig_regions$P_value[i]))
    }
  }
}

cat("\n================================================================================\n")
cat("ADNI Complete Validation Analysis (Part 1-11 + Part 13) Complete!\n")
cat("================================================================================\n")
