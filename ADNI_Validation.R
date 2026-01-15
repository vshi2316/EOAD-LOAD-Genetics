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
  sumstats = data.frame(chr = snp_info_eoad$chr, pos = snp_info_eoad$pos,
                        a0 = snp_info_eoad$a0, a1 = snp_info_eoad$a1,
                        beta = prs_weights$eoad_full),
  info_snp = adni_map, strand_flip = TRUE,
  join_by_pos = TRUE, remove_dups = TRUE
)

match_load <- snp_match(
  sumstats = data.frame(chr = snp_info_load$chr, pos = snp_info_load$pos,
                        a0 = snp_info_load$a0, a1 = snp_info_load$a1,
                        beta = prs_weights$load_full),
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
  
  fs_vars <- c("RID", "ST10CV", "ST128SV", "ST29SV", "ST88SV", 
               "ST24CV", "ST83CV", "ST37SV", "ST96SV", "ST44CV", "ST103CV")
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

## Section 4.1: Load CSF Data (AlzBio3 Platform)
## Per methods: CSF biomarkers restricted to AlzBio3 platform, ADNI phases 1/GO/2
csf_file <- "UPENNBIOMK_MASTER.csv"

if (file.exists(csf_file)) {
  csf_data <- fread(csf_file)
  
  # Filter to AlzBio3 platform and phases 1/GO/2
  csf_alzbio3 <- csf_data[csf_data$BATCH == "AlzBio3" | 
                            grepl("AlzBio3", csf_data$RUNDATE, ignore.case = TRUE), ]
  
  # Restrict to ADNI phases 1, GO, 2
  if ("PHASE" %in% colnames(csf_alzbio3)) {
    csf_alzbio3 <- csf_alzbio3[csf_alzbio3$PHASE %in% c("ADNI1", "ADNIGO", "ADNI2"), ]
  }
  
  # Extract baseline CSF
  csf_bl <- csf_alzbio3[csf_alzbio3$VISCODE2 == "bl" | csf_alzbio3$VISCODE == "bl", ]
  
  # Select CSF variables
  csf_vars <- c("RID", "ABETA", "TAU", "PTAU")
  csf_vars_available <- csf_vars[csf_vars %in% colnames(csf_bl)]
  
  if (length(csf_vars_available) > 1) {
    csf_subset <- csf_bl[, ..csf_vars_available]
    csf_subset <- csf_subset[!duplicated(csf_subset$RID), ]
    
    # Convert to numeric
    for (col in csf_vars_available[-1]) {
      csf_subset[[col]] <- as.numeric(as.character(csf_subset[[col]]))
    }
    
    # Calculate derived ratios per methods
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

## Section 4.2: CSF Amyloid Positivity (Abeta < 192 pg/mL)
if ("ABETA" %in% colnames(merged_data)) {
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
                    "Ventricle_Normalized", "CA1_Bilateral", "Subiculum_Bilateral")

## Section 5.2: PRS-Biomarker Association Function
## Covariates: age, sex, APOE4 count, education (PTEDUCAT) per methods
run_prs_biomarker_association <- function(data, prs_var, biomarker_var) {
  
  # Check variable availability
  required_vars <- c(prs_var, biomarker_var, "AGE", "Sex_Numeric", "APOE4_Count", "PTEDUCAT")
  if (!all(required_vars %in% colnames(data))) {
    return(NULL)
  }
  
  # Remove missing values
  analysis_data <- data[complete.cases(data[, required_vars]), ]
  
  if (nrow(analysis_data) < 50) return(NULL)
  
  # Linear regression with covariates
  formula_str <- paste0(biomarker_var, " ~ ", prs_var, 
                        " + AGE + Sex_Numeric + APOE4_Count + PTEDUCAT")
  model <- lm(as.formula(formula_str), data = analysis_data)
  
  coef_summary <- summary(model)$coefficients
  
  if (prs_var %in% rownames(coef_summary)) {
    result <- data.frame(
      PRS = prs_var,
      Biomarker = biomarker_var,
      Beta = coef_summary[prs_var, "Estimate"],
      SE = coef_summary[prs_var, "Std. Error"],
      t_value = coef_summary[prs_var, "t value"],
      P_value = coef_summary[prs_var, "Pr(>|t|)"],
      N = nrow(analysis_data),
      R2 = summary(model)$r.squared,
      stringsAsFactors = FALSE
    )
    return(result)
  }
  return(NULL)
}

## Section 5.3: Run All PRS-Biomarker Associations
prs_biomarker_results <- list()

for (prs in prs_vars) {
  for (biomarker in biomarker_vars) {
    if (prs %in% colnames(merged_data) && biomarker %in% colnames(merged_data)) {
      result <- run_prs_biomarker_association(merged_data, prs, biomarker)
      if (!is.null(result)) {
        prs_biomarker_results[[paste(prs, biomarker, sep = "_")]] <- result
      }
    }
  }
}

prs_biomarker_df <- do.call(rbind, prs_biomarker_results)

if (!is.null(prs_biomarker_df) && nrow(prs_biomarker_df) > 0) {
  prs_biomarker_df$FDR <- p.adjust(prs_biomarker_df$P_value, method = "BH")
  prs_biomarker_df <- prs_biomarker_df[order(prs_biomarker_df$P_value), ]
  fwrite(prs_biomarker_df, "results/tables/ADNI_PRS_Biomarker_Associations.csv")
}

# ==============================================================================
# Part 6: Age-Stratified Analysis (<70 vs ≥70)
# ==============================================================================

## Section 6.1: Age-Stratified PRS-Biomarker Association
run_age_stratified_analysis <- function(data, prs_var, biomarker_var, age_group) {
  
  subset_data <- data[data$Age_Group == age_group, ]
  required_vars <- c(prs_var, biomarker_var, "AGE", "Sex_Numeric", "APOE4_Count", "PTEDUCAT")
  
  if (!all(required_vars %in% colnames(subset_data))) return(NULL)
  
  analysis_data <- subset_data[complete.cases(subset_data[, required_vars]), ]
  
  if (nrow(analysis_data) < 30) return(NULL)
  
  formula_str <- paste0(biomarker_var, " ~ ", prs_var, 
                        " + AGE + Sex_Numeric + APOE4_Count + PTEDUCAT")
  model <- lm(as.formula(formula_str), data = analysis_data)
  
  coef_summary <- summary(model)$coefficients
  
  if (prs_var %in% rownames(coef_summary)) {
    result <- data.frame(
      Age_Group = age_group,
      PRS = prs_var,
      Biomarker = biomarker_var,
      Beta = coef_summary[prs_var, "Estimate"],
      SE = coef_summary[prs_var, "Std. Error"],
      P_value = coef_summary[prs_var, "Pr(>|t|)"],
      N = nrow(analysis_data),
      stringsAsFactors = FALSE
    )
    return(result)
  }
  return(NULL)
}

## Section 6.2: Run Age-Stratified Analyses
age_stratified_results <- list()
age_groups <- c("Younger (<70)", "Older (≥70)")

for (age_grp in age_groups) {
  for (prs in prs_vars) {
    for (biomarker in biomarker_vars) {
      if (all(c(prs, biomarker, "Age_Group") %in% colnames(merged_data))) {
        result <- run_age_stratified_analysis(merged_data, prs, biomarker, age_grp)
        if (!is.null(result)) {
          age_stratified_results[[paste(age_grp, prs, biomarker, sep = "_")]] <- result
        }
      }
    }
  }
}

age_stratified_df <- do.call(rbind, age_stratified_results)

if (!is.null(age_stratified_df) && nrow(age_stratified_df) > 0) {
  age_stratified_df$FDR <- p.adjust(age_stratified_df$P_value, method = "BH")
  fwrite(age_stratified_df, "results/tables/ADNI_Age_Stratified_Associations.csv")
}

## Section 6.3: Interaction Analysis (PRS × Age Group)
run_interaction_analysis <- function(data, prs_var, biomarker_var) {
  
  required_vars <- c(prs_var, biomarker_var, "Age_Group", "AGE", 
                     "Sex_Numeric", "APOE4_Count", "PTEDUCAT")
  
  if (!all(required_vars %in% colnames(data))) return(NULL)
  
  analysis_data <- data[complete.cases(data[, required_vars]), ]
  
  if (nrow(analysis_data) < 100) return(NULL)
  
  # Model with interaction term
  formula_str <- paste0(biomarker_var, " ~ ", prs_var, " * Age_Group + ",
                        "AGE + Sex_Numeric + APOE4_Count + PTEDUCAT")
  model <- lm(as.formula(formula_str), data = analysis_data)
  
  coef_summary <- summary(model)$coefficients
  interaction_term <- paste0(prs_var, ":Age_GroupOlder (≥70)")
  
  if (interaction_term %in% rownames(coef_summary)) {
    result <- data.frame(
      PRS = prs_var,
      Biomarker = biomarker_var,
      Interaction_Beta = coef_summary[interaction_term, "Estimate"],
      Interaction_SE = coef_summary[interaction_term, "Std. Error"],
      Interaction_P = coef_summary[interaction_term, "Pr(>|t|)"],
      N = nrow(analysis_data),
      stringsAsFactors = FALSE
    )
    return(result)
  }
  return(NULL)
}

interaction_results <- list()

for (prs in prs_vars) {
  for (biomarker in biomarker_vars) {
    if (all(c(prs, biomarker) %in% colnames(merged_data))) {
      result <- run_interaction_analysis(merged_data, prs, biomarker)
      if (!is.null(result)) {
        interaction_results[[paste(prs, biomarker, sep = "_")]] <- result
      }
    }
  }
}

interaction_df <- do.call(rbind, interaction_results)

if (!is.null(interaction_df) && nrow(interaction_df) > 0) {
  interaction_df$FDR <- p.adjust(interaction_df$Interaction_P, method = "BH")
  fwrite(interaction_df, "results/tables/ADNI_PRS_Age_Interaction.csv")
}

# ==============================================================================
# Part 7: MCI-to-Dementia Conversion (Cox Regression)
# ==============================================================================

## Section 7.1: Prepare Longitudinal Data for Survival Analysis
## Per methods: Cox proportional hazards with Surv(Time_Years, Event)
adni_long <- fread("ADNIMERGE.csv")

# Calculate time to conversion and event status
mci_subjects <- adni_long[adni_long$DX_bl == "LMCI" | adni_long$DX_bl == "EMCI", ]

# Get first dementia diagnosis date
conversion_data <- mci_subjects %>%
  group_by(RID) %>%
  arrange(EXAMDATE) %>%
  summarise(
    Baseline_Date = first(EXAMDATE),
    Baseline_DX = first(DX_bl),
    Last_Date = last(EXAMDATE),
    Last_DX = last(DX),
    Converted = any(DX == "Dementia", na.rm = TRUE),
    Conversion_Date = ifelse(any(DX == "Dementia", na.rm = TRUE),
                             EXAMDATE[which(DX == "Dementia")[1]], NA),
    .groups = "drop"
  ) %>%
  as.data.frame()

# Calculate time in years
conversion_data$Baseline_Date <- as.Date(conversion_data$Baseline_Date)
conversion_data$Last_Date <- as.Date(conversion_data$Last_Date)
conversion_data$Conversion_Date <- as.Date(conversion_data$Conversion_Date)

conversion_data$Time_Years <- ifelse(
  conversion_data$Converted,
  as.numeric(difftime(conversion_data$Conversion_Date, 
                      conversion_data$Baseline_Date, units = "days")) / 365.25,
  as.numeric(difftime(conversion_data$Last_Date, 
                      conversion_data$Baseline_Date, units = "days")) / 365.25
)

conversion_data$Event <- as.integer(conversion_data$Converted)

# Remove subjects with invalid follow-up
conversion_data <- conversion_data[conversion_data$Time_Years > 0, ]

## Section 7.2: Merge with PRS and Covariates
survival_data <- merge(conversion_data, merged_data, by = "RID", all.x = TRUE)

## Section 7.3: Cox Proportional Hazards Regression
## Covariates: age, sex, APOE4 count, education per methods
run_cox_regression <- function(data, prs_var) {
  
  required_vars <- c(prs_var, "Time_Years", "Event", "AGE", 
                     "Sex_Numeric", "APOE4_Count", "PTEDUCAT")
  
  if (!all(required_vars %in% colnames(data))) return(NULL)
  
  analysis_data <- data[complete.cases(data[, required_vars]), ]
  
  if (nrow(analysis_data) < 50 || sum(analysis_data$Event) < 10) return(NULL)
  
  # Cox regression with Surv(Time_Years, Event)
  formula_str <- paste0("Surv(Time_Years, Event) ~ ", prs_var,
                        " + AGE + Sex_Numeric + APOE4_Count + PTEDUCAT")
  
  cox_model <- coxph(as.formula(formula_str), data = analysis_data)
  cox_summary <- summary(cox_model)
  
  coef_df <- as.data.frame(cox_summary$coefficients)
  
  if (prs_var %in% rownames(coef_df)) {
    result <- data.frame(
      PRS = prs_var,
      HR = exp(coef_df[prs_var, "coef"]),
      HR_Lower = exp(coef_df[prs_var, "coef"] - 1.96 * coef_df[prs_var, "se(coef)"]),
      HR_Upper = exp(coef_df[prs_var, "coef"] + 1.96 * coef_df[prs_var, "se(coef)"]),
      Beta = coef_df[prs_var, "coef"],
      SE = coef_df[prs_var, "se(coef)"],
      Z = coef_df[prs_var, "z"],
      P_value = coef_df[prs_var, "Pr(>|z|)"],
      N = nrow(analysis_data),
      N_Events = sum(analysis_data$Event),
      Concordance = cox_summary$concordance[1],
      stringsAsFactors = FALSE
    )
    return(result)
  }
  return(NULL)
}

## Section 7.4: Run Cox Regression for All PRS
cox_results <- list()

for (prs in prs_vars) {
  if (prs %in% colnames(survival_data)) {
    result <- run_cox_regression(survival_data, prs)
    if (!is.null(result)) {
      cox_results[[prs]] <- result
    }
  }
}

cox_df <- do.call(rbind, cox_results)

if (!is.null(cox_df) && nrow(cox_df) > 0) {
  cox_df$FDR <- p.adjust(cox_df$P_value, method = "BH")
  cox_df <- cox_df[order(cox_df$P_value), ]
  fwrite(cox_df, "results/tables/ADNI_Cox_MCI_Conversion.csv")
}

## Section 7.5: Age-Stratified Cox Regression
run_cox_age_stratified <- function(data, prs_var, age_group) {
  
  subset_data <- data[data$Age_Group == age_group, ]
  required_vars <- c(prs_var, "Time_Years", "Event", "AGE", 
                     "Sex_Numeric", "APOE4_Count", "PTEDUCAT")
  
  if (!all(required_vars %in% colnames(subset_data))) return(NULL)
  
  analysis_data <- subset_data[complete.cases(subset_data[, required_vars]), ]
  
  if (nrow(analysis_data) < 30 || sum(analysis_data$Event) < 5) return(NULL)
  
  formula_str <- paste0("Surv(Time_Years, Event) ~ ", prs_var,
                        " + AGE + Sex_Numeric + APOE4_Count + PTEDUCAT")
  
  cox_model <- coxph(as.formula(formula_str), data = analysis_data)
  cox_summary <- summary(cox_model)
  coef_df <- as.data.frame(cox_summary$coefficients)
  
  if (prs_var %in% rownames(coef_df)) {
    result <- data.frame(
      Age_Group = age_group,
      PRS = prs_var,
      HR = exp(coef_df[prs_var, "coef"]),
      HR_Lower = exp(coef_df[prs_var, "coef"] - 1.96 * coef_df[prs_var, "se(coef)"]),
      HR_Upper = exp(coef_df[prs_var, "coef"] + 1.96 * coef_df[prs_var, "se(coef)"]),
      P_value = coef_df[prs_var, "Pr(>|z|)"],
      N = nrow(analysis_data),
      N_Events = sum(analysis_data$Event),
      stringsAsFactors = FALSE
    )
    return(result)
  }
  return(NULL)
}

cox_age_results <- list()

for (age_grp in age_groups) {
  for (prs in prs_vars) {
    if (all(c(prs, "Age_Group") %in% colnames(survival_data))) {
      result <- run_cox_age_stratified(survival_data, prs, age_grp)
      if (!is.null(result)) {
        cox_age_results[[paste(age_grp, prs, sep = "_")]] <- result
      }
    }
  }
}

cox_age_df <- do.call(rbind, cox_age_results)

if (!is.null(cox_age_df) && nrow(cox_age_df) > 0) {
  cox_age_df$FDR <- p.adjust(cox_age_df$P_value, method = "BH")
  fwrite(cox_age_df, "results/tables/ADNI_Cox_Age_Stratified.csv")
}

# ==============================================================================
# Part 8: Unsupervised K-means Clustering (k=3)
# ==============================================================================

## Section 8.1: Prepare Clustering Data
## Per methods: k-means with k=3, silhouette analysis for validation
clustering_vars <- c("PRS_EOAD_Global", "PRS_LOAD_Global", 
                     "PRS_EOAD_TCell", "PRS_EOAD_Microglia", 
                     "PRS_LOAD_TCell", "PRS_LOAD_Microglia")

clustering_vars_available <- clustering_vars[clustering_vars %in% colnames(merged_data)]

if (length(clustering_vars_available) >= 3) {
  
  cluster_data <- merged_data[, c("PTID", clustering_vars_available)]
  cluster_data <- cluster_data[complete.cases(cluster_data), ]
  
  # Scale features for clustering
  cluster_matrix <- scale(cluster_data[, clustering_vars_available])
  rownames(cluster_matrix) <- cluster_data$PTID
  
  ## Section 8.2: Determine Optimal k (Silhouette Method)
  silhouette_scores <- sapply(2:6, function(k) {
    km <- kmeans(cluster_matrix, centers = k, nstart = 25, iter.max = 100)
    sil <- silhouette(km$cluster, dist(cluster_matrix))
    mean(sil[, 3])
  })
  
  optimal_k <- which.max(silhouette_scores) + 1
  
  ## Section 8.3: K-means Clustering (k=3 per methods)
  set.seed(42)
  km_result <- kmeans(cluster_matrix, centers = 3, nstart = 25, iter.max = 100)
  
  cluster_data$Cluster <- factor(km_result$cluster)
  
  # Calculate silhouette for k=3
  sil_k3 <- silhouette(km_result$cluster, dist(cluster_matrix))
  avg_silhouette <- mean(sil_k3[, 3])
  
  ## Section 8.4: Cluster Characterization
  cluster_summary <- cluster_data %>%
    group_by(Cluster) %>%
    summarise(
      N = n(),
      across(all_of(clustering_vars_available), 
             list(Mean = ~mean(., na.rm = TRUE), SD = ~sd(., na.rm = TRUE))),
      .groups = "drop"
    )
  
  fwrite(cluster_summary, "results/tables/ADNI_Cluster_Summary.csv")
  
  ## Section 8.5: Merge Cluster Assignments
  merged_data <- merge(merged_data, 
                       cluster_data[, c("PTID", "Cluster")], 
                       by = "PTID", all.x = TRUE)
  
  ## Section 8.6: Cluster-Biomarker Associations
  cluster_biomarker_results <- list()
  
  for (biomarker in biomarker_vars) {
    if (biomarker %in% colnames(merged_data) && "Cluster" %in% colnames(merged_data)) {
      
      analysis_data <- merged_data[!is.na(merged_data$Cluster) & 
                                     !is.na(merged_data[[biomarker]]), ]
      
      if (nrow(analysis_data) >= 30) {
        # ANOVA for cluster differences
        anova_result <- aov(as.formula(paste(biomarker, "~ Cluster")), 
                            data = analysis_data)
        anova_summary <- summary(anova_result)
        
        cluster_biomarker_results[[biomarker]] <- data.frame(
          Biomarker = biomarker,
          F_value = anova_summary[[1]]["Cluster", "F value"],
          P_value = anova_summary[[1]]["Cluster", "Pr(>F)"],
          N = nrow(analysis_data),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  cluster_biomarker_df <- do.call(rbind, cluster_biomarker_results)
  
  if (!is.null(cluster_biomarker_df) && nrow(cluster_biomarker_df) > 0) {
    cluster_biomarker_df$FDR <- p.adjust(cluster_biomarker_df$P_value, method = "BH")
    fwrite(cluster_biomarker_df, "results/tables/ADNI_Cluster_Biomarker_ANOVA.csv")
  }
}

# ==============================================================================
# Part 9: Longitudinal Cognitive Trajectory (LMM)
# ==============================================================================

## Section 9.1: Prepare Longitudinal Cognitive Data
cognitive_vars <- c("MMSE", "CDRSB", "ADAS11", "ADAS13")

adni_long <- fread("ADNIMERGE.csv")

# Calculate time from baseline
adni_long <- adni_long %>%
  group_by(RID) %>%
  arrange(EXAMDATE) %>%
  mutate(
    Baseline_Date = first(EXAMDATE),
    Time_Years = as.numeric(difftime(as.Date(EXAMDATE), 
                                     as.Date(Baseline_Date), 
                                     units = "days")) / 365.25
  ) %>%
  ungroup() %>%
  as.data.frame()

# Merge with PRS
long_data <- merge(adni_long, 
                   merged_data[, c("RID", prs_vars, "Age_Group", "APOE4_Count", 
                                   "Sex_Numeric", "PTEDUCAT")],
                   by = "RID", all.x = TRUE)

## Section 9.2: Linear Mixed Model Function
## Per methods: LMM with random intercept and slope
run_lmm_trajectory <- function(data, prs_var, cognitive_var) {
  
  required_vars <- c(prs_var, cognitive_var, "Time_Years", "RID",
                     "AGE", "Sex_Numeric", "APOE4_Count", "PTEDUCAT")
  
  if (!all(required_vars %in% colnames(data))) return(NULL)
  
  analysis_data <- data[complete.cases(data[, required_vars]), ]
  
  # Require minimum observations
  n_subjects <- length(unique(analysis_data$RID))
  n_obs <- nrow(analysis_data)
  
  if (n_subjects < 50 || n_obs < 100) return(NULL)
  
  # LMM: cognitive ~ PRS * Time + covariates + (1 + Time | RID)
  formula_str <- paste0(cognitive_var, " ~ ", prs_var, " * Time_Years + ",
                        "AGE + Sex_Numeric + APOE4_Count + PTEDUCAT + ",
                        "(1 + Time_Years | RID)")
  
  tryCatch({
    lmm_model <- lmer(as.formula(formula_str), data = analysis_data,
                      control = lmerControl(optimizer = "bobyqa"))
    
    lmm_summary <- summary(lmm_model)
    coef_df <- as.data.frame(lmm_summary$coefficients)
    
    interaction_term <- paste0(prs_var, ":Time_Years")
    
    if (interaction_term %in% rownames(coef_df)) {
      result <- data.frame(
        PRS = prs_var,
        Cognitive = cognitive_var,
        PRS_Beta = coef_df[prs_var, "Estimate"],
        PRS_SE = coef_df[prs_var, "Std. Error"],
        PRS_P = 2 * pnorm(-abs(coef_df[prs_var, "t value"])),
        Interaction_Beta = coef_df[interaction_term, "Estimate"],
        Interaction_SE = coef_df[interaction_term, "Std. Error"],
        Interaction_P = 2 * pnorm(-abs(coef_df[interaction_term, "t value"])),
        N_Subjects = n_subjects,
        N_Observations = n_obs,
        stringsAsFactors = FALSE
      )
      return(result)
    }
  }, error = function(e) {
    return(NULL)
  })
  
  return(NULL)
}

## Section 9.3: Run LMM for All PRS-Cognitive Combinations
lmm_results <- list()

for (prs in prs_vars) {
  for (cog in cognitive_vars) {
    if (all(c(prs, cog) %in% colnames(long_data))) {
      result <- run_lmm_trajectory(long_data, prs, cog)
      if (!is.null(result)) {
        lmm_results[[paste(prs, cog, sep = "_")]] <- result
      }
    }
  }
}

lmm_df <- do.call(rbind, lmm_results)

if (!is.null(lmm_df) && nrow(lmm_df) > 0) {
  lmm_df$Interaction_FDR <- p.adjust(lmm_df$Interaction_P, method = "BH")
  lmm_df <- lmm_df[order(lmm_df$Interaction_P), ]
  fwrite(lmm_df, "results/tables/ADNI_LMM_Cognitive_Trajectory.csv")
}

## Section 9.4: Age-Stratified LMM
run_lmm_age_stratified <- function(data, prs_var, cognitive_var, age_group) {
  
  subset_data <- data[data$Age_Group == age_group, ]
  required_vars <- c(prs_var, cognitive_var, "Time_Years", "RID",
                     "AGE", "Sex_Numeric", "APOE4_Count", "PTEDUCAT")
  
  if (!all(required_vars %in% colnames(subset_data))) return(NULL)
  
  analysis_data <- subset_data[complete.cases(subset_data[, required_vars]), ]
  
  n_subjects <- length(unique(analysis_data$RID))
  n_obs <- nrow(analysis_data)
  
  if (n_subjects < 30 || n_obs < 60) return(NULL)
  
  formula_str <- paste0(cognitive_var, " ~ ", prs_var, " * Time_Years + ",
                        "AGE + Sex_Numeric + APOE4_Count + PTEDUCAT + ",
                        "(1 + Time_Years | RID)")
  
  tryCatch({
    lmm_model <- lmer(as.formula(formula_str), data = analysis_data,
                      control = lmerControl(optimizer = "bobyqa"))
    
    lmm_summary <- summary(lmm_model)
    coef_df <- as.data.frame(lmm_summary$coefficients)
    
    interaction_term <- paste0(prs_var, ":Time_Years")
    
    if (interaction_term %in% rownames(coef_df)) {
      result <- data.frame(
        Age_Group = age_group,
        PRS = prs_var,
        Cognitive = cognitive_var,
        Interaction_Beta = coef_df[interaction_term, "Estimate"],
        Interaction_SE = coef_df[interaction_term, "Std. Error"],
        Interaction_P = 2 * pnorm(-abs(coef_df[interaction_term, "t value"])),
        N_Subjects = n_subjects,
        stringsAsFactors = FALSE
      )
      return(result)
    }
  }, error = function(e) {
    return(NULL)
  })
  
  return(NULL)
}

lmm_age_results <- list()

for (age_grp in age_groups) {
  for (prs in prs_vars) {
    for (cog in cognitive_vars) {
      if (all(c(prs, cog, "Age_Group") %in% colnames(long_data))) {
        result <- run_lmm_age_stratified(long_data, prs, cog, age_grp)
        if (!is.null(result)) {
          lmm_age_results[[paste(age_grp, prs, cog, sep = "_")]] <- result
        }
      }
    }
  }
}

lmm_age_df <- do.call(rbind, lmm_age_results)

if (!is.null(lmm_age_df) && nrow(lmm_age_df) > 0) {
  lmm_age_df$FDR <- p.adjust(lmm_age_df$Interaction_P, method = "BH")
  fwrite(lmm_age_df, "results/tables/ADNI_LMM_Age_Stratified.csv")
}

# ==============================================================================
# Part 10: Visualization
# ==============================================================================

## Section 10.1: PRS-Biomarker Heatmap
if (exists("prs_biomarker_df") && nrow(prs_biomarker_df) > 0) {
  
  # Create effect size matrix
  heatmap_data <- prs_biomarker_df %>%
    select(PRS, Biomarker, Beta) %>%
    pivot_wider(names_from = Biomarker, values_from = Beta) %>%
    column_to_rownames("PRS")
  
  # Create significance matrix
  sig_data <- prs_biomarker_df %>%
    mutate(Sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    )) %>%
    select(PRS, Biomarker, Sig) %>%
    pivot_wider(names_from = Biomarker, values_from = Sig) %>%
    column_to_rownames("PRS")
  
  if (ncol(heatmap_data) > 0 && nrow(heatmap_data) > 0) {
    pdf("results/figures/ADNI_PRS_Biomarker_Heatmap.pdf", width = 12, height = 10)
    pheatmap(
      as.matrix(heatmap_data),
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      display_numbers = as.matrix(sig_data),
      fontsize_number = 12,
      main = "PRS-Biomarker Associations (Beta)",
      fontsize = 10
    )
    dev.off()
  }
}

## Section 10.2: Cox Regression Forest Plot
if (exists("cox_df") && nrow(cox_df) > 0) {
  
  cox_plot_data <- cox_df %>%
    filter(!is.na(HR)) %>%
    arrange(HR)
  
  cox_plot_data$PRS <- factor(cox_plot_data$PRS, levels = cox_plot_data$PRS)
  
  p_forest <- ggplot(cox_plot_data, aes(x = HR, y = PRS)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = HR_Lower, xmax = HR_Upper), height = 0.2) +
    geom_point(aes(color = FDR < 0.05), size = 3) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "#B2182B"),
                       name = "FDR < 0.05") +
    labs(x = "Hazard Ratio (95% CI)", y = "",
         title = "PRS Association with MCI-to-Dementia Conversion") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave("results/figures/ADNI_Cox_Forest_Plot.pdf", p_forest, 
         width = 10, height = 8)
}

## Section 10.3: Age-Stratified Effect Comparison
if (exists("age_stratified_df") && nrow(age_stratified_df) > 0) {
  
  # Select key PRS for visualization
  key_prs <- c("PRS_EOAD_Global", "PRS_LOAD_Global", 
               "PRS_EOAD_TCell", "PRS_LOAD_Microglia")
  
  plot_data <- age_stratified_df %>%
    filter(PRS %in% key_prs, Biomarker %in% c("Hippo_Normalized", "WMH_Normalized"))
  
  if (nrow(plot_data) > 0) {
    p_age <- ggplot(plot_data, aes(x = PRS, y = Beta, fill = Age_Group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_errorbar(aes(ymin = Beta - 1.96 * SE, ymax = Beta + 1.96 * SE),
                    position = position_dodge(width = 0.8), width = 0.2) +
      facet_wrap(~Biomarker, scales = "free_y") +
      scale_fill_manual(values = c("Younger (<70)" = "#4DAF4A", "Older (≥70)" = "#984EA3")) +
      labs(x = "", y = "Effect Size (Beta)", 
           title = "Age-Stratified PRS Effects on Neuroimaging") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")
    
    ggsave("results/figures/ADNI_Age_Stratified_Effects.pdf", p_age, 
           width = 12, height = 6)
  }
}

## Section 10.4: Cluster Visualization
if (exists("cluster_data") && "Cluster" %in% colnames(cluster_data)) {
  
  # PCA for visualization
  pca_result <- prcomp(cluster_matrix, scale. = FALSE)
  
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Cluster = cluster_data$Cluster
  )
  
  p_cluster <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point(alpha = 0.7, size = 2) +
    stat_ellipse(level = 0.95) +
    scale_color_brewer(palette = "Set1") +
    labs(title = paste0("K-means Clustering (k=3, Silhouette = ", 
                        round(avg_silhouette, 3), ")"),
         x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
    theme_minimal()
  
  ggsave("results/figures/ADNI_Cluster_PCA.pdf", p_cluster, width = 8, height = 6)
}

## Section 10.5: Kaplan-Meier Survival Curves by PRS Tertiles
if (exists("survival_data") && "PRS_EOAD_Global" %in% colnames(survival_data)) {
  
  # Create PRS tertiles
  survival_data$PRS_Tertile <- cut(
    survival_data$PRS_EOAD_Global,
    breaks = quantile(survival_data$PRS_EOAD_Global, c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE
  )
  
  km_data <- survival_data[!is.na(survival_data$PRS_Tertile) & 
                             !is.na(survival_data$Time_Years) &
                             !is.na(survival_data$Event), ]
  
  if (nrow(km_data) > 50) {
    km_fit <- survfit(Surv(Time_Years, Event) ~ PRS_Tertile, data = km_data)
    
    p_km <- ggsurvplot(
      km_fit,
      data = km_data,
      pval = TRUE,
      risk.table = TRUE,
      palette = c("#2166AC", "#F7F7F7", "#B2182B"),
      title = "MCI-to-Dementia Conversion by EOAD PRS Tertile",
      xlab = "Time (Years)",
      ylab = "Conversion-Free Probability",
      legend.title = "PRS Tertile",
      legend.labs = c("Low", "Medium", "High")
    )
    
    pdf("results/figures/ADNI_KM_Survival.pdf", width = 10, height = 8)
    print(p_km)
    dev.off()
  }
}

# ==============================================================================
# Part 11: Results Summary
# ==============================================================================

## Section 11.1: Compile Summary Statistics
summary_stats <- list()

# Sample sizes
summary_stats$N_Total <- nrow(merged_data)
summary_stats$N_with_PRS <- sum(!is.na(merged_data$PRS_EOAD_Global))
summary_stats$N_with_CSF <- sum(!is.na(merged_data$ABETA))
summary_stats$N_with_FreeSurfer <- sum(!is.na(merged_data$Hippo_Normalized))

# Age distribution
if ("AGE" %in% colnames(merged_data)) {
  summary_stats$Age_Mean <- mean(merged_data$AGE, na.rm = TRUE)
  summary_stats$Age_SD <- sd(merged_data$AGE, na.rm = TRUE)
  summary_stats$N_Younger <- sum(merged_data$Age_Group == "Younger (<70)", na.rm = TRUE)
  summary_stats$N_Older <- sum(merged_data$Age_Group == "Older (≥70)", na.rm = TRUE)
}

# Significant associations
if (exists("prs_biomarker_df") && nrow(prs_biomarker_df) > 0) {
  summary_stats$N_Sig_PRS_Biomarker <- sum(prs_biomarker_df$FDR < 0.05)
}

if (exists("cox_df") && nrow(cox_df) > 0) {
  summary_stats$N_Sig_Cox <- sum(cox_df$FDR < 0.05)
}

if (exists("lmm_df") && nrow(lmm_df) > 0) {
  summary_stats$N_Sig_LMM <- sum(lmm_df$Interaction_FDR < 0.05)
}

## Section 11.2: Export Summary
summary_df <- data.frame(
  Metric = names(summary_stats),
  Value = unlist(summary_stats)
)

fwrite(summary_df, "results/tables/ADNI_Analysis_Summary.csv")

## Section 11.3: Key Findings Table
key_findings <- data.frame(
  Analysis = c("PRS-Biomarker", "Cox Regression", "LMM Trajectory", "Clustering"),
  Description = c(
    "PRS associations with CSF and neuroimaging biomarkers",
    "MCI-to-dementia conversion hazard ratios",
    "PRS effects on cognitive decline trajectories",
    "Unsupervised clustering based on PRS profiles"
  ),
  N_Tests = c(
    ifelse(exists("prs_biomarker_df"), nrow(prs_biomarker_df), 0),
    ifelse(exists("cox_df"), nrow(cox_df), 0),
    ifelse(exists("lmm_df"), nrow(lmm_df), 0),
    ifelse(exists("cluster_biomarker_df"), nrow(cluster_biomarker_df), 0)
  ),
  N_Significant = c(
    ifelse(exists("prs_biomarker_df"), sum(prs_biomarker_df$FDR < 0.05), 0),
    ifelse(exists("cox_df"), sum(cox_df$FDR < 0.05), 0),
    ifelse(exists("lmm_df"), sum(lmm_df$Interaction_FDR < 0.05), 0),
    ifelse(exists("cluster_biomarker_df"), sum(cluster_biomarker_df$FDR < 0.05), 0)
  )
)

fwrite(key_findings, "results/tables/ADNI_Key_Findings.csv")

cat("\n")
cat("==============================================================================\n")
cat("ADNI Validation Analysis Complete\n")
cat("==============================================================================\n")
cat("Results saved to: results/tables/ and results/figures/\n")
cat("Total subjects:", summary_stats$N_Total, "\n")
cat("Subjects with PRS:", summary_stats$N_with_PRS, "\n")
cat("Subjects with CSF:", summary_stats$N_with_CSF, "\n")
cat("Subjects with FreeSurfer:", summary_stats$N_with_FreeSurfer, "\n")
cat("==============================================================================\n")



