# Genetic Heterogeneity of Early-Onset versus Late-Onset Alzheimer's Disease: From Polygenic Architecture to Cell-Type-Specific Mechanisms

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains analysis code for a study of genetic heterogeneity between early-onset Alzheimer's disease (EOAD; onset <65 years) and late-onset Alzheimer's disease (LOAD). The analyses integrate GWAS summary statistics, gene-based association testing, spatial transcriptomic enrichment, transcriptome-wide association studies (TWAS), pathway-specific polygenic risk scores (PRS), ADNI genotype-based validation and phenotype-level contextual analyses in A4, HABS and AIBL.

The central result is that EOAD and LOAD show substantial global genetic overlap but diverge at the pathway level. LOAD signals are dominated by highly polygenic microglial and immune-related architecture, whereas EOAD shows a more sample-size-limited but internally robust lipid, oligodendrocyte and myelination-related component. The repository also includes robustness and sensitivity analyses added to address EOAD discovery power, lack of independent EOAD GWAS replication, and ADNI MRI/WMH acquisition heterogeneity.

## Repository Structure

The current repository uses top-level scripts. Two supplemental robustness scripts have descriptive file names without `.R` extensions because they were uploaded with the current manuscript-support naming convention.

```text
EOAD-LOAD-Genetics/
├── EOAD_LOAD_Pathway_Discovery.R
├── PRS_LDpred2_Analysis.R
├── ADNI_Validation.R
├── A4_HABS_AIBL_Validation.R
├── EOAD Internal Robustness Analyses
├── WMH Sensitivity Analyses
├── LICENSE
└── README.md
```

## Analysis Modules

### 1. Pathway Discovery

**Script:** `EOAD_LOAD_Pathway_Discovery.R`

Gene-based and pathway-level dissection of EOAD versus LOAD genetic architecture.

Main analyses:

- MAGMA gene-based association using precomputed `.genes.out` files.
- Graham et al. glial module enrichment using one-sided Fisher's exact tests.
- GO Biological Process gene set enrichment analysis with `clusterProfiler`.
- MAGMA-style gene-property and cell-type marker analyses.
- Conditional regression assessing whether T-cell signal is independent of microglial gene expression.
- LOAD downsampling validation matching LOAD effective sample size to EOAD to distinguish power-driven and architecture-driven enrichment.
- Publication-quality pathway and enrichment visualizations.

Primary outputs include gene-level association tables, module enrichment tables, GSEA outputs, conditional regression results, downsampling summaries and figures for pathway-level comparison.

### 2. PRS Construction and TWAS

**Script:** `PRS_LDpred2_Analysis.R`

Bayesian PRS construction, pathway-specific risk partitioning and S-PrediXcan TWAS.

Main analyses:

- S-PrediXcan TWAS using GTEx v8 brain prediction models across Cortex, Anterior cingulate cortex BA24, Frontal Cortex BA9, Hippocampus and Putamen basal ganglia.
- Wilcoxon rank-sum tests comparing oligodendrocyte-myelin genes with background genes.
- LOAD and multivariate aging controls to test EOAD specificity.
- LDpred2-auto with multiple MCMC chains and chain bagging.
- Pathway-specific PRS for T-cell activation, activated microglia, Abeta clearance, APP metabolism, oligodendrocyte and myelination pathways.
- Cellular-origin burden quantification.
- APOE-excluded PRS sensitivity by masking chr19:44-46 Mb.

Primary outputs include LDpred2 weights, pathway-specific PRS weights, APOE-masked PRS weights, cellular-origin burden tables and TWAS summary tables.

### 3. ADNI Validation

**Script:** `ADNI_Validation.R`

Primary genotype-based validation in ADNI.

Main analyses:

- Integration of ADNI pathway-specific PRS, baseline clinical variables, FreeSurfer MRI measures, CSF biomarkers and longitudinal diagnosis data.
- PRS-biomarker associations with adjustment for age, sex, education and APOE epsilon 4 allele count.
- APOE-stratified PRS-biomarker analyses.
- MCI-to-dementia conversion models.
- PRS-neuroimaging associations and age interaction models.
- Sliding-window age-dependency mapping.
- Unsupervised k-means genetic subtyping using six pathway-specific PRS.
- Hierarchical clustering concordance checks.
- Regional cortical and subcortical analyses.
- Demographic and nested complete-case summaries.

Primary outputs include ADNI association tables, conversion models, age-stratified MRI results, sliding-window summaries, clustering outputs, subtype summaries and regional brain maps.

### 4. External Cohort Contextual Analyses

**Script:** `A4_HABS_AIBL_Validation.R`

Phenotype-level contextual analyses in A4, HABS and AIBL. These cohorts do not provide the same genotype-based pathway PRS framework used in ADNI, so the analyses are interpreted as phenotype-level context rather than genetic replication.

Main analyses:

- A4: WMH-cognition association using heteroscedasticity-consistent standard errors, adjusted for age, sex, education, APOE epsilon 4 and amyloid burden.
- HABS: WMH, plasma p-tau217 and cognition mediation analyses, including bootstrap mediation and SEM confirmation.
- HABS: age-stratified mediation and formal interaction testing.
- AIBL: longitudinal conversion analyses using Cox proportional hazards and clinical utility models.
- Participant-flow and missing-data assessments.
- Cross-cohort forest plot visualization.

Primary outputs include A4 WMH-cognition results, HABS mediation and SEM tables, AIBL survival results, missing-data assessments, clinical utility metrics and cross-cohort figures.

### 5. EOAD Internal Robustness

**Script:** `EOAD Internal Robustness Analyses`

Internal robustness analyses addressing the limited EOAD GWAS discovery sample size and absence of a directly downloadable independent EOAD GWAS replication resource.

Main analyses:

- EOAD power and minimum detectable effect calculations using EOAD effective sample size.
- Leave-one-gene robustness analyses for the EOAD oligodendrocyte-myelin S-PrediXcan signal.
- APOE locus exclusion sensitivity for chr19:44-46 Mb.
- MHC locus exclusion sensitivity for chr6:25-34 Mb.
- APOE-or-MHC combined exclusion sensitivity.
- Matched random gene-set permutations using 2,000 matched random gene sets per tissue.
- Matching based on TWAS model SNP count, predicted-expression variance and MAGMA gene-level SNP-count bins, with relaxed matching when strata are sparse.

Outputs correspond to Supplementary Tables 40-42:

- `S40_EOAD_power_minimum_detectable_effect.csv`
- `S41_leave_one_gene_summary.csv`
- `S41_leave_one_gene_detail.csv`
- `S41_locus_exclusion_sensitivity.csv`
- `S42_matched_random_gene_set_permutation.csv`
- `S42_permutation_diagnostics_first1000.csv`

### 6. ADNI WMH Scanner-Aware Sensitivity

**Script:** `WMH Sensitivity Analyses`

Scanner-aware MRI and WMH sensitivity analyses addressing residual acquisition heterogeneity in ADNI.

Main analyses:

- Exact matching of analytic FreeSurfer rows back to ADNI imaging metadata using RID, intracranial volume and WMH volume.
- Exact matching of analytic UCD FLAIR WMH rows back to UCD imaging metadata using RID and total WMH volume.
- Batch summaries for ADNI phase, field strength, visit code, manufacturer and scanner model.
- FreeSurfer MRI sensitivity models with field-strength adjustment and 1.5T-only restriction.
- Age-stratified FreeSurfer sensitivity models.
- UCD FLAIR WMH sensitivity models with manufacturer adjustment, scanner-model adjustment, manufacturer-effect residualization and manufacturer-specific restrictions.

Outputs correspond to Supplementary Table 39:

- `ADNI_FreeSurfer_batch_enriched.csv`
- `ADNI_UCD_WMH_batch_enriched.csv`
- `ADNI_scanner_batch_counts.csv`
- `ADNI_scanner_aware_sensitivity_results.csv`
- `ADNI_scanner_sensitivity_summary.txt`

## Requirements

### R Version

- R >= 4.2.0

### R Packages

```r
install.packages(c(
  "data.table", "dplyr", "tidyr", "ggplot2", "cowplot",
  "broom", "survival", "survminer", "lme4", "lmerTest",
  "cluster", "factoextra", "pheatmap", "corrplot",
  "pROC", "lavaan", "logistf", "PRROC", "sandwich",
  "lmtest", "mediation", "meta", "ggseg", "flexsurv",
  "dcurves", "patchwork", "viridis", "RColorBrewer",
  "ggrepel", "scales", "ComplexHeatmap", "circlize"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "bigsnpr", "clusterProfiler", "org.Hs.eg.db",
  "org.Mm.eg.db", "biomaRt",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "GenomicRanges", "enrichplot"
))
```

## Example Usage

The scripts are written as function-based analysis modules. Replace the example paths with local paths after obtaining the required controlled-access or public datasets.

### Pathway discovery

```r
source("EOAD_LOAD_Pathway_Discovery.R")

# Example only. Modify file paths and arguments to match local data.
# results <- run_pathway_discovery(
#   eoad_magma_file = "data/MAGMA/EOAD.genes.out",
#   load_magma_file = "data/MAGMA/LOAD.genes.out",
#   output_dir = "results/pathway_discovery"
# )
```

### PRS and TWAS

```r
source("PRS_LDpred2_Analysis.R")

# Example only. Modify file paths and arguments to match local data.
# results <- run_prs_analysis(
#   eoad_gwas_file  = "data/EOAD_GWAS_hg19.txt",
#   load_gwas_file  = "data/LOAD_GWAS_hg19.txt",
#   aging_gwas_file = "data/Aging_GWAS_hg19.txt",
#   ldref_dir       = "data/ldref_hm3_plus",
#   output_dir      = "results/prs",
#   eoad_n_case     = 1573,
#   eoad_n_control  = 199505,
#   load_n_case     = 85934,
#   load_n_control  = 401577
# )
```

### ADNI validation

```r
source("ADNI_Validation.R")

# Example only. Modify file paths and arguments to match local data.
# results <- run_adni_validation(
#   clinical_file = "data/ADNIMERGE.csv",
#   prs_file      = "data/ADNI_pathway_PRS.csv",
#   mri_file      = "data/UCSFFSX51_ADNI1GO2.csv",
#   csf_file      = "data/UPENNBIOMK_MASTER.csv",
#   strem2_file   = "data/ADNI_CSF_sTREM2.csv",
#   output_dir    = "results/adni_validation"
# )
```

### External cohort contextual analyses

```r
source("A4_HABS_AIBL_Validation.R")

# Example only. Modify file paths and arguments to match local data.
# results <- run_cross_cohort_validation(
#   a4_subjinfo   = "data/A4/SUBJINFO.csv",
#   a4_ptdemog    = "data/A4/PTDEMOG.csv",
#   a4_sppacc     = "data/A4/SPPACC.csv",
#   a4_petsuvr    = "data/A4/PETSUVR.csv",
#   a4_vmri       = "data/A4/VMRI.csv",
#   habs_genomics_xlsx  = "data/HABS/Genomics.xlsx",
#   habs_biomarkers_csv = "data/HABS/Biomarkers.csv",
#   habs_wmh_xlsx       = "data/HABS/WMH.xlsx",
#   habs_clinical_csv   = "data/HABS/Clinical.csv",
#   aibl_ptdemog  = "data/AIBL/ptdemog.csv",
#   aibl_apoeres  = "data/AIBL/apoeres.csv",
#   aibl_pdxconv  = "data/AIBL/pdxconv.csv",
#   output_dir    = "results/external_cohorts"
# )
```

### EOAD internal robustness

Because this file name contains spaces and has no `.R` extension, use the exact quoted file name.

```r
source("EOAD Internal Robustness Analyses")

twas_files <- c(
  "Cortex" =
    "data/SPrediXcan/EOAD_hg19__PM__mashr_Brain_Cortex.csv",
  "Anterior cingulate cortex BA24" =
    "data/SPrediXcan/EOAD_hg19__PM__mashr_Brain_Anterior_cingulate_cortex_BA24.csv",
  "Putamen basal ganglia" =
    "data/SPrediXcan/EOAD_hg19__PM__mashr_Brain_Putamen_basal_ganglia.csv",
  "Frontal Cortex BA9" =
    "data/SPrediXcan/EOAD_hg19__PM__mashr_Brain_Frontal_Cortex_BA9.csv",
  "Hippocampus" =
    "data/SPrediXcan/EOAD_hg19__PM__mashr_Brain_Hippocampus.csv"
)

# Example only. Modify file paths and arguments to match local data.
# robustness_results <- run_eoad_internal_robustness(
#   magma_gene_file = "data/MAGMA/Supplementary_Table_S1_EOAD_MAGMA_Complete.csv",
#   oligo_gene_file = "data/SPrediXcan/Supplementary_Table_2_Oligodendrocyte_Genes.csv",
#   twas_files = twas_files,
#   output_dir = "results/eoad_internal_robustness",
#   eoad_effective_n = 1573,
#   n_perm = 2000,
#   random_seed = 20260525
# )
```

### ADNI WMH scanner-aware sensitivity

Because this file name contains spaces and has no `.R` extension, use the exact quoted file name.

```r
source("WMH Sensitivity Analyses")

# Example only. Modify file paths and arguments to match local data.
# scanner_results <- run_adni_scanner_aware_sensitivity(
#   final_analytic_file =
#     "data/ADNI_Age_Stratified_Results/ADNI_Age_Stratified_Final.csv",
#   stratified_analytic_file =
#     "data/ADNI_Validation_Results_Stratified/ADNI_Merged_Data_Stratified.csv",
#   freesurfer_metadata_file =
#     "data/ADNI/data/MRI/MRI_FreeSurfer.csv",
#   ucd_wmh_metadata_file =
#     "data/ADNI/data/MRI/UCD_WMH_03Jan2026.csv",
#   output_dir = "results/adni_scanner_sensitivity"
# )
```

## Data Availability

### GWAS Summary Statistics

| Dataset | Source | Sample Size | Population | Access |
|---|---|---:|---|---|
| EOAD | FinnGen R11 AD_EO_EXMORE | 1,573 cases; 199,505 controls | Finnish / European | FinnGen public release R11 |
| LOAD | EADB consortium | 85,934 cases; 401,577 controls | European | EBI GWAS Catalog / EADB |
| Aging | multivariate aging GWAS | trait-specific sample sizes | European | Edinburgh DataShare and associated public repositories |

### Spatial and Transcriptomic Resources

The spatial transcriptomic and TWAS analyses use public or controlled resources including GTEx v8 prediction models, LIBD spatial transcriptomic data, Maynard et al. DLPFC spatial data and MOSTA mouse spatial transcriptomic data. Users should download these resources from the original data providers and preserve the file names used in the example paths or update the script arguments accordingly.

### Validation Cohorts

ADNI, A4, HABS and AIBL data require approval from their respective data access committees or study investigators. This repository does not redistribute individual-level clinical, imaging, genotype or biomarker data.

## Reproducibility Notes

1. Scripts assume that GWAS summary statistics have been harmonized to the column names expected by each module or supplied through the argument names described in the example usage.
2. ADNI, A4, HABS and AIBL individual-level data require approved data-use agreements and are not included in this repository.
3. File paths in all examples are placeholders and should be modified to match the local directory structure.
4. LDpred2 and genotype-based PRS construction require sufficient memory for LD reference matrices; at least 32 GB RAM is recommended.
5. MAGMA analyses require precomputed gene annotation files and an appropriate European LD reference panel.
6. S-PrediXcan analyses require GTEx v8 prediction model `.db` files and covariance files.
7. The EOAD internal robustness script uses 2,000 matched random gene-set permutations by default. Increase `n_perm` for more precise empirical P values if runtime allows.
8. The scanner-aware sensitivity script is not a full cross-cohort ComBat harmonization. It provides acquisition-aware sensitivity checks using the scanner metadata available for the ADNI MRI/WMH subsets.
9. The files `EOAD Internal Robustness Analyses` and `WMH Sensitivity Analyses` are valid R source files despite lacking `.R` extensions. Use their exact names in `source()` calls.

## External Tools

| Tool | Purpose |
|---|---|
| LDSC | Cross-trait genetic correlation and SNP heritability |
| HESS | Local SNP heritability partitioning |
| MAGMA | Gene-based association and gene-property analysis |
| LAVA | Local genetic correlation |
| conjFDR | Conditional and conjunctional FDR analyses |
| FUMA | Functional mapping and annotation |
| gsMap | Spatial transcriptomic enrichment mapping |
| S-PrediXcan / MetaXcan | Transcriptome-wide association analysis |
| LDpred2 | Bayesian PRS construction |
| g:Profiler | Functional enrichment of gene sets |

## Citation

If you use this repository, please cite the associated manuscript and the original data resources, software packages and cohorts listed in the manuscript.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contact

For questions about the code, please open an issue in this repository or contact the corresponding author.

## Acknowledgments

We thank the FinnGen Consortium, the European Alzheimer & Dementia Biobank Consortium, UK Biobank, ADNI, A4, HABS, AIBL, GTEx and the developers of the external tools used in this work. Data collection and sharing for ADNI was funded by the Alzheimer's Disease Neuroimaging Initiative and related supporting organizations. The A4 Study, HABS and AIBL are acknowledged according to their required data-use terms.
