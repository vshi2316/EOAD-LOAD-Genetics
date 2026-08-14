[Uploading README.md…]()
# EOAD and later-onset Alzheimer disease genetics

This repository contains the analysis code for genetic architecture, independent pathway replication, participant-level biomarker associations, multimodal pathology and structural analyses, and clinical expression across Alzheimer disease settings.

## Study workflow

The analysis begins with gene and pathway comparisons between early-onset Alzheimer disease (EOAD) and the European Alzheimer and Dementia Biobank (EADB) Alzheimer disease dataset. Two amyloid-processing pathways are tested in an independent age-at-onset genome-wide association study. Variant-level replication and an inverse-variance weighted analysis refine the replicated signal. White-matter traits, transcriptome-wide association results, and single-nucleus RNA sequencing characterize a parallel oligodendrocyte and myelin context. Participant-level analyses then connect an amyloid-processing polygenic score with continuous Centiloid burden in Alzheimer's Disease Neuroimaging Initiative (ADNI) and evaluate amyloid, tau, temporal cortical structure, and cognition in Health and Aging Brain Study: Health Disparities (HABS-HD). A4 and ADNI place pathology-associated cognitive expression across preclinical and clinical settings.

## Repository contents

| Script | Analysis | Principal manuscript output |
|---|---|---|
| `00_preflight.R` | Package, executable, and input checks | Reproducibility audit |
| `01_genetic_architecture.R` | EOAD and EADB gene, pathway, enrichment, and overlap analyses | Figure 1 and supporting tables |
| `02_polygenic_and_transcriptomic_analysis.R` | LDpred2, pathway scores, tissue models, and APOE-region sensitivity | Figure 1, Figure 3, and supporting tables |
| `03_independent_pathway_replication.R` | Prespecified pathway tests in Bradley age-at-onset datasets | Figure 2 and replication tables |
| `04_variant_replication.R` | Variant and locus comparison, including rs56368748 | Figure 2 and variant tables |
| `05_meta_enhanced_analysis.R` | FinnGen and Bradley inverse-variance weighted analysis followed by MAGMA | Figure 2 and meta-analysis tables |
| `06_genetic_robustness.R` | Gene-set restriction, permutation, colocalization, and sensitivity analyses | Supplementary results |
| `07_white_matter_analysis.R` | White-matter genetic associations, functional annotation, and imaging metadata utilities | Figure 3 and supporting tables |
| `08_single_nucleus_analysis.R` | GSE272082 cell-type annotation, pseudobulk testing, and pathway activity | Figure 3 and supplementary results |
| `09_adni_pathway_prs_analysis.R` | Entry point for ADNI pathway scoring and Centiloid models | Figure 4 |
| `09_1_adni_pathway_prs_data.R` | PLINK pathway scoring and phenotype merge | ADNI analysis dataset |
| `09_2_adni_pathway_prs_centiloid.R` | Continuous Centiloid models and center sensitivity | Figure 4 |
| `10_habs_multimodal_analysis.R` | Entry point for the HABS-HD modules | Figure 5 and supporting tables |
| `10_1_habs_multimodal_data.R` | Date matching and temporal cortical score construction | HABS-HD analysis dataset |
| `10_2_habs_longitudinal_cognition.R` | Cognitive orientation and longitudinal matching | Longitudinal analysis dataset |
| `10_3_habs_cognition_models.R` | Concurrent and longitudinal cognitive models | Secondary HABS-HD results |
| `10_4_habs_primary_models.R` | Three primary HC3 models, multiplicity control, and sensitivity analyses | Figure 5 and primary tables |
| `10_5_habs_regional_analysis.R` | Regional tau associations and temporal cortical surface display | Figure 5 source data |
| `11_clinical_expression.R` | A4 and ADNI pathology-to-cognition models with HABS-HD structural context | Figure 6 |
| `12_clinical_figures.R` | Vector and high-resolution raster exports for Figures 4 and 6 | Figure 4 and Figure 6 |
| `run_all.R` | Ordered execution of the analysis modules | Complete result tree |

## Software

The analyses use R 4.4 or later. MAGMA is required for gene and gene-set analysis. PLINK or PLINK 2 is required when pathway scores are calculated from genotype dosage files outside the repository.

CRAN packages include `data.table`, `dplyr`, `tidyr`, `readr`, `stringr`, `tibble`, `purrr`, `ggplot2`, `ggrepel`, `scales`, `forcats`, `broom`, `lmtest`, `sandwich`, `cowplot`, `patchwork`, `RColorBrewer`, `readxl`, `ggseg`, `bigsnpr`, `bigstatsr`, `svglite`, and `ragg`.

Bioconductor packages include `GenomicRanges`, `IRanges`, `rtracklayer`, `BSgenome`, `SNPlocs.Hsapiens.dbSNP155.GRCh37`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `ComplexHeatmap`, `edgeR`, `Seurat`, `SeuratObject`, and the genome annotation packages named in the scripts.

## Data access

Summary statistics and controlled cohort data retain their source licenses and access conditions. Participant-level data are excluded from this repository.

| Dataset | Access route | Use in this repository |
|---|---|---|
| FinnGen | FinnGen public results service | EOAD discovery and meta-analysis |
| EADB | Access conditions of the contributing study | Alzheimer disease comparison |
| Bradley age-at-onset GWAS | Study repository or author-designated archive | Independent pathway and variant replication |
| GSE272082 | NCBI Gene Expression Omnibus | Single-nucleus analysis |
| ADNI | LONI Image and Data Archive, subject to ADNI approval | Pathway PRS, Centiloid, and cognitive expression |
| A4 | LONI Image and Data Archive, subject to A4 approval | Preclinical Centiloid and PACC analysis |
| HABS-HD | LONI Image and Data Archive, subject to HABS-HD approval | FBB PET, PI-2620 PET, cortical thickness, cognition, and covariates |

MAGMA reference files and linkage disequilibrium reference panels are obtained under their respective distribution terms.

## Directory layout

The scripts use repository-relative directories by default.

```text
EOAD-LOAD-Genetics/
  data/
    A4/
    ADNI/
    Bradley/
    EADB/
    FinnGen/
    GSE272082/
    HABS_HD/
    reference/
  results/
  config.R
  run_all.R
```

`data/`, `results/`, and `config.R` are excluded by `.gitignore`.

## Configuration

Create `config.R` from `config_template.R`. Each setting may point to an authorized storage location outside the repository. Environment variables take precedence over default relative paths.

Core project variables:

```r
Sys.setenv(
  EOAD_PROJECT_DIR = getwd(),
  EOAD_DATA_DIR = file.path(getwd(), "data"),
  EOAD_RESULTS_DIR = file.path(getwd(), "results")
)
```

Genetic analyses also require the variables declared near the start of scripts `03`, `04`, and `05`, including `FINNGEN_EOAD_SUMSTATS`, `BRADLEY_RAW_DIR`, `MAGMA_EXE`, `MAGMA_BFILE`, `MAGMA_GENE_LOC`, and `MAGMA_SET_ANNOT`.

## Analysis-ready input contracts

### ADNI pathway PRS and Centiloid

`ADNI_PATHWAY_ANALYSIS_FILE` points to a tab-separated file with one row per participant and the following fields:

```text
RID
Centiloid
PRS_AmyloidBetaClearance_z
PRS_NegativeRegulationAPPCatabolism_z
age_z
sex
education_z
apoe4_z
baseline_dx
site
PC1_z ... PC10_z
```

The score columns contain standardized pathway scores constructed from the reported genetic weights after allele harmonization and genotype quality control. `09_1_adni_pathway_prs_data.R` calculates both scores. `09_2_adni_pathway_prs_centiloid.R` fits the primary joint model, center-clustered inference, center fixed-effect sensitivity, binary amyloid thresholds, winsorized Centiloid, ancestry restriction, and leave-one-center-out analyses.

### A4 clinical expression

`A4_CLINICAL_ANALYSIS_FILE` contains:

```text
subject_id
pathology
cognition
age
sex
education
apoe4
```

`pathology` is Centiloid. `cognition` is oriented so that higher values indicate greater cognitive impairment. For PACC, multiply the original score by minus one before export.

### ADNI clinical expression

`ADNI_CLINICAL_ANALYSIS_FILE` contains:

```text
subject_id
pathology
cognition
age
sex
education
apoe4
diagnosis
```

`pathology` is Centiloid and `cognition` is ADAS13. Each participant contributes one prespecified visit. The diagnosis field supports clinical-stage sensitivity analyses.

### HABS-HD files

The HABS-HD variables in `config_template.R` point to the authorized spreadsheets or comma-separated files for FBB PET, PI-2620 PET, cortical thickness, genomics, diagnosis, demographics, cognitive z scores, MMSE, and education.

The cortical thickness workbook must contain bilateral entorhinal, parahippocampal, inferior temporal, and middle temporal measures. Month-level visit dates are assigned to day 15. The primary matching window is 180 days, with 90-day and 365-day sensitivity analyses.

## Temporal structural score

The score uses eight cortical thickness measures:

```text
left and right entorhinal cortex
left and right parahippocampal cortex
left and right inferior temporal cortex
left and right middle temporal cortex
```

Participants require all eight measurements. Five deterministic folds are assigned before model fitting. For each held-out fold, every measure is standardized with the mean and standard deviation calculated from the other four folds. The score is the negative, equally weighted mean of the eight standardized values. Higher scores indicate thinner temporal cortex. Score construction uses the cortical thickness measurements alone.

## Statistical models

The independent pathway replication tests the five prespecified pathways carried from the FinnGen analysis. Benjamini-Hochberg control is applied across this family.

The ADNI genetic bridge uses continuous Centiloid as the primary outcome. The joint model contains both replicated pathway scores with age, sex, education, APOE dosage, baseline diagnosis, and ten ancestry principal components. Research center defines the cluster for robust inference.

The HABS-HD primary family contains three associations:

```text
FBB Centiloid to PI-2620 mesial-temporal tau
PI-2620 tau to the temporal structural score
temporal structural score to concurrent global cognition
```

Models use HC3 standard errors. Age, sex, education, ethnicity, APOE4, MRI scanner, and the applicable upstream pathology measurements are included according to the model definitions in `10_4_habs_primary_models.R`. Benjamini-Hochberg correction is applied jointly to the three primary coefficients. Regional tau associations form a separate family. Scanner-specific estimates require at least 100 complete observations.

Longitudinal models evaluate whether structural change predicts later cognition. Follow-up MRI occurs 365 to 1,825 days after baseline. Cognitive follow-up occurs 90 to 730 days after the selected MRI. These analyses report the observed null boundary alongside the cross-sectional results.

## Running the analyses

Run the complete workflow from the repository root:

```r
source("config.R", encoding = "UTF-8")
source("run_all.R", encoding = "UTF-8")
```

Individual modules can be run after their prerequisite outputs are present:

```r
source("03_independent_pathway_replication.R", encoding = "UTF-8")
source("09_adni_pathway_prs_analysis.R", encoding = "UTF-8")
source("10_habs_multimodal_analysis.R", encoding = "UTF-8")
source("11_clinical_expression.R", encoding = "UTF-8")
source("12_clinical_figures.R", encoding = "UTF-8")
```

Every analysis directory records model tables, sample counts, audit tables, and `session_info.txt`. Figure exports use PDF for vector output and TIFF at 600 dots per inch where raster output is requested.

## Result checks

Successful reproduction should recover the reported result pattern:

| Analysis | Expected pattern |
|---|---|
| Independent pathway replication | Two of five prespecified pathways pass false discovery rate control |
| Variant replication | rs56368748 has a directionally concordant association |
| ADNI APP-related pathway score | Positive association with continuous Centiloid and stable direction across center omissions |
| HABS-HD amyloid and tau | Positive association |
| HABS-HD tau and temporal structure | Positive association with the vulnerability-oriented score |
| HABS-HD temporal structure and cognition | Greater structural vulnerability accompanies poorer concurrent cognition |
| HABS-HD longitudinal structure and cognition | No supported association in the available longitudinal sample |

Exact coefficients and sample sizes are written by the scripts and can be checked against the article source-data tables.

## Reproducibility records

Analysis tables include the fitted formula, focal coefficient, standard error, confidence interval, P value, adjusted P value where applicable, and complete-case sample size. Date-matching audits report eligible pairs under each window. Scanner audits report complete observations before scanner-stratified models are fitted. Random procedures use fixed seeds declared in the corresponding script.

## License

Code is distributed under the MIT License. Source datasets remain governed by their original terms.

## Citation

Please cite the associated article when using this workflow. Dataset-specific acknowledgments and citations should follow the requirements of FinnGen, EADB, Bradley and colleagues, ADNI, A4, HABS-HD, and the Gene Expression Omnibus record for GSE272082.
