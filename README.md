[Uploading README.md…]()
# EOAD and later-onset Alzheimer disease genetics

This repository contains the analysis workflow for shared genetic architecture, independent pathway replication, variant-level follow-up, white-matter and oligodendrocyte context, participant-level amyloid associations, HABS-HD multimodal pathology models, age modification, and clinical expression across complementary Alzheimer disease settings.

## Study workflow

The workflow follows the order used for the manuscript evidence. Genome-wide architecture is characterized in FinnGen early-onset Alzheimer disease (EOAD) and the broader European Alzheimer and Dementia Biobank (EADB) comparator. Five prespecified amyloid-processing pathways are tested in an independent age-at-onset genome-wide association study. Variant-level replication and inverse-variance weighted meta-analysis refine the replicated signal. White-matter traits, tissue-level transcriptome-wide association study (TWAS) results, and single-nucleus RNA sequencing provide a parallel oligodendrocyte and myelin context. Participant-level analyses then test pathway polygenic scores against continuous amyloid burden in the Alzheimer Disease Neuroimaging Initiative (ADNI). HABS-HD analyses quantify amyloid, tau, temporal cortical structure, concurrent cognition, longitudinal boundary results, and prespecified age modification. A4 and ADNI provide clinical-expression context across preclinical and clinical settings.

## Repository contents

| Script | Analysis | Principal output |
|---|---|---|
| `00_preflight.R` | Package, executable, input, and directory checks | Reproducibility audit |
| `01_genetic_architecture.R` | Genome-wide correlation, regional architecture, gene, pathway, and cell-marker analyses | Figure 1 and supporting tables |
| `02_polygenic_and_transcriptomic_analysis.R` | LDpred2 models, pathway scores, tissue analyses, and APOE-region sensitivity | Figures 1 and 3 supporting results |
| `03_independent_pathway_replication.R` | Prespecified pathway tests in independent age-at-onset GWAS data | Figure 2 and replication tables |
| `04_variant_replication.R` | Variant harmonization, rs56368748 follow-up, and locus comparison | Figure 2 and variant tables |
| `05_meta_enhanced_analysis.R` | FinnGen and Bradley inverse-variance weighted meta-analysis followed by MAGMA | Figure 2 and meta-analysis tables |
| `06_genetic_robustness.R` | Gene-set restriction, matched permutation, colocalization, SMR, and boundary analyses | Supplementary results |
| `07_white_matter_analysis.R` | White-matter conjunctional FDR, precision attenuation, and functional annotation | Figure 3 and supporting tables |
| `08_single_nucleus_analysis.R` | GSE272082 cell annotation, donor-level pseudobulk models, and pathway activity | Figure 3 and supplementary results |
| `09_adni_pathway_prs_analysis.R` | Entry point for ADNI pathway scoring and Centiloid analysis | Figure 4 inputs |
| `09_1_adni_pathway_prs_data.R` | Genotype scoring, allele harmonization, and phenotype merge | ADNI analysis dataset |
| `09_2_adni_pathway_prs_centiloid.R` | Continuous Centiloid models, robust inference, and center sensitivity | Figure 4 source tables |
| `10_habs_multimodal_analysis.R` | Ordered entry point for the HABS-HD modules | Figure 5 inputs and HABS-HD results |
| `10_1_habs_multimodal_data.R` | MRI-anchored date matching and temporal structural score construction | HABS-HD analysis dataset |
| `10_2_habs_longitudinal_cognition.R` | Cognitive orientation, baseline matching, and longitudinal matching | Longitudinal analysis dataset |
| `10_3_habs_cognition_models.R` | Concurrent and longitudinal cognitive models | Secondary HABS-HD results |
| `10_4_habs_primary_models.R` | Three prespecified HC3 models, primary-family correction, and sensitivity analyses | Figure 5 and primary tables |
| `10_5_habs_regional_analysis.R` | Regional tau associations and temporal structural summaries | Figure 5 source data |
| `10_6_habs_age65_multimodal_interactions.R` | Prespecified age-group and continuous-age interactions for the three HABS-HD primary associations | Age-interaction tables and audit files |
| `11_clinical_expression.R` | A4 and ADNI pathology-to-cognition models with HABS-HD context | Figure 6 source tables |
| `13_manuscript_figures.R` | Submission figure assembly for Figures 3 to 6, including HABS-HD age panels | Manuscript Figures 3 to 6 |
| `run_all.R` | Ordered execution of the main workflow | Complete result tree |

The default workflow includes the prespecified age-interaction analysis through `10_6_habs_age65_multimodal_interactions.R`.

## Software

The workflow uses R 4.4 or later. MAGMA is required for gene and gene-set analysis. PLINK or PLINK 2 is required when pathway scores are calculated from genotype dosage files outside the repository.

CRAN packages include `data.table`, `dplyr`, `tidyr`, `readr`, `stringr`, `tibble`, `purrr`, `ggplot2`, `ggrepel`, `scales`, `forcats`, `broom`, `lmtest`, `sandwich`, `cowplot`, `patchwork`, `RColorBrewer`, `readxl`, `ggseg`, `bigsnpr`, `bigstatsr`, `svglite`, and `ragg`.

Bioconductor packages include `GenomicRanges`, `IRanges`, `rtracklayer`, `BSgenome`, `SNPlocs.Hsapiens.dbSNP155.GRCh37`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `ComplexHeatmap`, `edgeR`, `Seurat`, `SeuratObject`, and the annotation packages required by the installed scripts.

## Data access

Summary statistics and controlled cohort data retain their source licenses and access conditions. Participant-level data remain outside this repository.

| Dataset | Access route | Repository use |
|---|---|---|
| FinnGen | FinnGen public results service | EOAD discovery and meta-analysis |
| EADB | Access conditions of the contributing study | Alzheimer disease comparison |
| Bradley age-at-onset GWAS | Study repository or author-designated archive | Independent pathway and variant replication |
| GSE272082 | NCBI Gene Expression Omnibus | Single-nucleus analysis |
| ADNI | LONI Image and Data Archive, subject to ADNI approval | Pathway PRS, Centiloid, and clinical-expression analyses |
| A4 | LONI Image and Data Archive, subject to A4 approval | Preclinical Centiloid and PACC analysis |
| HABS-HD | LONI Image and Data Archive, subject to HABS-HD approval | Florbetaben PET, PI-2620 PET, cortical thickness, cognition, and covariates |

MAGMA reference files and linkage disequilibrium reference panels are obtained under their respective distribution terms.

## Directory layout

The scripts use repository-relative directories when the required environment variables are absent.

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

Controlled data, large reference files, and generated results are excluded by `.gitignore`.

## Configuration

Create `config.R` from `config_template.R`. Environment variables take precedence over repository-relative defaults and may point to authorized storage locations outside the repository.

```r
Sys.setenv(
  EOAD_PROJECT_DIR = getwd(),
  EOAD_DATA_DIR = file.path(getwd(), "data"),
  EOAD_RESULTS_DIR = file.path(getwd(), "results")
)
```

Genetic analyses also require the variables declared near the start of scripts `03`, `04`, and `05`, including `FINNGEN_EOAD_SUMSTATS`, `BRADLEY_RAW_DIR`, `MAGMA_EXE`, `MAGMA_BFILE`, `MAGMA_GENE_LOC`, and `MAGMA_SET_ANNOT`.

The figure assembly script accepts `EOAD_ANALYSIS_ROOT`, `EOAD_RESULTS_DIR`, `SUPPLEMENTARY_TABLE_FILE`, and `MANUSCRIPT_FIGURE_DIR`. When these variables are defined, all result and table inputs are resolved from those locations. When they are absent, the script searches configured result roots and available project-level directories for the required analysis markers.

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

The score columns contain standardized pathway scores constructed from the reported genetic weights after allele harmonization and genotype quality control. `09_1_adni_pathway_prs_data.R` calculates both scores. `09_2_adni_pathway_prs_centiloid.R` fits the primary joint model, center-robust inference, center fixed-effect sensitivity, binary amyloid thresholds, winsorized and trimmed Centiloid models, ancestry restriction, and leave-one-center-out analyses.

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

`pathology` is Centiloid and `cognition` is ADAS13. Each participant contributes the quality-controlled baseline record marked time = 0 in the longitudinal Centiloid export. The export contains one record per RID, representing the earliest eligible amyloid acquisition after quality control. Later acquisitions remain in the longitudinal files, and unresolved within-RID date ties are excluded from the baseline export.

### HABS-HD files

The HABS-HD variables in `config_template.R` point to authorized files for florbetaben PET, PI-2620 PET, cortical thickness, diagnosis, demographics, cognitive z scores, MMSE, and education.

The cortical thickness input contains bilateral entorhinal, parahippocampal, inferior temporal, and middle temporal measures. Month-level visit dates are assigned to day 15. The primary matching window is 180 days, with 90-day and 365-day direct PET-interval sensitivity analyses within the MRI-anchored participant set.

## Temporal structural score

The score uses eight cortical thickness measures:

```text
left and right entorhinal cortex
left and right parahippocampal cortex
left and right inferior temporal cortex
left and right middle temporal cortex
```

Participants require all eight measurements. Five deterministic folds are assigned before model fitting. For each held-out fold, every measure is standardized with the mean and standard deviation calculated from the other four folds. The score is the negative, equally weighted mean of the eight standardized values. Higher scores indicate thinner temporal cortex. Score construction uses cortical thickness measurements alone.

## HABS-HD age modification

The age analysis uses the frozen HABS-HD baseline dataset and the frozen temporal structural score. Age at the MRI anchor defines two descriptive groups: younger than 65 years and 65 years or older. The grouping describes age at assessment and does not assign age at disease onset.

Each of the three primary models adds an exposure by age-group interaction while retaining standardized continuous age and the original covariates. The three interaction terms form one Benjamini-Hochberg family. Continuous-age interactions form a separate three-test family. Age-stratified coefficients describe the direction and magnitude within each group. Diagnosis-adjusted models and analyses restricted to Vida1 scans provide prespecified sensitivity checks.

The three primary interaction models are:

```text
PI-2620 mesial-temporal tau ~ florbetaben Centiloid * age group + covariates
temporal structural score ~ PI-2620 mesial-temporal tau * age group + covariates
global cognition ~ temporal structural score * age group + covariates
```

The age-defined subgroup is interpreted as an age-modified clinical-expression analysis within HABS-HD. It is not treated as an EOAD case series or as an independent EOAD replication cohort.

## Statistical models

Independent pathway replication tests the five prespecified pathways carried from the FinnGen analysis. Benjamini-Hochberg control is applied across this family.

The ADNI genetic bridge uses continuous Centiloid as the primary outcome. The joint model contains both replicated pathway scores with age, sex, education, APOE dosage, baseline diagnosis, and ten ancestry principal components. Research center defines the cluster for robust inference.

The HABS-HD primary family contains three associations:

```text
florbetaben Centiloid to PI-2620 mesial-temporal tau
PI-2620 tau to the temporal structural score
temporal structural score to concurrent global cognition
```

Models use heteroscedasticity-consistent type 3 covariance. Age, sex, education, ethnicity, APOE4, MRI scanner, and applicable upstream pathology measurements are included according to the model definitions in `10_4_habs_primary_models.R`. Benjamini-Hochberg correction is applied jointly to the three primary coefficients. Regional tau associations form a separate family. Scanner-specific estimates require at least 100 complete observations.

Longitudinal models evaluate whether structural change predicts later cognition. Follow-up MRI occurs 365 to 1,825 days after baseline. Cognitive follow-up occurs 90 to 730 days after the selected MRI. These analyses report the observed longitudinal boundary alongside the cross-sectional results.

## Running the workflow

Run the main workflow from the repository root:

```r
source("config.R", encoding = "UTF-8")
source("run_all.R", encoding = "UTF-8")
```

The ordered workflow runs the genetic modules, ADNI modules, HABS-HD modules including age modification, clinical-expression models, and final manuscript figure assembly.

Individual modules can be run after prerequisite outputs are present:

```r
source("03_independent_pathway_replication.R", encoding = "UTF-8")
source("09_adni_pathway_prs_analysis.R", encoding = "UTF-8")
source("10_habs_multimodal_analysis.R", encoding = "UTF-8")
source("10_6_habs_age65_multimodal_interactions.R", encoding = "UTF-8")
source("11_clinical_expression.R", encoding = "UTF-8")
source("13_manuscript_figures.R", encoding = "UTF-8")
```

Every analysis directory records model tables, sample counts, audit tables, decision files where applicable, and `session_info.txt`. Figure exports use PDF for vector output and TIFF at 600 dots per inch where raster output is requested.

## Result checks

Successful reproduction should recover the following result pattern:

| Analysis | Expected pattern |
|---|---|
| Independent pathway replication | Two of five prespecified pathways pass false discovery rate control |
| Variant replication | rs56368748 has directionally concordant effects |
| ADNI APP-related pathway score | Positive continuous Centiloid association with stable direction across center omissions |
| HABS-HD amyloid and tau | Positive association |
| HABS-HD tau and temporal structure | Positive association with the vulnerability-oriented score |
| HABS-HD temporal structure and cognition | Greater structural vulnerability accompanies poorer concurrent cognition |
| HABS-HD age modification | Younger participants show weaker amyloid-to-tau and tau-to-structure coupling in the prespecified interaction models |
| HABS-HD longitudinal structure and cognition | No supported association in the available longitudinal sample |

Exact coefficients and sample sizes are written by the scripts and can be checked against the article source-data tables.

## Reproducibility records

Analysis tables include the fitted formula, focal coefficient, standard error, confidence interval, P value, adjusted P value where applicable, and complete-case sample size. Date-matching audits report eligible pairs under each window. Scanner audits report complete observations before scanner-stratified models are fitted. Age-modification outputs include group counts, continuous-age interactions, stratified estimates, diagnosis-adjusted estimates, Vida1 sensitivity estimates, and a decision table. Random procedures use fixed seeds declared in the corresponding script.

## License

Code is distributed under the MIT License. Source datasets remain governed by their original terms.

## Citation

Please cite the associated article when using this workflow. Dataset-specific acknowledgments and citations should follow the requirements of FinnGen, EADB, Bradley and colleagues, ADNI, A4, HABS-HD, and the Gene Expression Omnibus record for GSE272082.
