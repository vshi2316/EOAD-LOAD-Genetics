[Uploading README.md…]()
# Early-Onset Alzheimer Disease Within the Broader Alzheimer Disease Genetic Architecture

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Scope

This repository contains the reproducible analysis code for a multi-layer study of early-onset Alzheimer disease (EOAD) within the broader Alzheimer disease and related dementia genetic architecture. The study uses FinnGen EOAD summary statistics as the focused discovery resource, the EADB stage I AD and ADRD meta-analysis as the broader reference architecture, and Bradley and Pottier multi-cohort EOAD summary statistics as an independent genetic evaluation.

The analysis proceeds from genetic discovery to independent replication, meta-enhanced prioritization, molecular and cellular context, and individual-level clinical expression. The clinical expression layer uses harmonized analyses across A4, ADNI and HABS-HD. AIBL is retained as a longitudinal conversion benchmark when the available integrated file does not contain quantitative pathology and WMH variables compatible with the shared model.

The repository contains code and configuration templates only. Controlled-access participant data, genotype files, imaging files, GWAS summary statistics, transcriptomic matrices and software binaries are not redistributed.

## Evidence hierarchy

The manuscript distinguishes the following evidence levels.

1. The EADB AD and ADRD resource provides the broader genetic reference architecture.
2. FinnGen EOAD supplies the focused discovery analysis.
3. Bradley and Pottier Age65 supplies the primary independent EOAD summary-statistics evaluation. Age70 is a sensitivity definition.
4. FinnGen EOAD and Bradley Age65 are combined in a fixed-effect inverse-variance meta-analysis after the external evaluation.
5. Spatial transcriptomic, GTEx v8, and GSE272082 analyses provide molecular and cellular context.
6. ADNI provides genotype-enabled phenotype evaluation of pathway-specific scores.
7. A4, ADNI and HABS-HD provide harmonized pathology-to-cognition analyses across clinical settings. AIBL provides longitudinal clinical progression context.

These layers answer different questions. External EOAD data support genetic replication, meta-analysis improves statistical resolution, transcriptomic analyses localize candidate biology, and the clinical cohorts characterize how pathology relates to cognition across disease stages.

## Repository layout

```text
00_preflight.R
01. EOAD_EADB_Pathway_Discovery.R
02. PRS_LDpred2_Analysis.R
03. ADNI_Validation.R
04_1. A4_HABS-HD_AIBL_Validation.R
04_2. Harmonized Clinical Expression Across Alzheimer Disease Stages.R
05_1. EOAD Internal Robustness Analyses.R
05_2. WMH Sensitivity Analyses.R
06. Concordance and Phenotype-Level Contextualization Framework.R
07. Bradley 2025 EOAD GWAS Replication.R
08. Bradley 2025 SNP_Locus-Level Replication.R
09. FinnGen EOAD + Bradley Age65 IVW Meta-Enhanced MAGMA Analysis.R
10. GSE272082 Orthogonal Single-Nucleus Transcriptomic Validation.R
config_template.R
LICENSE
README.md
```

The numbered files define the manuscript-facing analytical order. Module `04_2` contains the complete harmonized clinical analysis in one script, including model construction, heterogeneity analysis, bounded sensitivity checks and Figure 6 generation.

## Manuscript-to-code correspondence

| Manuscript component | Code module | Main role |
|---|---|---|
| Shared EOAD and broader AD genetic architecture | `01` | FinnGen EOAD versus EADB pathway, gene-property and enrichment analyses |
| Polygenic architecture and transcriptome-wide analysis | `02` | LDpred2 scores, APOE-region sensitivity and S-PrediXcan analyses |
| ADNI pathway-score phenotype evaluation | `03` | Imaging, CSF, cognition and progression associations |
| Cohort-specific phenotype context | `04_1` | Original A4, HABS and AIBL context analyses retained for supplementary reporting |
| Harmonized clinical expression | `04_2` | Shared pathology-to-cognition model, WMH and age modifiers, stage heterogeneity and Figure 6 |
| EOAD internal robustness | `05_1` | Leave-one-gene, APOE and MHC exclusion, detectable-effect and matched gene-set analyses |
| Imaging sensitivity | `05_2` | Scanner, field-strength, manufacturer and acquisition sensitivity analyses |
| Concordance and boundary diagnostics | `06` | Published gene-set concordance, phenotype boundaries and targeted candidate-gene diagnostics |
| Independent EOAD summary-statistics replication | `07` | Bradley and Pottier Age65 primary analysis and Age70 sensitivity analysis |
| SNP and locus replication | `08` | Allele harmonization, LD clumping, directional testing and rs56368748 evaluation |
| Meta-enhanced gene and pathway analysis | `09` | FinnGen plus Bradley Age65 fixed-effect IVW meta-analysis |
| Single-nucleus transcriptomic context | `10` | GSE272082 cell-type reconstruction, pseudobulk analysis and candidate pathway context |

The manuscript methods describe harmonized clinical expression before independent replication because it is a prespecified phenotype layer. The results place the clinical expression results after genetic replication and molecular context so that the clinical analyses are interpreted as translational contextualization rather than as a substitute for EOAD genetic replication.

## Analysis order

Run the modules in this order after configuring local paths:

```r
source("config.R")
source("00_preflight.R")

source("01. EOAD_EADB_Pathway_Discovery.R")
source("02. PRS_LDpred2_Analysis.R")
source("03. ADNI_Validation.R")
source("04_1. A4_HABS-HD_AIBL_Validation.R")
source("04_2. Harmonized Clinical Expression Across Alzheimer Disease Stages.R")
source("05_1. EOAD Internal Robustness Analyses.R")
source("05_2. WMH Sensitivity Analyses.R")
source("06. Concordance and Phenotype-Level Contextualization Framework.R")
source("07. Bradley 2025 EOAD GWAS Replication.R")
source("08. Bradley 2025 SNP_Locus-Level Replication.R")
source("09. FinnGen EOAD + Bradley Age65 IVW Meta-Enhanced MAGMA Analysis.R")
source("10. GSE272082 Orthogonal Single-Nucleus Transcriptomic Validation.R")
```

Modules `01` to `06` contain functions and documented entry points. Modules `07` to `10` contain executable workflows. Module `04_2` is executable and runs the complete harmonized clinical-expression layer through one call.

## Configuration

Copy `config_template.R` to an untracked file named `config.R` and edit the local paths. Do not commit `config.R` when it contains controlled-access locations.

```r
source("config.R")
source("00_preflight.R")
```

The configuration uses the following variables.

| Variable | Purpose |
|---|---|
| `EOAD_PROJECT_DIR` | Root directory for local data, results and reference resources |
| `EOAD_DATA_DIR` | Local data directory |
| `EOAD_RESULTS_DIR` | Results directory |
| `FINNGEN_EOAD_SUMSTATS` | FinnGen EOAD summary statistics |
| `EADB_ADRD_SUMSTATS` | EADB stage I AD and ADRD summary statistics |
| `BRADLEY_RAW_DIR` | Bradley and Pottier input directory |
| `MAGMA_EXE` | MAGMA executable |
| `MAGMA_BFILE` | MAGMA reference genotype prefix |
| `MAGMA_GENE_LOC` | MAGMA gene-location file |
| `MAGMA_SET_ANNOT` | MAGMA gene-set annotation file |
| `HG38_TO_HG19_CHAIN` | Coordinate conversion chain when required |
| `GSE272082_DATA_DIR` | GSE272082 expression and metadata directory |
| `EOAD_A4_INTEGRATED_FILE` | Subject-level A4 integrated file for harmonized clinical expression |
| `EOAD_HABS_HD_INTEGRATED_FILE` | Subject-level HABS-HD integrated file |
| `EOAD_AIBL_INTEGRATED_FILE` | AIBL integrated longitudinal file |
| `EOAD_ADNI_CENTILOID_FILE` | ADNI Centiloid analysis file |
| `EOAD_ADNI_ADAS13_FILE` | ADNI ADAS13 analysis file |
| `EOAD_ADNI_WMH_FILE` | ADNI WMH derivative file |

The clinical module accepts CSV and TSV files. Values coded as empty strings, `NA`, `N/A`, `-999` or `-99` are treated as missing.

## Clinical input requirements

The harmonized clinical module requires the following columns in the integrated cohort files.

### A4

`BID`, `Centiloid`, `Log_WMH`, `PACC`, `Age`, `Gender`, `Education`, and `APOE4_Carrier`.

### HABS-HD

`Med_ID`, `pTau217`, `Log_WMH`, `MMSE`, `Age`, `Gender`, `Education`, and `APOE4_Carrier`.

The integrated HABS-HD `pTau217` field is treated as the supplied transformed or standardized scale. The script does not apply a second logarithmic transformation.

### ADNI

The Centiloid and ADAS13 files require `RID`, `date`, `Centiloid`, `baseline_age`, `sex`, `education`, `APOE4`, `baseline_dx`, and `time`. The WMH file requires `RID`, `EXAMDATE`, and `TOTAL_WMH`.

ADAS13 and WMH observations are matched to the closest compatible visit within the prespecified windows of 90 and 180 days, respectively. ADNI longitudinal models retain a participant-level random intercept and adjust for time and baseline diagnosis.

### AIBL

The AIBL audit expects an identifier column and the longitudinal `Time` and `Event` fields. AIBL is included in the shared pathology model only when the integrated file also contains a quantitative pathology field and a compatible WMH field. With the current integrated file, it remains a longitudinal conversion benchmark.

## Harmonized clinical-expression module

Run the complete layer with:

```r
source("04_2. Harmonized Clinical Expression Across Alzheimer Disease Stages.R")
```

The module executes four linked components.

### Component 1: harmonized primary models

The primary question is how pathology burden relates to cognitive impairment after adjustment for age, sex, education and APOE ε4 carrier status.

The three prespecified layers are:

```text
cognition impairment ~ pathology burden + covariates
cognition impairment ~ pathology burden * WMH + covariates
cognition impairment ~ pathology burden * age + covariates
```

All continuous variables are z-standardized within cohort and the cognitive outcome is oriented so that larger values indicate greater impairment. A4 and ADNI use Centiloid as the amyloid measure. HABS-HD uses p-tau217 as a downstream pathology marker and is reported separately from the Centiloid-only comparison.

The primary cross-sectional models use one complete observation per participant. ADNI observations are ordered by participant and pathology date before the first eligible observation is retained for the participant-level model.

HC3 robust standard errors are used for the primary linear models. The script writes the complete-case sample size, coefficient, standard error, confidence interval, nominal P value, model status and within-layer FDR-adjusted P value.

### Component 2: ADNI longitudinal sensitivity

ADNI longitudinal observations are matched to cognition and WMH measurements using the same visit windows. The random-intercept models include time, baseline age, baseline diagnosis, sex, education and APOE ε4. These models provide longitudinal sensitivity evidence and are not presented as EOAD-specific replication.

### Component 3: stage and common-support heterogeneity

The A4 to ADNI difference is evaluated using prespecified clinical-stage and age-support analyses. ADNI is summarized across baseline CN, EMCI, LMCI and AD groups, with a cognitively normal sensitivity analysis and a common 65 to 85 year range. The analysis reports cohort-specific slopes, pairwise contrasts and I-squared values. No sliding age windows, data-driven cohort removal or two-study meta-regression is used.

### Component 4: bounded diagnostic sensitivity and Figure 6

The validation layer evaluates baseline-diagnosis adjustment, stage-matched cognitively normal analyses, common-age analyses, quadratic pathology terms and standardized-residual influence checks. It reports whether the pathology main effect changes by more than 20 percent after influence filtering and whether the quadratic term is supported.

The Figure 6 script uses the covariate-aligned pathology main effects and partial-residual summaries. It writes PDF, SVG, PNG and TIFF files under `results/clinical_expression_figures`. It does not write to a manuscript directory.

## Clinical interpretation boundaries

The harmonized clinical layer supports the following interpretation.

- A4 and ADNI provide compatible Centiloid-based amyloid-to-cognition estimates.
- HABS-HD provides an independent community-sample p-tau217-to-cognition context using a different pathology scale.
- AIBL provides longitudinal clinical conversion context when quantitative pathology and WMH are unavailable.
- WMH interactions are reported with cohort-specific estimates and heterogeneity. They are not promoted to a universal modifier without consistent support in compatible cohorts.
- Age interactions are reported as prespecified continuous interactions and bounded support analyses. Sliding-window positives are not used as primary evidence.
- The clinical layer characterizes clinical expression of Alzheimer pathology. It does not constitute direct replication of the EOAD genetic pathways.

## Main outputs

Module `04_2` writes the following result directories below `EOAD_RESULTS_DIR`.

### `clinical_expression`

```text
primary_three_layers_HC3.tsv
ADNI_longitudinal_time_adjusted.tsv
meta_three_layers.tsv
clinical_evidence_gates.tsv
AIBL_role_audit.tsv
source_cohort_sizes.tsv
sessionInfo.txt
```

### `clinical_expression_heterogeneity`

```text
pathology_cognition_stage_estimates.tsv
prespecified_slope_contrasts.tsv
heterogeneity_stage_common_support.tsv
WMH_interaction_context.tsv
age_interaction_context.tsv
heterogeneity_decisions.tsv
sessionInfo.txt
```

### `clinical_expression_validation`

```text
covariate_aligned_models.tsv
ADNI_CN_stage_matched.tsv
primary_model_diagnostics.tsv
analysis_decisions.tsv
sessionInfo.txt
```

### `clinical_expression_figures`

```text
Figure_6.pdf
Figure 6.svg
Figure 6.png
Figure 6.tiff
Figure6_panel_A_source.tsv
Figure6_panels_B_D_source.tsv
Figure6_panel_E_source.tsv
Figure6_panel_F_source.tsv
sessionInfo.txt
```

## Existing supplementary analyses

Module `04_1` remains in the repository because the manuscript retains the original cohort-specific phenotype context in the supplementary material. It covers A4 WMH-cognition association, HABS-HD WMH and p-tau217 analyses, mediation and structural equation models, age-stratified context, clinical utility diagnostics, and AIBL APOE progression. Those analyses are interpreted as phenotype context and are separate from the harmonized primary model in `04_2`.

## Earlier modules

### `01. EOAD_EADB_Pathway_Discovery.R`

Reads FinnGen EOAD and EADB results and evaluates MAGMA gene-level associations, gene-set enrichment, cell-type gene-property models and effective-sample-size sensitivity analyses.

### `02. PRS_LDpred2_Analysis.R`

Fits LDpred2 scores, constructs pathway-specific scores, evaluates APOE-region exclusion and runs configured S-PrediXcan analyses. Pathway scores are burden metrics based on posterior SNP weights near pathway genes and must be interpreted with LD structure in mind.

### `03. ADNI_Validation.R`

Projects pathway-specific genetic scores into ADNI and tests imaging, CSF biomarker, cognition and progression phenotypes. These models are distinct from the harmonized pathology-expression models in `04_2`.

### `05_1. EOAD Internal Robustness Analyses.R`

Evaluates detectable effects, leave-one-gene sensitivity, APOE and MHC exclusion and matched random gene-set permutations.

### `05_2. WMH Sensitivity Analyses.R`

Evaluates scanner, field-strength, manufacturer and acquisition sensitivity within ADNI imaging analyses.

### `06. Concordance and Phenotype-Level Contextualization Framework.R`

Separates genome-wide discovery, published gene-set concordance, phenotype-level context and the targeted boundary analyses.

### `07. Bradley 2025 EOAD GWAS Replication.R`

Harmonizes coordinates and alleles and evaluates the prespecified FinnGen findings in Bradley and Pottier Age65, with Age70 sensitivity analysis.

### `08. Bradley 2025 SNP_Locus-Level Replication.R`

Performs LD-clumped locus comparisons, direction testing, exact-variant or proxy evaluation and separate APOE-region indexing.

### `09. FinnGen EOAD + Bradley Age65 IVW Meta-Enhanced MAGMA Analysis.R`

Performs fixed-effect inverse-variance meta-analysis followed by gene-level and pathway-level MAGMA analyses.

### `10. GSE272082 Orthogonal Single-Nucleus Transcriptomic Validation.R`

Reconstructs major cell classes, aggregates donor-level pseudobulk counts and evaluates candidate pathway expression context across region and cell type.

## Software requirements

Core R packages include:

```r
install.packages(c(
  "data.table", "dplyr", "tidyr", "readr", "stringr", "tibble", "purrr",
  "ggplot2", "ggrepel", "scales", "forcats", "broom", "survival",
  "lmtest", "sandwich", "pROC", "cowplot", "RColorBrewer", "bigsnpr",
  "bigstatsr", "nlme", "patchwork", "svglite", "ragg"
))
```

The discovery and transcriptomic modules also require the Bioconductor packages listed by `00_preflight.R`, together with MAGMA, S-PrediXcan and their reference resources where applicable.

## Reproducibility safeguards

- Freeze FinnGen-defined pathways before testing Bradley Age65.
- Record the number of variants and genes before and after coordinate harmonization.
- Use one APOE-region index variant for locus-level summaries.
- Keep pathology scales separate when their biological measurements differ.
- Report complete-case sample sizes for each model rather than treating them as the cohort total.
- Preserve cohort-specific estimates when heterogeneity is substantial.
- Treat unavailable AIBL pathology or WMH fields as unavailable rather than silently dropping the cohort.
- Record `sessionInfo()` in every clinical output directory.

## Data availability

FinnGen, EADB, Bradley and Pottier, GTEx, LIBD, GSE272082 and imaging-genetics resources are obtained from their original repositories or consortium portals. ADNI, A4, HABS-HD and AIBL require approval under their respective data-use procedures. No controlled-access participant data are included here.

## Citation

Please cite the associated manuscript and the original data resources used by each module, including FinnGen, EADB, Bradley and Pottier, ADNI, A4, HABS-HD, AIBL, GTEx, LIBD, Maynard et al. dorsolateral prefrontal cortex spatial transcriptomics and GSE272082.

## License

The analysis code is released under the MIT License. See `LICENSE` for details.
