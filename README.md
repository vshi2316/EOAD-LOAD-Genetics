[README.md](https://github.com/user-attachments/files/30106346/README.md)
# Early-Onset Alzheimer's Disease Within the Broader AD/ADRD Genetic Architecture

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

The repository contains analysis code for a multi-layer study of early-onset Alzheimer's disease (EOAD) within the broader genetic architecture of Alzheimer's disease and related dementia (AD/ADRD). The European Alzheimer & Dementia Biobank (EADB) stage I meta-analysis supplies the broad AD/ADRD reference architecture. FinnGen Release 11 EOAD summary statistics identify focused genetic, pathway, spatial, and transcriptome-wide enrichment patterns against that reference.

Independent EOAD genetic evaluation uses Bradley/Pottier 2025 multi-cohort summary statistics. Age65 is the primary external dataset, and Age70 is retained as a sensitivity definition. FinnGen EOAD and Bradley Age65 are subsequently combined through fixed-effect inverse-variance meta-analysis to increase gene-level and pathway-level resolution.

Spatial transcriptomic resources, GTEx v8 transcriptome-wide association analyses, and GSE272082 single-nucleus RNA sequencing provide anatomical and cellular context. The Alzheimer's Disease Neuroimaging Initiative (ADNI) provides individual-level phenotype evaluation of pathway-specific polygenic risk scores across imaging, cerebrospinal fluid biomarkers, and clinical progression. A4, HABS, and AIBL extend the phenotype-level context.

## Evidence Framework

The manuscript uses the following evidence hierarchy:

1. **Broader reference architecture:** EADB AD/ADRD summary statistics define the large-sample reference profile.
2. **Focused EOAD findings:** FinnGen EOAD identifies pathway, spatial, transcriptome-wide, and polygenic patterns within the shared architecture.
3. **External EOAD genetic replication:** Bradley/Pottier Age65 evaluates prespecified FinnGen loci and pathways. Age70 supplies a broader sensitivity analysis.
4. **Meta-enhanced integration:** FinnGen EOAD and Bradley Age65 meta-analysis increases statistical resolution after external evaluation.
5. **Transcriptomic context:** gsMap, GTEx v8, and GSE272082 localize genetic findings across brain regions and cell types.
6. **Individual-level phenotype evaluation:** ADNI evaluates pathway-specific scores against imaging, biomarker, and progression phenotypes. A4, HABS, and AIBL provide supporting phenotype-level analyses.

The analytical roles remain distinct. Bradley/Pottier supplies independent EOAD genetic replication. Meta-analysis increases power. GSE272082 supplies transcriptomic context. ADNI supplies genotype-enabled phenotype evaluation. A4, HABS, and AIBL extend the clinical context.

## Repository Structure

```text
01. EOAD_EADB_Pathway_Discovery.R
02. PRS_LDpred2_Analysis.R
03. ADNI_Validation.R
04. A4_HABS_AIBL_Validation.R
05_1. EOAD Internal Robustness Analyses.R
05_2. WMH Sensitivity Analyses.R
06. Concordance and Phenotype-Level Contextualization Framework.R
07. Bradley 2025 EOAD GWAS Replication.R
08. Bradley 2025 SNP_Locus-Level Replication.R
09. FinnGen EOAD + Bradley Age65 IVW Meta-Enhanced MAGMA Analysis.R
10. GSE272082 Orthogonal Single-Nucleus Transcriptomic Validation.R
LICENSE
README.md
```

All analysis scripts use the `.R` extension. Comparator identifiers, objects, environment variables, output names, and manuscript-facing labels use `EADB`, with `EADB AD/ADRD comparator` retained as the full scientific description.

## Data Inputs

Large genotype, GWAS, imaging, and transcriptomic files are not redistributed in this repository. Access remains subject to the terms of the original data providers.

### GWAS Summary Statistics

- FinnGen Release 11 `AD_EO_EXMORE` summary statistics for EOAD.
- EADB stage I AD/ADRD meta-analysis summary statistics.
- Bradley/Pottier 2025 non-Hispanic White Age65 and Age70 EOAD summary statistics.
- Multivariate healthspan, parental lifespan, and longevity summary statistics.
- UK Biobank Oxford BIG40 imaging-genetics summary statistics for white matter microstructure.

The Bradley/Pottier input files used by modules 07 to 09 are expected to follow the source-study base model:

```text
NHW_Age65_sumstats_base_model_SNP-SEX-10PCs-Array.txt
NHW_Age70_sumstats_base_model_SNP-SEX-10PCs-Array.txt
```

### Reference Resources

- MAGMA executable.
- 1000 Genomes Phase 3 European reference panel.
- MAGMA gene-location file matched to the target genomic build.
- Gene-set annotations from the selected MSigDB release.
- GRCh38-to-GRCh37 chain file when Bradley coordinates require liftOver.
- `SNPlocs.Hsapiens.dbSNP155.GRCh37` for coordinate-to-rsID mapping.
- HapMap3+ European linkage disequilibrium reference for LDpred2.

### Clinical and Transcriptomic Resources

- ADNI genotype, magnetic resonance imaging, white matter hyperintensity, cerebrospinal fluid biomarker, APOE, diagnosis, and longitudinal follow-up data.
- A4, HABS, and AIBL phenotype data used under the corresponding data-use agreements.
- GTEx v8 S-PrediXcan prediction models and brain-tissue outputs.
- LIBD human hippocampal spatial transcriptomic data.
- Maynard et al. dorsolateral prefrontal cortex spatial transcriptomic data.
- GSE272082 filtered single-nucleus expression matrices and sample metadata.

## Local Configuration

Local paths should be supplied through environment variables rather than hard-coded manuscript paths.

```r
Sys.setenv(
  EOAD_PROJECT_DIR = "H:/AD_project",
  EOAD_DATA_DIR = "H:/AD_project/data",
  EOAD_RESULTS_DIR = "H:/AD_project/results",
  FINNGEN_EOAD_SUMSTATS = "H:/AD_project/data/GWAS/finngen_R11_AD_EO_EXMORE.gz",
  EADB_ADRD_SUMSTATS = "H:/AD_project/data/GWAS/EADB_stage1_AD_ADRD.tsv.gz",
  BRADLEY_RAW_DIR = "H:/AD_project/data/Bradley_Pottier_2025",
  BRADLEY_BUILD = "hg38",
  TARGET_MAGMA_BUILD = "hg19",
  HG38_TO_HG19_CHAIN = "H:/AD_project/data/reference/hg38ToHg19.over.chain.gz",
  MAGMA_EXE = "C:/tools/magma/magma.exe",
  MAGMA_BFILE = "H:/AD_project/data/reference/magma/g1000_eur/g1000_eur",
  MAGMA_GENE_LOC = "H:/AD_project/data/reference/magma/ENSGv110.coding.genes.txt",
  MAGMA_SET_ANNOT = "H:/AD_project/data/reference/magma/MSigDB_20231Hs_MAGMA.txt"
)
```

## Software Dependencies

Core R packages include:

```r
install.packages(c(
  "data.table", "dplyr", "tidyr", "readr", "stringr", "tibble", "purrr",
  "ggplot2", "ggrepel", "scales", "forcats", "broom", "survival",
  "lmtest", "sandwich", "pROC", "cowplot", "RColorBrewer",
  "bigsnpr", "bigstatsr"
))
```

Selected Bioconductor packages include:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GenomicRanges", "IRanges", "rtracklayer", "BSgenome",
  "SNPlocs.Hsapiens.dbSNP155.GRCh37", "clusterProfiler",
  "org.Hs.eg.db", "enrichplot", "ComplexHeatmap", "edgeR"
))
```

Additional software includes MAGMA, S-PrediXcan, and the reference resources required by each program.

## Analysis Modules

### 01. EOAD Pathway Discovery Against the EADB AD/ADRD Reference

```r
source("01. EOAD_EADB_Pathway_Discovery.R")
```

Module 01 reads FinnGen EOAD and EADB AD/ADRD MAGMA results and performs:

- MAGMA gene-level association summaries.
- Glial module enrichment.
- Gene Ontology Biological Process gene set enrichment analysis.
- Cell-type gene-property analysis.
- Conditional immune-cell models adjusted for microglial expression.
- Effective-sample-size attenuation of EADB summary statistics as a descriptive power-sensitivity analysis.

The module compares enrichment profiles. The EADB resource represents a broader AD/ADRD reference because contributing cohorts did not apply a uniform age-at-onset threshold.

### 02. LDpred2 Polygenic Architecture and Transcriptome-Wide Association Analysis

```r
source("02. PRS_LDpred2_Analysis.R")
```

Module 02 performs:

- LDpred2-auto chain fitting and posterior effect aggregation.
- Genome-wide and pathway-specific score construction.
- APOE-region exclusion sensitivity analysis.
- S-PrediXcan analysis across GTEx v8 brain tissues.
- Leave-one-gene, locus-exclusion, and matched gene-set sensitivity analyses where configured.

Pathway-specific scores are approximate burden metrics derived by retaining posterior SNP weights near pathway genes. Interpretation should account for linkage disequilibrium between pathway and non-pathway variants.

### 03. ADNI Individual-Level Phenotype Evaluation

```r
source("03. ADNI_Validation.R")
```

Module 03 projects pathway-specific scores into ADNI and evaluates associations with:

- Structural magnetic resonance imaging measures.
- White matter hyperintensity volume.
- Cerebrospinal fluid amyloid, tau, phosphorylated tau, and sTREM2.
- Mild cognitive impairment-to-dementia conversion.
- Age-stratified and continuous age-interaction models.
- Scanner-aware imaging sensitivity analyses.

ADNI analyses focus on continuous pathway-score associations with imaging, biomarker, and progression phenotypes. Continuous score models carry greater inferential weight than exploratory PRS profiles or extreme-group contrasts.

### 04. A4, HABS, and AIBL Phenotype-Level Context

```r
source("04. A4_HABS_AIBL_Validation.R")
```

Module 04 evaluates:

- White matter hyperintensity volume and cognition in A4.
- White matter hyperintensity volume, plasma p-tau217, and cognition in HABS.
- APOE epsilon 4 and clinical progression in AIBL.
- Participant flow and sensitivity to influential observations.

A4, HABS, and AIBL extend phenotype-level context through white matter, biomarker, cognition, and progression analyses.

### 05_1. EOAD Internal Robustness Analyses

```r
source("05_1. EOAD Internal Robustness Analyses.R")
```

Analyses include minimum detectable effects, leave-one-gene sensitivity, APOE and MHC exclusion, and matched random gene-set permutation. The resulting estimates quantify uncertainty arising from the modest effective sample size of FinnGen EOAD.

### 05_2. White Matter Hyperintensity Sensitivity Analyses

```r
source("05_2. WMH Sensitivity Analyses.R")
```

Analyses include field-strength, manufacturer, scanner-model, and acquisition-related sensitivity checks that evaluate imaging robustness within ADNI.

### 06. Concordance, Phenotype Context, and Boundary Analyses

```r
source("06. Concordance and Phenotype-Level Contextualization Framework.R")
```

Module 06 separates:

- Genome-wide discovery in FinnGen EOAD and EADB AD/ADRD.
- Concordance with published Bradley/Pottier evidence classes.
- Phenotype-level analyses in ADNI, A4, HABS, and AIBL.
- Targeted 44-gene oligodendrocyte and myelin boundary analyses.

The restricted 44-gene results remain supplementary because matched permutation, single-cell eQTL colocalization, and SMR or HEIDI screening did not support a concentrated candidate-gene mechanism.

### 07. Bradley/Pottier EOAD External Genetic Replication

```r
source("07. Bradley 2025 EOAD GWAS Replication.R")
```

Age65 supplies the primary external EOAD dataset, and Age70 supplies sensitivity analysis. The module performs coordinate harmonization, rsID mapping, MAGMA gene analysis, genome-wide gene-set concordance, and same-method testing of the five pathways prespecified from FinnGen EOAD discovery.

Age65 supported amyloid-beta clearance and negative regulation of amyloid precursor protein catabolism after false discovery rate correction. Age70 supported all five pathways in sensitivity analysis. The primary replication claim is restricted to the Age65-supported findings.

### 08. Bradley/Pottier SNP and Locus Replication

```r
source("08. Bradley 2025 SNP_Locus-Level Replication.R")
```

Module 08 performs:

- Allele harmonization between FinnGen and Bradley/Pottier.
- Linkage disequilibrium clumping of FinnGen EOAD variants.
- Direction testing across independent non-APOE loci.
- Exact-variant or highest-linkage-disequilibrium-proxy evaluation.
- Separate summary of the APOE region through one index variant.

Among 25 represented non-APOE loci, 14 had concordant directions in Age65. The directional test was not significant. rs56368748 replicated after correction, and rs405509 supplied strong APOE-region support.

### 09. FinnGen EOAD and Bradley Age65 Meta-Enhanced Analysis

```r
source("09. FinnGen EOAD + Bradley Age65 IVW Meta-Enhanced MAGMA Analysis.R")
```

Module 09 performs allele harmonization, fixed-effect inverse-variance meta-analysis, Cochran Q and I-squared calculations, MAGMA gene analysis, and same-method evaluation of the five FinnGen-defined pathways.

Following external evaluation, the meta-analysis increases statistical resolution through combined FinnGen and Bradley evidence.

### 10. GSE272082 Single-Nucleus Transcriptomic Context

```r
source("10. GSE272082 Orthogonal Single-Nucleus Transcriptomic Validation.R")
```

Module 10 reads filtered 10x expression matrices, reconstructs major cell classes using the original-study marker framework, applies confidence filtering, aggregates donor-level pseudobulk expression, and tests candidate pathway shifts across brain region and cell type.

The resulting signals provide region-resolved and cell-type-resolved transcriptomic context, with inference restricted to pathway expression patterns across the available donors.

## Recommended Run Order

```r
source("01. EOAD_EADB_Pathway_Discovery.R")
source("02. PRS_LDpred2_Analysis.R")
source("03. ADNI_Validation.R")
source("04. A4_HABS_AIBL_Validation.R")
source("05_1. EOAD Internal Robustness Analyses.R")
source("05_2. WMH Sensitivity Analyses.R")
source("06. Concordance and Phenotype-Level Contextualization Framework.R")
source("07. Bradley 2025 EOAD GWAS Replication.R")
source("08. Bradley 2025 SNP_Locus-Level Replication.R")
source("09. FinnGen EOAD + Bradley Age65 IVW Meta-Enhanced MAGMA Analysis.R")
source("10. GSE272082 Orthogonal Single-Nucleus Transcriptomic Validation.R")
```

Module 08 uses harmonized Bradley inputs prepared by module 07. Module 09 uses the Bradley Age65 input generated or validated by module 07. Module 10 is independent of modules 07 to 09.

## EADB Comparator Labels

Internal objects and exported fields use `EADB`, including names such as `EADB_Microglia` and `PRS_EADB_*`. Plot labels, table titles, legends, captions, and exported manuscript-facing fields use the following mapping:

```r
trait_display <- c(
  EOAD = "FinnGen EOAD",
  EADB = "EADB AD/ADRD comparator",
  Aging = "Multivariate aging"
)

display_trait <- function(x) {
  dplyr::recode(as.character(x), !!!trait_display, .default = as.character(x))
}
```

Use `EOAD minus EADB comparator` for difference scores and `effective-sample-size-attenuated EADB comparator` for the deterministic power-sensitivity analysis.

## Interpretation Boundaries

Supported manuscript language includes:

- Bradley/Pottier Age65 supplied independent EOAD summary-statistics replication.
- Age70 supplied a broader sensitivity analysis.
- FinnGen plus Age65 meta-analysis increased power for gene and pathway analysis.
- GSE272082 supplied single-nucleus transcriptomic context.
- ADNI supplied individual-level phenotype evaluation of pathway-specific scores.
- A4, HABS, and AIBL extended phenotype-level context.

## Reproducibility Notes

- Set local paths with `Sys.setenv()`.
- Preserve the source genomic build and record each liftOver step.
- Report variant counts before and after coordinate and rsID harmonization.
- Use linkage disequilibrium-independent loci for directional replication tests.
- Represent the APOE region with a single index variant in locus-level summaries.
- Freeze FinnGen-defined replication pathways before testing Bradley Age65.
- Record the multiple-testing family used for each analysis.
- Report unavailable optional inputs as unavailable rather than omitting the corresponding output silently.
- Confirm pathway SNP coverage before interpreting a pathway-specific score.

## Data Availability

FinnGen, EADB, Bradley/Pottier, GTEx, LIBD, GSE272082, and the imaging-genetics resources are available through their original repositories or consortium portals. ADNI, A4, HABS, and AIBL data require approval under their respective access procedures. Controlled-access participant data are excluded from the repository.

## Citation

Please cite the associated manuscript and the original data resources used by each analysis, including FinnGen, EADB, Bradley/Pottier 2025, ADNI, A4, HABS, AIBL, GTEx, LIBD, Maynard et al. dorsolateral prefrontal cortex spatial transcriptomics, and GSE272082.

## License

The analysis code is released under the MIT License. See `LICENSE` for details.
