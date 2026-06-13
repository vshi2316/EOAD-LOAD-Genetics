# Genetic Heterogeneity of Early-Onset versus Late-Onset Alzheimer's Disease

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the analysis code for a multi-layer investigation of genetic heterogeneity between early-onset Alzheimer's disease (EOAD) and late-onset Alzheimer's disease (LOAD). The study integrates GWAS summary statistics, MAGMA gene and pathway analyses, spatial transcriptomic annotation, TWAS, pathway-specific polygenic risk scores, ADNI/A4/HABS/AIBL phenotype-level contextualization, independent EOAD GWAS replication using Bradley/Pottier 2025 summary statistics, meta-enhanced EOAD genetic analysis, and orthogonal single-nucleus transcriptomic contextualization using GSE272082.

The central model is that EOAD and LOAD can show substantial genome-wide overlap while still differing at the pathway and cellular-context levels. FinnGen EOAD and EADB LOAD are used as the primary discovery framework because they allow a parallel EOAD-versus-LOAD comparison. Bradley/Pottier 2025 EOAD GWAS summary statistics are used as independent external EOAD genetic replication, with Age65 as the primary replication definition and Age70 as sensitivity analysis. FinnGen plus Bradley Age65 inverse-variance meta-analysis is used only as a meta-enhanced power analysis after external replication. ADNI, A4, HABS, and AIBL are used for clinical and phenotype-level contextualization. GSE272082 is used for orthogonal single-nucleus transcriptomic contextualization and is not presented as an independent genetic replication dataset.

## Repository Structure

```text
01. EOAD_LOAD_Pathway_Discovery.R
02. PRS_LDpred2_Analysis.R
03. ADNI_Validation.R
04. A4_HABS_AIBL_Validation.R
05_1. EOAD Internal Robustness Analyses
05_2. WMH Sensitivity Analyses
06. Concordance and Phenotype-Level Contextualization
07. Bradley_Pottier_2025_EOAD_GWAS_Replication.R
08. Bradley_Pottier_2025_SNP_Locus_Replication.R
09. FinnGen_Bradley_Age65_IVW_Meta_MAGMA.R
10. GSE272082 Orthogonal Single-Nucleus Transcriptomic Validation
LICENSE
README.md
```

## Evidence Framework

The current manuscript uses a three-layer evidence chain:

1. Discovery: FinnGen EOAD and EADB LOAD define the primary EOAD-versus-LOAD genetic heterogeneity framework.
2. Independent EOAD genetic replication: Bradley/Pottier 2025 EOAD GWAS summary statistics test whether the main EOAD pathway architecture is externally supported.
3. Biological and clinical contextualization: meta-enhanced MAGMA, spatial transcriptomics, TWAS, ADNI/A4/HABS/AIBL, and GSE272082 place the replicated genetic signals into brain, cell-type, and phenotype contexts.

The repository deliberately separates independent genetic replication from contextualization. This distinction is important for interpretation:

- Bradley Age65 and Age70 are independent EOAD genetic replication resources.
- FinnGen plus Bradley Age65 meta-analysis is a power-enhancement analysis, not a new independent validation cohort.
- ADNI is a genotype-enabled clinical and neuroimaging contextualization cohort, not an independent EOAD GWAS replication cohort.
- A4, HABS, and AIBL provide phenotype-level context and sensitivity analyses, not GWAS replication.
- GSE272082 provides orthogonal single-nucleus transcriptomic support, not donor-level genetic replication.

## Data Inputs

The scripts expect local access to public or controlled-access data files. Large GWAS, genotype, and imaging data are not redistributed in this repository.

### GWAS Summary Statistics

Required or optional GWAS inputs include:

- FinnGen Release 11 EOAD summary statistics for the `AD_EO_EXMORE` endpoint.
- EADB LOAD summary statistics.
- Bradley/Pottier 2025 EOAD GWAS summary statistics:
  - `NHW_Age65_sumstats_base_model_SNP-SEX-10PCs-Array.txt`
  - `NHW_Age70_sumstats_base_model_SNP-SEX-10PCs-Array.txt`
- Multivariate healthspan/lifespan/longevity GWAS summary statistics.
- UK Biobank/Oxford BIG40 imaging genetics summary statistics for white-matter microstructure phenotypes.

### Reference Resources

MAGMA analyses require:

- MAGMA executable.
- 1000 Genomes Phase 3 EUR MAGMA binary reference.
- MAGMA gene-location file matching the coordinate build used for SNP annotation.
- Optional MSigDB MAGMA gene-set annotation file.

Bradley/Pottier coordinate processing may require:

- `hg38ToHg19.over.chain.gz` if the raw Bradley file is in GRCh38 and the MAGMA reference is GRCh37/hg19.
- Bioconductor `SNPlocs.Hsapiens.dbSNP155.GRCh37` for chr:pos-to-rsID mapping.

### Clinical and Transcriptomic Data

Clinical and transcriptomic contextualization modules use:

- ADNI genotype, MRI, CSF, APOE, diagnosis, and longitudinal clinical data.
- A4, HABS, and AIBL phenotype-level datasets where available.
- GTEx v8 S-PrediXcan brain tissue prediction outputs.
- LIBD hippocampus spatial transcriptomic data.
- Maynard et al. DLPFC spatial transcriptomic data.
- GSE272082 filtered single-nucleus gene-expression h5 matrices and sample metadata.

## Recommended Local Configuration

Set paths at runtime with environment variables rather than editing scripts. Example:

```r
Sys.setenv(
  EOAD_PROJECT_DIR = "H:/AD_project",
  EOAD_DATA_DIR = "H:/AD_project/data",
  EOAD_RESULTS_DIR = "H:/AD_project/results",
  FINNGEN_EOAD_SUMSTATS = "H:/AD_project/data/GWAS/finngen_R11_AD_EO_EXMORE.gz",
  BRADLEY_RAW_DIR = "H:/AD_project/data/Bradley_Pottier_2025",
  BRADLEY_BUILD = "hg38",
  TARGET_MAGMA_BUILD = "hg19",
  HG38_TO_HG19_CHAIN = "H:/AD_project/data/reference/hg38ToHg19.over.chain.gz",
  MAGMA_EXE = "C:/tools/magma/magma.exe",
  MAGMA_BFILE = "H:/AD_project/data/reference/magma/g1000_eur/g1000_eur",
  MAGMA_GENE_LOC = "H:/AD_project/data/reference/magma/ENSGv110.coding.genes.txt",
  MAGMA_SET_ANNOT = "H:/AD_project/data/reference/magma/MSigDB_20231Hs_MAGMA.txt",
  CANDIDATE_GENESET_FILE = "H:/AD_project/data/candidate_gene_sets/eoad_load_candidate_gene_sets.tsv"
)
```

All newly added Bradley/meta scripts use standard R, Bioconductor, and command-line MAGMA calls. They do not require local wrapper packages.

## Software Dependencies

Core R packages used across modules include:

```r
install.packages(c(
  "data.table", "dplyr", "tidyr", "readr", "stringr", "tibble", "purrr",
  "ggplot2", "ggrepel", "scales", "forcats", "broom", "survival",
  "lmtest", "sandwich", "pROC", "cowplot", "RColorBrewer"
))
```

Selected Bioconductor packages used by specific modules include:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "GenomicRanges", "IRanges", "rtracklayer", "BSgenome",
  "SNPlocs.Hsapiens.dbSNP155.GRCh37", "clusterProfiler",
  "org.Hs.eg.db", "enrichplot", "ComplexHeatmap", "edgeR"
))
```

Additional package requirements depend on the module:

- LDpred2 and genotype handling: `bigsnpr`, `bigstatsr`.
- Single-nucleus RNA-seq: `Seurat`, `SeuratObject`, `Matrix`, `edgeR`.
- TWAS and external command-line analyses: S-PrediXcan and MAGMA installed separately.

## Analysis Modules

### 01. EOAD/LOAD Pathway Discovery

Script:

```r
source("01. EOAD_LOAD_Pathway_Discovery.R")
```

Purpose:

- Read precomputed FinnGen EOAD and EADB LOAD MAGMA gene-level outputs.
- Identify EOAD and LOAD pathway architectures using MAGMA, curated module enrichment, GSEA, gene-property analysis, conditional regression, and LOAD downsampling.
- Establish the primary discovery framework for EOAD-versus-LOAD heterogeneity.

Main outputs:

- EOAD and LOAD MAGMA gene-level summaries.
- Graham et al. glial module enrichment tables.
- GO-BP GSEA outputs.
- Cell-type marker gene-property outputs.
- LOAD downsampling results.
- Publication-ready pathway and cell-type specificity figures.

Interpretation:

- This is the discovery module.
- LOAD downsampling helps distinguish power-dependent enrichment from disease-subtype architecture.

### 02. PRS Construction and TWAS

Script:

```r
source("02. PRS_LDpred2_Analysis.R")
```

Purpose:

- Construct genome-wide and pathway-specific PRS using LDpred2.
- Partition PRS into biologically defined pathway components.
- Perform S-PrediXcan TWAS across GTEx v8 brain tissues for selected EOAD-relevant genes.

Main analyses:

- LDpred2-auto chain fitting and posterior effect aggregation.
- Pathway-masked PRS construction.
- APOE-excluded PRS sensitivity.
- Brain-region TWAS of genetically predicted expression.

Interpretation:

- This module links GWAS architecture to pathway-specific genetic risk scores and genetically predicted brain expression.
- The output supports downstream ADNI clinical contextualization.

### 03. ADNI Clinical and Neuroimaging Contextualization

Script:

```r
source("03. ADNI_Validation.R")
```

Purpose:

- Project pathway-specific PRS into ADNI genotype-enabled participants.
- Test associations with MRI, WMH, CSF biomarkers, and MCI-to-dementia conversion.

Main analyses:

- Linear models for neuroimaging and CSF phenotypes.
- Logistic or Cox models for clinical conversion.
- Age-stratified and age-interaction analyses.
- Sliding-window age analyses.
- Exploratory subtype analyses when retained.

Interpretation:

- ADNI is used as clinical and neuroimaging contextualization.
- ADNI should not be described as independent EOAD GWAS replication.
- PRS results should be interpreted as phenotype-level translation of upstream genetic architecture.

### 04. A4, HABS, and AIBL Phenotype-Level Contextualization

Script:

```r
source("04. A4_HABS_AIBL_Validation.R")
```

Purpose:

- Extend phenotype-level context across preclinical and clinical cohorts.
- Evaluate WMH, cognition, amyloid/tau-related measures, and conversion-related outcomes where available.

Main analyses:

- HABS mediation and structural-equation style contextual analyses.
- A4 phenotype-level models of WMH and cognition.
- AIBL conversion and diagnostic-context analyses.
- Cross-cohort forest plots and participant-flow summaries.

Interpretation:

- These cohorts provide translational phenotype context.
- They should not be used to claim independent genetic replication unless cohort-specific genome-wide genotype and PRS analyses are explicitly available and adequately powered.

### 05_1. EOAD Internal Robustness Analyses

Script:

```r
source("05_1. EOAD Internal Robustness Analyses")
```

Purpose:

- Quantify limits of FinnGen EOAD discovery power.
- Test internal robustness of TWAS and pathway conclusions.

Main analyses:

- Minimum detectable effect/power calculations.
- Leave-one-gene TWAS sensitivity.
- APOE/MHC exclusion sensitivity.
- Matched random gene-set permutation.

Interpretation:

- This module documents discovery-stage uncertainty and guards against overinterpretation of small EOAD discovery sample size.

### 05_2. WMH Sensitivity Analyses

Script:

```r
source("05_2. WMH Sensitivity Analyses")
```

Purpose:

- Evaluate scanner and acquisition-related robustness of ADNI WMH/MRI results.

Main analyses:

- FreeSurfer and UCD FLAIR metadata matching.
- Field-strength sensitivity.
- Manufacturer/model sensitivity.
- Batch-effect residualization.

Interpretation:

- This module supports robustness of imaging contextualization.
- It does not create a new genetic replication layer.

### 06. Concordance and Phenotype-Level Contextualization Framework

Script:

```r
source("06. Concordance and Phenotype-Level Contextualization")
```

Purpose:

- Provide a late-stage support framework separating discovery, published EOAD gene-set concordance, phenotype-level context, and boundary/sensitivity analyses.

Main analyses:

- MAGMA gene-level overview.
- Unbiased GO/KEGG summary.
- Published EOAD gene-set concordance.
- APOE/chromosome 19 sensitivity.
- Targeted 44-gene oligodendrocyte/myelin boundary analysis.
- ADNI/HABS/A4/AIBL phenotype-level context summaries.

Interpretation:

- The targeted 44-gene oligodendrocyte/myelin result should be treated as a boundary or sensitivity analysis, not as the primary replication result.
- This distinction is important because the Bradley/Pottier independent EOAD GWAS replication provides the stronger genetic validation layer.

### 07. Bradley/Pottier 2025 EOAD GWAS Replication

Script:

```r
source("07. Bradley_Pottier_2025_EOAD_GWAS_Replication.R")
```

Purpose:

- Use Bradley/Pottier 2025 EOAD summary statistics as independent external EOAD genetic replication.
- Primary replication uses Age65.
- Age70 is retained as a sensitivity analysis.

Main analyses:

- Standardise Bradley Age65 and Age70 summary statistics.
- Optionally lift GRCh38 coordinates to GRCh37/hg19.
- Map chr:pos to rsID with standard Bioconductor dbSNP resources.
- Run MAGMA gene-level and gene-set analyses.
- Compare Bradley Age65/Age70 MAGMA outputs with FinnGen EOAD and LOAD MAGMA outputs.
- Test pre-specified candidate pathway gene sets with APOE-included and APOE-excluded sensitivity.

Main outputs:

```text
results/07_Bradley_Pottier_2025_EOAD_GWAS_replication/tables/Table_07_MAGMA_gene_level_overview.tsv
results/07_Bradley_Pottier_2025_EOAD_GWAS_replication/tables/Table_07_Bradley_vs_FinnGen_EOAD_gene_level_concordance.tsv
results/07_Bradley_Pottier_2025_EOAD_GWAS_replication/tables/Table_07_Bradley_vs_FinnGen_EOAD_MAGMA_gene_set_concordance.tsv
results/07_Bradley_Pottier_2025_EOAD_GWAS_replication/tables/Table_07_candidate_pathway_replication_by_MAGMA_rank.tsv
results/07_Bradley_Pottier_2025_EOAD_GWAS_replication/tables/Table_07_candidate_pathway_replication_wide.tsv
results/07_Bradley_Pottier_2025_EOAD_GWAS_replication/figures/Figure_07_Bradley_Age65_gene_level_concordance_APOE_excluded.png
results/07_Bradley_Pottier_2025_EOAD_GWAS_replication/figures/Figure_07_Bradley_candidate_pathway_replication_APOE_excluded.png
```

Interpretation:

- This is the primary independent EOAD genetic replication module.
- Results should be written as external support for the pathway-level EOAD genetic architecture.
- Non-replicating candidate sets, especially narrow oligodendrocyte/myelin boundary sets, should not be overstated.

### 08. Bradley/Pottier 2025 SNP/Locus-Level Replication

Script:

```r
source("08. Bradley_Pottier_2025_SNP_Locus_Replication.R")
```

Purpose:

- Test whether FinnGen EOAD SNP-level and locus-level signals show directional support in Bradley Age65/Age70.

Main analyses:

- Harmonise FinnGen EOAD and Bradley variants by rsID and alleles.
- Summarise SNP-level directional concordance across predefined FinnGen P-value thresholds.
- Define approximately independent FinnGen lead loci and test Bradley support within locus windows.
- Evaluate candidate pathway gene intervals for SNP-level support, using APOE-excluded primary reporting.

Main outputs:

```text
results/08_Bradley_Pottier_2025_SNP_locus_replication/tables/Table_08_SNP_harmonization_QC.tsv
results/08_Bradley_Pottier_2025_SNP_locus_replication/tables/Table_08_SNP_directional_replication_summary.tsv
results/08_Bradley_Pottier_2025_SNP_locus_replication/tables/Table_08_FinnGen_EOAD_independent_lead_SNPs.tsv
results/08_Bradley_Pottier_2025_SNP_locus_replication/tables/Table_08_Bradley_locus_level_replication_of_FinnGen_leads.tsv
results/08_Bradley_Pottier_2025_SNP_locus_replication/tables/Table_08_locus_replication_summary.tsv
results/08_Bradley_Pottier_2025_SNP_locus_replication/tables/Table_08_candidate_gene_locus_SNP_support_summary.tsv
results/08_Bradley_Pottier_2025_SNP_locus_replication/figures/Figure_08_SNP_effect_concordance_FinnGen_vs_Bradley_Age65_p1e5.png
results/08_Bradley_Pottier_2025_SNP_locus_replication/figures/Figure_08_candidate_gene_locus_SNP_support_APOE_excluded.png
```

Interpretation:

- This module complements script 07 by testing SNP and locus directionality.
- Locus support should be interpreted as directional external support, not as fine-mapped causality.

### 09. FinnGen EOAD + Bradley Age65 IVW Meta-Enhanced MAGMA

Script:

```r
source("09. FinnGen_Bradley_Age65_IVW_Meta_MAGMA.R")
```

Purpose:

- Combine FinnGen EOAD and Bradley Age65 EOAD summary statistics using fixed-effect inverse-variance weighted meta-analysis.
- Re-run MAGMA and candidate pathway support using the meta-enhanced EOAD summary statistics.
- Export optional meta-enhanced PRS weights.

Main analyses:

- Allele harmonisation between FinnGen EOAD and Bradley Age65.
- Fixed-effect IVW meta-analysis.
- Heterogeneity statistics, including Cochran's Q and I2.
- MAGMA gene-level and gene-set analyses.
- Candidate pathway support across FinnGen EOAD, Bradley Age65, meta-enhanced EOAD, and LOAD.

Main outputs:

```text
results/09_FinnGen_Bradley_Age65_IVW_meta_MAGMA/tables/FinnGen_Bradley_Age65_IVW_meta_full_summary_statistics.tsv
results/09_FinnGen_Bradley_Age65_IVW_meta_MAGMA/PRS_weights/FinnGen_Bradley_Age65_IVW_meta_PRS_weights.tsv
results/09_FinnGen_Bradley_Age65_IVW_meta_MAGMA/tables/Table_09_meta_harmonization_QC.tsv
results/09_FinnGen_Bradley_Age65_IVW_meta_MAGMA/tables/Table_09_meta_summary_QC.tsv
results/09_FinnGen_Bradley_Age65_IVW_meta_MAGMA/tables/Table_09_MAGMA_gene_level_overview.tsv
results/09_FinnGen_Bradley_Age65_IVW_meta_MAGMA/tables/Table_09_MAGMA_gene_set_overview.tsv
results/09_FinnGen_Bradley_Age65_IVW_meta_MAGMA/tables/Table_09_candidate_pathway_gene_set_support_by_MAGMA_rank.tsv
results/09_FinnGen_Bradley_Age65_IVW_meta_MAGMA/figures/Figure_09_meta_candidate_pathway_heatmap_APOE_excluded.png
results/09_FinnGen_Bradley_Age65_IVW_meta_MAGMA/figures/Figure_09_top1000_meta_SNPs.png
```

Interpretation:

- This is a meta-enhanced analysis after external replication.
- It should not be described as independent replication.
- If optional PRS weights are projected into ADNI, SNP coverage must be checked carefully. Low pathway SNP coverage should be treated as technical limitation rather than biological null evidence.

### 10. GSE272082 Orthogonal Single-Nucleus Transcriptomic Validation

Script:

```r
source("10. GSE272082 Orthogonal Single-Nucleus Transcriptomic Validation")
```

Purpose:

- Use public GSE272082 single-nucleus gene-expression data from sporadic EOAD and control brain samples to contextualize candidate genetic pathways at region and cell-type levels.

Main analyses:

- Read filtered 10x gene-expression h5 matrices.
- Reconstruct major cell classes using marker genes aligned with the original study framework.
- Apply marker-signal confidence filtering.
- Aggregate pseudobulk expression by brain region, reconstructed cell type, and sample.
- Test candidate pathway gene-set shifts in sEOAD versus control within region-cell-type strata.
- Export manuscript-facing FDR-significant pseudobulk candidate pathway signals.

Main outputs:

```text
results/10_GSE272082_orthogonal_validation/tables/candidate_gene_set_tests_significant.tsv
results/10_GSE272082_orthogonal_validation/tables/candidate_gene_set_tests_manuscript_focus.tsv
results/10_GSE272082_orthogonal_validation/figures/Figure_GSE272082_pseudobulk_gene_set_signals.png
results/10_GSE272082_orthogonal_validation/analysis_summary.txt
```

Interpretation:

- GSE272082 supports orthogonal biological contextualization.
- Donor-level module-score separation may be underpowered and should not be used as the main validation claim.
- Pseudobulk pathway-level results can be used as supportive orthogonal evidence when they are aligned with the genetic pathway architecture.

## Recommended Run Order

Run the modules in this order for the current manuscript structure:

```r
source("01. EOAD_LOAD_Pathway_Discovery.R")
source("02. PRS_LDpred2_Analysis.R")
source("03. ADNI_Validation.R")
source("04. A4_HABS_AIBL_Validation.R")
source("05_1. EOAD Internal Robustness Analyses")
source("05_2. WMH Sensitivity Analyses")
source("06. Concordance and Phenotype-Level Contextualization")
source("07. Bradley_Pottier_2025_EOAD_GWAS_Replication.R")
source("08. Bradley_Pottier_2025_SNP_Locus_Replication.R")
source("09. FinnGen_Bradley_Age65_IVW_Meta_MAGMA.R")
source("10. GSE272082 Orthogonal Single-Nucleus Transcriptomic Validation")
```

For the Bradley/meta section alone, run:

```r
source("07. Bradley_Pottier_2025_EOAD_GWAS_Replication.R")
source("08. Bradley_Pottier_2025_SNP_Locus_Replication.R")
source("09. FinnGen_Bradley_Age65_IVW_Meta_MAGMA.R")
```

Script 08 depends on the Bradley MAGMA input files generated by script 07. Script 09 also depends on the Bradley Age65 MAGMA input generated by script 07. Script 10 is independent of scripts 07 to 09 and uses GSE272082 single-nucleus RNA-seq data.

## Candidate Gene-Set File Format

Scripts 07 to 09 expect a pre-specified candidate gene-set file. The recommended format is tab-separated:

```text
set_name    gene
Adhesion    APP
Adhesion    NCAM1
Amyloid_Clearance    APOE
Immune_Microglia    TREM2
Glutamate_Synapse    GRIN2A
```

Accepted gene columns include `gene`, `GENE`, `gene_id`, `magma_gene_id`, `symbol`, or `gene_symbol`. Accepted set columns include `set_name`, `pathway`, `gene_set`, or `module`.

Use the same candidate gene sets across FinnGen, Bradley, meta-analysis, and GSE272082 to avoid post hoc pathway selection.

## Output Numbering and Manuscript Mapping

The manuscript-facing supplementary table order for the newly added evidence layer is:

- Bradley/Pottier independent EOAD GWAS replication tables first.
- FinnGen plus Bradley Age65 meta-enhanced tables next.
- GSE272082 orthogonal single-nucleus transcriptomic tables after the Bradley/meta tables.

The corresponding manuscript logic is:

- Methods: independent EOAD GWAS replication and meta-analysis are described before GSE272082.
- Results: Bradley/Pottier Age65/Age70 replication and meta-enhanced analysis are reported before GSE272082.
- Main figure: the most important Bradley/meta/GSE272082 panels can be combined as a final integrated figure, with all extended QC and sensitivity outputs in supplementary materials.

## Interpretation Guardrails

Use the following language boundaries when writing the manuscript:

- Appropriate: “Bradley/Pottier Age65 provided independent EOAD summary-statistics replication.”
- Appropriate: “Age70 was used as a sensitivity analysis across a broader early/younger-onset definition.”
- Appropriate: “FinnGen plus Bradley Age65 meta-analysis improved power for downstream gene and pathway analyses.”
- Appropriate: “GSE272082 provided orthogonal single-nucleus transcriptomic contextualization.”
- Appropriate: “ADNI provided clinical and neuroimaging contextualization of pathway-specific genetic risk.”
- Avoid: “GSE272082 independently replicated the genetic findings.”
- Avoid: “ADNI independently replicated the EOAD GWAS.”
- Avoid: “Meta-analysis is an independent validation dataset.”
- Avoid: “The targeted 44-gene oligodendrocyte/myelin set is the primary replicated genetic mechanism” unless the independent replication results support that claim.

## Reproducibility Notes

- Large data files are not included in this repository.
- Controlled-access or consortium data should be obtained under the appropriate data-use agreements.
- Path variables should be set with `Sys.setenv()` for local runs.
- MAGMA output file names can differ by MAGMA version. Scripts 07 and 09 attempt to detect common `.genes.out`, `.genes.out.txt`, `.gsa.out`, and `.gsa.out.txt` forms.
- If Bradley/Pottier files have already been converted to rsID/hg19 MAGMA input format, set `BRADLEY_AGE65_MAGMA_INPUT` and `BRADLEY_AGE70_MAGMA_INPUT` to skip raw coordinate conversion.
- If raw Bradley/Pottier files are used, install the required Bioconductor packages and provide the correct chain file when coordinate conversion is needed.
- If optional upstream inputs are unavailable, affected downstream summaries should be marked as not run rather than silently omitted.

## Citation

If you use this code, please cite the associated manuscript and the original data resources, including FinnGen, EADB, Bradley/Pottier 2025 EOAD GWAS, ADNI, A4, HABS, AIBL, GTEx, LIBD, Maynard et al. DLPFC spatial transcriptomics, and GSE272082.

## License

This repository is released under the MIT License. See `LICENSE` for details.
