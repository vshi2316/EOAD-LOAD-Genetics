# EOAD-LOAD-Genetics
# Genetic Heterogeneity of Early-Onset vs Late-Onset Alzheimer's Disease

## Repository Structure

```
EOAD-LOAD-Genetics/
├── 01_Pathway_Discovery/
│   └── EOAD_LOAD_Pathway_Discovery.R      # WGCNA, MAGMA, GSEA analysis
├── 02_PRS_Construction/
│   └── PRS_LDpred2_Analysis.R             # Polygenic risk score construction
├── 03_ADNI_Validation/
│   └── ADNI_Validation.R                  # ADNI cohort validation
├── 04_External_Cohort_Validation/
│   └── A4_HABS_AIBL_Validation.R          # A4, HABS, AIBL validation
├── LICENSE
└── README.md
```

---

## Analysis Descriptions

### 1. Data-Driven Pathway Discovery (`01_Pathway_Discovery/`)

**EOAD_LOAD_Pathway_Discovery.R** - Comprehensive pathway analysis integrating:

- **WGCNA**: Weighted gene co-expression network analysis on GSE137313 mouse hippocampus RNA-seq (N=48 samples across 4 age groups: 8, 16, 32, 72 weeks)
- **MAGMA v1.10**: Gene-based association with boundaries extended 35kb upstream/10kb downstream, Bonferroni correction (P < 0.05/19,000 genes)
- **GSEA**: Gene set enrichment using clusterProfiler (GO-BP, KEGG pathways, FDR < 0.05)
- **Cell Type Analysis**: MAGMA gene-property with curated marker sets (T-cell: 64 genes, Microglia: 58 genes, Oligodendrocyte: 62 genes)
- **Graham et al. Module Enrichment**: Validation using published human microglial (26 genes), human oligodendrocyte (29 genes), mouse ARM (17 genes), mouse phagolysosomal (45 genes), and mouse longevity (38 genes) signatures
- **Conditional Analysis**: Testing T-cell independence controlling for microglial expression (Z_MAGMA = β₀ + β₁×T_Cell + β₂×Microglia + ε)

### 2. Polygenic Risk Score Construction (`02_PRS_Construction/`)

**PRS_LDpred2_Analysis.R** - Bayesian PRS construction:

- **LDpred2-auto**: 30 Markov chains, 500 burn-in, 500 iterations, sparse=TRUE, shrink_corr=0.95
- **Chain Bagging**: Aggregating posterior estimates across converged chains
- **APOE Exclusion**: Extended APOE region (chr19: 44-46 Mb) masked
- **Pathway-specific PRS** (SNPs within 100kb of pathway genes):
  - T-cell activation (56 genes)
  - Activated microglia (58 genes)
  - Amyloid-beta clearance (52 genes)
  - APP metabolism (38 genes)
  - Oligodendrocyte function (62 genes)
  - Myelination (52 genes)

### 3. ADNI Validation (`03_ADNI_Validation/`)

**ADNI_Validation.R** - Primary validation cohort analysis:

- **Neuroimaging**: FreeSurfer UCSFFSX51 pipeline
  - Total WMH volume (ST128SV), ICV-normalized
  - Bilateral hippocampus (ST29SV + ST88SV)
  - Bilateral entorhinal cortex (ST24CV + ST83CV)
  - Hippocampal subfields (CA1: ST131HS + ST139HS, Subiculum: ST137HS + ST145HS)
- **CSF Biomarkers**: AlzBio3 platform, restricted to ADNI phases 1/GO/2
  - Aβ42, total tau, phosphorylated tau
  - Amyloid positivity threshold: Aβ42 < 192 pg/mL
- **Cox Regression**: MCI-to-dementia conversion using `coxph(Surv(Time_Years, Event) ~ PRS + covariates)`
- **K-means Clustering**: Genetic subtype identification (k=3, silhouette validation)
- **Linear Mixed Models**: Longitudinal cognitive trajectory (MMSE, CDRSB, ADAS11, ADAS13)
- **Age Stratification**: <70 vs ≥70 years
- **Covariates**: Age, sex, APOE ε4 count, education (PTEDUCAT)

### 4. External Cohort Validation (`04_External_Cohort_Validation/`)

**A4_HABS_AIBL_Validation.R** - Multi-cohort validation:

- **A4 Study**: Multi-group SEM using lavaan
  - Model: APOE4 → Centiloid → PACC
  - Estimator: MLR (maximum likelihood with robust standard errors)
  - Missing data: FIML (full information maximum likelihood)
  - Age groups: <70 vs ≥70 years
  - Fit indices: CFI, TLI, RMSEA, SRMR
  - Bootstrap CI: 1000 iterations for indirect effects

- **HABS**: Clinical utility assessment
  - Outcome: Cognitive impairment (CDR ≥ 0.5 OR MMSE < 24)
  - Firth-corrected logistic regression
  - AUC with DeLong 95% CI
  - AUPRC (area under precision-recall curve)
  - NRI with 3 threshold sets: high (0.20, 0.40, 0.60), standard (0.10, 0.20), fine (0.15, 0.25, 0.35, 0.50)
  - IDI with bootstrap (1000 iterations)
  - DCA (decision curve analysis, threshold 0.01-0.80)
  - Train/test split: 70/30 with fixed seed

- **AIBL**: Survival analysis
  - Cox proportional hazards regression
  - AFT models with Weibull distribution
  - Time ratio estimates for APOE ε4 effects

---

## Requirements

### R Version
- R >= 4.2.0

### R Packages

```r
# CRAN packages
install.packages(c(
  "data.table", "dplyr", "tidyr", "ggplot2", "cowplot",
  "survival", "survminer", "lme4", "lmerTest",
  "cluster", "factoextra", "pheatmap", "corrplot",
  "pROC", "lavaan", "logistf", "PRROC"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "bigsnpr", "WGCNA", "clusterProfiler", 
  "org.Hs.eg.db", "org.Mm.eg.db", "biomaRt",
  "TxDb.Hsapiens.UCSC.hg19.knownGene"
))
```

---

## Data Availability

### GWAS Summary Statistics

| Dataset | Source | Sample Size | Population | PubMed ID | Access |
|---------|--------|-------------|------------|-----------|--------|
| **EOAD** | FinnGen R11 | Cases: 1,573<br>Controls: 199,505 | European (Finnish) | 36653562 | [Download](https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_AD_EO_EXMORE.gz) |
| **LOAD** | EADB Consortium | Cases: 85,934<br>Controls: 401,577 | European (Multi-country) | 35379992 | [EBI GWAS](https://www.ebi.ac.uk/gwas/) (GCST90027158) |
| **Aging** | Timmers et al. | Healthspan: 300,477<br>Lifespan: 1,012,240<br>Longevity: 36,745 | European | 32678081<br>31413261 | [Edinburgh DataShare](https://datashare.ed.ac.uk/) |
| **TPMI** | Taiwan PMI | Cases: 2,828<br>Controls: 320,866 | East Asian (Han Chinese) | 41094136 | [TPMI PheWeb](https://pheweb.ibms.sinica.edu.tw/) |

#### Dataset Details

**Early-Onset Alzheimer's Disease (EOAD)**
- Definition: Age of onset < 65 years (ICD-10 code G30.0)
- Endpoint: AD_EO_EXMORE (strict exclusion criteria removing other dementia subtypes)
- Genotyping: Illumina Global Screening Array
- Imputation: SISu v3 Finnish reference panel

**Late-Onset Alzheimer's Disease (LOAD)**
- Composition: 39,106 clinically diagnosed + 46,828 proxy cases
- Diagnosis: NINCDS-ADRDA or DSM criteria
- Cohorts: 42 cohorts across Europe and North America
- SNPs: 21,101,114 SNPs post-QC
- Platforms: Affymetrix and Illumina arrays

**Multivariate Aging**
- Healthspan: Years free from major age-related diseases (28.3% unhealthy)
- Parental lifespan: Age at death or current age (60% deceased)
- Longevity: Survival to 90th percentile vs 60th percentile

**Trans-Ancestry AD (Han Chinese)**
- PheCode: 290.11
- Source: Taiwan Biobank electronic medical records
- Genotyping: Axiom Genome-Wide TWB 2.0 Array
- Imputation: 1000 Genomes EAS reference panel

### White Matter Microstructure (UK Biobank BIG40)

| Tract | UKB Field | Sample Size | PubMed ID | Description |
|-------|-----------|-------------|-----------|-------------|
| **Corpus Callosum Body** | 25059 | 33,224 | 30305740 | Primary inter-hemispheric motor/premotor tract |
| **Cingulum Hippocampus (R)** | 25092 | 33,224 | 30305740 | Limbic tract connecting cingulate-entorhinal-hippocampus |
| **Uncinate Fasciculus (R)** | 25100 | 33,224 | 30305740 | Temporal-orbitofrontal connection for semantic memory |
| **Fornix** | 25061 | 33,224 | 30305740 | Primary hippocampal efferent pathway (Papez circuit) |

- Phenotype: Mean fractional anisotropy (FA) on TBSS skeleton
- Method: Diffusion MRI (dMRI) with tract-based spatial statistics
- Access: [UK Biobank BIG40](https://open.win.ox.ac.uk/ukbiobank/big40/)

### Transcriptomic Data

**GSE137313 - Mouse Hippocampus Temporal RNA-seq**
- Species: C57BL/6J wild-type mice
- Timepoints: 8 weeks (young), 16 weeks (adult), 32 weeks (middle-aged), 72 weeks (aged)
- Samples: N=24 (pooled across 4 timepoints)
- Normalization: Transcripts per million (TPM)
- PubMed ID: 32274467
- Access: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137313) | [Mouseac](http://www.mouseac.org/)

### Validation Cohorts

| Cohort | Description | Access |
|--------|-------------|--------|
| **ADNI** | Alzheimer's Disease Neuroimaging Initiative<br>Genome-wide genotyping, neuroimaging, CSF, longitudinal | [adni.loni.usc.edu](https://adni.loni.usc.edu/) |
| **A4 Study** | Anti-Amyloid Treatment in Asymptomatic AD<br>Cognitively normal, elevated amyloid, Centiloid scale | [ida.loni.usc.edu](https://ida.loni.usc.edu/) |
| **HABS** | Harvard Aging Brain Study<br>Multimodal imaging, plasma pTau217 (ALZpath Simoa v2), WMH | [habs.mgh.harvard.edu](https://habs.mgh.harvard.edu/) |
| **AIBL** | Australian Imaging, Biomarker and Lifestyle<br>Non-North American replication, longitudinal conversion | [aibl.csiro.au](https://aibl.csiro.au/) |

---

## External Tools

Computational workflows employed publicly available implementations:

| Tool | Version | Purpose | URL |
|------|---------|---------|-----|
| **LDSC** | 1.0.1 | Genetic correlation, heritability | [GitHub](https://github.com/bulik/ldsc) |
| **LAVA** | 0.1.0 | Local genetic correlation (~1Mb segments) | [GitHub](https://github.com/josefin-werme/LAVA) |
| **conjFDR** | - | Conditional/conjunctional FDR for pleiotropic variants | [GitHub](https://github.com/precimed/pleiofdr) |
| **FUMA** | v1.5.4 | Functional mapping and annotation | [Web](https://fuma.ctglab.nl/) |
| **gsMap** | - | Spatial transcriptomic mapping (Cauchy combination test) | [GitHub](https://github.com/JianYang-Lab/gsMap) |
| **HESS** | 0.5.4 | Local SNP heritability (1,703 LD-independent regions) | [GitHub](https://github.com/huwenboshi/hess) |
| **Genomic SEM** | 0.0.5 | Structural equation modeling of genetic covariance | [GitHub](https://github.com/GenomicSEM/GenomicSEM) |
| **g:Profiler** | - | GO/KEGG pathway enrichment | [Web](https://biit.cs.ut.ee/gprofiler/gost) |

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---


## Acknowledgments

Data collection and sharing for this project was funded by:
- Alzheimer's Disease Neuroimaging Initiative (ADNI) (National Institutes of Health Grant U01 AG024904 and DOD ADNI W81XWH-12-2-0012)
- FinnGen Consortium
- European Alzheimer & Dementia Biobank (EADB) Consortium
- UK Biobank
- Taiwan Precision Medicine Initiative

We thank all participants and their families for their contributions to research.
