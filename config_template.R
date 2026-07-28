# Copy this file to config.R, replace the example paths, and source it before
# running analysis modules. Do not commit controlled-access paths or credentials.

project_dir <- normalizePath(
  Sys.getenv("EOAD_PROJECT_DIR", unset = getwd()),
  winslash = "/",
  mustWork = FALSE
)
magma_executable <- Sys.getenv(
  "MAGMA_EXE",
  unset = file.path(project_dir, "tools", "magma", "magma")
)

Sys.setenv(
  EOAD_PROJECT_DIR = project_dir,
  EOAD_DATA_DIR = file.path(project_dir, "data"),
  EOAD_RESULTS_DIR = file.path(project_dir, "results"),
  FINNGEN_EOAD_SUMSTATS = file.path(
    project_dir, "data", "GWAS", "finngen_R11_AD_EO_EXMORE.gz"
  ),
  EADB_ADRD_SUMSTATS = file.path(
    project_dir, "data", "GWAS", "EADB_stage1_AD_ADRD.tsv.gz"
  ),
  BRADLEY_RAW_DIR = file.path(
    project_dir, "data", "Bradley_Pottier_2025"
  ),
  BRADLEY_BUILD = "hg38",
  TARGET_MAGMA_BUILD = "hg19",
  HG38_TO_HG19_CHAIN = file.path(
    project_dir, "data", "reference", "hg38ToHg19.over.chain.gz"
  ),
  MAGMA_EXE = magma_executable,
  MAGMA_BFILE = file.path(
    project_dir, "data", "reference", "magma", "g1000_eur", "g1000_eur"
  ),
  MAGMA_GENE_LOC = file.path(
    project_dir, "data", "reference", "magma", "ENSGv110.coding.genes.txt"
  ),
  MAGMA_SET_ANNOT = file.path(
    project_dir, "data", "reference", "magma", "MSigDB_20231Hs_MAGMA.txt"
  ),
  GSE272082_DATA_DIR = file.path(project_dir, "data", "GSE272082")
)
