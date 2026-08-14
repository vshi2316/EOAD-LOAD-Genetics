# Copy this file to config.R and enter the locations of controlled-access data.

Sys.setenv(
  EOAD_DATA_DIR = file.path(getwd(), "data"),
  EOAD_RESULTS_DIR = file.path(getwd(), "results"),
  HABS_DATA_DIR = file.path(getwd(), "data", "HABS_HD"),
  HABS_FBB_FILE = file.path(getwd(), "data", "HABS_HD", "FBB_PET.xlsx"),
  HABS_TAU_FILE = file.path(getwd(), "data", "HABS_HD", "PI2620_PET.xlsx"),
  HABS_CORTICAL_THICKNESS_FILE = file.path(getwd(), "data", "HABS_HD", "cortical_thickness.xlsx"),
  HABS_GENOMICS_FILE = file.path(getwd(), "data", "HABS_HD", "genomics.xlsx"),
  HABS_DIAGNOSIS_FILE = file.path(getwd(), "data", "HABS_HD", "diagnosis.csv"),
  HABS_DEMOGRAPHICS_FILE = file.path(getwd(), "data", "HABS_HD", "demographics.csv"),
  HABS_COGNITION_FILE = file.path(getwd(), "data", "HABS_HD", "cognition_zscores.csv"),
  HABS_MMSE_FILE = file.path(getwd(), "data", "HABS_HD", "mmse.csv"),
  HABS_EDUCATION_FILE = file.path(getwd(), "data", "HABS_HD", "education.csv"),
  PLINK_EXE = "plink",
  ADNI_PLINK_PREFIX = file.path(getwd(), "data", "ADNI", "genotype", "adni_qc"),
  ADNI_PATHWAY_WEIGHT_FILE = file.path(getwd(), "data", "ADNI", "pathway_score_weights.tsv"),
  ADNI_CENTILOID_COVARIATE_FILE = file.path(getwd(), "data", "ADNI", "ADNI_Centiloid_covariates.tsv.gz"),
  ADNI_PATHWAY_ANALYSIS_FILE = file.path(getwd(), "data", "ADNI", "ADNI_pathway_PRS_Centiloid.tsv.gz"),
  A4_CLINICAL_ANALYSIS_FILE = file.path(getwd(), "data", "A4", "A4_Centiloid_PACC.tsv.gz"),
  ADNI_CLINICAL_ANALYSIS_FILE = file.path(getwd(), "data", "ADNI", "ADNI_Centiloid_ADAS13.tsv.gz")
)
