# ==============================================================================
# EOAD Internal Robustness Analyses
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})


# ==============================================================================
# Utility functions
# ==============================================================================

safe_qcut <- function(x, n_bins = 4) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_integer_, length(x))
  ok <- is.finite(x)

  if (sum(ok) < n_bins || length(unique(x[ok])) < n_bins) {
    return(out)
  }

  ranks <- rank(x[ok], ties.method = "first")
  breaks <- unique(quantile(ranks, probs = seq(0, 1, length.out = n_bins + 1),
                            na.rm = TRUE, type = 7))

  if (length(breaks) <= 2) {
    return(out)
  }

  out[ok] <- as.integer(cut(ranks, breaks = breaks, include.lowest = TRUE,
                           labels = FALSE)) - 1L
  out
}

wilcox_greater <- function(foreground, background) {
  foreground <- foreground[is.finite(foreground)]
  background <- background[is.finite(background)]

  if (length(foreground) < 3 || length(background) < 20) {
    return(NA_real_)
  }

  suppressWarnings(
    wilcox.test(foreground, background, alternative = "greater",
                exact = FALSE)$p.value
  )
}

signed_rank_greater <- function(values) {
  values <- values[is.finite(values)]

  if (length(values) < 5) {
    return(NA_real_)
  }

  suppressWarnings(
    tryCatch(
      wilcox.test(values, mu = 0, alternative = "greater",
                  exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  )
}

in_locus <- function(df, chr, start, stop) {
  chr_vec <- suppressWarnings(as.numeric(df$CHR))
  start_vec <- suppressWarnings(as.numeric(df$START))
  stop_vec <- suppressWarnings(as.numeric(df$STOP))

  chr_vec == chr &
    !is.na(start_vec) &
    !is.na(stop_vec) &
    start_vec <= stop &
    stop_vec >= start
}


# ==============================================================================
# 1. Power / minimum detectable effect analysis
# ==============================================================================

calculate_power_table <- function(eoad_effective_n = 1573,
                                  alpha_levels = list(
                                    "Genome-wide SNP threshold" = 5e-8,
                                    "MAGMA gene Bonferroni threshold" = 2.63e-6,
                                    "Nominal pathway screen" = 0.05
                                  ),
                                  mafs = c(0.01, 0.05, 0.10, 0.20, 0.40)) {
  rows <- list()

  for (label in names(alpha_levels)) {
    alpha <- alpha_levels[[label]]
    zcrit <- qnorm(1 - alpha / 2)
    z80 <- zcrit + qnorm(0.80)
    z90 <- zcrit + qnorm(0.90)

    for (maf in mafs) {
      denom <- sqrt(eoad_effective_n * 2 * maf * (1 - maf))
      beta80 <- z80 / denom
      beta90 <- z90 / denom

      rows[[length(rows) + 1]] <- data.frame(
        Analysis_level = label,
        Alpha = alpha,
        EOAD_effective_N = eoad_effective_n,
        MAF = maf,
        Z_for_80pct_power = z80,
        Approx_min_detectable_OR_80pct_power = exp(beta80),
        Z_for_90pct_power = z90,
        Approx_min_detectable_OR_90pct_power = exp(beta90),
        Interpretation = paste(
          "Approximate single-variant effect boundary;",
          "pathway and gene-set tests can aggregate weaker effects but remain",
          "constrained by EOAD discovery power."
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  bind_rows(rows)
}


# ==============================================================================
# 2. Data loading
# ==============================================================================

load_oligo_genes <- function(oligo_gene_file) {
  genes <- fread(oligo_gene_file)
  stopifnot("Gene" %in% colnames(genes))
  sort(unique(as.character(genes$Gene[!is.na(genes$Gene)])))
}

load_magma_coordinates <- function(magma_gene_file) {
  magma <- fread(magma_gene_file)
  required <- c("Gene_Symbol", "CHR", "START", "STOP", "NSNPS", "ZSTAT", "P")
  missing <- setdiff(required, colnames(magma))
  if (length(missing) > 0) {
    stop("MAGMA gene file is missing columns: ", paste(missing, collapse = ", "))
  }

  magma %>%
    select(all_of(required)) %>%
    distinct(Gene_Symbol, .keep_all = TRUE)
}

load_tissue_twas <- function(tissue_name, twas_file, magma_coords, oligo_genes) {
  twas <- fread(twas_file)
  required <- c("gene_name", "zscore", "pvalue")
  missing <- setdiff(required, colnames(twas))
  if (length(missing) > 0) {
    stop("S-PrediXcan file is missing columns: ", paste(missing, collapse = ", "))
  }

  twas <- twas %>%
    rename(Gene = gene_name, Z = zscore, P_TWAS = pvalue) %>%
    mutate(
      Tissue = tissue_name,
      Gene = as.character(Gene),
      Is_OligoMyelin = Gene %in% oligo_genes
    ) %>%
    left_join(magma_coords, by = c("Gene" = "Gene_Symbol")) %>%
    mutate(
      Gene_Length = suppressWarnings(as.numeric(STOP) - as.numeric(START) + 1),
      Z = suppressWarnings(as.numeric(Z)),
      P_TWAS = suppressWarnings(as.numeric(P_TWAS)),
      var_g = if ("var_g" %in% colnames(.)) suppressWarnings(as.numeric(var_g)) else NA_real_,
      n_snps_in_model = if ("n_snps_in_model" %in% colnames(.)) suppressWarnings(as.numeric(n_snps_in_model)) else NA_real_,
      n_snps_used = if ("n_snps_used" %in% colnames(.)) suppressWarnings(as.numeric(n_snps_used)) else NA_real_,
      NSNPS = suppressWarnings(as.numeric(NSNPS))
    )

  twas
}

add_matching_bins <- function(df) {
  df$n_snps_in_model_bin <- safe_qcut(df$n_snps_in_model, 4)
  df$var_g_bin <- safe_qcut(df$var_g, 4)
  df$NSNPS_bin <- safe_qcut(df$NSNPS, 4)
  df$Gene_Length_bin <- safe_qcut(df$Gene_Length, 4)
  df
}


# ==============================================================================
# 3. Robustness analyses
# ==============================================================================

build_candidate_index_list <- function(oligo_df, background_df) {
  bg <- background_df
  levels <- list(
    c("n_snps_in_model_bin", "var_g_bin", "NSNPS_bin"),
    c("n_snps_in_model_bin", "var_g_bin"),
    c("n_snps_in_model_bin"),
    c("NSNPS_bin")
  )

  candidates <- vector("list", nrow(oligo_df))

  for (i in seq_len(nrow(oligo_df))) {
    selected <- integer(0)

    for (cols in levels) {
      keep <- rep(TRUE, nrow(bg))
      for (col in cols) {
        if (col %in% colnames(bg) && !is.na(oligo_df[[col]][i])) {
          keep <- keep & bg[[col]] == oligo_df[[col]][i]
        }
      }
      idx <- which(keep)
      if (length(idx) > 0) {
        selected <- idx
        break
      }
    }

    if (length(selected) == 0) {
      selected <- seq_len(nrow(bg))
    }

    candidates[[i]] <- selected
  }

  candidates
}

sample_matched_z <- function(background_z, candidate_index_list) {
  vapply(candidate_index_list, function(idx) {
    background_z[sample(idx, size = 1)]
  }, numeric(1))
}

run_tissue_robustness <- function(tissue_name, twas_df, n_perm = 2000,
                                  random_seed = 20260525) {
  set.seed(random_seed)

  df <- twas_df %>%
    filter(is.finite(Z)) %>%
    add_matching_bins()

  oligo_df <- df %>% filter(Is_OligoMyelin)
  background_df <- df %>% filter(!Is_OligoMyelin)

  observed_median <- median(oligo_df$Z, na.rm = TRUE)
  observed_mean <- mean(oligo_df$Z, na.rm = TRUE)
  observed_abs_median <- median(abs(oligo_df$Z), na.rm = TRUE)

  observed_mw_p <- wilcox_greater(oligo_df$Z, background_df$Z)
  observed_signed_p <- signed_rank_greater(oligo_df$Z)

  leave_detail <- lapply(sort(unique(oligo_df$Gene)), function(gene) {
    loo <- oligo_df %>% filter(Gene != gene)
    data.frame(
      Tissue = tissue_name,
      Dropped_gene = gene,
      N_remaining_oligo_myelin_genes = nrow(loo),
      Median_Z_after_drop = median(loo$Z, na.rm = TRUE),
      Mean_Z_after_drop = mean(loo$Z, na.rm = TRUE),
      MannWhitney_greater_P_after_drop = wilcox_greater(loo$Z, background_df$Z),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  leave_summary <- data.frame(
    Tissue = tissue_name,
    N_oligo_myelin_genes = nrow(oligo_df),
    Observed_median_Z = observed_median,
    Observed_mean_Z = observed_mean,
    Observed_MannWhitney_greater_P = observed_mw_p,
    Observed_signed_rank_greater_P = observed_signed_p,
    Leave_one_gene_min_median_Z = min(leave_detail$Median_Z_after_drop, na.rm = TRUE),
    Leave_one_gene_max_median_Z = max(leave_detail$Median_Z_after_drop, na.rm = TRUE),
    Leave_one_gene_min_P = min(leave_detail$MannWhitney_greater_P_after_drop, na.rm = TRUE),
    Leave_one_gene_max_P = max(leave_detail$MannWhitney_greater_P_after_drop, na.rm = TRUE),
    Leave_one_gene_retained_nominal_P_lt_0.05 =
      sum(leave_detail$MannWhitney_greater_P_after_drop < 0.05, na.rm = TRUE),
    Worst_case_dropped_gene =
      leave_detail$Dropped_gene[which.max(leave_detail$MannWhitney_greater_P_after_drop)],
    Interpretation =
      ifelse(max(leave_detail$MannWhitney_greater_P_after_drop, na.rm = TRUE) < 0.05,
             "Robust", "Sensitive or not significant in full model"),
    stringsAsFactors = FALSE
  )

  apoe <- in_locus(df, chr = 19, start = 44000000, stop = 46000000)
  mhc <- in_locus(df, chr = 6, start = 25000000, stop = 34000000)
  locus_masks <- list(
    APOE_chr19_44_46Mb = apoe,
    MHC_chr6_25_34Mb = mhc,
    APOE_or_MHC = apoe | mhc
  )

  locus_results <- lapply(names(locus_masks), function(label) {
    mask <- locus_masks[[label]]
    d2 <- df[!mask, ]
    ol2 <- d2 %>% filter(Is_OligoMyelin)
    bg2 <- d2 %>% filter(!Is_OligoMyelin)

    data.frame(
      Tissue = tissue_name,
      Excluded_locus = label,
      N_total_genes_removed = sum(mask, na.rm = TRUE),
      N_oligo_myelin_genes_removed = sum(mask & df$Is_OligoMyelin, na.rm = TRUE),
      N_oligo_myelin_genes_remaining = nrow(ol2),
      Median_Z_after_locus_exclusion = median(ol2$Z, na.rm = TRUE),
      MannWhitney_greater_P_after_locus_exclusion =
        wilcox_greater(ol2$Z, bg2$Z),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  candidate_indices <- build_candidate_index_list(oligo_df, background_df)
  background_z <- background_df$Z

  random_medians <- numeric(n_perm)
  random_means <- numeric(n_perm)
  random_abs_medians <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    sampled_z <- sample_matched_z(background_z, candidate_indices)
    random_medians[i] <- median(sampled_z, na.rm = TRUE)
    random_means[i] <- mean(sampled_z, na.rm = TRUE)
    random_abs_medians[i] <- median(abs(sampled_z), na.rm = TRUE)
  }

  permutation_summary <- data.frame(
    Tissue = tissue_name,
    N_oligo_myelin_genes = nrow(oligo_df),
    N_background_genes = nrow(background_df),
    N_permutations = n_perm,
    Observed_median_Z = observed_median,
    Random_median_Z_mean = mean(random_medians),
    Random_median_Z_2.5pct = as.numeric(quantile(random_medians, 0.025)),
    Random_median_Z_97.5pct = as.numeric(quantile(random_medians, 0.975)),
    Empirical_P_median_Z_greater =
      (1 + sum(random_medians >= observed_median)) / (n_perm + 1),
    Observed_mean_Z = observed_mean,
    Random_mean_Z_mean = mean(random_means),
    Random_mean_Z_2.5pct = as.numeric(quantile(random_means, 0.025)),
    Random_mean_Z_97.5pct = as.numeric(quantile(random_means, 0.975)),
    Empirical_P_mean_Z_greater =
      (1 + sum(random_means >= observed_mean)) / (n_perm + 1),
    Observed_median_abs_Z = observed_abs_median,
    Random_median_abs_Z_mean = mean(random_abs_medians),
    Empirical_P_median_abs_Z_greater =
      (1 + sum(random_abs_medians >= observed_abs_median)) / (n_perm + 1),
    Observed_percentile_vs_random_median = mean(random_medians <= observed_median),
    Matching_variables =
      "n_snps_in_model, var_g, MAGMA NSNPS bins; relaxed matching used when sparse",
    stringsAsFactors = FALSE
  )

  permutation_diagnostics <- data.frame(
    Tissue = tissue_name,
    Random_median_Z = random_medians[seq_len(min(1000, n_perm))],
    Random_mean_Z = random_means[seq_len(min(1000, n_perm))],
    Random_median_abs_Z = random_abs_medians[seq_len(min(1000, n_perm))],
    stringsAsFactors = FALSE
  )

  list(
    leave_summary = leave_summary,
    leave_detail = leave_detail,
    locus_results = locus_results,
    permutation_summary = permutation_summary,
    permutation_diagnostics = permutation_diagnostics
  )
}


# ==============================================================================
# 4. Main execution wrapper
# ==============================================================================

run_eoad_internal_robustness <- function(
    magma_gene_file,
    oligo_gene_file,
    twas_files,
    output_dir = "results/eoad_internal_robustness",
    eoad_effective_n = 1573,
    n_perm = 2000,
    random_seed = 20260525) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  cat("=== EOAD Internal Robustness Analyses ===\n\n")
  cat(sprintf("Output directory: %s\n", output_dir))
  cat(sprintf("Matched random gene-set permutations per tissue: %d\n\n", n_perm))

  cat("[1/4] Calculating minimum detectable effects...\n")
  power_table <- calculate_power_table(eoad_effective_n = eoad_effective_n)
  fwrite(power_table,
         file.path(output_dir, "S40_EOAD_power_minimum_detectable_effect.csv"))

  cat("[2/4] Loading TWAS and MAGMA inputs...\n")
  oligo_genes <- load_oligo_genes(oligo_gene_file)
  magma_coords <- load_magma_coordinates(magma_gene_file)

  all_results <- list()

  cat("[3/4] Running leave-one-gene, locus-exclusion and permutation analyses...\n")
  for (tissue_name in names(twas_files)) {
    cat(sprintf("  %s\n", tissue_name))
    twas_df <- load_tissue_twas(
      tissue_name = tissue_name,
      twas_file = twas_files[[tissue_name]],
      magma_coords = magma_coords,
      oligo_genes = oligo_genes
    )

    all_results[[tissue_name]] <- run_tissue_robustness(
      tissue_name = tissue_name,
      twas_df = twas_df,
      n_perm = n_perm,
      random_seed = random_seed
    )
  }

  leave_summary <- bind_rows(lapply(all_results, `[[`, "leave_summary"))
  leave_detail <- bind_rows(lapply(all_results, `[[`, "leave_detail"))
  locus_results <- bind_rows(lapply(all_results, `[[`, "locus_results"))
  permutation_summary <- bind_rows(lapply(all_results, `[[`, "permutation_summary"))
  permutation_diagnostics <- bind_rows(lapply(all_results, `[[`, "permutation_diagnostics"))

  cat("[4/4] Writing results tables...\n")
  fwrite(leave_summary, file.path(output_dir, "S41_leave_one_gene_summary.csv"))
  fwrite(leave_detail, file.path(output_dir, "S41_leave_one_gene_detail.csv"))
  fwrite(locus_results, file.path(output_dir, "S41_locus_exclusion_sensitivity.csv"))
  fwrite(permutation_summary,
         file.path(output_dir, "S42_matched_random_gene_set_permutation.csv"))
  fwrite(permutation_diagnostics,
         file.path(output_dir, "S42_permutation_diagnostics_first1000.csv"))

  cat("\n=== EOAD Internal Robustness Complete ===\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  print(sessionInfo())

  invisible(list(
    power = power_table,
    leave_summary = leave_summary,
    leave_detail = leave_detail,
    locus_results = locus_results,
    permutation_summary = permutation_summary,
    permutation_diagnostics = permutation_diagnostics
  ))
}


# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================
#
# twas_files <- c(
#   "Cortex" =
#     "data/SPrediXcan/EOAD_hg19__PM__mashr_Brain_Cortex.csv",
#   "Anterior cingulate cortex BA24" =
#     "data/SPrediXcan/EOAD_hg19__PM__mashr_Brain_Anterior_cingulate_cortex_BA24.csv",
#   "Putamen basal ganglia" =
#     "data/SPrediXcan/EOAD_hg19__PM__mashr_Brain_Putamen_basal_ganglia.csv",
#   "Frontal Cortex BA9" =
#     "data/SPrediXcan/EOAD_hg19__PM__mashr_Brain_Frontal_Cortex_BA9.csv",
#   "Hippocampus" =
#     "data/SPrediXcan/EOAD_hg19__PM__mashr_Brain_Hippocampus.csv"
# )
#
# results <- run_eoad_internal_robustness(
#   magma_gene_file = "data/MAGMA/Supplementary_Table_S1_EOAD_MAGMA_Complete.csv",
#   oligo_gene_file = "data/SPrediXcan/Supplementary_Table_2_Oligodendrocyte_Genes.csv",
#   twas_files = twas_files,
#   output_dir = "results/eoad_internal_robustness",
#   eoad_effective_n = 1573,
#   n_perm = 2000,
#   random_seed = 20260525
# )
