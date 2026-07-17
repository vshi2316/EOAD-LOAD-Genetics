# ==============================================================================
# Concordance and phenotype-level contextualization framework
# ==============================================================================
#
# Purpose: Separate the analytical roles of:
#            1. FinnGen EOAD discovery against the EADB AD/ADRD reference,
#            2. cross-method concordance with published EOAD genetic evidence,
#            3. phenotype-level contextualization in ADNI/HABS/A4/AIBL, and
#            4. targeted negative sensitivity analyses for the 44-gene
#               oligodendrocyte/myelin candidate set.
#
# The implementation avoids proprietary or locally wrapped packages and
# uses CRAN/base R only (data.table + stats) and can summarize outputs generated
# by open-source tools such as MAGMA, clusterProfiler, coloc and SMR/HEIDI.
#
# Terminology guardrails:
#   - "Discovery" is reserved for FinnGen EOAD analyses evaluated against the
#     broader EADB AD/ADRD GWAS, MAGMA, and GSEA reference profile.
#   - "Concordance" describes comparison with published EOAD evidence, including
#     Bradley et al. 2025 gene sets, unless formal sample independence is proven.
#   - "Phenotype-level contextualization" describes ADNI/HABS/A4/AIBL analyses.
#   - The targeted 44-gene analyses are sensitivity tests, not proof of a
#     concentrated oligodendrocyte/myelin causal chain.
#
# Dependencies: data.table
# Optional upstream tools represented by input files: MAGMA, clusterProfiler,
#   coloc, SMR/HEIDI.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)


# ==============================================================================
# 0. Shared Helpers
# ==============================================================================

theme_message <- function(...) {
  message(sprintf(...))
}

safe_fread <- function(file, label = basename(file), required = FALSE, ...) {
  if (is.null(file) || is.na(file) || !nzchar(file) || !file.exists(file)) {
    msg <- sprintf("Missing %s: %s", label, as.character(file))
    if (required) {
      stop(msg, call. = FALSE)
    }
    warning(msg, call. = FALSE)
    return(NULL)
  }
  fread(file, ...)
}

require_columns <- function(dt, cols, label) {
  missing <- setdiff(cols, names(dt))
  if (length(missing) > 0) {
    stop(
      label, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

first_existing_col <- function(dt, candidates, label) {
  hit <- intersect(candidates, names(dt))
  if (length(hit) == 0) {
    stop(
      label, " requires one of: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  hit[1]
}

num <- function(x) suppressWarnings(as.numeric(x))
int <- function(x) suppressWarnings(as.integer(gsub("^chr", "", as.character(x), ignore.case = TRUE)))

normalize_gene <- function(x) {
  out <- toupper(trimws(as.character(x)))
  out[out %in% c("", "NA", "NAN")] <- NA_character_
  out
}

make_bins <- function(x, n_bins = 10) {
  x <- num(x)
  if (!any(is.finite(x))) {
    return(rep(1L, length(x)))
  }
  x[!is.finite(x)] <- stats::median(x[is.finite(x)], na.rm = TRUE)
  breaks <- unique(as.numeric(stats::quantile(
    x,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE,
    type = 8
  )))
  if (length(breaks) <= 2) {
    return(rep(1L, length(x)))
  }
  as.integer(cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE))
}

write_status_table <- function(file, status, reason, ...) {
  fwrite(
    data.table(status = status, reason = reason, ...),
    file
  )
}


# ==============================================================================
# 1. Input Discovery and Configuration
# ==============================================================================

find_project_base <- function(start_dirs = c("H:/AD_*", "F:/AD_*", "data")) {
  candidates <- unique(unlist(lapply(start_dirs, Sys.glob), use.names = FALSE))
  candidates <- candidates[file.exists(file.path(
    candidates,
    "Stage2_Results",
    "Supplementary_Tables",
    "Supplementary_Table_S1_EOAD_MAGMA_Complete.csv"
  ))]
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  candidates[1]
}

make_default_config <- function(base_dir = find_project_base()) {
  if (is.na(base_dir)) {
    base_dir <- "."
  }

  list(
    base_dir = base_dir,
    output_dir = file.path(base_dir, "Revised_AandD_concordance_contextualization"),
    input_dir = file.path(base_dir, "Revised_AandD_concordance_inputs"),
    eoad_magma_file = file.path(
      base_dir,
      "Stage2_Results",
      "Supplementary_Tables",
      "Supplementary_Table_S1_EOAD_MAGMA_Complete.csv"
    ),
    eadb_magma_file = file.path(
      base_dir,
      "Stage2_Results",
      "Supplementary_Tables",
      "Supplementary_Table_S2_EADB_MAGMA_Complete.csv"
    ),
    gsea_dir = file.path(base_dir, "Stage2_Results", "NC_Tables"),
    gsea_file_template = "Table_%s_GSEA_%s.csv",
    literature_gene_set_file = file.path(
      base_dir,
      "Revised_AandD_concordance_inputs",
      "published_EOAD_gene_sets.csv"
    ),
    phenotype_summary_file = file.path(
      base_dir,
      "Revised_AandD_concordance_inputs",
      "phenotype_level_contextualization.csv"
    ),
    targeted_magma_file = file.path(
      base_dir,
      "SMR",
      "Targeted_MAGMA_candidate_gene_set_support",
      "Supplementary_Table_candidate_gene_set_MAGMA_matched_permutation.csv"
    ),
    coloc_summary_file = file.path(
      base_dir,
      "SMR",
      "Targeted_sc_eQTL_coloc_ABF_AandD_ready",
      "Supplementary_Table_sc_eQTL_coloc_ABF_summary_by_cell_type.csv"
    ),
    coloc_detail_file = file.path(
      base_dir,
      "SMR",
      "Targeted_sc_eQTL_coloc_ABF_AandD_ready",
      "Supplementary_Table_sc_eQTL_coloc_ABF_complete_results.csv"
    ),
    smr_file = file.path(
      base_dir,
      "SMR",
      "Targeted_SMR_HEIDI_coloc_susie_results",
      "Supplementary_Table_Targeted_SMR_HEIDI_primary_p5e-8.csv"
    ),
    n_perm = 10000L,
    random_seed = 20260530L
  )
}


# ==============================================================================
# 2. MAGMA Standardization
# ==============================================================================

standardize_magma <- function(dt, trait) {
  if (is.null(dt)) {
    stop(trait, " MAGMA input is NULL.", call. = FALSE)
  }
  setDT(dt)

  gene_col <- first_existing_col(dt, c("Gene_Symbol", "SYMBOL", "GENE", "gene"), paste0(trait, " MAGMA"))
  chr_col <- first_existing_col(dt, c("CHR", "chr", "chromosome"), paste0(trait, " MAGMA"))
  start_col <- first_existing_col(dt, c("START", "start", "BP", "pos"), paste0(trait, " MAGMA"))
  stop_col <- first_existing_col(dt, c("STOP", "stop", "END", "end"), paste0(trait, " MAGMA"))
  nsnps_col <- first_existing_col(dt, c("NSNPS", "n_snps", "N_SNPS"), paste0(trait, " MAGMA"))
  p_col <- first_existing_col(dt, c("P", "PVALUE", "p", "p_value"), paste0(trait, " MAGMA"))

  out <- copy(dt)
  out[, trait := trait]
  out[, Gene_Symbol := as.character(.SD[[1]]), .SDcols = gene_col]
  out[, gene_upper := normalize_gene(Gene_Symbol)]
  out[, CHR := int(.SD[[1]]), .SDcols = chr_col]
  out[, START := int(.SD[[1]]), .SDcols = start_col]
  out[, STOP := int(.SD[[1]]), .SDcols = stop_col]
  out[, NSNPS := int(.SD[[1]]), .SDcols = nsnps_col]
  out[, P := num(.SD[[1]]), .SDcols = p_col]

  if ("ZSTAT" %in% names(out)) {
    out[, ZSTAT := num(ZSTAT)]
  } else {
    out[, ZSTAT := stats::qnorm(1 - pmin(pmax(P, .Machine$double.xmin), 1) / 2)]
  }

  out <- out[
    !is.na(gene_upper) &
      is.finite(CHR) & is.finite(START) & is.finite(STOP) &
      is.finite(NSNPS) & is.finite(P) & is.finite(ZSTAT) &
      P > 0 & P <= 1
  ]

  setorder(out, P)
  out <- out[!duplicated(gene_upper)]
  out[, gene_length := pmax(1L, STOP - START + 1L)]
  out[, neglog10p := -log10(pmax(P, .Machine$double.xmin))]
  out[, MAGMA_FDR := p.adjust(P, method = "BH")]
  out[, MAGMA_rank_P := frank(P, ties.method = "min")]
  out[, MAGMA_rank_percentile := MAGMA_rank_P / .N]
  out[, .row_id := as.character(seq_len(.N))]
  out[, NSNPS_bin := make_bins(NSNPS, 10)]
  out[, gene_len_bin := make_bins(gene_length, 10)]
  out[]
}

read_magma_pair <- function(config) {
  eoad <- safe_fread(config$eoad_magma_file, "EOAD MAGMA", required = TRUE)
  eadb <- safe_fread(config$eadb_magma_file, "EADB MAGMA", required = TRUE)
  list(
    EOAD = standardize_magma(eoad, "EOAD"),
    EADB = standardize_magma(eadb, "EADB")
  )
}

write_magma_overview <- function(magma, output_dir) {
  overview <- rbindlist(lapply(names(magma), function(trait) {
    magma[[trait]][, .(
      trait = trait,
      n_genes = .N,
      n_nominal_p05 = sum(P < 0.05, na.rm = TRUE),
      n_fdr05 = sum(MAGMA_FDR < 0.05, na.rm = TRUE),
      min_p = min(P, na.rm = TRUE),
      top_gene = Gene_Symbol[which.min(P)]
    )]
  }))
  fwrite(
    overview,
    file.path(output_dir, "Supplementary_Table_44_MAGMA_gene_level_overview.csv")
  )
  overview
}


# ==============================================================================
# 3. Unbiased GSEA Discovery Summary
# ==============================================================================

read_gsea_results <- function(config, traits = c("EOAD", "EADB"),
                              databases = c("GO_BP", "GO_MF", "KEGG")) {
  manifest <- CJ(trait = traits, database = databases)
  manifest[, file := file.path(
    config$gsea_dir,
    sprintf(config$gsea_file_template, trait, database)
  )]

  out <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
    z <- safe_fread(manifest$file[i], paste(manifest$trait[i], manifest$database[i]))
    if (is.null(z)) {
      return(NULL)
    }
    z[, trait := manifest$trait[i]]
    z[, database := manifest$database[i]]
    z
  }), fill = TRUE)

  if (nrow(out) == 0) {
    return(NULL)
  }
  require_columns(out, c("Description", "p.adjust"), "GSEA results")
  out[, p.adjust_num := num(p.adjust)]
  out[, pvalue_num := if ("pvalue" %in% names(out)) num(pvalue) else NA_real_]
  out[, significant_fdr05 := p.adjust_num < 0.05]
  setorder(out, trait, database, p.adjust_num)
  out[]
}

keyword_summary <- function(gsea, keyword_dt) {
  rbindlist(lapply(seq_len(nrow(keyword_dt)), function(i) {
    z <- gsea[grepl(keyword_dt$pattern[i], Description, ignore.case = TRUE)]
    if (nrow(z) == 0) {
      return(data.table(
        keyword_group = keyword_dt$group[i],
        trait = unique(gsea$trait),
        database = unique(gsea$database),
        n_matching_terms = 0L,
        n_fdr05_terms = 0L,
        min_fdr = NA_real_,
        top_term = NA_character_,
        top_fdr = NA_real_
      ))
    }
    setorder(z, p.adjust_num)
    data.table(
      keyword_group = keyword_dt$group[i],
      trait = unique(gsea$trait),
      database = unique(gsea$database),
      n_matching_terms = nrow(z),
      n_fdr05_terms = sum(z$p.adjust_num < 0.05, na.rm = TRUE),
      min_fdr = min(z$p.adjust_num, na.rm = TRUE),
      top_term = z$Description[1],
      top_fdr = z$p.adjust_num[1]
    )
  }), fill = TRUE)
}

write_gsea_summary <- function(gsea_all, output_dir) {
  if (is.null(gsea_all) || nrow(gsea_all) == 0) {
    write_status_table(
      file.path(output_dir, "Supplementary_Table_45_GSEA_summary_NOT_RUN.csv"),
      "not_run",
      "No GSEA files were detected."
    )
    return(NULL)
  }

  overview <- gsea_all[, .(
    n_terms = .N,
    n_fdr05 = sum(significant_fdr05, na.rm = TRUE),
    min_fdr = min(p.adjust_num, na.rm = TRUE),
    top_term = Description[which.min(p.adjust_num)]
  ), by = .(trait, database)]

  top30 <- gsea_all[, head(.SD, 30), by = .(trait, database)]
  top30[, panel := "top30"]
  overview[, panel := "overview"]

  keyword_groups <- data.table(
    group = c(
      "glutamate_synaptic",
      "immune_microglia",
      "intracellular_signaling",
      "amyloid_APP_clearance",
      "cell_adhesion",
      "lipid_cholesterol",
      "white_matter_myelin"
    ),
    pattern = c(
      "glutamate|glutamatergic|synap|neurotransmitter",
      "immune|microglia|leukocyte|lymphocyte|inflamm|cytokine|complement|phagocyt|T cell|B cell|natural killer",
      "signaling|phosphorylation|kinase|MAPK|GTPase|intracellular|second messenger",
      "amyloid|amyloid-beta|APP|amyloid precursor|beta clearance",
      "adhesion|cell-cell|plasma-membrane adhesion",
      "lipid|cholesterol|sterol|lipoprotein|phospholipid",
      "myelin|oligodendrocyte|axon ensheath|white matter|glial"
    )
  )
  keyword <- rbindlist(lapply(
    split(gsea_all, by = c("trait", "database"), keep.by = TRUE),
    keyword_summary,
    keyword_dt = keyword_groups
  ), fill = TRUE)

  fwrite(overview, file.path(output_dir, "Supplementary_Table_45_unbiased_GSEA_summary.csv"))
  fwrite(overview, file.path(output_dir, "Supplementary_Table_45A_unbiased_GSEA_overview.csv"))
  fwrite(top30, file.path(output_dir, "Supplementary_Table_45B_unbiased_GSEA_top30.csv"))
  fwrite(keyword, file.path(output_dir, "Supplementary_Table_45C_GSEA_keyword_concordance.csv"))

  list(overview = overview, top30 = top30, keyword = keyword)
}


# ==============================================================================
# 4. Published EOAD Gene-Set Concordance
# ==============================================================================

write_literature_template <- function(config) {
  dir.create(dirname(config$literature_gene_set_file), recursive = TRUE, showWarnings = FALSE)
  template <- sub("\\.csv$", "_TEMPLATE.csv", config$literature_gene_set_file)
  if (!file.exists(template)) {
    fwrite(
      data.table(
        source = c("Bradley2025_AandD", "Bradley2025_AandD", "Bradley2025_AandD"),
        set_name = c("S10_prioritized_genes", "MAGMA_evidence_genes", "MetaBrain_colocalization_genes"),
        gene = c("ADD_HGNC_SYMBOL", "ADD_HGNC_SYMBOL", "ADD_HGNC_SYMBOL"),
        notes = c(
          "Replace with curated genes from the published paper or supplement.",
          "The file supports cross-method concordance without constituting subtype-specific replication.",
          "Use direct HGNC symbols; one gene per row."
        )
      ),
      template
    )
  }
  template
}

compute_set_metrics <- function(x) {
  p <- pmax(num(x$P), .Machine$double.xmin)
  data.table(
    n_genes = nrow(x),
    mean_z = mean(num(x$ZSTAT), na.rm = TRUE),
    median_z = stats::median(num(x$ZSTAT), na.rm = TRUE),
    mean_neglog10p = mean(-log10(p), na.rm = TRUE),
    median_neglog10p = stats::median(-log10(p), na.rm = TRUE),
    min_p = min(p, na.rm = TRUE),
    n_nominal_p05 = sum(p < 0.05, na.rm = TRUE),
    prop_nominal_p05 = mean(p < 0.05, na.rm = TRUE),
    mean_rank_percentile = mean(num(x$MAGMA_rank_percentile), na.rm = TRUE)
  )
}

run_matched_rank_enrichment <- function(magma, genes, source, set_name,
                                        sensitivity = "all_regions",
                                        n_perm = 10000L, seed = 1L) {
  genes <- unique(normalize_gene(genes))
  genes <- genes[!is.na(genes)]

  cand <- magma[gene_upper %in% genes]
  bg <- magma[!gene_upper %in% genes]

  if (nrow(cand) < 3 || nrow(bg) < 100) {
    return(data.table(
      sensitivity = sensitivity,
      source = source,
      set_name = set_name,
      trait = unique(magma$trait)[1],
      status = "too_few_genes_in_MAGMA",
      n_input_genes = length(genes),
      n_in_MAGMA = nrow(cand),
      n_background = nrow(bg)
    ))
  }

  obs <- compute_set_metrics(cand)

  pools <- vector("list", nrow(cand))
  for (i in seq_len(nrow(cand))) {
    pool <- bg[NSNPS_bin == cand$NSNPS_bin[i] & gene_len_bin == cand$gene_len_bin[i]]
    if (nrow(pool) == 0) {
      pool <- bg[NSNPS_bin == cand$NSNPS_bin[i]]
    }
    if (nrow(pool) == 0) {
      pool <- bg
    }
    pools[[i]] <- pool$.row_id
  }

  set.seed(seed)
  setkey(bg, .row_id)
  null <- vector("list", n_perm)
  for (b in seq_len(n_perm)) {
    sampled_ids <- vapply(pools, function(pool) sample(pool, 1), character(1))
    null[[b]] <- compute_set_metrics(bg[sampled_ids])
  }
  null <- rbindlist(null)

  wilcox_z <- suppressWarnings(stats::wilcox.test(
    cand$ZSTAT,
    bg$ZSTAT,
    alternative = "greater",
    exact = FALSE
  ))
  wilcox_rank <- suppressWarnings(stats::wilcox.test(
    cand$MAGMA_rank_percentile,
    bg$MAGMA_rank_percentile,
    alternative = "less",
    exact = FALSE
  ))

  cbind(
    data.table(
      sensitivity = sensitivity,
      source = source,
      set_name = set_name,
      trait = unique(magma$trait)[1],
      status = "completed",
      n_input_genes = length(genes),
      n_in_MAGMA = nrow(cand),
      n_background = nrow(bg),
      matched_genes = paste(cand$Gene_Symbol, collapse = ";")
    ),
    obs,
    data.table(
      empirical_p_mean_z = (1 + sum(null$mean_z >= obs$mean_z)) / (n_perm + 1),
      empirical_p_mean_neglog10p =
        (1 + sum(null$mean_neglog10p >= obs$mean_neglog10p)) / (n_perm + 1),
      empirical_p_min_p = (1 + sum(null$min_p <= obs$min_p)) / (n_perm + 1),
      empirical_p_rank =
        (1 + sum(null$mean_rank_percentile <= obs$mean_rank_percentile)) /
        (n_perm + 1),
      wilcox_z_greater_p = wilcox_z$p.value,
      wilcox_rank_better_p = wilcox_rank$p.value,
      n_permutations = n_perm
    )
  )
}

filter_magma_by_sensitivity <- function(magma, sensitivity) {
  if (sensitivity == "all_regions") {
    return(magma)
  }
  if (sensitivity == "exclude_chr19_APOE_region") {
    return(magma[!(CHR == 19 & START <= 46000000 & STOP >= 44000000)])
  }
  if (sensitivity == "exclude_chr19_all") {
    return(magma[CHR != 19])
  }
  stop("Unknown sensitivity: ", sensitivity, call. = FALSE)
}

run_published_gene_set_concordance <- function(magma, config) {
  template <- write_literature_template(config)
  literature_sets <- safe_fread(config$literature_gene_set_file, "published EOAD gene sets")
  if (is.null(literature_sets)) {
    return(data.table(
      status = "not_run",
      reason = paste0("No curated published EOAD gene-set file found. Template: ", template),
      required_columns = "source,set_name,gene",
      terminology = "cross_method_concordance"
    ))
  }

  require_columns(literature_sets, c("source", "set_name", "gene"), "Published EOAD gene sets")
  literature_sets[, gene_upper := normalize_gene(gene)]
  literature_sets <- literature_sets[
    !is.na(gene_upper) &
      !gene_upper %in% c("ADD_HGNC_SYMBOL", "REPLACE_WITH_GENE")
  ]
  if (nrow(literature_sets) == 0) {
    return(data.table(
      status = "not_run",
      reason = "Published EOAD gene-set file contained no usable gene symbols.",
      required_columns = "source,set_name,gene",
      terminology = "cross_method_concordance"
    ))
  }

  setkey(literature_sets, source, set_name)
  set_defs <- unique(literature_sets[, .(source, set_name)])
  sensitivities <- c("all_regions", "exclude_chr19_APOE_region", "exclude_chr19_all")

  out <- rbindlist(lapply(seq_len(nrow(set_defs)), function(i) {
    src <- set_defs$source[i]
    set_nm <- set_defs$set_name[i]
    genes <- literature_sets[.(src, set_nm), gene_upper]

    rbindlist(lapply(sensitivities, function(sens) {
      rbindlist(lapply(names(magma), function(trait) {
        magma_sens <- filter_magma_by_sensitivity(magma[[trait]], sens)
        run_matched_rank_enrichment(
          magma = magma_sens,
          genes = genes,
          source = src,
          set_name = set_nm,
          sensitivity = sens,
          n_perm = config$n_perm,
          seed = config$random_seed + i * 1000 +
            match(trait, names(magma)) * 10 + match(sens, sensitivities)
        )
      }), fill = TRUE)
    }), fill = TRUE)
  }), fill = TRUE)

  for (pcol in c("empirical_p_mean_z", "empirical_p_rank", "wilcox_z_greater_p")) {
    fdr_col <- sub("^", "FDR_", pcol)
    out[status == "completed", (fdr_col) := p.adjust(.SD[[1]], method = "BH"),
        by = .(trait, sensitivity), .SDcols = pcol]
  }

  out[]
}

write_concordance_outputs <- function(concordance, output_dir) {
  if (is.null(concordance) || nrow(concordance) == 0 ||
      all(concordance$status == "not_run")) {
    fwrite(
      concordance,
      file.path(output_dir, "Supplementary_Table_46_Bradley_gene_set_concordance_NOT_RUN.csv")
    )
    fwrite(
      data.table(
        status = "not_run",
        reason = "Bradley concordance was not run, so APOE/chromosome 19 sensitivity could not be evaluated.",
        terminology = "cross_method_concordance_sensitivity"
      ),
      file.path(output_dir, "Supplementary_Table_47_APOE_chr19_concordance_sensitivity_NOT_RUN.csv")
    )
    return(invisible(NULL))
  }

  all_regions <- concordance[sensitivity == "all_regions"]
  sensitivity <- concordance[sensitivity != "all_regions"]

  fwrite(
    all_regions,
    file.path(output_dir, "Supplementary_Table_46_Bradley_gene_set_concordance.csv")
  )
  fwrite(
    sensitivity,
    file.path(output_dir, "Supplementary_Table_47_APOE_chr19_concordance_sensitivity.csv")
  )
}

write_concordance_not_run_outputs <- function(output_dir, reason) {
  not_run <- data.table(
    status = "not_run",
    reason = reason,
    terminology = "cross_method_concordance"
  )
  fwrite(
    not_run,
    file.path(output_dir, "Supplementary_Table_46_Bradley_gene_set_concordance_NOT_RUN.csv")
  )
  fwrite(
    data.table(
      status = "not_run",
      reason = paste(
        reason,
        "APOE/chromosome 19 sensitivity requires completed concordance results."
      ),
      terminology = "cross_method_concordance_sensitivity"
    ),
    file.path(output_dir, "Supplementary_Table_47_APOE_chr19_concordance_sensitivity_NOT_RUN.csv")
  )
  invisible(not_run)
}


# ==============================================================================
# 5. Targeted 44-Gene Sensitivity Summary
# ==============================================================================

summarize_targeted_sensitivity <- function(config, output_dir) {
  targeted_summary <- list()

  targeted_magma <- safe_fread(config$targeted_magma_file, "targeted 44-gene MAGMA")
  if (!is.null(targeted_magma)) {
    setDT(targeted_magma)
    targeted_summary[["magma"]] <- targeted_magma[, .(
      evidence_layer = "Targeted 44-gene MAGMA/permutation",
      test = paste(set_name, sensitivity, sep = " | "),
      n_genes = n_genes,
      primary_metric = "empirical_p_mean_z",
      primary_value = num(empirical_p_mean_z),
      secondary_metric = "empirical_p_mean_neglog10p",
      secondary_value = num(empirical_p_mean_neglog10p),
      interpretation = ifelse(
        num(empirical_p_mean_z) < 0.05 |
          num(empirical_p_mean_neglog10p) < 0.05,
        "nominal_positive_sensitivity_signal",
        "not_supportive"
      )
    )]
  }

  coloc_summary <- safe_fread(config$coloc_summary_file, "single-cell eQTL coloc summary")
  if (!is.null(coloc_summary)) {
    setDT(coloc_summary)
    targeted_summary[["coloc"]] <- coloc_summary[, .(
      evidence_layer = "Targeted single-cell eQTL coloc.abf",
      test = as.character(cell_type),
      n_genes = genes_tested,
      primary_metric = "max_PP_H4",
      primary_value = num(max_PP_H4),
      secondary_metric = "n_PP_H4_ge_075",
      secondary_value = num(n_PP_H4_ge_075),
      interpretation = ifelse(num(max_PP_H4) >= 0.75, "supportive", "not_supportive")
    )]
  }

  smr_res <- safe_fread(config$smr_file, "targeted SMR/HEIDI")
  if (!is.null(smr_res)) {
    setDT(smr_res)
    completed <- smr_res[is.finite(num(p_SMR))]
    targeted_summary[["smr"]] <- data.table(
      evidence_layer = "Targeted SMR/HEIDI primary threshold",
      test = "all_available_tests",
      n_genes = if (nrow(completed) > 0 && "gene_tested" %in% names(completed)) {
        uniqueN(completed$gene_tested)
      } else {
        0L
      },
      primary_metric = "min_p_SMR",
      primary_value = if (nrow(completed) > 0) min(num(completed$p_SMR), na.rm = TRUE) else NA_real_,
      secondary_metric = "n_p_SMR_lt_0.05",
      secondary_value = if (nrow(completed) > 0) sum(num(completed$p_SMR) < 0.05, na.rm = TRUE) else 0,
      interpretation = if (nrow(completed) > 0 &&
                           any(num(completed$p_SMR) < 0.05, na.rm = TRUE)) {
        "nominal_positive_sensitivity_signal"
      } else {
        "not_supportive"
      }
    )
  }

  out <- rbindlist(targeted_summary, fill = TRUE)
  if (nrow(out) == 0) {
    out <- data.table(
      evidence_layer = "Targeted 44-gene sensitivity",
      status = "not_found",
      interpretation = "No targeted MAGMA, coloc or SMR result files were detected."
    )
  }

  fwrite(
    out,
    file.path(output_dir, "Supplementary_Table_48_targeted_44_gene_sensitivity_summary.csv")
  )

  coloc_detail <- safe_fread(config$coloc_detail_file, "single-cell eQTL coloc detail")
  if (!is.null(coloc_detail)) {
    fwrite(
      coloc_detail,
      file.path(output_dir, "Supplementary_Table_49_sc_eQTL_coloc_detail.csv")
    )
  } else {
    write_status_table(
      file.path(output_dir, "Supplementary_Table_49_sc_eQTL_coloc_detail_NOT_RUN.csv"),
      "not_run",
      "No single-cell eQTL colocalization detail file was detected."
    )
  }

  out
}


# ==============================================================================
# 6. Orthogonal Phenotype-Level Contextualization
# ==============================================================================

write_phenotype_template <- function(config) {
  dir.create(dirname(config$phenotype_summary_file), recursive = TRUE, showWarnings = FALSE)
  template <- sub("\\.csv$", "_TEMPLATE.csv", config$phenotype_summary_file)
  if (!file.exists(template)) {
    fwrite(
      data.table(
        cohort = c("ADNI", "HABS", "A4", "AIBL"),
        modality = c(
          "genotype_imaging",
          "WMH_tau_cognition",
          "amyloid_enriched_WMH_cognition",
          "conversion_benchmark"
        ),
        phenotype = c(
          "Oligodendrocyte_PRS_x_age_on_WMH",
          "WMH_to_pTau217_to_cognition",
          "WMH_to_PACC",
          "APOE4_to_conversion"
        ),
        beta_or_effect = NA_real_,
        se = NA_real_,
        p_value = NA_real_,
        n = NA_integer_,
        interpretation = c(
          "Fill from existing ADNI result table.",
          "Fill from existing HABS association or mediation-style table.",
          "Fill from existing A4 result table.",
          "Fill from existing AIBL survival or logistic table."
        ),
        circularity_note = "Phenotype-level contextualization; not direct gene-level validation."
      ),
      template
    )
  }
  template
}

write_phenotype_context <- function(config, output_dir) {
  template <- write_phenotype_template(config)
  phenotype <- safe_fread(config$phenotype_summary_file, "phenotype-level contextualization")
  if (is.null(phenotype)) {
    write_status_table(
      file.path(output_dir, "Supplementary_Table_50_phenotype_contextualization_NOT_RUN.csv"),
      "not_run",
      paste0("No phenotype-level contextualization file found. Template: ", template),
      terminology = "phenotype_level_contextualization"
    )
    return(NULL)
  }
  fwrite(
    phenotype,
    file.path(output_dir, "Supplementary_Table_50_phenotype_level_contextualization.csv")
  )
  phenotype
}


# ==============================================================================
# 7. Language and Provenance Notes
# ==============================================================================

write_language_notes <- function(output_dir) {
  notes <- c(
    "Recommended manuscript framing",
    "",
    "The revised framework characterizes EOAD/EADB heterogeneity through",
    "unbiased pathway discovery, cross-method concordance with published EOAD",
    "genetic evidence, phenotype-level contextualization in external cohorts,",
    "and targeted negative sensitivity analyses for the 44-gene",
    "oligodendrocyte/myelin candidate set.",
    "",
    "Preferred terms:",
    "- discovery: current EOAD/EADB GWAS, MAGMA and GSEA outputs.",
    "- concordance: overlap or enrichment relative to published EOAD evidence.",
    "- phenotype-level contextualization: ADNI/HABS/A4/AIBL clinical and imaging results.",
    "- targeted sensitivity: 44-gene MAGMA, coloc and SMR/HEIDI analyses.",
    "",
    "Avoid:",
    "- clinical deployment claims.",
    "- subtype-specific genetic replication claims unless sample independence is proven.",
    "- causal oligodendrocyte/myelin claims from negative targeted sensitivity analyses.",
    "",
    "Conservative conclusion:",
    "EOAD and EADB show strong global genetic overlap with pathway-level",
    "heterogeneity. EADB shows broader immune and microglial enrichment, whereas",
    "EOAD shows a narrower architecture involving adhesion, amyloid-clearance,",
    "selected immune, lipid and white-matter-related signals. The targeted",
    "44-gene oligodendrocyte/myelin model remains exploratory and was not",
    "supported as a concentrated causal chain by the targeted sensitivity layers."
  )
  writeLines(
    notes,
    con = file.path(output_dir, "manuscript_language_guardrails.txt"),
    useBytes = TRUE
  )
}

write_supplementary_table_presence_check <- function(output_dir) {
  expected <- data.table(
    table_number = 44:50,
    expected_layer = c(
      "MAGMA gene-level overview",
      "Unbiased GSEA discovery summary",
      "Published EOAD gene-set concordance",
      "APOE/chromosome 19 concordance sensitivity",
      "Targeted 44-gene sensitivity summary",
      "Single-cell eQTL colocalization detail",
      "Phenotype-level contextualization"
    ),
    pattern = c(
      "^Supplementary_Table_44.*\\.csv$",
      "^Supplementary_Table_45.*\\.csv$",
      "^Supplementary_Table_46.*\\.csv$",
      "^Supplementary_Table_47.*\\.csv$",
      "^Supplementary_Table_48.*\\.csv$",
      "^Supplementary_Table_49.*\\.csv$",
      "^Supplementary_Table_50.*\\.csv$"
    )
  )

  files <- list.files(output_dir, pattern = "\\.csv$", full.names = FALSE)
  expected[, matched_files := vapply(
    pattern,
    function(x) paste(files[grepl(x, files)], collapse = ";"),
    character(1)
  )]
  expected[, present := nzchar(matched_files)]
  expected[, status := ifelse(present, "present", "missing")]

  fwrite(
    expected,
    file.path(output_dir, "Supplementary_Table_44_50_presence_check.csv")
  )
  expected[]
}


# ==============================================================================
# 8. Master Runner
# ==============================================================================

run_revised_aandd_concordance_contextualization <- function(config = make_default_config()) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$input_dir, recursive = TRUE, showWarnings = FALSE)

  theme_message("Base directory: %s", config$base_dir)
  theme_message("Output directory: %s", config$output_dir)
  theme_message("Permutations per concordance test: %s", config$n_perm)

  magma <- tryCatch(
    read_magma_pair(config),
    error = function(e) {
      warning(conditionMessage(e), call. = FALSE)
      NULL
    }
  )

  if (is.null(magma)) {
    magma_overview <- write_status_table(
      file.path(config$output_dir, "Supplementary_Table_44_MAGMA_gene_level_overview_NOT_RUN.csv"),
      "not_run",
      "EOAD and/or EADB MAGMA input files were not available; MAGMA overview could not be generated."
    )
  } else {
    magma_overview <- write_magma_overview(magma, config$output_dir)
  }

  gsea_all <- read_gsea_results(config)
  gsea_summary <- write_gsea_summary(gsea_all, config$output_dir)

  if (is.null(magma)) {
    concordance <- write_concordance_not_run_outputs(
      config$output_dir,
      "EOAD/EADB MAGMA inputs were unavailable; published gene-set concordance could not be evaluated."
    )
  } else {
    concordance <- run_published_gene_set_concordance(magma, config)
    write_concordance_outputs(concordance, config$output_dir)
  }

  targeted <- summarize_targeted_sensitivity(config, config$output_dir)
  phenotype <- write_phenotype_context(config, config$output_dir)
  write_language_notes(config$output_dir)
  table_presence <- write_supplementary_table_presence_check(config$output_dir)

  summary_file <- file.path(config$output_dir, "analysis_summary.txt")
  sink(summary_file)
  cat("=======================================================================\n")
  cat("EOAD/EADB concordance and phenotype-level contextualization framework\n")
  cat("=======================================================================\n")
  cat("Date:", as.character(Sys.Date()), "\n")
  cat("Base directory:", config$base_dir, "\n")
  cat("Output directory:", config$output_dir, "\n")
  cat("Permutations per gene-set concordance test:", config$n_perm, "\n\n")
  cat("Core logic:\n")
  cat("  1. Discovery uses current EOAD/EADB GWAS, MAGMA and GSEA outputs.\n")
  cat("  2. Bradley 2025 comparisons are cross-method concordance analyses.\n")
  cat("  3. ADNI/HABS/A4/AIBL layers are phenotype-level contextualization.\n")
  cat("  4. Targeted 44-gene analyses are sensitivity tests, not causal proof.\n\n")
  cat("Main output files:\n")
  print(data.table(
    file = list.files(config$output_dir, full.names = FALSE),
    full_path = file.path(config$output_dir, list.files(config$output_dir, full.names = FALSE))
  ))
  cat("=======================================================================\n")
  sink()

  list(
    config = config,
    magma_overview = magma_overview,
    gsea_summary = gsea_summary,
    concordance = concordance,
    targeted_sensitivity = targeted,
    phenotype_context = phenotype,
    table_presence = table_presence,
    summary_file = summary_file
  )
}


# ==============================================================================
# 9. Example
# ==============================================================================
# Uncomment and edit paths as needed:
#
# config <- make_default_config(base_dir = "H:/AD_project_directory")
# config$literature_gene_set_file <- "data/published_EOAD_gene_sets.csv"
# config$phenotype_summary_file <- "data/phenotype_level_contextualization.csv"
# results <- run_revised_aandd_concordance_contextualization(config)
#
# To run with automatically detected local paths:
# results <- run_revised_aandd_concordance_contextualization()
# ==============================================================================
