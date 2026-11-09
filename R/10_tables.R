# =============================================================================
# 10_tables.R —— 论文表格：Table 1（队列基线特征）与补充表汇总（Tables_all.xlsx，S1–S10）
# =============================================================================
source("R/00_config.R")
need(c("survival", "openxlsx"))
`%||%` <- function(a, b) if (is.null(a)) b else a
STEP <- "10_tables"; OUT <- res_dir(STEP)

ALL_RDS <- sub("\\.rds$", "", list.files(DIR_PRO, pattern = "\\.rds$"))
META_ALL <- data.frame(
  cohort   = c("TCGA_KIRC", "GSE29609", "GSE22541", "GSE167573", "GSE53757", "GSE73731", "GSE67501", "IMvigor210", "GSE176307"),
  label    = c("TCGA-KIRC", "GSE29609", "GSE22541", "GSE167573", "GSE53757", "GSE73731", "GSE67501", "IMvigor210", "GSE176307"),
  platform = c("RNA-seq (STAR, log2 TPM+1)", "Agilent Whole Human Genome (GPL1708)", "Affymetrix HG-U133 Plus 2.0 (GPL570)",
               "Transcriptome array / RNA-seq", "Affymetrix HG-U133 Plus 2.0 (GPL570)", "Affymetrix HG-U133 Plus 2.0 (GPL570)",
               "Illumina HumanHT-12 WG-DASL (GPL14951)", "RNA-seq (log2 TPM+1)", "RNA-seq (log2 normalized)"),
  setting  = c("Primary ccRCC (72 paired normals)", "Primary ccRCC, nephrectomy", "Primary & metastatic ccRCC",
               "Renal cell carcinoma cohort", "72 paired ccRCC tumor and normal tissues",
               "Primary ccRCC across grades and stages", "Metastatic RCC treated with Nivolumab (anti-PD-1)",
               "Advanced urothelial/renal-tract carcinoma, Atezolizumab", "Advanced carcinoma, anti-PD-1/PD-L1 blockade"),
  role     = c("Training", "Validation (OS)", "Validation (DFS)", "Validation", "Validation (Diagnostic & Stage)",
               "Validation (Grade & Stage)", "Validation (RCC anti-PD-1)", "Validation (anti-PD-L1)", "Validation (anti-PD-1/PD-L1)"),
  endpoint = c("OS, DSS, PFI", "OS", "DFS", "OS", "Tumor vs Normal, Stage", "Fuhrman Grade, Stage",
               "Objective response (ORR)", "OS, ORR", "OS, PFS, ORR"),
  stringsAsFactors = FALSE
)
META <- META_ALL[META_ALL$cohort %in% ALL_RDS, ]

fmt_med <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[!is.na(x)]
  if (!length(x)) return("NA")
  q <- stats::quantile(x, c(0.5, 0.25, 0.75))
  sprintf("%.0f (%.0f-%.0f)", q[1], q[2], q[3])
}
fmt_n <- function(k, n) if (is.na(k) || n == 0) "NA" else sprintf("%d (%.0f%%)", k, 100 * k / n)

rows <- lapply(seq_len(nrow(META)), function(i) {
  n <- META$cohort[i]
  cl <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$clin
  if (n == "TCGA_KIRC") cl <- cl[cl$group == "Tumor", ]
  N <- nrow(cl)
  male <- if ("gender" %in% names(cl)) {
    g <- tolower(substr(as.character(cl$gender), 1, 1))
    fmt_n(sum(g == "m", na.rm = TRUE), sum(!is.na(g)))
  } else "NA"
  age <- if ("age" %in% names(cl)) fmt_med(cl$age) else "NA"
  adv_st <- if ("stage" %in% names(cl)) {
    s <- toupper(as.character(cl$stage))
    ok <- !is.na(s) & nzchar(s)
    fmt_n(sum(grepl("III|IV|3|4", s[ok])), sum(ok))
  } else "NA"
  hg <- if ("grade" %in% names(cl)) {
    g <- toupper(as.character(cl$grade))
    ok <- !is.na(g) & nzchar(g) & !g %in% c("GX", "NA")
    fmt_n(sum(grepl("3|4|HIGH", g[ok])), sum(ok))
  } else "NA"
  has_os <- !is.na(cl$OS.time) & !is.na(cl$OS) & cl$OS.time > 0
  ev <- if (sum(has_os) > 0) as.character(sum(cl$OS[has_os] == 1)) else "NA"
  fu <- if (sum(has_os) >= 10) {
    km <- tryCatch(survival::survfit(survival::Surv(OS.time, 1 - OS) ~ 1, data = cl[has_os, ]), error = function(e) NULL)
    med_d <- if (!is.null(km)) as.numeric(summary(km)$table["median"]) else NA_real_
    if (is.na(med_d)) med_d <- stats::median(cl$OS.time[has_os], na.rm = TRUE)
    sprintf("%.1f", med_d / 30.4375)
  } else "NA"
  resp <- if ("response" %in% names(cl)) {
    r <- cl$response
    fmt_n(sum(r %in% c("CR/PR", "Responder", "R"), na.rm = TRUE), sum(r %in% c("CR/PR", "SD/PD", "Responder", "Non-responder"), na.rm = TRUE))
  } else "NA"

  data.frame(
    Cohort = META$label[i], Platform = META$platform[i], Setting = META$setting[i], Role = META$role[i],
    Patients = N, `Age, median (IQR)` = age, `Male, n (%)` = male,
    `Stage III-IV, n (%)` = adv_st, `Grade G3-4, n (%)` = hg,
    `End point` = META$endpoint[i], Events = ev, `Median follow-up, months` = fu,
    `Responders, n (%)` = resp, check.names = FALSE
  )
})
tab1 <- do.call(rbind, rows)
print(tab1, row.names = FALSE)
utils::write.csv(tab1, file.path(OUT, "Table1_cohorts.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 补充表汇总工作簿 Tables_all.xlsx
# -----------------------------------------------------------------------------
R1 <- function(...) file.path(DIR_RES, ...)
SUPP <- list(
  "S1a MTRG catalog"           = R1("01_geneset", "Table_S1_MTRG_catalog.csv"),
  "S1b Gene sets"              = file.path(DIR_GS, "MTRG_genesets_table.csv"),
  "S1c Gene set coverage"      = R1("04_bulk_score", "geneset_coverage_by_cohort.csv"),
  "S2a scRNA samples"          = R1("03_scRNA", "scRNA_samples.csv"),
  "S2b Cluster annotation"     = R1("03_scRNA", "cluster_annotation_auto.csv"),
  "S2c copyKAT reference"      = R1("03_scRNA", "copykat_reference_cells.csv"),
  "S2d copyKAT by sample"      = R1("03_scRNA", "copykat_by_sample.csv"),
  "S2e copyKAT normal control" = R1("03_scRNA", "copykat_normal_control.csv"),
  "S2f Score per sample"       = R1("03_scRNA", "MT_score_per_sample_class.csv"),
  "S2g Pseudobulk DE"          = R1("03_scRNA", "pseudobulk_DE_High_vs_Low.csv"),
  "S2h Hallmark paired"        = R1("03_scRNA", "Hallmark_High_vs_Low_paired.csv"),
  "S2i Hallmark G1 only"       = R1("03_scRNA", "Hallmark_all_vs_G1only.csv"),
  "S2j Cell cycle by sample"   = R1("03_scRNA", "cellcycle_by_sample.csv"),
  "S2k CellChat"               = R1("03_scRNA", "CellChat_count_weight_from_malignant.csv"),
  "S3a Spearman TCGA"          = R1("04_bulk_score", "signature_spearman_TCGA.csv"),
  "S3b Prolif by cohort"       = R1("04_bulk_score", "proliferation_correlation_by_cohort.csv"),
  "S3c Prolif overlap"         = R1("04_bulk_score", "proliferation_overlap_summary.csv"),
  "S3d Per-gene rho Prolif"    = R1("04_bulk_score", "proliferation_overlap_per_gene_rho.csv"),
  "S3e Mutation frequency"     = R1("04_bulk_score", "TCGA_mutation_frequency.csv"),
  "S3f Score by mutation"      = R1("04_bulk_score", "MT_score_by_mutation.csv"),
  "S3g TCGA Cox"               = R1("04_bulk_score", "TCGA_Cox_MT_scores.csv"),
  "S3h External Cox MT_core"   = R1("04_bulk_score", "External_Cox_MT_core.csv"),
  "S4 Candidate genes"         = R1("05_candidates", "Table_S3_candidates_unicox.csv"),
  "S5a Model combinations"     = R1("06_ml_signature", "Table_S4_Cindex_101_models.csv"),
  "S5b MTRS coefficients"      = R1("09_hub_gene_spatial", "RLRS_coefficients.csv"),
  "S5c MTRS by cohort"         = R1("06_ml_signature", "RLRS_Cox_by_cohort.csv"),
  "S5d MTRS timeROC"           = R1("06_ml_signature", "RLRS_timeROC_AUC.csv"),
  "S5e MTRS genes prolif"      = R1("09_hub_gene_spatial", "RLRS_genes_proliferation.csv"),
  "S5f MTRS gene coverage"     = R1("09_hub_gene_spatial", "RLRS_gene_coverage.csv"),
  "S6a Random gene sets"       = R1("06b_specificity", "random_signature_benchmark.csv"),
  "S6b Random by cohort"       = R1("06b_specificity", "random_signature_benchmark_by_cohort.csv"),
  "S6c Proliferation-adjusted" = R1("06b_specificity", "Table_S6_proliferation_adjusted.csv"),
  "S6d Published signatures"   = R1("06_ml_signature", "Table_S5_published_models_Cindex.csv"),
  "S6e Signature sources"      = file.path(DIR_GS, "published_KIRC_signatures_sources.csv"),
  "S6f MTRS rank"              = R1("06_ml_signature", "RLRS_rank_among_published.csv"),
  "S7a TCGA Cox"               = R1("07_clinical_nomogram", "TCGA_uni_multi_Cox.csv"),
  "S7b External multivariable" = R1("07_clinical_nomogram", "External_multivariate_Cox.csv"),
  "S7c C-index nomogram"       = R1("07_clinical_nomogram", "Cindex_nomogram_RLRS_stage.csv"),
  "S8a Immune infiltration"    = R1("08_tme_therapy", "TME_High_vs_Low_all_methods.csv"),
  "S8b Checkpoints"            = R1("08_tme_therapy", "Checkpoints_High_vs_Low.csv"),
  "S8c TIDE"                   = R1("08_tme_therapy", "TIDE_High_vs_Low.csv"),
  "S8d Immunotherapy"          = R1("08_tme_therapy", "ICI_response_summary.csv"),
  "S8e Drug sensitivity"       = R1("08_tme_therapy", "oncoPredict_High_vs_Low.csv"),
  "S9a SHAP ranking"           = R1("09_hub_gene_spatial", "SHAP_gene_ranking.csv"),
  "S9b Model genes scRNA"      = R1("09_hub_gene_spatial", "hub_scRNA_by_cell_class.csv"),
  "S9c Spatial programmes"     = R1("09_hub_gene_spatial", "Spatial_RLcore_vs_programs_by_slice.csv"),
  "S9d Spatial model genes"    = R1("09_hub_gene_spatial", "Spatial_RLcore_vs_hub_spearman.csv"),
  "S9e Spatial depth check"    = R1("09_hub_gene_spatial", "Spatial_depth_check_by_slice.csv"),
  "S9f Spatial genes adjusted" = R1("09_hub_gene_spatial", "Spatial_hub_partial_by_slice.csv")
)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Table 1"); openxlsx::writeData(wb, "Table 1", tab1)
clean_df_labels <- function(df) {
  fix_str <- function(x) {
    x <- gsub("RLRS", "MTRS", x, fixed = TRUE)
    x <- gsub("RL_core", "MT_core", x, fixed = TRUE)
    x <- gsub("RLcore", "MT_core", x, fixed = TRUE)
    x <- gsub("RL_all", "MT_all", x, fixed = TRUE)
    x <- gsub("RHI", "MTI", x, fixed = TRUE)
    x
  }
  names(df) <- fix_str(names(df))
  for (j in seq_along(df)) {
    if (is.character(df[[j]])) df[[j]] <- fix_str(df[[j]])
    else if (is.factor(df[[j]])) df[[j]] <- fix_str(as.character(df[[j]]))
  }
  df
}
for (nm in names(SUPP)) {
  if (file.exists(SUPP[[nm]])) {
    openxlsx::addWorksheet(wb, nm)
    df_in <- clean_df_labels(utils::read.csv(SUPP[[nm]], check.names = FALSE, stringsAsFactors = FALSE))
    openxlsx::writeData(wb, nm, df_in)
  }
}

pk <- c("Seurat", "harmony", "copykat", "UCell", "CellChat", "edgeR", "DESeq2", "WGCNA", "GSVA", "msigdbr",
        "Mime1", "glmnet", "randomForestSRC", "survival", "timeROC", "rms", "dcurves", "IOBR", "oncoPredict",
        "shapviz", "TCGAbiolinks", "GEOquery")
S10 <- data.frame(
  package = pk,
  version = vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), "")
)
S10 <- rbind(data.frame(package = "R", version = paste(R.version$major, R.version$minor, sep = ".")), S10)
utils::write.csv(S10, file.path(OUT, "Table_S10_package_versions.csv"), row.names = FALSE)
openxlsx::addWorksheet(wb, "S10 Package versions"); openxlsx::writeData(wb, "S10 Package versions", S10)
openxlsx::saveWorkbook(wb, file.path(OUT, "Tables_all.xlsx"), overwrite = TRUE)
save_session(STEP)
log_msg("10_tables 完成")
