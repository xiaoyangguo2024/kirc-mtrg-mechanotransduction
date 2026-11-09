# =============================================================================
# 08_tme_therapy.R —— 肿瘤免疫微环境（TME）与抗 PD-1/PD-L1 免疫治疗及肾癌靶向药敏感性
#   1. TCGA-KIRC：IOBR 7 种免疫反卷积算法（CIBERSORT / EPIC / MCPcounter / quanTIseq / TIMER / xCell / ESTIMATE）
#   2. 免疫检查点与血管生成/髓系炎症特征（IMmotion150 Signatures）
#   3. 本地 TIDE（复用 BLCA_DIR/conda_py/bin/python 运行 tidepy）
#   4. 真实免疫治疗队列应答验证：GSE67501（RCC Nivolumab anti-PD-1）、IMvigor210（anti-PD-L1）、GSE176307（anti-PD-1/PD-L1）
#   5. GDSC2 靶向药敏预测（oncoPredict，复用 BLCA_DIR/data/raw/oncoPredict，聚焦 ccRCC 一线 TKI 与 mTOR 抑制剂）
# =============================================================================
source("R/00_config.R")
need(c("survival", "survminer", "pROC", "ggplot2", "ggpubr"))
STEP <- "08_tme_therapy"; OUT <- res_dir(STEP)

tcga <- readRDS(file.path(DIR_PRO, "TCGA_KIRC.rds"))
tum  <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
rlrs <- readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
d <- merge(tum, rlrs[rlrs$cohort == "TCGA_KIRC", c("ID", "RLRS", "RLRS_z", "RLRS_group")], by = "ID")
d$group <- factor(d$RLRS_group, levels = c("High", "Low"))
PAL <- c(High = "#D7301F", Low = "#2C7FB8")

wilcox_table <- function(df, vars, group = "group") {
  out <- do.call(rbind, lapply(vars, function(v) {
    x <- df[[v]]; g <- df[[group]]
    data.frame(
      variable = v,
      median_High = stats::median(x[g == "High"], na.rm = TRUE),
      median_Low  = stats::median(x[g == "Low"], na.rm = TRUE),
      rho_RLRS    = suppressWarnings(stats::cor(x, df$RLRS, method = "spearman", use = "complete.obs")),
      p           = tryCatch(stats::wilcox.test(x ~ g)$p.value, error = function(e) NA_real_)
    )
  }))
  out$FDR <- stats::p.adjust(out$p, "BH")
  out
}

# =============================================================================
# 1. IOBR 7 种免疫浸润反卷积（输入非 log 的 TPM，TIMER group_list = "kirc"）
# =============================================================================
f_tme <- file.path(DIR_PRO, "TCGA_IOBR_tme.rds")
if (!file.exists(f_tme) && requireNamespace("IOBR", quietly = TRUE)) {
  tpm <- tcga$tpm[, d$ID]
  meths <- c("cibersort", "mcpcounter", "xcell", "epic", "timer", "quantiseq", "estimate")
  tme <- lapply(meths, function(m) tryCatch({
    args <- list(eset = tpm, method = m, arrays = FALSE)
    if (m == "cibersort") args$perm <- 200
    if (m == "timer") args$group_list <- rep("kirc", ncol(tpm))
    if (m == "quantiseq") { args$tumor <- TRUE; args$scale_mrna <- TRUE }
    log_msg("IOBR 反卷积：", m)
    do.call(IOBR::deconvo_tme, args)
  }, error = function(e) { log_msg("IOBR ", m, " 失败：", conditionMessage(e)); NULL }))
  tme <- Reduce(function(a, b) merge(a, b, by = "ID"), Filter(Negate(is.null), tme))
  saveRDS(tme, f_tme)
}
if (file.exists(f_tme)) {
  tme <- readRDS(f_tme)
  dt <- merge(d[, c("ID", "RLRS", "group")], tme, by = "ID")
  cells <- setdiff(colnames(tme), "ID")
  cells <- cells[vapply(dt[, cells], function(x) is.numeric(x) && stats::sd(x, na.rm = TRUE) > 0, logical(1))]
  wt <- wilcox_table(dt, cells)
  utils::write.csv(wt[order(wt$p), ], file.path(OUT, "TME_High_vs_Low_all_methods.csv"), row.names = FALSE)
}

# =============================================================================
# 2. 免疫检查点
# =============================================================================
cp <- intersect(CHECKPOINTS, rownames(tcga$expr))
dc <- data.frame(d[, c("ID", "RLRS", "group")], t(tcga$expr[cp, d$ID]), check.names = FALSE)
utils::write.csv(wilcox_table(dc, cp), file.path(OUT, "Checkpoints_High_vs_Low.csv"), row.names = FALSE)

# =============================================================================
# 3. 本地 TIDE（调用 BLCA_DIR/conda_py/bin/python）
# =============================================================================
tide_in <- tcga$expr[, d$ID]
tide_in <- tide_in - rowMeans(tide_in)
utils::write.table(round(tide_in, 4), file.path(OUT, "TIDE_input.txt"), sep = "\t", quote = FALSE, col.names = NA)
f_tide <- file.path(OUT, "TIDE_output.csv")
py <- file.path(BLCA_DIR, "conda_py", "bin", "python")
if (!file.exists(py)) py <- file.path(PROJ, "conda_py", "bin", "python")
if (!file.exists(f_tide) && file.exists(py)) {
  log_msg("用本地 tidepy 计算 TCGA-KIRC TIDE 评分...")
  system2(py, c("-I", shQuote(file.path(PROJ, "R", "tide_run.py")),
                shQuote(file.path(OUT, "TIDE_input.txt")), shQuote(f_tide)))
}
if (file.exists(f_tide)) {
  ti <- utils::read.csv(f_tide, check.names = FALSE)
  dt <- merge(d[, c("ID", "RLRS", "group")], ti, by.x = "ID", by.y = "Patient")
  vars <- intersect(c("TIDE", "Dysfunction", "Exclusion", "MDSC", "CAF", "TAM M2", "IFNG", "CD8"), colnames(dt))
  utils::write.csv(wilcox_table(dt, vars), file.path(OUT, "TIDE_High_vs_Low.csv"), row.names = FALSE)
}

# =============================================================================
# 4. 免疫治疗真实队列验证：GSE67501（RCC anti-PD-1 Nivolumab）、IMvigor210、GSE176307
# =============================================================================
ici_rows <- list()
for (n in c("GSE67501", "IMvigor210", "GSE176307")) {
  f_rds <- file.path(DIR_PRO, paste0(n, ".rds"))
  if (!file.exists(f_rds) || !n %in% unique(rlrs$cohort)) next
  cl <- readRDS(f_rds)$clin
  e <- merge(cl, rlrs[rlrs$cohort == n, c("ID", "RLRS", "RLRS_z", "RLRS_group")], by = "ID")
  e$group <- factor(e$RLRS_group, levels = c("High", "Low"))
  r <- e[e$response %in% c("CR/PR", "SD/PD"), ]
  if (nrow(r) < 5 || length(unique(r$response)) < 2) next
  r$response <- factor(r$response, levels = c("SD/PD", "CR/PR"))
  roc <- pROC::roc(r$response, r$RLRS, levels = c("SD/PD", "CR/PR"), direction = "<", quiet = TRUE)
  ci  <- tryCatch(as.numeric(pROC::ci.auc(roc)), error = function(err) c(NA_real_, as.numeric(roc$auc), NA_real_))
  tb  <- table(r$group, r$response)
  ici_rows[[n]] <- data.frame(
    cohort = n, n_evaluable = nrow(r), n_responders = sum(r$response == "CR/PR"),
    AUC = ci[2], AUC_lower = ci[1], AUC_upper = ci[3],
    p_wilcox = tryCatch(stats::wilcox.test(RLRS ~ response, r)$p.value, error = function(err) NA_real_),
    ORR_High = mean(r$response[r$group == "High"] == "CR/PR"),
    ORR_Low  = mean(r$response[r$group == "Low"] == "CR/PR"),
    p_fisher = tryCatch(stats::fisher.test(tb)$p.value, error = function(err) NA_real_)
  )
}
if (length(ici_rows)) {
  ici_tab <- do.call(rbind, ici_rows); print(ici_tab, row.names = FALSE)
  utils::write.csv(ici_tab, file.path(OUT, "ICI_response_summary.csv"), row.names = FALSE)
}

# =============================================================================
# 5. GDSC2 药物敏感性预测（oncoPredict，直接复用 BLCA_DIR/data/raw/oncoPredict）
# =============================================================================
op_dir <- file.path(BLCA_DIR, "data", "raw", "oncoPredict")
if (!dir.exists(op_dir)) op_dir <- file.path(DIR_RAW, "oncoPredict")
f_gx <- file.path(op_dir, "GDSC2_Expr (RMA Normalized and Log Transformed).rds")
f_gr <- file.path(op_dir, "GDSC2_Res.rds")
f_pred <- file.path(OUT, "oncoPredict", "calcPhenotype_Output", "DrugPredictions.csv")
if (!file.exists(f_pred) && file.exists(f_gx) && file.exists(f_gr) && requireNamespace("oncoPredict", quietly = TRUE)) {
  dir.create(file.path(OUT, "oncoPredict"), showWarnings = FALSE, recursive = TRUE)
  old <- setwd(file.path(OUT, "oncoPredict"))
  oncoPredict::calcPhenotype(
    trainingExprData = readRDS(f_gx), trainingPtype = exp(readRDS(f_gr)),
    testExprData = as.matrix(tcga$expr[, d$ID]), batchCorrect = "eb",
    powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.2, minNumSamples = 10,
    printOutput = TRUE, removeLowVaringGenesFrom = "rawData", folder = TRUE
  )
  setwd(old)
}
if (file.exists(f_pred)) {
  pred <- utils::read.csv(f_pred, row.names = 1, check.names = FALSE)
  dp <- data.frame(d[match(rownames(pred), d$ID), c("ID", "RLRS", "group")], pred, check.names = FALSE)
  drugs <- colnames(pred)
  wt <- wilcox_table(dp, drugs)
  utils::write.csv(wt[order(wt$p), ], file.path(OUT, "oncoPredict_High_vs_Low.csv"), row.names = FALSE)
}

save_session(STEP)
log_msg("08_tme_therapy 完成")
