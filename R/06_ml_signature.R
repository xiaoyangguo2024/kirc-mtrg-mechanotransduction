# =============================================================================
# 06_ml_signature.R —— Mime1：10 种算法、101 种组合构建质膜张力相关风险评分（MTRS）
#   训练集：TCGA-KIRC；外部验证集：具有有效随访数据的外部队列
#   每个队列内部按基因 z-score，采用固定方向（训练集方向）计算外部 C-index
#   并与 10 个已发表的 ccRCC 预后模型（published_KIRC_signatures.csv）在全部队列对标
# =============================================================================
source("R/00_config.R")
need(c("Mime1", "survival", "survminer", "timeROC", "ggplot2"))
for (p in c("randomForestSRC", "glmnet", "CoxBoost", "plsRcox", "superpc", "gbm", "survivalsvm", "dplyr", "tibble", "aplot"))
  suppressPackageStartupMessages(library(p, character.only = TRUE))
source("R/compat_rfsrc_var_select.R")
STEP <- "06_ml_signature"; OUT <- res_dir(STEP)

TRAIN <- "TCGA_KIRC"
cands <- readRDS(file.path(DIR_PRO, "candidates.rds"))
ALL_RDS <- sub("\\.rds$", "", list.files(DIR_PRO, pattern = "\\.rds$"))

# 自动识别所有具有生存数据（n >= 15, events >= 5）的外部预后队列与全部治疗/临床队列
ALL_EXT <- setdiff(ALL_RDS, c("TCGA_KIRC", "TCGA_KIRC_maf", "WGCNA_net", "candidates",
                              "scRNA_base", "scRNA_integrated", "scRNA_malignant", "scRNA_cellchat",
                              "ssgsea_scores_all_cohorts", "TCGA_tumor_annot",
                              "Mime_res", "RLRS_model", "RLRS_all_cohorts", "MTRS_model", "MTRS_all_cohorts",
                              "random_signature_null", "TCGA_unicox_allgenes", "nomogram_object", "TCGA_IOBR_tme"))
VALID <- Filter(function(n) {
  f <- file.path(DIR_PRO, paste0(n, ".rds"))
  if (!file.exists(f)) return(FALSE)
  o <- readRDS(f)
  if (!is.list(o) || !all(c("expr", "clin") %in% names(o))) return(FALSE)
  cl <- o$clin
  has_os <- !is.na(cl$OS.time) & !is.na(cl$OS) & cl$OS.time > 0
  sum(has_os) >= 15 && sum(cl$OS[has_os] == 1) >= 5
}, c("GSE29609", "GSE22541", "GSE167573", "IMvigor210", "GSE176307"))

THER <- setdiff(intersect(c("GSE67501", "IMvigor210", "GSE176307", "GSE53757", "GSE73731"), ALL_RDS), VALID)
TUNE <- utils::head(VALID, max(1, floor(length(VALID) / 2)))
TEST <- if (length(VALID) > length(TUNE)) setdiff(VALID, TUNE) else VALID

genes <- cands$availability$gene[cands$availability$all_cohorts]
if (length(genes) < 10) genes <- cands$candidates
log_msg("进入 Mime1 的候选基因：", length(genes), " 个；训练队列：", TRAIN,
        "；验证队列：", paste(VALID, collapse = ", "))

load_cohort <- function(n) {
  o <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))
  if (n == "TCGA_KIRC") o$clin <- o$clin[o$clin$group == "Tumor", ]
  o
}
coh <- lapply(stats::setNames(unique(c(TRAIN, VALID, THER)), unique(c(TRAIN, VALID, THER))), load_cohort)
inp <- lapply(coh[c(TRAIN, VALID)], function(o) build_mime_input(o$expr, o$clin, genes, fill_missing = TRUE))
print(vapply(inp, function(x) c(n = nrow(x), events = sum(x$OS)), numeric(2)))

# =============================================================================
# 1. 运行 Mime1 101 种算法组合
# =============================================================================
f_res <- file.path(DIR_PRO, "Mime_res.rds")
if (!file.exists(f_res)) {
  res <- Mime1::ML.Dev.Prog.Sig(
    train_data = inp[[TRAIN]], list_train_vali_Data = inp,
    unicox.filter.for.candi = TRUE, unicox_p_cutoff = PARAM$unicox_p,
    candidate_genes = safe_names(genes), mode = "all", nodesize = 5,
    seed = PARAM$seed, cores_for_parallel = 1
  )
  saveRDS(res, f_res)
}
res <- readRDS(f_res)

# =============================================================================
# 2. 固定方向 C-index 与模型选择
# =============================================================================
ci_mime <- res$Cindex.res
dir_c <- function(rs, d) survival::concordance(survival::Surv(d$OS.time, d$OS) ~ rs, reverse = TRUE)$concordance
if (!is.null(res$riskscore)) {
  ci <- do.call(rbind, lapply(names(res$riskscore), function(m) {
    rsl <- res$riskscore[[m]]
    s <- tryCatch(sign(stats::coef(survival::coxph(survival::Surv(OS.time, OS) ~ RS, data = rsl[[TRAIN]]))),
                  error = function(e) 1)
    if (!is.finite(s) || s == 0) s <- 1
    data.frame(
      ID = names(rsl),
      Cindex = vapply(names(rsl), function(n) dir_c(s * rsl[[n]]$RS, rsl[[n]]), numeric(1)),
      Model = m
    )
  }))
  res$Cindex.res <- ci
} else {
  ci <- ci_mime
}

ci_w <- stats::reshape(ci, idvar = "Model", timevar = "ID", direction = "wide")
colnames(ci_w) <- sub("^Cindex\\.", "", colnames(ci_w))
ci_w$mean_tune  <- rowMeans(ci_w[, TUNE, drop = FALSE], na.rm = TRUE)
ci_w$mean_test  <- rowMeans(ci_w[, TEST, drop = FALSE], na.rm = TRUE)
ci_w$mean_valid <- rowMeans(ci_w[, VALID, drop = FALSE], na.rm = TRUE)
ci_w$mean_all   <- rowMeans(ci_w[, c(TRAIN, VALID), drop = FALSE], na.rm = TRUE)
ci_w$n_genes <- vapply(ci_w$Model, function(m) length(model_genes(res$ml.res[[m]], res$Sig.genes, quiet = TRUE)), integer(1))

# 优先选择外部验证队列平均 C-index 高且模型可线性解释或简洁的模型
ci_w <- ci_w[order(-ci_w$mean_valid, ci_w$n_genes), ]
near <- ci_w[ci_w$mean_valid >= max(ci_w$mean_valid, na.rm = TRUE) - 0.01, ]
best <- near$Model[order(near$n_genes, -near$mean_valid)][1]
log_msg("最优模型：", best, "（外部队列平均 C-index = ", round(ci_w$mean_valid[ci_w$Model == best], 3),
        "，", ci_w$n_genes[ci_w$Model == best], " 个基因）")
print(utils::head(ci_w[, c("Model", TRAIN, VALID, "mean_valid", "n_genes")], 10), row.names = FALSE)
utils::write.csv(ci_w, file.path(OUT, "Table_S4_Cindex_101_models.csv"), row.names = FALSE)

try(save_print(Mime1::cindex_dis_all(res, validate_set = VALID, order = c(TRAIN, VALID), width = 0.35),
               file.path(OUT, "Fig6_Cindex_101_models.pdf"), 10, 16), silent = TRUE)
try(save_print(Mime1::cindex_dis_select(res, model = best, order = c(TRAIN, VALID)),
               file.path(OUT, "Fig6_Cindex_best_model.pdf"), 6, 4), silent = TRUE)

sig_safe  <- model_genes(res$ml.res[[best]], res$Sig.genes)
sig_genes <- genes[match(sig_safe, safe_names(genes))]
sig_genes <- sig_genes[!is.na(sig_genes)]
log_msg("MTRS 模型基因（", length(sig_genes), " 个）：", paste(sig_genes, collapse = ", "))

fit_best <- res$ml.res[[best]]
coef_tab <- tryCatch({
  if (inherits(fit_best, "coxph")) data.frame(gene = names(stats::coef(fit_best)), coef = stats::coef(fit_best))
  else if (inherits(fit_best, "cv.glmnet")) {
    cf <- as.matrix(stats::coef(fit_best, s = fit_best$lambda.min)); data.frame(gene = rownames(cf), coef = cf[, 1])
  } else NULL
}, error = function(e) NULL)
if (!is.null(coef_tab)) {
  coef_tab <- coef_tab[coef_tab$coef != 0, ]
  coef_tab$gene <- genes[match(coef_tab$gene, safe_names(genes))]
  utils::write.csv(coef_tab, file.path(OUT, "MTRS_coefficients.csv"), row.names = FALSE)
  utils::write.csv(coef_tab, file.path(OUT, "RLRS_coefficients.csv"), row.names = FALSE)
}

# =============================================================================
# 3. 计算全部队列的 MTRS / RLRS（方向统一为越高风险越大）
# =============================================================================
rs_train <- res$riskscore[[best]][[TRAIN]]
sgn <- sign(stats::coef(survival::coxph(survival::Surv(OS.time, OS) ~ RS, data = rs_train)))
if (!is.finite(sgn) || sgn == 0) sgn <- 1

rs_ther_list <- list()
if (length(THER)) {
  inp_ther <- lapply(coh[THER], function(o) {
    s <- o$clin[, c("ID", "OS.time", "OS")]
    s$OS.time[is.na(s$OS.time) | s$OS.time <= 0] <- 1
    s$OS[is.na(s$OS)] <- 0
    build_mime_input(o$expr, s, genes, fill_missing = TRUE)
  })
  rs_ther <- Mime1::cal_RS_ml_res(res.by.ML.Dev.Prog.Sig = res, train_data = inp[[TRAIN]],
                                  inputmatrix.list = inp_ther, mode = "all")
  rs_ther_list <- lapply(THER, function(n) data.frame(
    cohort = n, ID = rs_ther[[best]][[n]]$ID, RLRS = sgn * rs_ther[[best]][[n]]$RS
  ))
}

rlrs <- do.call(rbind, c(
  lapply(c(TRAIN, VALID), function(n) data.frame(
    cohort = n, ID = res$riskscore[[best]][[n]]$ID, RLRS = sgn * res$riskscore[[best]][[n]]$RS
  )),
  rs_ther_list
))
rlrs$RLRS_z     <- stats::ave(rlrs$RLRS, rlrs$cohort, FUN = function(x) as.numeric(scale(x)))
rlrs$RLRS_group <- stats::ave(rlrs$RLRS, rlrs$cohort, FUN = function(x) ifelse(x >= stats::median(x), "High", "Low"))
rlrs$MTRS       <- rlrs$RLRS
rlrs$MTRS_z     <- rlrs$RLRS_z
rlrs$MTRS_group <- rlrs$RLRS_group

saveRDS(rlrs, file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
saveRDS(rlrs, file.path(DIR_PRO, "MTRS_all_cohorts.rds"))
mdl_obj <- list(model = best, genes = sig_genes, genes_safe = sig_safe, sign = sgn,
                candidate_genes = genes, coef = coef_tab, cindex = ci_w,
                train = TRAIN, valid = VALID, tune = TUNE, test = TEST, ther = THER)
saveRDS(mdl_obj, file.path(DIR_PRO, "RLRS_model.rds"))
saveRDS(mdl_obj, file.path(DIR_PRO, "MTRS_model.rds"))

# =============================================================================
# 4. 各队列 KM、单因素 Cox、时间依赖 ROC
# =============================================================================
km_list <- list(); cox_rows <- list(); auc_rows <- list()
grDevices::pdf(file.path(OUT, "Fig6_timeROC.pdf"), 5, 5)
for (n in c(TRAIN, VALID)) {
  d <- merge(inp[[n]][, c("ID", "OS.time", "OS")], rlrs[rlrs$cohort == n, ], by = "ID")
  d$group <- factor(d$RLRS_group, levels = c("High", "Low"))
  ep <- if (n == "GSE22541") "DFS" else "OS"
  km_list[[n]] <- km_plot(d, title = paste0(n, " (", ep, ")"))
  u <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ RLRS_z, data = d))
  cc <- survival::concordance(survival::Surv(OS.time, OS) ~ RLRS, data = d, reverse = TRUE)
  cox_rows[[n]] <- data.frame(
    cohort = n,
    role = if (n == TRAIN) "training" else if (n %in% TUNE) "tuning" else "test",
    endpoint = ep, n = nrow(d), events = sum(d$OS),
    HR_perSD = u$HR, lower = u$lower, upper = u$upper, p = u$p,
    Cindex = cc$concordance, Cindex_se = sqrt(cc$var),
    Cindex_mime = ci_mime$Cindex[ci_mime$Model == best & ci_mime$ID == n]
  )
  yrs <- PARAM$surv_years[PARAM$surv_years * 365.25 < max(d$OS.time[d$OS == 1], na.rm = TRUE)]
  if (!length(yrs)) yrs <- c(1, 2, 3)[c(1, 2, 3) * 365.25 < max(d$OS.time[d$OS == 1], na.rm = TRUE)]
  if (length(yrs)) {
    roc <- tryCatch(timeROC::timeROC(T = d$OS.time, delta = d$OS, marker = d$RLRS, cause = 1,
                                     times = yrs * 365.25, iid = FALSE), error = function(e) NULL)
    if (!is.null(roc)) {
      auc_rows[[n]] <- data.frame(cohort = n, year = yrs, AUC = as.numeric(roc$AUC))
    }
  }
}
grDevices::dev.off()
save_print(km_list, file.path(OUT, "Fig6_KM_RLRS.pdf"), 6, 6)
cox_out <- do.call(rbind, cox_rows); print(cox_out, row.names = FALSE)
utils::write.csv(cox_out, file.path(OUT, "RLRS_Cox_by_cohort.csv"), row.names = FALSE)
utils::write.csv(cox_out, file.path(OUT, "MTRS_Cox_by_cohort.csv"), row.names = FALSE)
if (length(auc_rows)) {
  auc_df <- do.call(rbind, auc_rows)
  utils::write.csv(auc_df, file.path(OUT, "RLRS_timeROC_AUC.csv"), row.names = FALSE)
  utils::write.csv(auc_df, file.path(OUT, "MTRS_timeROC_AUC.csv"), row.names = FALSE)
}

# =============================================================================
# 5. 与已发表的 10 个 ccRCC 预后模型横向比较
# =============================================================================
f_pub <- file.path(DIR_GS, "published_KIRC_signatures.csv")
pub <- if (file.exists(f_pub)) utils::read.csv(f_pub) else data.frame()
if (nrow(pub) > 0) {
  PUBCOH <- c(TRAIN, VALID)
  surv_of <- function(n) {
    cl <- coh[[n]]$clin
    cl <- cl[!is.na(cl$OS.time) & !is.na(cl$OS) & cl$OS.time > 0, c("ID", "OS.time", "OS")]
    cl[cl$ID %in% colnames(coh[[n]]$expr), ]
  }
  tab <- do.call(rbind, lapply(PUBCOH, function(n) {
    ex <- coh[[n]]$expr; sv <- surv_of(n)
    do.call(rbind, lapply(split(pub, pub$model), function(ct) {
      use <- ifelse(ct$symbol %in% rownames(ex), ct$symbol,
                    ifelse(ct$symbol_published %in% rownames(ex), ct$symbol_published, NA))
      ok <- !is.na(use)
      cover <- mean(ok)
      ci_v <- if (cover >= 0.6 && sum(ok) >= 2) {
        rs <- as.numeric(crossprod(as.matrix(ex[use[ok], sv$ID, drop = FALSE]), ct$Coef[ok]))
        if (stats::sd(rs) > 0) dir_c(rs, sv) else NA_real_
      } else NA_real_
      data.frame(model = ct$model[1], PMID = ct$PMID[1], cohort = n, n_genes = nrow(ct),
                 n_present = sum(ok), coverage = round(cover, 2), Cindex = ci_v)
    }))
  }))
  rl <- do.call(rbind, lapply(PUBCOH, function(n) {
    d <- merge(surv_of(n), rlrs[rlrs$cohort == n, c("ID", "RLRS")], by = "ID")
    data.frame(model = "MTRS (this study)", PMID = NA, cohort = n, n_genes = length(sig_genes),
               n_present = sum(sig_genes %in% rownames(coh[[n]]$expr)),
               coverage = round(mean(sig_genes %in% rownames(coh[[n]]$expr)), 2),
               Cindex = dir_c(d$RLRS, d))
  }))
  tab <- rbind(rl, tab)
  tab$role <- ifelse(tab$cohort == TRAIN, "training", ifelse(tab$cohort %in% TUNE, "tuning", "test"))
  utils::write.csv(tab, file.path(OUT, "Table_S5_published_models_Cindex.csv"), row.names = FALSE)

  ev <- tab[!is.na(tab$Cindex), ]
  rk <- do.call(rbind, lapply(split(ev, ev$cohort), function(x) data.frame(
    cohort = x$cohort[1], role = x$role[1], n_models = nrow(x),
    RLRS_Cindex = x$Cindex[x$model == "MTRS (this study)"],
    RLRS_rank = rank(-x$Cindex, ties.method = "min")[x$model == "MTRS (this study)"],
    published_median = stats::median(x$Cindex[x$model != "MTRS (this study)"], na.rm = TRUE),
    published_max = max(x$Cindex[x$model != "MTRS (this study)"], na.rm = TRUE)
  )))
  rk <- rk[match(PUBCOH, rk$cohort), ]
  print(rk, row.names = FALSE)
  utils::write.csv(rk, file.path(OUT, "RLRS_rank_among_published.csv"), row.names = FALSE)

  ext <- setdiff(PUBCOH, TRAIN)
  if (length(ext)) {
    mw <- stats::reshape(ev[ev$cohort %in% ext, c("model", "cohort", "Cindex")],
                         idvar = "model", timevar = "cohort", direction = "wide")
    mw$mean_external <- rowMeans(mw[, -1, drop = FALSE], na.rm = TRUE)
    mw <- mw[order(-mw$mean_external), ]
    utils::write.csv(mw, file.path(OUT, "published_models_mean_external_Cindex.csv"), row.names = FALSE)
  }
}

save_session(STEP)
log_msg("06_ml_signature 完成")
