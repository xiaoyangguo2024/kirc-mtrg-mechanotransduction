# =============================================================================
# 06b_specificity.R —— MTRS 特异性检验：
#   1. 1,000 次随机基因集置换基准（Null A：全转录组随机表达基因；Null B：TCGA 单因素预后显著基因）
#   2. 增殖评分（Prolif = E2F + G2M）与分期（AJCC Stage）多变量校正及似然比检验（LRT）
# =============================================================================
source("R/00_config.R")
need(c("survival", "glmnet", "ggplot2"))
STEP <- "06b_specificity"; OUT <- res_dir(STEP)

mdl    <- readRDS(file.path(DIR_PRO, "RLRS_model.rds"))
rlrs   <- readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
scores <- readRDS(file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds"))
TRAIN  <- mdl$train
VALID  <- mdl$valid
OSCOH  <- c(TRAIN, VALID)
N_RAND <- 1000

load_cohort <- function(n) {
  o <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))
  if (n == TRAIN) o$clin <- o$clin[o$clin$group == "Tumor", ]
  o
}
coh <- lapply(stats::setNames(OSCOH, OSCOH), load_cohort)
sig <- mdl$genes; k <- length(sig)
log_msg("MTRS 基因数 k = ", k, "：", paste(sig, collapse = ", "))

# =============================================================================
# 1. 1,000 次随机基因集基准检验
# =============================================================================
cidx <- function(lp, d) survival::concordance(survival::Surv(OS.time, OS) ~ lp,
                                              data = data.frame(d, lp = as.numeric(lp)), reverse = TRUE)$concordance
ridge_cindex <- function(genes, seed) {
  inp <- lapply(coh[c(TRAIN, VALID)], function(o) build_mime_input(o$expr, o$clin, genes, fill_missing = TRUE))
  x <- as.matrix(inp[[TRAIN]][, -(1:3)]); y <- survival::Surv(inp[[TRAIN]]$OS.time, inp[[TRAIN]]$OS)
  set.seed(seed)
  cv <- glmnet::cv.glmnet(x, y, family = "cox", alpha = 0, nfolds = 10)
  b <- as.matrix(stats::coef(cv, s = "lambda.min"))[, 1]
  vapply(VALID, function(n) cidx(as.matrix(inp[[n]][, names(b), drop = FALSE]) %*% b, inp[[n]]), numeric(1))
}

univ <- rownames(coh[[TRAIN]]$expr)
univ <- univ[rowMeans(coh[[TRAIN]]$expr[univ, ]) >= 1]
# 要求至少在多数验证队列中存在
in_val <- rowSums(sapply(VALID, function(n) univ %in% rownames(coh[[n]]$expr)))
univ <- univ[in_val >= max(1, length(VALID) - 1)]
log_msg("随机抽样的基因全集：", length(univ), " 个")

f_uc <- file.path(DIR_PRO, "TCGA_unicox_allgenes.rds")
if (!file.exists(f_uc)) {
  tr <- coh[[TRAIN]]; ok <- !is.na(tr$clin$OS.time) & tr$clin$OS.time > 0 & !is.na(tr$clin$OS)
  ex <- zscore_rows(tr$expr[univ, tr$clin$ID[ok]])
  sv <- survival::Surv(tr$clin$OS.time[ok], tr$clin$OS[ok])
  p_uc <- vapply(univ, function(g) summary(survival::coxph(sv ~ ex[g, ]))$coefficients[1, 5], numeric(1))
  saveRDS(p_uc, f_uc)
}
p_uc <- readRDS(f_uc)
pool_B <- names(p_uc)[!is.na(p_uc) & p_uc < PARAM$unicox_p]
log_msg("零分布 B 基因池（TCGA 单因素 Cox P < ", PARAM$unicox_p, "）：", length(pool_B), " 个")

f_null <- file.path(DIR_PRO, "random_signature_null.rds")
if (!file.exists(f_null)) {
  set.seed(PARAM$seed)
  setsA <- replicate(N_RAND, sample(univ, k), simplify = FALSE)
  setsB <- replicate(N_RAND, sample(pool_B, min(k, length(pool_B))), simplify = FALSE)
  run <- function(sets, tag) {
    out <- matrix(NA_real_, length(sets), length(VALID), dimnames = list(NULL, VALID))
    for (i in seq_along(sets)) {
      out[i, ] <- tryCatch(ridge_cindex(sets[[i]], PARAM$seed + i), error = function(e) rep(NA_real_, length(VALID)))
      if (i %% 200 == 0) log_msg(tag, "：", i, "/", length(sets))
    }
    out
  }
  obs <- ridge_cindex(sig, PARAM$seed)
  nullA <- run(setsA, "零分布 A"); nullB <- run(setsB, "零分布 B")
  saveRDS(list(setsA = setsA, setsB = setsB, nullA = nullA, nullB = nullB, obs = obs), f_null)
}
nl <- readRDS(f_null)

ci_best <- unlist(mdl$cindex[mdl$cindex$Model == mdl$model, VALID])
summ <- function(null, label) {
  m <- rowMeans(null, na.rm = TRUE); m <- m[!is.na(m)]
  data.frame(
    null = label, n_valid_draws = length(m), null_mean = mean(m), null_sd = stats::sd(m),
    null_q95 = stats::quantile(m, 0.95),
    RLRS_ridge_mean = mean(nl$obs, na.rm = TRUE),
    p_RLRS_ridge = (1 + sum(m >= mean(nl$obs, na.rm = TRUE))) / (length(m) + 1),
    RLRS_best_mean = mean(ci_best, na.rm = TRUE),
    p_RLRS_best = (1 + sum(m >= mean(ci_best, na.rm = TRUE))) / (length(m) + 1),
    row.names = NULL
  )
}
bench <- rbind(summ(nl$nullA, "A: random expressed genes"), summ(nl$nullB, "B: random TCGA-prognostic genes"))
print(bench, row.names = FALSE)
utils::write.csv(bench, file.path(OUT, "random_signature_benchmark.csv"), row.names = FALSE)

per_coh <- do.call(rbind, lapply(VALID, function(n) data.frame(
  cohort = n, RLRS_ridge = nl$obs[n], RLRS_best = ci_best[n],
  pctl_A = mean(nl$nullA[, n] < nl$obs[n], na.rm = TRUE),
  pctl_B = mean(nl$nullB[, n] < nl$obs[n], na.rm = TRUE)
)))
utils::write.csv(per_coh, file.path(OUT, "random_signature_benchmark_by_cohort.csv"), row.names = FALSE)

# =============================================================================
# 2. 增殖评分（Prolif）与分期校正
# =============================================================================
adv_stage <- function(n, c0) {
  if ("stage" %in% colnames(c0)) {
    s <- toupper(as.character(c0$stage))
    return(ifelse(grepl("III|IV|3|4", s), 1L, ifelse(grepl("I|II|1|2", s), 0L, NA_integer_)))
  }
  rep(NA_integer_, nrow(c0))
}
rows <- list()
for (n in OSCOH) {
  c0 <- coh[[n]]$clin
  d <- merge(data.frame(ID = c0$ID, OS.time = c0$OS.time, OS = c0$OS, adv = adv_stage(n, c0)),
             rlrs[rlrs$cohort == n, c("ID", "RLRS")], by = "ID")
  d <- merge(d, scores[[n]][, c("ID", "Prolif", "MT_core", "RL_core")], by = "ID")
  d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
  d$RLRS_z   <- as.numeric(scale(d$RLRS))
  d$Prolif_z <- as.numeric(scale(d$Prolif))
  S <- survival::Surv(d$OS.time, d$OS)
  m0 <- survival::coxph(S ~ Prolif_z, data = d)
  m1 <- survival::coxph(S ~ Prolif_z + RLRS_z, data = d)
  mr <- survival::coxph(S ~ RLRS_z, data = d)
  r <- data.frame(
    cohort = n, n = nrow(d), events = sum(d$OS),
    rho_RLRS_Prolif = stats::cor(d$RLRS, d$Prolif, method = "spearman"),
    rho_RLRS_RLcore = stats::cor(d$RLRS, d$MT_core, method = "spearman"),
    C_Prolif = summary(m0)$concordance[1],
    C_RLRS = summary(mr)$concordance[1],
    C_Prolif_RLRS = summary(m1)$concordance[1]
  )
  a1 <- summary(m1)$coefficients
  r$HR_RLRS_adjProlif <- exp(a1["RLRS_z", 1]); r$p_RLRS_adjProlif <- a1["RLRS_z", 5]
  r$LRT_p_add_RLRS <- stats::anova(m0, m1)[2, "Pr(>|Chi|)"]
  if (sum(!is.na(d$adv)) > 20 && length(unique(stats::na.omit(d$adv))) == 2) {
    m2 <- survival::coxph(S ~ RLRS_z + Prolif_z + adv, data = d)
    a2 <- summary(m2)$coefficients
    r$HR_RLRS_adjProlifStage <- exp(a2["RLRS_z", 1]); r$p_RLRS_adjProlifStage <- a2["RLRS_z", 5]
  } else {
    r$HR_RLRS_adjProlifStage <- NA_real_; r$p_RLRS_adjProlifStage <- NA_real_
  }
  rows[[n]] <- r
}
adj <- do.call(rbind, rows); print(adj, row.names = FALSE, digits = 3)
utils::write.csv(adj, file.path(OUT, "Table_S6_proliferation_adjusted.csv"), row.names = FALSE)
save_session(STEP)
log_msg("06b_specificity 完成")
