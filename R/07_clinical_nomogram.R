# =============================================================================
# 07_clinical_nomogram.R —— MTRS 的独立预后价值与临床列线图（Nomogram）
#   单/多因素 Cox（TCGA-KIRC + 外部验证队列）、rms 列线图、1/3/5 年 1000 次 bootstrap 校准曲线、
#   决策曲线分析（DCA）、时间依赖 ROC 对比与临床亚组 KM
# =============================================================================
source("R/00_config.R")
need(c("survival", "survminer", "rms", "dcurves", "timeROC", "ggplot2"))
STEP <- "07_clinical_nomogram"; OUT <- res_dir(STEP)
NOMO_VARS <- c("RLRS", "age", "stage_bin", "grade_bin")

tum  <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
rlrs <- readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
mdl  <- readRDS(file.path(DIR_PRO, "RLRS_model.rds"))
d <- merge(tum, rlrs[rlrs$cohort == "TCGA_KIRC", c("ID", "RLRS", "RLRS_z", "RLRS_group")], by = "ID")
d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
d$time_y    <- d$OS.time / 365.25
d$age10     <- d$age / 10
d$age_group <- factor(ifelse(d$age > 60, ">60", "<=60"), levels = c("<=60", ">60"))
d$gender    <- factor(d$gender, levels = c("female", "male"))
d$stage_bin <- factor(ifelse(d$stage %in% c("III", "IV"), "III-IV",
                      ifelse(d$stage %in% c("I", "II"), "I-II", NA)), levels = c("I-II", "III-IV"))
d$grade_bin <- factor(ifelse(d$grade %in% c("G3", "G4"), "G3-4",
                      ifelse(d$grade %in% c("G1", "G2"), "G1-2", NA)), levels = c("G1-2", "G3-4"))
d$group     <- factor(d$RLRS_group, levels = c("High", "Low"))
log_msg("TCGA-KIRC 纳入 ", nrow(d), " 例（OS 事件 ", sum(d$OS), " 例）")

# =============================================================================
# 1. 单因素与多因素 Cox（TCGA-KIRC）
# =============================================================================
vars <- c(RLRS_z = "MTRS (per SD)", age10 = "Age (per 10 years)", gender = "Gender (male vs female)",
          stage_bin = "AJCC Stage (III-IV vs I-II)", grade_bin = "Fuhrman Grade (G3-4 vs G1-2)")
uni <- do.call(rbind, lapply(names(vars), function(v) {
  r <- cox_table(survival::coxph(stats::as.formula(paste("survival::Surv(OS.time, OS) ~", v)), data = d))
  data.frame(variable = vars[[v]], r[, -1])
}))
fit_multi <- survival::coxph(survival::Surv(OS.time, OS) ~ RLRS_z + age10 + gender + stage_bin + grade_bin, data = d)
multi <- data.frame(variable = unname(vars), cox_table(fit_multi)[, -1])
cox_all <- rbind(data.frame(model = "Univariate", uni), data.frame(model = "Multivariate", multi))
print(cox_all, row.names = FALSE)
utils::write.csv(cox_all, file.path(OUT, "TCGA_uni_multi_Cox.csv"), row.names = FALSE)

# 外部队列多因素 Cox
ext_rows <- list()
for (n in mdl$valid) {
  cl <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$clin
  e <- merge(cl, rlrs[rlrs$cohort == n, c("ID", "RLRS_z")], by = "ID")
  e <- e[!is.na(e$OS.time) & e$OS.time > 0 & !is.na(e$OS), ]
  if ("stage" %in% colnames(e)) {
    st_raw <- toupper(as.character(e$stage))
    e$stage_bin <- factor(ifelse(grepl("III|IV|3|4", st_raw), "III-IV",
                          ifelse(grepl("I|II|1|2", st_raw), "I-II", NA)), levels = c("I-II", "III-IV"))
  }
  if ("grade" %in% colnames(e)) {
    gr_raw <- toupper(as.character(e$grade))
    e$grade_bin <- factor(ifelse(grepl("3|4|HIGH", gr_raw), "G3-4",
                          ifelse(grepl("1|2|LOW", gr_raw), "G1-2", NA)), levels = c("G1-2", "G3-4"))
  }
  if ("age" %in% colnames(e)) e$age10 <- suppressWarnings(as.numeric(e$age)) / 10
  cov <- intersect(c("age10", "gender", "stage_bin", "grade_bin"), colnames(e))
  cov <- cov[vapply(cov, function(v) sum(!is.na(e[[v]])) >= 15 && length(unique(stats::na.omit(e[[v]]))) > 1, logical(1))]
  e_cc <- e[stats::complete.cases(e[, c("OS.time", "OS", "RLRS_z", cov), drop = FALSE]), ]
  if (nrow(e_cc) < 15) next
  f <- stats::as.formula(paste("survival::Surv(OS.time, OS) ~", paste(c("RLRS_z", cov), collapse = " + ")))
  m_full <- survival::coxph(f, data = e_cc)
  r <- cox_table(m_full)
  m_clin <- if (length(cov)) survival::coxph(stats::as.formula(paste("survival::Surv(OS.time, OS) ~",
                                                                     paste(cov, collapse = " + "))), data = e_cc) else NULL
  ext_rows[[n]] <- data.frame(
    cohort = n, n = nrow(e_cc), events = sum(e_cc$OS),
    covariates = if (length(cov)) paste(cov, collapse = "+") else "unadjusted",
    r,
    LRT_p_add_RLRS = if (is.null(m_clin)) NA_real_ else stats::anova(m_clin, m_full)[2, "Pr(>|Chi|)"],
    C_clinical = if (is.null(m_clin)) NA_real_ else survival::concordance(m_clin)$concordance,
    C_clinical_RLRS = survival::concordance(m_full)$concordance
  )
}
if (length(ext_rows)) {
  ext_tab <- do.call(rbind, ext_rows); print(ext_tab, row.names = FALSE, digits = 3)
  utils::write.csv(ext_tab, file.path(OUT, "External_multivariate_Cox.csv"), row.names = FALSE)
}

# =============================================================================
# 2. 列线图（rms）与 1/3/5 年校准曲线（1000 次 bootstrap）
# =============================================================================
dn <- d[stats::complete.cases(d[, c("time_y", "OS", NOMO_VARS)]), ]
dd <- rms::datadist(dn); options(datadist = "dd")
f_nomo <- stats::as.formula(paste("Surv(time_y, OS) ~", paste(NOMO_VARS, collapse = " + ")))
fit_cph <- rms::cph(f_nomo, data = dn, x = TRUE, y = TRUE, surv = TRUE)
S <- rms::Survival(fit_cph)
nom <- rms::nomogram(fit_cph, lp = FALSE,
                     fun = list(function(x) S(1, x), function(x) S(3, x), function(x) S(5, x)),
                     funlabel = c("1-year OS", "3-year OS", "5-year OS"),
                     fun.at = c(0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.1))
grDevices::pdf(file.path(OUT, "Fig7_nomogram.pdf"), 9, 6); plot(nom); grDevices::dev.off()
saveRDS(nom, file.path(DIR_PRO, "nomogram_object.rds"))
f_nomo_cox <- stats::as.formula(paste("survival::Surv(time_y, OS) ~", paste(NOMO_VARS, collapse = " + ")))
utils::write.csv(cox_table(survival::coxph(f_nomo_cox, data = dn)), file.path(OUT, "nomogram_model_Cox.csv"), row.names = FALSE)

grDevices::pdf(file.path(OUT, "Fig7_calibration.pdf"), 5, 5)
cols <- c("#E41A1C", "#377EB8", "#4DAF4A")
cal_rows <- list()
for (k in seq_along(PARAM$surv_years)) {
  t <- PARAM$surv_years[k]
  fit_t <- rms::cph(f_nomo, data = dn, x = TRUE, y = TRUE, surv = TRUE, time.inc = t)
  set.seed(PARAM$seed)
  cal <- rms::calibrate(fit_t, cmethod = "KM", method = "boot", u = t, m = floor(nrow(dn) / 4), B = 1000)
  cal_rows[[k]] <- data.frame(year = t, as.data.frame(unclass(cal)[, , drop = FALSE]), check.names = FALSE)
  plot(cal, lwd = 2, lty = 1, errbar.col = cols[k], col = cols[k], add = k > 1,
       xlim = c(0, 1), ylim = c(0, 1), subtitles = FALSE,
       xlab = "Nomogram-predicted OS", ylab = "Observed OS (Kaplan-Meier)")
}
graphics::abline(0, 1, lty = 2, col = "grey50")
grDevices::dev.off()
utils::write.csv(do.call(rbind, cal_rows), file.path(OUT, "calibration_bootstrap.csv"), row.names = FALSE)

# =============================================================================
# 3. C-index、时间依赖 ROC 与 DCA（Nomogram vs MTRS vs Stage）
# =============================================================================
fit_rlrs  <- survival::coxph(survival::Surv(time_y, OS) ~ RLRS, data = dn)
fit_stage <- survival::coxph(survival::Surv(time_y, OS) ~ stage_bin + grade_bin, data = dn)
fit_nomo  <- survival::coxph(f_nomo_cox, data = dn)
lp <- data.frame(Nomogram = stats::predict(fit_nomo, type = "lp"),
                 RLRS = stats::predict(fit_rlrs, type = "lp"),
                 Stage = stats::predict(fit_stage, type = "lp"))
cidx <- vapply(lp, function(x) survival::concordance(survival::Surv(dn$time_y, dn$OS) ~ x, reverse = TRUE)$concordance, numeric(1))
print(round(cidx, 3))
utils::write.csv(data.frame(model = names(cidx), Cindex = cidx), file.path(OUT, "Cindex_nomogram_RLRS_stage.csv"), row.names = FALSE)

yrs <- PARAM$surv_years[PARAM$surv_years < max(dn$time_y[dn$OS == 1])]
aucs <- sapply(lp, function(x) timeROC::timeROC(T = dn$time_y, delta = dn$OS, marker = x, cause = 1,
                                                times = yrs, iid = FALSE)$AUC)
utils::write.csv(data.frame(year = yrs, matrix(aucs, nrow = length(yrs), dimnames = list(NULL, colnames(lp)))),
                 file.path(OUT, "timeROC_nomogram_RLRS_stage.csv"), row.names = FALSE)

risk_at <- function(fit, t) {
  sf <- survival::survfit(fit, newdata = dn)
  1 - as.numeric(summary(sf, times = t, extend = TRUE)$surv)
}
dca_obj <- lapply(c(1, 3, 5)[c(1, 3, 5) < max(dn$time_y)], function(t) {
  dd_t <- data.frame(time_y = dn$time_y, OS = dn$OS, Nomogram = risk_at(fit_nomo, t),
                     RLRS = risk_at(fit_rlrs, t), Stage = risk_at(fit_stage, t))
  list(t = t, dca = dcurves::dca(Surv(time_y, OS) ~ Nomogram + RLRS + Stage, data = dd_t, time = t,
                                 thresholds = seq(0, 0.7, 0.01)))
})
utils::write.csv(do.call(rbind, lapply(dca_obj, function(o) data.frame(year = o$t, as.data.frame(o$dca$dca)))),
                 file.path(OUT, "DCA_net_benefit.csv"), row.names = FALSE)

# 4. 临床亚组 KM
sub_km <- list()
for (v in c("stage_bin", "grade_bin", "age_group", "gender")) {
  for (lv in levels(d[[v]])) {
    s <- d[!is.na(d[[v]]) & d[[v]] == lv, ]
    if (nrow(s) < 25 || length(unique(s$group)) < 2) next
    sub_km[[paste(v, lv)]] <- km_plot(s, title = paste0(v, ": ", lv, " (n = ", nrow(s), ")"))
  }
}
save_print(sub_km, file.path(OUT, "Fig7_subgroup_KM.pdf"), 6, 6)
save_session(STEP)
log_msg("07_clinical_nomogram 完成")
