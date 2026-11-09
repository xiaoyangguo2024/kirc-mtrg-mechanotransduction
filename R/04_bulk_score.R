# =============================================================================
# 04_bulk_score.R —— Bulk 队列 ssGSEA 质膜张力评分；
#   TCGA-KIRC 中与肿瘤/正常（配对及非配对）、AJCC 分期、Fuhrman 分级、
#   ccRCC 关键驱动突变（VHL / PBRM1 / SETD2 / BAP1 / KDM5C / MTOR / TP53）、TMB、
#   血管生成 / 缺氧 / EMT / IFN-γ / 髓系炎症 / 增殖的关联与生存分析
# =============================================================================
source("R/00_config.R")
need(c("GSVA", "msigdbr", "survival", "survminer", "ggplot2", "ggpubr"))
STEP <- "04_bulk_score"; OUT <- res_dir(STEP)
GS <- readRDS(file.path(DIR_GS, "MTRG_genesets.rds"))

# 自动发现 data/processed 下所有已整理的 bulk 队列
ALL_RDS <- sub("\\.rds$", "", list.files(DIR_PRO, pattern = "\\.rds$"))
COHORTS <- intersect(c("TCGA_KIRC", "GSE29609", "GSE22541", "GSE167573", "GSE53757", "GSE73731",
                       "GSE67501", "IMvigor210", "GSE176307"), ALL_RDS)
log_msg("纳入 04_bulk_score 的队列：", paste(COHORTS, collapse = ", "))
coh <- lapply(stats::setNames(COHORTS, COHORTS), function(n) readRDS(file.path(DIR_PRO, paste0(n, ".rds"))))

# =============================================================================
# 1. 基因集在各队列中的覆盖率
# =============================================================================
cov <- do.call(rbind, lapply(COHORTS, function(n) data.frame(
  cohort = n, geneset = names(GS), n_genes = lengths(GS),
  n_found = vapply(GS, function(g) sum(g %in% rownames(coh[[n]]$expr)), integer(1))
)))
cov$pct <- round(100 * cov$n_found / cov$n_genes, 1)
utils::write.csv(cov, file.path(OUT, "geneset_coverage_by_cohort.csv"), row.names = FALSE)
print(stats::xtabs(pct ~ cohort + geneset, cov))

# =============================================================================
# 2. ssGSEA：质膜张力基因集 + 肾癌微环境/机制通路（血管生成、缺氧、EMT、IFN、cGAS-STING、增殖）
# =============================================================================
H  <- get_msig("H")
C2 <- get_msig("C2")
pick <- function(lst, nm) { if (!nm %in% names(lst)) log_msg("【核对】MSigDB 中找不到 ", nm); lst[[nm]] }
SIG <- c(GS, list(
  Hypoxia         = pick(H, "HALLMARK_HYPOXIA"),
  Angiogenesis    = pick(H, "HALLMARK_ANGIOGENESIS"),
  EMT             = pick(H, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  TGF_beta        = pick(H, "HALLMARK_TGF_BETA_SIGNALING"),
  Glycolysis      = pick(H, "HALLMARK_GLYCOLYSIS"),
  IFN_alpha       = pick(H, "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
  IFN_gamma       = pick(H, "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  Inflammatory    = pick(H, "HALLMARK_INFLAMMATORY_RESPONSE"),
  IL6_JAK_STAT3   = pick(H, "HALLMARK_IL6_JAK_STAT3_SIGNALING"),
  E2F             = pick(H, "HALLMARK_E2F_TARGETS"),
  G2M             = pick(H, "HALLMARK_G2M_CHECKPOINT"),
  DNA_repair      = pick(H, "HALLMARK_DNA_REPAIR"),
  cGAS_STING_KEGG = pick(C2, "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY"),
  GEP18           = GEP18,
  IMmotion_Angio  = ANGIO_SIG,
  IMmotion_Myeloid = MYELOID_INFLAM
))
SIG <- SIG[lengths(SIG) > 0]

f_sc <- file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds")
if (!file.exists(f_sc)) {
  scores <- lapply(coh, function(o) {
    s <- as.data.frame(t(run_ssgsea(o$expr, SIG)))
    s$MTI    <- as.numeric(scale(s$MT_stiff) - scale(s$MT_buffer))
    s$Prolif <- as.numeric(scale(s$E2F) + scale(s$G2M)) / 2
    # 兼容别名
    s$RL_core <- s$MT_core
    s$RL_all  <- s$MT_all
    s$RHI     <- s$MTI
    s$ID <- rownames(s)
    s
  })
  saveRDS(scores, f_sc)
}
scores <- readRDS(f_sc)

# =============================================================================
# 3. TCGA-KIRC：肿瘤 vs 正常（全部 + 72 对配对）、分期（Stage I-IV）、Fuhrman 分级（G1-G4）
# =============================================================================
tc <- merge(coh$TCGA_KIRC$clin, scores$TCGA_KIRC, by = "ID")
tc$group <- factor(tc$group, levels = c("Normal", "Tumor"))
tum <- tc[tc$group == "Tumor", ]
tum$stage_bin <- factor(ifelse(tum$stage %in% c("III", "IV"), "III-IV",
                        ifelse(tum$stage %in% c("I", "II"), "I-II", NA)), levels = c("I-II", "III-IV"))
tum$grade_bin <- factor(ifelse(tum$grade %in% c("G3", "G4"), "G3-4",
                        ifelse(tum$grade %in% c("G1", "G2"), "G1-2", NA)), levels = c("G1-2", "G3-4"))

plots <- list()
for (s in c("MT_core", "MTI")) {
  plots[[paste0("TvN_", s)]] <- ggpubr::ggboxplot(tc, "group", s, color = "group", add = "jitter",
                                                  palette = c("#2C7FB8", "#D7301F")) +
    ggpubr::stat_compare_means(method = "wilcox.test")
  pat <- intersect(tc$patient[tc$group == "Tumor"], tc$patient[tc$group == "Normal"])
  plots[[paste0("paired_", s)]] <- ggpubr::ggpaired(tc[tc$patient %in% pat, ], x = "group", y = s, id = "patient",
                                                   color = "group", line.color = "grey70",
                                                   palette = c("#2C7FB8", "#D7301F")) +
    ggpubr::stat_compare_means(paired = TRUE, method = "wilcox.test") +
    ggplot2::labs(title = paste0("Paired KIRC (n = ", length(pat), ")"))
  plots[[paste0("stage_", s)]] <- ggpubr::ggboxplot(tum[!is.na(tum$stage), ], "stage", s, color = "stage", add = "jitter") +
    ggpubr::stat_compare_means(method = "kruskal.test")
  plots[[paste0("grade_", s)]] <- ggpubr::ggboxplot(tum[!is.na(tum$grade), ], "grade", s, color = "grade", add = "jitter") +
    ggpubr::stat_compare_means(method = "kruskal.test")
}

# 同步检验外部配对/分期队列 GSE53757（72 对配对肿瘤/正常）与 GSE73731（265 例分级/分期）
ext_clin_summary <- list()
if ("GSE53757" %in% names(coh)) {
  d53 <- merge(coh$GSE53757$clin, scores$GSE53757, by = "ID")
  p_tvn <- tryCatch(stats::wilcox.test(MT_core ~ group, data = d53)$p.value, error = function(e) NA_real_)
  p_st  <- tryCatch(stats::kruskal.test(MT_core ~ stage, data = d53[d53$group == "Tumor", ])$p.value, error = function(e) NA_real_)
  ext_clin_summary$GSE53757 <- data.frame(cohort = "GSE53757", test = c("Tumor_vs_Normal", "Stage_I_to_IV"),
                                          p_MT_core = c(p_tvn, p_st))
}
if ("GSE73731" %in% names(coh)) {
  d73 <- merge(coh$GSE73731$clin, scores$GSE73731, by = "ID")
  p_gr <- tryCatch(stats::kruskal.test(MT_core ~ grade, data = d73)$p.value, error = function(e) NA_real_)
  p_st <- tryCatch(stats::kruskal.test(MT_core ~ stage, data = d73)$p.value, error = function(e) NA_real_)
  ext_clin_summary$GSE73731 <- data.frame(cohort = "GSE73731", test = c("Grade_G1_to_G4", "Stage_I_to_IV"),
                                          p_MT_core = c(p_gr, p_st))
}
if (length(ext_clin_summary)) {
  utils::write.csv(do.call(rbind, ext_clin_summary), file.path(OUT, "External_stage_grade_validation.csv"), row.names = FALSE)
}

# =============================================================================
# 4. 体细胞突变：VHL / PBRM1 / SETD2 / BAP1 / KDM5C / MTOR / TP53 / PTEN / ATM 与 TMB
# =============================================================================
maf <- readRDS(file.path(DIR_PRO, "TCGA_KIRC_maf.rds"))
maf$sample <- substr(maf$Tumor_Sample_Barcode, 1, 15)
NONSYN <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation",
            "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
ns <- maf[maf$Variant_Classification %in% NONSYN, ]
mut_samples <- unique(maf$sample)
tum$TMB <- ifelse(tum$ID %in% mut_samples,
                  as.numeric(table(factor(ns$sample, levels = mut_samples))[tum$ID]) / 38, NA_real_)

MUT_GENES <- c("VHL", "PBRM1", "SETD2", "BAP1", "KDM5C", "MTOR", "TP53", "PTEN", "ATM", "ARID1A", "PIK3CA", "SMARCA4")
for (g in MUT_GENES) {
  tum[[paste0("mut_", g)]] <- ifelse(tum$ID %in% mut_samples,
                                     ifelse(tum$ID %in% ns$sample[ns$Hugo_Symbol == g], "Mut", "WT"), NA_character_)
}
mut_freq <- data.frame(
  gene = MUT_GENES,
  n_mut = vapply(MUT_GENES, function(g) sum(tum[[paste0("mut_", g)]] == "Mut", na.rm = TRUE), integer(1)),
  n_tested = sum(tum$ID %in% mut_samples)
)
mut_freq$freq <- round(mut_freq$n_mut / mut_freq$n_tested, 3)
print(mut_freq)
utils::write.csv(mut_freq, file.path(OUT, "TCGA_mutation_frequency.csv"), row.names = FALSE)

mut_long <- do.call(rbind, lapply(MUT_GENES, function(g)
  data.frame(gene = g, status = tum[[paste0("mut_", g)]], MT_core = tum$MT_core, MTI = tum$MTI)))
mut_long <- mut_long[!is.na(mut_long$status), ]
mut_test <- do.call(rbind, lapply(MUT_GENES, function(g) {
  d <- mut_long[mut_long$gene == g, ]
  wp <- function(f) tryCatch(stats::wilcox.test(f, d)$p.value, error = function(e) NA_real_)
  data.frame(
    gene = g, n_mut = sum(d$status == "Mut"),
    median_diff_MT_core = stats::median(d$MT_core[d$status == "Mut"]) - stats::median(d$MT_core[d$status == "WT"]),
    p_MT_core = wp(MT_core ~ status), p_MTI = wp(MTI ~ status),
    p_RL_core = wp(MT_core ~ status), p_RHI = wp(MTI ~ status)
  )
}))
mut_test$FDR_MT_core <- stats::p.adjust(mut_test$p_MT_core, "BH")
mut_test$FDR_MTI     <- stats::p.adjust(mut_test$p_MTI, "BH")
mut_test$FDR_RL_core <- mut_test$FDR_MT_core
mut_test$FDR_RHI     <- mut_test$FDR_MTI
utils::write.csv(mut_test, file.path(OUT, "MT_score_by_mutation.csv"), row.names = FALSE)
utils::write.csv(mut_test, file.path(OUT, "Rloop_score_by_mutation.csv"), row.names = FALSE)

# =============================================================================
# 5. 通路相关性矩阵
# =============================================================================
cor_vars <- intersect(
  c("MT_core", "MT_all", "MT_cortex", "MT_BAR", "MT_channel", "MT_ECM", "MTI",
    "Hypoxia", "Angiogenesis", "IMmotion_Angio", "EMT", "TGF_beta", "IFN_gamma", "GEP18",
    "IMmotion_Myeloid", "IL6_JAK_STAT3", "cGAS_STING_KEGG", "E2F", "G2M", "Prolif"),
  colnames(tum)
)
cm_r <- stats::cor(tum[, cor_vars], method = "spearman", use = "pairwise.complete.obs")
utils::write.csv(cm_r, file.path(OUT, "signature_spearman_TCGA.csv"))
save_print(plots, file.path(OUT, "Fig4_TCGA_associations.pdf"), 7, 6)

# =============================================================================
# 6. TCGA-KIRC 生存分析（OS / DSS / PFI）：未校正 vs 校正年龄、分期、分级、增殖评分
# =============================================================================
EPS <- list(OS = c("OS.time", "OS"), DSS = c("DSS.time", "DSS"), PFI = c("PFI.time", "PFI"))
km_list <- list(); cox_rows <- list()
for (s in c("MT_core", "MTI", "MT_all", "MT_cortex", "MT_BAR", "MT_channel", "MT_ECM")) {
  for (ep in names(EPS)) {
    tt <- EPS[[ep]][1]; ee <- EPS[[ep]][2]
    d <- tum[!is.na(tum[[tt]]) & tum[[tt]] > 0 & !is.na(tum[[ee]]), ]
    d$group <- median_group(d[[s]]); d$z <- as.numeric(scale(d[[s]]))
    km_list[[paste(s, ep)]] <- km_plot(d, tt, ee, "group", title = paste("TCGA-KIRC", s, ep))
    f_uni <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ z", tt, ee))
    f_adj <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ z + age + stage_bin + Prolif", tt, ee))
    u <- cox_table(survival::coxph(f_uni, data = d))
    a <- cox_table(survival::coxph(f_adj, data = d))
    cox_rows[[paste(s, ep)]] <- data.frame(
      score = s, endpoint = ep, n = nrow(d), events = sum(d[[ee]] == 1),
      HR_uni = u$HR[1], lower_uni = u$lower[1], upper_uni = u$upper[1], p_uni = u$p[1],
      HR_adj = a$HR[a$term == "z"], lower_adj = a$lower[a$term == "z"],
      upper_adj = a$upper[a$term == "z"], p_adj = a$p[a$term == "z"]
    )
  }
}
cox_tcga <- do.call(rbind, cox_rows)
print(cox_tcga[cox_tcga$endpoint == "OS", ], row.names = FALSE)
utils::write.csv(cox_tcga, file.path(OUT, "TCGA_Cox_MT_scores.csv"), row.names = FALSE)
utils::write.csv(cox_tcga, file.path(OUT, "TCGA_Cox_Rloop_scores.csv"), row.names = FALSE)
save_print(km_list, file.path(OUT, "Fig4_TCGA_KM.pdf"), 6, 6)

# 外部有随访数据的队列
ext <- list()
for (n in setdiff(COHORTS, "TCGA_KIRC")) {
  d <- merge(coh[[n]]$clin, scores[[n]], by = "ID")
  d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
  if (nrow(d) < 15 || sum(d$OS == 1, na.rm = TRUE) < 5) next
  d$z <- as.numeric(scale(d$MT_core))
  u <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ z, data = d))
  ext[[n]] <- data.frame(cohort = n, n = nrow(d), events = sum(d$OS == 1),
                         HR_perSD = u$HR, lower = u$lower, upper = u$upper, p = u$p)
}
if (length(ext)) {
  ext_df <- do.call(rbind, ext)
  utils::write.csv(ext_df, file.path(OUT, "External_Cox_MT_core.csv"), row.names = FALSE)
  utils::write.csv(ext_df, file.path(OUT, "External_Cox_RL_core.csv"), row.names = FALSE)
}

saveRDS(tum, file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
save_session(STEP)
log_msg("04_bulk_score 完成")
