# =============================================================================
# 04b_proliferation_overlap.R —— 检验 MT_core / MTI 与增殖（Prolif）及 EMT 的重叠与混杂效应
# =============================================================================
source("R/00_config.R")
need(c("GSVA", "msigdbr", "survival"))
OUT <- res_dir("04_bulk_score")
GS <- readRDS(file.path(DIR_GS, "MTRG_genesets.rds"))
H <- get_msig("H")
prolif_genes <- union(H[["HALLMARK_E2F_TARGETS"]], H[["HALLMARK_G2M_CHECKPOINT"]])
emt_genes    <- H[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]]
ov_prolif <- intersect(GS$MT_core, prolif_genes)
ov_emt    <- intersect(GS$MT_core, emt_genes)
log_msg("MT_core ", length(GS$MT_core), " 个基因中：属于 E2F/G2M 有 ", length(ov_prolif),
        " 个；属于 Hallmark EMT 有 ", length(ov_emt), " 个")

tcga <- readRDS(file.path(DIR_PRO, "TCGA_KIRC.rds"))
tum  <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
t1   <- tum[tum$group == "Tumor", ]
ex   <- tcga$expr[, t1$ID]
s_sub <- run_ssgsea(ex, list(
  MT_core_noProlif = setdiff(GS$MT_core, prolif_genes),
  MT_core_noEMT    = setdiff(GS$MT_core, emt_genes)
))
t1$MT_core_noProlif <- as.numeric(s_sub["MT_core_noProlif", t1$ID])
t1$MT_core_noEMT    <- as.numeric(s_sub["MT_core_noEMT", t1$ID])

ct <- function(a, b) stats::cor.test(a, b, method = "spearman", exact = FALSE)
r0 <- ct(t1$MT_core, t1$Prolif)
r1 <- ct(t1$MT_core_noProlif, t1$Prolif)
r2 <- ct(t1$MT_core, t1$EMT)
r3 <- ct(t1$MT_core_noEMT, t1$EMT)

d <- t1[!is.na(t1$OS.time) & t1$OS.time > 0 & !is.na(t1$OS), ]
d$z_np <- as.numeric(scale(d$MT_core_noProlif))
d$z_ne <- as.numeric(scale(d$MT_core_noEMT))
u_np <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ z_np, data = d))
a_np <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ z_np + age + stage_bin + Prolif, data = d))
u_ne <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ z_ne, data = d))
a_ne <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ z_ne + age + stage_bin + EMT, data = d))

summ <- data.frame(
  item = c("MT_core genes", "MT_core genes in Hallmark E2F/G2M", "MT_core genes in Hallmark EMT",
           "rho MT_core vs Prolif", "rho MT_core_noProlif vs Prolif",
           "rho MT_core vs EMT", "rho MT_core_noEMT vs EMT",
           "OS HR per SD MT_core_noProlif (unadjusted)", "OS HR per SD MT_core_noProlif (adjusted age, stage, Prolif)",
           "OS HR per SD MT_core_noEMT (unadjusted)", "OS HR per SD MT_core_noEMT (adjusted age, stage, EMT)"),
  value = c(length(GS$MT_core), length(ov_prolif), length(ov_emt),
            round(r0$estimate, 3), round(r1$estimate, 3), round(r2$estimate, 3), round(r3$estimate, 3),
            sprintf("%.2f (%.2f-%.2f)", u_np$HR[1], u_np$lower[1], u_np$upper[1]),
            sprintf("%.2f (%.2f-%.2f)", a_np$HR[a_np$term == "z_np"], a_np$lower[a_np$term == "z_np"], a_np$upper[a_np$term == "z_np"]),
            sprintf("%.2f (%.2f-%.2f)", u_ne$HR[1], u_ne$lower[1], u_ne$upper[1]),
            sprintf("%.2f (%.2f-%.2f)", a_ne$HR[a_ne$term == "z_ne"], a_ne$lower[a_ne$term == "z_ne"], a_ne$upper[a_ne$term == "z_ne"])),
  p = c(NA, NA, NA, r0$p.value, r1$p.value, r2$p.value, r3$p.value,
        u_np$p[1], a_np$p[a_np$term == "z_np"], u_ne$p[1], a_ne$p[a_ne$term == "z_ne"])
)
print(summ, row.names = FALSE)
utils::write.csv(summ, file.path(OUT, "proliferation_overlap_summary.csv"), row.names = FALSE)
utils::write.csv(data.frame(gene = sort(ov_prolif), set = "E2F_G2M"), file.path(OUT, "proliferation_overlap_genes.csv"), row.names = FALSE)

g <- intersect(GS$MT_core, rownames(ex))
rho_p <- vapply(g, function(x) suppressWarnings(stats::cor(as.numeric(ex[x, ]), t1$Prolif, method = "spearman")), numeric(1))
rho_e <- vapply(g, function(x) suppressWarnings(stats::cor(as.numeric(ex[x, ]), t1$EMT, method = "spearman")), numeric(1))
pg <- data.frame(gene = g, rho_Prolif = rho_p, rho_EMT = rho_e, in_E2F_G2M = g %in% prolif_genes, in_EMT = g %in% emt_genes)
utils::write.csv(pg[order(-pg$rho_Prolif), ], file.path(OUT, "proliferation_overlap_per_gene_rho.csv"), row.names = FALSE)

sc_all <- readRDS(file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds"))
rc <- do.call(rbind, lapply(names(sc_all), function(n) {
  x <- sc_all[[n]]
  if (n == "TCGA_KIRC") x <- x[x$ID %in% t1$ID, ]
  r <- ct(x$MT_core, x$Prolif); r2 <- ct(x$MTI, x$Prolif); r3 <- ct(x$MT_core, x$EMT)
  data.frame(cohort = n, n = nrow(x),
             rho_MTcore_Prolif = unname(r$estimate), p_MTcore_Prolif = r$p.value,
             rho_MTI_Prolif = unname(r2$estimate), p_MTI_Prolif = r2$p.value,
             rho_MTcore_EMT = unname(r3$estimate), p_MTcore_EMT = r3$p.value)
}))
utils::write.csv(rc, file.path(OUT, "proliferation_correlation_by_cohort.csv"), row.names = FALSE)
log_msg("04b_proliferation_overlap 完成")
