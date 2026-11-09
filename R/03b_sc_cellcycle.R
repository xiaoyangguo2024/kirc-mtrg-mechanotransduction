# =============================================================================
# 03b_sc_cellcycle.R —— 单细胞层面检验：恶性细胞膜张力评分（MT_core）是否受细胞周期混杂
#   1. 每个样本内 MT_core_UCell 与 S.Score / G2M.Score 的 Spearman 相关
#   2. MT 高/低组中处于 S/G2M 期细胞比例（样本为单位，配对 Wilcoxon）
#   3. 仅保留 G1 期恶性细胞，样本内重新按中位数分组，检验 Hallmark 差异是否稳健存在
#   4. 敏感性分析：仅保留 copyKAT 判定为 aneuploid 的肿瘤上皮细胞（results/03_scRNA/sensitivity_copykat_rule）
# =============================================================================
source("R/00_config.R")
need(c("Seurat", "UCell", "edgeR", "ggplot2", "ggpubr"))
OUT <- res_dir("03_scRNA")
mal <- readRDS(file.path(DIR_PRO, "scRNA_malignant.rds"))
cc <- Seurat::cc.genes.updated.2019
mal <- CellCycleScoring(mal, s.features = intersect(cc$s.genes, rownames(mal)),
                        g2m.features = intersect(cc$g2m.genes, rownames(mal)), set.ident = FALSE)
md <- mal@meta.data
md$cycling <- md$Phase %in% c("S", "G2M")
log_msg("恶性细胞 ", nrow(md), " 个；S/G2M 期比例 ", round(100 * mean(md$cycling), 1), "%")

per_s <- do.call(rbind, lapply(split(md, md$sample), function(x) data.frame(
  sample = x$sample[1], n_cells = nrow(x),
  rho_MTcore_S = stats::cor(x$MT_core_UCell, x$S.Score, method = "spearman"),
  rho_MTcore_G2M = stats::cor(x$MT_core_UCell, x$G2M.Score, method = "spearman"),
  rho_RLcore_S = stats::cor(x$MT_core_UCell, x$S.Score, method = "spearman"),
  rho_RLcore_G2M = stats::cor(x$MT_core_UCell, x$G2M.Score, method = "spearman"),
  pct_cycling_High = 100 * mean(x$cycling[x$MT_group == "High"]),
  pct_cycling_Low = 100 * mean(x$cycling[x$MT_group == "Low"])
)))
utils::write.csv(per_s, file.path(OUT, "cellcycle_by_sample.csv"), row.names = FALSE)
p_cyc <- tryCatch(stats::wilcox.test(per_s$pct_cycling_High, per_s$pct_cycling_Low, paired = TRUE)$p.value, error = function(e) NA_real_)
log_msg("样本内 MT_core 与 G2M 评分相关中位数 rho = ", round(stats::median(per_s$rho_MTcore_G2M, na.rm = TRUE), 2),
        "；S/G2M 期比例 High vs Low P = ", signif(p_cyc, 3))

# 仅用 G1 期细胞重新检验 Hallmark（每个肿瘤样本至少 20 个 G1 期恶性细胞以保留全部 7 例肿瘤配对）
g1 <- md[md$Phase == "G1", ]
keep <- names(which(table(g1$sample) >= 20))
if (length(keep) < 3) keep <- names(which(table(g1$sample) >= 10))
g1 <- g1[g1$sample %in% keep, ]
g1$G1_group <- stats::ave(g1$MT_core_UCell, g1$sample, FUN = function(v) ifelse(v >= stats::median(v), "High", "Low"))
hcols <- grep("^HALLMARK_.*_UCell$", colnames(g1), value = TRUE)
if (!length(hcols)) {
  H <- get_msig("H")
  mal <- UCell::AddModuleScore_UCell(mal, features = H, name = "_UCell", BPPARAM = BiocParallel::SerialParam())
  hcols <- paste0(names(H), "_UCell")
  g1 <- cbind(g1, mal@meta.data[rownames(g1), hcols])
}
hm <- stats::aggregate(g1[, hcols], by = list(sample = g1$sample, group = g1$G1_group), FUN = mean)
hall_g1 <- do.call(rbind, lapply(hcols, function(h) {
  w <- stats::reshape(hm[, c("sample", "group", h)], idvar = "sample", timevar = "group", direction = "wide")
  d <- w[[paste0(h, ".High")]] - w[[paste0(h, ".Low")]]
  data.frame(pathway = sub("^HALLMARK_", "", sub("_UCell$", "", h)), mean_diff = mean(d, na.rm = TRUE),
             p = tryCatch(stats::wilcox.test(d, mu = 0)$p.value, error = function(e) NA_real_))
}))
hall_g1$FDR <- stats::p.adjust(hall_g1$p, "BH")
hall_g1 <- hall_g1[order(hall_g1$mean_diff), ]
utils::write.csv(hall_g1, file.path(OUT, "Hallmark_High_vs_Low_G1only.csv"), row.names = FALSE)

f_all <- file.path(OUT, "Hallmark_High_vs_Low_paired.csv")
if (file.exists(f_all)) {
  ha <- utils::read.csv(f_all)
  both <- merge(ha[, c("pathway", "mean_diff", "FDR")], hall_g1[, c("pathway", "mean_diff", "FDR")],
                by = "pathway", suffixes = c("_all", "_G1"))
  utils::write.csv(both, file.path(OUT, "Hallmark_all_vs_G1only.csv"), row.names = FALSE)
  log_msg("G1 期细胞中仍 FDR < 0.05 的 Hallmark 通路数：", sum(both$FDR_G1 < 0.05, na.rm = TRUE),
          "（全部恶性细胞中：", sum(both$FDR_all < 0.05, na.rm = TRUE), "）")
}

# =============================================================================
# 4. 敏感性分析：仅保留 copyKAT 判定为 aneuploid 的肿瘤上皮细胞（Figure S2）
# =============================================================================
OUT_SUB <- file.path(OUT, "sensitivity_copykat_rule")
dir.create(OUT_SUB, recursive = TRUE, showWarnings = FALSE)
mal_an <- subset(mal, subset = copykat == "aneuploid")
n_an <- table(mal_an$sample)
keep_an <- names(n_an)[n_an >= 20]
mal_an <- subset(mal_an, subset = sample %in% keep_an)
grp_an <- character(ncol(mal_an))
for (s in keep_an) {
  i <- mal_an$sample == s
  grp_an[i] <- ifelse(mal_an$MT_core_UCell[i] >= stats::median(mal_an$MT_core_UCell[i]), "High", "Low")
}
mal_an$MT_group <- factor(grp_an, levels = c("Low", "High"))
md_an <- mal_an@meta.data
md_an$cycling <- md_an$Phase %in% c("S", "G2M")

pb_an <- AggregateExpression(mal_an, assays = "RNA", group.by = c("sample", "MT_group"), return.seurat = FALSE)$RNA
pb_an <- as.matrix(pb_an)
pb_smp <- sub("_[^_]+$", "", colnames(pb_an))
pb_grp <- factor(sub("^.*_", "", colnames(pb_an)), levels = c("Low", "High"))
y_an <- edgeR::DGEList(pb_an, group = pb_grp)
y_an <- y_an[edgeR::filterByExpr(y_an, group = pb_grp), , keep.lib.sizes = FALSE]
y_an <- edgeR::calcNormFactors(y_an)
des_an <- stats::model.matrix(~ factor(pb_smp) + pb_grp)
y_an <- edgeR::estimateDisp(y_an, des_an)
fit_an <- edgeR::glmQLFit(y_an, des_an)
deg_an <- edgeR::topTags(edgeR::glmQLFTest(fit_an, coef = "pb_grpHigh"), n = Inf)$table
deg_an$gene <- rownames(deg_an)
utils::write.csv(deg_an, file.path(OUT_SUB, "pseudobulk_DE_High_vs_Low.csv"), row.names = FALSE)

per_an <- do.call(rbind, lapply(split(md_an, md_an$sample), function(x) data.frame(
  sample = x$sample[1], n_cells = nrow(x),
  rho_MTcore_S = stats::cor(x$MT_core_UCell, x$S.Score, method = "spearman"),
  rho_MTcore_G2M = stats::cor(x$MT_core_UCell, x$G2M.Score, method = "spearman"),
  rho_RLcore_S = stats::cor(x$MT_core_UCell, x$S.Score, method = "spearman"),
  rho_RLcore_G2M = stats::cor(x$MT_core_UCell, x$G2M.Score, method = "spearman"),
  pct_cycling_High = 100 * mean(x$cycling[x$MT_group == "High"]),
  pct_cycling_Low = 100 * mean(x$cycling[x$MT_group == "Low"])
)))
utils::write.csv(per_an, file.path(OUT_SUB, "cellcycle_by_sample.csv"), row.names = FALSE)

hm_an_all <- stats::aggregate(md_an[, hcols], by = list(sample = md_an$sample, group = md_an$MT_group), FUN = mean)
hall_an_all <- do.call(rbind, lapply(hcols, function(h) {
  w <- stats::reshape(hm_an_all[, c("sample", "group", h)], idvar = "sample", timevar = "group", direction = "wide")
  d <- w[[paste0(h, ".High")]] - w[[paste0(h, ".Low")]]
  data.frame(pathway = sub("^HALLMARK_", "", sub("_UCell$", "", h)), mean_diff = mean(d, na.rm = TRUE),
             p = tryCatch(stats::wilcox.test(d, mu = 0)$p.value, error = function(e) NA_real_))
}))
hall_an_all$FDR <- stats::p.adjust(hall_an_all$p, "BH")

g1_an <- md_an[md_an$Phase == "G1", ]
keep_g1_an <- names(which(table(g1_an$sample) >= 10))
g1_an <- g1_an[g1_an$sample %in% keep_g1_an, ]
g1_an$G1_group <- stats::ave(g1_an$MT_core_UCell, g1_an$sample, FUN = function(v) ifelse(v >= stats::median(v), "High", "Low"))
hm_an_g1 <- stats::aggregate(g1_an[, hcols], by = list(sample = g1_an$sample, group = g1_an$G1_group), FUN = mean)
hall_an_g1 <- do.call(rbind, lapply(hcols, function(h) {
  w <- stats::reshape(hm_an_g1[, c("sample", "group", h)], idvar = "sample", timevar = "group", direction = "wide")
  d <- w[[paste0(h, ".High")]] - w[[paste0(h, ".Low")]]
  data.frame(pathway = sub("^HALLMARK_", "", sub("_UCell$", "", h)), mean_diff = mean(d, na.rm = TRUE),
             p = tryCatch(stats::wilcox.test(d, mu = 0)$p.value, error = function(e) NA_real_))
}))
hall_an_g1$FDR <- stats::p.adjust(hall_an_g1$p, "BH")
both_an <- merge(hall_an_all[, c("pathway", "mean_diff", "FDR")], hall_an_g1[, c("pathway", "mean_diff", "FDR")],
                 by = "pathway", suffixes = c("_all", "_G1"))
utils::write.csv(both_an, file.path(OUT_SUB, "Hallmark_all_vs_G1only.csv"), row.names = FALSE)
if (file.exists(file.path(OUT, "CellChat_count_weight_from_malignant.csv"))) {
  file.copy(file.path(OUT, "CellChat_count_weight_from_malignant.csv"),
            file.path(OUT_SUB, "CellChat_count_weight_from_malignant.csv"), overwrite = TRUE)
}

save_session("03b_sc_cellcycle")
log_msg("03b_sc_cellcycle 完成")
