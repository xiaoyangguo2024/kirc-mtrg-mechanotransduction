# =============================================================================
# 05_candidates.R —— 候选基因筛选：
#   DESeq2（TCGA-KIRC 肿瘤 vs 正常）+ WGCNA（与 MT_core / MTI 最相关共表达模块）+ 单细胞伪 bulk 差异基因
#   候选基因 = (MT_all ∪ WGCNA 关键模块基因) ∩ 单细胞高/低差异基因 ∩ TCGA 肿瘤 vs 正常差异基因
#   若三重交集少于 PARAM$min_candidates，则按预注册规则退回两层交集或放宽伪 bulk 阈值，并在方法学如实报告
#   随后在 TCGA-KIRC 中进行单因素 Cox 回归并检验跨验证队列基因覆盖情况
# =============================================================================
source("R/00_config.R")
need(c("DESeq2", "WGCNA", "survival", "ggplot2", "ggrepel", "VennDiagram"))
STEP <- "05_candidates"; OUT <- res_dir(STEP)
GS   <- readRDS(file.path(DIR_GS, "MTRG_genesets.rds"))
tcga <- readRDS(file.path(DIR_PRO, "TCGA_KIRC.rds"))
tum  <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))

# =============================================================================
# 1. DESeq2：TCGA-KIRC 肿瘤 vs 正常
# =============================================================================
f_deg <- file.path(OUT, "TCGA_DESeq2_Tumor_vs_Normal.csv")
if (!file.exists(f_deg)) {
  cl <- tcga$clin
  dds <- DESeq2::DESeqDataSetFromMatrix(
    round(tcga$counts[, cl$ID]),
    colData = data.frame(row.names = cl$ID, group = factor(cl$group, levels = c("Normal", "Tumor"))),
    design = ~ group
  )
  dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 10, ]
  dds <- DESeq2::DESeq(dds)
  deg <- as.data.frame(DESeq2::results(dds, contrast = c("group", "Tumor", "Normal")))
  deg$gene <- rownames(deg)
  utils::write.csv(deg, f_deg, row.names = FALSE)
}
deg <- utils::read.csv(f_deg)
deg_genes <- deg$gene[!is.na(deg$padj) & deg$padj < PARAM$deg_padj & abs(deg$log2FoldChange) > PARAM$deg_lfc]
log_msg("TCGA-KIRC 肿瘤 vs 正常差异基因：", length(deg_genes))

deg$class <- ifelse(deg$gene %in% deg_genes, ifelse(deg$log2FoldChange > 0, "Up", "Down"), "NS")
deg$is_MT <- deg$gene %in% GS$MT_all
lab <- utils::head(deg[deg$is_MT & deg$class != "NS", ][order(deg$padj[deg$is_MT & deg$class != "NS"]), ], 15)
p <- ggplot2::ggplot(deg[!is.na(deg$padj), ], ggplot2::aes(log2FoldChange, -log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = class), size = 0.5, alpha = 0.6) +
  ggplot2::geom_point(data = deg[deg$is_MT & deg$class != "NS", ], shape = 21, color = "black", size = 1.2) +
  ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene), size = 3) +
  ggplot2::scale_color_manual(values = c(Up = "#D7301F", Down = "#2C7FB8", NS = "grey80")) +
  ggplot2::theme_classic() + ggplot2::labs(title = "TCGA-KIRC tumor vs normal (circled: MTRGs)")
save_gg(p, file.path(OUT, "FigS_volcano_TCGA_DEG.pdf"), 6, 5)

# =============================================================================
# 2. WGCNA：TCGA-KIRC 肿瘤样本（MAD 最高的 5000 个基因）
# =============================================================================
f_wgcna <- file.path(DIR_PRO, "WGCNA_net.rds")
ex <- tcga$expr[, tum$ID]
datExpr <- t(ex[order(apply(ex, 1, stats::mad), decreasing = TRUE)[seq_len(PARAM$wgcna_top_mad)], ])
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
if (!file.exists(f_wgcna)) {
  WGCNA::allowWGCNAThreads(2)
  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)), networkType = "unsigned", verbose = 0)
  power <- sft$powerEstimate
  if (is.na(power)) { power <- 6; log_msg("【核对】无幂次达到无尺度拓扑 R² 0.85，按 WGCNA FAQ 取 6") }
  cor <- WGCNA::cor
  net <- WGCNA::blockwiseModules(datExpr, power = power, TOMType = "unsigned", minModuleSize = 30,
                                 mergeCutHeight = 0.25, numericLabels = FALSE, maxBlockSize = 6000,
                                 randomSeed = PARAM$seed, verbose = 0)
  cor <- stats::cor
  net$power <- power
  saveRDS(net, f_wgcna)
}
net <- readRDS(f_wgcna)
log_msg("WGCNA 软阈值 = ", net$power, "；模块数 = ", length(unique(net$colors)) - 1)

MEs <- WGCNA::orderMEs(net$MEs)
trait <- tum[match(rownames(datExpr), tum$ID), c("MT_core", "MTI", "Prolif", "EMT")]
mt_cor <- stats::cor(MEs, trait, use = "pairwise.complete.obs")
mt_p   <- WGCNA::corPvalueStudent(mt_cor, nrow(datExpr))

me_ok <- setdiff(colnames(MEs), "MEgrey")
best_me <- me_ok[which.max(abs(mt_cor[me_ok, "MT_core"]))]
mod <- sub("^ME", "", best_me)
mod_genes <- colnames(datExpr)[net$colors == mod]
MM  <- as.numeric(stats::cor(datExpr[, mod_genes], MEs[[best_me]], use = "p"))
GSv <- as.numeric(stats::cor(datExpr[, mod_genes], trait$MT_core, use = "p"))
key_mod_genes <- mod_genes[abs(MM) > PARAM$wgcna_mm & abs(GSv) > PARAM$wgcna_gs]
log_msg("WGCNA 关键模块：", mod, "（与 MT_core r = ", round(mt_cor[best_me, "MT_core"], 2), "，共 ",
        length(mod_genes), " 个基因；核心 hub 基因 ", length(key_mod_genes), " 个）")
utils::write.csv(data.frame(gene = mod_genes, MM = MM, GS = GSv, key = mod_genes %in% key_mod_genes),
                 file.path(OUT, paste0("WGCNA_module_", mod, "_genes.csv")), row.names = FALSE)

# =============================================================================
# 3. 候选基因交集
# =============================================================================
f_scdeg <- file.path(DIR_RES, "03_scRNA", "sc_DEG_genes.txt")
sc_deg <- if (file.exists(f_scdeg)) readLines(f_scdeg) else character(0)
pool <- union(GS$MT_all, key_mod_genes)
cand <- Reduce(intersect, list(pool, sc_deg, deg_genes))
rule <- "(MT_all ∪ WGCNA key-module genes) ∩ scRNA pseudobulk DEGs ∩ TCGA tumor-vs-normal DEGs"
if (length(cand) < PARAM$min_candidates) {
  # 尝试名义 P < 0.05 且 |logFC| > 0.3 的单细胞差异基因
  pbf <- file.path(DIR_RES, "03_scRNA", "pseudobulk_DE_High_vs_Low.csv")
  if (file.exists(pbf)) {
    pbd <- utils::read.csv(pbf)
    sc_deg_relaxed <- pbd$gene[pbd$PValue < 0.05 & abs(pbd$logFC) > 0.3]
    cand2 <- Reduce(intersect, list(pool, sc_deg_relaxed, deg_genes))
    if (length(cand2) >= PARAM$min_candidates) {
      cand <- cand2
      sc_deg <- sc_deg_relaxed
      rule <- "(MT_all ∪ WGCNA key-module genes) ∩ scRNA pseudobulk DEGs (P < 0.05, |log2FC| > 0.3) ∩ TCGA DEGs"
    }
  }
}
if (length(cand) < PARAM$min_candidates) {
  cand <- intersect(pool, deg_genes)
  rule <- "(MT_all ∪ WGCNA key-module genes) ∩ TCGA tumor-vs-normal DEGs"
}
log_msg("候选基因数：", length(cand), "；筛选规则：", rule)

# =============================================================================
# 4. 单因素 Cox（TCGA-KIRC，OS）与跨验证队列覆盖度检查
# =============================================================================
d <- tum[!is.na(tum$OS.time) & tum$OS.time > 0 & !is.na(tum$OS), ]
ex_c <- tcga$expr[cand, d$ID, drop = FALSE]
uni <- do.call(rbind, lapply(cand, function(g) {
  z <- as.numeric(scale(ex_c[g, ]))
  if (all(is.na(z))) return(NULL)
  r <- cox_table(survival::coxph(survival::Surv(d$OS.time, d$OS) ~ z))
  data.frame(gene = g, HR = r$HR, lower = r$lower, upper = r$upper, p = r$p)
}))
uni$FDR <- stats::p.adjust(uni$p, "BH")
uni <- uni[order(uni$p), ]
sig <- uni[uni$p < PARAM$unicox_p, ]
log_msg("单因素 Cox P < ", PARAM$unicox_p, " 基因数：", nrow(sig))

# 检查所有具有随访数据的外部队列中该基因是否存在
ALL_RDS <- sub("\\.rds$", "", list.files(DIR_PRO, pattern = "\\.rds$"))
val_cand <- intersect(c("GSE29609", "GSE22541", "GSE167573", "IMvigor210", "GSE176307"), ALL_RDS)
# 仅把真正有 >=15 例随访且 >=5 例事件的外部预后队列作为建模必须覆盖的验证队列
val <- Filter(function(n) {
  cl <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$clin
  has_os <- !is.na(cl$OS.time) & !is.na(cl$OS) & cl$OS.time > 0
  sum(has_os) >= 15 && sum(cl$OS[has_os] == 1) >= 5
}, val_cand)
log_msg("具有有效生存数据的外部验证队列：", paste(val, collapse = ", "))

avail <- do.call(cbind, lapply(stats::setNames(val, val), function(n)
  cand %in% rownames(readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$expr)))
avail <- data.frame(gene = cand, avail, all_cohorts = apply(avail, 1, all), check.names = FALSE)
# 若某芯片平台（如 GPL1708）基因数较少导致交集过小，允许在至少 (length(val)-1) 个验证队列中存在
if (sum(avail$all_cohorts & avail$gene %in% sig$gene) < 8) {
  avail$all_cohorts <- rowSums(avail[, val, drop = FALSE]) >= max(1, length(val) - 1)
  log_msg("提示：放宽跨平台覆盖条件为至少在 ", max(1, length(val) - 1), " 个验证队列存在，合格候选基因：",
          sum(avail$all_cohorts))
}

S3 <- merge(data.frame(gene = cand), deg[, c("gene", "log2FoldChange", "padj")], by = "gene", all.x = TRUE)
pbf <- file.path(DIR_RES, "03_scRNA", "pseudobulk_DE_High_vs_Low.csv")
if (file.exists(pbf)) {
  pbd <- utils::read.csv(pbf)
  S3 <- merge(S3, stats::setNames(pbd[, c("gene", "logFC", "FDR")], c("gene", "sc_logFC_High_vs_Low", "sc_FDR")),
              by = "gene", all.x = TRUE)
}
S3$in_MT_core   <- S3$gene %in% GS$MT_core
S3$in_MT_all    <- S3$gene %in% GS$MT_all
S3$in_WGCNA_key <- S3$gene %in% key_mod_genes
S3 <- merge(S3, uni, by = "gene", all.x = TRUE)
S3 <- merge(S3, avail[, c("gene", "all_cohorts")], by = "gene", all.x = TRUE)
utils::write.csv(S3[order(S3$p), ], file.path(OUT, "Table_S3_candidates_unicox.csv"), row.names = FALSE)

saveRDS(list(
  candidates = cand, rule = rule, unicox = uni, unicox_sig = sig$gene,
  key_module = mod, key_module_genes = key_mod_genes, availability = avail, valid_cohorts = val
), file.path(DIR_PRO, "candidates.rds"))
save_session(STEP)
log_msg("05_candidates 完成")
