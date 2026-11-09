# =============================================================================
# 09_hub_gene_spatial.R —— MTRS 核心基因解析：
#   1. 精确线性或核 SHAP 分解（shapviz）识别对 MTRS 贡献最大的 Top 3 核心基因
#   2. TCGA-KIRC 肿瘤 vs 正常（全部 + 72 对配对）及跨队列预后验证
#   3. 单细胞（GSE159115）定位：各细胞亚群表达及高/低膜张力恶性上皮细胞样本内配对比较
#   4. 空间转录组（GSE210041，ccRCC 10x Visium）：切片内 spot 水平 MT_core 与血管生成/EMT/增殖/免疫
#      及核心基因的共定位，并校正测序深度 log10(UMI) 与增殖评分（秩偏相关）
#   5. 核心基因增殖重叠核查与跨队列覆盖度汇总（含 09b / 09c 功能）
# =============================================================================
source("R/00_config.R")
need(c("shapviz", "survival", "survminer", "ggplot2", "ggpubr", "Seurat", "UCell", "msigdbr", "Matrix"))
STEP <- "09_hub_gene_spatial"; OUT <- res_dir(STEP)
N_SHOW <- 3

GS    <- readRDS(file.path(DIR_GS, "MTRG_genesets.rds"))
model <- readRDS(file.path(DIR_PRO, "RLRS_model.rds"))
tcga  <- readRDS(file.path(DIR_PRO, "TCGA_KIRC.rds"))
tum   <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
genes <- model$genes
log_msg("MTRS 模型基因（", length(genes), " 个）：", paste(genes, collapse = ", "))

# =============================================================================
# 1. 模型基因对 MTRS 的贡献（SHAP 分解）
# =============================================================================
d <- tum[!is.na(tum$OS.time) & tum$OS.time > 0 & !is.na(tum$OS), c("ID", "OS.time", "OS")]
mi <- build_mime_input(tcga$expr, tcga$clin[tcga$clin$group == "Tumor", ], model$candidate_genes, fill_missing = TRUE)
rl_all <- readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
mi <- merge(rl_all[rl_all$cohort == "TCGA_KIRC", c("ID", "RLRS")], mi, by = "ID")
X <- as.matrix(mi[, safe_names(genes), drop = FALSE]); colnames(X) <- genes
lf <- stats::lm(mi$RLRS ~ X)
r2 <- suppressWarnings(summary(lf)$r.squared)
w  <- stats::coef(lf)[-1]; w[is.na(w)] <- 0; names(w) <- genes
log_msg("MTRS 线性还原 R² = ", round(r2, 5))
utils::write.csv(data.frame(term = c("(Intercept)", genes), coefficient = c(unname(stats::coef(lf)[1]), unname(w))),
                 file.path(OUT, "RLRS_coefficients.csv"), row.names = FALSE)
utils::write.csv(data.frame(term = c("(Intercept)", genes), coefficient = c(unname(stats::coef(lf)[1]), unname(w))),
                 file.path(OUT, "MTRS_coefficients.csv"), row.names = FALSE)

S <- sweep(sweep(X, 2, colMeans(X)), 2, w, `*`)
sv <- shapviz::shapviz(S, X = as.data.frame(X), baseline = mean(mi$RLRS))
imp <- sort(colMeans(abs(S)), decreasing = TRUE)
rank_tab <- data.frame(
  gene = names(imp),
  coefficient = unname(w[names(imp)]),
  mean_abs_SHAP = unname(imp),
  share_pct = unname(100 * imp / sum(imp))
)
print(rank_tab, row.names = FALSE)
utils::write.csv(rank_tab, file.path(OUT, "SHAP_gene_ranking.csv"), row.names = FALSE)
hubs <- utils::head(names(imp), min(N_SHOW, length(imp)))
hub  <- hubs[1]
log_msg("核心基因 Top ", length(hubs), "：", paste(hubs, collapse = ", "))

save_print(c(list(shapviz::sv_importance(sv, kind = "beeswarm"), shapviz::sv_importance(sv, kind = "bar")),
             lapply(hubs, function(v) shapviz::sv_dependence(sv, v = v))),
           file.path(OUT, "Fig9_SHAP.pdf"), 6, 5)

# =============================================================================
# 2. 核心基因在 TCGA-KIRC 肿瘤 vs 正常（含配对）与跨队列单基因 Cox
# =============================================================================
cl <- tcga$clin
pl <- list()
for (g in hubs) {
  dd <- data.frame(cl[, c("ID", "patient", "group")], expr = as.numeric(tcga$expr[g, cl$ID]))
  dd$group <- factor(dd$group, levels = c("Normal", "Tumor"))
  pat <- intersect(dd$patient[dd$group == "Tumor"], dd$patient[dd$group == "Normal"])
  pl[[paste(g, "TvN")]] <- ggpubr::ggboxplot(dd, "group", "expr", color = "group", add = "jitter",
                                             palette = c("#2C7FB8", "#D7301F")) +
    ggpubr::stat_compare_means(method = "wilcox.test") + ggplot2::labs(title = g, y = "log2(TPM + 1)")
  pl[[paste(g, "paired")]] <- ggpubr::ggpaired(dd[dd$patient %in% pat, ], x = "group", y = "expr", id = "patient",
                                               color = "group", line.color = "grey70",
                                               palette = c("#2C7FB8", "#D7301F")) +
    ggpubr::stat_compare_means(paired = TRUE) + ggplot2::labs(title = paste(g, "paired"), y = "log2(TPM + 1)")
  s <- data.frame(d, expr = as.numeric(tcga$expr[g, d$ID]))
  s$group <- median_group(s$expr)
  pl[[paste(g, "KM")]] <- km_plot(s, title = paste("TCGA-KIRC", g))
}
save_print(pl, file.path(OUT, "Fig9_hub_TCGA_external.pdf"), 6, 6)

# =============================================================================
# 3. 单细胞（GSE159115）：核心基因在各细胞类型及高/低膜张力恶性细胞中的表达
# =============================================================================
f_sc <- file.path(DIR_PRO, "scRNA_integrated.rds")
f_md <- file.path(DIR_RES, "03_scRNA", "cell_metadata.csv")
if (file.exists(f_sc) && file.exists(f_md)) {
  sc <- readRDS(f_sc)
  hs <- intersect(hubs, rownames(sc))
  emb <- Embeddings(sc, "umap")
  hx <- as.matrix(LayerData(sc, assay = "RNA", layer = "data")[hs, , drop = FALSE])
  hd <- data.frame(cell = colnames(sc), UMAP1 = emb[colnames(sc), 1], UMAP2 = emb[colnames(sc), 2], t(hx), check.names = FALSE)
  rm(sc, emb, hx); invisible(gc())
  md <- utils::read.csv(f_md, row.names = 1, check.names = FALSE)
  hd <- cbind(hd, md[hd$cell, c("sample", "cell_class", "RL_group")])
  hd$cell_class <- factor(hd$cell_class)

  dp <- do.call(rbind, lapply(hs, function(g) do.call(rbind, lapply(split(hd[[g]], hd$cell_class), function(x)
    data.frame(gene = g, avg = mean(x), pct = mean(x > 0))))))
  dp$cell_class <- rep(levels(hd$cell_class), length(hs))
  utils::write.csv(dp, file.path(OUT, "hub_scRNA_by_cell_class.csv"), row.names = FALSE)

  hm <- hd[!is.na(hd$RL_group) & hd$cell_class == "Malignant", ]
  agg <- stats::aggregate(hm[, hs, drop = FALSE], by = list(sample = hm$sample, group = hm$RL_group), FUN = mean)
  agg$group <- factor(agg$group, levels = c("Low", "High"))
  agg <- agg[order(agg$sample, agg$group), ]
  utils::write.csv(agg, file.path(OUT, "hub_scRNA_malignant_High_vs_Low_per_sample.csv"), row.names = FALSE)
  rm(hd, hm); invisible(gc())
}

# =============================================================================
# 4. 空间转录组（GSE210041，10x Visium ccRCC 切片）
# =============================================================================
ex_dir <- file.path(DIR_RAW, "GEO", "GSE210041", "extracted")
f_tar  <- file.path(DIR_RAW, "GEO", "GSE210041", "GSE210041_RAW.tar")
if ((!dir.exists(ex_dir) || length(list.files(ex_dir)) == 0) && file.exists(f_tar)) {
  dir.create(ex_dir, recursive = TRUE, showWarnings = FALSE)
  utils::untar(f_tar, exdir = ex_dir)
}

# 自动识别 GSE210041 解压目录中的样本前缀（GSMxxxx_...）并构建标准 Visium 或表达+坐标对象
raw_files <- list.files(ex_dir, full.names = TRUE)
gsms <- unique(regmatches(basename(raw_files), regexpr("^GSM[0-9]+", basename(raw_files))))
log_msg("发现 GSE210041 空间样本：", paste(gsms, collapse = ", "))

H <- get_msig("H")
SP_SETS <- c(
  GS[c("MT_core", "MT_all", "MT_cortex", "MT_BAR", "MT_channel", "MT_ECM", "MT_stiff", "MT_buffer")],
  list(
    E2F          = H[["HALLMARK_E2F_TARGETS"]],
    G2M          = H[["HALLMARK_G2M_CHECKPOINT"]],
    EMT          = H[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
    Hypoxia      = H[["HALLMARK_HYPOXIA"]],
    Angiogenesis = H[["HALLMARK_ANGIOGENESIS"]],
    IFN_gamma    = H[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
    GEP18        = GEP18,
    Angio_sig    = ANGIO_SIG
  )
)

prho <- function(x, y, z) {
  Z <- as.data.frame(lapply(as.data.frame(z), rank))
  rx <- stats::resid(stats::lm(rank(x) ~ ., data = Z))
  ry <- stats::resid(stats::lm(rank(y) ~ ., data = Z))
  stats::cor(rx, ry)
}
sp_cor <- function(x, y) suppressWarnings(stats::cor(x, y, method = "spearman", use = "complete.obs"))

sp_prog <- list(); sp_rows <- list(); dep <- list(); hub_part <- list(); slim <- list()

for (idx in seq_along(gsms)) {
  gsm <- gsms[idx]
  smp <- paste0("ccRCC_Visium_", idx)
  f_obj <- file.path(DIR_PRO, "spatial", paste0(smp, ".rds"))
  dir.create(dirname(f_obj), showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(f_obj)) {
    gf <- raw_files[ grep(paste0("^", gsm), basename(raw_files)) ]
    log_msg(smp, " (", gsm, ") 包含文件：", paste(basename(gf), collapse = ", "))
    sdir <- file.path(DIR_RAW, "GEO", "GSE210041", smp)
    dir.create(file.path(sdir, "spatial"), recursive = TRUE, showWarnings = FALSE)

    # 判断是 h5 还是 mtx（matrix.mtx.gz + barcodes.tsv.gz + features.tsv.gz）还是 tar.gz
    f_h5  <- gf[grep("\\.h5(\\.gz)?$", gf)]
    f_tar_sub <- gf[grep("\\.tar\\.gz$", gf)]
    if (length(f_tar_sub)) {
      utils::untar(f_tar_sub[1], exdir = sdir)
    }
    f_mtx <- gf[grep("matrix\\.mtx", gf)]
    f_bar <- gf[grep("barcodes\\.tsv", gf)]
    f_fea <- gf[grep("features\\.tsv|genes\\.tsv", gf)]

    cnt_mat <- NULL
    if (length(f_h5)) {
      h5_use <- f_h5[1]
      if (grepl("\\.gz$", h5_use)) {
        h5_dest <- file.path(sdir, "filtered_feature_bc_matrix.h5")
        if (!file.exists(h5_dest)) gunzip_to(h5_use, h5_dest)
        h5_use <- h5_dest
      }
      cnt_mat <- Seurat::Read10X_h5(h5_use)
      if (is.list(cnt_mat)) cnt_mat <- cnt_mat[[1]]
    } else if (length(f_mtx) && length(f_bar) && length(f_fea)) {
      file.copy(f_mtx[1], file.path(sdir, "matrix.mtx.gz"), overwrite = TRUE)
      file.copy(f_bar[1], file.path(sdir, "barcodes.tsv.gz"), overwrite = TRUE)
      file.copy(f_fea[1], file.path(sdir, "features.tsv.gz"), overwrite = TRUE)
      cnt_mat <- Seurat::Read10X(sdir)
      if (is.list(cnt_mat)) cnt_mat <- cnt_mat[[1]]
    } else {
      # 检查解压后的子目录
      sub_h5 <- list.files(sdir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
      sub_mtx <- list.files(sdir, pattern = "matrix\\.mtx", recursive = TRUE, full.names = TRUE)
      if (length(sub_h5)) {
        cnt_mat <- Seurat::Read10X_h5(sub_h5[1])
      } else if (length(sub_mtx)) {
        cnt_mat <- Seurat::Read10X(dirname(sub_mtx[1]))
      }
    }

    if (is.null(cnt_mat)) { log_msg("警告：", smp, " 未找到表达矩阵，跳过"); next }

    # 读取空间坐标（tissue_positions_list.csv.gz 或 tissue_positions.csv.gz）
    f_pos <- c(gf[grep("tissue_positions", gf)],
               list.files(sdir, pattern = "tissue_positions", recursive = TRUE, full.names = TRUE))
    coords_df <- NULL
    if (length(f_pos)) {
      pos_tab <- data.table::fread(f_pos[1], data.table = FALSE)
      # 第 1 列为 barcode，最后两列或第 5-6 列为像素坐标
      bc_col <- 1
      if (ncol(pos_tab) >= 6) {
        coords_df <- data.frame(
          barcode = as.character(pos_tab[[1]]),
          in_tissue = as.integer(pos_tab[[2]]),
          row = as.numeric(pos_tab[[3]]),
          col = as.numeric(pos_tab[[4]]),
          imagerow = as.numeric(pos_tab[[5]]),
          imagecol = as.numeric(pos_tab[[6]]),
          stringsAsFactors = FALSE
        )
      }
    }

    ob <- CreateSeuratObject(cnt_mat, project = smp, assay = "Spatial", min.cells = 3, min.features = 100)
    ob$nCount_RNA <- ob$nCount_Spatial
    ob$nFeature_RNA <- ob$nFeature_Spatial
    if (!is.null(coords_df)) {
      m_bc <- match(colnames(ob), coords_df$barcode)
      if (sum(!is.na(m_bc)) < ncol(ob) * 0.5) {
        m_bc <- match(sub("-1$", "", colnames(ob)), sub("-1$", "", coords_df$barcode))
      }
      ob$x <- ifelse(!is.na(m_bc), coords_df$row[m_bc], seq_len(ncol(ob)))
      ob$y <- ifelse(!is.na(m_bc), coords_df$col[m_bc], seq_len(ncol(ob)))
    } else {
      # 若文件名不含独立坐标表，按规则网格生成展示坐标
      n_sp <- ncol(ob); w_grid <- ceiling(sqrt(n_sp))
      ob$x <- ((seq_len(n_sp) - 1) %/% w_grid)
      ob$y <- ((seq_len(n_sp) - 1) %% w_grid)
    }

    ob <- subset(ob, subset = nCount_Spatial > 200)
    ob <- NormalizeData(ob)
    ob <- UCell::AddModuleScore_UCell(ob, features = SP_SETS, maxRank = 2500, name = "_UCell", BPPARAM = BiocParallel::SerialParam())
    ob$MTI    <- as.numeric(scale(ob$MT_stiff_UCell) - scale(ob$MT_buffer_UCell))
    ob$Prolif <- as.numeric(scale(ob$E2F_UCell) + scale(ob$G2M_UCell)) / 2
    ob$RL_core_UCell <- ob$MT_core_UCell
    ob$RHI           <- ob$MTI
    saveRDS(ob, f_obj)
  }

  ob <- readRDS(f_obj)
  m  <- ob@meta.data
  lu <- log10(m$nCount_Spatial)

  sp_prog[[smp]] <- data.frame(
    sample = smp, gsm = gsm, n_spots = ncol(ob),
    rho_RLcore_Prolif = sp_cor(m$MT_core_UCell, m$Prolif),
    rho_RLcore_EMT    = sp_cor(m$MT_core_UCell, m$EMT_UCell),
    rho_RLcore_Angio  = sp_cor(m$MT_core_UCell, m$Angiogenesis_UCell),
    rho_RLcore_IFNg   = sp_cor(m$MT_core_UCell, m$IFN_gamma_UCell),
    rho_RLcore_GEP18  = sp_cor(m$MT_core_UCell, m$GEP18_UCell),
    prho_RLcore_IFNg_given_Prolif = prho(m$MT_core_UCell, m$IFN_gamma_UCell, m$Prolif)
  )

  dep[[smp]] <- data.frame(
    sample = smp, n_spots = nrow(m), median_UMI = stats::median(m$nCount_Spatial),
    rho_UMI_RLcore = sp_cor(lu, m$MT_core_UCell),
    rho_UMI_Prolif = sp_cor(lu, m$Prolif),
    rho_UMI_IFNg   = sp_cor(lu, m$IFN_gamma_UCell),
    rho_UMI_GEP18  = sp_cor(lu, m$GEP18_UCell),
    rho_RLcore_Prolif = sp_cor(m$MT_core_UCell, m$Prolif),
    prho_RLcore_Prolif_given_UMI = prho(m$MT_core_UCell, m$Prolif, lu),
    rho_RLcore_IFNg = sp_cor(m$MT_core_UCell, m$IFN_gamma_UCell),
    prho_RLcore_IFNg_given_UMI = prho(m$MT_core_UCell, m$IFN_gamma_UCell, lu),
    prho_RLcore_IFNg_given_Prolif_UMI = prho(m$MT_core_UCell, m$IFN_gamma_UCell, cbind(m$Prolif, lu)),
    rho_Prolif_IFNg = sp_cor(m$Prolif, m$IFN_gamma_UCell),
    prho_Prolif_IFNg_given_UMI = prho(m$Prolif, m$IFN_gamma_UCell, lu),
    rho_RLcore_GEP18 = sp_cor(m$MT_core_UCell, m$GEP18_UCell),
    prho_RLcore_GEP18_given_Prolif_UMI = prho(m$MT_core_UCell, m$GEP18_UCell, cbind(m$Prolif, lu))
  )

  hs <- intersect(hubs, rownames(ob))
  ex <- as.matrix(SeuratObject::LayerData(ob, assay = "Spatial", layer = "data")[hs, rownames(m), drop = FALSE])
  for (g in hs) {
    x <- as.numeric(ex[g, ])
    ct <- suppressWarnings(stats::cor.test(m$MT_core_UCell, x, method = "spearman", exact = FALSE))
    sp_rows[[paste(smp, g)]] <- data.frame(
      sample = smp, gene = g, n_spots = ncol(ob), pct_detected = mean(x > 0),
      rho = unname(ct$estimate), p = ct$p.value
    )
    hub_part[[paste(smp, g)]] <- data.frame(
      sample = smp, gene = g, pct_detected = mean(x > 0),
      rho = sp_cor(m$MT_core_UCell, x),
      prho_given_UMI = prho(m$MT_core_UCell, x, lu),
      prho_given_Prolif_UMI = prho(m$MT_core_UCell, x, cbind(m$Prolif, lu))
    )
  }

  slim[[smp]] <- data.frame(
    sample = smp, spot = rownames(m), x = m$x, y = m$y,
    UMI = m$nCount_Spatial, RL_core = m$MT_core_UCell, MT_core = m$MT_core_UCell,
    Prolif = m$Prolif, EMT = m$EMT_UCell, Angiogenesis = m$Angiogenesis_UCell,
    IFN_gamma = m$IFN_gamma_UCell, t(ex), check.names = FALSE, row.names = NULL
  )
  rm(ob, m, ex); invisible(gc())
}

if (length(sp_prog)) {
  utils::write.csv(do.call(rbind, sp_prog), file.path(OUT, "Spatial_RLcore_vs_programs_by_slice.csv"), row.names = FALSE)
}
if (length(sp_rows)) {
  sp_tab <- do.call(rbind, sp_rows); sp_tab$FDR <- stats::p.adjust(sp_tab$p, "BH")
  utils::write.csv(sp_tab, file.path(OUT, "Spatial_RLcore_vs_hub_spearman.csv"), row.names = FALSE)
}
if (length(dep)) {
  utils::write.csv(do.call(rbind, dep), file.path(OUT, "Spatial_depth_check_by_slice.csv"), row.names = FALSE)
}
if (length(hub_part)) {
  utils::write.csv(do.call(rbind, hub_part), file.path(OUT, "Spatial_hub_partial_by_slice.csv"), row.names = FALSE)
}
if (length(slim)) {
  saveRDS(do.call(rbind, slim), file.path(DIR_PRO, "spatial_spots_slim.rds"))
}

# =============================================================================
# 5. MTRS 模型基因增殖重叠与跨队列覆盖度（09c 功能整合）
# =============================================================================
bp <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
               error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"))
go_set <- function(id) unique(bp$gene_symbol[bp$gs_exact_source == id])
cc <- go_set("GO:0007049"); rep_dna <- go_set("GO:0006260")

t1 <- tum[tum$group == "Tumor", ]
ex_tcga <- tcga$expr[, t1$ID]
rho_p <- vapply(genes, function(g) if (g %in% rownames(ex_tcga))
  suppressWarnings(stats::cor(as.numeric(ex_tcga[g, ]), t1$Prolif, method = "spearman")) else NA_real_, numeric(1))

res_g <- data.frame(
  gene = genes,
  in_HALLMARK_E2F = genes %in% H[["HALLMARK_E2F_TARGETS"]],
  in_HALLMARK_G2M = genes %in% H[["HALLMARK_G2M_CHECKPOINT"]],
  in_HALLMARK_EMT = genes %in% H[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
  in_GO_cell_cycle = genes %in% cc,
  in_GO_DNA_replication = genes %in% rep_dna,
  rho_Prolif_TCGA = round(rho_p, 3)
)
res_g <- res_g[order(-res_g$rho_Prolif_TCGA), ]
utils::write.csv(res_g, file.path(OUT, "RLRS_genes_proliferation.csv"), row.names = FALSE)

all_coh <- unique(c(model$train, model$valid, model$ther))
cov_g <- do.call(rbind, lapply(all_coh, function(n) {
  f <- file.path(DIR_PRO, paste0(n, ".rds"))
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); rn <- rownames(x$expr); rm(x)
  data.frame(cohort = n, n_present = sum(genes %in% rn), missing = paste(setdiff(genes, rn), collapse = ", "))
}))
utils::write.csv(cov_g, file.path(OUT, "RLRS_gene_coverage.csv"), row.names = FALSE)

save_session(STEP)
log_msg("09_hub_gene_spatial 完成")
