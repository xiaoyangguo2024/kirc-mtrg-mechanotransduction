# =============================================================================
# 03_scRNA.R —— ccRCC 单细胞转录组（GSE159115）：
#   质控 → Harmony 整合 → 细胞注释 → copyKAT 非整倍体分析（含正常肾上皮阴性对照）
#   → UCell 质膜张力（MT_core / MT_all / MTI）评分 → 样本内中位数分组
#   → 伪 bulk 差异分析（edgeR 配对模型）→ Hallmark 富集 → CellChat 细胞通讯
# =============================================================================
source("R/00_config.R")
need(c("Seurat", "harmony", "UCell", "edgeR", "data.table", "Matrix", "ggplot2", "patchwork", "ggpubr", "ggrepel", "hdf5r"))
STEP <- "03_scRNA"; OUT <- res_dir(STEP)
GS <- readRDS(file.path(DIR_GS, "MTRG_genesets.rds"))

MALIGNANT_RULE <- "all_tumor_epi"  # 肿瘤样本中全部肾小管上皮细胞定义为肿瘤上皮（Malignant），同时输出 copyKAT 比例与敏感性对比

# =============================================================================
# 1. 构建 GSE159115 样本表（排除非 ccRCC 病理类型，仅保留 ccRCC 肿瘤与正常肾组织）
# =============================================================================
sc_dir <- file.path(DIR_RAW, "GEO", "GSE159115")
pd_file <- file.path(sc_dir, "GSE159115_pData.rds")
if (!file.exists(pd_file)) {
  es_sc <- GEOquery::getGEO("GSE159115", GSEMatrix = TRUE, getGPL = FALSE, destdir = file.path(DIR_RAW, "GEO"))[[1]]
  pd_sc <- Biobase::pData(es_sc)
  saveRDS(pd_sc, pd_file)
} else {
  pd_sc <- readRDS(pd_file)
}

supp_col <- grep("^supplementary_file", colnames(pd_sc), value = TRUE)[1]
all_info <- apply(pd_sc, 1, paste, collapse = " ")
is_non_ccrcc <- grepl("papillary|chromophobe|pRCC|chRCC", all_info, ignore.case = TRUE) &
                !grepl("clear cell|ccRCC|normal|benign", pd_sc$title, ignore.case = TRUE)

SC_SAMPLES <- data.frame(
  dataset = "GSE159115",
  gsm     = pd_sc$geo_accession,
  title   = as.character(pd_sc$title),
  file    = basename(as.character(pd_sc[[supp_col]])),
  tissue  = ifelse(grepl("Benign Adjacent|normal|benign|NAT", all_info, ignore.case = TRUE), "Normal", "Tumor"),
  stringsAsFactors = FALSE
)
SC_SAMPLES <- SC_SAMPLES[!is_non_ccrcc & grepl("\\.h5$", SC_SAMPLES$file), ]
# 生成不含下划线的简洁样本编号（Seurat AggregateExpression 以 "_" 分隔分组）
SC_SAMPLES$sample <- paste0(ifelse(SC_SAMPLES$tissue == "Tumor", "ccRCC", "Norm"),
                            ave(seq_len(nrow(SC_SAMPLES)), SC_SAMPLES$tissue, FUN = seq_along))
print(SC_SAMPLES)
utils::write.csv(SC_SAMPLES, file.path(OUT, "scRNA_samples.csv"), row.names = FALSE)

# =============================================================================
# 2. 读取 10x HDF5、质控、Harmony 整合、聚类
# =============================================================================
f_base <- file.path(DIR_PRO, "scRNA_base.rds")
if (!file.exists(f_base)) {
  objs <- lapply(seq_len(nrow(SC_SAMPLES)), function(i) {
    s <- SC_SAMPLES[i, ]
    f_h5 <- file.path(sc_dir, s$file)
    if (!file.exists(f_h5)) f_h5 <- gsm_suppl(s$gsm, s$file, "GSE159115")
    m <- Seurat::Read10X_h5(f_h5, use.names = TRUE, unique.features = TRUE)
    if (is.list(m)) m <- m[[1]]
    if (ncol(m) > PARAM$sc_max_cells_per_sample) {
      set.seed(PARAM$seed + i)
      m <- m[, sort(sample(ncol(m), PARAM$sc_max_cells_per_sample))]
    }
    colnames(m) <- paste0(s$sample, "_", make.unique(sub("[.-][0-9]+$", "", colnames(m))))
    o <- CreateSeuratObject(m, project = s$sample, min.cells = 3, min.features = PARAM$sc_min_features)
    o$sample  <- s$sample
    o$gsm     <- s$gsm
    o$dataset <- s$dataset
    o$tissue  <- s$tissue
    log_msg(s$sample, " (", s$gsm, ")：", ncol(o), " 个细胞")
    o
  })
  sc <- merge(objs[[1]], objs[-1]); rm(objs); invisible(gc())
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  save_gg(VlnPlot(sc, c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", pt.size = 0, ncol = 1),
          file.path(OUT, "QC_before_filter.pdf"), 12, 10)
  sc <- subset(sc, subset = nFeature_RNA > PARAM$sc_min_features &
                            nFeature_RNA < PARAM$sc_max_features &
                            percent.mt < PARAM$sc_max_mt)
  sc <- NormalizeData(sc)
  sc <- FindVariableFeatures(sc, nfeatures = 2000)
  sc <- ScaleData(sc)
  sc <- RunPCA(sc, npcs = PARAM$sc_npcs, verbose = FALSE)
  sc <- harmony::RunHarmony(sc, group.by.vars = "sample")
  sc <- RunUMAP(sc, reduction = "harmony", dims = 1:PARAM$sc_npcs)
  sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:PARAM$sc_npcs)
  sc <- FindClusters(sc, resolution = PARAM$sc_resolution)
  sc <- JoinLayers(sc)
  saveRDS(sc, f_base)
}
sc <- readRDS(f_base)
if ("ccRCC13" %in% sc$sample) {
  old_s <- paste0("ccRCC", seq_len(nrow(SC_SAMPLES)))
  new_s <- SC_SAMPLES$sample
  old_cells <- colnames(sc)
  new_cells <- old_cells
  for (i in seq_len(nrow(SC_SAMPLES))) {
    idx <- sc$gsm == SC_SAMPLES$gsm[i]
    new_cells[idx] <- sub(paste0("^", old_s[i], "_"), paste0(new_s[i], "_"), old_cells[idx])
  }
  sc <- RenameCells(sc, new.names = new_cells)
  sc$sample <- SC_SAMPLES$sample[match(sc$gsm, SC_SAMPLES$gsm)]
  sc$tissue <- SC_SAMPLES$tissue[match(sc$gsm, SC_SAMPLES$gsm)]
  saveRDS(sc, f_base)
  ck_dir0 <- file.path(OUT, "copykat")
  unlink(list.files(ck_dir0, pattern = "(gene_by_cell|bin_by_cell|input_counts)", full.names = TRUE))
  tmp_d <- file.path(ck_dir0, "tmp_mig"); dir.create(tmp_d, showWarnings = FALSE)
  for (i in seq_len(nrow(SC_SAMPLES))) {
    f_old <- file.path(ck_dir0, paste0(old_s[i], "_copykat_prediction.txt"))
    if (file.exists(f_old)) {
      ln <- readLines(f_old)
      ln <- gsub(paste0(old_s[i], "_"), paste0(new_s[i], "_"), ln, fixed = TRUE)
      writeLines(ln, file.path(tmp_d, paste0(new_s[i], "_copykat_prediction.txt")))
    }
    j_old <- file.path(ck_dir0, paste0(old_s[i], "_copykat_bin_by_cell_heatmap.jpeg"))
    if (file.exists(j_old)) {
      file.copy(j_old, file.path(tmp_d, paste0(new_s[i], "_copykat_bin_by_cell_heatmap.jpeg")), overwrite = TRUE)
    }
  }
  unlink(list.files(ck_dir0, pattern = "^ccRCC", full.names = TRUE))
  for (ff in list.files(tmp_d, full.names = TRUE)) file.rename(ff, file.path(ck_dir0, basename(ff)))
  unlink(tmp_d, recursive = TRUE)
}
if ("scale.data" %in% SeuratObject::Layers(sc[["RNA"]])) sc[["RNA"]]$scale.data <- NULL
print(table(sc$sample, sc$tissue))

# =============================================================================
# 3. 细胞类型注释（结合 ccRCC / 肾小管经典 marker 打分）
# =============================================================================
MARKERS$Epithelial <- unique(c(MARKERS$Epithelial, "LRP2", "CUBN", "UMOD", "SLC12A1", "AQP2", "ATP6V1B1", "KRT19"))
mk <- lapply(MARKERS, function(g) intersect(g, rownames(sc)))
mk <- mk[lengths(mk) >= 2]
sc <- AddModuleScore(sc, features = mk, name = "mk")
mk_cols <- paste0("mk", seq_along(mk))
cl_mean <- stats::aggregate(sc@meta.data[, mk_cols], by = list(cluster = sc$seurat_clusters), FUN = mean)
ann <- data.frame(cluster = cl_mean$cluster, celltype_auto = names(mk)[apply(as.matrix(cl_mean[, -1]), 1, which.max)])
f_manual <- file.path(OUT, "cluster_annotation_manual.csv")
if (file.exists(f_manual)) {
  man <- utils::read.csv(f_manual)
  ann$celltype <- man$celltype[match(ann$cluster, man$cluster)]
} else {
  ann$celltype <- ann$celltype_auto
  utils::write.csv(ann, file.path(OUT, "cluster_annotation_auto.csv"), row.names = FALSE)
}
sc$celltype <- ann$celltype[match(as.character(sc$seurat_clusters), as.character(ann$cluster))]
for (cc in mk_cols) sc[[cc]] <- NULL

save_gg(DotPlot(sc, features = unique(unlist(mk)), group.by = "seurat_clusters") + RotatedAxis(),
        file.path(OUT, "DotPlot_markers_by_cluster.pdf"), 18, 10)
save_gg(DotPlot(sc, features = unique(unlist(mk)), group.by = "celltype") + RotatedAxis(),
        file.path(OUT, "Fig2_DotPlot_markers_by_celltype.pdf"), 16, 6)
save_gg(DimPlot(sc, group.by = "celltype", label = TRUE, repel = TRUE) +
          DimPlot(sc, group.by = "sample") + DimPlot(sc, group.by = "tissue"),
        file.path(OUT, "Fig2_UMAP_celltype_sample_tissue.pdf"), 21, 6)

# =============================================================================
# 4. copyKAT：逐个样本串行运行（释放 Seurat 大对象后运行，每样本上皮细胞上限 350 个，兼顾精度与 8GB 内存）
# =============================================================================
REF_TYPES <- c("T_cell", "NK", "B_cell", "Plasma", "Myeloid")
N_REF <- 150
N_EPI_MAX <- 350
ck_dir <- file.path(OUT, "copykat"); dir.create(ck_dir, showWarnings = FALSE)
ck_file <- function(smp) file.path(ck_dir, paste0(smp, "_copykat_prediction.txt"))

ck_plan <- list()
for (smp in SC_SAMPLES$sample) {
  epi <- colnames(sc)[sc$sample == smp & sc$celltype == "Epithelial"]
  if (length(epi) < 30) { log_msg(smp, "：上皮细胞不足 30 个（", length(epi), "），跳过 copyKAT"); next }
  set.seed(PARAM$seed)
  epi_use <- if (length(epi) > N_EPI_MAX) sample(epi, N_EPI_MAX) else epi
  own <- colnames(sc)[sc$sample == smp & sc$celltype %in% REF_TYPES]
  ref <- if (length(own) >= N_REF) sample(own, N_REF) else {
    pool <- setdiff(colnames(sc)[sc$celltype %in% REF_TYPES], own)
    c(own, sample(pool, min(length(pool), N_REF - length(own))))
  }
  ck_plan[[smp]] <- list(epi = epi_use, all_epi = epi, ref = ref, n_own = length(own),
                         tissue = SC_SAMPLES$tissue[SC_SAMPLES$sample == smp])
}

utils::write.csv(data.frame(
  sample = names(ck_plan),
  tissue = vapply(ck_plan, function(p) p$tissue, ""),
  n_epithelial = vapply(ck_plan, function(p) length(p$epi), 1L),
  n_reference = vapply(ck_plan, function(p) length(p$ref), 1L),
  n_reference_same_sample = vapply(ck_plan, function(p) p$n_own, 1L)
), file.path(OUT, "copykat_reference_cells.csv"), row.names = FALSE)

todo <- names(ck_plan)[!file.exists(vapply(names(ck_plan), ck_file, ""))]
if (length(todo) && requireNamespace("copykat", quietly = TRUE)) {
  need("copykat")
  cnt <- LayerData(sc, assay = "RNA", layer = "counts")
  f_in <- stats::setNames(file.path(ck_dir, paste0(todo, "_input_counts.rds")), todo)
  for (smp in todo) {
    m <- cnt[, c(ck_plan[[smp]]$epi, ck_plan[[smp]]$ref)]
    saveRDS(m[Matrix::rowSums(m) > 0, ], f_in[[smp]])
  }
  rm(cnt, m)
  f_tmp <- file.path(DIR_PRO, "scRNA_annotated_tmp.rds")
  saveRDS(sc, f_tmp, compress = FALSE); rm(sc); invisible(gc())

  for (smp in todo) {
    raw <- as.matrix(readRDS(f_in[[smp]]))
    colnames(raw) <- make.names(colnames(raw))
    ref <- ck_plan[[smp]]$ref
    log_msg("copyKAT：", smp, "（", length(ck_plan[[smp]]$epi), " 个上皮 + ", length(ref), " 个参照）")
    old <- setwd(ck_dir)
    ok <- tryCatch({
      invisible(copykat::copykat(
        rawmat = raw, id.type = "S", ngene.chr = 5, win.size = 25, KS.cut = 0.1, sam.name = smp,
        distance = "euclidean", norm.cell.names = make.names(ref),
        output.seg = "FALSE", plot.genes = "FALSE", genome = "hg20",
        n.cores = 1L
      ))
      TRUE
    }, error = function(e) { log_msg("copyKAT 失败：", smp, "：", conditionMessage(e)); FALSE })
    setwd(old)
    rm(raw); invisible(gc())
    if (ok && file.exists(ck_file(smp))) {
      unlink(f_in[[smp]])
      big <- list.files(ck_dir, pattern = paste0("^", smp, "_copykat_.*(gene_by_cell|bin_by_cell)\\.txt$"), full.names = TRUE)
      unlink(big)
    }
  }
  sc <- readRDS(f_tmp); unlink(f_tmp)
}

ck <- do.call(rbind, lapply(names(ck_plan), function(smp) {
  if (!file.exists(ck_file(smp))) return(NULL)
  pr <- utils::read.delim(ck_file(smp))
  epi <- ck_plan[[smp]]$epi
  out <- data.frame(cell = epi[match(pr$cell.names, make.names(epi))], copykat = pr$copykat.pred)
  out[!is.na(out$cell), ]
}))
sc$copykat <- if (is.null(ck)) NA_character_ else ck$copykat[match(colnames(sc), ck$cell)]

is_te <- sc$celltype == "Epithelial" & sc$tissue == "Tumor"
is_ne <- sc$celltype == "Epithelial" & sc$tissue == "Normal"
sc$cell_class <- sc$celltype
sc$cell_class[is_ne] <- "Normal_epi"
sc$cell_class[is_te] <- if (MALIGNANT_RULE == "all_tumor_epi") "Malignant" else
  ifelse(sc$copykat[is_te] %in% "aneuploid", "Malignant",
         ifelse(sc$copykat[is_te] %in% "diploid", "Tumor_diploid_epi", "Tumor_epi_undefined"))

# 输出肿瘤与正常上皮 copyKAT 统计表
ck_lab_t <- ifelse(sc$copykat[is_te] %in% "aneuploid", "aneuploid",
                   ifelse(sc$copykat[is_te] %in% "diploid", "diploid", "not_analysed_or_undefined"))
ck_tab <- as.data.frame.matrix(table(sc$sample[is_te],
                                     factor(ck_lab_t, levels = c("aneuploid", "diploid", "not_analysed_or_undefined"))))
ck_tab$aneuploid_frac <- round(ck_tab$aneuploid / pmax(1, ck_tab$aneuploid + ck_tab$diploid), 3)
utils::write.csv(ck_tab, file.path(OUT, "copykat_by_sample.csv"))

if (any(is_ne)) {
  ck_lab_n <- ifelse(sc$copykat[is_ne] %in% "aneuploid", "aneuploid",
                     ifelse(sc$copykat[is_ne] %in% "diploid", "diploid", "not_analysed_or_undefined"))
  ck_norm <- as.data.frame.matrix(table(sc$sample[is_ne],
                                        factor(ck_lab_n, levels = c("aneuploid", "diploid", "not_analysed_or_undefined"))))
  ck_norm_df <- data.frame(
    sample = rownames(ck_norm),
    n_epithelial = rowSums(ck_norm),
    n_aneuploid = ck_norm$aneuploid,
    n_diploid = ck_norm$diploid,
    n_undefined = ck_norm$not_analysed_or_undefined,
    aneuploid_frac = round(ck_norm$aneuploid / pmax(1, ck_norm$aneuploid + ck_norm$diploid), 3)
  )
  utils::write.csv(ck_norm_df, file.path(OUT, "copykat_normal_control.csv"), row.names = FALSE)
}

# =============================================================================
# 5. UCell 质膜张力评分（MT_core / MT_all / MT_stiff / MT_buffer / MTI）
# =============================================================================
sc <- UCell::AddModuleScore_UCell(
  sc,
  features = GS[c("MT_core", "MT_all", "MT_cortex", "MT_BAR", "MT_channel", "MT_ECM", "MT_stiff", "MT_buffer")],
  maxRank = 2500, name = "_UCell", BPPARAM = BiocParallel::SerialParam()
)
sc$MTI <- as.numeric(scale(sc$MT_stiff_UCell) - scale(sc$MT_buffer_UCell))
# 同步保存兼容别名字段（便于下游模块统一调用）
sc$RL_core_UCell <- sc$MT_core_UCell
sc$RL_all_UCell  <- sc$MT_all_UCell
sc$RHI           <- sc$MTI

save_gg(FeaturePlot(sc, c("MT_core_UCell", "MTI"), order = TRUE),
        file.path(OUT, "Fig2_UMAP_MT_score.pdf"), 12, 5)
save_gg(VlnPlot(sc, "MT_core_UCell", group.by = "cell_class", pt.size = 0) + NoLegend(),
        file.path(OUT, "Fig2_Vln_MT_by_cellclass.pdf"), 10, 5)

md <- sc@meta.data
agg <- stats::aggregate(
  cbind(MT_core_UCell, MTI) ~ sample + cell_class,
  data = md[md$cell_class %in% c("Malignant", "Normal_epi"), ], FUN = mean
)
utils::write.csv(agg, file.path(OUT, "MT_score_per_sample_class.csv"), row.names = FALSE)

# =============================================================================
# 6. 恶性上皮细胞样本内中位数分组 + 伪 bulk edgeR 配对差异分析
# =============================================================================
mal <- subset(sc, subset = cell_class == "Malignant")
n_mal <- table(mal$sample)
keep_s <- names(n_mal)[n_mal >= 2 * PARAM$sc_min_cells_pb]
log_msg("纳入高/低分组的肿瘤样本：", paste(keep_s, collapse = ", "))
mal <- subset(mal, subset = sample %in% keep_s)
grp <- character(ncol(mal))
for (s in keep_s) {
  i <- mal$sample == s
  grp[i] <- ifelse(mal$MT_core_UCell[i] >= stats::median(mal$MT_core_UCell[i]), "High", "Low")
}
mal$MT_group <- factor(grp, levels = c("Low", "High"))
mal$RL_group <- mal$MT_group

pb <- AggregateExpression(mal, assays = "RNA", group.by = c("sample", "MT_group"), return.seurat = FALSE)$RNA
pb <- as.matrix(pb)
pb_smp <- sub("_[^_]+$", "", colnames(pb))
pb_grp <- factor(sub("^.*_", "", colnames(pb)), levels = c("Low", "High"))
y <- edgeR::DGEList(pb, group = pb_grp)
y <- y[edgeR::filterByExpr(y, group = pb_grp), , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y)
design <- stats::model.matrix(~ factor(pb_smp) + pb_grp)
y <- edgeR::estimateDisp(y, design)
fit <- edgeR::glmQLFit(y, design)
deg_pb <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = "pb_grpHigh"), n = Inf)$table
deg_pb$gene <- rownames(deg_pb)
utils::write.csv(deg_pb, file.path(OUT, "pseudobulk_DE_High_vs_Low.csv"), row.names = FALSE)

sc_deg <- deg_pb$gene[deg_pb$FDR < PARAM$sc_deg_fdr & abs(deg_pb$logFC) > PARAM$sc_deg_lfc]
writeLines(sc_deg, file.path(OUT, "sc_DEG_genes.txt"))
log_msg("伪 bulk 差异基因：", length(sc_deg), " 个（FDR < ", PARAM$sc_deg_fdr, "，|log2FC| > ", PARAM$sc_deg_lfc, "）")

deg_pb$sig <- ifelse(deg_pb$FDR < PARAM$sc_deg_fdr & deg_pb$logFC > PARAM$sc_deg_lfc, "Up",
                     ifelse(deg_pb$FDR < PARAM$sc_deg_fdr & deg_pb$logFC < -PARAM$sc_deg_lfc, "Down", "NS"))
top_lab <- utils::head(deg_pb[deg_pb$sig != "NS", ][order(deg_pb$FDR[deg_pb$sig != "NS"]), ], 20)
p <- ggplot2::ggplot(deg_pb, ggplot2::aes(logFC, -log10(FDR), color = sig)) +
  ggplot2::geom_point(size = 0.6) +
  ggplot2::scale_color_manual(values = c(Up = "#D7301F", Down = "#2C7FB8", NS = "grey75")) +
  ggrepel::geom_text_repel(data = top_lab, ggplot2::aes(label = gene), size = 3, max.overlaps = 30) +
  ggplot2::theme_classic() + ggplot2::labs(title = "Malignant ccRCC cells: MT-high vs MT-low (pseudobulk)")
save_gg(p, file.path(OUT, "Fig3_volcano_pseudobulk.pdf"), 6, 5)

# =============================================================================
# 7. Hallmark 通路（样本内 High - Low 配对 Wilcoxon）
# =============================================================================
H <- get_msig("H")
mal <- UCell::AddModuleScore_UCell(mal, features = H, name = "_UCell", BPPARAM = BiocParallel::SerialParam())
hcols <- paste0(names(H), "_UCell")
hm <- stats::aggregate(mal@meta.data[, hcols], by = list(sample = mal$sample, group = mal$MT_group), FUN = mean)
hall <- do.call(rbind, lapply(hcols, function(h) {
  w <- stats::reshape(hm[, c("sample", "group", h)], idvar = "sample", timevar = "group", direction = "wide")
  d <- w[[paste0(h, ".High")]] - w[[paste0(h, ".Low")]]
  data.frame(pathway = sub("^HALLMARK_", "", sub("_UCell$", "", h)), mean_diff = mean(d, na.rm = TRUE),
             p = tryCatch(stats::wilcox.test(d, mu = 0)$p.value, error = function(e) NA_real_))
}))
hall$FDR <- stats::p.adjust(hall$p, "BH")
hall <- hall[order(hall$mean_diff), ]
utils::write.csv(hall, file.path(OUT, "Hallmark_High_vs_Low_paired.csv"), row.names = FALSE)

# =============================================================================
# 8. CellChat：高/低膜张力恶性细胞与肿瘤微环境基质/血管/免疫细胞互作
# =============================================================================
if (requireNamespace("CellChat", quietly = TRUE)) {
  others <- colnames(sc)[sc$tissue == "Tumor" &
                           sc$cell_class %in% c("T_cell", "NK", "B_cell", "Plasma", "Myeloid", "Mast",
                                                "Fibroblast", "Myofib_Peri", "Endothelial")]
  cc_cells <- c(colnames(mal), others)
  meta <- data.frame(
    labels = c(paste0("Malig_", mal$MT_group), sc$cell_class[match(others, colnames(sc))]),
    samples = factor(sc$sample[match(cc_cells, colnames(sc))]),
    row.names = cc_cells
  )
  cellchat <- CellChat::createCellChat(object = LayerData(sc, assay = "RNA", layer = "data")[, cc_cells],
                                       meta = meta, group.by = "labels")
  cellchat@DB <- CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean")
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  saveRDS(cellchat, file.path(DIR_PRO, "scRNA_cellchat.rds"))

  comm <- CellChat::subsetCommunication(cellchat)
  utils::write.csv(comm[comm$source %in% c("Malig_High", "Malig_Low"), ],
                   file.path(OUT, "CellChat_from_malignant.csv"), row.names = FALSE)
  net <- data.frame(
    source = rep(rownames(cellchat@net$count), ncol(cellchat@net$count)),
    target = rep(colnames(cellchat@net$count), each = nrow(cellchat@net$count)),
    count  = as.vector(cellchat@net$count),
    weight = as.vector(cellchat@net$weight)
  )
  utils::write.csv(net[net$source %in% c("Malig_High", "Malig_Low"), ],
                   file.path(OUT, "CellChat_count_weight_from_malignant.csv"), row.names = FALSE)
}

# =============================================================================
# 9. 保存单细胞整合对象与元数据
# =============================================================================
saveRDS(sc, file.path(DIR_PRO, "scRNA_integrated.rds"))
saveRDS(mal, file.path(DIR_PRO, "scRNA_malignant.rds"))
meta_out <- sc@meta.data[, intersect(
  c("sample", "gsm", "dataset", "tissue", "seurat_clusters", "celltype", "copykat", "cell_class",
    "MT_core_UCell", "MT_all_UCell", "MT_cortex_UCell", "MT_BAR_UCell", "MT_channel_UCell", "MT_ECM_UCell",
    "MTI", "RL_core_UCell", "RL_all_UCell", "RHI"),
  colnames(sc@meta.data)
)]
meta_out$MT_group <- mal$MT_group[match(rownames(meta_out), colnames(mal))]
meta_out$RL_group <- meta_out$MT_group
utils::write.csv(meta_out, file.path(OUT, "cell_metadata.csv"))
save_session(STEP)
log_msg("03_scRNA 完成")
