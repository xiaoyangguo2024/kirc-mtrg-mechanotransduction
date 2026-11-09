# =============================================================================
# 11_figures.R —— KIRC_MTRG 论文主图与补充图：统一出版级主题重绘
#   每张图输出 600 dpi PNG（results/figures/FigureN.png）和矢量 PDF（FigureN.pdf）
#   图中标注的统计量同步写入 results/figures/FigureN_stats.csv，确保正文与图表完全一致
# =============================================================================
source("R/00_config.R")
need(c("ggplot2", "patchwork", "ragg", "survival", "scales", "ggrepel"))
FIG <- res_dir("figures")
args <- commandArgs(trailingOnly = TRUE)
WHICH <- if (length(args)) args else c(1:9, paste0("S", 1:3))

FONT <- "Arial"
BASE <- 7
theme_pub <- function(base = BASE) {
  ggplot2::theme_classic(base_size = base, base_family = FONT) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = "black", size = base - 0.5),
      axis.line = ggplot2::element_line(linewidth = 0.3),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      plot.title = ggplot2::element_text(size = base + 0.5, hjust = 0.5),
      plot.tag = ggplot2::element_text(size = base + 4, face = "bold", family = FONT),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base, face = "bold"),
      legend.key.size = ggplot2::unit(3, "mm"),
      legend.text = ggplot2::element_text(size = base - 0.5),
      legend.title = ggplot2::element_text(size = base - 0.5)
    )
}
ggplot2::theme_set(theme_pub())
C_HI <- "#C0392B"; C_LO <- "#2E6DA4"; C_GREY <- "grey55"
PAL2 <- c(High = C_HI, Low = C_LO)

fmt_p <- function(p, digits = 2) {
  ifelse(is.na(p), "NA", ifelse(p < 1e-4, "<0.0001", formatC(signif(p, digits), format = "fg", digits = digits)))
}
lab_p <- function(p) ifelse(!is.na(p) & p < 1e-4, "P < 0.0001", paste0("P = ", fmt_p(p)))
exact_p <- function(p) paste0("P = ", signif(p, 2))
lab_q <- function(q) ifelse(!is.na(q) & q < 1e-4, "FDR < 0.0001", paste0("FDR = ", fmt_p(q)))

STATS <- list()
add_stat <- function(fig, panel, what, value, detail = "") {
  STATS[[length(STATS) + 1]] <<- data.frame(
    figure = fig, panel = panel, statistic = what,
    value = as.character(value), detail = detail
  )
}

save_fig <- function(p, n, w = 180, h = 200) {
  f <- file.path(FIG, paste0("Figure", n))
  tmp <- tempfile(fileext = ".png")
  ragg::agg_png(tmp, width = w, height = h, units = "mm", res = 600, background = "white")
  print(p); grDevices::dev.off()
  stopifnot(file.copy(tmp, paste0(f, ".png"), overwrite = TRUE)); unlink(tmp)
  grDevices::cairo_pdf(paste0(f, ".pdf"), width = w / 25.4, height = h / 25.4, family = FONT)
  print(p); grDevices::dev.off()
  st <- do.call(rbind, Filter(function(x) x$figure == n, STATS))
  if (!is.null(st)) utils::write.csv(st, paste0(f, "_stats.csv"), row.names = FALSE)
  log_msg("Figure ", n, " 已保存（", w, " × ", h, " mm）")
}

km_gg <- function(time, status, group, title = NULL, xlab = "Years", ylab = "Survival probability",
                  cols = PAL2, years = TRUE) {
  t <- if (years) time / 365.25 else time
  d <- data.frame(t = t, s = status, g = factor(group, levels = names(cols)))
  d <- d[stats::complete.cases(d), ]
  fit <- survival::survfit(survival::Surv(t, s) ~ g, data = d)
  grp <- sub("^g=", "", names(fit$strata))
  curve <- do.call(rbind, lapply(seq_along(grp), function(i) {
    ss <- fit[i]
    data.frame(g = grp[i], t = c(0, ss$time), surv = c(1, ss$surv), cens = c(0, ss$n.censor))
  }))
  lr <- survival::survdiff(survival::Surv(t, s) ~ g, data = d)
  p <- stats::pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  cx <- summary(survival::coxph(survival::Surv(t, s) ~ relevel(g, ref = names(cols)[2]), data = d))$conf.int
  ntab <- table(d$g)
  lab <- sprintf("HR %.2f (%.2f-%.2f)\nLog-rank %s", cx[1, 1], cx[1, 3], cx[1, 4], lab_p(p))
  curve$g <- factor(curve$g, levels = names(cols),
                    labels = sprintf("%s (n = %d)", names(cols), as.integer(ntab[names(cols)])))
  cl <- stats::setNames(cols, levels(curve$g))
  late <- curve$surv[curve$t >= 2 / 3 * max(curve$t)]
  hr_top <- length(late) > 0 && max(late) < 0.75
  gg <- ggplot2::ggplot(curve, ggplot2::aes(t, surv, colour = g)) +
    ggplot2::geom_step(linewidth = 0.45) +
    ggplot2::geom_point(data = curve[curve$cens > 0, ], shape = 3, size = 0.6, stroke = 0.3) +
    ggplot2::annotate("text", x = Inf, y = if (hr_top) 1 else 0.02, hjust = 1.03, vjust = if (hr_top) 1 else 0, label = lab,
                      size = (BASE - 1) / ggplot2::.pt, family = FONT, lineheight = 0.9) +
    ggplot2::scale_colour_manual(values = cl, name = NULL) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme(
      legend.position = "inside", legend.position.inside = c(0.01, 0.01), legend.justification = c(0, 0),
      legend.background = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = BASE - 1.5),
      legend.key.height = ggplot2::unit(2.5, "mm")
    )
  attr(gg, "stats") <- c(HR = cx[1, 1], lower = cx[1, 3], upper = cx[1, 4], p_logrank = p)
  gg
}

forest_gg <- function(df, label = "label", group = NULL, xlab = "Hazard ratio per SD (95% CI)", cols = NULL) {
  df$label <- factor(df[[label]], levels = rev(unique(df[[label]])))
  aes_c <- if (is.null(group)) ggplot2::aes(HR, label) else ggplot2::aes(HR, label, colour = .data[[group]])
  g <- ggplot2::ggplot(df, aes_c) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), width = 0.25, linewidth = 0.4, orientation = "y",
                           position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(size = 1.3, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_x_log10() + ggplot2::labs(x = xlab, y = NULL)
  if (!is.null(cols)) g <- g + ggplot2::scale_colour_manual(values = cols, name = NULL)
  g
}

# =============================================================================
# Figure 1：研究设计流程图
# =============================================================================
fig1 <- function() {
  W <- 180; H <- 121
  FILL <- c(sc = "#F7E1DD", bulk = "#DDE9F5", sp = "#E0F0E1", gs = "#EDEDED", mod = "#FBF0D5", ev = "#FFFFFF")
  box <- function(id, x0, x1, y1, h, fill, title, body)
    data.frame(id = id, x0 = x0, x1 = x1, y1 = y1, y0 = y1 - h, fill = fill, title = title, body = body)
  bx <- rbind(
    box("gs", 2, 178, 120, 9.5, "gs", "Plasma membrane tension regulator gene sets (MTRGs)",
        paste0("MT_core: 97 regulators (Cortex 37, BAR/caveolae 27, Mechanochannels 13, ECM-FA 20)  ·  ",
               "MT_all: 125  ·  MTI = stiffening (70) − buffering (27)")),
    box("sc", 2, 58, 105.5, 15, "sc", "Single-cell RNA-seq", "GSE159115\n7 ccRCC tumours, 6 normal kidneys\n(14,208 quality-controlled cells)"),
    box("bulk", 62, 118, 105.5, 15, "bulk", "Bulk transcriptomes",
        "TCGA-KIRC (533 tumours, 72 normal)\n+ 7 GEO and ICI cohorts\n1,505 primary/metastatic samples"),
    box("sp", 122, 178, 105.5, 15, "sp", "Spatial transcriptomics",
        "GSE210041\n2 ccRCC tumour sections\n(10x Genomics Visium)"),
    box("sc2", 2, 58, 85.5, 17.5, "sc", "Tumour epithelial cells",
        "Harmony integration, annotation\ncopyKAT with normal-kidney control\nUCell scores; MT-high vs -low cells\n(pseudobulk DE, G1 check, CellChat)"),
    box("bulk2", 62, 118, 85.5, 17.5, "bulk", "Membrane tension scores",
        "ssGSEA scores vs EMT, hypoxia,\nproliferation, VHL/PBRM1/SETD2/BAP1\nmutations and survival across cohorts"),
    box("sp2", 122, 178, 85.5, 17.5, "sp", "Spot-level programmes",
        "MT_core score vs EMT, angiogenesis,\nhypoxia and immune programmes\n(adjusted for UMI and proliferation)"),
    box("cand", 2, 118, 63, 12.5, "mod", "Candidate genes",
        "(MT_all or WGCNA key-module genes) also differentially expressed in scRNA-seq\nand in TCGA tumour vs normal kidney; univariable Cox regression in TCGA"),
    box("mime", 2, 170, 45.5, 24, "mod", "Mime: 10 survival algorithms, 101 combinations",
        "Within-cohort z-scores; fixed-direction C-index → Membrane Tension-Related Risk Score (MTRS)"),
    box("train", 4, 56.67, 36, 12.5, "ev", "Training", "TCGA-KIRC (OS, n = 531)"),
    box("tune", 59.67, 112.33, 36, 12.5, "ev", "Tuning (model selection)", "GSE29609 (OS), GSE22541 (DFS)"),
    box("test", 115.33, 168, 36, 12.5, "ev", "External & clinical cohorts",
        "IMvigor210, GSE176307 (OS & ICI)\nGSE67501 (RCC Nivolumab)\nGSE53757, GSE73731 (Stage/Grade)"),
    box("spec", 2, 41.75, 16.5, 15, "ev", "Specificity",
        "1,000 random gene sets × 2\nproliferation adjustment\n10 published ccRCC signatures"),
    box("clin", 44.75, 84.5, 16.5, 15, "ev", "Clinical model", "Multivariable Cox\nnomogram, calibration\ndecision curves"),
    box("ther", 87.5, 127.25, 16.5, 15, "ev", "TME & therapy",
        "ICI: GSE67501, IMvigor210, GSE176307\nIOBR (7 methods), TIDE\noncoPredict (TKI & mTORi)"),
    box("genes", 130.25, 170, 16.5, 15, "mod", "Model genes",
        "SHAP decomposition\n72-pair tumour vs normal\nsingle-cell and spatial Visium")
  )
  bx$xm <- (bx$x0 + bx$x1) / 2
  fx <- bx$xm[match(c("spec", "clin", "ther", "genes"), bx$id)]
  ar <- data.frame(
    x    = c(30, 90, 150, 30, 90, 150, 30, 90, 60, fx),
    y    = c(110.5, 110.5, 110.5, 90.5, 90.5, 90.5, 68, 68, 50.5, rep(21.5, 4)),
    xend = c(30, 90, 150, 30, 90, 150, 30, 90, 60, fx),
    yend = c(105.5, 105.5, 105.5, 85.5, 85.5, 85.5, 63, 63, 45.5, rep(16.5, 4))
  )
  sp_line <- data.frame(x = c(175, 175), y = c(68, 9))
  sp_arr  <- data.frame(x = 175, y = 9, xend = 170, yend = 9)
  pt <- ggplot2::.pt
  arw <- ggplot2::arrow(length = ggplot2::unit(1.4, "mm"), type = "closed")
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = bx, ggplot2::aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill),
                       colour = "grey35", linewidth = 0.3) +
    ggplot2::geom_text(data = bx, ggplot2::aes(x = xm, y = y1 - 1.6, label = title), vjust = 1, fontface = "bold",
                       size = BASE / pt, family = FONT) +
    ggplot2::geom_text(data = bx, ggplot2::aes(x = xm, y = y1 - 5, label = body), vjust = 1, lineheight = 0.95,
                       size = (BASE - 0.5) / pt, family = FONT) +
    ggplot2::geom_segment(data = ar, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.35,
                          colour = "grey30", arrow = arw) +
    ggplot2::geom_path(data = sp_line, ggplot2::aes(x, y), linewidth = 0.35, colour = "grey30", linetype = "22") +
    ggplot2::geom_segment(data = sp_arr, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.35,
                          colour = "grey30", linetype = "22", arrow = arw) +
    ggplot2::scale_fill_manual(values = FILL, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, H), expand = c(0, 0)) +
    ggplot2::coord_fixed() + ggplot2::theme_void()
  save_fig(p, 1, W, H)
}

# =============================================================================
# 单细胞数据缓存与配色（GSE159115）
# =============================================================================
CT_NAMES <- c(
  Malignant   = "Malignant epithelial",
  Normal_epi  = "Normal epithelial",
  Endothelial = "Endothelial",
  Myofib_Peri = "Myofibroblast/pericyte",
  T_cell      = "T/NK cell",
  Myeloid     = "Myeloid",
  Mast        = "Mast cell"
)
CT_COLS <- c(
  `Malignant epithelial`   = "#C0392B",
  `Normal epithelial`      = "#E59866",
  Endothelial              = "#B7950B",
  `Myofibroblast/pericyte` = "#27AE60",
  `T/NK cell`              = "#2E6DA4",
  Myeloid                  = "#8E44AD",
  `Mast cell`              = "#A04000"
)
MK_SHOW <- list(
  `Malignant epithelial`   = c("CA9", "NDUFA4L2", "EPCAM", "KRT8", "PAX8"),
  `Normal epithelial`      = c("ALDOB", "GPX3", "UMOD", "SLC12A1"),
  Endothelial              = c("PECAM1", "VWF", "KDR", "PLVAP"),
  `Myofibroblast/pericyte` = c("ACTA2", "PDGFRB", "RGS5", "TAGLN"),
  `T/NK cell`              = c("CD3D", "CD3E", "NKG7", "GNLY"),
  Myeloid                  = c("CD68", "CD14", "LYZ", "C1QA"),
  `Mast cell`              = c("TPSAB1", "KIT", "CPA3")
)
SC <- NULL
load_sc <- function() {
  if (is.null(SC)) {
    need(c("Seurat", "Matrix"))
    sc <- readRDS(file.path(DIR_PRO, "scRNA_integrated.rds"))
    md <- sc@meta.data
    emb <- SeuratObject::Embeddings(sc, "umap")
    md$UMAP1 <- emb[rownames(md), 1]; md$UMAP2 <- emb[rownames(md), 2]
    md$ct <- factor(CT_NAMES[md$cell_class], levels = intersect(CT_NAMES, CT_NAMES[unique(md$cell_class)]))
    g <- intersect(unique(unlist(MK_SHOW)), rownames(sc))
    ex <- SeuratObject::LayerData(sc, assay = "RNA", layer = "data")[g, rownames(md), drop = FALSE]
    rm(sc); invisible(gc())
    SC <<- list(md = md, mk = ex)
  }
  SC
}

# =============================================================================
# Figure 2：单细胞全景、copyKAT 与 MT_core 评分（GSE159115）
# =============================================================================
fig2 <- function() {
  s <- load_sc(); md <- s$md
  set.seed(PARAM$seed); md <- md[sample(nrow(md)), ]
  cen <- stats::aggregate(cbind(UMAP1, UMAP2) ~ ct, data = md, FUN = stats::median)
  add_stat(2, "A", "cells after QC", nrow(md),
           sprintf("%d samples (%d tumour, %d normal)", length(unique(md$sample)),
                   length(unique(md$sample[md$tissue == "Tumor"])),
                   length(unique(md$sample[md$tissue == "Normal"]))))
  for (k in levels(md$ct)) add_stat(2, "A", paste("cells:", k), sum(md$ct == k), sprintf("%.1f%%", 100 * mean(md$ct == k)))
  pA <- ggplot2::ggplot(md, ggplot2::aes(UMAP1, UMAP2, colour = ct)) +
    ggplot2::geom_point(size = 0.25, stroke = 0) +
    ggrepel::geom_text_repel(data = cen, ggplot2::aes(label = ct), colour = "black", size = BASE / ggplot2::.pt, family = FONT,
                             seed = 1, box.padding = 0.5, min.segment.length = 0, segment.size = 0.2,
                             bg.color = "white", bg.r = 0.12) +
    ggplot2::scale_colour_manual(values = CT_COLS, guide = "none") + ggplot2::coord_fixed() +
    ggplot2::labs(title = sprintf("%s cells, %d samples (GSE159115)", format(nrow(md), big.mark = ","), length(unique(md$sample)))) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

  md$tissue_lab <- factor(ifelse(md$tissue == "Tumor", "ccRCC tumour", "Benign adjacent"), levels = c("ccRCC tumour", "Benign adjacent"))
  pB <- ggplot2::ggplot(md[order(md$tissue_lab == "Benign adjacent"), ], ggplot2::aes(UMAP1, UMAP2, colour = tissue_lab)) +
    ggplot2::geom_point(size = 0.25, stroke = 0) +
    ggplot2::scale_colour_manual(values = c(`ccRCC tumour` = C_HI, `Benign adjacent` = C_LO), name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2))) + ggplot2::coord_fixed() +
    ggplot2::labs(title = "Tissue of origin") +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0))

  ex <- s$mk; ct <- s$md$ct
  mk <- unlist(MK_SHOW); mk <- mk[mk %in% rownames(ex)]
  dp <- do.call(rbind, lapply(levels(ct), function(k) {
    m <- ex[mk, which(ct == k), drop = FALSE]
    data.frame(ct = k, gene = mk, avg = Matrix::rowMeans(m), pct = Matrix::rowMeans(m > 0))
  }))
  dp$avg_s <- stats::ave(dp$avg, dp$gene, FUN = function(v) (v - min(v)) / max(diff(range(v)), 1e-9))
  dp$gene <- factor(dp$gene, levels = unique(mk)); dp$ct <- factor(dp$ct, levels = rev(levels(ct)))
  pC <- ggplot2::ggplot(dp, ggplot2::aes(gene, ct)) +
    ggplot2::geom_point(ggplot2::aes(size = pct, colour = avg_s)) +
    ggplot2::scale_size_area(max_size = 2.6, name = "Fraction\nexpressing", breaks = c(0.25, 0.5, 0.75)) +
    ggplot2::scale_colour_gradient(low = "grey90", high = C_HI, name = "Scaled\nexpression", breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
                   legend.box = "horizontal")

  te <- md[md$cell_class %in% c("Malignant", "Normal_epi") & md$tissue == "Tumor", ]
  te$call <- ifelse(te$copykat %in% "aneuploid", "Aneuploid",
                    ifelse(te$copykat %in% "diploid", "Diploid", "Not analysed/undefined"))
  cb <- as.data.frame(table(sample = te$sample, call = te$call))
  cb$set <- "Tumour"
  f_ctl <- file.path(DIR_RES, "03_scRNA", "copykat_normal_control.csv")
  if (file.exists(f_ctl)) {
    ctl <- utils::read.csv(f_ctl)
    cb <- rbind(cb, data.frame(
      sample = rep(ctl$sample, 3),
      call = rep(c("Aneuploid", "Diploid", "Not analysed/undefined"), each = nrow(ctl)),
      Freq = c(ctl$n_aneuploid, ctl$n_diploid, ctl$n_undefined), set = "Normal"
    ))
    for (i in seq_len(nrow(ctl))) add_stat(2, "D", paste("copyKAT negative control, aneuploid fraction:", ctl$sample[i]),
                                           ctl$aneuploid_frac[i], paste(ctl$n_epithelial[i], "normal epithelial cells"))
  }
  tot <- stats::aggregate(Freq ~ sample + set, data = cb, FUN = sum)
  cb$frac <- cb$Freq / tot$Freq[match(paste(cb$sample, cb$set), paste(tot$sample, tot$set))]
  an <- cb[cb$call == "Aneuploid" & cb$set == "Tumour", ]
  for (i in seq_len(nrow(an))) add_stat(2, "D", paste("aneuploid fraction of tumour epithelial cells:", an$sample[i]),
                                        round(an$frac[i], 3), paste(tot$Freq[tot$sample == an$sample[i] & tot$set == "Tumour"], "epithelial cells"))
  add_stat(2, "D", "tumour epithelial cells called aneuploid (all samples)", sum(te$call == "Aneuploid"),
           sprintf("of %d (%.1f%%)", nrow(te), 100 * mean(te$call == "Aneuploid")))
  ord <- tot$sample[order(tot$set != "Tumour", tot$sample)]
  cb$sample <- factor(cb$sample, levels = unique(ord))
  cb$set <- factor(cb$set, levels = c("Tumour", "Normal"))
  pD <- ggplot2::ggplot(cb, ggplot2::aes(sample, frac, fill = call)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_grid(~ set, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = c(Aneuploid = C_HI, Diploid = C_LO, `Not analysed/undefined` = "grey80"), name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Epithelial cells (copyKAT)", title = "copyKAT calls") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 1.5),
                   legend.position = "top", strip.clip = "off")

  pE <- ggplot2::ggplot(md[order(md$MT_core_UCell), ], ggplot2::aes(UMAP1, UMAP2, colour = MT_core_UCell)) +
    ggplot2::geom_point(size = 0.25, stroke = 0) +
    ggplot2::scale_colour_gradientn(colours = c("grey88", "#F5B7B1", C_HI, "#641E16"), name = "MT_core\n(UCell)") +
    ggplot2::coord_fixed() + ggplot2::labs(title = "Membrane tension regulator activity") +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

  cls <- c(Malignant = "Tumour", Normal_epi = "Normal")
  sub_md <- md[md$cell_class %in% names(cls), ]
  agg <- stats::aggregate(MT_core_UCell ~ sample + cell_class, data = sub_md, FUN = mean)
  nn  <- stats::aggregate(MT_core_UCell ~ sample + cell_class, data = sub_md, FUN = length)
  agg <- agg[nn$MT_core_UCell[match(paste(agg$sample, agg$cell_class), paste(nn$sample, nn$cell_class))] >= 10, ]
  agg$cls <- factor(cls[agg$cell_class], levels = cls)
  p_mn <- stats::wilcox.test(agg$MT_core_UCell[agg$cell_class == "Malignant"],
                             agg$MT_core_UCell[agg$cell_class == "Normal_epi"])$p.value
  add_stat(2, "F", "per-sample mean MT_core, tumour vs normal epithelium, Wilcoxon P", signif(p_mn, 3),
           sprintf("%d tumour samples vs %d normal samples", sum(agg$cell_class == "Malignant"), sum(agg$cell_class == "Normal_epi")))
  add_stat(2, "F", "median of per-sample mean MT_core: tumour / normal epithelium",
           sprintf("%.3f / %.3f", stats::median(agg$MT_core_UCell[agg$cell_class == "Malignant"]),
                   stats::median(agg$MT_core_UCell[agg$cell_class == "Normal_epi"])))
  pF <- ggplot2::ggplot(agg, ggplot2::aes(cls, MT_core_UCell)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(colour = cls), size = 1.2, position = ggplot2::position_jitter(width = 0.08, seed = 1)) +
    ggplot2::annotate("text", x = 1.5, y = Inf, vjust = 1.3, label = lab_p(p_mn), size = (BASE - 0.5) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(C_HI, C_LO), cls), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::labs(x = "Epithelial cells", y = "MT_core (per-sample mean)", title = "Samples as the unit")

  pl <- Map(function(g, t) g + ggplot2::labs(tag = t), list(pA, pB, pC, pD, pE, pF), LETTERS[1:6])
  p <- patchwork::wrap_plots(
    patchwork::wrap_elements(full = pl[[1]] | pl[[2]]),
    patchwork::wrap_elements(full = pl[[3]]),
    patchwork::wrap_elements(full = (pl[[4]] | pl[[5]] | pl[[6]]) + patchwork::plot_layout(widths = c(1.35, 1, 0.7))),
    ncol = 1, heights = c(1.05, 0.62, 1)
  )
  save_fig(p, 2, 180, 215)
}

# =============================================================================
# Figure 3：膜张力高/低肿瘤上皮细胞（伪 bulk 差异、Hallmark、细胞周期、CellChat）
# =============================================================================
fig3 <- function(d3 = file.path(DIR_RES, "03_scRNA"), id = 3, cells = "Tumour epithelial cells",
                 all_lab = "All tumour epithelial cells",
                 lab_genes = c("CAV2", "COL4A1", "COL4A2", "NTN4", "ETS1", "GNG7", "RNF152", "ARC",
                               "RAPGEF3", "SUCNR1", "BCAM", "JUP", "SPOCK2", "CTNNAL1", "PHLDA1", "CDKN2B")) {
  GS <- readRDS(file.path(DIR_GS, "MTRG_genesets.rds"))
  de <- utils::read.csv(file.path(d3, "pseudobulk_DE_High_vs_Low.csv"))
  de$sig <- ifelse(de$PValue < 0.05 & de$logFC > 0.3, "Up in MT-high",
                   ifelse(de$PValue < 0.05 & de$logFC < -0.3, "Up in MT-low", "NS"))
  add_stat(id, "A", "pseudobulk DEGs up in MT-high", sum(de$sig == "Up in MT-high"), "P < 0.05, log2FC > 0.3")
  add_stat(id, "A", "pseudobulk DEGs up in MT-low", sum(de$sig == "Up in MT-low"), "P < 0.05, log2FC < -0.3")
  add_stat(id, "A", "DEGs that are MTRGs (MT_all)", sum(de$sig != "NS" & de$gene %in% GS$MT_all),
           sprintf("of %d DEGs; MT_core: %d", sum(de$sig != "NS"), sum(de$sig != "NS" & de$gene %in% GS$MT_core)))
  lab <- de[de$sig != "NS" & de$gene %in% lab_genes, ]
  if (nrow(lab) < 10) {
    rest <- de[de$sig != "NS" & !de$gene %in% lab$gene, ]
    rest <- rest[order(rest$PValue), ]
    lab <- rbind(lab, utils::head(rest, 10 - nrow(lab)))
  }
  pA <- ggplot2::ggplot(de, ggplot2::aes(logFC, -log10(PValue), colour = sig)) +
    ggplot2::geom_point(size = 0.35) +
    ggplot2::geom_vline(xintercept = c(-0.3, 0.3), linetype = 2, linewidth = 0.25, colour = "grey60") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25, colour = "grey60") +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene), size = (BASE - 1.5) / ggplot2::.pt, family = FONT,
                             max.overlaps = Inf, segment.size = 0.2, fontface = "italic", show.legend = FALSE,
                             box.padding = 0.35, point.padding = 0.15, min.segment.length = 0, force = 3, seed = 1) +
    ggplot2::scale_colour_manual(values = c(`Up in MT-high` = C_HI, `Up in MT-low` = C_LO, NS = "grey80"), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.2))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.14))) +
    ggplot2::labs(x = "log2 fold change (MT-high vs -low)", y = "-log10 P value", title = paste0(cells, ", pseudobulk")) +
    ggplot2::theme(legend.position = "bottom") + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2)))

  hb <- utils::read.csv(file.path(d3, "Hallmark_all_vs_G1only.csv"))
  ord <- hb[order(-abs(hb$mean_diff_all)), ]
  keep <- unique(c(ord$pathway[ord$FDR_all < 0.05 | ord$FDR_G1 < 0.05], ord$pathway))
  top <- hb[match(keep[seq_len(min(18, length(keep)))], hb$pathway), ]
  hl <- rbind(
    data.frame(pathway = top$pathway, set = all_lab, diff = top$mean_diff_all, FDR = top$FDR_all),
    data.frame(pathway = top$pathway, set = "G1-phase cells only", diff = top$mean_diff_G1, FDR = top$FDR_G1)
  )
  hl$pathway <- factor(gsub("_", " ", hl$pathway), levels = gsub("_", " ", top$pathway[order(top$mean_diff_all)]))
  hl$sig <- ifelse(!is.na(hl$FDR) & hl$FDR < 0.05, "FDR < 0.05", "NS")
  for (i in seq_len(nrow(hb))) if ((!is.na(hb$FDR_all[i]) && hb$FDR_all[i] < 0.05) || (!is.na(hb$FDR_G1[i]) && hb$FDR_G1[i] < 0.05))
    add_stat(id, "B", paste("Hallmark High-Low difference (all / G1):", hb$pathway[i]),
             sprintf("%.4f / %.4f", hb$mean_diff_all[i], hb$mean_diff_G1[i]),
             sprintf("FDR %s / %s", signif(hb$FDR_all[i], 2), signif(hb$FDR_G1[i], 2)))
  add_stat(id, "B", "Hallmark pathways with FDR < 0.05 (all tumour epithelial / G1 only)",
           sprintf("%d / %d", sum(hb$FDR_all < 0.05, na.rm = TRUE), sum(hb$FDR_G1 < 0.05, na.rm = TRUE)))
  pB <- ggplot2::ggplot(hl, ggplot2::aes(diff, pathway, fill = set, alpha = sig)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = stats::setNames(c(C_HI, "#7F8C8D"), c(all_lab, "G1-phase cells only")), name = NULL) +
    ggplot2::scale_alpha_manual(values = c(`FDR < 0.05` = 1, NS = 0.35), name = NULL) +
    ggplot2::labs(x = "UCell difference, MT-high minus -low\n(per-sample mean)", y = NULL) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical", axis.text.y = ggplot2::element_text(size = BASE - 1.5))

  cc <- utils::read.csv(file.path(d3, "cellcycle_by_sample.csv"))
  pc <- stats::wilcox.test(cc$pct_cycling_High, cc$pct_cycling_Low, paired = TRUE)$p.value
  add_stat(id, "C", "S/G2M cells, MT-high vs -low, median %",
           sprintf("%.1f vs %.1f", stats::median(cc$pct_cycling_High), stats::median(cc$pct_cycling_Low)),
           paste0(exact_p(pc), " (paired Wilcoxon, ", nrow(cc), " samples)"))
  add_stat(id, "C", "within-sample Spearman rho MT_core vs G2M score, median (range)",
           sprintf("%.2f (%.2f to %.2f)", stats::median(cc$rho_RLcore_G2M), min(cc$rho_RLcore_G2M), max(cc$rho_RLcore_G2M)))
  add_stat(id, "C", "within-sample Spearman rho MT_core vs S score, median (range)",
           sprintf("%.2f (%.2f to %.2f)", stats::median(cc$rho_RLcore_S), min(cc$rho_RLcore_S), max(cc$rho_RLcore_S)))
  cl <- rbind(
    data.frame(sample = cc$sample, grp = "Low", pct = cc$pct_cycling_Low),
    data.frame(sample = cc$sample, grp = "High", pct = cc$pct_cycling_High)
  )
  cl$grp <- factor(cl$grp, levels = c("Low", "High"))
  pC <- ggplot2::ggplot(cl, ggplot2::aes(grp, pct)) +
    ggplot2::geom_line(ggplot2::aes(group = sample), colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = grp), size = 1.1) +
    ggplot2::annotate("text", x = 1.5, y = Inf, vjust = 1.3, label = lab_p(pc), size = (BASE - 0.5) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_colour_manual(values = PAL2, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::labs(x = "MT_core group", y = "Cells in S or G2/M phase (%)")
  rl <- rbind(
    data.frame(sample = cc$sample, score = "S phase", rho = cc$rho_RLcore_S),
    data.frame(sample = cc$sample, score = "G2/M phase", rho = cc$rho_RLcore_G2M)
  )
  pD <- ggplot2::ggplot(rl, ggplot2::aes(score, rho)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.35) +
    ggplot2::geom_jitter(width = 0.1, size = 1, colour = C_HI) +
    ggplot2::labs(x = NULL, y = "Spearman ρ with MT_core\n(within sample)")

  f_cc <- file.path(d3, "CellChat_count_weight_from_malignant.csv")
  pE <- NULL
  if (file.exists(f_cc)) {
    cw <- utils::read.csv(f_cc)
    cw <- cw[!cw$target %in% c("Malig_High", "Malig_Low"), ]
    cw$grp <- factor(ifelse(cw$source == "Malig_High", "High", "Low"), levels = c("Low", "High"))
    cw$target <- factor(ifelse(cw$target %in% names(CT_NAMES), CT_NAMES[cw$target], cw$target))
    tot <- stats::aggregate(cbind(count, weight) ~ grp, data = cw, FUN = sum)
    add_stat(id, "E", "CellChat interactions sent by tumour epithelial cells to non-epithelial cells (Low / High)",
             sprintf("%d / %d", tot$count[tot$grp == "Low"], tot$count[tot$grp == "High"]),
             sprintf("total weight %.3f / %.3f", tot$weight[tot$grp == "Low"], tot$weight[tot$grp == "High"]))
    pE <- ggplot2::ggplot(cw, ggplot2::aes(count, target, fill = grp)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::scale_fill_manual(values = PAL2, name = "MT_core group") +
      ggplot2::labs(x = "Inferred interactions (CellChat)", y = NULL, title = paste(cells, "as senders")) +
      ggplot2::theme(legend.position = "bottom")
  }
  pl <- Filter(Negate(is.null), list(pA, pB, pC, pD, pE))
  pl <- Map(function(g, t) g + ggplot2::labs(tag = t), pl, LETTERS[seq_along(pl)])
  top <- patchwork::wrap_elements(full = (pl[[1]] | pl[[2]]) + patchwork::plot_layout(widths = c(1, 1.05)))
  mid <- patchwork::wrap_elements(full = pl[[3]] | pl[[4]])
  p <- if (length(pl) == 5) patchwork::wrap_plots(top, mid, patchwork::wrap_elements(full = pl[[5]]), ncol = 1, heights = c(1.15, 0.6, 0.55)) else
    patchwork::wrap_plots(top, mid, ncol = 1, heights = c(1.15, 0.6))
  save_fig(p, id, 180, 200)
}
figS2 <- function() {
  d_sub <- file.path(DIR_RES, "03_scRNA", "sensitivity_copykat_rule")
  if (file.exists(file.path(d_sub, "pseudobulk_DE_High_vs_Low.csv")) &&
      file.exists(file.path(d_sub, "Hallmark_all_vs_G1only.csv")) &&
      file.exists(file.path(d_sub, "cellcycle_by_sample.csv"))) {
    fig3(d_sub, "S2", cells = "Aneuploid epithelial cells", all_lab = "All aneuploid epithelial cells")
  } else {
    fig3(file.path(DIR_RES, "03_scRNA"), "S2", cells = "Malignant epithelial cells (copyKAT sensitivity)",
         all_lab = "All malignant epithelial cells")
  }
}

# =============================================================================
# Figure 4：Bulk 队列中的质膜张力评分（TCGA-KIRC 与外部队列）
# =============================================================================
fig4 <- function() {
  tcga <- readRDS(file.path(DIR_PRO, "TCGA_KIRC.rds"))
  sc_all <- readRDS(file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds"))
  tum <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
  tc <- merge(tcga$clin, sc_all$TCGA_KIRC, by = "ID")

  # A：TCGA-KIRC 72 对配对肿瘤与癌旁正常肾组织（MT_core 与 MTI）
  pat <- intersect(tc$patient[tc$group == "Tumor"], tc$patient[tc$group == "Normal"])
  pr <- tc[tc$patient %in% pat, c("patient", "group", "MT_core", "MTI")]
  pr <- pr[!duplicated(pr[, c("patient", "group")]), ]
  long <- rbind(
    data.frame(pr[, c("patient", "group")], score = "MT_core", value = pr$MT_core),
    data.frame(pr[, c("patient", "group")], score = "MTI", value = pr$MTI)
  )
  long$group <- factor(long$group, levels = c("Normal", "Tumor"), labels = c("Normal", "Tumour"))
  long$score <- factor(long$score, levels = c("MT_core", "MTI"))
  pv <- vapply(levels(long$score), function(s) {
    w <- stats::reshape(long[long$score == s, c("patient", "group", "value")], idvar = "patient", timevar = "group",
                        direction = "wide")
    stats::wilcox.test(w$value.Tumour, w$value.Normal, paired = TRUE)$p.value
  }, numeric(1))
  unp <- vapply(c("MT_core", "MTI"), function(s) stats::wilcox.test(tc[[s]] ~ tc$group)$p.value, numeric(1))
  for (s in names(pv)) {
    add_stat(4, "A", paste(s, "paired tumour vs normal, Wilcoxon signed-rank P"), signif(pv[[s]], 3), paste(length(pat), "pairs"))
    add_stat(4, "A", paste(s, "all tumour vs normal, Wilcoxon P"), signif(unp[[s]], 3),
             sprintf("%d tumours, %d normals", sum(tc$group == "Tumor"), sum(tc$group == "Normal")))
  }
  pA <- ggplot2::ggplot(long, ggplot2::aes(group, value)) +
    ggplot2::geom_line(ggplot2::aes(group = patient), colour = "grey75", linewidth = 0.2) +
    ggplot2::geom_point(ggplot2::aes(colour = group), size = 0.7) +
    ggplot2::geom_boxplot(ggplot2::aes(colour = group), fill = NA, width = 0.35, outlier.shape = NA, linewidth = 0.35) +
    ggplot2::facet_wrap(~ score, scales = "free_y", nrow = 1) +
    ggplot2::geom_text(data = data.frame(score = factor(names(pv), levels = levels(long$score)), lab = lab_p(pv)),
                       ggplot2::aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, size = (BASE - 0.5) / ggplot2::.pt,
                       inherit.aes = FALSE, family = FONT) +
    ggplot2::scale_colour_manual(values = c(Normal = C_LO, Tumour = C_HI), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
    ggplot2::labs(x = NULL, y = "ssGSEA score", title = sprintf("TCGA-KIRC paired samples (n = %d)", length(pat)))

  # B：相关热图（MT 子模块 vs 机械转导/缺氧/血管生成/代谢/免疫/增殖）
  vars <- c(
    MT_core = "MT_core", MT_cortex = "MT_cortex", MT_BAR = "MT_BAR", MT_channel = "MT_channel",
    MT_ECM = "MT_ECM", MTI = "MTI", EMT = "EMT", TGF_beta = "TGF-β signalling",
    Angiogenesis = "Angiogenesis", Hypoxia = "Hypoxia", YAP_TAZ = "YAP/TAZ signature",
    IFN_gamma = "IFN-γ response", GEP18 = "T cell-inflamed GEP", Prolif = "Proliferation"
  )
  t1 <- tum[tum$group == "Tumor", ]
  vars <- vars[names(vars) %in% colnames(t1)]
  cm <- stats::cor(t1[, names(vars)], method = "spearman", use = "pairwise.complete.obs")
  rows <- c("MT_core", "MT_cortex", "MT_BAR", "MT_channel", "MT_ECM", "MTI")
  rows <- intersect(rows, rownames(cm))
  hm <- data.frame(x = rep(colnames(cm), each = length(rows)), y = rep(rows, ncol(cm)),
                   r = as.vector(cm[rows, ]))
  hm$x <- factor(vars[hm$x], levels = vars); hm$y <- factor(hm$y, levels = rev(rows))
  for (v in intersect(c("Prolif", "EMT", "TGF_beta", "Angiogenesis", "Hypoxia", "IFN_gamma", "GEP18"), colnames(cm)))
    for (s in c("MT_core", "MT_cortex", "MT_ECM", "MTI"))
      add_stat(4, "B", sprintf("Spearman rho %s vs %s (TCGA tumours)", s, v), round(cm[s, v], 2), paste(nrow(t1), "tumours"))
  pB <- ggplot2::ggplot(hm, ggplot2::aes(x, y, fill = r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = (BASE - 1.8) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_fill_gradient2(low = C_LO, mid = "white", high = C_HI, limits = c(-1, 1), name = "Spearman ρ") +
    ggplot2::labs(x = NULL, y = NULL) + ggplot2::coord_fixed() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1), axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), legend.position = "right")

  # C：MT_core 与 EMT 及增殖评分关系
  rho_emt <- stats::cor.test(t1$MT_core, t1$EMT, method = "spearman", exact = FALSE)
  rho_pro <- stats::cor.test(t1$MT_core, t1$Prolif, method = "spearman", exact = FALSE)
  add_stat(4, "C", "Spearman rho MT_core vs EMT score", round(rho_emt$estimate, 2),
           paste0(exact_p(rho_emt$p.value), "; n = ", nrow(t1)))
  add_stat(4, "C", "Spearman rho MT_core vs proliferation score", round(rho_pro$estimate, 2),
           paste0(exact_p(rho_pro$p.value), "; n = ", nrow(t1)))
  pC <- ggplot2::ggplot(t1, ggplot2::aes(EMT, MT_core)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6, colour = "grey30") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = C_HI, linewidth = 0.5, se = FALSE) +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.3, family = FONT, size = (BASE - 0.5) / ggplot2::.pt,
                      label = sprintf("EMT: ρ = %.2f, %s\nProlif: ρ = %.2f, %s",
                                      rho_emt$estimate, lab_p(rho_emt$p.value), rho_pro$estimate, lab_p(rho_pro$p.value))) +
    ggplot2::labs(x = "EMT score (ssGSEA)", y = "MT_core score", title = "TCGA-KIRC tumours")

  # D：ccRCC 核心体细胞突变（VHL、PBRM1、SETD2、BAP1）
  mt <- utils::read.csv(file.path(DIR_RES, "04_bulk_score", "MT_score_by_mutation.csv"))
  genes_m <- intersect(c("VHL", "PBRM1", "SETD2", "BAP1"), mt$gene)
  ml <- do.call(rbind, lapply(genes_m, function(g) data.frame(gene = g, status = t1[[paste0("mut_", g)]], MT_core = t1$MT_core)))
  ml <- ml[!is.na(ml$status), ]
  ml$status <- factor(ml$status, levels = c("WT", "Mut")); ml$gene <- factor(ml$gene, levels = genes_m)
  labd <- data.frame(gene = factor(genes_m, levels = genes_m),
                     lab = lab_q(mt$FDR_MT_core[match(genes_m, mt$gene)]))
  for (g in genes_m) add_stat(4, "D", paste("MT_core mutant vs wild-type, Wilcoxon FDR:", g),
                              signif(mt$FDR_MT_core[mt$gene == g], 3),
                              sprintf("%d mutant / %d tested", mt$n_mut[mt$gene == g], sum(!is.na(t1[[paste0("mut_", g)]]))))
  pD <- ggplot2::ggplot(ml, ggplot2::aes(status, MT_core, colour = status)) +
    ggplot2::geom_jitter(width = 0.2, size = 0.35, alpha = 0.6) +
    ggplot2::geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, linewidth = 0.35, colour = "black") +
    ggplot2::facet_wrap(~ gene, nrow = 1) +
    ggplot2::geom_text(data = labd, ggplot2::aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, inherit.aes = FALSE,
                       size = (BASE - 1) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::scale_colour_manual(values = c(WT = C_LO, Mut = C_HI), guide = "none") +
    ggplot2::labs(x = NULL, y = "MT_core score")

  # E：GSE53757 72 对配对肿瘤与正常肾组织验证
  g53 <- readRDS(file.path(DIR_PRO, "GSE53757.rds"))$clin
  g53 <- merge(g53, sc_all$GSE53757[, c("ID", "MT_core", "MTI")], by = "ID")
  g53$group <- factor(g53$group, levels = c("Normal", "Tumor"), labels = c("Normal", "Tumour"))
  p_g53 <- stats::wilcox.test(MT_core ~ group, data = g53)$p.value
  add_stat(4, "E", "GSE53757 tumour vs normal MT_core Wilcoxon P", signif(p_g53, 3),
           sprintf("n = %d (%d pairs)", nrow(g53), nrow(g53) / 2))
  pE <- ggplot2::ggplot(g53, ggplot2::aes(group, MT_core, colour = group)) +
    ggplot2::geom_jitter(width = 0.18, size = 0.45, alpha = 0.65) +
    ggplot2::geom_boxplot(fill = NA, outlier.shape = NA, width = 0.45, linewidth = 0.35, colour = "black") +
    ggplot2::annotate("text", x = 1.5, y = Inf, vjust = 1.3, family = FONT, size = (BASE - 0.5) / ggplot2::.pt,
                      label = lab_p(p_g53)) +
    ggplot2::scale_colour_manual(values = c(Normal = C_LO, Tumour = C_HI), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
    ggplot2::labs(x = NULL, y = "MT_core score", title = "GSE53757 (n = 144)")

  # F：TCGA-KIRC 中各 MT 子模块评分的 Cox（每 SD）
  cx <- utils::read.csv(file.path(DIR_RES, "04_bulk_score", "TCGA_Cox_MT_scores.csv"))
  cx <- cx[cx$endpoint == "OS", ]
  d_os <- t1[!is.na(t1$OS.time) & t1$OS.time > 0 & !is.na(t1$OS), ]
  fr <- do.call(rbind, lapply(cx$score, function(s) {
    u <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ scale(d_os[[s]]), data = d_os))
    rbind(
      data.frame(score = s, model = "Unadjusted", HR = u$HR, lower = u$lower, upper = u$upper, p = u$p),
      data.frame(score = s, model = "Adjusted", HR = cx$HR_adj[cx$score == s], lower = cx$lower_adj[cx$score == s],
                 upper = cx$upper_adj[cx$score == s], p = cx$p_adj[cx$score == s])
    )
  }))
  for (i in seq_len(nrow(fr))) add_stat(4, "F", sprintf("TCGA OS HR per SD, %s, %s", fr$score[i], fr$model[i]),
                                        sprintf("%.2f (%.2f-%.2f)", fr$HR[i], fr$lower[i], fr$upper[i]), exact_p(fr$p[i]))
  fr$label <- fr$score
  fr$model <- factor(fr$model, levels = c("Adjusted", "Unadjusted"))
  pF <- forest_gg(fr, group = "model") +
    ggplot2::geom_text(ggplot2::aes(x = max(fr$upper) * 1.15, label = fmt_p(p)), hjust = 0, show.legend = FALSE,
                       size = (BASE - 1) / ggplot2::.pt, family = FONT, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_colour_manual(values = c(Unadjusted = "grey35", Adjusted = C_HI), name = NULL,
                                 breaks = c("Unadjusted", "Adjusted"),
                                 labels = c("Unadjusted", "Adjusted for age, stage, proliferation")) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::theme(plot.margin = ggplot2::margin(5, 22, 5, 5)) +
    ggplot2::labs(title = "TCGA-KIRC, overall survival") +
    ggplot2::theme(legend.position = "bottom")

  # G：各队列中 MT_core 的 Cox
  ex <- utils::read.csv(file.path(DIR_RES, "04_bulk_score", "External_Cox_MT_core.csv"))
  ex <- rbind(
    data.frame(cohort = "TCGA-KIRC", n = nrow(d_os), events = sum(d_os$OS),
               fr[fr$score == "MT_core" & fr$model == "Unadjusted", c("HR", "lower", "upper", "p")]),
    data.frame(ex[, c("cohort", "n", "events")], HR = ex$HR_perSD, ex[, c("lower", "upper", "p")])
  )
  ex$label <- sprintf("%s (%d/%d)", ex$cohort, ex$events, ex$n)
  ex$label[ex$cohort == "GSE22541"] <- paste0(ex$label[ex$cohort == "GSE22541"], " DFS")
  for (i in seq_len(nrow(ex))) add_stat(4, "G", paste("MT_core HR per SD,", ex$cohort[i]),
                                        sprintf("%.2f (%.2f-%.2f)", ex$HR[i], ex$lower[i], ex$upper[i]), exact_p(ex$p[i]))
  pG <- forest_gg(ex) + ggplot2::geom_point(size = 1.3, colour = C_HI) +
    ggplot2::geom_text(ggplot2::aes(x = max(ex$upper) * 1.15, label = fmt_p(p)), hjust = 0, size = (BASE - 1) / ggplot2::.pt,
                       family = FONT) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::theme(plot.margin = ggplot2::margin(5, 22, 5, 5)) +
    ggplot2::labs(title = "MT_core score across cohorts (events/patients)")

  top <- (pA | pC) + patchwork::plot_layout(widths = c(1.15, 1))
  mid <- pB
  low <- (pD | pE) + patchwork::plot_layout(widths = c(2.2, 1))
  bot <- (pF | pG) + patchwork::plot_layout(widths = c(1, 1.1))
  p <- (top / mid / low / bot) + patchwork::plot_layout(heights = c(1, 0.65, 0.9, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
  save_fig(p, 4, 180, 235)
}

# =============================================================================
# 共用：各队列的 MTRS 与临床数据
# =============================================================================
COH_LAB  <- c(TCGA_KIRC = "TCGA-KIRC", GSE29609 = "GSE29609", GSE22541 = "GSE22541",
              IMvigor210 = "IMvigor210", GSE176307 = "GSE176307", GSE67501 = "GSE67501",
              GSE53757 = "GSE53757", GSE73731 = "GSE73731")
COH_ROLE <- c(TCGA_KIRC = "training", GSE29609 = "tuning, OS", GSE22541 = "tuning, DFS",
              IMvigor210 = "test, ICI", GSE176307 = "test, ICI", GSE67501 = "RCC Nivolumab")
COH_COLS <- c(`TCGA-KIRC` = "grey30", GSE29609 = "#2E6DA4", GSE22541 = "#C0392B",
              IMvigor210 = "#8E44AD", GSE176307 = "#E67E22")
SURV_COH <- c("TCGA_KIRC", "GSE29609", "GSE22541", "IMvigor210", "GSE176307")

coh_clin <- function(n) {
  if (n == "TCGA_KIRC") readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds")) else readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$clin
}
rlrs_all <- function() readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
rlrs_cohort <- function(n) {
  r <- rlrs_all(); r <- r[r$cohort == n, c("ID", "RLRS", "RLRS_z", "RLRS_group")]
  merge(r, coh_clin(n), by = "ID")
}
txt <- function(...) ggplot2::annotate("text", size = (BASE - 0.5) / ggplot2::.pt, family = FONT, ...)

# =============================================================================
# Figure 5：候选基因与 MTRS 的构建、多队列验证
# =============================================================================
fig5 <- function() {
  GS <- readRDS(file.path(DIR_GS, "MTRG_genesets.rds"))
  cands <- readRDS(file.path(DIR_PRO, "candidates.rds"))
  model <- readRDS(file.path(DIR_PRO, "RLRS_model.rds"))
  d6 <- file.path(DIR_RES, "06_ml_signature")
  sc_deg <- readLines(file.path(DIR_RES, "03_scRNA", "sc_DEG_genes.txt"))
  if (!length(sc_deg)) {
    pbd <- utils::read.csv(file.path(DIR_RES, "03_scRNA", "pseudobulk_DE_High_vs_Low.csv"))
    sc_deg <- pbd$gene[pbd$PValue < 0.05 & abs(pbd$logFC) > 0.3]
  }
  uc <- utils::read.csv(file.path(DIR_RES, "05_candidates", "Table_S3_candidates_unicox.csv"))

  pool <- union(GS$MT_all, cands$key_module_genes)
  n_steps <- c(length(pool), length(intersect(pool, sc_deg)), length(cands$candidates), sum(uc$all_cohorts),
               sum(uc$all_cohorts & uc$p < PARAM$unicox_p), length(model$genes))
  lab_steps <- c("MT_all or WGCNA key-module genes", "+ DEG, MT-high vs -low cells", "+ DEG, tumour vs normal (TCGA)",
                 "+ measured in all validation cohorts", "+ univariate Cox P < 0.05 (TCGA)", "Genes in the selected MTRS model")
  for (i in seq_along(n_steps)) add_stat(5, "A", paste("candidate selection:", lab_steps[i]), n_steps[i])
  st <- data.frame(step = factor(lab_steps, levels = rev(lab_steps)), n = n_steps)
  pA <- ggplot2::ggplot(st, ggplot2::aes(y = step)) +
    ggplot2::geom_segment(ggplot2::aes(x = 1, xend = n, yend = step), linewidth = 3.2, colour = "#D98880") +
    ggplot2::geom_text(ggplot2::aes(x = n, label = format(n, big.mark = ",")), hjust = -0.2,
                       size = (BASE - 0.5) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0, 0.22))) +
    ggplot2::labs(x = "Genes (log scale)", y = NULL, title = "Candidate genes")

  ug <- uc[uc$gene %in% model$genes, ]
  ug <- ug[order(-ug$HR), ]
  ug$dir <- ifelse(ug$HR > 1, "Higher risk", "Lower risk")
  add_stat(5, "B", "model genes with HR > 1 in TCGA (univariate, per SD)", sum(ug$HR > 1), paste("of", nrow(ug)))
  pB <- forest_gg(ug, label = "gene", group = "dir", xlab = "HR per SD (95% CI)",
                  cols = c(`Higher risk` = C_HI, `Lower risk` = C_LO)) +
    ggplot2::labs(title = "Model genes, TCGA-KIRC OS") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"), legend.position = "none")

  ci <- utils::read.csv(file.path(d6, "Table_S4_Cindex_101_models.csv"), check.names = FALSE)
  ci$sel <- ci$Model == model$model
  prof <- apply(round(ci[, SURV_COH], 6), 1, paste, collapse = "_")
  rk <- rank(-ci$mean_valid)[ci$sel]
  rt <- stats::cor.test(ci$mean_tune, ci$mean_test, method = "spearman", exact = FALSE)
  add_stat(5, "C", "algorithm combinations evaluated", nrow(ci), sprintf("%d distinct C-index profiles", length(unique(prof))))
  add_stat(5, "C", "selected model", model$model,
           sprintf("mean C-index validation %.3f (tuning %.3f, test %.3f); %d genes",
                   ci$mean_valid[ci$sel], ci$mean_tune[ci$sel], ci$mean_test[ci$sel], ci$n_genes[ci$sel]))
  add_stat(5, "C", "rank of the selected model by mean validation C-index", rk,
           sprintf("of %d; mean validation C-index range %.3f-%.3f", nrow(ci), min(ci$mean_valid), max(ci$mean_valid)))
  pC <- ggplot2::ggplot(ci, ggplot2::aes(mean_tune, mean_test)) +
    ggplot2::geom_point(data = ci[!ci$sel, ], colour = "grey60", size = 0.8, alpha = 0.8) +
    ggplot2::geom_point(data = ci[ci$sel, ], colour = C_HI, size = 1.8) +
    txt(x = ci$mean_tune[ci$sel], y = ci$mean_test[ci$sel], label = "Selected", hjust = 1.1, vjust = 1.4, colour = C_HI) +
    ggplot2::labs(x = "Mean C-index, ccRCC validation cohorts", y = "Mean C-index, ICI test cohorts",
                  title = sprintf("%d algorithm combinations", nrow(ci)))

  cx <- utils::read.csv(file.path(d6, "RLRS_Cox_by_cohort.csv"))
  cd <- data.frame(cohort = cx$cohort, C = cx$Cindex, lo = cx$Cindex - 1.96 * cx$Cindex_se, hi = cx$Cindex + 1.96 * cx$Cindex_se)
  for (i in seq_len(nrow(cd))) add_stat(5, "D", paste("C-index (fixed direction):", cd$cohort[i]),
                                        sprintf("%.3f (%.3f-%.3f)", cd$C[i], cd$lo[i], cd$hi[i]), COH_ROLE[cd$cohort[i]])
  cd$lab <- factor(COH_LAB[cd$cohort], levels = rev(COH_LAB[cd$cohort]))
  cd$role <- factor(sub(",.*", "", COH_ROLE[cd$cohort]), levels = c("training", "tuning", "test"))
  pD <- ggplot2::ggplot(cd, ggplot2::aes(C, lab, colour = role)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.4, orientation = "y") +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(values = c(training = "grey30", tuning = C_LO, test = C_HI), name = NULL) +
    ggplot2::labs(x = "C-index (95% CI)", y = NULL, title = "Selected MTRS model") +
    ggplot2::theme(legend.position = "bottom")

  au <- utils::read.csv(file.path(d6, "RLRS_timeROC_AUC.csv"))
  for (i in seq_len(nrow(au))) add_stat(5, "J", sprintf("time-dependent AUC, %d year(s): %s", au$year[i], au$cohort[i]), sprintf("%.3f", au$AUC[i]))
  au$lab <- factor(COH_LAB[au$cohort], levels = COH_LAB[SURV_COH])
  pE <- ggplot2::ggplot(au, ggplot2::aes(factor(year), AUC, colour = lab, group = lab)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.4) + ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_colour_manual(values = COH_COLS, name = NULL) +
    ggplot2::coord_cartesian(ylim = c(0.45, 0.85)) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::labs(x = "Years", y = "Time-dependent AUC", title = "Discrimination over time") +
    ggplot2::theme(legend.position = "bottom", legend.key.size = ggplot2::unit(2.5, "mm"),
                   legend.text = ggplot2::element_text(size = BASE - 1.5), legend.margin = ggplot2::margin(0, 0, 0, 0))

  kms <- lapply(SURV_COH, function(n) {
    d <- rlrs_cohort(n); d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
    g <- km_gg(d$OS.time, d$OS, d$RLRS_group, title = sprintf("%s (%s)", COH_LAB[n], COH_ROLE[n]),
               ylab = if (n == "GSE22541") "Disease-free survival" else "Overall survival")
    s <- attr(g, "stats")
    add_stat(5, "E-I", paste("KM MTRS high vs low:", n), sprintf("HR %.2f (%.2f-%.2f)", s["HR"], s["lower"], s["upper"]),
             paste0("log-rank ", exact_p(s["p_logrank"]), "; n = ", nrow(d), ", events = ", sum(d$OS)))
    g
  })
  tagged <- Map(function(g, t) g + ggplot2::labs(tag = t), c(list(pA, pB, pC, pD), kms, list(pE)), LETTERS[1:10])
  top <- (tagged[[1]] | tagged[[2]]) / (patchwork::free(tagged[[3]]) | tagged[[4]])
  bottom <- patchwork::wrap_elements(full = patchwork::wrap_plots(tagged[5:10], ncol = 3))
  p <- patchwork::wrap_plots(list(patchwork::wrap_elements(full = top), bottom), ncol = 1, heights = c(1.05, 1))
  save_fig(p, 5, 180, 225)
}

# =============================================================================
# Figure 6：特异性：随机基因集置换、增殖校正、与 10 个已发表 ccRCC 模型对标
# =============================================================================
adv_stage <- function(n, c0) {
  if ("stage" %in% colnames(c0)) {
    s <- toupper(as.character(c0$stage))
    return(ifelse(grepl("III|IV|3|4", s), 1L, ifelse(grepl("I|II|1|2", s), 0L, NA_integer_)))
  }
  rep(NA_integer_, nrow(c0))
}
fig6 <- function() {
  d6b <- file.path(DIR_RES, "06b_specificity"); d6 <- file.path(DIR_RES, "06_ml_signature")
  nl <- readRDS(file.path(DIR_PRO, "random_signature_null.rds"))
  bm <- utils::read.csv(file.path(d6b, "random_signature_benchmark.csv"))
  obs <- bm$RLRS_best_mean[1]
  nA <- rowMeans(nl$nullA, na.rm = TRUE); nB <- rowMeans(nl$nullB, na.rm = TRUE)
  pv <- bm$p_RLRS_best
  nm <- c("Random expressed genes", "Random prognostic genes")
  dn <- rbind(data.frame(null = nm[1], C = nA), data.frame(null = nm[2], C = nB))
  lb <- data.frame(null = nm, p = pv)
  for (i in 1:2) add_stat(6, "A", paste("random-signature null:", nm[i]),
                          sprintf("mean %.3f, 95th percentile %.3f", bm$null_mean[i], bm$null_q95[i]),
                          sprintf("MTRS selected model %.3f, empirical P = %.3g; ridge %.3f, P = %.3g",
                                  bm$RLRS_best_mean[i], bm$p_RLRS_best[i], bm$RLRS_ridge_mean[i], bm$p_RLRS_ridge[i]))
  pA <- ggplot2::ggplot(dn, ggplot2::aes(C)) +
    ggplot2::geom_histogram(bins = 35, fill = "grey75", colour = NA) +
    ggplot2::geom_vline(xintercept = obs, colour = C_HI, linewidth = 0.5) +
    ggplot2::geom_text(data = lb, ggplot2::aes(x = obs, y = Inf, label = paste0("MTRS\nempirical ", lab_p(p))),
                       hjust = 1.05, vjust = 1.3, size = (BASE - 1) / ggplot2::.pt, family = FONT, lineheight = 0.9) +
    ggplot2::facet_wrap(~ null, ncol = 1) +
    ggplot2::labs(x = "Mean C-index, four validation cohorts", y = "Random 16-gene sets", title = "Random-signature benchmark")

  sc <- readRDS(file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds"))
  r <- rlrs_all()
  t1 <- merge(r[r$cohort == "TCGA_KIRC", ], sc$TCGA_KIRC[, c("ID", "Prolif")], by = "ID")
  rr <- stats::cor.test(t1$RLRS, t1$Prolif, method = "spearman", exact = FALSE)
  pB <- ggplot2::ggplot(t1, ggplot2::aes(Prolif, RLRS_z)) +
    ggplot2::geom_point(size = 0.45, colour = "grey45", alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = C_HI, linewidth = 0.5) +
    txt(x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, label = sprintf("ρ = %.2f", rr$estimate)) +
    ggplot2::labs(x = "Proliferation score", y = "MTRS (z-score)", title = "TCGA-KIRC")

  adj <- utils::read.csv(file.path(d6b, "Table_S6_proliferation_adjusted.csv"))
  for (i in seq_len(nrow(adj))) add_stat(6, "C", paste("Spearman rho MTRS vs proliferation / MT_core:", adj$cohort[i]),
                                         sprintf("%.2f / %.2f", adj$rho_RLRS_Prolif[i], adj$rho_RLRS_RLcore[i]))
  rl <- rbind(
    data.frame(cohort = adj$cohort, what = "Proliferation", rho = abs(adj$rho_RLRS_Prolif)),
    data.frame(cohort = adj$cohort, what = "|MT_core|", rho = abs(adj$rho_RLRS_RLcore))
  )
  rl$lab <- factor(COH_LAB[rl$cohort], levels = rev(COH_LAB[adj$cohort]))
  pC <- ggplot2::ggplot(rl, ggplot2::aes(rho, lab, fill = what)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
    ggplot2::scale_fill_manual(values = c(Proliferation = "grey60", `|MT_core|` = C_HI), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(x = "|Spearman ρ| with MTRS", y = NULL, title = "Correlation with MTRS") +
    ggplot2::theme(legend.position = "bottom")

  cc <- rbind(
    data.frame(cohort = adj$cohort, model = "Proliferation", C = adj$C_Prolif),
    data.frame(cohort = adj$cohort, model = "MTRS", C = adj$C_RLRS),
    data.frame(cohort = adj$cohort, model = "Proliferation + MTRS", C = adj$C_Prolif_RLRS)
  )
  for (i in seq_len(nrow(adj))) add_stat(6, "D", paste("Cox C-index proliferation / MTRS / both:", adj$cohort[i]),
                                         sprintf("%.3f / %.3f / %.3f", adj$C_Prolif[i], adj$C_RLRS[i], adj$C_Prolif_RLRS[i]),
                                         paste("LRT for adding MTRS", exact_p(adj$LRT_p_add_RLRS[i])))
  cc$lab <- factor(COH_LAB[cc$cohort], levels = rev(COH_LAB[adj$cohort]))
  cc$model <- factor(cc$model, levels = c("Proliferation", "MTRS", "Proliferation + MTRS"))
  pD <- ggplot2::ggplot(cc, ggplot2::aes(C, lab, colour = model, shape = model)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_point(size = 1.5, position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::scale_colour_manual(values = c(Proliferation = "grey50", MTRS = C_HI, `Proliferation + MTRS` = "#641E16"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(Proliferation = 16, MTRS = 17, `Proliferation + MTRS` = 15), name = NULL) +
    ggplot2::labs(x = "C-index", y = NULL, title = "Added value over proliferation") +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")

  fr <- list()
  for (n in adj$cohort) {
    c0 <- coh_clin(n)
    d <- merge(data.frame(ID = c0$ID, OS.time = c0$OS.time, OS = c0$OS, adv = adv_stage(n, c0)),
               r[r$cohort == n, c("ID", "RLRS")], by = "ID")
    d <- merge(d, sc[[n]][, c("ID", "Prolif")], by = "ID")
    d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
    d$RLRS_z <- as.numeric(scale(d$RLRS)); d$Prolif_z <- as.numeric(scale(d$Prolif))
    one <- function(f, what) {
      m <- survival::coxph(f, data = d); s <- summary(m)
      data.frame(cohort = n, model = what, HR = s$conf.int["RLRS_z", 1], lower = s$conf.int["RLRS_z", 3],
                 upper = s$conf.int["RLRS_z", 4], p = s$coefficients["RLRS_z", 5])
    }
    o <- rbind(
      one(survival::Surv(OS.time, OS) ~ RLRS_z, "Unadjusted"),
      one(survival::Surv(OS.time, OS) ~ RLRS_z + Prolif_z, "Adjusted for proliferation")
    )
    if (sum(!is.na(d$adv)) > 20 && length(unique(stats::na.omit(d$adv))) == 2)
      o <- rbind(o, one(survival::Surv(OS.time, OS) ~ RLRS_z + Prolif_z + adv, "Adjusted for proliferation and stage"))
    fr[[n]] <- o
  }
  fr <- do.call(rbind, fr)
  for (i in seq_len(nrow(fr))) add_stat(6, "E", paste("MTRS HR per SD,", fr$model[i], ":", fr$cohort[i]),
                                        sprintf("%.2f (%.2f-%.2f)", fr$HR[i], fr$lower[i], fr$upper[i]), exact_p(fr$p[i]))
  fr$lab <- COH_LAB[fr$cohort]
  fr$model <- factor(fr$model, levels = c("Unadjusted", "Adjusted for proliferation", "Adjusted for proliferation and stage"))
  pE <- forest_gg(fr, label = "lab", group = "model", xlab = "HR per SD of MTRS (95% CI)",
                  cols = c(Unadjusted = "grey55", `Adjusted for proliferation` = C_HI, `Adjusted for proliferation and stage` = "#641E16")) +
    ggplot2::labs(title = "Survival models") +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")

  pub <- utils::read.csv(file.path(d6, "Table_S5_published_models_Cindex.csv"))
  rk <- utils::read.csv(file.path(d6, "RLRS_rank_among_published.csv"))
  ext <- c("TCGA_KIRC", "GSE29609", "GSE22541", "IMvigor210", "GSE176307")
  pub <- pub[pub$cohort %in% ext & !is.na(pub$Cindex), ]
  rk <- rk[match(ext, rk$cohort), ]
  ylab_rk <- stats::setNames(sprintf("%s (%d/%d)", COH_LAB[ext], rk$RLRS_rank, rk$n_models), ext)
  pub$lab <- factor(ylab_rk[pub$cohort], levels = rev(ylab_rk))
  add_stat(6, "F", "published signatures evaluated", length(unique(pub$model[pub$model != "MTRS (this study)"])),
           "across five cohorts")
  for (i in seq_len(nrow(rk))) add_stat(6, "F", paste("MTRS rank among models by C-index:", rk$cohort[i]),
                                        sprintf("%d of %d", rk$RLRS_rank[i], rk$n_models[i]),
                                        sprintf("MTRS %.3f; published median %.3f, max %.3f", rk$RLRS_Cindex[i], rk$published_median[i], rk$published_max[i]))
  pF <- ggplot2::ggplot(pub[pub$model != "MTRS (this study)", ], ggplot2::aes(Cindex, lab)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3, colour = "grey45") +
    ggplot2::geom_point(position = ggplot2::position_jitter(height = 0.15, seed = 1), size = 0.4, colour = "grey70") +
    ggplot2::geom_point(data = pub[pub$model == "MTRS (this study)", ], colour = C_HI, size = 2.2, shape = 18) +
    ggplot2::labs(x = "C-index", y = NULL, title = "MTRS (red) vs published signatures")

  row1 <- (pA | pB | pC) + patchwork::plot_layout(widths = c(1.05, 0.9, 1))
  row2 <- (pD | pE | pF) + patchwork::plot_layout(widths = c(0.9, 1.05, 1.1))
  p <- (row1 / row2) + patchwork::plot_layout(heights = c(1, 1.1)) + patchwork::plot_annotation(tag_levels = "A")
  save_fig(p, 6, 180, 165)
}

# =============================================================================
# Figure 7：临床意义（多因素 Cox、列线图、校准、时间 ROC、DCA）
# =============================================================================
fig7 <- function() {
  need(c("rms", "gridGraphics"))
  d7 <- file.path(DIR_RES, "07_clinical_nomogram")
  um <- utils::read.csv(file.path(d7, "TCGA_uni_multi_Cox.csv"))
  um$model <- factor(ifelse(um$model == "Univariate", "Univariable", "Multivariable"), levels = c("Univariable", "Multivariable"))
  um$variable <- sub("^Gender", "Sex", um$variable)
  for (i in seq_len(nrow(um))) add_stat(7, "A", paste(um$model[i], "Cox, TCGA OS:", um$variable[i]),
                                        sprintf("%.2f (%.2f-%.2f)", um$HR[i], um$lower[i], um$upper[i]), exact_p(um$p[i]))
  pA <- forest_gg(um, label = "variable", group = "model", xlab = "HR (95% CI)",
                  cols = c(Univariable = "grey55", Multivariable = C_HI)) +
    ggplot2::labs(title = "TCGA-KIRC, overall survival") + ggplot2::theme(legend.position = "bottom")

  ex <- utils::read.csv(file.path(d7, "External_multivariate_Cox.csv"))
  ex <- ex[ex$term == "RLRS_z", ]
  cv_map <- c(
    `stage_bin+grade_bin` = "stage, grade",
    `gender+stage_bin+grade_bin` = "sex, stage, grade",
    `gender` = "sex",
    `age10+gender` = "age, sex"
  )
  tc <- um[um$model == "Multivariable" & grepl("^MTRS|^RLRS", um$variable), ]
  eb <- rbind(
    data.frame(lab = "TCGA-KIRC\n(age, sex, stage, grade)", HR = tc$HR[1], lower = tc$lower[1], upper = tc$upper[1], p = tc$p[1], LRT = NA),
    data.frame(lab = paste0(COH_LAB[ex$cohort], "\n(", ifelse(ex$covariates %in% names(cv_map), cv_map[ex$covariates], ex$covariates), ")"),
               HR = ex$HR, lower = ex$lower, upper = ex$upper, p = ex$p, LRT = ex$LRT_p_add_RLRS)
  )
  for (i in seq_len(nrow(ex))) add_stat(7, "B", paste("MTRS HR per SD adjusted for clinical factors:", ex$cohort[i]),
                                        sprintf("%.2f (%.2f-%.2f)", ex$HR[i], ex$lower[i], ex$upper[i]),
                                        sprintf("%s; LRT %s; C-index clinical %.3f vs clinical + MTRS %.3f", exact_p(ex$p[i]),
                                                exact_p(ex$LRT_p_add_RLRS[i]), ex$C_clinical[i], ex$C_clinical_RLRS[i]))
  eb$plab <- lab_p(eb$p)
  pB <- forest_gg(eb, label = "lab", xlab = "HR per SD of MTRS (95% CI)") +
    ggplot2::geom_text(ggplot2::aes(x = upper, label = plab), hjust = -0.15, size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.05, 0.45))) +
    ggplot2::labs(title = "Multivariable models")

  nom <- readRDS(file.path(DIR_PRO, "nomogram_object.rds"))
  nm_idx <- match(c("age", "stage_bin", "grade_bin", "RLRS_z", "RLRS"), names(nom))
  nm_new <- c("Age (years)", "Stage", "Grade", "MTRS (z)", "MTRS")
  for (ii in seq_along(nm_idx)) if (!is.na(nm_idx[ii])) names(nom)[nm_idx[ii]] <- nm_new[ii]
  pC <- patchwork::wrap_elements(full = ~ {
    graphics::par(mar = c(0.5, 0.5, 0.5, 0.5), family = FONT, cex = 0.56)
    plot(nom, xfrac = 0.28, lmgp = 0.18, cex.axis = 0.85, cex.var = 1, col.grid = grDevices::gray(c(0.85, 0.95)))
  })

  cal <- utils::read.csv(file.path(d7, "calibration_bootstrap.csv"))
  cal$yr <- factor(paste0(cal$year, "-year"), levels = c("1-year", "3-year", "5-year"))
  yr_cols <- c(`1-year` = "#E59866", `3-year` = C_HI, `5-year` = "#641E16")
  pD <- ggplot2::ggplot(cal, ggplot2::aes(mean.predicted, colour = yr)) +
    ggplot2::geom_abline(linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_line(ggplot2::aes(y = KM.corrected), linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(y = KM.corrected), size = 1) +
    ggplot2::scale_colour_manual(values = yr_cols, name = NULL) +
    ggplot2::coord_equal(xlim = c(0.2, 1), ylim = c(0.2, 1)) +
    ggplot2::labs(x = "Predicted survival (nomogram)", y = "Observed survival\n(bias-corrected)", title = "Calibration") +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
                   legend.background = ggplot2::element_blank())

  cn <- utils::read.csv(file.path(d7, "Cindex_nomogram_RLRS_stage.csv"))
  for (i in seq_len(nrow(cn))) add_stat(7, "E", paste("C-index, TCGA:", cn$model[i]), sprintf("%.3f", cn$Cindex[i]))
  ta <- utils::read.csv(file.path(d7, "timeROC_nomogram_RLRS_stage.csv"))
  tl <- stats::reshape(ta, direction = "long", varying = c("Nomogram", "RLRS", "Stage"), v.names = "AUC",
                       timevar = "model", times = c("Nomogram", "MTRS", "Stage"), idvar = "year")
  for (i in seq_len(nrow(tl))) add_stat(7, "E", sprintf("time-dependent AUC, %d year(s): %s", tl$year[i], tl$model[i]), sprintf("%.3f", tl$AUC[i]))
  m_cols <- c(Nomogram = "#641E16", MTRS = C_HI, Stage = "grey50")
  pE <- ggplot2::ggplot(tl, ggplot2::aes(factor(year), AUC, colour = model, group = model)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.4) + ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_colour_manual(values = m_cols, name = NULL) + ggplot2::coord_cartesian(ylim = c(0.45, 0.85)) +
    ggplot2::labs(x = "Years", y = "Time-dependent AUC", title = "TCGA-KIRC") +
    ggplot2::theme(legend.position = "bottom")

  dc <- utils::read.csv(file.path(d7, "DCA_net_benefit.csv"))
  dc <- dc[dc$year == 3 & dc$threshold <= 0.6, ]
  dc$label[dc$label == "RLRS"] <- "MTRS"
  dc$label <- factor(dc$label, levels = c("Treat All", "Treat None", "Nomogram", "MTRS", "Stage"))
  pF <- ggplot2::ggplot(dc, ggplot2::aes(threshold, net_benefit, colour = label, linetype = label)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = c(`Treat All` = "grey65", `Treat None` = "black", m_cols), name = NULL) +
    ggplot2::scale_linetype_manual(values = c(`Treat All` = 2, `Treat None` = 1, Nomogram = 1, MTRS = 1, Stage = 1), name = NULL) +
    ggplot2::coord_cartesian(ylim = c(-0.05, max(dc$net_benefit[dc$label != "Treat All"], na.rm = TRUE) * 1.15)) +
    ggplot2::labs(x = "Threshold probability (death within 3 years)", y = "Net benefit", title = "Decision curves, 3 years") +
    ggplot2::theme(legend.position = "bottom")

  p <- ((pA | pB) + patchwork::plot_layout(widths = c(1, 1))) /
    ((pC | pD) + patchwork::plot_layout(widths = c(1.45, 1))) /
    ((patchwork::free(pE) | pF) + patchwork::plot_layout(widths = c(1, 1.3))) +
    patchwork::plot_layout(heights = c(1, 1.15, 1)) + patchwork::plot_annotation(tag_levels = "A")
  save_fig(p, 7, 180, 205)
}

# =============================================================================
# Figure 8：免疫微环境与治疗反应（IOBR、检查点、TIDE、ICI、分期/分级递进、oncoPredict）
# =============================================================================
TME_MAP <- data.frame(
  variable = c(
    "T_cells_CD8_CIBERSORT", "CD8_Tcells_EPIC", "CD8_T_cells_MCPcounter", "T_cells_CD8_quantiseq", "T_cell_CD8_TIMER", "CD8+_T-cells_xCell",
    "T_cells_CD4_memory_activated_CIBERSORT", "CD4_Tcells_EPIC", "T_cells_CD4_quantiseq", "T_cell_CD4_TIMER", "CD4+_T-cells_xCell",
    "T_cells_regulatory_(Tregs)_CIBERSORT", "Tregs_quantiseq", "Tregs_xCell",
    "NK_cells_activated_CIBERSORT", "NKcells_EPIC", "NK_cells_MCPcounter", "NK_cells_quantiseq", "NK_cells_xCell",
    "B_cells_naive_CIBERSORT", "Bcells_EPIC", "B_lineage_MCPcounter", "B_cells_quantiseq", "B_cell_TIMER", "B-cells_xCell",
    "Macrophages_M1_CIBERSORT", "Macrophages_M1_quantiseq", "Macrophages_M1_xCell",
    "Macrophages_M2_CIBERSORT", "Macrophages_M2_quantiseq", "Macrophages_M2_xCell",
    "Macrophages_EPIC", "Monocytic_lineage_MCPcounter", "Macrophage_TIMER", "Macrophages_xCell",
    "Dendritic_cells_activated_CIBERSORT", "Myeloid_dendritic_cells_MCPcounter", "Dendritic_cells_quantiseq", "DC_TIMER", "DC_xCell",
    "Neutrophils_CIBERSORT", "Neutrophils_MCPcounter", "Neutrophils_quantiseq", "Neutrophil_TIMER", "Neutrophils_xCell",
    "CAFs_EPIC", "Fibroblasts_MCPcounter", "Fibroblasts_xCell",
    "Endothelial_EPIC", "Endothelial_cells_MCPcounter", "Endothelial_cells_xCell",
    "ImmuneScore_estimate", "StromalScore_estimate", "TumorPurity_estimate"
  ),
  cell = c(
    rep("CD8 T cell", 6), rep("CD4 T cell", 5), rep("Treg", 3), rep("NK cell", 5), rep("B cell", 6),
    rep("M1 macrophage", 3), rep("M2 macrophage", 3), rep("Macrophage/monocyte", 4), rep("Dendritic cell", 5),
    rep("Neutrophil", 5), rep("Fibroblast", 3), rep("Endothelial cell", 3), "Immune score", "Stromal score", "Tumour purity"
  ),
  method = c(
    "CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell", "CIBERSORT", "EPIC", "quanTIseq", "TIMER", "xCell",
    "CIBERSORT", "quanTIseq", "xCell", "CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "xCell",
    "CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell", "CIBERSORT", "quanTIseq", "xCell",
    "CIBERSORT", "quanTIseq", "xCell", "EPIC", "MCPcounter", "TIMER", "xCell",
    "CIBERSORT", "MCPcounter", "quanTIseq", "TIMER", "xCell", "CIBERSORT", "MCPcounter", "quanTIseq", "TIMER", "xCell",
    "EPIC", "MCPcounter", "xCell", "EPIC", "MCPcounter", "xCell", "ESTIMATE", "ESTIMATE", "ESTIMATE"
  ),
  check.names = FALSE
)
TME_NONCELL <- c("P-value_CIBERSORT", "Correlation_CIBERSORT", "RMSE_CIBERSORT", "Other_quantiseq", "otherCells_EPIC")

fig8 <- function() {
  need("pROC")
  d8 <- file.path(DIR_RES, "08_tme_therapy")
  tme <- utils::read.csv(file.path(d8, "TME_High_vs_Low_all_methods.csv"), check.names = FALSE)
  cellv <- tme[!tme$variable %in% TME_NONCELL & !grepl("_estimate$", tme$variable) &
                 !tme$variable %in% c("ImmuneScore_xCell", "StromaScore_xCell", "MicroenvironmentScore_xCell"), ]
  add_stat(8, "A", "cell-type estimates (7 methods) with FDR < 0.05: higher / lower in MTRS-high",
           sprintf("%d / %d", sum(cellv$FDR < 0.05 & cellv$median_High > cellv$median_Low),
                   sum(cellv$FDR < 0.05 & cellv$median_High < cellv$median_Low)), paste("of", nrow(cellv), "estimates"))
  ta <- merge(TME_MAP, tme, by = "variable")
  for (i in seq_len(nrow(ta))) if (ta$method[i] == "ESTIMATE" || ta$FDR[i] < 0.05)
    add_stat(8, "A", paste("rho with MTRS:", ta$cell[i], "/", ta$method[i]), sprintf("%.2f", ta$rho_RLRS[i]), lab_q(ta$FDR[i]))
  ta$cell <- factor(ta$cell, levels = rev(unique(TME_MAP$cell)))
  ta$method <- factor(ta$method, levels = c("CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell", "ESTIMATE"))
  ta$sig <- ifelse(ta$FDR < 0.05, "FDR < 0.05", "NS")
  pA <- ggplot2::ggplot(ta, ggplot2::aes(method, cell)) +
    ggplot2::geom_point(ggplot2::aes(size = pmin(-log10(FDR), 10), fill = rho_RLRS, colour = sig), shape = 21, stroke = 0.35) +
    ggplot2::scale_fill_gradient2(low = C_LO, mid = "white", high = C_HI, limits = c(-0.7, 0.7), oob = scales::squish,
                                  name = "Spearman ρ\nwith MTRS") +
    ggplot2::scale_size_area(max_size = 3.2, name = "-log10 FDR", breaks = c(2, 5, 10), labels = c("2", "5", "≥10")) +
    ggplot2::scale_colour_manual(values = c(`FDR < 0.05` = "black", NS = "grey80"), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "Immune and stromal cells, TCGA-KIRC") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.box = "vertical")

  ck <- utils::read.csv(file.path(d8, "Checkpoints_High_vs_Low.csv"))
  for (i in seq_len(nrow(ck))) add_stat(8, "B", paste("checkpoint gene rho with MTRS:", ck$variable[i]),
                                        sprintf("%.2f", ck$rho_RLRS[i]), lab_q(ck$FDR[i]))
  ck <- ck[order(ck$rho_RLRS), ]; ck$variable <- factor(ck$variable, levels = ck$variable)
  ck$sig <- ifelse(ck$FDR < 0.05, ifelse(ck$rho_RLRS > 0, "Positive", "Negative"), "NS")
  pB <- ggplot2::ggplot(ck, ggplot2::aes(rho_RLRS, variable, fill = sig)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c(Positive = C_HI, Negative = C_LO, NS = "grey75"), guide = "none") +
    ggplot2::labs(x = "Spearman ρ with MTRS", y = NULL, title = "Immune checkpoints") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))

  r <- rlrs_all(); rt <- r[r$cohort == "TCGA_KIRC", ]
  td <- utils::read.csv(file.path(d8, "TIDE_output.csv"), check.names = FALSE)
  td <- merge(rt, td, by.x = "ID", by.y = "Patient")
  tl <- rbind(
    data.frame(score = "TIDE", v = td$TIDE, g = td$RLRS_group),
    data.frame(score = "Dysfunction", v = td$Dysfunction, g = td$RLRS_group),
    data.frame(score = "Exclusion", v = td$Exclusion, g = td$RLRS_group)
  )
  tl$score <- factor(tl$score, levels = c("TIDE", "Dysfunction", "Exclusion")); tl$g <- factor(tl$g, levels = c("Low", "High"))
  tp <- do.call(rbind, lapply(levels(tl$score), function(s) data.frame(score = s, p = stats::wilcox.test(v ~ g, data = tl[tl$score == s, ])$p.value)))
  for (i in seq_len(nrow(tp))) add_stat(8, "C", paste("TIDE component, MTRS high vs low, Wilcoxon:", tp$score[i]), exact_p(tp$p[i]),
                                        sprintf("median high %.2f vs low %.2f; n = %d", stats::median(tl$v[tl$score == tp$score[i] & tl$g == "High"]),
                                                stats::median(tl$v[tl$score == tp$score[i] & tl$g == "Low"]), nrow(td)))
  tp$score <- factor(tp$score, levels = levels(tl$score))
  pC <- ggplot2::ggplot(tl, ggplot2::aes(g, v, fill = g)) +
    ggplot2::geom_boxplot(outlier.size = 0.2, width = 0.6, linewidth = 0.3) +
    ggplot2::geom_text(data = tp, ggplot2::aes(x = 1.5, y = Inf, label = lab_p(p)), vjust = 1.3, inherit.aes = FALSE,
                       size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::facet_wrap(~ score, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = PAL2, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::labs(x = NULL, y = "Score", title = "TIDE, TCGA-KIRC")

  ici <- utils::read.csv(file.path(d8, "ICI_response_summary.csv"))
  for (i in seq_len(nrow(ici))) add_stat(8, "D-E", paste("ICI response:", ici$cohort[i]),
                                         sprintf("ORR high %.1f%% vs low %.1f%%, Fisher %s", 100 * ici$ORR_High[i], 100 * ici$ORR_Low[i], exact_p(ici$p_fisher[i])),
                                         sprintf("AUC %.3f (%.3f-%.3f), Wilcoxon %s; %d evaluable, %d responders", ici$AUC[i], ici$AUC_lower[i],
                                                 ici$AUC_upper[i], exact_p(ici$p_wilcox[i]), ici$n_evaluable[i], ici$n_responders[i]))
  orr <- rbind(
    data.frame(cohort = ici$cohort, g = "High", orr = ici$ORR_High),
    data.frame(cohort = ici$cohort, g = "Low", orr = ici$ORR_Low)
  )
  orr$g <- factor(orr$g, levels = c("Low", "High")); orr$cohort <- factor(orr$cohort, levels = ici$cohort)
  ici$cohort <- factor(ici$cohort, levels = ici$cohort)
  pD <- ggplot2::ggplot(orr, ggplot2::aes(g, orr, fill = g)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * orr)), vjust = -0.4, size = (BASE - 1) / ggplot2::.pt, family = FONT) +
    ggplot2::geom_text(data = ici, ggplot2::aes(x = 1.5, y = Inf, label = paste("Fisher", lab_p(p_fisher))), vjust = 1.3,
                       inherit.aes = FALSE, size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::facet_wrap(~ cohort, nrow = 1) +
    ggplot2::scale_fill_manual(values = PAL2, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 0.55), expand = c(0, 0)) +
    ggplot2::labs(x = "MTRS group", y = "Objective response rate", title = "Anti-PD-1/PD-L1 therapy")

  roc_df <- do.call(rbind, lapply(as.character(ici$cohort), function(n) {
    e <- rlrs_cohort(n); e <- e[e$response %in% c("CR/PR", "SD/PD"), ]
    ro <- pROC::roc(factor(e$response, levels = c("SD/PD", "CR/PR")), e$RLRS, levels = c("SD/PD", "CR/PR"), direction = "<", quiet = TRUE)
    data.frame(cohort = n, fpr = 1 - ro$specificities, tpr = ro$sensitivities, auc = as.numeric(ro$auc))
  }))
  roc_df <- roc_df[order(roc_df$cohort, roc_df$fpr, roc_df$tpr), ]
  al <- stats::aggregate(auc ~ cohort, data = roc_df, FUN = function(x) x[1])
  al$lab <- sprintf("%s AUC %.2f", al$cohort, al$auc)
  coh_order <- as.character(ici$cohort)
  roc_df$cohort <- factor(roc_df$cohort, levels = coh_order, labels = al$lab[match(coh_order, al$cohort)])
  pE <- ggplot2::ggplot(roc_df, ggplot2::aes(fpr, tpr, colour = cohort)) +
    ggplot2::geom_abline(linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_path(linewidth = 0.45) +
    ggplot2::scale_colour_manual(values = stats::setNames(c("#C0392B", "#8E44AD", "#2E6DA4")[seq_along(levels(roc_df$cohort))], levels(roc_df$cohort)), name = NULL) +
    ggplot2::coord_equal() + ggplot2::labs(x = "1 - specificity", y = "Sensitivity", title = "Response (CR/PR)") +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0),
                   legend.background = ggplot2::element_blank())

  # F：GSE73731（n = 265）ccRCC Fuhrman 分级与分期中的 MTRS
  g73 <- rlrs_cohort("GSE73731")
  gr_raw <- toupper(trimws(as.character(g73$grade)))
  gr_clean <- ifelse(grepl("1|I\\b", gr_raw) & !grepl("2|3|4|II|III|IV", gr_raw), "G1",
              ifelse(grepl("2|II\\b", gr_raw) & !grepl("3|4|III|IV", gr_raw), "G2",
              ifelse(grepl("3|III\\b", gr_raw) & !grepl("4|IV", gr_raw), "G3",
              ifelse(grepl("4|IV\\b", gr_raw), "G4", NA_character_))))
  g73$grade_clean <- factor(gr_clean, levels = c("G1", "G2", "G3", "G4"))
  g73_sub <- g73[!is.na(g73$grade_clean), ]
  kw_g73 <- stats::kruskal.test(RLRS_z ~ grade_clean, data = g73_sub)$p.value
  add_stat(8, "F", "GSE73731 MTRS across Fuhrman grades G1-G4, Kruskal-Wallis P", signif(kw_g73, 3),
           sprintf("n = %d", nrow(g73_sub)))
  pF <- ggplot2::ggplot(g73_sub, ggplot2::aes(grade_clean, RLRS_z, colour = grade_clean)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
    ggplot2::geom_jitter(width = 0.15, size = 0.35, alpha = 0.6) +
    ggplot2::annotate("text", x = 2.5, y = Inf, vjust = 1.3, label = paste("KW", lab_p(kw_g73)),
                      size = (BASE - 1) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_colour_manual(values = c(G1 = "#2E6DA4", G2 = "#5DADE2", G3 = "#E67E22", G4 = "#C0392B"), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
    ggplot2::labs(x = "Fuhrman grade", y = "MTRS (z-score)", title = "GSE73731 (n = 256)")

  # G：GSE53757（n = 72 tumors）AJCC 分期中的 MTRS
  g53 <- rlrs_cohort("GSE53757")
  g53_t <- g53[g53$group == "Tumor" & !is.na(g53$stage), ]
  st_raw <- toupper(trimws(as.character(g53_t$stage)))
  st_clean <- ifelse(grepl("IV|4", st_raw), "IV",
              ifelse(grepl("III|3", st_raw), "III",
              ifelse(grepl("II|2", st_raw), "II",
              ifelse(grepl("I|1", st_raw), "I", NA_character_))))
  g53_t$stage <- factor(st_clean, levels = c("I", "II", "III", "IV"))
  g53_t <- g53_t[!is.na(g53_t$stage), ]
  kw_g53 <- stats::kruskal.test(RLRS_z ~ stage, data = g53_t)$p.value
  add_stat(8, "G", "GSE53757 MTRS across AJCC stages I-IV, Kruskal-Wallis P", signif(kw_g53, 3),
           sprintf("n = %d", nrow(g53_t)))
  pG <- ggplot2::ggplot(g53_t, ggplot2::aes(stage, RLRS_z, colour = stage)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
    ggplot2::geom_jitter(width = 0.15, size = 0.5, alpha = 0.7) +
    ggplot2::annotate("text", x = 2.5, y = Inf, vjust = 1.3, label = paste("KW", lab_p(kw_g53)),
                      size = (BASE - 1) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_colour_manual(values = c(I = "#2E6DA4", II = "#5DADE2", III = "#E67E22", IV = "#C0392B"), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
    ggplot2::labs(x = "AJCC stage", y = "MTRS (z-score)", title = "GSE53757 tumours (n = 72)")

  # H：oncoPredict 预测的药物 IC50（聚焦肾癌靶向药与显著差异药物）
  f_op <- file.path(d8, "oncoPredict_High_vs_Low.csv")
  pH <- NULL
  if (file.exists(f_op)) {
    op <- utils::read.csv(f_op, check.names = FALSE)
    pred <- utils::read.csv(file.path(d8, "oncoPredict", "calcPhenotype_Output", "DrugPredictions.csv"), row.names = 1, check.names = FALSE)
    add_stat(8, "H", "drugs with predicted IC50 differing between MTRS groups (FDR < 0.05): lower / higher in MTRS-high",
             sprintf("%d / %d", sum(op$FDR < 0.05 & op$median_High < op$median_Low), sum(op$FDR < 0.05 & op$median_High > op$median_Low)),
             paste("of", nrow(op), "drugs"))
    cand_drugs <- grep("^(Sorafenib|Axitinib|Rapamycin|Temsirolimus|AZD8055|Dactolisib|Olaparib|Cisplatin|Gemcitabine|Gefitinib|Erlotinib)_",
                       colnames(pred), value = TRUE, ignore.case = TRUE)
    cand_drugs <- cand_drugs[!duplicated(sub("_.*", "", cand_drugs))]
    sig_cand <- cand_drugs[cand_drugs %in% op$variable[op$FDR < 0.05]]
    if (length(sig_cand) < 6) {
      top_op <- op$variable[!op$variable %in% sig_cand][seq_len(6 - length(sig_cand))]
      focus <- c(sig_cand, top_op)[1:6]
    } else {
      focus <- sig_cand[1:6]
    }
    for (v in focus) add_stat(8, "H", paste("predicted IC50, MTRS high vs low:", v),
                              sprintf("median %.3g vs %.3g", op$median_High[op$variable == v], op$median_Low[op$variable == v]),
                              lab_q(op$FDR[op$variable == v]))
    dl <- do.call(rbind, lapply(focus, function(v) data.frame(drug = sub("_.*", "", v), id = rownames(pred), v = log10(pred[[v]]))))
    dl$g <- factor(rt$RLRS_group[match(dl$id, rt$ID)], levels = c("Low", "High"))
    dl <- dl[!is.na(dl$g), ]
    qq <- data.frame(drug = sub("_.*", "", focus), q = op$FDR[match(focus, op$variable)])
    pH <- ggplot2::ggplot(dl, ggplot2::aes(g, v, fill = g)) +
      ggplot2::geom_boxplot(outlier.size = 0.2, width = 0.6, linewidth = 0.3) +
      ggplot2::geom_text(data = qq, ggplot2::aes(x = 1.5, y = Inf, label = lab_q(q)), vjust = 1.3, inherit.aes = FALSE,
                         size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
      ggplot2::facet_wrap(~ drug, ncol = 3, scales = "free_y") +
      ggplot2::scale_fill_manual(values = PAL2, guide = "none") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
      ggplot2::labs(x = "MTRS group", y = "Predicted log10 IC50", title = "Predicted drug sensitivity (GDSC2, oncoPredict)")
  }
  pl <- Filter(Negate(is.null), list(pA, pB, pC, pD, pE, pF, pG, pH))
  pl <- Map(function(g, t) g + ggplot2::labs(tag = t), pl, LETTERS[seq_along(pl)])
  top <- patchwork::wrap_elements(full = (pl[[1]] | (pl[[2]] / pl[[3]] + patchwork::plot_layout(heights = c(1, 0.8)))) +
                                    patchwork::plot_layout(widths = c(1.15, 1)))
  mid <- patchwork::wrap_elements(full = (pl[[4]] | pl[[5]] | pl[[6]]) + patchwork::plot_layout(widths = c(1.15, 0.9, 0.95)))
  bot <- patchwork::wrap_elements(full = if (length(pl) == 8) (pl[[7]] | pl[[8]]) + patchwork::plot_layout(widths = c(0.75, 1.25)) else pl[[7]])
  p <- patchwork::wrap_plots(top, mid, bot, ncol = 1, heights = c(1.3, 0.85, 1.0))
  save_fig(p, 8, 180, 245)
}

# =============================================================================
# Figure 9：核心基因（SHAP）→ TCGA 配对 → 跨队列预后 → 单细胞 → 空间转录组（GSE210041）
# =============================================================================
fig9 <- function() {
  d9 <- file.path(DIR_RES, "09_hub_gene_spatial")
  rk <- utils::read.csv(file.path(d9, "SHAP_gene_ranking.csv"))
  cf <- utils::read.csv(file.path(d9, "RLRS_coefficients.csv"))
  genes <- rk$gene; hubs <- utils::head(genes, 3)
  HUB_COLS <- stats::setNames(c(C_HI, "#E67E22", "#2E6DA4"), hubs)

  tcga <- readRDS(file.path(DIR_PRO, "TCGA_KIRC.rds"))
  mi <- build_mime_input(tcga$expr, tcga$clin[tcga$clin$group == "Tumor", ], genes)
  r0 <- rlrs_all(); mi <- merge(r0[r0$cohort == "TCGA_KIRC", c("ID", "RLRS")], mi, by = "ID")
  X <- as.matrix(mi[, safe_names(genes)]); colnames(X) <- genes
  w <- stats::setNames(cf$coefficient[match(genes, cf$term)], genes)
  S <- sweep(sweep(X, 2, colMeans(X)), 2, w, `*`)
  for (i in seq_len(nrow(rk))) add_stat(9, "A", paste("mean |SHAP|, share of total:", rk$gene[i]),
                                        sprintf("%.2f (%.1f%%)", rk$mean_abs_SHAP[i], rk$share_pct[i]),
                                        sprintf("coefficient %.2f", rk$coefficient[i]))
  ylab <- stats::setNames(sprintf("%s (%.1f%%)", rk$gene, rk$share_pct), rk$gene)
  sl <- data.frame(
    gene = factor(rep(genes, each = nrow(S)), levels = rev(genes)),
    shap = as.vector(S),
    z = pmax(pmin(as.vector(X), 2.5), -2.5)
  )
  pA <- ggplot2::ggplot(sl, ggplot2::aes(shap, gene, colour = z)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_point(position = ggplot2::position_jitter(height = 0.28, width = 0, seed = 1), size = 0.2, stroke = 0) +
    ggplot2::scale_colour_gradient2(low = C_LO, mid = "grey85", high = C_HI, name = "Expression\n(z-score)", limits = c(-2.5, 2.5)) +
    ggplot2::scale_y_discrete(labels = ylab) +
    ggplot2::labs(x = "SHAP value (contribution to MTRS)", y = NULL, title = "Gene contributions, TCGA-KIRC") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"), legend.position = "bottom",
                   legend.key.width = ggplot2::unit(5, "mm"), legend.key.height = ggplot2::unit(2, "mm"))

  cl <- tcga$clin
  pat <- intersect(cl$patient[cl$group == "Tumor"], cl$patient[cl$group == "Normal"])
  pr <- cl[cl$patient %in% pat, ]; pr <- pr[!duplicated(pr[, c("patient", "group")]), ]
  bl <- do.call(rbind, lapply(hubs, function(g) data.frame(gene = g, patient = pr$patient, group = pr$group,
                                                          v = as.numeric(tcga$expr[g, pr$ID]))))
  bl$group <- factor(bl$group, levels = c("Normal", "Tumor"), labels = c("Normal", "Tumour")); bl$gene <- factor(bl$gene, levels = hubs)
  pvB <- vapply(hubs, function(g) {
    x <- stats::reshape(bl[bl$gene == g, c("patient", "group", "v")], idvar = "patient", timevar = "group", direction = "wide")
    stats::wilcox.test(x$v.Tumour, x$v.Normal, paired = TRUE)$p.value
  }, numeric(1))
  for (g in hubs) {
    xa <- as.numeric(tcga$expr[g, cl$ID[cl$group == "Tumor"]]); xb <- as.numeric(tcga$expr[g, cl$ID[cl$group == "Normal"]])
    add_stat(9, "B", paste("tumour vs normal, paired Wilcoxon P (all-sample Wilcoxon P):", g), signif(pvB[[g]], 3),
             sprintf("%d pairs; all samples P = %s; median tumour %.2f vs normal %.2f", length(pat),
                     signif(stats::wilcox.test(xa, xb)$p.value, 3), stats::median(xa), stats::median(xb)))
  }
  pB <- ggplot2::ggplot(bl, ggplot2::aes(group, v)) +
    ggplot2::geom_line(ggplot2::aes(group = patient), colour = "grey75", linewidth = 0.2) +
    ggplot2::geom_point(ggplot2::aes(colour = group), size = 0.6) +
    ggplot2::geom_text(data = data.frame(gene = factor(hubs, levels = hubs), lab = lab_p(pvB)),
                       ggplot2::aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, inherit.aes = FALSE,
                       size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    ggplot2::scale_colour_manual(values = c(Normal = C_LO, Tumour = C_HI), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.2))) +
    ggplot2::labs(x = NULL, y = "log2(TPM + 1)", title = sprintf("TCGA-KIRC paired samples (n = %d)", length(pat))) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold.italic"),
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  hr <- list()
  for (n in SURV_COH) {
    o <- if (n == "TCGA_KIRC") tcga else readRDS(file.path(DIR_PRO, paste0(n, ".rds")))
    c0 <- coh_clin(n); c0 <- c0[!is.na(c0$OS.time) & c0$OS.time > 0 & !is.na(c0$OS) & c0$ID %in% colnames(o$expr), ]
    for (g in intersect(hubs, rownames(o$expr))) {
      z <- as.numeric(scale(as.numeric(o$expr[g, c0$ID])))
      cs <- summary(survival::coxph(survival::Surv(c0$OS.time, c0$OS) ~ z))
      hr[[paste(n, g)]] <- data.frame(cohort = n, gene = g, HR = cs$conf.int[1, 1], lower = cs$conf.int[1, 3],
                                      upper = cs$conf.int[1, 4], p = cs$coefficients[1, 5], n = nrow(c0), events = sum(c0$OS))
    }
  }
  hr <- do.call(rbind, hr)
  for (i in seq_len(nrow(hr))) add_stat(9, "C", paste("HR per SD:", hr$gene[i], hr$cohort[i]),
                                        sprintf("%.2f (%.2f-%.2f)", hr$HR[i], hr$lower[i], hr$upper[i]), exact_p(hr$p[i]))
  hr$lab <- factor(COH_LAB[hr$cohort], levels = COH_LAB[SURV_COH]); hr$gene <- factor(hr$gene, levels = hubs)
  pC <- forest_gg(hr, label = "lab", group = "gene", cols = HUB_COLS) +
    ggplot2::labs(title = "Overall survival (DFS in GSE22541)") +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(face = "italic"))

  cc <- utils::read.csv(file.path(d9, "hub_scRNA_by_cell_class.csv"))
  CLS <- c(Malignant = "Malignant epithelial", Normal_epi = "Normal epithelial", Endothelial = "Endothelial",
           Myofib_Peri = "Myofibroblast/pericyte", T_cell = "T/NK cell", Myeloid = "Myeloid", Mast = "Mast cell")
  cc$cls <- factor(ifelse(cc$cell_class %in% names(CLS), CLS[cc$cell_class], cc$cell_class),
                   levels = rev(CLS[CLS %in% CLS[cc$cell_class]]))
  cc$gene <- factor(cc$gene, levels = hubs)
  cc$avg_s <- stats::ave(cc$avg, cc$gene, FUN = function(v) (v - min(v)) / max(diff(range(v)), 1e-9))
  for (i in seq_len(nrow(cc))) add_stat(9, "D", paste("fraction expressing / mean expression:", cc$gene[i], cc$cell_class[i]),
                                        sprintf("%.3f / %.3f", cc$pct[i], cc$avg[i]))
  pD <- ggplot2::ggplot(cc, ggplot2::aes(gene, cls)) +
    ggplot2::geom_point(ggplot2::aes(size = pct, colour = avg_s)) +
    ggplot2::scale_size_area(max_size = 3.2, name = "Fraction\nexpressing", breaks = c(0.05, 0.25, 0.5)) +
    ggplot2::scale_colour_gradient(low = "grey90", high = C_HI, name = "Scaled\nexpression", breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = "Single-cell expression") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"), legend.box = "horizontal", legend.position = "right")

  hv <- utils::read.csv(file.path(d9, "hub_scRNA_malignant_High_vs_Low_per_sample.csv"), check.names = FALSE)
  el <- do.call(rbind, lapply(hubs, function(g) data.frame(gene = g, sample = hv$sample, grp = hv$group, v = hv[[g]])))
  el$grp <- factor(el$grp, levels = c("Low", "High")); el$gene <- factor(el$gene, levels = hubs)
  pvE <- vapply(hubs, function(g) {
    x <- stats::reshape(el[el$gene == g, c("sample", "grp", "v")], idvar = "sample", timevar = "grp", direction = "wide")
    suppressWarnings(stats::wilcox.test(x$v.High, x$v.Low, paired = TRUE, exact = FALSE)$p.value)
  }, numeric(1))
  for (g in hubs) {
    x <- stats::reshape(el[el$gene == g, c("sample", "grp", "v")], idvar = "sample", timevar = "grp", direction = "wide")
    add_stat(9, "E", paste("MT-high vs -low tumour epithelial cells, paired Wilcoxon P:", g), signif(pvE[[g]], 3),
             sprintf("higher in MT-high in %d of %d samples", sum(x$v.High > x$v.Low), nrow(x)))
  }
  pE <- ggplot2::ggplot(el, ggplot2::aes(grp, v)) +
    ggplot2::geom_line(ggplot2::aes(group = sample), colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = grp), size = 0.9) +
    ggplot2::geom_text(data = data.frame(gene = factor(hubs, levels = hubs), lab = lab_p(pvE)),
                       ggplot2::aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, inherit.aes = FALSE,
                       size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    ggplot2::scale_colour_manual(values = PAL2, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.2))) +
    ggplot2::labs(x = "MT_core group", y = "Mean expression per sample", title = "Tumour epithelial cells, samples as the unit") +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold.italic"))

  # F：空间分布（GSE210041 Visium 代表切片）
  slim <- readRDS(file.path(DIR_PRO, "spatial_spots_slim.rds"))
  show <- unique(slim$sample)[1]
  feats <- c(MT_core = "MT_core", EMT = "EMT", Angiogenesis = "Angiogenesis", stats::setNames(hubs[1], hubs[1]))
  feats <- feats[names(feats) %in% colnames(slim)]
  d_sp <- slim[slim$sample == show, ]
  fl <- do.call(rbind, lapply(names(feats), function(f) {
    data.frame(sample = show, feature = feats[[f]], x = d_sp$x, y = d_sp$y,
               v = rank(d_sp[[f]], ties.method = "average") / nrow(d_sp))
  }))
  fl$feature <- factor(fl$feature, levels = feats)
  sp_rho <- function(a, b) stats::cor(d_sp[[a]], d_sp[[b]], method = "spearman", use = "complete.obs")
  add_stat(9, "F", paste("Visium section", show, "spot-level rho MT_core vs EMT / Angiogenesis /", hubs[1]),
           sprintf("%.2f / %.2f / %.2f", sp_rho("MT_core", "EMT"), sp_rho("MT_core", "Angiogenesis"), sp_rho("MT_core", hubs[1])),
           paste(nrow(d_sp), "spots"))
  pF <- ggplot2::ggplot(fl, ggplot2::aes(y, -x, colour = v)) +
    ggplot2::geom_point(size = 0.4, stroke = 0, shape = 16) +
    ggplot2::facet_wrap(~ feature, nrow = 1) + ggplot2::coord_fixed() +
    ggplot2::scale_colour_viridis_c(option = "inferno", name = "Percentile\nwithin section",
                                    breaks = c(0, 0.5, 1), labels = c("0", "50", "100")) +
    ggplot2::labs(x = NULL, y = NULL, title = sprintf("ccRCC Visium section %s (%s spots, GSE210041)", show,
                                                       format(nrow(d_sp), big.mark = ","))) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())

  # G：空间转录组切片中 MT_core 与各功能通路的相关及偏相关
  sp_p <- utils::read.csv(file.path(d9, "Spatial_RLcore_vs_programs_by_slice.csv"))
  dc   <- utils::read.csv(file.path(d9, "Spatial_depth_check_by_slice.csv"))
  gv <- c(
    rho_RLcore_EMT = "EMT",
    rho_RLcore_Angio = "Angiogenesis",
    rho_RLcore_IFNg = "IFN-γ response",
    prho_RLcore_IFNg_given_Prolif_UMI = "IFN-γ | proliferation, UMI",
    rho_RLcore_GEP18 = "T cell-inflamed GEP",
    prho_RLcore_GEP18_given_Prolif_UMI = "GEP | proliferation, UMI",
    rho_RLcore_Prolif = "Proliferation"
  )
  m_sp <- merge(sp_p, dc, by = c("sample", "n_spots"), suffixes = c("", "_dc"))
  gl <- do.call(rbind, lapply(names(gv), function(v) {
    if (!v %in% colnames(m_sp)) return(NULL)
    data.frame(sample = m_sp$sample, what = gv[[v]], rho = m_sp[[v]], adj = grepl("^prho", v))
  }))
  gl$what <- factor(gl$what, levels = rev(gv[names(gv) %in% colnames(m_sp)]))
  for (v in intersect(names(gv), colnames(m_sp))) {
    add_stat(9, "G", paste("spot-level correlation with MT_core across Visium sections:", gv[[v]]),
             paste(sprintf("%s: %.2f", m_sp$sample, m_sp[[v]]), collapse = "; "))
  }
  pG <- ggplot2::ggplot(gl, ggplot2::aes(rho, what, colour = sample, shape = adj)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_point(size = 1.8, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::scale_colour_manual(values = c(ccRCC_Visium_1 = C_HI, ccRCC_Visium_2 = C_LO), name = "Visium section") +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    ggplot2::labs(x = "Spot-level correlation with MT_core", y = NULL, title = "Programmes (GSE210041)") +
    ggplot2::theme(legend.position = "bottom")

  # H：核心基因与 MT_core 的 spot 水平相关（未校正 vs 校正增殖和深度）
  hp <- utils::read.csv(file.path(d9, "Spatial_hub_partial_by_slice.csv"))
  hl <- rbind(
    data.frame(sample = hp$sample, gene = hp$gene, what = "Unadjusted", rho = hp$rho),
    data.frame(sample = hp$sample, gene = hp$gene, what = "Adjusted for proliferation & UMI", rho = hp$prho_given_Prolif_UMI)
  )
  hl <- hl[!is.na(hl$rho), ]
  hl$gene <- factor(hl$gene, levels = rev(hubs))
  hl$what <- factor(hl$what, levels = c("Unadjusted", "Adjusted for proliferation & UMI"))
  for (g in hubs) {
    h <- hp[hp$gene == g, ]
    add_stat(9, "H", paste("spot-level rho with MT_core (unadjusted / adjusted):", g),
             paste(sprintf("%s: %.2f / %.2f (detected %.1f%%)", h$sample, h$rho, h$prho_given_Prolif_UMI, 100 * h$pct_detected), collapse = "; "))
  }
  pH <- ggplot2::ggplot(hl, ggplot2::aes(rho, gene, fill = what)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65,
                      data = stats::aggregate(rho ~ gene + what, data = hl, FUN = mean)) +
    ggplot2::geom_point(ggplot2::aes(group = what), position = ggplot2::position_dodge(width = 0.75), size = 1.2, colour = "black") +
    ggplot2::scale_fill_manual(values = c(Unadjusted = "grey75", `Adjusted for proliferation & UMI` = "#F5B7B1"), name = NULL) +
    ggplot2::labs(x = "Spot-level correlation with MT_core", y = NULL, title = "Model genes in Visium") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"), legend.position = "bottom", legend.direction = "vertical")

  r1 <- (pA | (pB / pC + patchwork::plot_layout(heights = c(0.8, 1)))) + patchwork::plot_layout(widths = c(1, 1.1))
  r2 <- (pD | pE) + patchwork::plot_layout(widths = c(1, 1.15))
  r4 <- (pG | pH) + patchwork::plot_layout(widths = c(1.2, 1))
  p <- (r1 / r2 / pF / r4) + patchwork::plot_layout(heights = c(1.45, 0.8, 0.62, 0.8)) +
    patchwork::plot_annotation(tag_levels = "A")
  save_fig(p, 9, 180, 240)
}

# =============================================================================
# 补充图 S1：单细胞聚类、样本分布、各簇 marker、各样本细胞组成与质控指标
# =============================================================================
figS1 <- function() {
  s <- load_sc(); md0 <- s$md
  lev <- as.character(sort(unique(as.integer(as.character(md0$seurat_clusters)))))
  md0$cluster <- factor(as.character(md0$seurat_clusters), levels = lev)
  set.seed(PARAM$seed); md <- md0[sample(nrow(md0)), ]
  cen <- stats::aggregate(cbind(UMAP1, UMAP2) ~ cluster, data = md, FUN = stats::median)
  add_stat("S1", "A", "clusters", nlevels(md$cluster), "Harmony integration, resolution 0.8")
  noax <- ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  pA <- ggplot2::ggplot(md, ggplot2::aes(UMAP1, UMAP2, colour = cluster)) +
    ggplot2::geom_point(size = 0.25, stroke = 0) +
    ggrepel::geom_text_repel(data = cen, ggplot2::aes(label = cluster), colour = "black", size = BASE / ggplot2::.pt, family = FONT,
                             seed = 1, box.padding = 0.15, min.segment.length = 0.3, segment.size = 0.2,
                             bg.color = "white", bg.r = 0.12) +
    ggplot2::scale_colour_manual(values = grDevices::hcl.colors(nlevels(md$cluster), "Dark 3"), guide = "none") +
    ggplot2::coord_fixed() + ggplot2::labs(title = sprintf("%d clusters", nlevels(md$cluster))) + noax

  pB <- ggplot2::ggplot(md, ggplot2::aes(UMAP1, UMAP2, colour = sample)) +
    ggplot2::geom_point(size = 0.25, stroke = 0) +
    ggplot2::scale_colour_manual(values = grDevices::hcl.colors(length(unique(md$sample)), "Spectral"), name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 1.8), ncol = 2)) +
    ggplot2::coord_fixed() + ggplot2::labs(title = "Samples (7 ccRCC, 6 normal)") + noax +
    ggplot2::theme(legend.position = "right", legend.key.size = ggplot2::unit(2.2, "mm"))

  ex <- s$mk; mk <- unlist(MK_SHOW); mk <- mk[mk %in% rownames(ex)]
  ann <- tapply(as.character(md0$ct), md0$cluster, function(v) names(sort(table(v), decreasing = TRUE))[1])
  dp <- do.call(rbind, lapply(lev, function(k) {
    m <- ex[mk, md0$cluster == k, drop = FALSE]
    data.frame(cluster = k, gene = mk, avg = Matrix::rowMeans(m), pct = Matrix::rowMeans(m > 0))
  }))
  dp$avg_s <- stats::ave(dp$avg, dp$gene, FUN = function(v) (v - min(v)) / max(diff(range(v)), 1e-9))
  ylab <- paste0(lev, "  ", ann[lev])
  dp$y <- factor(paste0(dp$cluster, "  ", ann[dp$cluster]), levels = rev(ylab))
  dp$gene <- factor(dp$gene, levels = unique(mk))
  for (k in lev) add_stat("S1", "C", paste("cluster annotation:", k), ann[[k]], paste(sum(md0$cluster == k), "cells"))
  pC <- ggplot2::ggplot(dp, ggplot2::aes(gene, y)) +
    ggplot2::geom_point(ggplot2::aes(size = pct, colour = avg_s)) +
    ggplot2::scale_size_area(max_size = 2.4, name = "Fraction\nexpressing", breaks = c(0.25, 0.5, 0.75)) +
    ggplot2::scale_colour_gradient(low = "grey90", high = C_HI, name = "Scaled\nexpression", breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = "Marker genes by cluster") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
                   axis.text.y = ggplot2::element_text(size = BASE - 1.5), legend.box = "horizontal")

  cp <- as.data.frame(table(sample = md0$sample, ct = md0$ct))
  tot <- tapply(cp$Freq, cp$sample, sum); cp$frac <- cp$Freq / tot[as.character(cp$sample)]
  tis <- tapply(as.character(md0$tissue), md0$sample, function(v) v[1])
  cp$tissue <- factor(ifelse(tis[as.character(cp$sample)] == "Tumor", "Tumour", "Normal"), levels = c("Tumour", "Normal"))
  for (k in names(tot)) add_stat("S1", "D", paste("cells after QC:", k), tot[[k]])
  pD <- ggplot2::ggplot(cp, ggplot2::aes(sample, frac, fill = ct)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::facet_grid(~ tissue, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = CT_COLS, name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Cells", title = "Cell-type composition") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 1.5), strip.clip = "off",
                   legend.key.size = ggplot2::unit(2.4, "mm"), legend.text = ggplot2::element_text(size = BASE - 1.5))

  mtc <- grep("percent.?mt|pct.?mt|mito", names(md0), value = TRUE, ignore.case = TRUE)[1]
  qv <- c(nFeature_RNA = "Detected genes", stats::setNames("Mitochondrial reads (%)", mtc))
  qv <- qv[names(qv) %in% names(md0)]
  ql <- do.call(rbind, lapply(names(qv), function(v) data.frame(sample = md0$sample, metric = qv[[v]], value = md0[[v]])))
  ql$sample <- factor(ql$sample, levels = levels(factor(cp$sample)))
  pE <- ggplot2::ggplot(ql, ggplot2::aes(sample, value)) +
    ggplot2::geom_violin(scale = "width", fill = "grey80", colour = "grey40", linewidth = 0.2) +
    ggplot2::facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, title = "Quality control") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 1.5))

  pl <- Map(function(g, t) g + ggplot2::labs(tag = t), list(pA, pB, pC, pD, pE), LETTERS[1:5])
  p <- patchwork::wrap_plots(
    patchwork::wrap_elements(full = pl[[1]] | pl[[2]]),
    patchwork::wrap_elements(full = pl[[3]]),
    patchwork::wrap_elements(full = (pl[[4]] | pl[[5]]) + patchwork::plot_layout(widths = c(1.25, 1))),
    ncol = 1, heights = c(1, 1.1, 0.95)
  )
  save_fig(p, "S1", 180, 235)
}

# =============================================================================
# 补充图 S3：Bulk 队列中 MT_core 与 MTI 按分期、分级和体细胞突变的分布
# =============================================================================
figS3 <- function() {
  tum <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
  t1 <- tum[tum$group == "Tumor", ]
  mt <- utils::read.csv(file.path(DIR_RES, "04_bulk_score", "MT_score_by_mutation.csv"))
  sc2 <- c(MT_core = "MT_core", MTI = "MTI")
  box_by <- function(d, var, lev, test = c("kw", "wilcox"), title) {
    test <- match.arg(test)
    d <- d[!is.na(d[[var]]) & d[[var]] %in% lev, ]
    l <- do.call(rbind, lapply(names(sc2), function(s) data.frame(score = sc2[[s]], grp = factor(d[[var]], levels = lev), v = d[[s]])))
    l$score <- factor(l$score, levels = sc2)
    pv <- vapply(sc2, function(s) {
      x <- l[l$score == s, ]
      if (test == "kw") stats::kruskal.test(v ~ grp, data = x)$p.value else stats::wilcox.test(v ~ grp, data = x)$p.value
    }, numeric(1))
    for (s in sc2) add_stat("S3", title, paste(s, "by", var, if (test == "kw") "Kruskal-Wallis P" else "Wilcoxon P"),
                            signif(pv[[s]], 3), paste(sprintf("%s n = %d", lev, as.integer(table(d[[var]])[lev])), collapse = "; "))
    ggplot2::ggplot(l, ggplot2::aes(grp, v)) +
      ggplot2::geom_boxplot(ggplot2::aes(colour = grp), outlier.shape = NA, width = 0.6, linewidth = 0.35) +
      ggplot2::geom_jitter(ggplot2::aes(colour = grp), width = 0.15, size = 0.25, alpha = 0.6) +
      ggplot2::geom_text(data = data.frame(score = factor(sc2, levels = sc2), lab = lab_p(pv)),
                         ggplot2::aes(x = (length(lev) + 1) / 2, y = Inf, label = lab), vjust = 1.3, inherit.aes = FALSE,
                         size = (BASE - 1) / ggplot2::.pt, family = FONT) +
      ggplot2::facet_wrap(~ score, scales = "free_y", nrow = 1) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
      ggplot2::labs(x = NULL, y = "ssGSEA score") + ggplot2::theme(legend.position = "none")
  }
  pA <- box_by(t1, "stage", c("I", "II", "III", "IV"), "kw", "A") +
    ggplot2::scale_colour_manual(values = c(I = "#2E6DA4", II = "#5DADE2", III = "#E67E22", IV = C_HI)) +
    ggplot2::labs(title = "TCGA-KIRC AJCC stage (I-IV)")
  pB <- box_by(t1, "stage_bin", c("I-II", "III-IV"), "wilcox", "B") +
    ggplot2::scale_colour_manual(values = c(`I-II` = C_LO, `III-IV` = C_HI)) +
    ggplot2::labs(title = "AJCC stage (binary)")
  pC <- box_by(t1, "grade_bin", c("G1-2", "G3-4"), "wilcox", "C") +
    ggplot2::scale_colour_manual(values = c(`G1-2` = C_LO, `G3-4` = C_HI)) +
    ggplot2::labs(title = "Fuhrman grade (binary)")

  gm <- intersect(c("VHL", "PBRM1", "SETD2", "BAP1", "MTOR", "TP53"), mt$gene)
  ml <- do.call(rbind, lapply(gm, function(g) data.frame(gene = g, status = t1[[paste0("mut_", g)]], MTI = t1$MTI)))
  ml <- ml[!is.na(ml$status), ]; ml$status <- factor(ml$status, levels = c("WT", "Mut")); ml$gene <- factor(ml$gene, levels = gm)
  for (g in gm) add_stat("S3", "D", paste("MTI mutant vs wild-type, Wilcoxon FDR:", g), signif(mt$FDR_MTI[mt$gene == g], 3),
                         sprintf("%d mutant / %d tested", sum(ml$gene == g & ml$status == "Mut"), sum(ml$gene == g)))
  pD <- ggplot2::ggplot(ml, ggplot2::aes(status, MTI)) +
    ggplot2::geom_boxplot(ggplot2::aes(colour = status), outlier.shape = NA, width = 0.6, linewidth = 0.35) +
    ggplot2::geom_jitter(ggplot2::aes(colour = status), width = 0.15, size = 0.25, alpha = 0.6) +
    ggplot2::geom_text(data = data.frame(gene = factor(gm, levels = gm), lab = lab_q(mt$FDR_MTI[match(gm, mt$gene)])),
                       ggplot2::aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, inherit.aes = FALSE,
                       size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::facet_wrap(~ gene, nrow = 1) +
    ggplot2::scale_colour_manual(values = c(WT = "grey50", Mut = C_HI), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
    ggplot2::labs(x = NULL, y = "MTI (ssGSEA)", title = "MTI by somatic mutation")

  p <- (pA / (pB | pC) / pD) + patchwork::plot_layout(heights = c(1.0, 1.0, 0.9)) + patchwork::plot_annotation(tag_levels = "A")
  save_fig(p, "S3", 180, 190)
}

for (n in WHICH) {
  f <- get0(paste0("fig", n))
  if (is.null(f)) { log_msg("Figure ", n, " 的绘图函数未定义，跳过"); next }
  f()
}
log_msg("11_figures 完成")
