# =============================================================================
# 00_config.R —— KIRC_MTRG 全局路径、参数与通用函数
# =============================================================================
options(stringsAsFactors = FALSE, timeout = 3600, scipen = 999)

N_CORES <- 2L
options(mc.cores = N_CORES)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_FALLBACK = N_CORES, OMP_NUM_THREADS = N_CORES)

if (!file.exists("KIRC_MTRG.Rproj"))
  stop("请先把工作目录设到项目根目录 KIRC_MTRG/")

PROJ     <- normalizePath(".")
BLCA_DIR <- "../BLCA_Rloop"
DIR_RAW  <- file.path(PROJ, "data", "raw")
DIR_PRO  <- file.path(PROJ, "data", "processed")
DIR_GS   <- file.path(PROJ, "data", "geneset")
DIR_RES  <- file.path(PROJ, "results")
for (d in c(DIR_RAW, DIR_PRO, DIR_GS, DIR_RES)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

iobr_shared <- file.path(BLCA_DIR, "data", "raw", "IOBR_cache")
options(IOBR.cache_dir = if (dir.exists(iobr_shared)) iobr_shared else file.path(DIR_RAW, "IOBR_cache"))

if (requireNamespace("BiocParallel", quietly = TRUE))
  try(suppressWarnings(suppressMessages(BiocParallel::register(BiocParallel::SerialParam()))), silent = TRUE)

res_dir <- function(step) {
  d <- file.path(DIR_RES, step)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

PARAM <- list(
  seed            = 20261108,
  sc_min_features = 200,
  sc_max_features = 7000,
  sc_max_mt       = 25,
  sc_max_cells_per_sample = 1500,
  sc_npcs         = 30,
  sc_resolution   = 0.8,
  sc_min_cells_pb = 20,
  sc_deg_lfc      = 0.5,
  sc_deg_fdr      = 0.05,
  copykat_cores   = 2L,
  deg_lfc         = 1,
  deg_padj        = 0.05,
  wgcna_top_mad   = 5000,
  wgcna_mm        = 0.5,
  wgcna_gs        = 0.25,
  min_candidates  = 15,
  unicox_p        = 0.05,
  surv_years      = c(1, 3, 5)
)
set.seed(PARAM$seed)

MARKERS <- list(
  Epithelial  = c("EPCAM", "KRT8", "KRT18", "PAX8", "CA9", "NDUFA4L2", "ALDOB", "GPX3", "SLC17A3"),
  T_cell      = c("CD3D", "CD3E", "CD2", "CD8A", "CD4", "IL7R"),
  NK          = c("NKG7", "GNLY", "KLRD1", "KLRF1", "FGFBP2"),
  B_cell      = c("CD79A", "MS4A1", "CD19"),
  Plasma      = c("MZB1", "JCHAIN", "IGHG1"),
  Myeloid     = c("CD68", "CD14", "LYZ", "C1QA", "C1QB", "CD163", "APOE"),
  Mast        = c("TPSAB1", "CPA3", "KIT"),
  Fibroblast  = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "FAP"),
  Myofib_Peri = c("ACTA2", "RGS5", "TAGLN", "MCAM"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "PLVAP", "KDR", "FLT1")
)

GEP18 <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
           "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")

CHECKPOINTS <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1",
                 "SIGLEC15", "CD276", "VSIR", "BTLA")

ANGIO_SIG <- c("VEGFA", "KDR", "ESM1", "PECAM1", "ANGPTL4", "CD34")
MYELOID_INFLAM <- c("IL6", "CXCL1", "CXCL2", "CXCL3", "CXCL8", "PTGS2")

COL2 <- c(High = "#D7301F", Low = "#2C7FB8")

log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("缺少 R 包：", paste(miss, collapse = ", "))
  invisible(suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE)))
}

save_session <- function(step) {
  writeLines(capture.output(sessionInfo()), file.path(res_dir(step), "sessionInfo.txt"))
}

save_gg <- function(p, file, w = 6, h = 5) ggplot2::ggsave(file, p, width = w, height = h, limitsize = FALSE)

save_print <- function(obj, file, w = 7, h = 6) {
  grDevices::pdf(file, width = w, height = h, onefile = TRUE)
  on.exit(grDevices::dev.off())
  if (is.list(obj) && !inherits(obj, c("gg", "ggsurvplot"))) for (o in obj) print(o) else print(obj)
}

auto_log2 <- function(m) {
  qx <- as.numeric(stats::quantile(m, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
  need_log <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (need_log) {
    m[m < 0] <- 0
    m <- log2(m + 1)
    log_msg("检测到非 log 数据，已做 log2(x+1)")
  }
  m
}

gpl_symbol_map <- function(fd) {
  cand <- c("Symbol", "ILMN_Gene", "Gene Symbol", "GENE_SYMBOL", "Gene symbol", "gene_assignment", "ORF", "gene_symbol")
  hit <- cand[cand %in% colnames(fd)][1]
  if (is.na(hit)) stop("GPL 注释中找不到基因符号列，现有列：", paste(colnames(fd), collapse = ", "))
  sym <- as.character(fd[[hit]])
  if (hit == "gene_assignment")
    sym <- vapply(strsplit(sym, " // ", fixed = TRUE),
                  function(x) if (length(x) >= 2) trimws(x[2]) else NA_character_, character(1))
  data.frame(probe = as.character(fd$ID), symbol = trimws(sub("\\s*///.*$", "", sym)))
}

collapse_to_symbol <- function(m, symbols) {
  ok <- !is.na(symbols) & nzchar(symbols) & symbols != "---"
  m <- m[ok, , drop = FALSE]; symbols <- symbols[ok]
  ord <- order(rowMeans(m, na.rm = TRUE), decreasing = TRUE)
  m <- m[ord, , drop = FALSE]; symbols <- symbols[ord]
  keep <- !duplicated(symbols)
  out <- m[keep, , drop = FALSE]
  rownames(out) <- symbols[keep]
  out
}

probe_to_symbol <- function(expr, p2s) {
  p2s <- p2s[p2s$probe %in% rownames(expr), ]
  collapse_to_symbol(expr[p2s$probe, , drop = FALSE], p2s$symbol)
}

geo_pheno <- function(eset) {
  pd <- Biobase::pData(eset)
  ch <- grep(":ch1$", colnames(pd), value = TRUE)
  out <- pd[, c("geo_accession", "title", "source_name_ch1", ch), drop = FALSE]
  colnames(out) <- sub(":ch1$", "", colnames(out))
  out
}

MONTH2DAY <- 30.4375

zscore_rows <- function(m) {
  m <- as.matrix(m)
  s <- t(scale(t(m)))
  s[is.na(s)] <- 0
  s
}

safe_names <- function(x) gsub("[^A-Za-z0-9.]", ".", x)

build_mime_input <- function(expr, surv, genes, fill_missing = FALSE, scale = TRUE) {
  ids <- intersect(colnames(expr), surv$ID)
  g_have <- intersect(genes, rownames(expr))
  z <- if (scale) zscore_rows(expr[g_have, ids, drop = FALSE]) else as.matrix(expr[g_have, ids, drop = FALSE])
  miss <- setdiff(genes, g_have)
  if (length(miss)) {
    if (!fill_missing) stop("队列缺少基因：", paste(head(miss, 20), collapse = ", "))
    log_msg("提示：缺失 ", length(miss), " 个基因，已补 0：", paste(head(miss, 10), collapse = ", "))
    z <- rbind(z, matrix(0, nrow = length(miss), ncol = ncol(z), dimnames = list(miss, colnames(z))))
  }
  z <- z[genes, , drop = FALSE]
  s <- surv[match(ids, surv$ID), c("ID", "OS.time", "OS")]
  df <- data.frame(s, t(z), check.names = FALSE)
  colnames(df)[-(1:3)] <- safe_names(genes)
  df <- df[!is.na(df$OS.time) & !is.na(df$OS) & df$OS.time > 0, ]
  rownames(df) <- NULL
  df
}

model_genes <- function(fit, fallback, quiet = FALSE) {
  g <- tryCatch({
    if (inherits(fit, "rfsrc")) fit$xvar.names
    else if (inherits(fit, "coxph")) names(stats::coef(fit))
    else if (inherits(fit, "cv.glmnet")) {
      cf <- as.matrix(stats::coef(fit, s = fit$lambda.min)); rownames(cf)[cf[, 1] != 0]
    }
    else if (inherits(fit, "CoxBoost")) { cf <- stats::coef(fit); names(cf)[cf != 0] }
    else if (inherits(fit, "plsRcoxmodel")) colnames(fit$dataX)
    else if (inherits(fit, "survivalsvm")) fit$var.names
    else if (is.list(fit) && inherits(fit[[1]], "gbm")) fit[[1]]$var.names
    else if (is.list(fit) && inherits(fit[[1]], "glmnet")) rownames(fit[[1]]$beta)
    else NULL
  }, error = function(e) NULL)
  if (!length(g)) {
    if (!quiet) log_msg("【核对】无法从模型对象（", class(fit)[1], "）自动提取基因，暂用全部候选基因")
    g <- fallback
  }
  unique(g)
}

run_ssgsea <- function(expr, gsets) {
  expr <- as.matrix(expr)
  gsets <- lapply(gsets, function(g) intersect(g, rownames(expr)))
  gsets <- gsets[lengths(gsets) >= 5]
  if (utils::packageVersion("GSVA") >= "1.50.0") {
    GSVA::gsva(GSVA::ssgseaParam(expr, gsets), verbose = FALSE)
  } else {
    GSVA::gsva(expr, gsets, method = "ssgsea", verbose = FALSE)
  }
}

get_msig <- function(collection, subcollection = NULL) {
  df <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection),
    error = function(e) {
      if (is.null(subcollection)) msigdbr::msigdbr(species = "Homo sapiens", category = collection)
      else msigdbr::msigdbr(species = "Homo sapiens", category = collection, subcategory = subcollection)
    })
  split(df$gene_symbol, df$gs_name)
}

median_group <- function(x) factor(ifelse(x >= stats::median(x, na.rm = TRUE), "High", "Low"),
                                   levels = c("High", "Low"))

km_plot <- function(df, time = "OS.time", status = "OS", group = "group", title = "",
                    xlab = "Time (days)", palette = unname(COL2)) {
  fml <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time, status, group))
  fit <- survminer::surv_fit(fml, data = df)
  survminer::ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                        title = title, xlab = xlab, palette = palette, legend.title = "")
}

cox_table <- function(fit) {
  s <- summary(fit)
  data.frame(term = rownames(s$coefficients),
             HR = s$conf.int[, "exp(coef)"], lower = s$conf.int[, "lower .95"],
             upper = s$conf.int[, "upper .95"], p = s$coefficients[, ncol(s$coefficients)],
             row.names = NULL)
}

gunzip_to <- function(src, dest) {
  zz <- gzfile(src, "rb"); on.exit(close(zz))
  out <- file(dest, "wb"); on.exit(close(out), add = TRUE)
  repeat {
    b <- readBin(zz, "raw", n = 4e6)
    if (!length(b)) break
    writeBin(b, out)
  }
  invisible(dest)
}

get_geo_eset <- function(gse, gpl = NULL) {
  dest <- file.path(DIR_RAW, "GEO"); dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  es <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = TRUE, destdir = dest)
  if (length(es) > 1 && !is.null(gpl)) es <- es[grep(gpl, names(es))]
  es[[1]]
}

pcol <- function(pd, ...) {
  cand <- c(...)
  hit <- cand[cand %in% colnames(pd)][1]
  if (is.na(hit)) {
    hit <- colnames(pd)[match(tolower(cand), tolower(colnames(pd)))]
    hit <- hit[!is.na(hit)][1]
  }
  if (is.na(hit)) {
    log_msg("【核对】找不到字段：", paste(cand, collapse = " / "))
    return(rep(NA_character_, nrow(pd)))
  }
  trimws(as.character(pd[[hit]]))
}

# 继承 BLCA_Rloop/R/00_config.R 的代理配置与补充文件下载函数（geo_suppl / ensure_symbols / gsm_suppl）
local({
  ref_cfg <- file.path(BLCA_DIR, "R", "00_config.R")
  if (file.exists(ref_cfg)) {
    exprs <- parse(file = ref_cfg)
    for (ex in as.list(exprs)) {
      txt <- paste(deparse(ex), collapse = " ")
      if (grepl("proxyauth|^geo_suppl\\s*<-|^ensure_symbols\\s*<-|^gsm_suppl\\s*<-", txt)) {
        eval(ex, envir = .GlobalEnv)
      }
    }
  }
})
