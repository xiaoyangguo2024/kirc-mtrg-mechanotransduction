# =============================================================================
# 02_download_bulk.R —— 下载并整理 KIRC_MTRG 全部队列与单细胞/空转原始文件
#   训练队列：TCGA-KIRC（STAR counts/TPM + TCGA-CDR OS/DSS/PFI + 突变 MAF）
#   外部预后与临床验证队列：GSE29609、GSE22541、GSE167573、GSE53757、GSE73731
#   免疫治疗队列：GSE67501（RCC anti-PD-1 Nivolumab）、IMvigor210、GSE176307
#   单细胞与空间转录组预下载：GSE159115（14 例样本 h5 + 官方注释）、GSE210041（Visium）
# =============================================================================
source("R/00_config.R")
need(c("TCGAbiolinks", "SummarizedExperiment", "GEOquery", "Biobase", "data.table", "org.Hs.eg.db", "AnnotationDbi"))
STEP <- "02_download_bulk"; OUT <- res_dir(STEP)
GDC_DIR <- file.path(DIR_RAW, "GDCdata")
dir.create(GDC_DIR, recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))

COHORTS <- list()
save_cohort <- function(name, expr, clin, role, platform, endpoint) {
  expr <- as.matrix(expr); storage.mode(expr) <- "double"
  ids <- intersect(colnames(expr), clin$ID)
  if (length(ids) < ncol(expr) * 0.8)
    log_msg("【核对】", name, "：表达矩阵与临床表匹配 ", length(ids), "/", ncol(expr), " 个样本")
  clin <- clin[match(ids, clin$ID), , drop = FALSE]; rownames(clin) <- NULL
  obj <- list(expr = expr[, ids, drop = FALSE], clin = clin)
  saveRDS(obj, file.path(DIR_PRO, paste0(name, ".rds")))
  has_os <- !is.na(clin$OS.time) & !is.na(clin$OS) & clin$OS.time > 0
  COHORTS[[name]] <<- data.frame(
    cohort = name, role = role, platform = platform, endpoint = endpoint,
    n_samples = ncol(obj$expr), n_genes = nrow(obj$expr),
    n_with_survival = sum(has_os), n_events = sum(clin$OS[has_os] == 1),
    n_responders = if ("response" %in% colnames(clin)) sum(clin$response %in% c("CR/PR", "Responder", "R"), na.rm = TRUE) else NA)
  log_msg(name, "：", nrow(obj$expr), " 个基因 × ", ncol(obj$expr), " 个样本；有随访 ", sum(has_os),
          " 例，事件 ", sum(clin$OS[has_os] == 1), " 例")
  invisible(obj)
}

geo_expr_symbol <- function(es) {
  m <- auto_log2(Biobase::exprs(es))
  fd <- Biobase::fData(es)
  p2s <- tryCatch(gpl_symbol_map(fd), error = function(e) {
    for (col in colnames(fd)) {
      v <- trimws(as.character(fd[[col]]))
      if (sum(v %in% c("TP53", "VHL", "EGFR", "ACTB", "GAPDH", "MYC", "EZR", "MSN", "CAV1")) >= 3) {
        log_msg("自动识别 GPL 基因符号列：", col)
        return(data.frame(probe = as.character(fd$ID), symbol = sub("\\s*///.*$", "", v)))
      }
    }
    stop(conditionMessage(e))
  })
  probe_to_symbol(m, p2s)
}

# =============================================================================
# 1. TCGA-KIRC（训练集：表达 + 生存 + 体细胞突变）
# =============================================================================
f_tcga <- file.path(DIR_PRO, "TCGA_KIRC.rds")
if (!file.exists(f_tcga)) {
  log_msg("下载 TCGA-KIRC 表达数据（STAR - Counts）...")
  q <- TCGAbiolinks::GDCquery(project = "TCGA-KIRC", data.category = "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
  TCGAbiolinks::GDCdownload(q, directory = GDC_DIR, files.per.chunk = 50)
  se <- TCGAbiolinks::GDCprepare(q, directory = GDC_DIR)

  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  keep <- rd$gene_type == "protein_coding"
  cnt <- SummarizedExperiment::assay(se, "unstranded")[keep, ]
  tpm <- SummarizedExperiment::assay(se, "tpm_unstrand")[keep, ]
  sym <- rd$gene_name[keep]
  ord <- order(rowMeans(tpm), decreasing = TRUE)
  pick <- ord[!duplicated(sym[ord])]
  cnt <- cnt[pick, ]; tpm <- tpm[pick, ]
  rownames(cnt) <- rownames(tpm) <- sym[pick]

  bc <- sort(colnames(cnt))
  bc <- bc[substr(bc, 14, 15) %in% c("01", "11")]
  bc <- bc[!duplicated(substr(bc, 1, 15))]
  cnt <- cnt[, bc]; tpm <- tpm[, bc]
  colnames(cnt) <- colnames(tpm) <- substr(bc, 1, 15)

  cdr_file <- file.path(BLCA_DIR, "data", "raw", "xena", "Survival_SupplementalTable_S1_20171025_xena_sp")
  if (!file.exists(cdr_file)) {
    need("UCSCXenaTools")
    xe <- UCSCXenaTools::XenaGenerate(subset = XenaHostNames == "pancanAtlasHub")
    xe <- UCSCXenaTools::XenaFilter(xe, filterDatasets = "Survival_SupplementalTable_S1_20171025_xena_sp")
    xd <- UCSCXenaTools::XenaDownload(UCSCXenaTools::XenaQuery(xe), destdir = file.path(DIR_RAW, "xena"))
    cdr <- as.data.frame(UCSCXenaTools::XenaPrepare(xd))
  } else {
    cdr <- as.data.frame(data.table::fread(cdr_file))
  }
  names(cdr) <- make.names(names(cdr))
  cdr <- cdr[cdr$cancer.type.abbreviation == "KIRC", ]
  if (!nrow(cdr)) stop("TCGA-CDR 中未找到 KIRC 记录")

  stage <- sub("^Stage ([IV]+)[ABC]?$", "\\1", trimws(cdr$ajcc_pathologic_tumor_stage))
  stage[!stage %in% c("I", "II", "III", "IV")] <- NA
  grade <- trimws(cdr$histological_grade)
  grade[!grade %in% c("G1", "G2", "G3", "G4")] <- NA

  cdr_tab <- data.frame(
    ID = cdr$sample,
    OS = num(cdr$OS), OS.time = num(cdr$OS.time),
    DSS = num(cdr$DSS), DSS.time = num(cdr$DSS.time),
    PFI = num(cdr$PFI), PFI.time = num(cdr$PFI.time),
    age = num(cdr$age_at_initial_pathologic_diagnosis),
    gender = tolower(trimws(cdr$gender)),
    stage = stage,
    grade = grade
  )

  clin <- data.frame(
    ID = colnames(tpm),
    patient = substr(colnames(tpm), 1, 12),
    group = ifelse(substr(colnames(tpm), 14, 15) == "01", "Tumor", "Normal")
  )
  clin <- merge(clin, cdr_tab, by = "ID", all.x = TRUE, sort = FALSE)
  clin <- clin[match(colnames(tpm), clin$ID), ]
  clin[clin$group == "Normal", c("OS", "OS.time", "DSS", "DSS.time", "PFI", "PFI.time")] <- NA

  saveRDS(list(counts = cnt, tpm = tpm, expr = log2(tpm + 1), clin = clin), f_tcga)
}
tcga <- readRDS(f_tcga)
log_msg("TCGA-KIRC 样本构成：")
print(table(tcga$clin$group))
tc <- tcga$clin[tcga$clin$group == "Tumor", ]
has_os <- !is.na(tc$OS.time) & !is.na(tc$OS) & tc$OS.time > 0
COHORTS[["TCGA_KIRC"]] <- data.frame(
  cohort = "TCGA_KIRC", role = "training", platform = "RNA-seq (STAR, log2 TPM+1)",
  endpoint = "OS / DSS / PFI", n_samples = nrow(tc), n_genes = nrow(tcga$expr),
  n_with_survival = sum(has_os), n_events = sum(tc$OS[has_os] == 1),
  n_responders = NA
)

f_maf <- file.path(DIR_PRO, "TCGA_KIRC_maf.rds")
if (!file.exists(f_maf)) {
  log_msg("下载 TCGA-KIRC Masked Somatic Mutation MAF...")
  qm <- TCGAbiolinks::GDCquery(
    project = "TCGA-KIRC", data.category = "Simple Nucleotide Variation",
    access = "open", data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  TCGAbiolinks::GDCdownload(qm, directory = GDC_DIR)
  saveRDS(TCGAbiolinks::GDCprepare(qm, directory = GDC_DIR), f_maf)
}

# =============================================================================
# 2. GEO 独立验证队列
# =============================================================================
extract_kv <- function(es) {
  pd <- Biobase::pData(es)
  out <- data.frame(ID = pd$geo_accession, title = as.character(pd$title),
                    source = as.character(pd$source_name_ch1), stringsAsFactors = FALSE)
  ch_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
  for (cc in ch_cols) {
    vals <- as.character(pd[[cc]])
    for (i in seq_along(vals)) {
      parts <- unlist(strsplit(vals[i], ";\\s*"))
      for (p in parts) {
        if (grepl(":", p, fixed = TRUE)) {
          kv <- trimws(unlist(strsplit(p, ":", fixed = TRUE)))
          k <- tolower(gsub("[^a-zA-Z0-9]+", "_", kv[1]))
          v <- paste(kv[-1], collapse = ":")
          if (!k %in% colnames(out)) out[[k]] <- NA_character_
          out[[k]][i] <- trimws(v)
        }
      }
    }
  }
  out
}

# ---- (a) GSE29609（GPL1708）：39 例 ccRCC，含 OS 与 Fuhrman 分级、分期 ---------
tryCatch({
  log_msg("处理 GSE29609...")
  es <- get_geo_eset("GSE29609")
  kv <- extract_kv(es)
  log_msg("GSE29609 字段：", paste(colnames(kv), collapse = ", "))
  for (cc in colnames(kv)[-c(1:3)]) {
    log_msg("  ", cc, ": ", paste(head(unique(kv[[cc]]), 5), collapse = " | "))
  }
  m <- geo_expr_symbol(es)

  t_col <- grep("month|time|follow|surv", colnames(kv), value = TRUE)
  s_col <- grep("status|death|dead|alive|state|censor|event|progress", colnames(kv), value = TRUE)
  os_mo <- rep(NA_real_, nrow(kv))
  for (tc_name in t_col) {
    v <- num(gsub("[^0-9.]", "", kv[[tc_name]]))
    if (sum(!is.na(v) & v > 0) > sum(!is.na(os_mo) & os_mo > 0)) os_mo <- v
  }
  os_ev <- rep(NA_integer_, nrow(kv))
  for (sc_name in unique(c(s_col, colnames(kv)[-c(1:3)]))) {
    raw_v <- tolower(as.character(kv[[sc_name]]))
    if (any(grepl("dead|death|deceased|died|^1$", raw_v)) && any(grepl("alive|living|censored|^0$", raw_v))) {
      os_ev <- ifelse(grepl("dead|death|deceased|died|^1$", raw_v), 1L,
                      ifelse(grepl("alive|living|censored|^0$", raw_v), 0L, NA_integer_))
      break
    }
  }
  fg_col <- grep("fuhrman|grade", colnames(kv), value = TRUE)[1]
  st_col <- grep("stage|tnm", colnames(kv), value = TRUE)[1]
  clin <- data.frame(
    ID = kv$ID,
    OS.time = os_mo * MONTH2DAY,
    OS = os_ev,
    grade = if (!is.na(fg_col)) kv[[fg_col]] else NA_character_,
    stage = if (!is.na(st_col)) kv[[st_col]] else NA_character_
  )
  save_cohort("GSE29609", m, clin, "validation (OS)", "GPL1708 Agilent", "OS")
}, error = function(e) log_msg("GSE29609 处理警告：", conditionMessage(e)))

# ---- (b) GSE22541（GPL570）：ccRCC 原发与转移队列（含 DFS/随访）----------------
tryCatch({
  log_msg("处理 GSE22541...")
  es <- get_geo_eset("GSE22541", "GPL570")
  raw_pd <- Biobase::pData(es)
  kv <- extract_kv(es)
  log_msg("GSE22541 kv 列：", paste(colnames(kv), collapse = ", "))
  for (cc in colnames(kv)[-c(1:3)]) {
    log_msg("  ", cc, ": ", paste(head(unique(kv[[cc]]), 5), collapse = " | "))
  }
  ch_all <- apply(raw_pd[, grep("^characteristics_ch1|description|title", colnames(raw_pd), value = TRUE), drop = FALSE],
                  1, paste, collapse = " ; ")
  log_msg("GSE22541 样本描述示例：", head(ch_all, 3))

  m <- geo_expr_symbol(es)
  os_mo <- rep(NA_real_, nrow(raw_pd))
  os_ev <- rep(NA_integer_, nrow(raw_pd))
  for (cc in colnames(kv)[-c(1:3)]) {
    v <- kv[[cc]]
    vn <- num(gsub("[^0-9.]", "", v))
    if (sum(!is.na(vn) & vn > 0) >= 15 && max(vn, na.rm = TRUE) > 5) os_mo <- vn
    lv <- tolower(v)
    if (any(grepl("yes|dead|relapse|recur|metast|event|^1$", lv)) && any(grepl("no|alive|free|none|^0$", lv))) {
      os_ev <- ifelse(grepl("yes|dead|relapse|recur|metast|event|^1$", lv), 1L,
                      ifelse(grepl("no|alive|free|none|^0$", lv), 0L, NA_integer_))
    }
  }
  clin <- data.frame(
    ID = raw_pd$geo_accession,
    OS.time = ifelse(os_mo > 250, os_mo, os_mo * MONTH2DAY),
    OS = os_ev,
    title = as.character(raw_pd$title),
    characteristics = ch_all
  )
  save_cohort("GSE22541", m, clin, "validation (DFS/progression)", "GPL570 Affymetrix", "DFS")
}, error = function(e) log_msg("GSE22541 处理警告：", conditionMessage(e)))

# ---- (c) GSE167573（RNA-seq / 芯片 ccRCC 队列，检查是否含 OS/随访）------------
tryCatch({
  log_msg("处理 GSE167573...")
  es <- get_geo_eset("GSE167573")
  raw_pd <- Biobase::pData(es)
  kv <- extract_kv(es)
  log_msg("GSE167573 kv 列：", paste(colnames(kv), collapse = ", "))
  for (cc in colnames(kv)[-c(1:3)]) {
    log_msg("  ", cc, ": ", paste(head(unique(kv[[cc]]), 5), collapse = " | "))
  }
  if (nrow(Biobase::exprs(es)) > 100) {
    m <- geo_expr_symbol(es)
    t_col <- grep("month|time|follow|surv|os|dfs", colnames(kv), value = TRUE)
    s_col <- grep("status|death|dead|alive|state|censor|event", colnames(kv), value = TRUE)
    os_mo <- if (length(t_col)) num(gsub("[^0-9.]", "", kv[[t_col[1]]])) else rep(NA_real_, nrow(kv))
    os_ev <- rep(NA_integer_, nrow(kv))
    if (length(s_col)) {
      lv <- tolower(kv[[s_col[1]]])
      os_ev <- ifelse(grepl("dead|death|deceased|1|yes", lv), 1L,
                      ifelse(grepl("alive|living|0|no", lv), 0L, NA_integer_))
    }
    clin <- data.frame(ID = kv$ID, OS.time = ifelse(os_mo > 250, os_mo, os_mo * MONTH2DAY), OS = os_ev)
    save_cohort("GSE167573", m, clin, "validation", "GEO", "OS")
  }
}, error = function(e) log_msg("GSE167573 提示：", conditionMessage(e)))

# ---- (d) GSE53757（GPL570）：72 对配对 ccRCC 肿瘤与癌旁正常（共 144 例，含 Stage I-IV）----
tryCatch({
  log_msg("处理 GSE53757（72 对配对肿瘤/正常 + Stage I-IV）...")
  es <- get_geo_eset("GSE53757", "GPL570")
  kv <- extract_kv(es)
  log_msg("GSE53757 kv 列：", paste(colnames(kv), collapse = ", "))
  for (cc in colnames(kv)[-c(1:3)]) {
    log_msg("  ", cc, ": ", paste(head(unique(kv[[cc]]), 6), collapse = " | "))
  }
  m <- geo_expr_symbol(es)
  all_txt <- tolower(paste(kv$title, kv$source, apply(kv[, -c(1:3), drop = FALSE], 1, paste, collapse = " ")))
  grp <- ifelse(grepl("normal|adjacent|non-tumor|nontumor", all_txt), "Normal", "Tumor")
  st <- rep(NA_character_, nrow(kv))
  st[grepl("stage\\s*iv|stage\\s*4", all_txt)] <- "IV"
  st[grepl("stage\\s*iii|stage\\s*3", all_txt) & is.na(st)] <- "III"
  st[grepl("stage\\s*ii|stage\\s*2", all_txt) & is.na(st)] <- "II"
  st[grepl("stage\\s*i\\b|stage\\s*1", all_txt) & is.na(st)] <- "I"
  clin <- data.frame(ID = kv$ID, OS.time = NA_real_, OS = NA_integer_, group = grp, stage = st, title = kv$title)
  save_cohort("GSE53757", m, clin, "validation (paired Tumor/Normal & Stage)", "GPL570 Affymetrix", "Diagnostic / Stage")
}, error = function(e) log_msg("GSE53757 处理警告：", conditionMessage(e)))

# ---- (e) GSE73731（GPL570）：265 例 ccRCC 肿瘤（含 Fuhrman G1-G4 与 Stage I-IV）----
tryCatch({
  log_msg("处理 GSE73731（265 例 ccRCC 肿瘤，含 Fuhrman grade 与 stage）...")
  es <- get_geo_eset("GSE73731", "GPL570")
  kv <- extract_kv(es)
  log_msg("GSE73731 kv 列：", paste(colnames(kv), collapse = ", "))
  for (cc in colnames(kv)[-c(1:3)]) {
    log_msg("  ", cc, ": ", paste(head(unique(kv[[cc]]), 6), collapse = " | "))
  }
  m <- geo_expr_symbol(es)
  gr_col <- grep("grade", colnames(kv), value = TRUE)[1]
  st_col <- grep("stage", colnames(kv), value = TRUE)[1]
  gd_col <- grep("gender|sex", colnames(kv), value = TRUE)[1]
  clin <- data.frame(
    ID = kv$ID, OS.time = NA_real_, OS = NA_integer_,
    grade = if (!is.na(gr_col)) toupper(kv[[gr_col]]) else NA_character_,
    stage = if (!is.na(st_col)) toupper(kv[[st_col]]) else NA_character_,
    gender = if (!is.na(gd_col)) tolower(kv[[gd_col]]) else NA_character_
  )
  save_cohort("GSE73731", m, clin, "validation (Grade & Stage & TME)", "GPL570 Affymetrix", "Grade / Stage")
}, error = function(e) log_msg("GSE73731 处理警告：", conditionMessage(e)))

# =============================================================================
# 3. 免疫治疗队列（GSE67501 肾癌 Nivolumab + IMvigor210 + GSE176307）
# =============================================================================
tryCatch({
  log_msg("处理 GSE67501（RCC anti-PD-1 Nivolumab）...")
  es <- get_geo_eset("GSE67501")
  kv <- extract_kv(es)
  log_msg("GSE67501 kv 列：", paste(colnames(kv), collapse = ", "))
  for (cc in colnames(kv)[-c(1:3)]) {
    log_msg("  ", cc, ": ", paste(unique(kv[[cc]]), collapse = " | "))
  }
  m <- geo_expr_symbol(es)
  all_txt <- tolower(apply(kv, 1, paste, collapse = " "))
  resp <- ifelse(grepl("non-responder|nonresponder|nr|sd|pd|no response", all_txt), "SD/PD",
                 ifelse(grepl("responder|cr|pr|complete|partial", all_txt), "CR/PR", NA_character_))
  clin <- data.frame(ID = kv$ID, OS.time = NA_real_, OS = NA_integer_,
                     response = resp, title = kv$title)
  print(table(GSE67501_response = clin$response, useNA = "ifany"))
  save_cohort("GSE67501", m, clin, "immunotherapy (RCC Nivolumab anti-PD-1)", "GPL14951 Illumina", "ORR")
}, error = function(e) log_msg("GSE67501 处理警告：", conditionMessage(e)))

for (ic_name in c("IMvigor210", "GSE176307")) {
  src_rds <- file.path(BLCA_DIR, "data", "processed", paste0(ic_name, ".rds"))
  dst_rds <- file.path(DIR_PRO, paste0(ic_name, ".rds"))
  if (file.exists(src_rds)) {
    if (!file.exists(dst_rds)) file.copy(src_rds, dst_rds)
    obj <- readRDS(dst_rds)
    has_os <- !is.na(obj$clin$OS.time) & !is.na(obj$clin$OS) & obj$clin$OS.time > 0
    COHORTS[[ic_name]] <- data.frame(
      cohort = ic_name,
      role = if (ic_name == "IMvigor210") "immunotherapy (anti-PD-L1)" else "immunotherapy (anti-PD-1/PD-L1)",
      platform = "RNA-seq (log2 TPM+1)", endpoint = "ORR / OS",
      n_samples = ncol(obj$expr), n_genes = nrow(obj$expr),
      n_with_survival = sum(has_os), n_events = sum(obj$clin$OS[has_os] == 1),
      n_responders = sum(obj$clin$response == "CR/PR", na.rm = TRUE)
    )
    log_msg("已加载免疫治疗队列 ", ic_name, "：", ncol(obj$expr), " 例样本，应答者 ",
            sum(obj$clin$response == "CR/PR", na.rm = TRUE), " 例")
  }
}

# =============================================================================
# 4. 预下载单细胞（GSE159115）与空间转录组（GSE210041）原始文件，确保后续步骤零联网
# =============================================================================
tryCatch({
  log_msg("预下载单细胞数据集 GSE159115（ccRCC + 正常肾组织 10x h5 与官方注释）...")
  geo_suppl("GSE159115", "GSE159115_ccRCC_anno.csv.gz", "GSE159115")
  geo_suppl("GSE159115", "GSE159115_normal_anno.csv.gz", "GSE159115")
  es_sc <- GEOquery::getGEO("GSE159115", GSEMatrix = TRUE, getGPL = FALSE, destdir = file.path(DIR_RAW, "GEO"))[[1]]
  pd_sc <- Biobase::pData(es_sc)
  saveRDS(pd_sc, file.path(DIR_RAW, "GEO", "GSE159115", "GSE159115_pData.rds"))
  supp_col <- grep("^supplementary_file", colnames(pd_sc), value = TRUE)[1]
  for (i in seq_len(nrow(pd_sc))) {
    gsm <- pd_sc$geo_accession[i]
    surl <- as.character(pd_sc[[supp_col]][i])
    if (!is.na(surl) && grepl("\\.h5$", surl)) {
      fn <- basename(surl)
      gsm_suppl(gsm, fn, "GSE159115")
    }
  }
  log_msg("GSE159115 单细胞文件预下载完成")
}, error = function(e) log_msg("GSE159115 预下载提示：", conditionMessage(e)))

tryCatch({
  log_msg("预下载空间转录组数据集 GSE210041（10x Visium ccRCC）...")
  f_tar <- geo_suppl("GSE210041", "GSE210041_RAW.tar", "GSE210041")
  ex_dir <- file.path(DIR_RAW, "GEO", "GSE210041", "extracted")
  if (!dir.exists(ex_dir) || length(list.files(ex_dir)) == 0) {
    dir.create(ex_dir, recursive = TRUE, showWarnings = FALSE)
    utils::untar(f_tar, exdir = ex_dir)
  }
  log_msg("GSE210041 空间转录组文件预下载并解压完成，包含文件：",
          paste(head(list.files(ex_dir), 10), collapse = ", "))
}, error = function(e) log_msg("GSE210041 预下载提示：", conditionMessage(e)))

# =============================================================================
# 汇总输出 Table S2
# =============================================================================
S2 <- do.call(rbind, COHORTS)
print(S2, row.names = FALSE)
utils::write.csv(S2, file.path(OUT, "Table_S2_cohorts.csv"), row.names = FALSE)
save_session(STEP)
log_msg("02_download_bulk 全部完成")
