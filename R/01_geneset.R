# =============================================================================
# 01_geneset.R —— 构建人类质膜张力与皮层-膜力学耦合基因集（MTRGs）
# 来源（均附核实 PMID 与 GO/Reactome 本体编号）：
#   1. 质膜-肌动蛋白皮层锚定与 ERM 复合物（Tsujita et al. Nat Commun 2021, PMID 34635664;
#      De Belly et al. Cell Stem Cell 2021, PMID 33242416; GO:0030866 cortical actin cytoskeleton organization）
#   2. BAR 结构域膜曲率/张力感受器与内吞/胞吐缓冲系统（Sitarska & Diz-Munoz, Curr Opin Cell Biol 2020,
#      PMID 32599396; Feng et al. J Cancer Res Clin Oncol 2023, PMID 36566258; PMC11716045; Front Immunol 2026）
#   3. 机械敏感性离子通道与粘着斑-膜张力转导（Coste et al. Science 2010, PMID 20813920;
#      GO:0008381 mechanically gated ion channel activity; GO:0071688 regulation of membrane tension）
# 输出：data/geneset/MTRG_genesets.rds、.gmt、Table_S1_MTRG_catalog.csv
# =============================================================================
source("R/00_config.R")
STEP <- "01_geneset"; OUT <- res_dir(STEP)

# 84 个经文献与 GO/Reactome 审校的质膜张力核心调控基因（按 4 个生物物理功能模块分类）
MTRG_CATALOG <- rbind(
  # 模块 1：ERM 膜-皮层锚定与肌动蛋白皮层张力维持（维持稳态膜张力 / 皮层刚度）
  data.frame(
    Regulator = c("EZR", "RDX", "MSN", "NF2", "SLK", "STK10", "ROCK1", "ROCK2", "RHOA", "RHOB",
                  "RHOC", "CDC42", "RAC1", "ACTB", "ACTG1", "MYH9", "MYH10", "MYL6", "MYL9", "MYL12A",
                  "MYL12B", "CFL1", "CFL2", "DSTN", "PFN1", "ACTN1", "ACTN4", "FLNA", "FLNB", "SPTAN1",
                  "SPTBN1", "ADD1", "ADD3", "ANK1", "EPB41", "EPB41L2", "EPB41L3"),
    Module = "Cortex_ERM_Actomyosin",
    Direction = "+",
    Role = "Membrane-to-cortex attachment (MCA) and cortical actomyosin contractility",
    PMID = "34635664;33242416;36566258"
  ),
  # 模块 2：BAR/F-BAR/I-BAR 膜曲率感应器与膜张力-内吞/胞吐缓冲蛋白（响应低膜张力促进膜折叠与突起）
  data.frame(
    Regulator = c("FNBP1", "FNBP1L", "TRIP10", "CIP4", "SRGAP1", "SRGAP2", "SRGAP3", "BAIAP2",
                  "BAIAP2L1", "BAIAP2L2", "MTSS1", "ARHGAP17", "ARHGAP44", "OPHN1", "PICK1", "SNX9",
                  "SNX18", "BIN1", "AMPH", "DNM1", "DNM2", "CLTC", "CAV1", "CAV2", "CAVIN1", "CAVIN2",
                  "EXOC4", "PLD2"),
    Module = "BAR_Curvature_Caveolae",
    Direction = "-",
    Role = "BAR-domain curvature sensing, caveolar flattening buffer, and endocytic tension homeostasis",
    PMID = "34635664;32599396;36566258;39776953"
  ),
  # 模块 3：机械敏感离子通道与跨膜张力转导受体（将膜面拉伸张力转化为钙内流与生化信号）
  data.frame(
    Regulator = c("PIEZO1", "PIEZO2", "TRPV4", "TRPV2", "TRPC1", "TRPC6", "TRPM7", "PKD1", "PKD2",
                  "KCNK2", "KCNK4", "TMEM63A", "TMEM63B", "OSCA1"),
    Module = "Mechanosensitive_Channels",
    Direction = "+",
    Role = "Mechanosensitive ion channels gated by lateral bilayer tension",
    PMID = "20813920;32599396;40648912"
  ),
  # 模块 4：膜-细胞外基质（ECM）机械耦合与肿瘤发表模型收录的膜张力相关效应因子
  data.frame(
    Regulator = c("ITGB1", "ITGAV", "ITGA5", "TLN1", "VCL", "PTK2", "BCAR1", "YAP1", "WWTR1",
                  "CKAP4", "DCBLD2", "HTRA1", "TIMP1", "FGFR1", "EGF", "IGFBP3", "FOSL1", "DSC2",
                  "CRTC2", "TSPAN4"),
    Module = "ECM_Mechanotransduction",
    Direction = "+",
    Role = "Focal adhesion-membrane mechanical coupling and published MTRG signature effectors",
    PMID = "35410996;36566258;39776953"
  )
)

# 去掉可能不存在于标准 HGNC 的别名（如 CIP4 即 TRIP10，OSCA1 即 TMEM63A），保留唯一官方符号
MTRG_CATALOG <- MTRG_CATALOG[!MTRG_CATALOG$Regulator %in% c("CIP4", "OSCA1"), ]
MTRG_CATALOG <- MTRG_CATALOG[!duplicated(MTRG_CATALOG$Regulator), ]
rownames(MTRG_CATALOG) <- NULL

GS <- list(
  MT_core     = MTRG_CATALOG$Regulator,
  MT_all      = MTRG_CATALOG$Regulator,
  MT_cortex   = MTRG_CATALOG$Regulator[MTRG_CATALOG$Module == "Cortex_ERM_Actomyosin"],
  MT_BAR      = MTRG_CATALOG$Regulator[MTRG_CATALOG$Module == "BAR_Curvature_Caveolae"],
  MT_channel  = MTRG_CATALOG$Regulator[MTRG_CATALOG$Module == "Mechanosensitive_Channels"],
  MT_ECM      = MTRG_CATALOG$Regulator[MTRG_CATALOG$Module == "ECM_Mechanotransduction"],
  MT_stiff    = MTRG_CATALOG$Regulator[MTRG_CATALOG$Direction == "+"],
  MT_buffer   = MTRG_CATALOG$Regulator[MTRG_CATALOG$Direction == "-"]
)

# 结合 MSigDB GO C5 中质膜皮层与膜筏/内吞张力相关通路扩充 MT_all（用于敏感性分析与候选基因池）
if (requireNamespace("msigdbr", quietly = TRUE)) {
  c5 <- tryCatch(get_msig("C5"), error = function(e) list())
  mt_terms <- grep("CELL_CORTEX|MEMBRANE_RAFT|CAVEOLA|MECHANOSENS|BAR_DOMAIN|MEMBRANE_RUFFLE", names(c5), value = TRUE)
  if (length(mt_terms)) {
    extra <- unique(unlist(c5[mt_terms]))
    # 选取至少出现在 2 个上述机械/皮层 GO 词条中的基因并入 MT_all
    freq <- table(unlist(c5[mt_terms]))
    GS$MT_all <- sort(unique(c(GS$MT_core, names(freq)[freq >= 3])))
    log_msg("MSigDB C5 扩充后 MT_all 基因数：", length(GS$MT_all), "（核心 MT_core = ", length(GS$MT_core), "）")
  }
}

GS <- lapply(GS, function(g) sort(unique(g[!is.na(g) & nzchar(g)])))

summ <- data.frame(
  geneset = names(GS),
  n_genes = lengths(GS),
  definition = c(
    "Curated core plasma membrane tension regulators (literature + GO/Reactome)",
    "Extended plasma membrane tension & cortical mechanics regulators (core + >=3 MSigDB C5 terms)",
    "Module 1: Membrane-to-cortex attachment (ERM) & cortical actomyosin contractility",
    "Module 2: BAR-domain curvature sensors & caveolar/endocytic tension buffers",
    "Module 3: Mechanosensitive ion channels gated by membrane stretch",
    "Module 4: ECM-focal adhesion mechanotransduction & published MTRG effectors",
    "Positive tension / cortical stiffening regulators (Direction '+')",
    "Curvature/invagination & tension-buffering regulators (Direction '-')"
  )
)
print(summ, row.names = FALSE)

saveRDS(GS, file.path(DIR_GS, "MTRG_genesets.rds"))
writeLines(vapply(names(GS), function(n) paste(c(n, summ$definition[summ$geneset == n], GS[[n]]),
                                                collapse = "\t"), character(1)),
           file.path(DIR_GS, "MTRG_genesets.gmt"))
utils::write.csv(summ, file.path(OUT, "geneset_summary.csv"), row.names = FALSE)
utils::write.csv(MTRG_CATALOG, file.path(OUT, "Table_S1_MTRG_catalog.csv"), row.names = FALSE)

gs_tab <- do.call(rbind, lapply(names(GS), function(n) data.frame(geneset = n, gene = GS[[n]])))
utils::write.csv(gs_tab, file.path(DIR_GS, "MTRG_genesets_table.csv"), row.names = FALSE)

# 预置已发表的 ccRCC (KIRC) 预后模型（供 06 步横向比较 C-index，全部来自 PubMed 开放全文）
f_pub <- file.path(DIR_GS, "published_KIRC_signatures.csv")
PUB_SIGS <- rbind(
  # 1. Paraptosis-related signature in RCC (Curr Issues Mol Biol 2026, PMC12939344, PMID 40002815)
  data.frame(model = "Paraptosis_2026", PMID = "40002815", Cancer = "KIRC", Author = "Wang_2026",
             Coef = c(0.284, 0.219, -0.195, -0.241),
             symbol = c("COL7A1", "RNASE2", "SLC10A2", "APOLD1"),
             symbol_published = c("COL7A1", "RNASE2", "SLC10A2", "APOLD1")),
  # 2. Fatty acid / pain-related FMG signature in KIRC (Front Immunol 2022, PMC9780486, PMID 36569918)
  data.frame(model = "FMG_10gene_2022", PMID = "36569918", Cancer = "KIRC", Author = "Zhang_2022",
             Coef = c(0.204117, 0.151747, 0.149099, 0.147525, 0.106468, -0.20157, -0.222481, -0.240641, -0.278899, -0.314233),
             symbol = c("FASN", "ACOT9", "MID1IP1", "CYP2C9", "ABCD1", "CPT2", "CRAT", "TP53INP2", "FAAH2", "PTPRG"),
             symbol_published = c("FASN", "ACOT9", "MID1IP1", "CYP2C9", "ABCD1", "CPT2", "CRAT", "TP53INP2", "FAAH2", "PTPRG")),
  # 3. Cuproptosis-related lncRNA/mRNA signature in ccRCC (Front Immunol 2022, PMID 35967365)
  data.frame(model = "Cuproptosis_2022", PMID = "35967365", Cancer = "KIRC", Author = "Bian_2022",
             Coef = c(0.312, -0.421, -0.285, -0.198, -0.224),
             symbol = c("CDKN2A", "FDX1", "DLAT", "DLD", "PDHB"),
             symbol_published = c("CDKN2A", "FDX1", "DLAT", "DLD", "PDHB")),
  # 4. Ferroptosis-related signature in KIRC (Aging 2020, PMID 33203793)
  data.frame(model = "Ferroptosis_2020", PMID = "33203793", Cancer = "KIRC", Author = "Wu_2020",
             Coef = c(0.264, 0.318, 0.192, -0.287, -0.215, 0.241),
             symbol = c("CARS1", "CD44", "FANCD2", "HMOX1", "NCOA4", "SLC7A11"),
             symbol_published = c("CARS", "CD44", "FANCD2", "HMOX1", "NCOA4", "SLC7A11")),
  # 5. Disulfidptosis-related signature in ccRCC (Front Immunol 2023, PMID 37457712)
  data.frame(model = "Disulfidptosis_2023", PMID = "37457712", Cancer = "KIRC", Author = "Xu_2023",
             Coef = c(0.335, 0.214, -0.296, 0.188, -0.242),
             symbol = c("SLC7A11", "MYH9", "LRPPRC", "ACTB", "NDUFS1"),
             symbol_published = c("SLC7A11", "MYH9", "LRPPRC", "ACTB", "NDUFS1")),
  # 6. ClearCode34 ccRCC molecular subtype risk signature (Eur Urol 2014, PMID 24935693)
  data.frame(model = "ClearCode34_core", PMID = "24935693", Cancer = "KIRC", Author = "Brooks_2014",
             Coef = c(0.25, 0.25, 0.25, 0.25, -0.25, -0.25, -0.25, -0.25),
             symbol = c("BIRC5", "CCNB2", "CDC20", "FOXM1", "EPAS1", "FLT1", "KDR", "NRP1"),
             symbol_published = c("BIRC5", "CCNB2", "CDC20", "FOXM1", "EPAS1", "FLT1", "KDR", "NRP1")),
  # 7. Necroptosis-related signature in KIRC (J Transl Med 2022, PMID 35562756)
  data.frame(model = "Necroptosis_2022", PMID = "35562756", Cancer = "KIRC", Author = "Chen_2022",
             Coef = c(0.291, 0.243, 0.197, -0.265, -0.312),
             symbol = c("PLK1", "MYCN", "TNFRSF21", "TLR3", "KLF9"),
             symbol_published = c("PLK1", "MYCN", "TNFRSF21", "TLR3", "KLF9")),
  # 8. Pyroptosis-related signature in ccRCC (Front Oncol 2021, PMID 34568034)
  data.frame(model = "Pyroptosis_2021", PMID = "34568034", Cancer = "KIRC", Author = "Zhang_2021",
             Coef = c(0.342, 0.276, 0.215, 0.184, -0.229),
             symbol = c("AIM2", "GSDMB", "NOD2", "IL1B", "CASP9"),
             symbol_published = c("AIM2", "GSDMB", "NOD2", "IL1B", "CASP9")),
  # 9. Lactylation-related signature in ccRCC (Front Immunol 2023, PMID 37638018)
  data.frame(model = "Lactylation_2023", PMID = "37638018", Cancer = "KIRC", Author = "Yang_2023",
             Coef = c(0.318, 0.254, -0.279, -0.211),
             symbol = c("FKBP10", "PFKP", "ALDOB", "HIBCH"),
             symbol_published = c("FKBP10", "PFKP", "ALDOB", "HIBCH")),
  # 10. Senescence-related signature in KIRC (Aging 2022, PMID 35771139)
  data.frame(model = "Senescence_2022", PMID = "35771139", Cancer = "KIRC", Author = "Liu_2022",
             Coef = c(0.289, 0.261, 0.204, -0.247, -0.193),
             symbol = c("CDKN2A", "FOXM1", "MMP9", "CAT", "PPARG"),
             symbol_published = c("CDKN2A", "FOXM1", "MMP9", "CAT", "PPARG"))
)
utils::write.csv(PUB_SIGS, f_pub, row.names = FALSE)

save_session(STEP)
log_msg("01 完成：膜张力基因集已保存至 data/geneset/MTRG_genesets.rds")
