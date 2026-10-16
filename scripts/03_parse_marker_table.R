# ============================================================
# Script : 03_parse_marker_table.R
# Purpose: 解析 mmc2.xlsx，提取文献/snRNA/spRNA marker
# Input  : data/mmc2.xlsx
# Output : data/sn_root_markers.rds
# ============================================================

library(readxl)
library(dplyr)

file_path <- "data/mmc2.xlsx"
lit_marker <- read_excel(file_path, sheet = 1, skip = 2)

# ---------- 1. 文献 marker（左半） ----------
lit_s5 <- lit_marker[, c("geneID", "Expected cell type", "Name",
                         "Root", "Early_nodule", "Hypocotyl", "Leaf", "Pod")]
lit_s5 <- lit_s5[!is.na(lit_s5$geneID), ]

# ---------- 2. snRNA marker（中间） ----------
sn_s7 <- lit_marker[, c("Gene ID...19", "Avgerage log2FC...20", "Pct.1...21",
                        "Pct.2...22", "P value...23", "P value adjusted...24",
                        "Cluster...25", "Tissue...26")]
colnames(sn_s7) <- c("geneID", "avg_log2FC", "pct1", "pct2",
                     "p_val", "p_val_adj", "cluster", "tissue")
sn_s7 <- sn_s7[!is.na(sn_s7$cluster), ]

cluster_map <- unique(sn_s7[, c("cluster", "tissue")])

# ---------- 3. spRNA marker（右半） ----------
sp_s8 <- lit_marker[, c("Gene ID...27", "Avgerage log2FC...28", "Pct.1...29",
                        "Pct.2...30", "P value...31", "P value adjusted...32",
                        "Cluster...33", "Tissue...34")]
colnames(sp_s8) <- c("geneID", "avg_log2FC", "pct1", "pct2",
                     "p_val", "p_val_adj", "cluster", "tissue")
sp_s8 <- sp_s8[!is.na(sp_s8$cluster), ]

# ---------- 4. 筛 Root ----------
sn_root <- sn_s7 %>% filter(tissue == "Root")

# ---------- 5. 加 ann1. 前缀 ----------
sn_root <- sn_root %>% mutate(geneID_ann1 = paste0("ann1.", geneID))

saveRDS(sn_root, "data/sn_root_markers.rds")
cat("Root marker:", nrow(sn_root), "\n")
