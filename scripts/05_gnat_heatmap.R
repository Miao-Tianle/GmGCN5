# ============================================================
# Script : 05_gnat_heatmap.R
# Purpose: GNAT 家族基因在 12 类细胞的表达热图
# Input  : results/tables/seurat_obj_annotated.rds,
#          data/GNATjaizu.txt
# Output : results/figures/05_GNAT_heatmap.pdf
# ============================================================

library(Seurat)
library(pheatmap)

seurat_obj <- readRDS("results/tables/seurat_obj_annotated.rds")
DefaultAssay(seurat_obj) <- "SCT"

# ---------- 1. 读基因 ----------
genes <- readLines("data/GNATjaizu.txt")
genes <- trimws(genes[genes != ""])
genes_ann1 <- paste0("ann1.", genes)
genes_use  <- intersect(genes_ann1, rownames(seurat_obj))
cat("匹配:", length(genes_use), "/", length(genes), "\n")

# ---------- 2. 每类平均表达 ----------
avg <- AverageExpression(seurat_obj, features = genes_use,
                         group.by = "celltype_named",
                         assays = "SCT", layer = "data")$SCT
rownames(avg) <- sub("^ann1\\.", "", rownames(avg))

# ---------- 3. z-score + 热图 ----------
avg_z <- t(scale(t(avg)))
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
pdf("results/figures/05_GNAT_heatmap.pdf", width = 8, height = 10)
pheatmap(avg_z, cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize_row = 6, fontsize_col = 10, main = "GNAT family genes")
dev.off()
