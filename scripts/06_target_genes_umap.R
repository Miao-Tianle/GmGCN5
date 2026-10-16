# ============================================================
# Script : 06_target_genes_umap.R
# Purpose: 两个目标基因的 UMAP 表达
# Input  : results/tables/seurat_obj_annotated.rds
# Output : results/figures/06_target_genes_umap.pdf
# ============================================================

library(Seurat); library(patchwork); library(ggplot2)

seurat_obj <- readRDS("results/tables/seurat_obj_annotated.rds")
DefaultAssay(seurat_obj) <- "SCT"

g1 <- "ann1.Glyma.19G160700"
g2 <- "ann1.Glyma.03G158300"

p1 <- FeaturePlot(seurat_obj, features = g1) + ggtitle("Glyma.19G160700")
p2 <- FeaturePlot(seurat_obj, features = g2) + ggtitle("Glyma.03G158300")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/06_target_genes_umap.pdf", p1 | p2, width = 10, height = 5)
