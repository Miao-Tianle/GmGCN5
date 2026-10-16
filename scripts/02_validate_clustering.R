# ============================================================
# Script : 02_validate_clustering.R
# Purpose: 查看作者聚类信息，验证聚类可复现性
# Input  : data/GSE270392_Gm_atlas_Root.rna.seurat.obj.rds
# Output : results/tables/cluster_summary.txt
# ============================================================

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

seurat_obj <- readRDS("data/GSE270392_Gm_atlas_Root.rna.seurat.obj.rds")

# ---------- 1. 聚类列 ----------
colnames(seurat_obj@meta.data)

# ---------- 2. 各聚类分布 ----------
cat("--- SCT_snn_res.0.4 ---\n"); print(table(seurat_obj$SCT_snn_res.0.4))
cat("--- seurat_clusters ---\n"); print(table(seurat_obj$seurat_clusters))
cat("--- sub_clst1 ---\n");       print(table(seurat_obj$sub_clst1))
cat("--- sub_clst1_8 ---\n");     print(table(seurat_obj$sub_clst1_8))
cat("--- celltype ---\n");        print(table(seurat_obj$celltype))

# ---------- 3. 聚类关系 ----------
print(table(SCT_snn_res.0.4 = seurat_obj$SCT_snn_res.0.4,
            seurat_clusters  = seurat_obj$seurat_clusters))
print(table(seurat_clusters = seurat_obj$seurat_clusters,
            celltype        = seurat_obj$celltype))

# ---------- 4. 可视化 ----------
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "SCT_snn_res.0.4",
              label = TRUE, repel = TRUE) + ggtitle("First-round clustering")
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE) + ggtitle("Final clustering (12)")
p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype",
              label = TRUE, repel = TRUE) + ggtitle("Annotated celltype")
p1 | p2 | p3

# ---------- 5. 保存图 ----------
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/02_clustering_overview.pdf", p1 | p2 | p3,
       width = 18, height = 6)
