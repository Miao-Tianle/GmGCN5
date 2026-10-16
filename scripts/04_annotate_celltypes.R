# ============================================================
# Script : 04_annotate_celltypes.R
# Purpose: 用 Root marker 注释 12 类 RNA_x 细胞
# Input  : data/GSE270392_...rds, data/sn_root_markers.rds
# Output : results/tables/seurat_obj_annotated.rds
# ============================================================

library(Seurat)
library(dplyr)
library(pheatmap)

seurat_obj <- readRDS("data/GSE270392_Gm_atlas_Root.rna.seurat.obj.rds")
sn_root    <- readRDS("data/sn_root_markers.rds")

# ---------- 1. 匹配 marker ----------
markers_use <- sn_root %>% filter(geneID_ann1 %in% rownames(seurat_obj))

# ---------- 2. 每类平均表达 ----------
avg <- AverageExpression(seurat_obj,
                         features = unique(markers_use$geneID_ann1),
                         group.by = "celltype",
                         assays = "SCT", layer = "data")$SCT

# ---------- 3. 打分矩阵 ----------
score <- sapply(unique(markers_use$cluster), function(ct) {
  genes <- intersect(markers_use$geneID_ann1[markers_use$cluster == ct],
                     rownames(avg))
  if (length(genes) == 0) return(rep(NA, ncol(avg)))
  Matrix::colMeans(avg[genes, , drop = FALSE])
})
rownames(score) <- colnames(avg)
colnames(score) <- unique(markers_use$cluster)
print(round(score, 3))

# ---------- 4. 最优匹配 ----------
best <- colnames(score)[apply(score, 1, which.max)]
result <- data.frame(RNA_cluster = rownames(score),
                     best_root_celltype = best,
                     score = round(apply(score, 1, max), 3))
print(result)

# ---------- 5. 映射 ----------
map <- c(
  "RNA_0" = "Cortex_0", "RNA_1_0" = "Procambium_1_0",
  "RNA_1_1" = "Pericycle_1_1", "RNA_2" = "Epidermis_2",
  "RNA_3" = "Epidermis_3", "RNA_4" = "Endodermis_4",
  "RNA_5" = "Cortex_5", "RNA_6" = "Root_hair_6",
  "RNA_7" = "Cortex_7", "RNA_8_0" = "Phloem_8_0",
  "RNA_8_1" = "Xylem_8_1", "RNA_8_2" = "Phloem_8_2"
)
seurat_obj$celltype_named <- unname(map[as.character(seurat_obj$celltype)])

# ---------- 6. 验证 + 保存 ----------
table(seurat_obj$celltype_named, useNA = "ifany")
sum(is.na(seurat_obj$celltype_named))

saveRDS(seurat_obj, "results/tables/seurat_obj_annotated.rds")
