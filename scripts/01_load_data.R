# ============================================================
# Script : 01_load_data.R
# Purpose: 读取单细胞 Seurat 对象，查看基本信息
# Input  : data/GSE270392_Gm_atlas_Root.rna.seurat.obj.rds
# Output : 无（查看信息）
# ============================================================

library(Seurat)
library(dplyr)

# ---------- 1. 读取对象 ----------
seurat_obj <- readRDS("data/GSE270392_Gm_atlas_Root.rna.seurat.obj.rds")

# ---------- 2. 基本信息 ----------
seurat_obj
Assays(seurat_obj)
DefaultAssay(seurat_obj)

# ---------- 3. 细胞与基因 ----------
ncol(seurat_obj)
nrow(seurat_obj)
head(Cells(seurat_obj))
head(Features(seurat_obj))

# ---------- 4. 降维 ----------
Reductions(seurat_obj)
head(Embeddings(seurat_obj, "umap"))
head(Embeddings(seurat_obj, "pca"))

# ---------- 5. metadata ----------
head(seurat_obj@meta.data)
colnames(seurat_obj@meta.data)

# ---------- 6. 表达矩阵 ----------
GetAssayData(seurat_obj, assay = "SCT", layer = "counts")[1:5, 1:5]
GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[1:5, 1:5]
