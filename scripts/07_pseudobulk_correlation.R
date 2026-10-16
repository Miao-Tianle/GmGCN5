# ============================================================
# Script : 07_pseudobulk_correlation.R
# Purpose: 伪 bulk 斯皮尔曼相关，筛 |r| > 0.75
# Input  : results/tables/seurat_obj_annotated.rds
# Output : results/tables/only_target1_0.75.csv
#          results/tables/only_target2_0.75.csv
# ============================================================

library(Seurat); library(dplyr)

seurat_obj <- readRDS("results/tables/seurat_obj_annotated.rds")
DefaultAssay(seurat_obj) <- "SCT"

# ---------- 1. 全基因组伪 bulk ----------
avg_all <- AverageExpression(seurat_obj, group.by = "celltype_named",
                             assays = "SCT", layer = "data")$SCT
pb_mat  <- as.matrix(avg_all)
pb_filt <- pb_mat[Matrix::rowMeans(pb_mat > 0) >= 0.5, ]
cat("过滤后:", nrow(pb_filt), "\n")

# ---------- 2. 目标基因 ----------
target1 <- "ann1.Glyma.19G160700"
target2 <- "ann1.Glyma.03G158300"

# ---------- 3. 相关 ----------
cor1 <- apply(pb_filt, 1, function(x) cor(x, pb_filt[target1, ], method = "spearman"))
cor2 <- apply(pb_filt, 1, function(x) cor(x, pb_filt[target2, ], method = "spearman"))

df1 <- data.frame(gene = names(cor1), cor_target1 = as.numeric(cor1))
df2 <- data.frame(gene = names(cor2), cor_target2 = as.numeric(cor2))
merged <- inner_join(df1, df2, by = "gene") %>%
  filter(gene != target1, gene != target2)

# ---------- 4. 筛 >0.75 ----------
thr <- 0.75
only1 <- merged %>% filter(abs(cor_target1) > thr, abs(cor_target2) <= thr)
only2 <- merged %>% filter(abs(cor_target2) > thr, abs(cor_target1) <= thr)
both  <- merged %>% filter(abs(cor_target1) > thr, abs(cor_target2) > thr)

cat("only1:", nrow(only1), "| only2:", nrow(only2), "| both:", nrow(both), "\n")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(only1, "results/tables/only_target1_0.75.csv", row.names = FALSE)
write.csv(only2, "results/tables/only_target2_0.75.csv", row.names = FALSE)
write.csv(both,  "results/tables/both_0.75.csv",         row.names = FALSE)
