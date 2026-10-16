# ============================================================
# Script : 08_coexpression_cells.R
# Purpose: 单细胞内共表达细胞计数
# Input  : results/tables/seurat_obj_annotated.rds,
#          results/tables/only_target1_0.75.csv
# Output : results/tables/coexpr_cells_res1.csv
#          results/tables/coexpr_cells_res2.csv
# ============================================================

library(Seurat); library(dplyr)

seurat_obj <- readRDS("results/tables/seurat_obj_annotated.rds")
DefaultAssay(seurat_obj) <- "SCT"
expr_mat <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")

only1 <- read.csv("results/tables/only_target1_0.75.csv")
only2 <- read.csv("results/tables/only_target2_0.75.csv")

target1 <- "ann1.Glyma.19G160700"
target2 <- "ann1.Glyma.03G158300"

count_coexpr <- function(target, candidates) {
  celltypes <- unique(seurat_obj$celltype_named)
  lapply(celltypes, function(ct) {
    cells_ct <- colnames(seurat_obj)[seurat_obj$celltype_named == ct]
    tg <- expr_mat[target, cells_ct]
    n  <- sapply(candidates, function(g) sum(tg > 0 & expr_mat[g, cells_ct] > 0))
    data.frame(celltype = ct, gene = candidates,
               coexpr_cells = n, total_cells = length(cells_ct),
               pct = round(n / length(cells_ct) * 100, 1))
  }) %>% bind_rows()
}

res1 <- count_coexpr(target1, only1$gene)
res2 <- count_coexpr(target2, only2$gene)

write.csv(res1, "results/tables/coexpr_cells_res1.csv", row.names = FALSE)
write.csv(res2, "results/tables/coexpr_cells_res2.csv", row.names = FALSE)
