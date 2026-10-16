# ============================================================
# Script : 09_enrichment_filter.R
# Purpose: 共表达富集倍数，筛最终候选
# Input  : results/tables/seurat_obj_annotated.rds,
#          results/tables/only_target1_0.75.csv
# Output : results/tables/final_target1_coexpr5.csv
#          results/tables/final_target2_coexpr5.csv
# ============================================================

library(Seurat); library(dplyr)

seurat_obj <- readRDS("results/tables/seurat_obj_annotated.rds")
DefaultAssay(seurat_obj) <- "SCT"
expr_mat <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")

only1 <- read.csv("results/tables/only_target1_0.75.csv")
only2 <- read.csv("results/tables/only_target2_0.75.csv")
target1 <- "ann1.Glyma.19G160700"
target2 <- "ann1.Glyma.03G158300"

enrich_coexpr <- function(target, candidates) {
  celltypes <- unique(seurat_obj$celltype_named)
  lapply(celltypes, function(ct) {
    cells_ct <- colnames(seurat_obj)[seurat_obj$celltype_named == ct]
    tg  <- expr_mat[target, cells_ct] > 0
    n   <- length(cells_ct)
    ptg <- mean(tg)
    sapply(candidates, function(g) {
      gpos <- expr_mat[g, cells_ct] > 0
      pg   <- mean(gpos)
      obs  <- sum(tg & gpos)
      exp  <- ptg * pg * n
      data.frame(celltype = ct, gene = g, coexpr = obs,
                 expected = round(exp, 2),
                 enrich = round(obs / max(exp, 0.01), 2),
                 pct = round(obs / n * 100, 1))
    }, simplify = FALSE) %>% bind_rows()
  }) %>% bind_rows()
}

res1_enrich <- enrich_coexpr(target1, only1$gene)
res2_enrich <- enrich_coexpr(target2, only2$gene)

final1 <- res1_enrich %>% filter(coexpr >= 5, enrich > 2) %>%
  group_by(gene) %>% slice_max(enrich, n = 1) %>% ungroup() %>% arrange(desc(enrich))
final2 <- res2_enrich %>% filter(coexpr >= 5, enrich > 2) %>%
  group_by(gene) %>% slice_max(enrich, n = 1) %>% ungroup() %>% arrange(desc(enrich))

cat("target1:", nrow(final1), "| target2:", nrow(final2), "\n")

write.csv(final1, "results/tables/final_target1_coexpr5.csv", row.names = FALSE)
write.csv(final2, "results/tables/final_target2_coexpr5.csv", row.names = FALSE)
