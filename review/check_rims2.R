suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("stage8/annotated_cellcycle.rds")
x <- GetAssayData(obj, assay = "RNA", layer = "data")["RIMS2", ] > 0
r <- tapply(x, obj$sample, function(z) round(100 * mean(z), 2))
st <- tapply(obj$disease_stage, obj$sample, `[`, 1); sx <- tapply(obj$sex, obj$sample, `[`, 1)
print(data.frame(stage = st, sex = sx, pct_RIMS2_pos = r)[order(-r), ])
t <- obj$cell_type %in% c("T_CD4_naive","T_CD4_memory","T_memory","T_CD8_GZMK")
cat("\n% RIMS2+ among T cells, top samples:\n"); print(head(sort(tapply(x[t], obj$sample[t], function(z) round(100*mean(z),2)), decreasing = TRUE), 5))
