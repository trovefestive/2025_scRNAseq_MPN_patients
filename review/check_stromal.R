suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("stage8/annotated_cellcycle.rds")
e <- GetAssayData(obj, assay = "RNA", layer = "data")
strom <- colSums(e[intersect(c("COL1A1","COL1A2","COL3A1","DCN","CXCL12","PDGFRB"), rownames(e)), ] > 0) >= 3
cat("Cells with stromal programme (>=3 of COL1A1/COL1A2/COL3A1/DCN/CXCL12/PDGFRB):", sum(strom), "\n")
cat("  inside MSC_stromal (res 0.8 cluster 25):", sum(strom & obj$cell_type == "MSC_stromal"),
    "| elsewhere:", sum(strom & obj$cell_type != "MSC_stromal"), "\n")
print(head(sort(table(obj$cell_type[strom & obj$cell_type != "MSC_stromal"]), decreasing = TRUE), 6))
for (r in c("RNA_snn_res.0.8","RNA_snn_res.1","RNA_snn_res.1.2")) {
  cl <- obj[[r]][, 1]; tb <- table(cl[strom]); best <- names(which.max(tb))
  cat(sprintf("%s: best stromal cluster %s has %d cells, %d with programme (%.0f%%)\n", r, best,
      sum(cl == best), max(tb), 100 * max(tb) / sum(cl == best)))
}
cat("\nStromal-programme cells by stage:\n"); print(table(obj$disease_stage[strom]))
