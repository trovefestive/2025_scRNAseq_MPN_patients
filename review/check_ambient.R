suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("stage8/annotated_cellcycle.rds")
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
genes <- c("PRAME","ITGB3","THBS1","EGF","DKK1","PF4","PPBP","GP9")
k <- c("T_CD4_naive","T_CD8_GZMK","T_memory","T_CD4_memory")
i <- obj$cell_type %in% k & obj$disease_stage != "Reactive"
pts <- names(which(table(obj$sample[i]) >= 20)); st <- tapply(obj$disease_stage, obj$sample, `[`, 1)[pts]
cpm <- sapply(pts, function(p) { s <- Matrix::rowSums(cnt[, i & obj$sample == p, drop = FALSE]); 1e6 * s[genes] / sum(s) })
mf3 <- st == "PMF(MF3)"
cat("T cells (control), CPM: median Early vs median MF3, and MF3 patients above Early max\n")
for (g in genes) cat(sprintf("  %-6s Early %7.1f | MF3 %7.1f | above Early max %d/%d\n", g,
  median(cpm[g, !mf3]), median(cpm[g, mf3]), sum(cpm[g, mf3] > max(cpm[g, !mf3])), sum(mf3)))
