suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("stage8/annotated_cellcycle.rds")
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
genes <- c("PRAME","AXL","CH25H","THBS1","LUZP2","ITGB3","EGF","DKK1","DOK5","SEMA3A","DNASE1L3")
for (k in c("HSC_MPP","MEP")) {
  cat("\n==", k, ": CPM per patient (>=20 cells), stages R/P/M2/M3 ==\n")
  i <- obj$cell_type == k & obj$disease_stage != "Reactive"
  pts <- names(which(table(obj$sample[i]) >= 20))
  st <- tapply(obj$disease_stage, obj$sample, `[`, 1)[pts]
  cpm <- sapply(pts, function(p) { s <- Matrix::rowSums(cnt[, i & obj$sample == p, drop = FALSE]); 1e6 * s[genes] / sum(s) })
  ord <- order(st == "PMF(MF3)", pts)
  m <- round(cpm[, ord], 1); colnames(m) <- paste0(substr(pts[ord], 9, 11), ifelse(st[ord] == "PMF(MF3)", "*", ""))
  print(m)
  mf3 <- st == "PMF(MF3)"
  for (g in genes) cat(sprintf("  %-9s MF3 patients above Early max: %d/%d\n", g,
                               sum(cpm[g, mf3] > max(cpm[g, !mf3])), sum(mf3)))
}
