# Review check: is cluster 22 (Megakaryocyte_mature) real MKs, or MK-endothelial/stromal doublets?
suppressPackageStartupMessages({ library(Seurat) })
obj <- readRDS("stage8/annotated_cellcycle.rds")
ct  <- obj$cell_type
e   <- GetAssayData(obj, assay = "RNA", layer = "data")
foc <- c("Megakaryocyte_mature","Megakaryocyte","Megakaryocyte_prog","MSC_stromal","HSC_MPP","Monocyte_classical")

cat("== QC and Scrublet score by cell type ==\n")
q <- do.call(rbind, lapply(foc, function(k) { i <- ct == k
  data.frame(cell_type = k, n = sum(i), med_nFeature = median(obj$nFeature_RNA[i]),
             med_nCount = median(obj$nCount_RNA[i]), med_mt = round(median(obj$percent.mt[i]), 2),
             med_doublet = round(median(obj$doublet_score[i]), 3),
             q90_doublet = round(quantile(obj$doublet_score[i], 0.9), 3)) }))
print(q, row.names = FALSE)

genes <- list(
  MK_specific   = c("ITGA2B","GP9","GP1BA","PF4","PPBP","TUBB1","GP6"),
  Endothelial   = c("CDH5","EMCN","KDR","PLVAP","CLDN5","FLT1","ACKR1","TEK"),
  Stromal       = c("COL1A1","COL3A1","LEPR","CXCL12","DCN","PDGFRB"),
  Cluster22_sig = c("SERPINE1","VEGFC","PDGFB","PDGFA","CXCL1","THBS1","PDE3A"))
cat("\n== % of cells expressing (>0) ==\n")
for (g in names(genes)) {
  gg <- intersect(genes[[g]], rownames(e)); cat("--", g, "\n")
  m <- sapply(foc, function(k) round(100 * rowMeans(e[gg, ct == k, drop = FALSE] > 0), 1))
  colnames(m) <- c("MKmat","MK","MKprog","MSC","HSC","Mono"); print(m)
}

cat("\n== Co-expression inside cluster 22 (per cell) ==\n")
i22 <- ct == "Megakaryocyte_mature"
mk  <- colSums(e[intersect(c("ITGA2B","GP9","PF4","PPBP"), rownames(e)), i22] > 0) >= 2
endo<- colSums(e[intersect(c("CDH5","EMCN","KDR","PLVAP","CLDN5"), rownames(e)), i22] > 0) >= 2
str <- colSums(e[intersect(c("COL1A1","COL3A1","DCN","LEPR"), rownames(e)), i22] > 0) >= 2
cat(sprintf("MK program (>=2 of ITGA2B/GP9/PF4/PPBP): %.1f%%\n", 100*mean(mk)))
cat(sprintf("Endothelial program (>=2 of CDH5/EMCN/KDR/PLVAP/CLDN5): %.1f%%\n", 100*mean(endo)))
cat(sprintf("Stromal program (>=2 of COL1A1/COL3A1/DCN/LEPR): %.1f%%\n", 100*mean(str)))
cat(sprintf("MK AND endothelial: %.1f%% | MK AND stromal: %.1f%%\n", 100*mean(mk & endo), 100*mean(mk & str)))
sig <- colSums(e[intersect(c("SERPINE1","VEGFC","PDGFB","CXCL1"), rownames(e)), i22] > 0) >= 2
cat(sprintf("Cluster-22 signature (>=2 of SERPINE1/VEGFC/PDGFB/CXCL1): %.1f%% | of those, MK program: %.1f%% | endothelial: %.1f%%\n",
            100*mean(sig), 100*mean(mk[sig]), 100*mean(endo[sig])))
