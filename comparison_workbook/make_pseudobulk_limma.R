#!/usr/bin/env Rscript
# Whole-marrow pseudobulk per patient + limma-voom, laid out like the
# "Rampal MPN Analysis.xlsx" workbook (limma topTable columns) for side-by-side use.
# Counts from all QC-passed cells of each patient are summed into one bulk-like
# profile. Design ~ 0 + group + sex (MF3 is 8M/1F, so sex is adjusted).
# SRR27748964 (erythroid-only, ~76% late erythroid) is kept in the expression
# matrix but excluded from the model, as in CellChat (Stage 13).
suppressPackageStartupMessages({ library(Seurat); library(edgeR); library(limma) })
setwd(Sys.getenv("PROJ_DIR", getwd()))
od <- "comparison_workbook"
obj <- readRDS("stage8/annotated_cellcycle.rds")
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
smp <- sort(unique(obj$sample))
pb  <- sapply(smp, function(s) Matrix::rowSums(cnt[, obj$sample == s, drop = FALSE]))
grp_map <- c(Reactive = "Reactive", prefibrotic = "PrePMF", "PMF(MF2)" = "MF2", "PMF(MF3)" = "MF3")
first <- function(v) tapply(v, obj$sample, `[`, 1)[smp]
mk <- obj$cell_type %in% c("Megakaryocyte","Megakaryocyte_prog","Megakaryocyte_mature")
sra <- read.csv("SraRunTable.csv", check.names = FALSE)
meta <- data.frame(sample = smp, patient = first(obj$patient), disease_stage = first(obj$disease_stage),
  group = unname(grp_map[first(obj$disease_stage)]), sex = first(obj$sex), age = first(obj$age),
  collection_date = sra$Collection_Date[match(smp, sra$Run)],
  n_cells = as.integer(table(obj$sample)[smp]),
  pct_megakaryocyte_lineage = round(100 * tapply(mk, obj$sample, mean)[smp], 2),
  qc_flag = first(obj$qc_flag), in_model = smp != "SRR27748964", row.names = NULL)
rm(obj, cnt); invisible(gc())

y <- DGEList(pb, samples = meta)
mm <- meta[meta$in_model, ]
group <- factor(mm$group, levels = c("Reactive","PrePMF","MF2","MF3")); sex <- factor(mm$sex)
design <- model.matrix(~ 0 + group + sex); colnames(design) <- sub("^group", "", colnames(design))
ym <- y[, meta$in_model]
keep <- filterByExpr(ym, design = design)
ym <- calcNormFactors(ym[keep, , keep.lib.sizes = FALSE])
v <- voom(ym, design)
fit <- lmFit(v, design)
n <- table(group)
w_early <- sprintf("(%d*PrePMF + %d*MF2)/%d", n["PrePMF"], n["MF2"], n["PrePMF"] + n["MF2"])
w_mpn <- sprintf("(%d*PrePMF + %d*MF2 + %d*MF3)/%d", n["PrePMF"], n["MF2"], n["MF3"], sum(n[-1]))
cons <- c(PrePMF_vs_Reactive = "PrePMF - Reactive", MF2_vs_Reactive = "MF2 - Reactive",
          MF3_vs_Reactive = "MF3 - Reactive", MF3_vs_PrePMF = "MF3 - PrePMF",
          MF3_vs_Early = paste0("MF3 - ", w_early), MPNvsReactive = paste0(w_mpn, " - Reactive"))
cm <- makeContrasts(contrasts = cons, levels = design); colnames(cm) <- names(cons)
fit2 <- eBayes(contrasts.fit(fit, cm))
chrom <- read.delim("genesets/gene_chrom.tsv", header = FALSE, col.names = c("chr","gene"))
for (k in names(cons)) {
  tt <- topTable(fit2, coef = k, number = Inf, sort.by = "P")
  tt <- data.frame(symbol = rownames(tt), chromosome = chrom$chr[match(rownames(tt), chrom$gene)],
                   tt[, c("logFC","AveExpr","t","P.Value","adj.P.Val","B")], row.names = NULL)
  write.csv(tt, file.path(od, paste0("de_", k, ".csv")), row.names = FALSE)
  cat(sprintf("%-20s %-45s up %4d down %4d (adj.P<0.05, |logFC|>0.585)\n", k, cons[k],
      sum(tt$adj.P.Val < 0.05 & tt$logFC > 0.585), sum(tt$adj.P.Val < 0.05 & tt$logFC < -0.585)))
}
# Expression matrix for ALL 19 samples (same genes), log2 CPM
yall <- calcNormFactors(y[rownames(ym), , keep.lib.sizes = FALSE])
lc <- cpm(yall, log = TRUE, prior.count = 2)
write.csv(data.frame(symbol = rownames(lc), round(lc, 4), check.names = FALSE, row.names = NULL),
          file.path(od, "log2cpm_matrix.csv"), row.names = FALSE)
write.csv(meta, file.path(od, "metadata.csv"), row.names = FALSE)
writeLines(c(paste("genes_tested", nrow(ym)), paste0("n_", names(n), " ", n),
             paste("contrast", names(cons), cons)), file.path(od, "model_info.txt"))
cat("genes tested:", nrow(ym), "| samples in model:", ncol(ym), "\n")
