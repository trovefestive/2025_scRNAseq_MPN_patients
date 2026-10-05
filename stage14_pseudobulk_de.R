#!/usr/bin/env Rscript
# Stage 14 (added after review): within-cell-type disease comparison by
# per-patient pseudobulk, edgeR quasi-likelihood.
#
# Contrast: PMF-MF3 vs early disease (prefibrotic + PMF-MF2). Reactive (n = 2,
# non-MPN) is excluded: too few patients for a within-cell-type test.
# Design ~ sex + group: MF3 is 8 male / 1 female, early is 2 male / 6 female, so
# sex is nearly confounded with grade. Adjusting for it costs power; X/Y genes
# are flagged and excluded from the pathway ranking.
# Unit of replication = patient. A patient contributes a cell type only with
# >= 20 cells; a cell type is tested only with >= 3 patients per group.
# Label correction from the review: cluster 13 ("Megakaryocyte") is reported
# as Platelet_like_MK; cluster 22 ("Megakaryocyte_mature") as Megakaryocyte_intact.

suppressPackageStartupMessages({ library(Seurat); library(edgeR); library(fgsea) })
proj_dir <- Sys.getenv("PROJ_DIR", getwd()); setwd(proj_dir)
out_dir <- "stage14"; dir.create(out_dir, showWarnings = FALSE)
MIN_CELLS <- 20; MIN_PAT <- 3

obj <- readRDS("stage8/annotated_cellcycle.rds")
ct <- as.character(obj$cell_type)
ct[ct %in% c("Erythroid_cycling", "Erythroid_S_phase")] <- "Erythroblast"
ct[ct == "Megakaryocyte"] <- "Platelet_like_MK"
ct[ct == "Megakaryocyte_mature"] <- "Megakaryocyte_intact"
grp <- ifelse(obj$disease_stage == "PMF(MF3)", "MF3",
       ifelse(obj$disease_stage %in% c("prefibrotic", "PMF(MF2)"), "Early", NA))
keep <- !is.na(grp)
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[, keep]
meta <- data.frame(sample = obj$sample[keep], ct = ct[keep], group = grp[keep], sex = obj$sex[keep])
rm(obj); invisible(gc())

chrom <- read.delim("genesets/gene_chrom.tsv", header = FALSE, col.names = c("chr", "gene"))
sexchr <- unique(chrom$gene[chrom$chr %in% c("chrX", "chrY")])
hall <- readRDS("genesets/hallmark_human.rds")
focus_genes <- c("SERPINE1","VEGFC","PDGFB","PDGFA","CXCL1","THBS1","TGFB1","PF4","LOX","TIMP1")

pd <- meta[!duplicated(meta$sample), c("sample", "group", "sex")]
cat("Patients:\n"); print(table(group = pd$group, sex = pd$sex))

summary_rows <- list(); de_all <- list(); gsea_all <- list()
for (k in sort(unique(meta$ct))) {
  m <- meta[meta$ct == k, ]
  n_by_pt <- table(m$sample); pts <- names(n_by_pt)[n_by_pt >= MIN_CELLS]
  pinfo <- pd[match(pts, pd$sample), ]
  n_mf3 <- sum(pinfo$group == "MF3"); n_early <- sum(pinfo$group == "Early")
  if (n_mf3 < MIN_PAT || n_early < MIN_PAT) {
    summary_rows[[k]] <- data.frame(cell_type = k, n_MF3 = n_mf3, n_Early = n_early, tested = FALSE,
                                    n_genes = NA, n_FDR05 = NA, n_FDR05_autosomal = NA, sex_model = NA)
    next
  }
  # Pseudobulk: sum raw counts per patient
  pb <- sapply(pts, function(p) Matrix::rowSums(counts[, which(meta$ct == k & meta$sample == p), drop = FALSE]))
  grp_f <- factor(pinfo$group, levels = c("Early", "MF3"))
  sex_f <- factor(pinfo$sex)
  # Sex is only adjustable if both sexes appear and sex is not identical to group
  # Adjust for sex only if the sex + group design is full rank. If sex perfectly
  # separates the groups for this cell type's patients, grade and sex cannot be
  # separated at all -- fall back to ~ group and label the result as confounded.
  d_sex <- if (nlevels(sex_f) == 2) model.matrix(~ sex_f + grp_f) else NULL
  use_sex <- !is.null(d_sex) && qr(d_sex)$rank == ncol(d_sex)
  sex_note <- if (use_sex) "~ sex + group" else if (nlevels(sex_f) == 2)
    "~ group; SEX FULLY CONFOUNDED" else "~ group (one sex only)"
  design <- if (use_sex) d_sex else model.matrix(~ grp_f)
  y <- DGEList(pb, group = grp_f)
  y <- y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  res <- topTags(glmQLFTest(fit, coef = "grp_fMF3"), n = Inf)$table
  res$gene <- rownames(res); res$cell_type <- k; res$sex_chr <- res$gene %in% sexchr
  de_all[[k]] <- res
  summary_rows[[k]] <- data.frame(cell_type = k, n_MF3 = n_mf3, n_Early = n_early, tested = TRUE,
    n_genes = nrow(res), n_FDR05 = sum(res$FDR < 0.05),
    n_FDR05_autosomal = sum(res$FDR < 0.05 & !res$sex_chr),
    sex_model = sex_note)
  # GSEA on signed -log10 p, autosomal genes only
  a <- res[!res$sex_chr, ]
  st <- sign(a$logFC) * -log10(a$PValue); names(st) <- a$gene
  st <- sort(st[is.finite(st)], decreasing = TRUE)
  g <- fgsea(hall, st, minSize = 10, maxSize = 500, eps = 0)
  g$cell_type <- k; g$leadingEdge <- sapply(g$leadingEdge, paste, collapse = ";")
  gsea_all[[k]] <- g
}
summ <- do.call(rbind, summary_rows)
de   <- do.call(rbind, de_all)
gs   <- do.call(rbind, gsea_all)
write.csv(summ, file.path(out_dir, "pseudobulk_summary.csv"), row.names = FALSE)
write.csv(de,   file.path(out_dir, "pseudobulk_de_MF3_vs_early.csv"), row.names = FALSE)
write.csv(gs,   file.path(out_dir, "pseudobulk_fgsea_hallmark.csv"), row.names = FALSE)
cat("\n== Cell types tested (MF3 vs Early, per-patient pseudobulk) ==\n"); print(summ, row.names = FALSE)
cat("\n== Top autosomal DEGs per tested cell type (FDR < 0.05) ==\n")
for (k in names(de_all)) { r <- de_all[[k]]; r <- r[r$FDR < 0.05 & !r$sex_chr, ]
  if (nrow(r)) cat(sprintf("%-22s up: %s\n%-22s down: %s\n", k,
      paste(head(r$gene[r$logFC > 0], 10), collapse = ", "), "",
      paste(head(r$gene[r$logFC < 0], 10), collapse = ", "))) }
cat("\n== Focus genes in MK compartments ==\n")
fk <- de[de$gene %in% focus_genes & grepl("MK|Megakaryocyte", de$cell_type),
         c("cell_type","gene","logFC","PValue","FDR")]
fk[, 3:5] <- signif(fk[, 3:5], 3); print(fk[order(fk$cell_type, fk$gene), ], row.names = FALSE)
cat("\n== HALLMARK (padj < 0.05) ==\n")
s <- gs[!is.na(gs$padj) & gs$padj < 0.05, c("cell_type","pathway","NES","padj")]
s$pathway <- sub("HALLMARK_", "", s$pathway); s$NES <- round(s$NES, 2); s$padj <- signif(s$padj, 3)
print(s[order(s$cell_type, -abs(s$NES)), ], row.names = FALSE)
