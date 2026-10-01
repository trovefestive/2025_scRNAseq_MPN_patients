#!/usr/bin/env Rscript
# Stage 9 (downstream-analysis-plan.md): DEG per cluster (cluster vs all others),
# on ALL cells. Stage 7's FindAllMarkers was subsampled annotation-support; this is
# the real result.
#
# Plan thresholds: avg_log2FC > 1 and FDR < 0.2.
#   IMPORTANT: Seurat's `p_val_adj` is BONFERRONI (p_val * n_genes), not FDR. The
#   plan says FDR, so Benjamini-Hochberg is computed here from raw p_val and used
#   for the hit calling. Both columns are kept so either can be reported.
#
# Cell-cycle phase-split merge (Stage 8 decision 1): clusters 7 (Erythroid_cycling)
# and 8 (Erythroid_S_phase) are folded into Erythroblast -- both are erythroid by
# SingleR (85.9% / 99.1%) and by GYPA, differing only in cycle phase.
#   NOT merged: clusters 18 (Progenitor_G2M) and 26 (Cycling_progenitor). Stage 8
#   proposed folding them into a progenitor parent, but on inspection their lineage
#   is not supportable -- SingleR calls 26 "NK_cell" at 60%, and 26 is only 93 cells.
#   Merging would assert a lineage the data doesn't establish, so they stay separate
#   and are flagged as cell-cycle-defined, not interpretable as cell types.
#
# Set STAGE9_TEST_N=<n> to dry-run on n cells per sample -> stage9_test/.

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)
options(future.globals.maxSize = 100 * 1024^3)

LFC_CUT <- 1.0
FDR_CUT <- 0.2

test_n  <- Sys.getenv("STAGE9_TEST_N", "")
out_dir <- if (nzchar(test_n)) "stage9_test" else "stage9"
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

message("=== Loading stage8/annotated_cellcycle.rds ===")
obj <- readRDS("stage8/annotated_cellcycle.rds")
if (nzchar(test_n)) {
  n <- as.integer(test_n); set.seed(1)
  keep <- unlist(lapply(split(colnames(obj), obj$sample),
                        function(x) head(sample(x), min(n, length(x)))), use.names = FALSE)
  obj <- subset(obj, cells = keep); message("*** TEST MODE: ", ncol(obj), " cells ***")
}

ct <- as.character(obj$cell_type)
ct[ct %in% c("Erythroid_cycling", "Erythroid_S_phase")] <- "Erythroblast"
obj$cell_type_deg <- ct
Idents(obj) <- "cell_type_deg"
message("Groups for DEG: ", length(unique(ct)))
print(sort(table(ct), decreasing = TRUE))

message("\n=== FindAllMarkers (presto: ", requireNamespace("presto", quietly = TRUE), ") ===")
t0 <- Sys.time()
mk <- FindAllMarkers(obj, only.pos = FALSE, min.pct = 0.1,
                     logfc.threshold = 0.25, verbose = FALSE)
message("Elapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

# BH-FDR across the whole table (Seurat's p_val_adj is Bonferroni).
mk$FDR <- p.adjust(mk$p_val, method = "BH")
mk <- mk[order(mk$cluster, -mk$avg_log2FC), ]
write.csv(mk, file.path(out_dir, "deg_all_genes.csv"), row.names = FALSE)
message("Total tested rows: ", nrow(mk))

hits <- mk[mk$avg_log2FC > LFC_CUT & mk$FDR < FDR_CUT, ]
write.csv(hits, file.path(out_dir, "deg_hits_lfc1_fdr0.2.csv"), row.names = FALSE)
message("Hits (avg_log2FC > ", LFC_CUT, " & FDR < ", FDR_CUT, "): ", nrow(hits))

n_tbl <- hits %>% count(cluster, name = "n_hits") %>%
  left_join(as.data.frame(table(ct), responseName = "n_cells"),
            by = c("cluster" = "ct")) %>% arrange(desc(n_hits))
write.csv(n_tbl, file.path(out_dir, "deg_hit_counts.csv"), row.names = FALSE)
message("\n=== Hits per cell type ===")
print(as.data.frame(n_tbl), row.names = FALSE)

top <- hits %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 15) %>%
  summarise(top_genes = paste(gene, collapse = ", "), .groups = "drop")
write.csv(top, file.path(out_dir, "deg_top15_per_cluster.csv"), row.names = FALSE)
message("\n=== Top 15 by log2FC per cell type ===")
print(as.data.frame(top), row.names = FALSE)

# Also report how many hits survive the stricter Bonferroni, so the permissive
# FDR<0.2 threshold from the paper can be judged in context.
message("\nHits also passing Seurat's Bonferroni p_val_adj < 0.05: ",
        sum(hits$p_val_adj < 0.05), " of ", nrow(hits))

top5 <- hits %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>% pull(gene) %>% unique()
ggsave(file.path(out_dir, "plots", "deg_dotplot_top5.png"),
       DotPlot(obj, features = top5) + RotatedAxis() +
         theme(axis.text.x = element_text(size = 5)) +
         ggtitle("Top 5 DEG per cell type (log2FC > 1, FDR < 0.2)"),
       width = 26, height = 8, dpi = 150, limitsize = FALSE)

message("\n=== Stage 9 complete ===")
