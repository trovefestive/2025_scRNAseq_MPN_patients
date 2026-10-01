#!/usr/bin/env Rscript
# Stage 7 completion (consensus cell type labels) + Stage 8 (cell cycle scoring).
# Combined into one job because the labelling is a lookup and CellCycleScoring is
# cheap; both were agreed at the Stage 7 check-in.
#
# Motivating observation from Stage 7 markers: clusters 7, 8, 18 and 26 are defined
# almost entirely by proliferation genes (AURKA/CCNB1, replication histones, CDC20/
# PLK1, E2F1/MELK). That means cell cycle is partly DRIVING the clustering, so this
# stage also has to answer: does cell cycle need regressing out before Stage 9?
#
# Set STAGE8_TEST_N=<n> to dry-run on n cells per sample -> stage8_test/.

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)
options(future.globals.maxSize = 100 * 1024^3)
source("stage8_celltypes.R")

RES_COL <- "RNA_snn_res.0.8"
test_n  <- Sys.getenv("STAGE8_TEST_N", "")
out_dir <- if (nzchar(test_n)) "stage8_test" else "stage8"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Loading stage7/annotated.rds ===")
obj <- readRDS("stage7/annotated.rds")
if (nzchar(test_n)) {
  n <- as.integer(test_n); set.seed(1)
  keep <- unlist(lapply(split(colnames(obj), obj$sample),
                        function(x) head(sample(x), min(n, length(x)))), use.names = FALSE)
  obj <- subset(obj, cells = keep)
  message("*** TEST MODE: ", ncol(obj), " cells ***")
}

# ---- Consensus labels -----------------------------------------------------
cl <- as.character(obj[[RES_COL]][, 1])
stopifnot(all(unique(cl) %in% names(CELLTYPE)))
obj$cell_type <- unname(CELLTYPE[cl])
obj$lineage   <- unname(LINEAGE[obj$cell_type])
stopifnot(!any(is.na(obj$cell_type)), !any(is.na(obj$lineage)))
Idents(obj) <- "cell_type"
message("Assigned ", length(unique(obj$cell_type)), " cell types in ",
        length(unique(obj$lineage)), " lineages")

# ---- Cell cycle scoring ---------------------------------------------------
message("\n=== CellCycleScoring ===")
s.g   <- intersect(cc.genes.updated.2019$s.genes,   rownames(obj))
g2m.g <- intersect(cc.genes.updated.2019$g2m.genes, rownames(obj))
message("S genes: ", length(s.g), "/43 | G2M genes: ", length(g2m.g), "/54 present")
obj <- CellCycleScoring(obj, s.features = s.g, g2m.features = g2m.g, set.ident = FALSE)
obj$Phase <- factor(obj$Phase, levels = c("G1", "S", "G2M"))
message("\nOverall phase distribution:")
print(round(100 * prop.table(table(obj$Phase)), 2))
saveRDS(obj, file.path(out_dir, "annotated_cellcycle.rds"))
message("Saved: ", file.path(out_dir, "annotated_cellcycle.rds"))

# ---- Phase by cluster: does cell cycle drive the clustering? --------------
message("\n=== Phase composition per cluster ===")
ph <- as.data.frame.matrix(table(cl, obj$Phase))
ph$n_cells <- rowSums(ph)
ph$cell_type <- unname(CELLTYPE[rownames(ph)])
ph$pct_G1  <- round(100 * ph$G1  / ph$n_cells, 1)
ph$pct_S   <- round(100 * ph$S   / ph$n_cells, 1)
ph$pct_G2M <- round(100 * ph$G2M / ph$n_cells, 1)
ph$pct_cycling <- ph$pct_S + ph$pct_G2M
ph <- ph[order(-ph$pct_cycling), c("cell_type","n_cells","pct_G1","pct_S","pct_G2M","pct_cycling")]
write.csv(cbind(cluster = rownames(ph), ph), file.path(out_dir, "phase_per_cluster.csv"),
          row.names = FALSE)
print(ph)

message("\n>>> Clusters >50% cycling (cell cycle is a defining feature, not a nuisance):")
print(ph[ph$pct_cycling > 50, c("cell_type","n_cells","pct_cycling")])

# ---- Composition: cell type x disease stage ------------------------------
message("\n=== Cell type composition by disease stage (% within stage) ===")
comp <- as.data.frame.matrix(table(obj$cell_type, obj$disease_stage))
comp_pct <- round(100 * sweep(comp, 2, colSums(comp), "/"), 2)
write.csv(cbind(cell_type = rownames(comp), comp, comp_pct),
          file.path(out_dir, "composition_by_stage.csv"), row.names = FALSE)
print(comp_pct)

lin <- as.data.frame.matrix(table(obj$lineage, obj$disease_stage))
lin_pct <- round(100 * sweep(lin, 2, colSums(lin), "/"), 2)
write.csv(cbind(lineage = rownames(lin), lin, lin_pct),
          file.path(out_dir, "composition_by_stage_lineage.csv"), row.names = FALSE)
message("\n=== Lineage composition by disease stage (%) ===")
print(lin_pct)

# Per-sample composition, so any stage-level difference can be checked against
# per-patient spread rather than being read off pooled cells (n=19 patients, not
# 111,861 independent observations).
persamp <- as.data.frame.matrix(table(obj$sample, obj$cell_type))
persamp_pct <- round(100 * sweep(persamp, 1, rowSums(persamp), "/"), 2)
write.csv(cbind(sample = rownames(persamp_pct), persamp_pct),
          file.path(out_dir, "composition_per_sample.csv"), row.names = FALSE)

# ---- Plots ----------------------------------------------------------------
message("\n=== Plotting ===")
ggsave(file.path(plot_dir, "umap_cell_type.png"),
       DimPlot(obj, reduction = "umap", group.by = "cell_type", label = TRUE,
               label.size = 3, repel = TRUE, raster = TRUE) + NoLegend() +
         ggtitle("Consensus cell types (res 0.8)"), width = 10, height = 9, dpi = 150)
ggsave(file.path(plot_dir, "umap_lineage.png"),
       DimPlot(obj, reduction = "umap", group.by = "lineage", label = TRUE,
               label.size = 4, repel = TRUE, raster = TRUE) +
         ggtitle("Lineage"), width = 10, height = 8, dpi = 150)
ggsave(file.path(plot_dir, "umap_phase.png"),
       DimPlot(obj, reduction = "umap", group.by = "Phase", raster = TRUE) +
         ggtitle("Cell cycle phase"), width = 9, height = 7, dpi = 150)
ggsave(file.path(plot_dir, "umap_cc_scores.png"),
       FeaturePlot(obj, reduction = "umap", features = c("S.Score", "G2M.Score"),
                   raster = TRUE, ncol = 2), width = 12, height = 5, dpi = 150)

pd <- as.data.frame(table(obj$cell_type, obj$Phase))
names(pd) <- c("cell_type", "Phase", "n")
ggsave(file.path(plot_dir, "phase_by_celltype.png"),
       ggplot(pd, aes(cell_type, n, fill = Phase)) +
         geom_col(position = "fill") + coord_flip() +
         labs(y = "fraction of cells", x = NULL, title = "Cell cycle phase by cell type") +
         theme_classic(),
       width = 8, height = 8, dpi = 150)

cd <- as.data.frame(table(obj$lineage, obj$disease_stage))
names(cd) <- c("lineage", "stage", "n")
ggsave(file.path(plot_dir, "composition_by_stage.png"),
       ggplot(cd, aes(stage, n, fill = lineage)) +
         geom_col(position = "fill") +
         labs(y = "fraction of cells", x = NULL, title = "Lineage composition by disease stage") +
         theme_classic(),
       width = 7, height = 6, dpi = 150)

message("\n=== Stage 7 labels + Stage 8 cell cycle complete ===")
