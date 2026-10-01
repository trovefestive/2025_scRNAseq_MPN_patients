#!/usr/bin/env Rscript
# Stage 10 (downstream-analysis-plan.md): gene signature scores via UCell.
# The plan left gene sets "TBD"; chosen here = 50 MSigDB HALLMARK + 5 literature
# PMF sets (stage10_genesets.R). All 50 HALLMARK are scored rather than a
# hand-picked few, so the fibrosis/angiogenesis result can't be cherry-picking.
#
# Gene sets are read from genesets/hallmark_human.rds, pre-built on the LOGIN node:
# msigdbr 25.x fetches over the network and compute nodes have no internet (same
# constraint as celldex in Stage 7).
#
# Hypothesis under test (from Stage 9): cluster 22 Megakaryocyte_mature carries a
# pro-fibrotic/angiogenic program (SERPINE1/INHBA/VEGFC/PDGFA). Scoring independent
# literature sets is a non-circular check of that.
#
# STATISTICS: all stage-level comparisons are PER PATIENT (mean score per sample
# per cell type, n=19 patients), never pooled cells. Pooling 111,861 cells would
# manufacture tiny p-values from what is a 19-sample study with n=2 controls.
#
# Set STAGE10_TEST_N=<n> to dry-run -> stage10_test/.

suppressPackageStartupMessages({ library(Seurat); library(UCell); library(dplyr); library(ggplot2) })

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)
options(future.globals.maxSize = 100 * 1024^3)
source("stage10_genesets.R")

NCORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
test_n  <- Sys.getenv("STAGE10_TEST_N", "")
out_dir <- if (nzchar(test_n)) "stage10_test" else "stage10"
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

message("=== Loading stage8/annotated_cellcycle.rds ===")
obj <- readRDS("stage8/annotated_cellcycle.rds")
if (nzchar(test_n)) {
  n <- as.integer(test_n); set.seed(1)
  keep <- unlist(lapply(split(colnames(obj), obj$sample),
                        function(x) head(sample(x), min(n, length(x)))), use.names = FALSE)
  obj <- subset(obj, cells = keep); message("*** TEST MODE: ", ncol(obj), " cells ***")
}

hall <- readRDS("genesets/hallmark_human.rds")
sets <- c(hall, CUSTOM_SETS)
sets <- lapply(sets, function(g) intersect(g, rownames(obj)))
sets <- sets[lengths(sets) >= 10]
message("Scoring ", length(sets), " gene sets (", length(hall), " HALLMARK + ",
        length(CUSTOM_SETS), " custom, dropped those with <10 genes present)")
for (k in names(CUSTOM_SETS))
  message("  ", k, ": ", length(sets[[k]]), "/", length(CUSTOM_SETS[[k]]), " genes present")

message("\n=== AddModuleScore_UCell (", NCORES, " cores) ===")
t0 <- Sys.time()
obj <- AddModuleScore_UCell(obj, features = sets, ncores = NCORES, name = "")
message("Elapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
saveRDS(obj, file.path(out_dir, "ucell_scored.rds"))

# ---- Per-cell-type score summary -----------------------------------------
set_names <- names(sets)
md <- obj@meta.data[, c("sample", "disease_stage", "cell_type", set_names)]

message("\n=== Mean score per cell type (all cells) ===")
by_ct <- md %>% group_by(cell_type) %>%
  summarise(across(all_of(set_names), mean), n = n(), .groups = "drop")
write.csv(by_ct, file.path(out_dir, "score_by_celltype.csv"), row.names = FALSE)

key <- intersect(c("PMF_fibrosis","MK_profibrotic_secretome","MPN_inflammation",
                   "JAK_STAT_targets","Stromal_niche","HALLMARK_TGF_BETA_SIGNALING",
                   "HALLMARK_ANGIOGENESIS","HALLMARK_INFLAMMATORY_RESPONSE",
                   "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), set_names)
message("\nKey sets by cell type (top 8 cell types for PMF_fibrosis):")
print(as.data.frame(by_ct %>% select(cell_type, n, all_of(key)) %>%
      arrange(desc(PMF_fibrosis)) %>% head(8) %>% mutate(across(where(is.numeric), ~round(.x, 4)))),
      row.names = FALSE)

# ---- PER-PATIENT stage comparison ----------------------------------------
# One value per (sample, cell type): the mean score over that patient's cells.
# Stage-level statistics are then computed over PATIENTS, not cells.
message("\n=== Per-patient scores by cell type and stage ===")
per_pt <- md %>% group_by(sample, disease_stage, cell_type) %>%
  summarise(n_cells = n(), across(all_of(set_names), mean), .groups = "drop") %>%
  filter(n_cells >= 20)   # a patient needs >=20 cells of a type to contribute
write.csv(per_pt, file.path(out_dir, "score_per_patient_celltype.csv"), row.names = FALSE)

stage_sum <- per_pt %>% group_by(cell_type, disease_stage) %>%
  summarise(n_patients = n(), across(all_of(key), ~round(median(.x), 4)), .groups = "drop")
write.csv(stage_sum, file.path(out_dir, "score_by_stage_per_patient.csv"), row.names = FALSE)

message("\n>>> PMF_fibrosis, per-patient median by stage, MK/stromal cell types:")
foc <- c("Megakaryocyte_mature","Megakaryocyte","Megakaryocyte_prog","MSC_stromal","MEP","HSC_MPP")
print(as.data.frame(stage_sum %>% filter(cell_type %in% foc) %>%
      select(cell_type, disease_stage, n_patients, PMF_fibrosis,
             MK_profibrotic_secretome, HALLMARK_ANGIOGENESIS, HALLMARK_TGF_BETA_SIGNALING) %>%
      arrange(cell_type, disease_stage)), row.names = FALSE)

# Cross-cell-type ranking: is cluster 22 actually the top fibrosis scorer?
message("\n>>> Cross-check of the Stage 9 hypothesis -- PMF_fibrosis rank across cell types:")
print(as.data.frame(by_ct %>% select(cell_type, n, PMF_fibrosis, MK_profibrotic_secretome) %>%
      arrange(desc(PMF_fibrosis)) %>% mutate(across(where(is.numeric), ~round(.x, 4)))),
      row.names = FALSE)

# ---- Plots ----------------------------------------------------------------
message("\n=== Plotting ===")
for (s in intersect(c("PMF_fibrosis","MK_profibrotic_secretome","MPN_inflammation",
                      "HALLMARK_ANGIOGENESIS","HALLMARK_TGF_BETA_SIGNALING"), set_names)) {
  ggsave(file.path(out_dir, "plots", paste0("umap_", s, ".png")),
         FeaturePlot(obj, features = s, reduction = "umap", raster = TRUE) +
           ggtitle(paste("UCell:", s)), width = 7, height = 6, dpi = 150)
}

ggsave(file.path(out_dir, "plots", "vln_PMF_fibrosis_by_celltype.png"),
       VlnPlot(obj, features = "PMF_fibrosis", group.by = "cell_type", pt.size = 0) +
         NoLegend() + theme(axis.text.x = element_text(size = 7)) +
         ggtitle("PMF fibrosis signature by cell type"),
       width = 11, height = 6, dpi = 150)

# Per-patient dots by stage for the MK compartment -- shows the actual n behind
# every stage-level claim instead of hiding it behind pooled cells.
pp <- per_pt %>% filter(cell_type %in% foc)
pp$disease_stage <- factor(pp$disease_stage,
                           levels = c("Reactive","prefibrotic","PMF(MF2)","PMF(MF3)"))
ggsave(file.path(out_dir, "plots", "per_patient_fibrosis_by_stage.png"),
       ggplot(pp, aes(disease_stage, PMF_fibrosis)) +
         geom_boxplot(outlier.shape = NA, fill = "grey92") +
         geom_jitter(width = 0.15, size = 1.8, alpha = 0.85) +
         facet_wrap(~cell_type, scales = "free_y") +
         labs(x = NULL, y = "PMF_fibrosis (per-patient mean)",
              title = "Fibrosis signature by stage - one point per patient",
              subtitle = "n = 2 Reactive, 5 prefibrotic, 3 MF2, 9 MF3") +
         theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1)),
       width = 11, height = 7, dpi = 150)

message("\n=== Stage 10 complete ===")
