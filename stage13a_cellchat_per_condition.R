#!/usr/bin/env Rscript
# Stage 13a: one CellChat object per condition (Reactive / PrePMF / PMF).
# Comparison (liftCellChat + mergeCellChat) is Stage 13b, kept separate so a
# failure in the slow per-condition step doesn't force re-running everything.
#
# DB: CellChatDB.human via subsetDB() default = protein signaling only
#     (Secreted, ECM-Receptor, Cell-Cell Contact; Non-protein Signaling excluded).
# computeCommunProb: type = "triMean" (CellChat default, robust to dropout);
#     population.size = FALSE (default) -- abundance changes between conditions
#     are reported from composition (Stage 8), not folded into edge weights.
# filterCommunication(min.cells = 10): plan's "exclude groups with <10 cells".
#
# Usage: Rscript stage13a_cellchat_per_condition.R <Reactive|PrePMF|PMF>
# Set STAGE13_TEST_N=<n> to subsample n cells/sample -> stage13_test/.

suppressPackageStartupMessages({ library(Seurat); library(CellChat); library(future) })

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)
source("stage13_config.R")

cond <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(cond %in% unique(CONDITION))
NCORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
test_n  <- Sys.getenv("STAGE13_TEST_N", "")
out_dir <- if (nzchar(test_n)) "stage13_test" else "stage13"
dir.create(out_dir, showWarnings = FALSE)
options(future.globals.maxSize = 100 * 1024^3)
# CellChat's permutation test (computeCommunProb, nboot = 100) runs inside
# future.apply calls that do not request parallel-safe RNG (no future.seed), so
# future warns "UNRELIABLE VALUE". Workers in a multisession plan are fresh R
# sessions with independent default seeds, so the permutation nulls are still
# random, but p-values are NOT bit-reproducible across reruns. Accepted and
# documented rather than running sequentially (hours on the 85k-cell PMF object).
options(future.rng.onMisuse = "ignore")
plan("multisession", workers = NCORES)
set.seed(42)

message("=== Stage 13a: ", cond, " ===")
obj <- readRDS("stage8/annotated_cellcycle.rds")
obj$cc_group  <- unname(CC_GROUP[as.character(obj$cell_type)])
obj$condition <- unname(CONDITION[as.character(obj$disease_stage)])
stopifnot(!any(is.na(obj$condition)))
keep <- !is.na(obj$cc_group) & !(obj$sample %in% EXCLUDE_SAMPLES) & obj$condition == cond
obj <- subset(obj, cells = colnames(obj)[keep])

if (nzchar(test_n)) {
  n <- as.integer(test_n)
  keep <- unlist(lapply(split(colnames(obj), obj$sample),
                        function(x) head(sample(x), min(n, length(x)))), use.names = FALSE)
  obj <- subset(obj, cells = keep)
  message("*** TEST MODE ***")
}
message("Cells: ", ncol(obj), " | patients: ", length(unique(obj$sample)))
print(sort(table(obj$cc_group), decreasing = TRUE))

data_in <- GetAssayData(obj, assay = "RNA", layer = "data")
meta    <- data.frame(group = obj$cc_group, samples = obj$sample, row.names = colnames(obj))
cc <- createCellChat(object = data_in, meta = meta, group.by = "group")
cc@DB <- subsetDB(CellChatDB.human)
rm(obj, data_in); invisible(gc())

t0 <- Sys.time()
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- computeCommunProb(cc, type = "triMean")
cc <- filterCommunication(cc, min.cells = 10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)
cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
message("CellChat elapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

df <- subsetCommunication(cc)
df$category <- categorize(df$pathway_name)
write.csv(df, file.path(out_dir, paste0("lr_edges_", cond, ".csv")), row.names = FALSE)
message("Significant LR edges: ", nrow(df), " | pathways: ", length(cc@netP$pathways))
print(table(df$category))
saveRDS(cc, file.path(out_dir, paste0("cellchat_", cond, ".rds")))
message("=== done: ", cond, " ===")
