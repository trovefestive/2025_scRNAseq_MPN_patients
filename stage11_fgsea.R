#!/usr/bin/env Rscript
# Stage 11 (downstream-analysis-plan.md): fgsea + MSigDB HALLMARK on per-cell-type
# DEG rankings.
#
# METHOD NOTE -- why this does NOT reuse stage9/deg_all_genes.csv:
# Stage 9 ran FindAllMarkers with logfc.threshold = 0.25 and min.pct = 0.1, so its
# gene list is TRUNCATED. GSEA assumes a complete ranked list over all measured
# genes; feeding it a pre-filtered list inflates enrichment for the retained tail
# and is a real bias, not a technicality. So rankings are recomputed here over ALL
# genes with presto::wilcoxauc (the same Wilcoxon engine Seurat uses underneath),
# which returns every gene for every group cheaply.
#
# Ranking statistic = AUC - 0.5. AUC is bounded in [0,1] and immune to the p-value
# underflow that hits -log10(p) at n = 111,861 cells (many p_val are exactly 0).
#
# Gene sets: genesets/hallmark_human.rds, pre-built on the login node (msigdbr 25.x
# needs network; compute nodes have none).
#
# Set STAGE11_TEST_N=<n> to dry-run -> stage11_test/.

suppressPackageStartupMessages({
  library(Seurat); library(presto); library(fgsea); library(dplyr); library(ggplot2)
})

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)
options(future.globals.maxSize = 100 * 1024^3)
set.seed(42)

test_n  <- Sys.getenv("STAGE11_TEST_N", "")
out_dir <- if (nzchar(test_n)) "stage11_test" else "stage11"
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

message("=== Loading stage8/annotated_cellcycle.rds ===")
obj <- readRDS("stage8/annotated_cellcycle.rds")
if (nzchar(test_n)) {
  n <- as.integer(test_n); set.seed(1)
  keep <- unlist(lapply(split(colnames(obj), obj$sample),
                        function(x) head(sample(x), min(n, length(x)))), use.names = FALSE)
  obj <- subset(obj, cells = keep); message("*** TEST MODE: ", ncol(obj), " cells ***")
  set.seed(42)
}

# Same grouping as Stage 9 so the two stages are comparable.
ct <- as.character(obj$cell_type)
ct[ct %in% c("Erythroid_cycling", "Erythroid_S_phase")] <- "Erythroblast"
obj$cell_type_deg <- ct
message("Groups: ", length(unique(ct)))

message("\n=== presto::wilcoxauc over ALL genes ===")
t0 <- Sys.time()
wa <- wilcoxauc(obj, group_by = "cell_type_deg", assay = "data", seurat_assay = "RNA")
message("Elapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min | rows: ", nrow(wa))
message("Genes per group: ", round(nrow(wa) / length(unique(wa$group))))
write.csv(wa, file.path(out_dir, "wilcoxauc_full_rankings.csv.gz"), row.names = FALSE)

hall <- readRDS("genesets/hallmark_human.rds")
message("HALLMARK sets: ", length(hall))

message("\n=== fgsea per cell type ===")
res <- list()
for (g in sort(unique(wa$group))) {
  d <- wa[wa$group == g, ]
  stats <- d$auc - 0.5
  names(stats) <- d$feature
  stats <- stats[is.finite(stats)]
  stats <- sort(stats[!duplicated(names(stats))], decreasing = TRUE)
  # eps = 0 lets fgseaMultilevel estimate arbitrarily small p-values. With the
  # default eps (1e-10) the most extreme enrichments (e.g. COAGULATION in canonical
  # megakaryocytes, ES = 0.60) returned NA p-values -- a convergence limit, NOT an
  # absence of signal, and easily misread as a null result.
  fg <- fgsea(pathways = hall, stats = stats, minSize = 10, maxSize = 500,
              nPermSimple = 10000, eps = 0)
  if (nrow(fg)) { fg$cell_type <- g; res[[g]] <- fg }
  message("  ", g, ": ", nrow(fg), " pathways, ", sum(fg$padj < 0.05, na.rm = TRUE), " at padj<0.05")
}
gs <- bind_rows(res) %>% arrange(cell_type, padj)
gs$leadingEdge <- sapply(gs$leadingEdge, paste, collapse = ";")
write.csv(gs, file.path(out_dir, "fgsea_hallmark_all.csv"), row.names = FALSE)

sig <- gs[!is.na(gs$padj) & gs$padj < 0.05, ]
write.csv(sig, file.path(out_dir, "fgsea_hallmark_sig.csv"), row.names = FALSE)
message("\nSignificant (padj<0.05): ", nrow(sig), " cell-type x pathway pairs")

message("\n=== Top enriched pathway per cell type (by NES) ===")
top <- sig %>% group_by(cell_type) %>% slice_max(NES, n = 3) %>%
  summarise(top_pathways = paste0(sub("HALLMARK_", "", pathway),
                                  " (NES=", round(NES, 2), ")", collapse = "; "), .groups = "drop")
print(as.data.frame(top), row.names = FALSE)
write.csv(top, file.path(out_dir, "fgsea_top3_per_celltype.csv"), row.names = FALSE)

# Independent check of the Stage 9/10 fibrosis result: does an unsupervised
# HALLMARK scan flag fibrogenic/angiogenic programmes in the MK compartment
# without being told to look for them?
foc_paths <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_ANGIOGENESIS",
               "HALLMARK_TGF_BETA_SIGNALING","HALLMARK_COAGULATION",
               "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_HYPOXIA")
foc_ct <- c("Megakaryocyte_mature","Megakaryocyte","Megakaryocyte_prog","MSC_stromal")
message("\n>>> Fibrogenic/angiogenic HALLMARK in MK + stromal compartments:")
chk <- gs %>% filter(cell_type %in% foc_ct, pathway %in% foc_paths) %>%
  mutate(pathway = sub("HALLMARK_", "", pathway), NES = round(NES, 2),
         padj = signif(padj, 3)) %>%
  select(cell_type, pathway, NES, padj, size) %>% arrange(cell_type, desc(NES))
print(as.data.frame(chk), row.names = FALSE)
write.csv(chk, file.path(out_dir, "fgsea_fibrosis_check.csv"), row.names = FALSE)

# ---- NES heatmap ----------------------------------------------------------
hm <- gs %>% mutate(pathway = sub("HALLMARK_", "", pathway)) %>%
  group_by(pathway) %>% filter(any(padj < 0.05, na.rm = TRUE)) %>% ungroup()
ggsave(file.path(out_dir, "plots", "fgsea_nes_heatmap.png"),
       ggplot(hm, aes(cell_type, pathway, fill = NES)) +
         geom_tile() +
         geom_point(data = subset(hm, padj < 0.05), size = 0.5, colour = "black") +
         scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
         labs(x = NULL, y = NULL,
              title = "HALLMARK enrichment (NES) by cell type",
              subtitle = "dot = padj < 0.05") +
         theme_minimal(base_size = 8) +
         theme(axis.text.x = element_text(angle = 60, hjust = 1)),
       width = 12, height = 10, dpi = 150)

ggsave(file.path(out_dir, "plots", "fgsea_mk_fibrosis.png"),
       ggplot(chk, aes(reorder(pathway, NES), NES, fill = padj < 0.05)) +
         geom_col() + coord_flip() + facet_wrap(~cell_type) +
         scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = "#B2182B"),
                           name = "padj < 0.05") +
         labs(x = NULL, title = "Fibrogenic/angiogenic HALLMARK in MK + stromal") +
         theme_bw(base_size = 9),
       width = 10, height = 6, dpi = 150)

message("\n=== Stage 11 complete ===")
