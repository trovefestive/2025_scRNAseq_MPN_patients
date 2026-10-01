#!/usr/bin/env Rscript
# Stage 13b: cross-condition CellChat comparison (Reactive / PrePMF / PMF).
# Inputs: stage13/cellchat_<cond>.rds from Stage 13a.
#
# liftCellChat is required because the conditions do not share the same cell
# groups (MK_profibrotic has 0 cells in Reactive). Lifted groups carry zero
# signaling, so "absent in Reactive" shows as 0, not as a missing row.
#
# INTERPRETATION LIMIT: CellChat pools cells within a condition. There is no
# patient-level replication here (Reactive = 2 patients). Differences below are
# descriptive -- rankNet's do.stat p-values test LR-pair probabilities, not
# patients, and must not be read as patient-level significance.

suppressPackageStartupMessages({ library(CellChat); library(ggplot2); library(patchwork) })
proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)
source("stage13_config.R")
out_dir <- "stage13"; plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
safe_plot <- function(file, expr, w = 9, h = 7) {
  tryCatch({ png(file.path(plot_dir, file), width = w, height = h, units = "in", res = 150)
             print(expr); dev.off() },
           error = function(e) { try(dev.off(), silent = TRUE)
                                 message("  plot failed (", file, "): ", conditionMessage(e)) })
}

conds <- c("Reactive", "PrePMF", "PMF")
cl <- lapply(conds, function(c) readRDS(file.path(out_dir, paste0("cellchat_", c, ".rds"))))
names(cl) <- conds
for (c in conds) message(c, ": ", nlevels(cl[[c]]@idents), " groups, ",
                         length(cl[[c]]@netP$pathways), " pathways")

# ---- Totals per condition ----------------------------------------------------
tot <- data.frame(condition = conds,
                  n_groups       = sapply(cl, function(x) nlevels(x@idents)),
                  n_pathways     = sapply(cl, function(x) length(x@netP$pathways)),
                  n_interactions = sapply(cl, function(x) sum(x@net$count)),
                  total_strength = sapply(cl, function(x) round(sum(x@net$weight), 3)))
write.csv(tot, file.path(out_dir, "totals_by_condition.csv"), row.names = FALSE)
message("\n=== Totals ==="); print(tot, row.names = FALSE)

# ---- Edges by LR category (from Stage 13a tables) -------------------------
edges <- do.call(rbind, lapply(conds, function(c) {
  d <- read.csv(file.path(out_dir, paste0("lr_edges_", c, ".csv"))); d$condition <- c; d }))
cat_tab <- as.data.frame.matrix(table(edges$category, factor(edges$condition, levels = conds)))
write.csv(cbind(category = rownames(cat_tab), cat_tab), file.path(out_dir, "edges_by_category.csv"),
          row.names = FALSE)
message("\n=== Significant LR edges by category ==="); print(cat_tab)

# ---- Lift to a common group set and merge -----------------------------------
group.new <- sort(unique(unlist(lapply(cl, function(x) levels(x@idents)))))
cl <- lapply(cl, liftCellChat, group.new = group.new)
merged <- mergeCellChat(cl, add.names = names(cl))
saveRDS(merged, file.path(out_dir, "cellchat_merged.rds"))

# ---- Information flow per pathway, pairwise ---------------------------------
# rankNet(mode = "comparison") over a pair of conditions; return.data gives the
# per-pathway information flow (sum of communication probability) per condition.
pairs <- list(PMF_vs_Reactive = c(1, 3), PMF_vs_PrePMF = c(2, 3), PrePMF_vs_Reactive = c(1, 2))
flow_all <- list()
for (nm in names(pairs)) {
  r <- tryCatch(rankNet(merged, mode = "comparison", comparison = pairs[[nm]], stacked = FALSE,
                        do.stat = TRUE, return.data = TRUE),
                error = function(e) { message("rankNet ", nm, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  d <- r$signaling.contribution
  d$comparison <- nm
  d$category <- categorize(as.character(d$name))
  flow_all[[nm]] <- d
  safe_plot(paste0("infoflow_", nm, ".png"),
            rankNet(merged, mode = "comparison", comparison = pairs[[nm]], stacked = TRUE,
                    do.stat = TRUE), w = 8, h = 12)
}
flow <- do.call(rbind, flow_all)
write.csv(flow, file.path(out_dir, "information_flow_pairwise.csv"), row.names = FALSE)

# ---- Fibrosis-relevant senders: MK_profibrotic and MSC_stromal --------------
# Built from the per-condition edge tables (not the lifted object) so the
# numbers are exactly what each condition's CellChat inferred.
foc_pw <- c("TGFb","PDGF","VEGF","ACTIVIN","BMP","FGF","THBS","SPP1","FN1","COLLAGEN",
            "LAMININ","ANGPT","CXCL","CCL","MIF","GALECTIN","TENASCIN","PERIOSTIN")
for (src in c("MK_profibrotic", "MSC_stromal", "MK")) {
  d <- edges[edges$source == src, c("condition","source","target","ligand","receptor",
                                     "prob","pval","pathway_name","category")]
  d <- d[order(d$condition, -d$prob), ]
  write.csv(d, file.path(out_dir, paste0("outgoing_", src, ".csv")), row.names = FALSE)
  message("\n=== ", src, " outgoing: edges per condition ===")
  print(table(factor(d$condition, levels = conds)))
  s <- aggregate(prob ~ condition + pathway_name, data = d[d$pathway_name %in% foc_pw, ], FUN = sum)
  if (nrow(s)) {
    s <- s[order(s$condition, -s$prob), ]
    message("  fibrosis-relevant pathways, summed probability:")
    print(s, row.names = FALSE)
  }
}

# ---- Plots -------------------------------------------------------------------
safe_plot("compare_interactions.png",
          compareInteractions(merged, show.legend = FALSE, group = 1:3) +
            compareInteractions(merged, show.legend = FALSE, group = 1:3, measure = "weight"),
          w = 8, h = 4)
for (nm in c("PMF_vs_Reactive", "PMF_vs_PrePMF")) {
  safe_plot(paste0("diff_heatmap_", nm, ".png"),
            netVisual_heatmap(merged, comparison = pairs[[nm]], measure = "weight"), w = 9, h = 8)
}
safe_plot("PMF_bubble_MK_profibrotic.png",
          netVisual_bubble(cl[["PMF"]], sources.use = "MK_profibrotic", remove.isolate = TRUE,
                           signaling = intersect(foc_pw, cl[["PMF"]]@netP$pathways)),
          w = 10, h = 8)
message("\n=== Stage 13b complete ===")
