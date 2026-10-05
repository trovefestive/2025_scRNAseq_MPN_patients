# Is cluster 22's programme disease-associated, or just what intact MKs look like?
# Find intact-MK-like cells in EVERY sample (independent of cluster labels) and
# compare the programme across stages.
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("stage8/annotated_cellcycle.rds")
e <- GetAssayData(obj, assay = "RNA", layer = "data")
pos <- function(g) if (g %in% rownames(e)) e[g, ] > 0 else rep(FALSE, ncol(e))
mkprog <- (pos("ITGA2B") + pos("GP9") + pos("PF4") + pos("PPBP")) >= 2
mature <- pos("GP6") | pos("GP1BA")
# STRICT: all of ITGA2B, GP9, TUBB1 plus GP6 or GP1BA -- excludes cells carrying only ambient PF4/PPBP
intact <- pos("ITGA2B") & pos("GP9") & pos("TUBB1") & mature & obj$nCount_RNA >= 5000
sig    <- (pos("SERPINE1") + pos("VEGFC") + pos("PDGFB") + pos("CXCL1")) >= 2
stage  <- factor(obj$disease_stage, levels = c("Reactive","prefibrotic","PMF(MF2)","PMF(MF3)"))

cat("== STRICT intact MKs (ITGA2B+GP9+TUBB1 + GP6/GP1BA + nCount >= 5000) ==\n")
cat("total:", sum(intact), "| of which in cluster 22:", sum(intact & obj$cell_type == "Megakaryocyte_mature"), "\n")
print(table(cell_type = obj$cell_type[intact])[table(obj$cell_type[intact]) > 0])

cat("\n== By stage: intact MKs and how many carry the cluster-22 programme ==\n")
d <- data.frame(stage = stage, intact, sig)
r <- do.call(rbind, lapply(levels(stage), function(s) { x <- d[d$stage == s & d$intact, ]
  data.frame(stage = s, n_intact = nrow(x), n_with_programme = sum(x$sig),
             pct_with_programme = if (nrow(x)) round(100 * mean(x$sig), 1) else NA) }))
print(r, row.names = FALSE)

cat("\n== Per sample (intact MKs / with programme) ==\n")
ps <- tapply(seq_len(ncol(obj)), obj$sample, function(i)
  sprintf("%s %d/%d", obj$disease_stage[i[1]], sum(intact[i]), sum(intact[i] & sig[i])))
print(noquote(ps))

cat("\n== Per-gene % in intact MKs, by stage ==\n")
for (g in c("SERPINE1","VEGFC","PDGFB","PDGFA","CXCL1","TGFB1","THBS1","PF4")) {
  v <- tapply(pos(g)[intact], stage[intact], function(z) round(100 * mean(z), 1))
  cat(sprintf("%-9s %s\n", g, paste(sprintf("%s=%s", names(v), v), collapse = "  ")))
}

cat("\n== Sensitivity: MF3 without SRR27748963 ==\n")
k <- stage == "PMF(MF3)" & obj$sample != "SRR27748963" & intact
cat(sprintf("MF3 excl. 963: %d intact MKs, %d with programme (%.1f%%)\n", sum(k), sum(k & sig), 100*mean(sig[k])))
cat("\n== Per-patient % of intact MKs carrying programme (patients with >=10 intact MKs) ==\n")
pp <- do.call(rbind, lapply(split(seq_len(ncol(obj)), obj$sample), function(i) {
  n <- sum(intact[i]); data.frame(sample = obj$sample[i[1]], stage = obj$disease_stage[i[1]], n_intact = n,
    pct_prog = if (n) round(100*sum(intact[i] & sig[i])/n, 1) else NA) }))
pp <- pp[pp$n_intact >= 10, ]; print(pp[order(pp$stage, -pp$pct_prog), ], row.names = FALSE)
nm <- pp$stage == "PMF(MF3)"
if (sum(nm) >= 2 && sum(!nm) >= 2) {
  w <- wilcox.test(pp$pct_prog[nm], pp$pct_prog[!nm], exact = FALSE)
  cat(sprintf("\nPatient-level Wilcoxon, MF3 (n=%d) vs non-MF3 (n=%d): p = %.3g\n", sum(nm), sum(!nm), w$p.value))
}
