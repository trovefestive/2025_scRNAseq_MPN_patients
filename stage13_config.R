# Stage 13 configuration: cell-group merges + LR pathway categories.

# ---- Cell groups for CellChat (plan Stage 13 merges + earlier decisions) ----
# Plan: "merge naive T subsets -> naive T, NK1-4 -> NK, mono subsets -> monocyte".
# Stage 9: clusters 7/8 are cycle-phase splits of erythroblast -> Erythroblast.
# EXCLUDED (NA): Cycling_progenitor + Progenitor_G2M -- defined by cell-cycle
# genes, lineage not supportable (Stage 9), so their edges would be uninterpretable.
CC_GROUP <- c(
  Monocyte_classical = "Monocyte", Monocyte_nonclassical = "Monocyte",
  Monocyte_S100high = "Monocyte",
  T_CD4_naive = "T_naive", T_CD4_memory = "T_CD4_memory", T_memory = "T_memory",
  T_CD8_GZMK = "T_CD8_GZMK", NK = "NK",
  B_naive = "B_naive", Plasma_cell = "Plasma_cell", pDC = "pDC", cDC2 = "cDC2",
  Neutrophil = "Neutrophil", GMP_promyelocyte = "GMP_promyelocyte",
  HSC_MPP = "HSC_MPP", MEP = "MEP",
  Erythroblast = "Erythroblast", Erythroid_cycling = "Erythroblast",
  Erythroid_S_phase = "Erythroblast",
  Erythroid_late = "Erythroid_late", Erythroid_late_HBG = "Erythroid_late_HBG",
  Megakaryocyte_prog = "MK_prog", Megakaryocyte = "MK",
  Megakaryocyte_mature = "MK_profibrotic",
  MSC_stromal = "MSC_stromal",
  Cycling_progenitor = NA, Progenitor_G2M = NA
)

# Conditions per plan: Reactive vs Pre-PMF vs PMF (MF2 + MF3 pooled).
CONDITION <- c("Reactive" = "Reactive", "prefibrotic" = "PrePMF",
               "PMF(MF2)" = "PMF", "PMF(MF3)" = "PMF")

# Erythroid-only sample (Stage 7/8): ~76% late erythroid, would skew PMF
# erythroid averages -- excluded per standing decision 2.
EXCLUDE_SAMPLES <- c("SRR27748964")

# ---- LR pathway categories (plan: cytokine/chemokine, immune checkpoint,
# growth factor, other). CellChatDB has no such field, so this is an explicit
# pathway_name map; anything unlisted -> "Other".
# Co-stimulatory TNFSF ligands (CD40, CD70, CD137, OX40, GITRL, CD30) are filed
# under immune checkpoint (stimulatory checkpoints), not cytokine.
LR_CATEGORY <- c(
  setNames(rep("Cytokine_chemokine", 26), c(
    "CCL","CXCL","CX3C","XCR","IL1","IL2","IL4","IL6","IL10","IL12","IL16","IL17",
    "IFN-I","IFN-II","TNF","LT","LIGHT","TWEAK","TRAIL","FASLG","VEGI","CSF","CSF3",
    "LIFR","OSM","MIF")),
  setNames(rep("Cytokine_chemokine", 4), c("CHEMERIN","BAFF","APRIL","RANKL")),
  setNames(rep("Immune_checkpoint", 21), c(
    "PD-L1","PDL2","CD80","CD86","PVR","CD96","LAIR1","BTLA","CD160","VISTA",
    "CD276","CD200","SIRP","GALECTIN","ICOS","CD40","CD70","CD137","OX40","GITRL","CD30")),
  setNames(rep("Growth_factor", 31), c(
    "TGFb","BMP","BMP10","GDF","ACTIVIN","NODAL","AMH","MSTN","PDGF","VEGF","FGF",
    "EGF","EPGN","NRG","IGF","IGFBP","HGF","NGF","NT","GDNF","KIT","FLT3","THPO",
    "EPO","ANGPT","ANGPTL","PTN","MK","GAS","PROS","GRN"))
)
categorize <- function(p) ifelse(p %in% names(LR_CATEGORY), LR_CATEGORY[p], "Other")
