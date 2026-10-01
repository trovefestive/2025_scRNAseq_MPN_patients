# Custom gene sets for Stage 10, on top of the 50 MSigDB HALLMARK sets.
#
# DELIBERATELY LITERATURE-BASED, NOT DERIVED FROM OUR OWN STAGE 9 DEG.
# Scoring a signature built from cluster 22's own marker list against cluster 22
# would be circular and guarantee a "hit". These are canonical PMF / marrow-fibrosis
# effectors from the myelofibrosis literature, so cluster 22 scoring high on them is
# an independent test of the Stage 9 observation rather than a restatement of it.
CUSTOM_SETS <- list(
  # Core fibrogenic effectors driving marrow fibrosis in PMF
  PMF_fibrosis = c("TGFB1","TGFB2","TGFB3","PDGFA","PDGFB","FGF2","SERPINE1","THBS1",
                   "LOX","LOXL2","TIMP1","MMP9","SPP1","CTGF","CCN2","FN1",
                   "COL1A1","COL1A2","COL3A1","INHBA","BMP4"),
  # MK-derived secretory factors implicated in fibrotic niche remodelling
  MK_profibrotic_secretome = c("PF4","PPBP","TGFB1","PDGFA","PDGFB","VEGFA","VEGFC",
                               "SERPINE1","THBS1","CXCL1","CXCL5","IGF1","LIF","INHBA"),
  # Inflammatory cytokine milieu characteristic of MPN
  MPN_inflammation = c("IL1B","IL6","TNF","CXCL8","CXCL1","CXCL2","CCL2","CCL3","CCL4",
                       "OSM","LIF","IL1RN","S100A8","S100A9","NLRP3","PTGS2"),
  # JAK-STAT output -- the pathway constitutively active in MPN
  JAK_STAT_targets = c("SOCS1","SOCS2","SOCS3","CISH","PIM1","MYC","BCL2L1","STAT1",
                       "STAT3","STAT5A","STAT5B","IRF1","OSMR","IL2RA"),
  # Niche/stromal support -- CXCL12-abundant reticular (CAR) cell program
  Stromal_niche = c("CXCL12","KITLG","ANGPT1","LEPR","VCAM1","SPP1","NES","GREM1",
                    "WNT5A","IGF1","COL1A1","LIFR")
)
