# Stage 7 marker panel, curated for MPN bone marrow.
# HPCA is a bulk primary-cell reference: good for broad immune lineages, weak on
# BM-specific compartments (MK, erythroid precursors, HSPC subsets), which is
# exactly where this project's biology lives -- hence the manual panel.
MARKERS <- list(
  HSPC        = c("CD34","CRHBP","AVP","HLF","MEIS1","ZNF521"),
  MEP         = c("GATA1","KLF1","TFRC","CNRIP1"),
  Erythroid   = c("HBB","HBD","AHSP","ALAS2","GYPA","CA1","SLC4A1","PRDX2"),
  Megakaryo   = c("PF4","PPBP","ITGA2B","GP9","GP1BB","VWF","TUBB1","CMTM5","NRGN"),
  GranProg    = c("MPO","ELANE","PRTN3","AZU1","CTSG"),
  Monocyte    = c("LYZ","CD14","VCAN","FCN1","S100A8","S100A9","MNDA","CTSS"),
  Mono_CD16   = c("FCGR3A","MS4A7","CDKN1C"),
  cDC         = c("CD1C","FCER1A","CLEC9A"),
  pDC         = c("LILRA4","IL3RA","CLEC4C"),
  B_cell      = c("MS4A1","CD79A","CD79B","PAX5","EBF1","BANK1"),
  ProB        = c("VPREB1","DNTT","IGLL1"),
  Plasma      = c("MZB1","JCHAIN","XBP1","SDC1"),
  T_cell      = c("CD3D","CD3E","CD2","IL7R"),
  T_naive     = c("CCR7","LEF1","SELL","TCF7"),
  T_CD8       = c("CD8A","CD8B","GZMK"),
  NK          = c("NKG7","GNLY","KLRD1","KLRF1","PRF1"),
  Cycling     = c("MKI67","TOP2A","TPX2","CDK1"),
  Baso_Mast   = c("MS4A2","CPA3","TPSAB1","HDC"),
  Stromal     = c("COL1A1","COL3A1","LUM","PDGFRB","CXCL12","LEPR"),
  Endothelial = c("PECAM1","CDH5")
)
