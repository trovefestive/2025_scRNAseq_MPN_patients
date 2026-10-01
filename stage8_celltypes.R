# Consensus cell type labels for RNA_snn_res.0.8 (27 clusters).
# Source: SingleR (HPCA, cluster + per-cell) reconciled against the curated marker
# panel. Where they disagree, MARKERS WIN -- HPCA is a bulk primary-cell reference
# and is unreliable on BM-specific compartments. Disagreements are noted inline.
CELLTYPE <- c(
  "0"  = "Monocyte_classical",   # VCAN, CD14, CD163, NLRP3
  "1"  = "T_CD4_naive",          # LEF1, NELL2, TSHZ2
  "2"  = "B_naive",              # MS4A1, PAX5, FCRL1, TNFRSF13C
  "3"  = "T_CD8_GZMK",           # GZMK, CD8A/B, GZMH
  "4"  = "NK",                   # KLRF1, NCR1, NCAM1, GNLY
  "5"  = "HSC_MPP",              # AVP, HLF, CRHBP, BAALC -- SingleR said "CMP";
                                 # AVP/HLF/CRHBP are canonical HSC, marker wins
  "6"  = "Erythroid_late",       # HBA1/2, HBB, ALAS2, HBM -- SingleR said "BM"
  "7"  = "Erythroid_cycling",    # AURKA, CCNB1/2, PTTG1 + GYPA
  "8"  = "Erythroid_S_phase",    # replication histones H2AC/H3C, E2F7/E2F8
  "9"  = "Megakaryocyte_prog",   # ITGA2B, CMTM5, MYL9, PDLIM1 -- SingleR "MEP" 32.7%
  "10" = "MEP",                  # CNRIP1, XACT, ST6GAL2, RYR3
  "11" = "T_memory",             # THEMIS, CAMK4, CD6, ANK3
  "12" = "Erythroblast",         # SOX6, TMCC2, ABCC13, BEST3
  "13" = "Megakaryocyte",        # PF4, PPBP, GP9, TUBB1, CLEC1B, ACRBP
  "14" = "Erythroid_late_HBG",   # HBG1/HBG2, TRIM58, BPGM, ALAS2 (fetal Hb)
  "15" = "Monocyte_nonclassical",# CDKN1C, HES4, PELATON (CD16+)
  "16" = "T_CD4_memory",         # ICOS, CD28, IL7R, TRAT1
  "17" = "Neutrophil",           # CYP4F3, MGAM, AQP9, PADI4, FCAR
  "18" = "Progenitor_G2M",       # CDC20, PLK1, UBE2C, AURKB -- cycling, lineage unclear
  "19" = "Monocyte_S100high",    # S100A8/A9/A12, CSTA, SERPINA1
  "20" = "GMP_promyelocyte",     # MPO, ELANE, PRTN3, AZU1, CTSG, MS4A3
  "21" = "cDC2",                 # CD1C, CD1E, FCER1A, CLEC10A
  "22" = "Megakaryocyte_mature", # PDE3A, SERPINE1, VEGFC, CXCL1 -- LEAST CONFIDENT
  "23" = "pDC",                  # LILRA4, CLEC4C, SCT -- SingleR said "Pro-B_cell_CD34+",
                                 # which is plainly wrong; these are canonical pDC genes
  "24" = "Plasma_cell",          # IGHG1, IGHA1, SDC1, IGLL5
  "25" = "MSC_stromal",          # COL1A2, COL3A1, CXCL12, DCN, CDH11, LIFR
  "26" = "Cycling_progenitor"    # E2F1, MELK, CLSPN, PKMYT1
)
# Coarse grouping for composition stats and CellChat (Stage 13).
LINEAGE <- c(
  Monocyte_classical="Myeloid", Monocyte_nonclassical="Myeloid", Monocyte_S100high="Myeloid",
  Neutrophil="Myeloid", GMP_promyelocyte="Myeloid", cDC2="Myeloid", pDC="Myeloid",
  Erythroid_late="Erythroid", Erythroid_cycling="Erythroid", Erythroid_S_phase="Erythroid",
  Erythroblast="Erythroid", Erythroid_late_HBG="Erythroid", MEP="Erythroid",
  Megakaryocyte="Megakaryocyte", Megakaryocyte_prog="Megakaryocyte",
  Megakaryocyte_mature="Megakaryocyte",
  T_CD4_naive="T_NK", T_CD8_GZMK="T_NK", T_memory="T_NK", T_CD4_memory="T_NK", NK="T_NK",
  B_naive="B_Plasma", Plasma_cell="B_Plasma",
  HSC_MPP="HSPC", Progenitor_G2M="HSPC", Cycling_progenitor="HSPC",
  MSC_stromal="Stromal"
)
