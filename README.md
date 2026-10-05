# scRNA-seq of MPN bone marrow: myelofibrosis (PMF) cohort

Reanalysis of public single-cell RNA-seq from bone marrow aspirates of patients with
primary myelofibrosis (PMF), pre-fibrotic PMF and reactive (non-MPN) marrow.
Methods follow Jung et al. ([PMC11959246](https://pmc.ncbi.nlm.nih.gov/articles/PMC11959246/)),
adapted to this cohort.

- **Data:** BioProject PRJNA1070224, 19 SRR runs (10x 3' GEX, NovaSeq 6000)
- **Cohort:** PMF-MF3 n=9, PMF-MF2 n=3, prefibrotic n=5, reactive BM n=2 (one sample per patient)
- **Cells after QC:** 111,861 across 19 patients
- **Compute:** SLURM HPC cluster, R 4.3.3 / Seurat 5.3.0

The cohort has no ET/PV samples, so the paper's cross-MPN-subtype megakaryocyte
projection is out of scope, and the two reactive marrows stand in for normal BM.

## Methods

| Stage | What | Script | Output |
|---|---|---|---|
| 0 | Download + FastQC/MultiQC | `fastq_download.slurm`, `run_fastqc_multiqc.slurm` | not tracked |
| 1 | Cell Ranger 10.0.0 count (GRCh38-2024-A) | `cellranger_count.slurm` | not tracked |
| 2 | Seurat QC: nFeature >= 200, percent.mt <= 20 | `seurat_qc.R` | `seurat_qc/qc_summary.csv` |
| 3 | Scrublet doublet removal, per sample | `scrublet_doublets.py`, `apply_doublet_filter.R` | `seurat_qc/*summary.csv` |
| 4 | Merge, LogNormalize, 2,000 HVGs (per sample), PCA | `stage4_merge_normalize.R` | `stage4/` |
| 5 | Harmony on sample, dims 1:30 | `stage5_harmony.R` | `stage5/` |
| 6 | Clustering (6 resolutions), UMAP | `stage6_cluster_umap.R` | `stage6/` |
| 7 | SingleR (HPCA) + curated BM marker panel | `stage7_annotate.R` | `stage7/` |
| 8 | Consensus cell types, cell-cycle scoring, composition | `stage8_cellcycle.R` | `stage8/` |
| 9 | Markers per cell type (log2FC > 1, BH-FDR < 0.2) | `stage9_deg.R` | `stage9/` |
| 10 | UCell signatures (50 HALLMARK + 5 literature PMF sets) | `stage10_ucell.R` | `stage10/` |
| 11 | fgsea HALLMARK on full per-cell-type rankings | `stage11_fgsea.R` | `stage11/` |
| 12 | MK cross-subtype projection | skipped: no ET/PV in cohort | |
| 13 | CellChat per condition, then cross-condition comparison | `stage13a_cellchat_per_condition.R`, `stage13b_compare.R` | `stage13/` |
| 14 | Per-cell-type pseudobulk DE, MF3 vs prefibrotic + MF2 (edgeR QL, `~ sex + group`, patient as replicate) | `stage14_pseudobulk_de.R` | `stage14/` |

## Main findings

Statistics below are per patient (n = 19 marrows).

**1. Megakaryocytes expand with disease stage.** The megakaryocyte-lineage fraction
rises with stage (Spearman rho 0.65, p = 0.003): median 1.5% reactive, 3.3%
prefibrotic, 4.0% MF2, 9.9% MF3. It stays significant without the most MK-rich
marrow (p = 0.006).

**2. MF3 marrow contains intact megakaryocytes with a fibrogenic programme.**
One cluster (723 cells) is intact megakaryocytes: about 19,000 UMIs per cell, with
the mature MK glycoproteins GP6 and GP1BA in about 80% of cells. It isn't doublets
(0.4% co-express endothelial genes, 0% stromal). Its markers are *SERPINE1, VEGFC,
CXCL1, PDGFB, PDGFA* and *THBS1*, and three independent analyses call it fibrogenic:

- marker genes (Stage 9);
- literature fibrosis signatures not derived from this data: ranked 1st of 26
  cell types on both `PMF_fibrosis` and `MK_profibrotic_secretome`, while
  unremarkable on inflammation (Stage 10);
- unsupervised GSEA over all 50 HALLMARK sets: COAGULATION (padj 2.3e-4),
  TGF-beta (0.013), ANGIOGENESIS (0.044), EMT (0.049) (Stage 11).

Intact megakaryocytes were captured almost only from MF3 marrow (1,056 of 1,075
cells; none from reactive marrow). In MF3, 10-65% of each patient's intact
megakaryocytes carry the programme; 1 of the 19 captured from earlier stages does.
The cluster's share of cells is higher in MF3 (p = 0.003).

**3. Megakaryocyte progenitors shift toward the same programme in MF3.** MF3 MK
progenitors score higher than earlier stages on fibrosis (p = 0.002), MK secretome
(p = 0.02) and angiogenesis (p = 0.05) signatures. Platelet-like MK particles
(about 1,000 UMIs per cell) do not change across stages (0.093 / 0.107 / 0.106).

**4. MF3 progenitors express PRAME.** In the pseudobulk comparison, HSC/MPP and MEP
express the cancer-testis antigen *PRAME* in 8 of 8 MF3 patients and none of 6
earlier-stage patients; T cells from the same marrows don't express it. MF3 MEPs
also express *EGF* and *DKK1* above the highest earlier-stage level in 7 of 8
patients.

**5. Fibroblast-like stromal cells carry the strongest fibrotic programme** (EMT
NES 2.28, padj 2.7e-19). Of the 96 cells with a fibroblast programme, 94 come from
prefibrotic marrow and none from MF3.

**6. Cell-cell communication (CellChat; Reactive / Pre-PMF / PMF = MF2 + MF3).**
Pre-PMF and PMF infer similar totals (3,976 and 4,086 interactions).
Relative to Pre-PMF, PMF shows more MIF, MHC-I, ADGRE5, cyclophilin A and galectin
signaling and less MHC-II. In PMF, intact megakaryocytes send mainly
THBS1 -> CD47/CD36, MIF -> CD74 and TGFB1 -> TGF-beta receptor signals.
