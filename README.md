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

## Main findings

All stage-level comparisons are made **per patient** (n = 19), not over pooled cells.
With 2 reactive and 3 MF2 marrows, these are descriptive results, not tests.

**1. A pro-fibrotic megakaryocyte subset that tracks fibrosis grade.**
One megakaryocyte cluster (723 cells; MK marker score 2.48 vs 0.90 endothelial)
is defined by a secretory programme: *SERPINE1, VEGFC, CXCL1, PDGFB, PDGFA, THBS1*.
Three independent analyses agree it is fibrogenic:

- marker genes (Stage 9);
- literature fibrosis signatures not derived from this data: ranked 1st of 26
  cell types on both `PMF_fibrosis` and `MK_profibrotic_secretome`, while
  unremarkable on inflammation (Stage 10);
- unsupervised GSEA over all 50 HALLMARK sets: COAGULATION (padj 2.3e-4),
  TGF-beta (0.013), ANGIOGENESIS (0.044), EMT (0.049) (Stage 11).

It is present in 8 of 9 MF3 marrows (median 0.64% of cells) and essentially absent
at MF2 (<= 0.06%), pre-fibrotic (<= 0.19%) and reactive (0%). MF2 and MF3 are both
overt PMF, so this follows fibrosis grade, not diagnosis alone.

**2. Megakaryocyte expansion with disease stage.** Median MK-lineage fraction per
patient: reactive 1.5%, prefibrotic 3.3%, MF2 4.0%, MF3 9.9%. One MF3 marrow is at
39%; excluding it the MF3 median is 9.5%.

**3. The fibrogenic shift is confined to MK subsets.** Within MK progenitors, MF3
patients score about 3x higher on fibrosis and 4-5x on the MK secretome than earlier
stages. Canonical megakaryocytes are flat across stages (0.093 / 0.107 / 0.106).

**4. Stromal cells carry the strongest fibrotic programme** (EMT NES 2.28,
padj 2.7e-19, 200 cells). Their *abundance* is lowest at MF3, but fibrotic marrow
aspirates poorly (dry tap), so that drop cannot be separated from sampling.

**5. Cell-cell communication (CellChat; Reactive / Pre-PMF / PMF = MF2 + MF3).**
Pre-PMF and PMF infer similar totals (3,976 and 4,086 interactions); Reactive infers
904, which mostly reflects 2 patients and very small cell groups, not less signaling.
Relative to Pre-PMF, PMF shows more MIF, MHC-I, ADGRE5, cyclophilin A and galectin
signaling and less MHC-II. In PMF the pro-fibrotic MK subset sends mainly
THBS1 -> CD47/CD36, MIF -> CD74 and TGFB1 -> TGF-beta receptor signals.

Two limits shape these results:

- The apparent drop in collagen, fibronectin, VCAM and CXCL signaling in PMF comes
  from stromal cells and follows how many were captured (135 Pre-PMF vs 54 PMF), so
  it reflects poor aspiration of fibrotic marrow, not less matrix signaling.
- The MK-to-stroma axis can't be tested here. The pro-fibrotic MK subset expresses
  *PDGFA*, *PDGFB* and *VEGFC* in roughly half its cells, but PMF stroma is
  under-sampled and the data contain no endothelial cells, so CellChat finds almost
  no receivers for those ligands.

CellChat pools cells within each condition; its p-values compare ligand-receptor
probabilities, not patients.
