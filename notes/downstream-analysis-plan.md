---
name: mpn-downstream-analysis-plan
description: "Staged downstream scRNA-seq analysis plan for the MPN PMF project, mirroring the methods of Jung et al. (PMC11959246)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 72c42200-abb9-4445-b438-3c6461464b83
  modified: 2026-08-26T18:30:48.216Z
---

Source paper methods (Jung et al., PMC11959246, MPN bone marrow scRNA-seq) — full methods text was pasted into the 2026-08-26 session. This is the adapted execution plan for the project directory. See [[mpn-scrnaseq-project-status]] for live progress against this plan, and [[scrnaseq-sra-pipeline]] for the upstream download/QC pipeline already completed.

## Cohort scope caveat
This BioProject (PRJNA1070224, 19 SRR runs) contains only **PMF/pre-PMF/reactive-BM** samples:
- PMF(MF2) — Overt PMF — n=3
- PMF(MF3) — Overt PMF — n=9
- prefibrotic — Pre-PMF — n=5
- reactive (health_state="Reactive", disease_stage="-") — control — n=2

No ET/PV/post-ET-MF/post-PV-MF samples exist in this cohort, so **the paper's MK-subset cross-MPN-subtype projection step (TransferData) is out of scope** here. For the CellChat normal-vs-disease comparison, the 2 Reactive-BM samples substitute for the paper's external "normal BM" reference.

## Environment (as of 2026-08-26)
- `cellranger/10.0.0` module available (paper used 3.1.0 — version noted, not expected to be material).
- `refdata-gex-GRCh38-2024-A` (10x GRCh38 3' GEX reference) downloaded + extracted to `reference/` — done.
- `scrnaseq-r` conda env (R 4.3.3) built with Seurat, Harmony, SingleR + celldex (Human Primary Cell Atlas reference), UCell, fgsea, msigdbr — done. **CellChat still needs installing** (GitHub-only package, not on conda/bioconda) before Stage 13.
- Scrublet: reused the existing `rnaseq` conda env (already had scrublet 0.2.3) rather than building a new one.

## Pipeline stages
1. **Cell Ranger count** — per-SRR array job against `refdata-gex-GRCh38-2024-A`.
2. **Seurat QC** — filter `nFeature_RNA >= 200`, `percent.mt <= 20`.
3. **Doublet removal** — Scrublet per sample, drop predicted doublets pre-merge.
4. **Normalize** — LogNormalize + `FindVariableFeatures(vst)`.
5. **Integration** — Harmony (batch = sample/patient).
6. **Clustering/UMAP** — PCA → graph clustering → UMAP on harmony embeddings.
7. **Annotation** — SingleR (Human Primary Cell Atlas via celldex) + manual marker curation for MPN BM cell types.
8. **Cell cycle scoring** — Seurat `CellCycleScoring`.
9. **DEG per cluster** — `FindMarkers`, keep `avg_log2FC > 1` & FDR < 0.2 (paper's thresholds).
10. **Gene signature scores** — UCell, gene sets TBD (MPN/fibrosis/inflammation-relevant).
11. **GSEA** — fgsea + MSigDB HALLMARK on DEG rankings per cluster.
12. ~~MK cross-subtype projection~~ — skipped, cohort doesn't include ET/PV/post-ET/post-PV-MF.
13. **CellChat** — per condition (Reactive vs Pre-PMF vs PMF); LR pairs grouped into cytokine/chemokine, immune checkpoint, growth factor, other; merge naive T subsets → naive T, NK1-4 → NK, mono subsets → monocyte; exclude groups with <10 cells; `liftCellChat` for cross-condition comparison.

## Execution mode
**Stage-by-stage with a check-in after each stage** (explicit user preference, 2026-08-26) — do not chain multiple stages' jobs without pausing to report results and get confirmation.
