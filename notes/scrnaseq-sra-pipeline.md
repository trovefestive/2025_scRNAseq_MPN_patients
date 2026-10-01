---
name: scrnaseq-sra-pipeline
description: "Reusable HPC Slurm pipeline pattern for downloading + QC'ing 10x-style scRNA-seq FASTQs from SRA, applicable to other MPN/scRNA-seq SRA datasets"
metadata: 
  node_type: memory
  type: project
  originSessionId: 440ee9b8-b066-4c89-bca6-8eaea7c4d9cb
  modified: 2026-08-26T18:41:47.532Z
---

Pipeline built for the 2025 MPN patients scRNA-seq project (the project directory), PRJNA1070224, 19 runs of 10x-style single-cell bone marrow aspirate data. This pattern is meant to be reused for future similar datasets (other SRA-hosted 10x scRNA-seq cohorts).

**Source data:** `SraRunTable.csv` downloaded from SRA Run Selector, placed in the project root. Contains one row per SRR run with metadata (disease_stage, health_state, sex, tissue, LibrarySource=TRANSCRIPTOMIC SINGLE CELL, etc.) — this is the metadata file all downstream scripts key off of.

**Step 1 — Download (`fastq_download.slurm`):**
- Slurm array job, 1 task per SRR run (`--array=1-N`, N = row count in SraRunTable.csv), not a single sequential-loop job — lets runs download in parallel and retry individually on failure.
- `runs.txt` (list of SRR accessions) must be pre-generated before submission (`tail -n +2 SraRunTable.csv | cut -d ',' -f 1 > runs.txt`) to avoid a race condition if generated inside the array job itself.
- Uses `fasterq-dump --split-files --include-technical`. The `--include-technical` flag is essential for 10x-style single-cell data — plain `--split-files` can silently drop the barcode/UMI read, making the data unusable for Cell Ranger/STARsolo. Produces 3 files per run, **but `_1`/`_2`/`_3` do NOT map to R1/R2/I1 in the assumed order** — verify actual read lengths per dataset before trusting any fixed mapping (see correction below).

**⚠️ `_1`/`_2`/`_3` ≠ R1/R2/I1 — always verify by read length, don't assume.** For this PRJNA1070224 dataset (confirmed via FastQC on all 19 runs, consistent): `_1`=8bp=**I1** (sample index), `_2`=28bp=**R1** (16bp CB + 12bp UMI, 10x v3 chemistry), `_3`=91bp=**R2** (cDNA). This is the *opposite* of a naive `_1→R1, _2→R2, _3→I1` assumption — `fasterq-dump --include-technical`'s output file order is apparently determined by the SRA submission's spot layout, not a fixed convention, so it can vary by dataset/submitter. **Before staging files for Cell Ranger, check read lengths in the FastQC/MultiQC `multiqc_general_stats.txt` (avg_sequence_length column) for each `_N` suffix and confirm which is 8bp (I1), ~26-28bp (R1, barcode+UMI), and ~90bp+ (R2, cDNA) — never hardcode the mapping.** Cell Ranger requires exact `<sample>_S1_L001_R1_001.fastq.gz` / `_R2_` / `_I1_` naming, so build a symlink directory (e.g. `fastq_cellranger/`) mapping the *verified* suffixes to the correct Cell Ranger read-type names rather than renaming/touching the original downloaded files.
- Raw uncompressed FASTQ output is much larger than the SRA `Bytes` column suggests (~5-8x) — size the scratch quota accordingly before starting.
- Each array task gzips only its own run's files after extraction (`pigz fastq/${srr}_*.fastq`), not a bulk compress-everything-at-the-end step.
- Env: `module load miniforge; conda activate rnaseq` (env has sra-tools + pigz).

**Step 2 — QC (`run_fastqc_multiqc.slurm`):**
- Single job (not an array — FastQC parallelizes internally via `--threads`, one thread per file), submitted with `sbatch --dependency=afterok:<download_array_jobid>` so it auto-starts only once every download task succeeds, no manual babysitting required.
- FastQC on all `fastq/*.fastq.gz`, then MultiQC to aggregate into one report.
- **No trimming step** — Cell Ranger/STARsolo handle adapter/poly-A trimming internally, so QC-only is the right scope for 10x data at this stage.
- Expected/normal FastQC findings for 10x reads: R1 (barcode/UMI) and I1 (sample index) will fail "Per base sequence content" and "Sequence Duplication Levels" checks because they're short structured reads, not random biological sequence — this is not a data quality problem. Only judge library/sequencing quality from R2 (cDNA) metrics.

**Reference scripts:** both scripts were adapted from an earlier project's templates (`2018_U2AF1_S34F_GSE112174`), which used simpler single-job/sequential patterns appropriate to a smaller dataset — the array-based scaling and `--include-technical`/no-trim decisions above are the deltas that matter for scRNA-seq-scale SRA datasets specifically.

See [[mpn-scrnaseq-project-status]] for the live status of this specific run.
