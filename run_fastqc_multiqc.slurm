#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=standard
#SBATCH --job-name=fastqc_multiqc
#SBATCH -o logs/fastqc_multiqc_%j.out
#SBATCH -e logs/fastqc_multiqc_%j.err

# Load required modules – FastQC for per-file analysis, MultiQC for summary
module purge
module load fastqc          # verify exact version with `module avail fastqc`
module load multiqc           # use the exact version available on your HPC

# Project root and directories
PROJ_DIR="${SLURM_SUBMIT_DIR}"
FASTQ_DIR=$PROJ_DIR/fastq
QC_DIR=$PROJ_DIR/fastqc_reports
MULTIQC_DIR=$PROJ_DIR/multiqc_report

mkdir -p "$QC_DIR" "$MULTIQC_DIR" "$PROJ_DIR/logs"

# NOTE: this is 10x-style single-cell data (--include-technical was used when
# downloading), so each run has 3 fastq.gz files. IMPORTANT: the _1/_2/_3
# suffixes do NOT map to R1/R2/I1 in the naive order -- verified via FastQC
# read length for this dataset: _1 = I1 (8bp sample index), _2 = R1 (28bp
# cell barcode + UMI, structured), _3 = R2 (91bp cDNA, biological read).
# See fastq_cellranger/ for correctly-renamed symlinks used by Cell Ranger.
# FastQC will very likely flag I1/R1 (_1/_2) for "Per base sequence content"
# and "Sequence Duplication Levels" — that's expected for structured technical
# reads, not a sign of bad data. Judge library quality from R2 (_3) metrics.
# No adapter/quality trimming is run here: Cell Ranger / STARsolo handle
# trimming internally during alignment/quantification.

# Step 1: FastQC on all gzipped FASTQ files
echo "Running FastQC on all .fastq.gz files..."
fastqc \
  --threads 16 \
  --outdir "$QC_DIR" \
  "$FASTQ_DIR"/*.fastq.gz

# Step 2: MultiQC to combine the FastQC output
echo "Generating MultiQC report..."
multiqc \
  --force \
  --outdir "$MULTIQC_DIR" \
  "$QC_DIR"

echo "QC finished. Open $MULTIQC_DIR/multiqc_report.html in a browser."
