#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=standard
#SBATCH --mem=96G
#SBATCH -J "cellranger_count"
#SBATCH --array=1-19
#SBATCH -o logs/cellranger_count_%A_%a.out
#SBATCH -e logs/cellranger_count_%A_%a.err

# Load required modules
module purge
module load cellranger/10.0.0

# --- Move to the project directory ---
PROJ_DIR="${SLURM_SUBMIT_DIR}"
cd "$PROJ_DIR"

FASTQ_DIR="$PROJ_DIR/fastq_cellranger"
REF_DIR="$PROJ_DIR/reference/refdata-gex-GRCh38-2024-A"
OUT_DIR="$PROJ_DIR/cellranger_count"

# --- runs.txt must already exist (same file used by fastq_download.slurm) ---
if [ ! -f runs.txt ]; then
    echo "runs.txt not found. Generate it before submitting, e.g.:"
    echo "  tail -n +2 SraRunTable.csv | cut -d ',' -f 1 > runs.txt"
    exit 1
fi

mkdir -p "$OUT_DIR" logs

# --- Pick this array task's SRR accession ---
srr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" runs.txt)

if [ -z "$srr" ]; then
    echo "No run found for array index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Running cellranger count for ${srr} ..."

# NOTE on FASTQ naming: fasterq-dump's _1/_2/_3 suffixes do NOT map to
# R1/R2/I1 in the naive order for this dataset (confirmed via FastQC read
# lengths on all 19 runs) -- _1=I1 (8bp), _2=R1 (28bp barcode+UMI),
# _3=R2 (91bp cDNA). fastq_cellranger/ holds symlinks with correct
# Cell Ranger-compliant names (<SRR>_S1_L001_{R1,R2,I1}_001.fastq.gz)
# pointing at the original fastq/ files -- see [[scrnaseq-sra-pipeline]].

cd "$OUT_DIR"

cellranger count \
    --id="${srr}" \
    --transcriptome="$REF_DIR" \
    --fastqs="$FASTQ_DIR" \
    --sample="${srr}" \
    --create-bam=true \
    --localcores="${SLURM_CPUS_PER_TASK}" \
    --localmem=90

echo "Done with ${srr}."
