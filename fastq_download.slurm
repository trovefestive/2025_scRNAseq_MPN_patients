#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=standard
#SBATCH --mem=128G
#SBATCH -J "fastq_download"
#SBATCH --array=1-19
#SBATCH -o logs/fastq_download_%A_%a.out
#SBATCH -e logs/fastq_download_%A_%a.err

# Load required modules
module purge
module load miniforge
conda activate rnaseq

# --- Move to the project directory ---
cd "${SLURM_SUBMIT_DIR}"

# --- runs.txt must already exist (generate it once before submitting, see below) ---
if [ ! -f runs.txt ]; then
    echo "runs.txt not found. Generate it before submitting, e.g.:"
    echo "  tail -n +2 SraRunTable.csv | cut -d ',' -f 1 > runs.txt"
    exit 1
fi

# --- Create output/log directories if they don't already exist ---
mkdir -p fastq logs

# --- Pick this array task's SRR accession ---
srr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" runs.txt)

if [ -z "$srr" ]; then
    echo "No run found for array index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Downloading ${srr} ..."

# --include-technical keeps the barcode/UMI (and index) reads that 10x-style
# single-cell libraries need for Cell Ranger / STARsolo processing, which
# plain --split-files can silently drop.
fasterq-dump "$srr" \
    --split-files \
    --include-technical \
    -e "${SLURM_CPUS_PER_TASK}" \
    -p \
    -O fastq

# --- Compress this run's FASTQ files ---
echo "Compressing FASTQ files for ${srr} ..."
pigz -p "${SLURM_CPUS_PER_TASK}" fastq/"${srr}"_*.fastq

echo "Done with ${srr}."
