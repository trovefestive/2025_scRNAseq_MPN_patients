#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --partition=standard
#SBATCH -J "stage8_cellcycle"
#SBATCH -o logs/stage8_cellcycle_%j.out
#SBATCH -e logs/stage8_cellcycle_%j.err

# Lighter than Stage 7 (42G, no SingleR here) -- 96G is headroom, not a guess.
module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs
echo "=== Stage 7 labels + Stage 8 cell cycle ==="
conda activate scrnaseq-r
Rscript stage8_cellcycle.R
