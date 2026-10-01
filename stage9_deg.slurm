#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=standard
#SBATCH -J "stage9_deg"
#SBATCH -o logs/stage9_deg_%j.out
#SBATCH -e logs/stage9_deg_%j.err

# presto is present, so Wilcoxon is fast (0.5 min on a 5.7k-cell subsample).
module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs
echo "=== Stage 9: DEG per cell type (all cells) ==="
conda activate scrnaseq-r
Rscript stage9_deg.R
