#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=standard
#SBATCH -J "stage11_fgsea"
#SBATCH -o logs/stage11_fgsea_%j.out
#SBATCH -e logs/stage11_fgsea_%j.err

module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs
echo "=== Stage 11: fgsea HALLMARK on full per-cell-type rankings ==="
conda activate scrnaseq-r
Rscript stage11_fgsea.R
