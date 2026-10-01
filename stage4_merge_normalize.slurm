#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --partition=standard
#SBATCH -J "stage4_merge_normalize"
#SBATCH -o logs/stage4_merge_normalize_%j.out
#SBATCH -e logs/stage4_merge_normalize_%j.err

module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"

PROJ_DIR="${SLURM_SUBMIT_DIR}"
cd "$PROJ_DIR"
mkdir -p logs

echo "=== Stage 4: merge + normalize + HVG + scale + PCA (scrnaseq-r env) ==="
conda activate scrnaseq-r
Rscript stage4_merge_normalize.R
