#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --partition=standard
#SBATCH -J "stage6_cluster_umap"
#SBATCH -o logs/stage6_cluster_umap_%j.out
#SBATCH -e logs/stage6_cluster_umap_%j.err

# Stage 5 peaked at 34 GB; the SNN graph over 111,861 cells plus a 6-resolution
# Louvain sweep is the memory-heavy part here, hence the larger request.

module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"

PROJ_DIR="${SLURM_SUBMIT_DIR}"
cd "$PROJ_DIR"
mkdir -p logs

echo "=== Stage 6: clustering + UMAP (scrnaseq-r env) ==="
conda activate scrnaseq-r
Rscript stage6_cluster_umap.R
