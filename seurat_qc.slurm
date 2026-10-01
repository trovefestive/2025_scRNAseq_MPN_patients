#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=standard
#SBATCH --mem=128G
#SBATCH -J "seurat_qc"
#SBATCH -o logs/seurat_qc_%j.out
#SBATCH -e logs/seurat_qc_%j.err

module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scrnaseq-r

PROJ_DIR="${SLURM_SUBMIT_DIR}"
cd "$PROJ_DIR"
mkdir -p logs seurat_qc

Rscript seurat_qc.R
