#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=standard
#SBATCH -J "stage10_ucell"
#SBATCH -o logs/stage10_ucell_%j.out
#SBATCH -e logs/stage10_ucell_%j.err

# Gene sets come from genesets/hallmark_human.rds, pre-built on the login node --
# msigdbr 25.x downloads over the network and compute nodes have no internet.
module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs
echo "=== Stage 10: UCell signature scoring ==="
conda activate scrnaseq-r
Rscript stage10_ucell.R
