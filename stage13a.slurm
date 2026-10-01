#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --partition=standard
#SBATCH -o logs/stage13a_%x_%j.out
#SBATCH -e logs/stage13a_%x_%j.err
# Usage: sbatch -J cc_<COND> stage13a.slurm <Reactive|PrePMF|PMF>
# 8 multisession workers each hold a copy of the data matrix -> 192G headroom.
module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"
cd "${SLURM_SUBMIT_DIR}"
conda activate scrnaseq-r
Rscript stage13a_cellchat_per_condition.R "$1"
