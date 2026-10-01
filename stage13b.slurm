#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --partition=standard
#SBATCH -J stage13b_compare
#SBATCH -o logs/stage13b_%j.out
#SBATCH -e logs/stage13b_%j.err
module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"
cd "${SLURM_SUBMIT_DIR}"
conda activate scrnaseq-r
Rscript stage13b_compare.R
