#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=standard
#SBATCH -J "stage7_annotate"
#SBATCH -o logs/stage7_annotate_%j.out
#SBATCH -e logs/stage7_annotate_%j.err

# Per-cell SingleR over 111,861 cells against HPCA is the long pole here, hence
# the 12h wall time. HPCA comes from the ExperimentHub cache in ~/.cache/R --
# compute nodes have no outbound internet, so the cache must already be populated
# (done on the login node 2026-09-09).

module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs

echo "=== Stage 7: annotation (SingleR + curated markers) ==="
conda activate scrnaseq-r
Rscript stage7_annotate.R
