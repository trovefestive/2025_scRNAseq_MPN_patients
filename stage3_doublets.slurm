#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=standard
#SBATCH --mem=64G
#SBATCH -J "stage3_doublets"
#SBATCH -o logs/stage3_doublets_%j.out
#SBATCH -e logs/stage3_doublets_%j.err

module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"

PROJ_DIR="${SLURM_SUBMIT_DIR}"
cd "$PROJ_DIR"
mkdir -p logs

echo "=== Stage 3a: Scrublet doublet detection (rnaseq env) ==="
conda activate rnaseq
python3 scrublet_doublets.py
rc=$?
if [ $rc -ne 0 ]; then
    echo "scrublet_doublets.py failed with exit code $rc"
    exit $rc
fi
conda deactivate

echo "=== Stage 3b: apply doublet filter + QC flag (scrnaseq-r env) ==="
conda activate scrnaseq-r
Rscript apply_doublet_filter.R
