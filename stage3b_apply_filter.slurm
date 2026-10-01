#!/bin/bash

# Stage 3b only. Stage 3a (scrublet_doublets.py) already completed successfully
# in job 18932577 -- all 19 seurat_qc/<SRR>_scrublet.csv files are on disk and
# their calls are recorded in seurat_qc/scrublet_summary.csv. Scrublet simulates
# doublets stochastically, so re-running 3a would change those calls; this script
# reruns only the R filter step that failed.

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=standard
#SBATCH -J "stage3b_apply_filter"
#SBATCH -o logs/stage3b_apply_filter_%j.out
#SBATCH -e logs/stage3b_apply_filter_%j.err

module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"

PROJ_DIR="${SLURM_SUBMIT_DIR}"
cd "$PROJ_DIR"
mkdir -p logs

echo "=== Stage 3b: apply doublet filter + QC flag (scrnaseq-r env) ==="
conda activate scrnaseq-r
Rscript apply_doublet_filter.R
