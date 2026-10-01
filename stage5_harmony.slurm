#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --partition=standard
#SBATCH -J "stage5_harmony"
#SBATCH -o logs/stage5_harmony_%j.out
#SBATCH -e logs/stage5_harmony_%j.err

# Stage 4 peaked at 17.8 GB on the same object, so 96G is headroom for
# RunHarmony + JoinLayers rather than a guess.

module purge
module load miniforge/26.3.2
source "$(conda info --base)/etc/profile.d/conda.sh"

PROJ_DIR="${SLURM_SUBMIT_DIR}"
cd "$PROJ_DIR"
mkdir -p logs

echo "=== Stage 5: Harmony integration (scrnaseq-r env) ==="
conda activate scrnaseq-r
Rscript stage5_harmony.R
