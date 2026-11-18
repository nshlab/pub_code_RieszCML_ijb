#!/bin/bash
#SBATCH --job-name=RieszCML-comparison
#SBATCH --account=nhejazi_lab
#SBATCH --partition=hsph
#SBATCH --output=comparison_of_methods/slurm/logs/%x_%A_%a.out
#SBATCH --error=comparison_of_methods/slurm/logs/%x_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=24G
#SBATCH --array=1-100

# Submit from the repository root, so that .Rprofile activates renv
# and the log paths above resolve.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module load R/4.5.3

unset _R_CHECK_LIMIT_CORES_

Rscript comparison_of_methods/R/01_run_on_cluster/05_run_experiment.R "$SLURM_ARRAY_TASK_ID"
