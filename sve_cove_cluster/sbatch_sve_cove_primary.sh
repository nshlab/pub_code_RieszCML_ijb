#!/bin/bash
## Primary COVE-calibrated SVE analysis (single simulated trial, 11 deltas).
##
## SEQUENTIAL BY DESIGN (SVE_WORKERS=1): with one worker the script never
## touches {future}, so there is nothing to serialize, fork, or interrupt.
## Measured sequential footprint: ~1.5 GB peak, ~15-20 min -- this is a
## one-off job and does not need parallelism. (Empirical cross-check: the
## Monte Carlo array ran this same analysis sequentially in < 6 GB per task
## on this cluster, while multisession runs OOM'd 64-120 GB allocations --
## the blow-up lives in the parallel path, not the analysis.)
##
## If you later want the ~3x speedup and have diagnosed the multisession
## issue: raise SVE_WORKERS (max useful = 5, one per outer fold), budget
## ~2 GB + workers x 1.5 GB, and leave future.globals.maxSize at its
## default -- if the 500 MB globals guard fires, that error is the
## diagnosis, not an obstacle.
#SBATCH --job-name=RieszCML_sve_primary
#SBATCH --partition=hsph
#SBATCH --account=nhejazi_lab
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

module load R/4.5.3

export SVE_WORKERS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript sve_cove_cluster/sve_cove_rieszcml.R
