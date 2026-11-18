#!/bin/bash
## Monte Carlo as a job array, M = 2000: one replicate per task
## (requeue-safe; finished ids exit in seconds). Delta grid:
## seq(-1, 1.5, 0.5), set inside sve_cove_mc_rep.R. Measured per task:
## ~3 min, ~0.8 GB peak on one core; at %100 concurrency the full array
## finishes in roughly an hour. Coverage MC error at M = 2000: about
## +/- 0.5 percentage points.
##
## REUSE: this points at the SAME mc_out_d6 directory as the M = 200 run.
## Replicate ids are the seeds (seed_base + id), so the existing
## rep_0001..rep_0200 files are exactly what tasks 1-200 would recompute --
## the requeue-skip keeps them and only ids 201-2000 do work. (Only valid
## because the settings in sve_cove_mc_rep.R are unchanged; if you change
## the grid / n / folds, use a fresh outdir.)
##
## If submission fails with "Invalid job array specification", your
## cluster's MaxArraySize is below 2001 (check:
##   scontrol show config | grep MaxArraySize
## ); submit in chunks with the offset instead, e.g.
##   sbatch --array=1-1000 sbatch_sve_cove_mc_array.sh
##   sbatch --array=1-1000 --export=ALL,SVE_REP_OFFSET=1000 sbatch_sve_cove_mc_array.sh
##
## After completion:  Rscript sve_cove_mc_collect.R mc_out_d6
## (the collector reports any missing ids in a paste-ready --array= list)
#SBATCH --job-name=RieszCML_sve_mc
#SBATCH --partition=hsph
#SBATCH --account=nhejazi_lab
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=00:30:00
## one core per task BY DESIGN: parallelism is across array tasks (see the
## plan(sequential) note in sve_cove_mc_rep.R); with one core, even an
## accidental multicore plan degenerates safely to sequential
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --array=1-2000%100

module load R/4.5.3

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript sve_cove_cluster/sve_cove_mc_rep.R "${SLURM_ARRAY_TASK_ID}" mc_out_d6
