## =============================================================================
## Run a single Monte Carlo replicate of the COVE-calibrated SVE simulation, for
## use with a SLURM job array (see sbatch_sve_cove_mc_array.sh).
##
##   Usage: Rscript sve_cove_mc_rep.R [rep_id] [outdir]
##
## =============================================================================

library(here)

if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  stop("renv is not active: run from the repository root so .Rprofile is read.")
}

outdir <- if (length(args) >= 2) {
  here("sve_cove_cluster", args[2])
} else {
  here("sve_cove_cluster", "mc_out_d6")
}

args <- commandArgs(trailingOnly = TRUE)
rep_id <- if (length(args) >= 1) as.integer(args[1]) else
  as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (!isTRUE(is.finite(rep_id))) {
  stop("Provide a replicate id (argument 1) or run under a SLURM job array.")
}
## Optional offset for clusters whose MaxArraySize is smaller than the target
## number of replicates: submit several --array=1-K chunks with
## SVE_REP_OFFSET=0, K, 2K, ... and the effective replicate id (and hence the
## seed and output file) becomes task_id + offset.
offset <- suppressWarnings(as.integer(Sys.getenv("SVE_REP_OFFSET", "0")))
if (isTRUE(is.finite(offset)) && offset > 0L) rep_id <- rep_id + offset
outdir <- if (length(args) >= 2) args[2] else "mc_out_d6"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

outfile <- file.path(outdir, sprintf("rep_%04d.rds", rep_id))
if (file.exists(outfile)) {
  message("Output for replicate ", rep_id, " already exists; nothing to do.")
  quit(save = "no")
}

## ---- simulation settings (edit here, identically for all array tasks) -------
n              <- 28000
scenario       <- "cove"
delta_grid     <- seq(-1, 1.5, by = 0.5)
n_folds        <- 3
inner_n_folds  <- 3
seed_base      <- 3000

# function definitions
source(here("sve_cove_cluster", "sve_cove_rieszcml.R"))

options(future.globals.maxSize = +Inf)

## IMPORTANT: keep this `sequential`. In the job-array design, parallelism
## lives ACROSS array tasks (the --array=...%K throttle), not within them:
## one task = one replicate = one core. Do NOT set plan(multicore) here --
## forked R workers un-share copy-on-write memory during garbage collection,
## transiently multiplying resident memory several-fold. The cgroup OOM
## signature of that mistake: tasks killed in ~40 s with TotalCPU > Elapsed
## while sacct's sampled MaxRSS still reads below the request. Under the
## sequential plan a replicate peaks at ~0.75 GB, and future.globals.maxSize
## is irrelevant (no globals are exported anywhere).
future::plan(future::sequential)

set.seed(seed_base + rep_id)
df  <- sim_data(n, scenario)
res <- sve_analysis(df, delta_grid = delta_grid,
                    n_folds = n_folds, inner_n_folds = inner_n_folds)
attr(res, "support") <- NULL      # keep per-replicate files small
res$rep <- rep_id

## write atomically so an interrupted task never leaves a partial file
tmp <- paste0(outfile, ".tmp")
saveRDS(res, tmp)
file.rename(tmp, outfile)
cat("replicate", rep_id, "written to", outfile, "\n")
