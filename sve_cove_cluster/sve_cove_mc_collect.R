## =============================================================================
## Collect and summarize job-array Monte Carlo output.
##   Usage (from the repo root):
##     Rscript sve_cove_cluster/sve_cove_mc_collect.R [outdir] [expected_reps]
## Writes sve_mc_results_cove.csv, sve_mc_summary_cove.csv, sve_mc_table.tex
## to sve_cove_cluster/results/, and reports any missing replicate ids.
## =============================================================================

library(here)

args <- commandArgs(trailingOnly = TRUE)
outdir        <- if (length(args) >= 1) {
  here("sve_cove_cluster", args[1])
} else {
  here("sve_cove_cluster", "mc_out_d6")
}
expected_reps <- if (length(args) >= 2) as.integer(args[2]) else 2000L

res_dir <- here("sve_cove_cluster", "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

source(here("sve_cove_cluster", "sve_cove_rieszcml.R"))

files <- sort(list.files(outdir, pattern = "^rep_\\d+\\.rds$", full.names = TRUE))
if (length(files) == 0L) stop("No replicate files found in ", outdir)
results <- do.call(rbind, lapply(files, readRDS))

ids <- sort(as.integer(sub("^rep_(\\d+)\\.rds$", "\\1", basename(files))))
missing <- setdiff(seq_len(expected_reps), ids)
cat(length(ids), "of", expected_reps, "replicates collected.\n")
if (length(missing) > 0) {
  cat("Missing replicate ids (rerun with --array=",
      paste(missing, collapse = ","), ")\n", sep = "")
}

summ <- summarize_sve_mc(results, scenario = "cove")
summ <- dplyr::arrange(summ, dplyr::desc(delta))
print(summ, digits = 3, row.names = FALSE)

write.csv(results, file.path(res_dir, "sve_mc_results_cove.csv"), row.names = FALSE)
write.csv(summ,    file.path(res_dir, "sve_mc_summary_cove.csv"), row.names = FALSE)
writeLines(as.character(mc_table(summ)), file.path(res_dir, "sve_mc_table.tex"))
cat("Wrote results to", res_dir, "\n")
