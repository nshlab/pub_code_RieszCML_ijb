# R/05_run_experiment.R

library(here)

args <- commandArgs(trailingOnly = TRUE)
job_id <- as.integer(args[1])

prefix_dir <- here("comparison_of_methods")
prefix_code_dir <- here("comparison_of_methods/R/01_run_on_cluster")

source(here(prefix_code_dir, "00_dependencies.R"))

n_workers <- future::availableCores()
message(paste0("available cores detected: ", n_workers))

switch(Sys.info()[["sysname"]],
       Windows = { future::plan(future::multisession, workers = n_workers) },
       Darwin  = { future::plan(future::multisession, workers = n_workers-1) },
       Linux   = { future::plan(future::multicore, workers = n_workers) })

source(here(prefix_code_dir, "01_job_array_config.R"))
source(here(prefix_code_dir, "02_dgp.R"))
source(here(prefix_code_dir, "03_methods.R"))
source(here(prefix_code_dir, "04_experiment.R"))

# distinct reproducible rng per job
set.seed(20261026L + job_id)

out_dir <- file.path(
  here(prefix_dir, "results"),
  "sim_results",
  sprintf("job_%04d", job_id)
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

options(future.globals.maxSize = Inf)

fit_results <- simChef::fit_experiment(
  experiment,
  n_reps = n_reps_job,

  future.globals = c(
    "expit",
    "f_g",
    "f_q",
    "g_fun",
    "q_fun",
    ".compute_truth",
    "truth_vec",
    "simulate_nonlinear_binary_data",
    "nonlinear_binary_dgp",

    ".formula_Q",
    ".formula_g",
    ".g_bound",
    ".truncate_g",
    ".fit_nadir_sl",
    ".Q_preds",
    ".make_nadir_binary_learners",
    ".make_nadir_mean_only_learner",
    ".make_sl3_binary_stack",
    ".make_sl3_superlearner",
    ".make_sl3_mean_learner",

    ".make_nuis_ate_nadir_sl",
    ".make_nuis_ate_nadir_q_misspecified",
    ".make_nuis_ate_nadir_g_misspecified",

    ".soft_try",
    ".result_from_fit",
    ".run_tmle3_ate_sl",
    ".run_tmle3_ate_custom",

    "riesz_onestep_nadir_sl_fun",
    "riesz_tmle_nadir_sl_fun",
    "tmle3_sl_fun",
    'fit_tmle3',

    "riesz_onestep_nadir_q_misspecified_fun",
    "riesz_onestep_nadir_g_misspecified_fun",
    "riesz_tmle_nadir_q_misspecified_fun",
    "riesz_tmle_nadir_g_misspecified_fun",
    "tmle3_q_misspecified_fun",
    "tmle3_g_misspecified_fun",

    "m_riesz_onestep_nadir_sl",
    "m_riesz_tmle_nadir_sl",
    "m_tmle3_sl",

    "m_riesz_onestep_nadir_q_misspecified",
    "m_riesz_onestep_nadir_g_misspecified",
    "m_riesz_tmle_nadir_q_misspecified",
    "m_riesz_tmle_nadir_g_misspecified",
    "m_tmle3_q_misspecified",
    "m_tmle3_g_misspecified",

    "RieszCurve",
    "RieszFit",
    "riesz_curve_catalog",
    "riesz_estimate",
    "riesz_tmle",
    "variable_type",

    "lubridate::time_length",
    "lubridate::interval"
  )
)

saveRDS(
  list(
    job_id = job_id,
    n_reps_job = n_reps_job,
    fit_results = fit_results
  ),
  file = file.path(out_dir, "fit_results.rds")
)
