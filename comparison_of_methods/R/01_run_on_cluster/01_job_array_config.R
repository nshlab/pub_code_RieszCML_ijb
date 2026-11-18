# R/01_job_array_config.R

n_jobs <- 100

n_grid <- round(cumsum(rep(sqrt(30), 10))^2, 0)[-c(3, 5, 7, 8, 9)]

# new sample size range:
n_grid <- round(cumsum(rep(sqrt(100), 10))^2, 0)[-c(3, 5, 7, 9, 10)]


sim_array_config <- data.frame(
  job_id = 1:n_jobs,
  dgp = rep("nonlinear_binary", n_jobs),
  sample_size = rep(n_grid, length.out = 100),
  n_reps_job = rep(100, n_jobs)
)

dgp_choice  <- sim_array_config[job_id, "dgp"]
sample_size <- sim_array_config[job_id, "sample_size"]
n_reps_job  <- sim_array_config[job_id, "n_reps_job"]
