library(here)

# the purpose of a meal is to digest its ingredients.
#
# we will consume the fit results from each of our jobs in a job array,
# and then run evaluation and visualization.

# dependencies
source(here("comparison_of_methods/R/01_run_on_cluster/00_dependencies.R"))

# assign job_id to be an NA_integer_
job_id <- NA_integer_
# define the sim_array_config
source(here("comparison_of_methods/R/01_run_on_cluster/01_job_array_config.R"))

# read in fit results
result_dirs <- list.files(here("comparison_of_methods/results/sim_results/"))
fit_results <- lapply(result_dirs, \(job_dir) {
  if (file.exists(here('comparison_of_methods/results/sim_results', job_dir, 'fit_results.rds'))) {
    return(readRDS(here('comparison_of_methods/results/sim_results', job_dir, 'fit_results.rds')))
  } else {
    return(NULL)
  }
})

# extract fit_results -- and expand into one giant dataframe
fit_results <- purrr::map_dfr(
  fit_results,
  \(x) {
    # append the job_id from each fit_results
    bind_cols(tibble(
      job_id     = x$job_id
    ),
    as_tibble(x$fit_results))
  }
)

fit_results$.rep <- as.integer(fit_results$.rep)

# bring in our setup data describing each job
# sim_array_config has important information describing every job.
#
# we left_join our fit_results into it by job_id so that we have
# that information available to us in evaluation and visualization.
fit_results <- dplyr::left_join(sim_array_config, fit_results, by = c('job_id' = 'job_id'))

# construct experiment ----------------------------------------------------

source(here("comparison_of_methods/R/01_run_on_cluster/02_dgp.R"))
source(here("comparison_of_methods/R/01_run_on_cluster/03_methods.R"))

# NOTE that the experiment defined here may differ from the one defined in
# R/01_run_on_cluster/04_experiment.R
#
# The reason for that is that the experiment that has run_experiment() called
# on it on the cluster may be just a subset of the whole experiment.
#
# The experiment we define here should include all of the dgps and methods
# used.
#
# This experiment is for:
#  - running the evaluators and visualizers, and
#  - rendering the documentation
#
experiment <-
  simChef::create_experiment(
    name = "RieszCML nadir::super_learner vs tmle3 sl3 Super Learner"
  ) |>
  simChef::add_dgp(dgp = nonlinear_binary_dgp) |>
  simChef::add_method(m_riesz_onestep_nadir_sl) |>
  simChef::add_method(m_riesz_tmle_nadir_sl) |>
  simChef::add_method(m_tmle3_sl)

# load evaluators and visualizers -----------------------------------------

source(here("comparison_of_methods/R/02_evaluators_and_visualizers/05_evaluators.R"))
source(here("comparison_of_methods/R/02_evaluators_and_visualizers/06_visualizers.R"))

# add evaluators and visualizers to recipe --------------------------------

# add the evaluator and visualizations to the experiment
experiment <- experiment |>
  add_evaluator(eval_ate_summary) |>
  add_visualizer(viz_coverage) |>
  add_visualizer(viz_mean_se_over_empirical_sd)

# we have to call `save = TRUE` in order for documentation to render properly
eval_results <- experiment$evaluate(fit_results)
#
# # setup and render documentation ------------------------------------------
#
# # render the experiment's documentation
# render_docs(experiment)
#
# # View the rendered file now!
# browseURL(here("results/ATE methods/ATE methods.html"))


# sample sizes appearing in your summary table
sample_size_grid <- sort(unique(eval_results$summary$sample_size))

# reference lines for panels
ref_lines_df <- tibble::tribble(
  ~measure,               ~y,
  "Coverage",             0.95,
  "Bias",                 0,
  "nMSE / Eff. Bound",    1
) |>
  tidyr::crossing(
    framework = c("{RieszCML}", "{tmle3}")
  ) |>
  mutate(
    measure = factor(
      measure,
      levels = c("Bias", "Coverage", "Mean Squared Error", "nMSE / Eff. Bound")
    ),
    framework = factor(framework, levels = c("{RieszCML}", "{tmle3}"))
  )

plot_df <- eval_results$summary |>
  filter(!stringr::str_detect(.method_name, "onestep")) |>
  mutate(
    estimator = ifelse(
      stringr::str_detect(.method_name, "tmle"),
      "TMLE",
      "1-step"
    ),
    framework = ifelse(
      stringr::str_detect(.method_name, "riesz"),
      "{RieszCML}",
      "{tmle3}"
    ),
    framework = factor(framework, levels = c("{RieszCML}", "{tmle3}")),
    nuisances = case_when(
      stringr::str_detect(.method_name, "g_misspecified") ~ "PS misspecified",
      stringr::str_detect(.method_name, "q_misspecified") ~ "OR misspecified",
      TRUE ~ "Correctly specified"
    ),
    `nMSE / Eff. Bound` = n_mse / eff_bound$var_eif
  ) |>
  tidyr::pivot_longer(
    cols = c("bias", "mse", "coverage", "nMSE / Eff. Bound"),
    names_to = "measure",
    values_to = "value"
  ) |>
  mutate(
    measure = factor(
      measure,
      levels = c("bias", "coverage", "mse", "nMSE / Eff. Bound"),
      labels = c("Bias", "Coverage", "Mean Squared Error", "nMSE / Eff. Bound")
    )
  )

plot_df |>
  ggplot(aes(
    x = sample_size,
    y = value,
    color = nuisances,
    linetype = nuisances
  )) +
  geom_hline(
    data = ref_lines_df,
    aes(yintercept = y),
    inherit.aes = FALSE,
    color = "black",
    linetype = "solid",
    linewidth = 0.8
  ) +
  geom_point(position = position_dodge(width = 400), alpha = 0.75) +
  geom_line(position = position_dodge(width = 400), alpha = 0.75) +
  theme_minimal() +
  ggh4x::facet_grid2(
    rows = vars(measure),      # metrics down the rows
    cols = vars(framework),    # {RieszCML} / {tmle3} across the columns
    scales = "free_y",
    switch = "y"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      measure == "Coverage" ~ scale_y_continuous(
        limits = c(0.8, 1),
        breaks = c(0.8, 0.85, 0.9, 0.95, 1),
        labels = label_number(accuracy = 0.01)
      ),
      measure == "Bias" ~ scale_y_continuous(
        limits = c(-.005, 0.04),
        labels = label_number(accuracy = 0.01)
      ),
      measure == "Mean Squared Error" ~ scale_y_continuous(
        n.breaks = 4,
        labels = label_number(accuracy = 0.001)
      ),
      measure == "nMSE / Eff. Bound" ~ scale_y_continuous(
        limits = c(.75, 1.26),
        labels = label_number(accuracy = 0.01)
      )
    )
  ) +
  scale_x_continuous(breaks = c(100, 1600, 3600, 6400)) +
  see::scale_color_okabeito() +
  labs(
    title = "Convergence diagnostics for {RieszCML} and Standard TMLE",
    subtitle = "Black lines indicate ideal bias, nominal coverage, and efficiency ratio = 1",
    x = "sample size",
    y = "",
    caption = paste0(
      "All experiments were run for ",
      comma_format()(unique(eval_results$summary$n_reps)[1]),
      " replications."
    )
  ) +
  theme(
    plot.caption.position = "plot",
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey90"),
    strip.placement = "outside",
    strip.text.x = element_text(size = 11, face = "bold"),
    plot.title.position = "plot"
  )

dir.create(here("comparison_of_methods/results/figures"), showWarnings = FALSE)

# portrait, letter-width: fills a text column at 6.5in wide
ggsave(here("comparison_of_methods/results/figures/convergence_diagnostics.png"),  width = 6.5, height = 8, scale = 1, dpi = 300)
ggsave(here("comparison_of_methods/results/figures/convergence_diagnostics.pdf"),  width = 6.5, height = 8, scale = 1)
ggsave(here("comparison_of_methods/results/figures/convergence_diagnostics.jpeg"), width = 6.5, height = 8, scale = 1, dpi = 300)
