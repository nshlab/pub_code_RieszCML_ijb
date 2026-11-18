# R/02_evaluators_and_visualizers/05_evaluators.R

# Evaluator: summarize coverage, bias, sqrt(n)*bias, etc.
eval_ate_summary <- create_evaluator(
  .name = "summary",
  .eval_fun = function(fit_results, vary_params = NULL, ...) {

    fit_results %>%
      group_by(.dgp_name, .method_name, sample_size) %>%
      summarize(
        n_reps = dplyr::n(),
        coverage = mean(truth <= ci_high & truth >= ci_low, na.rm = TRUE),
        bias = mean(estimate - truth, na.rm = TRUE),
        abs_bias = mean(abs(estimate - truth), na.rm=TRUE),
        est_var = mean(se^2, na.rm = TRUE),
        mean_error = mean(estimate - truth, na.rm=TRUE),
        sqrt_n_bias = mean(sqrt(sample_size) * (estimate - truth), na.rm = TRUE),
        mse   = mean((estimate - truth)^2, na.rm = TRUE),
        n_mse = first(sample_size) * mse,
        .groups = "drop"
      )
  }
)
