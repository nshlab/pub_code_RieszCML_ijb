# R/06_visualizers.R

# visualization specific dependencies
library(ggh4x)
library(scales)
library(ggplot2)

visualize_coverage_by_n <- function(fit_results, eval_results, ...) {
  fit_results |>
    filter(! is.na(method)) |>
    # filter(method %in% c('riesz_onestep',
    #                      'riesz_tmle',
    #                      'tmle3')) |>
    mutate(covers = truth >= ci_low & truth <= ci_high) |>
    group_by(method, sample_size) |>
    summarize(
      coverage = mean(covers, na.rm=TRUE)) |>
    ggplot(aes(x = sample_size, y = coverage)) +
    geom_hline(yintercept = 0.95, linetype = 'solid', color = 'firebrick') +
    geom_point(alpha= .7) +
    geom_line() +
    theme_bw() +
    facet_wrap(~ method) +
    ylim(c(0, 1)) +
    labs(
      x = "Sample Size",
      y = "Coverage",
      title = 'Coverage of estimators vs. nominal 95% level',
      subtitle = 'Comparing the {RieszCML} one-step and TMLE vs. {tmle3}',
      caption = paste0('Each experiment (per sample size) was repeated ', unique(eval_results$summary$n_reps)[1], ' times.'))
}


viz_coverage <- create_visualizer(
  .viz_fun = visualize_coverage_by_n,
  .name = 'coverage')

# 1) Bias vs sample size, faceted by DGP
viz_bias_by_n <- create_visualizer(
  .name = "bias_by_n",
  .viz_fun = function(fit_results, eval_results, vary_params = NULL, ...) {

    dat <- NULL
    if (!is.null(eval_results) && "summary" %in% names(eval_results)) {
      dat <- eval_results[["summary"]]
    }
    if (is.null(dat) || nrow(dat) == 0) {
      dat <- fit_results %>%
        group_by(.dgp_name, .method_name, sample_size) %>%
        summarize(
          bias = mean(estimate - truth, na.rm = TRUE),
          .groups = "drop"
        )
    }

    x_breaks <- sort(unique(dat$sample_size))

    dat |>
      filter(sample_size > 30) |>
      filter(abs(bias) < 2) |>
      filter(! is.na(.method_name)) |>
      ggplot(aes(x = sample_size, y = bias,
                 color = .method_name)) +
      geom_hline(yintercept = 0, linetype = "solid", alpha = 0.8) +
      geom_point(size = 2, shape=7) +
      geom_line(aes(group = .method_name), linewidth = .75) +
      facet_wrap(~ .dgp_name, scales = 'free') +
      theme_bw() +
      theme(panel.grid.minor = element_blank()) +
      facet_wrap(~.method_name) +
      labs(
        x = "Sample Size",
        y = "Bias",
        color = "Method",
        shape = "Method",
        title = "Bias of ATE estimators across sample sizes and DGPs"
      )
  },
  .description = "
### Bias vs sample size

Average bias \\(\\mathbb{E}[\\hat\\tau - \\tau]\\) by method and sample size, faceted by DGP.
Dashed horizontal line at 0.
",
  .doc_show = TRUE,
  .doc_options = list(width = 10, height = 6)
)


# 2) sqrt(n) * bias vs sample size, faceted by DGP
viz_sqrtn_bias_by_n <- create_visualizer(
  .name = "sqrtn_bias_by_n",
  .viz_fun = function(fit_results, eval_results, vary_params = NULL, ...) {

    dat <- NULL
    if (!is.null(eval_results) && "summary" %in% names(eval_results)) {
      dat <- eval_results[["summary"]]
    }
    if (is.null(dat) || nrow(dat) == 0) {
      dat <- fit_results %>%
        filter(! is.na(.method_name)) |>
        group_by(.dgp_name, .method_name, sample_size) %>%
        summarize(
          sqrt_n_bias = mean(sqrt(sample_size) * (estimate - truth), na.rm = TRUE),
          .groups = "drop"
        )
    }

    x_breaks <- sort(unique(dat$sample_size))

    dat |>
      filter(sample_size > 30) |>
    ggplot(aes(x = sample_size, y = sqrt_n_bias,
                    color = .method_name)) +
      geom_hline(yintercept = 0, linetype = "solid", alpha = .8) +
      geom_point(size = 2, shape = 7) +
      geom_line(aes(group = .method_name), linewidth = 0.75) +
      facet_wrap(~ .method_name) +
      theme_bw() +
      labs(
        x = "Sample Size",
        y = expression(sqrt(n) %*% Bias),
        color = "Method",
        shape = "Method",
        title = expression(paste("Average ", sqrt(n), " * error across sample sizes and DGPs"))
      )
  },
  .description = "
### \\(\\sqrt{n}\\)·bias vs sample size

Average scaled bias \\(\\mathbb{E}[\\sqrt{n}(\\hat\\tau - \\tau)]\\) by method and sample size, faceted by DGP.
Dashed horizontal line at 0.
",
  .doc_show = TRUE,
  .doc_options = list(width = 10, height = 6)
)


# sqrt(n) * MSE vs sample size, faceted by method
viz_sqrtn_mse_by_n <- create_visualizer(
  .name = "sqrtn_bias_by_n",
  .viz_fun = function(fit_results, eval_results, vary_params = NULL, ...) {

    dat <- NULL
    if (!is.null(eval_results) && "summary" %in% names(eval_results)) {
      dat <- eval_results[["summary"]]
    }
    if (is.null(dat) || nrow(dat) == 0) {
      dat <- fit_results %>%
        filter(! is.na(.method_name)) |>
        group_by(.dgp_name, .method_name, sample_size) %>%
        summarize(
          sqrt_n_mse = mean(sqrt(sample_size) * (estimate - truth)^2, na.rm = TRUE),
          .groups = "drop"
        )
    }

    dat |>
      filter(sample_size > 30) |>
      ggplot(aes(x = sample_size, y = sqrt_n_mse)) +
      geom_point() +
      geom_line(aes(group = .method_name), linewidth = 0.75) +
      facet_wrap(~ .method_name) +
      ylim(c(0, NA)) +
      theme_bw() +
      labs(
        x = "Sample Size",
        y = expression(sqrt(n) %*% MSE),
        title = expression(paste(sqrt(n) %*% MSE, " across sample sizes")),
        caption = paste0('Each experiment (per sample size) was repeated ', unique(eval_results$ate_summary$n_reps)[1], ' times.')
      )
  },
  .description = "
### \\(\\sqrt{n}\\)·MSE vs sample size

Average scaled MSE \\(\\mathbb{E}[\\sqrt{n}(\\hat\\tau - \\tau)^2]\\) by method and sample size, faceted by DGP.
Dashed horizontal line at 0.
",
  .doc_show = TRUE,
  .doc_options = list(width = 10, height = 6)
)




viz_mse_by_n <- create_visualizer(
  .name = "mse_by_n",
  .viz_fun = function(fit_results, eval_results, vary_params = NULL, ...) {

    dat <- NULL
    if (!is.null(eval_results) && "summary" %in% names(eval_results)) {
      dat <- eval_results[["summary"]]
    }
    if (is.null(dat) || nrow(dat) == 0) {
      dat <- fit_results %>%
        filter(! is.na(.method_name)) |>
        group_by(.dgp_name, .method_name, sample_size) %>%
        summarize(
          mse = mean((estimate - truth)^2, na.rm = TRUE),
          .groups = "drop"
        )
    }

    dat |>
      filter(sample_size > 30) |>
      ggplot(aes(x = sample_size, y = mse)) +
      geom_point(size = 2) +
      geom_line(aes(group = .method_name), linewidth = 0.75) +
      facet_wrap(~ .method_name) +
      ylim(c(0, NA)) +
      theme_bw() +
      labs(
        x = "Sample Size",
        y = expression( MSE),
        title = expression(paste(MSE, " across sample sizes")),
        caption = paste0('Each experiment (per sample size) was repeated ', unique(eval_results$ate_summary$n_reps)[1], ' times.')
      )
  },
  .description = "
### MSE vs sample size

MSE \\(\\mathbb{E}[(\\hat\\tau - \\tau)^2]\\) by method and sample size, faceted by DGP.
",
  .doc_show = TRUE,
  .doc_options = list(width = 10, height = 6)
)




visualize_sqrtn_error_histograms <- function(fit_results, eval_results, vary_params = NULL,
                                             dgp_to_plot = "Simple DGP", ...) {

  stopifnot(is.data.frame(fit_results))
  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    stop("Package 'ggh4x' is required for this visualizer (facet_grid2).")
  }

  dat <- fit_results %>%
    mutate(
      sqrt_n_error = sqrt(sample_size) * (estimate - truth)
    ) %>%
    filter(!is.na(sqrt_n_error)) |>
    filter(sample_size > 30)

  # Make sample_size an ordered factor for nicer facet ordering
  ss_levels <- sort(unique(dat$sample_size))
  dat <- dat %>%
    mutate(
      sample_size = factor(
        sample_size,
        levels = ss_levels,
        labels = paste0("sample size: ", ss_levels)
      )
    )

  dat |>
    ggplot(aes(x = sqrt_n_error, fill = .method_name)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
    geom_density(alpha = 0.25) +
    ggh4x::facet_grid2(
      .method_name ~ sample_size,
      scales = "free",
      independent = "x",
    ) +
    scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
    theme_bw() +
    labs(
      x = expression(sqrt(n) %*% Error),
      y = "Density",
      fill = "Method/Estimator",
      title = expression(paste("Limiting distributions of ", sqrt(n) %*% Error)),
      subtitle = paste0(
        "Rows correspond to estimators; columns to increasing sample sizes. "
      )
    )



  dat |>
    filter(! is.na(eps)) |>
    ggplot(aes(x = eps, fill = .method_name)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
    geom_density(alpha = 0.25) +
    ggh4x::facet_grid2(
      .method_name ~ sample_size,
      scales = "free",
      independent = "x",
    ) +
    scale_x_continuous(labels = scales::number_format(accuracy = .01)) +
    theme_bw() +
    labs(
      x = "Epsilon",
      y = "Density",
      fill = "Method/Estimator",
      title = expression(paste("Limiting distributions of Epsilon from Fluctuation Step")),
      subtitle = paste0(
        "Rows correspond to estimators; columns to increasing sample sizes. "
      )
    )
}


viz_mean_se_over_empirical_sd <- create_visualizer(
  .viz_fun = function(fit_results, eval_results, ...) {
    empirical_sd_and_mean_se <- fit_results |>
      group_by(method, sample_size) |>
      summarize(
        mean_se = mean(se, na.rm=TRUE),
        empirical_estimate_sd = sd(estimate, na.rm=TRUE))

    empirical_sd_and_mean_se |>
      filter(! is.na(method)) |>
      filter(sample_size > 30) |>
      tidyr::pivot_longer(
        cols = c(mean_se, empirical_estimate_sd),
        names_to = 'measure') |>
      ggplot(aes(x = sample_size, y = value, linetype = measure, shape = measure)) +
      geom_point() +
      geom_line() +
      facet_wrap(~method) +
      theme_minimal() +
      labs(title = 'Comparing Empirical SD of Estimates vs. Mean(Standard Error)')
  },
  .name = 'Mean SE / Empirical SD(Estimates)'
)



# multipanel figure of convergence diagnostics ----------------------------


viz_multipanel_convergence_diagnostics <- create_visualizer(
  .viz_fun = function(fit_results, eval_results, ...) {

    # sample sizes appearing in your summary table
    sample_size_grid <- sort(unique(eval_results$summary$sample_size))

    eff_bound_df <- tibble::tibble(
      sample_size = sample_size_grid
    ) |>
      dplyr::mutate(
        # `Mean Squared Error` = eff_bound$var_eif / sample_size,
        `n * MSE`      = eff_bound$var_eif
      ) |>
      tidyr::pivot_longer(
        cols = c(`n * MSE`),
        names_to = "measure",
        values_to = "value"
      ) |>
      tidyr::crossing(
        framework = c("{RieszCML}", "Standard")
      )

    eval_results$summary |>
      filter() |>
      filter( ! stringr::str_detect(.method_name, 'onestep')) |>
      tidyr::pivot_longer(
        cols = c(
          'bias', 'mse', 'n_mse',
          'coverage'),
        values_to = 'value',
        names_to = 'measure') |>
      mutate(
        estimator = ifelse(
          stringr::str_detect(.method_name, 'tmle'), 'TMLE', '1-step'),
        framework = ifelse(
          stringr::str_detect(.method_name, 'riesz'), '{RieszCML}', 'Standard'),
        nuisances = case_when(
          stringr::str_detect(.method_name, 'g_misspecified') ~ 'PS misspecified',
          stringr::str_detect(.method_name, 'q_misspecified') ~ 'OR misspecified',
          TRUE ~ 'Correctly specified')) |>
      mutate(
        measure = factor(
          measure,
          levels = c("bias", "coverage", "mse", "n_mse"),
          labels = c("Bias", "Coverage", 'Mean Squared Error', 'n * MSE')
        )) |>
      ggplot(aes(
        y = value, x = sample_size,
        # shape = .method_name,
        color = nuisances, # interaction(nuisances, framework),
        # linewidth = framework,
        # shape = framework,
        linetype = nuisances
        )) +
      # ref line for coverage
      geom_hline(
        data = data.frame(
          measure = 'Coverage',
          y = 0.95),
        mapping = aes(yintercept = y),
        linetype = 'solid',
        color = 'black') +
      # ref line for bias
      geom_hline(
        data = data.frame(
          measure = 'Bias',
          y = 0.),
        mapping = aes(yintercept = y),
        linetype = 'solid',
        color = 'black') +
      # efficiency-bound lines for MSE panels
      geom_line(
        data = eff_bound_df,
        aes(x = sample_size, y = value),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.8,
        linetype = "solid"
      ) +
      geom_point(alpha = 0.75) +
      geom_line(alpha = 0.75) +
      theme_minimal() +
      # facet_grid(framework~measure, scales = 'free') +
      ggh4x::facet_grid2(
        cols = vars(measure),
        rows = vars(framework),
        scales = "free_y",
        independent = "y",
        switch = 'y'
      ) +
      facetted_pos_scales(
        y = list(
          measure == "Coverage" ~ scale_y_continuous(
            limits = c(0.8, 1),
            breaks = c(0.8, 0.85, 0.9, 0.95, 1),
            labels = label_number(accuracy = 0.01)
          ),
          measure == "Bias" ~ scale_y_continuous(
            limits = c(
              0,
              .055
            ),
            labels = label_number(accuracy = 0.01)
          ),
          measure == "Mean Squared Error" ~ scale_y_continuous(
            # limits = c(0, 0.05),
            n.breaks = 4,
            labels = label_number(accuracy = 0.01)
          ),
          measure == "n * MSE" ~ scale_y_continuous(
            # limits = c(.95, 1.55),
            labels = label_number(accuracy = 0.01)
          )
        )
      ) +
      see::scale_color_okabeito() +
      # guides(color = guide_none()) +
      scale_linewidth_manual(values = c(`Standard` =.5, `{RieszCML}` = .75)) +
      labs(
        title = "Convergence diagnostics for {RieszCML} TMLE and standard software",
        subtitle = "Black lines indicate ideal bias, coverage, and efficiency bound.",
        x = 'sample size',
        y = '',
        caption = paste0('All experiments were run for ', comma_format()(unique(eval_results$summary$n_reps)[1]), ' replications.')
      ) +
      theme(plot.caption.position = 'plot') +
      # guides(color = guide_legend(nrow = 2),
      #        linetype = guide_legend(nrow = 2)) +
      theme(legend.position = 'bottom',
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = 'grey90'),
            # strip.text.y.left = element_text(angle = 0),
            strip.placement = "outside",
            plot.title.position = 'plot')


# new vers ----------------------------------------------------------------

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
        framework = c("{RieszCML}", "Standard")
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
          "Standard"
        ),
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
      geom_point(alpha = 0.75) +
      geom_line(alpha = 0.75) +
      theme_minimal() +
      ggh4x::facet_grid2(
        cols = vars(measure),
        rows = vars(framework),
        scales = "free_y",
        independent = "y",
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
            limits = c(0, 0.055),
            labels = label_number(accuracy = 0.01)
          ),
          measure == "Mean Squared Error" ~ scale_y_continuous(
            n.breaks = 4,
            labels = label_number(accuracy = 0.01)
          ),
          measure == "nMSE / Eff. Bound" ~ scale_y_continuous(
            labels = label_number(accuracy = 0.01)
          )
        )
      ) +
      scale_x_continuous(breaks = c(100, 1600, 3600, 6400)) +
      see::scale_color_okabeito() +
      labs(
        title = "Convergence diagnostics for {RieszCML} TMLE and standard software",
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
        plot.title.position = "plot"
      )



    ggsave(here("comparison_of_methods/results/figures/convergence_diagnostics.png"), width=8, height=4, scale=1)
    ggsave(here("comparison_of_methods/results/figures/convergence_diagnostics.pdf"), width=8, height=4, scale=1)
    ggsave(here("comparison_of_methods/results/figures/convergence_diagnostics.jpeg"), width=8, height=4, scale=1)

    },
  .name = 'Mean SE / Empirical SD(Estimates)'
)




viz_sqrtn_error_hist_simple_dgp <- create_visualizer(
  .name = "sqrtn_error_hist_simple_dgp",
  .viz_fun = visualize_sqrtn_error_histograms,
  .doc_options = list(width = 12, height = 8)
)

viz_sqrtn_error_hist_complex_dgp <- create_visualizer(
  .name = "sqrtn_error_hist_complex_dgp",
  .viz_fun = visualize_sqrtn_error_histograms,
  dgp_to_plot = "Complex (nonlinear) DGP",
  .doc_options = list(width = 12, height = 8)
)
