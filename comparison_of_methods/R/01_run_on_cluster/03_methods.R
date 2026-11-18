# R/03_methods.R

# Learner libraries --------------------------------------------------------

.make_nadir_binary_learners <- function() {
  list(
    mean = nadir::lnr_mean,
    glm  = nadir::lnr_glm,
    mars = nadir::lnr_earth,
    glmnet = nadir::lnr_glmnet
  )
}

.make_sl3_binary_stack <- function(formula) {
  sl3::Stack$new(
    sl3::Lrnr_mean$new(),
    sl3::Lrnr_glm$new(formula = formula),
    sl3::Lrnr_earth$new(formula = formula),
    sl3::Lrnr_glmnet$new(formula = formula)
  )
}

.make_sl3_superlearner <- function(formula = formula) {
  sl3::Lrnr_sl$new(
    learners = .make_sl3_binary_stack(formula = formula),
    metalearner = sl3::Lrnr_nnls$new()
  )
}

.make_sl3_mean_learner <- function() {
  sl3::Lrnr_mean$new()
}

# RieszCML nuisance fits with nadir ---------------------------------------

.formula_Q <- Y ~ A + L2 + L4*L5 + A*L1 + A*L2
.formula_g <- A ~ L2 + L3 + L4*L5 + L1*L2

.g_bound <- 0.025   # matches tmle3's default propensity bound
.truncate_g <- function(g) pmin(pmax(g, .g_bound), 1 - .g_bound)

.fit_nadir_sl <- function(data, formula) {
  nadir::super_learner(
    data = data, formula = formula,
    learners = .make_nadir_binary_learners(),
    outcome_type = "binary", n_folds = 10,
    determine_super_learner_weights = nadir::determine_super_learner_weights_nnls
  )
}

.Q_preds <- function(m_fit, data) {
  d1 <- data; d1$A <- 1
  d0 <- data; d0$A <- 0
  list(m  = m_fit$oof_predict(data),
       m1 = m_fit$oof_predict(d1),
       m0 = m_fit$oof_predict(d0))
}

.make_nuis_ate_nadir_sl <- function() {
  function(data) {
    c(.Q_preds(.fit_nadir_sl(data, .formula_Q), data),
      list(g = .truncate_g(.fit_nadir_sl(data, .formula_g)$oof_predict(data))))
  }
}

.make_nuis_ate_nadir_q_misspecified <- function() {
  function(data) {
    m_const <- rep(mean(data$Y), nrow(data))
    list(m = m_const, m1 = m_const, m0 = m_const,
         g = .truncate_g(.fit_nadir_sl(data, .formula_g)$oof_predict(data)))
  }
}

.make_nuis_ate_nadir_g_misspecified <- function() {
  function(data) {
    c(.Q_preds(.fit_nadir_sl(data, .formula_Q), data),
      list(g = .truncate_g(rep(mean(data$A), nrow(data)))))
  }
}

.soft_try <- function(expr) {
  tryCatch(
    list(value = eval.parent(substitute(expr)), ok = TRUE, err = NA_character_),
    error = function(e) list(value = NULL, ok = FALSE, err = conditionMessage(e))
  )
}

.result_from_fit <- function(fit, truth, method_name) {
  list(
    estimate = fit$estimate,
    truth = truth$ate,
    se = fit$se,
    ci_low = fit$ci_low,
    ci_high = fit$ci_high,
    method = method_name,
    ic_sd = if (!is.null(fit$ic)) stats::sd(fit$ic - fit$estimate) else NA_real_,
    eps = if (!is.null(fit$eps)) {
      if (length(fit$eps) == 1L) as.numeric(fit$eps) else mean(as.numeric(fit$eps))
    } else {
      NA_real_
    }
  )
}

# -------------------------------------------------------------------------
# RieszCML nuisance fits with one nuisance intentionally misspecified
# -------------------------------------------------------------------------


# Riesz one-step and TMLE --------------------------------------------------

riesz_onestep_nadir_sl_fun <- function(df, truth) {
  time_a <- Sys.time()
  rc_ate <- RieszCML::riesz_curve_catalog$ate(
    nuis = .make_nuis_ate_nadir_sl()
  )

  res <- .soft_try({
    fit <- RieszCML::riesz_estimate(
      data = df,
      rc = rc_ate
    )

    .result_from_fit(fit, truth, "riesz_onestep_nadir_sl")
  })
  time_b <- Sys.time()
  runtime <- lubridate::time_length(lubridate::interval(time_a, time_b), 'seconds')

  if (!res$ok) {
    return(list(
      estimate = NA_real_,
      truth = truth$ate,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      method = "riesz_onestep_nadir_sl",
      ic_sd = NA_real_,
      eps = NA_real_,
      runtime = runtime,
      err = res$err
    ))
  }
  res$value$runtime <- runtime
  return(res$value)
}

riesz_tmle_nadir_sl_fun <- function(df, truth) {
  time_a <- Sys.time()
  rc_ate <- RieszCML::riesz_curve_catalog$ate(
    nuis = .make_nuis_ate_nadir_sl()
  )

  res <- .soft_try({
    fit <- RieszCML::riesz_tmle(
      data = df,
      rc = rc_ate,
      fluctuation_type = "logistic",
      fluctuation_method = "covariate",
      bounds = c(0, 1),
      outcome_col = "Y"
    )

    .result_from_fit(fit, truth, "riesz_tmle_nadir_sl")
  })
  time_b <- Sys.time()
  runtime <- lubridate::time_length(lubridate::interval(time_a, time_b), 'seconds')

  if (!res$ok) {
    return(list(
      estimate = NA_real_,
      truth = truth$ate,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      method = "riesz_tmle_nadir_sl",
      ic_sd = NA_real_,
      eps = NA_real_,
      runtime = runtime,
      err = res$err
    ))
  }
  res$value$runtime <- runtime
  return(res$value)
}


# -------------------------------------------------------------------------
# Riesz one-step / TMLE with one nuisance misspecified
# -------------------------------------------------------------------------

riesz_onestep_nadir_q_misspecified_fun <- function(df, truth) {
  time_a <- Sys.time()
  rc_ate <- RieszCML::riesz_curve_catalog$ate(
    nuis = .make_nuis_ate_nadir_q_misspecified()
  )

  res <- .soft_try({
    fit <- RieszCML::riesz_estimate(
      data = df,
      rc = rc_ate
    )

    .result_from_fit(fit, truth, "riesz_onestep_nadir_q_misspecified")
  })

  time_b <- Sys.time()
  runtime <- lubridate::time_length(
    lubridate::interval(time_a, time_b), "seconds"
  )

  if (!res$ok) {
    return(list(
      estimate = NA_real_,
      truth = truth$ate,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      method = "riesz_onestep_nadir_q_misspecified",
      ic_sd = NA_real_,
      eps = NA_real_,
      runtime = runtime,
      err = res$err
    ))
  }

  res$value$runtime <- runtime
  res$value
}

riesz_onestep_nadir_g_misspecified_fun <- function(df, truth) {
  time_a <- Sys.time()
  rc_ate <- RieszCML::riesz_curve_catalog$ate(
    nuis = .make_nuis_ate_nadir_g_misspecified()
  )

  res <- .soft_try({
    fit <- RieszCML::riesz_estimate(
      data = df,
      rc = rc_ate
    )

    .result_from_fit(fit, truth, "riesz_onestep_nadir_g_misspecified")
  })

  time_b <- Sys.time()
  runtime <- lubridate::time_length(
    lubridate::interval(time_a, time_b), "seconds"
  )

  if (!res$ok) {
    return(list(
      estimate = NA_real_,
      truth = truth$ate,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      method = "riesz_onestep_nadir_g_misspecified",
      ic_sd = NA_real_,
      eps = NA_real_,
      runtime = runtime,
      err = res$err
    ))
  }

  res$value$runtime <- runtime
  res$value
}

riesz_tmle_nadir_q_misspecified_fun <- function(df, truth) {
  time_a <- Sys.time()
  rc_ate <- RieszCML::riesz_curve_catalog$ate(
    nuis = .make_nuis_ate_nadir_q_misspecified()
  )

  res <- .soft_try({
    fit <- RieszCML::riesz_tmle(
      data = df,
      rc = rc_ate,
      fluctuation_type = "logistic",
      fluctuation_method = "covariate",
      bounds = c(0, 1),
      outcome_col = "Y"
    )

    .result_from_fit(fit, truth, "riesz_tmle_nadir_q_misspecified")
  })

  time_b <- Sys.time()
  runtime <- lubridate::time_length(
    lubridate::interval(time_a, time_b), "seconds"
  )

  if (!res$ok) {
    return(list(
      estimate = NA_real_,
      truth = truth$ate,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      method = "riesz_tmle_nadir_q_misspecified",
      ic_sd = NA_real_,
      eps = NA_real_,
      runtime = runtime,
      err = res$err
    ))
  }

  res$value$runtime <- runtime
  res$value
}

riesz_tmle_nadir_g_misspecified_fun <- function(df, truth) {
  time_a <- Sys.time()
  rc_ate <- RieszCML::riesz_curve_catalog$ate(
    nuis = .make_nuis_ate_nadir_g_misspecified()
  )

  res <- .soft_try({
    fit <- RieszCML::riesz_tmle(
      data = df,
      rc = rc_ate,
      fluctuation_type = "logistic",
      fluctuation_method = "covariate",
      bounds = c(0, 1),
      outcome_col = "Y"
    )

    .result_from_fit(fit, truth, "riesz_tmle_nadir_g_misspecified")
  })

  time_b <- Sys.time()
  runtime <- lubridate::time_length(
    lubridate::interval(time_a, time_b), "seconds"
  )

  if (!res$ok) {
    return(list(
      estimate = NA_real_,
      truth = truth$ate,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      method = "riesz_tmle_nadir_g_misspecified",
      ic_sd = NA_real_,
      eps = NA_real_,
      runtime = runtime,
      err = res$err
    ))
  }

  res$value$runtime <- runtime
  res$value
}

# -------------------------------------------------------------------------
# tmle3 helper with one nuisance misspecified
# -------------------------------------------------------------------------

.run_tmle3_ate_custom <- function(df, truth, learner_Y, learner_A, method_name) {
  time_a <- Sys.time()

  res <- .soft_try({
    df_tmle <- data.table::as.data.table(df)

    node_list <- list(
      W = c("L1", "L2", "L3", "L4", "L5"),
      A = "A",
      Y = "Y"
    )

    ate_spec <- tmle3::tmle_ATE(
      treatment_level = 1,
      control_level = 0
    )

    learner_list <- list(
      A = learner_A,
      Y = learner_Y
    )

    fit <- tmle3::tmle3(
      tmle_spec = ate_spec,
      data = df_tmle,
      node_list = node_list,
      learner_list = learner_list
    )

    summ <- fit$summary
    ate_row <- summ[summ$param_type == "ATE", ]
    if (nrow(ate_row) == 0) {
      ate_row <- summ[1, ]
    }

    list(
      estimate = as.numeric(ate_row$tmle_est),
      truth = truth$ate,
      se = as.numeric(ate_row$se),
      ci_low = as.numeric(ate_row$lower),
      ci_high = as.numeric(ate_row$upper),
      method = method_name,
      ic_sd = NA_real_,
      eps = NA_real_
    )
  })

  time_b <- Sys.time()
  runtime <- lubridate::time_length(
    lubridate::interval(time_a, time_b), "seconds"
  )

  if (!res$ok) {
    return(list(
      estimate = NA_real_,
      truth = truth$ate,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      method = method_name,
      ic_sd = NA_real_,
      eps = NA_real_,
      runtime = runtime,
      err = res$err
    ))
  }

  res$value$runtime <- runtime
  res$value
}

tmle3_q_misspecified_fun <- function(df, truth) {
  .run_tmle3_ate_custom(
    df = df,
    truth = truth,
    learner_Y = .make_sl3_mean_learner(),
    learner_A = .make_sl3_superlearner(formula = "~ L2 + L3 + L4*L5 + L1 * L2"),
    method_name = "tmle3_q_misspecified"
  )
}

tmle3_g_misspecified_fun <- function(df, truth) {
  .run_tmle3_ate_custom(
    df = df,
    truth = truth,
    learner_Y = .make_sl3_superlearner(formula = "~ A + L2 + L4*L5 + A*L1 + A*L2"),
    learner_A = .make_sl3_mean_learner(),
    method_name = "tmle3_g_misspecified"
  )
}

# tmle3 with sl3 Super Learner --------------------------------------------

.run_tmle3_ate_sl <- function(df, truth, method_name) {
  time_a <- Sys.time()
  res <- .soft_try({
    df_tmle <- data.table::as.data.table(df)

    node_list <- list(
      W = c("L1", "L2", "L3", "L4", "L5"),
      A = "A",
      Y = "Y"
    )

    ate_spec <- tmle3::tmle_ATE(
      treatment_level = 1,
      control_level = 0
    )

    learner_Y <- .make_sl3_superlearner(formula = "~ A + L2 + L4 * L5 + A * L1 + A*L2")
    learner_A <- .make_sl3_superlearner(formula = "~ L2 + L3 + L4*L5 + L1 * L2")

    learner_list <- list(
      A = learner_A,
      Y = learner_Y
    )

    folds <- origami::make_folds(
      n = nrow(df_tmle),
      fold_fun = origami::folds_vfold,
      V = 10
    )

    tmle_task <- ate_spec$make_tmle_task(
      data = df_tmle,
      node_list = node_list,
      folds = folds
    )

    # fit <- tmle3::tmle3(
    #   tmle_spec = ate_spec,
    #   data = df_tmle,
    #   node_list = node_list,
    #   learner_list = learner_list
    # )

    initial_likelihood <- ate_spec$make_initial_likelihood(tmle_task, learner_list)

    updater <- ate_spec$make_updater()
    targeted_likelihood <- ate_spec$make_targeted_likelihood(initial_likelihood, updater)

    tmle_params <- ate_spec$make_params(tmle_task, targeted_likelihood)
    updater$tmle_params <- tmle_params

    fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)



    # tmle3_Fit objects expose a summary table
    summ <- fit$summary

    # pull the ATE row; usually there is only one row here
    ate_row <- summ[summ$param_type == "ATE", ]
    if (nrow(ate_row) == 0) {
      ate_row <- summ[1, ]
    }

    # results
    estimate <- as.numeric(ate_row$tmle_est)
    se_hat   <- as.numeric(ate_row$se)
    ci_low   <- as.numeric(ate_row$lower)
    ci_high  <- as.numeric(ate_row$upper)

    list(
      estimate = estimate,
      truth = truth$ate,
      se = se_hat,
      ci_low = ci_low,
      ci_high = ci_high,
      method = method_name,
      ic_sd = sd(fit$estimates[[1]]$IC),
      eps = fit$updater$epsilons[[1]]$Y,
      steps = fit$steps
    )
  })
  time_b <- Sys.time()
  runtime <- lubridate::time_length(lubridate::interval(time_a, time_b))

  if (!res$ok) {
    return(list(
      estimate = NA_real_,
      truth = truth$ate,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      method = method_name,
      ic_sd = NA_real_,
      eps = NA_real_,
      runtime = runtime,
      err = res$err
    ))
  }

  res$value$runtime <- runtime
  return(res$value)
}

tmle3_sl_fun <- function(df, truth) {
  .run_tmle3_ate_sl(
    df = df,
    truth = truth,
    method_name = "tmle3_sl"
  )
}


# simChef methods ----------------------------------------------------------

m_riesz_onestep_nadir_sl <- simChef::create_method(
  .method_fun = riesz_onestep_nadir_sl_fun,
  .name = "riesz_onestep_nadir_sl"
)

m_riesz_tmle_nadir_sl <- simChef::create_method(
  .method_fun = riesz_tmle_nadir_sl_fun,
  .name = "riesz_tmle_nadir_sl"
)

m_tmle3_sl <- simChef::create_method(
  .method_fun = tmle3_sl_fun,
  .name = "tmle3_sl"
)

# -------------------------------------------------------------------------
# simChef methods for DR checks
# -------------------------------------------------------------------------

m_riesz_onestep_nadir_q_misspecified <- simChef::create_method(
  .method_fun = riesz_onestep_nadir_q_misspecified_fun,
  .name = "riesz_onestep_nadir_q_misspecified"
)

m_riesz_onestep_nadir_g_misspecified <- simChef::create_method(
  .method_fun = riesz_onestep_nadir_g_misspecified_fun,
  .name = "riesz_onestep_nadir_g_misspecified"
)

m_riesz_tmle_nadir_q_misspecified <- simChef::create_method(
  .method_fun = riesz_tmle_nadir_q_misspecified_fun,
  .name = "riesz_tmle_nadir_q_misspecified"
)

m_riesz_tmle_nadir_g_misspecified <- simChef::create_method(
  .method_fun = riesz_tmle_nadir_g_misspecified_fun,
  .name = "riesz_tmle_nadir_g_misspecified"
)

m_tmle3_q_misspecified <- simChef::create_method(
  .method_fun = tmle3_q_misspecified_fun,
  .name = "tmle3_q_misspecified"
)

m_tmle3_g_misspecified <- simChef::create_method(
  .method_fun = tmle3_g_misspecified_fun,
  .name = "tmle3_g_misspecified"
)
