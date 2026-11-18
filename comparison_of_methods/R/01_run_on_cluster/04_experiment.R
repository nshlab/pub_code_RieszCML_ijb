# R/04_experiment.R

experiment <-
  simChef::create_experiment(
    name = "RieszCML nadir::super_learner vs tmle3 sl3 Super Learner"
  ) |>
  simChef::add_dgp(dgp = nonlinear_binary_dgp) |>
  simChef::add_method(m_riesz_onestep_nadir_sl) |>
  simChef::add_method(m_riesz_tmle_nadir_sl) |>
  simChef::add_method(m_tmle3_sl) |>
  simChef::add_method(m_riesz_onestep_nadir_q_misspecified) |>
  simChef::add_method(m_riesz_onestep_nadir_g_misspecified) |>
  simChef::add_method(m_riesz_tmle_nadir_q_misspecified) |>
  simChef::add_method(m_riesz_tmle_nadir_g_misspecified) |>
  simChef::add_method(m_tmle3_q_misspecified) |>
  simChef::add_method(m_tmle3_g_misspecified)
