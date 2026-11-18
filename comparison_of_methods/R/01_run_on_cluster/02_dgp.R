# R/02_dgp.R

# Single-timepoint binary-treatment, binary-outcome DGP.
# Interaction-heavy so that Super Learner is useful.
# Target parameter: ATE = E[Y^1 - Y^0].

.truth_mc_n <- 2e6

expit <- plogis

f_g <- function(L1, L2, L3, L4, L5) {
  -0.4 +
    0.6 * L2 -
    0.5 * L3 +
    0.5 * L4 * L5 -
    0.4 * L1 * L2
}

f_q <- function(a, L1, L2, L3, L4, L5) {
  -0.8 +
    0.9 * a +
    0.6 * L2  +
    0.8 * L4 * L5 +
    0.7 * a * L1  -
    0.6 * a * L2
}

g_fun <- function(L1, L2, L3, L4, L5) {
  expit(f_g(L1, L2, L3, L4, L5))
}

q_fun <- function(a, L1, L2, L3, L4, L5) {
  expit(f_q(a, L1, L2, L3, L4, L5))
}

.compute_truth <- function(n = .truth_mc_n) {
  L1 <- rnorm(n)
  L2 <- rnorm(n)
  L3 <- rnorm(n)
  L4 <- rbinom(n, 1, 0.5)
  L5 <- runif(n, -1, 1)

  Q1 <- q_fun(1, L1, L2, L3, L4, L5)
  Q0 <- q_fun(0, L1, L2, L3, L4, L5)

  list(
    mu1 = mean(Q1),
    mu0 = mean(Q0),
    ate = mean(Q1 - Q0)
  )
}

truth_vec <- .compute_truth()

simulate_nonlinear_binary_data <- function(n = 1000) {
  L1 <- rnorm(n)
  L2 <- rnorm(n)
  L3 <- rnorm(n)
  L4 <- rbinom(n, 1, 0.5)
  L5 <- runif(n, -1, 1)

  g <- g_fun(L1, L2, L3, L4, L5)
  A <- rbinom(n, size = 1, prob = g)

  Q <- q_fun(A, L1, L2, L3, L4, L5)
  Y <- rbinom(n, size = 1, prob = Q)

  df <- data.frame(
    L1 = L1,
    L2 = L2,
    L3 = L3,
    L4 = L4,
    L5 = L5,
    A = A,
    Y = Y
  )

  list(df = df, truth = truth_vec)
}

nonlinear_binary_dgp <- simChef::create_dgp(
  .dgp_fun = function(n = sample_size) simulate_nonlinear_binary_data(n = n),
  .name = "Nonlinear binary-treatment DGP for Super Learner comparison",
  n = sample_size
)





# efficiency bound from DGP --------------------------------------------------------

# ------------------------------------------------------------------
# Monte Carlo estimate of the semiparametric efficiency bound
# for the binary-treatment ATE under the known DGP
# ------------------------------------------------------------------

estimate_ate_efficiency_bound <- function(n_mc = 2e6, seed = 123) {
  set.seed(seed)

  # Draw W from the true marginal law used in the DGP
  L1 <- rnorm(n_mc)
  L2 <- rnorm(n_mc)
  L3 <- rnorm(n_mc)
  L4 <- rbinom(n_mc, 1, 0.5)
  L5 <- runif(n_mc, -1, 1)

  # True nuisance functions from the DGP
  g  <- g_fun(L1, L2, L3, L4, L5)
  Q1 <- q_fun(1, L1, L2, L3, L4, L5)
  Q0 <- q_fun(0, L1, L2, L3, L4, L5)

  psi0 <- mean(Q1 - Q0)

  # Generate observed data O = (W, A, Y) from the true law
  A <- rbinom(n_mc, 1, g)
  QAW <- ifelse(A == 1, Q1, Q0)
  Y <- rbinom(n_mc, 1, QAW)

  # True EIF for the ATE
  eif <- (Q1 - Q0) +
    A / g * (Y - Q1) -
    (1 - A) / (1 - g) * (Y - Q0) -
    psi0

  list(
    psi0 = psi0,
    var_eif = mean(eif^2),   # mean(eif^2) since E{EIF}=0
    sd_eif = sqrt(mean(eif^2)),
    mc_mean_eif = mean(eif)
  )
}

eff_bound <- estimate_ate_efficiency_bound(n_mc = 2e6, seed = 1)

# eff_bound
#   $psi0
#   $var_eif      <-- asymptotic variance constant
#   $sd_eif
#   $mc_mean_eif  <-- should be near 0


