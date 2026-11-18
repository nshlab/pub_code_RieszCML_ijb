## =============================================================================
## Stochastic-interventional vaccine efficacy, SVE(delta), under two-phase
## sampling of the immune marker and right-censoring of the outcome,
## estimated with {RieszCML} + cross-fitted super-learner nuisances ({nadir}).
##
## COVE-calibrated demonstration analysis.
##
## Target parameter (Hejazi et al. 2023, IJID; Gilbert et al. 2024, Vaccine):
##
##   SVE(delta) = 1 - psi1(delta) / psi0 ,
##
##   psi1(delta) = E_{S,L | A=1}[ Qbar(S + delta, L) ],
##                 Qbar(s,l) = E[ Y | S=s, L=l, A=1, C=1 ]
##   psi0        = E_L[ E( Y | A=0, L, C=1 ) ]
##
## Observed data:  O = (L, A, RS, C, CY, R)  --  see sim_data() below.
##
## What is new relative to the previous version of this script:
##
##   1. NUISANCES VIA nadir::crossfit_super_learner() WITH OUT-OF-FOLD
##      PREDICTIONS. Every nuisance regression is now a cross-fitted super
##      learner; the influence-curve machinery consumes only out-of-fold
##      (OOF) predictions via `$oof_predictions` and
##      `$oof_predict_modified(modify)` (e.g., predictions at the shifted
##      marker S + delta), the property required of cross-fitted nuisance
##      estimators in semiparametric workflows.
##
##   2. SUBSET FITS VIA ZERO WEIGHTS. RieszCML requires length-n nuisance
##      vectors, but several nuisances are fit on subsets (phase-two sample,
##      uncensored units) with Horvitz-Thompson reweighting. We pass the full
##      data with observation weights that are ZERO off the relevant subset
##      (e.g. R*C/pw for the outcome regression). nadir propagates these
##      weights to the base learners AND to the ensemble-weight determination,
##      so zero-weight rows influence nothing -- yet every row still receives
##      an OOF prediction. This replaces the old subset-fit-then-predict
##      plumbing entirely.
##
##   3. NUISANCES FIT ONCE, SHARED EVERYWHERE. The delta-independent
##      cross-fits (outcome regression, censoring mechanism, marker-density
##      working model, denominator models) are fit a single time; each delta
##      on the grid only requires cheap `$oof_predict_modified()` evaluations
##      plus one small stage-2 cross-fit. The one-step estimator and the TMLE
##      are computed from the SAME precomputed RieszCurve objects (RieszCurve
##      accepts a precomputed named list for `nuis`), so adding an estimator
##      is free.
##
##   4. CASE-STRATIFIED CROSS-FITTING FOLDS. With ~46 vaccine breakthrough
##      cases among ~14,000 vaccinees (COVE calibration), random folds can
##      starve a fold of cases. cv_strat_CY_schema() assigns folds randomly
##      WITHIN case strata; because it reads the stratum off the data, the
##      same schema serves the outer cross-fitting split and the inner
##      ensemble-weight split.
##
##   5. LEARNER LIBRARY: mean + GLM + GAM + glmnet. Deliberately light and
##      fast. NOTE: MARS ({earth}) was removed from the library because case
##      weights disable earth's fast update algorithm, inflating fit times by
##      roughly two orders of magnitude (~0.02s -> ~7.5s per fit at n = 3000
##      in our profiling) -- and every fit in this analysis is HT-weighted.
##      mgcv::gam() with per-learner smooth formulas (nadir's formula-list
##      interface) provides the flexible learner at no weighted-fit penalty.
##
##   6. POSITIVITY / STABILITY NOTE. In the rare-outcome COVE calibration,
##      downward shifts beyond delta ~ -1 SD push the marker density ratios
##      g(S - delta | L)/g(S | L) into a near-positivity-violation regime and
##      the TMLE fluctuation can destabilize (we observed psi1 collapse at
##      delta = -1.5). The default grid below, delta in [-1, 1.5], stays on
##      well-supported territory; widen it deliberately, not by default.
##
## Estimation strategy (paper: Theorem 2 / Algorithm 1): the numerator is a
## nested linear functional estimated with a ComposedRieszCurve built from
## THREE stages, each a primitive whose Riesz representer is already known
## from the literature. No representer for the combined problem (an MTP with
## right censoring under two-phase sampling) is derived anywhere -- to our
## knowledge none appears in the literature -- because none is needed: the
## composition assembles it. The nesting is
##
##   psi1 = E[ Q_1(V) ],       Q_1(V)     = E[ h_2(S, L) | R = 1, V ],
##                             h_2(S, L)  = Q_2(S + delta, L),
##                             Q_2(S, L)  = E[ h_1 | S, L ],
##                             h_1(S, L)  = Qbar(S, L),
##                             Qbar(S, L) = E[ CY | S, L, C = 1 ],
##
## with V = (CY, C, L) the always-observed variables:
##   * stage 1 (innermost, acts on CY): the censoring/missingness primitive
##     (cf. `riesz_curve_catalog$missing_mean`; Rose & van der Laan 2011),
##       alpha_1 = C / pC(S, L),          alpha*_1 = 1 / pC(S, L).
##   * stage 2: the modified-treatment-policy primitive
##     (cf. `riesz_curve_catalog$mtp`; Diaz et al. 2021),
##       alpha_2  = g(S - delta | L) / g(S | L),
##       alpha*_2 = g(S | L) / g(S + delta | L).
##     Because h_1 is a function of (S, L), its "regression" Q_2 =
##     E[h_1 | S, L] conditions on every argument of h_1 and is available in
##     CLOSED FORM: Q_2 = h_1. We therefore declare f as the `h[j=1]` column
##     itself, so this stage's residual term in the influence curve is
##     exactly zero -- yet the stage still contributes its representer to
##     every suffix product, which is its entire job.
##   * stage 3 (outermost): the two-phase sampling primitive (Rose & van der
##     Laan 2011; Hejazi et al. 2020),
##       alpha_3 = R / pi(V),             alpha*_3 = 1 / pi(V).
##   ComposedRieszCurve assembles the Theorem-2 EIF with suffix products
##       omega_1 = alpha_1 * alpha_2 * alpha_3,
##       omega_2 = alpha_2 * alpha_3,     omega_3 = alpha_3,
##   so the familiar combined weight ratio * (C/pC) * (R/pi) multiplying the
##   innermost residual (CY - Qbar) EMERGES from the composition rather than
##   being derived by hand, and the Horvitz-Thompson reweighting of the
##   phase-two fluctuations rides inside the omegas automatically (no
##   `fluctuation_weights` anywhere).
##   * denominator: single RieszCurve with the product representer
##       [1(A=0)/P(A=0)] * [C/pC(A=0, L)].
##   * SVE(delta) = 1 - psi1/psi0 via riesz_delta() (delta method).
##
## Two-phase weights pw are the DESIGN sampling probabilities (known in a
## designed two-phase sample). Estimating pi flexibly instead requires the
## additional targeting step of Hejazi et al. (2020) (clever covariate
## G/g_C), which is not implemented here.
##
## Runtime and memory guidance (measured at the defaults below):
##   * primary analysis (n = 28,000, 5 x 5 folds, 11 deltas):
##       ~15-20 min single-core, peak RSS ~1.5 GB SEQUENTIAL. Cross-fitting
##       parallelizes over the 5 outer folds via {future}; the run block
##       below picks this up from SVE_WORKERS or SLURM_CPUS_PER_TASK (workers
##       beyond the fold count buy nothing). Under multisession each worker
##       is a separate R session peaking at ~1-1.5 GB, so budget memory as
##       ~2 GB (parent) + workers x 1.5 GB: request ~16 GB for 5 workers, or
##       set SVE_WORKERS to fewer than the CPU count on a tight allocation.
##       A worker OOM-killed by the kernel surfaces as a {future}
##       "FutureInterruptError ... MultisessionFuture interrupted" in the
##       parent -- the job need not be flagged OUT_OF_MEMORY by SLURM.
##   * one Monte Carlo replicate (n = 28,000, 3 x 3 folds, 3 deltas):
##       ~2 min single-core, peak RSS ~0.75 GB. On a cluster, run the Monte
##       Carlo as a SLURM JOB ARRAY (one replicate per task; see
##       `sve_cove_mc_rep.R` / `sbatch_sve_cove_mc_array.sh` and collect with
##       `sve_cove_mc_collect.R`) rather than via in-session {future}
##       workers: memory then never exceeds the single-replicate footprint,
##       tasks are requeue-safe, and the scheduler backfills small tasks
##       quickly. If you do use `run_simulation()` with
##       plan(multisession, workers = K), budget roughly K x 1 GB of memory
##       and M / K x 2 min of walltime, and pin BLAS threading
##       (OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1) to avoid
##       oversubscription. Note: {future}'s "total size of the globals ...
##       exceeds the maximum allowed size" error is an option limit
##       (future.globals.maxSize), not an out-of-memory condition.
## =============================================================================

invisible(suppressMessages({
  library(RieszCML)
  library(nadir)
  library(ggplot2)
  ## Register S3 predict methods needed by learner closures (required in
  ## fresh sessions if cross-fit objects are restored from .rds files).
  loadNamespace("glmnet"); loadNamespace("mgcv")
}))

options(future.globals.maxSize = +Inf)
mem.maxVSize(vsize = +Inf)

## =============================================================================
## 1. Data-generating process
## =============================================================================

sim_data <- function(n, scenario = c("cove", "base")) {
  scenario <- match.arg(scenario)
  L1 <- rnorm(n)
  L2 <- rbinom(n, 1, 0.5)
  A  <- rbinom(n, 1, 0.5)

  ## post-vaccination immune marker (vaccine arm only), Gaussian given L, on
  ## a standardized (sd ~ 1) scale so delta is in SD units, as in the
  ## HVTN 505 / COVE immune-correlates analyses
  S_full <- 0.5 + 0.3 * L1 + 0.2 * L2 + rnorm(n, sd = 1)
  S_full[A == 0] <- NA

  if (scenario == "base") {
    pC <- plogis(2 + 0.2 * L1 - 0.3 * L2 +
                   ifelse(A == 1, 0.3 * (S_full - 0.5), 0))
    C  <- rbinom(n, 1, pC)
    pY <- ifelse(A == 1,
                 plogis(-2.0 - 0.8 * S_full + 0.4 * L1 + 0.3 * L2),
                 plogis(-1.2 + 0.4 * L1 + 0.3 * L2))
    Y  <- rbinom(n, 1, pY)
    pi0 <- plogis(-0.8 + 0.5 * L1)
  } else {
    ## COVE-calibrated (Gilbert et al. 2022, Science): rare outcome
    ## (~0.004 vaccine, ~0.061 placebo risk through 100 days post day 57),
    ## light censoring, case-cohort sampling: all observed vaccine
    ## breakthrough cases + an ~7% immunogenicity subcohort (1010/14064)
    pC <- plogis(3 + 0.2 * L1 - 0.1 * L2)
    C  <- rbinom(n, 1, pC)
    pY <- ifelse(A == 1,
                 plogis(-5.26 - 0.9 * S_full + 0.3 * L1 + 0.2 * L2),
                 plogis(-2.80 + 0.3 * L1 + 0.2 * L2))
    Y  <- rbinom(n, 1, pY)
    pi0 <- plogis(-2.5 + 0.3 * L1)
  }
  CY <- C * Y

  pi_R <- ifelse(C == 1 & Y == 1, 1, pi0)
  R    <- ifelse(A == 1, rbinom(n, 1, pi_R), 0L)

  data.frame(L1 = L1, L2 = L2, A = A,
             S  = ifelse(R == 1, S_full, NA),
             C  = C, CY = CY, R = R,
             pw = ifelse(A == 1, pi_R, NA))   # design sampling probability
}

## Monte Carlo truth -----------------------------------------------------------
get_truth <- function(delta_grid, scenario = c("cove", "base"),
                      n_mc = 2e6, seed = 777) {
  scenario <- match.arg(scenario)
  set.seed(seed)
  L1 <- rnorm(n_mc); L2 <- rbinom(n_mc, 1, 0.5)
  S  <- 0.5 + 0.3 * L1 + 0.2 * L2 + rnorm(n_mc, sd = 1)
  if (scenario == "base") {
    psi0 <- mean(plogis(-1.2 + 0.4 * L1 + 0.3 * L2))
    psi1 <- vapply(delta_grid, function(d)
      mean(plogis(-2.0 - 0.8 * (S + d) + 0.4 * L1 + 0.3 * L2)), numeric(1))
  } else {
    psi0 <- mean(plogis(-2.80 + 0.3 * L1 + 0.2 * L2))
    psi1 <- vapply(delta_grid, function(d)
      mean(plogis(-5.26 - 0.9 * (S + d) + 0.3 * L1 + 0.2 * L2)), numeric(1))
  }
  data.frame(delta = delta_grid, psi1 = psi1, psi0 = psi0,
             sve = 1 - psi1 / psi0)
}

## =============================================================================
## 2. Cross-fitted nuisance estimation ({nadir})
## =============================================================================

## ---- fold schema stratified on case status ----------------------------------
## Random fold assignment within CY strata, so the rare cases are spread
## evenly across folds. Reads the stratum off the data itself, so the same
## schema serves the outer cross-fitting split and the inner ensemble-weight
## split (which sees only fold-subsetted data).
cv_strat_CY_schema <- function(data, n_folds) {
  fold <- integer(nrow(data))
  for (s in unique(data$CY)) {
    idx <- which(data$CY == s)
    fold[idx] <- sample(rep_len(seq_len(n_folds), length(idx)))
  }
  list(
    training_data   = lapply(seq_len(n_folds),
                             function(k) data[fold != k, , drop = FALSE]),
    validation_data = lapply(seq_len(n_folds),
                             function(k) data[fold == k, , drop = FALSE])
  )
}

## ---- light-and-fast learner library -----------------------------------------
## mean + GLM + GAM + glmnet: all four take Horvitz-Thompson observation
## weights without a speed penalty. MARS/{earth} is deliberately avoided:
## case weights disable earth's fast update algorithm, inflating fit times by
## ~two orders of magnitude, and every fit here is weighted.
sve_learners <- list(mean   = nadir::lnr_mean,
                     glm    = nadir::lnr_glm,
                     gam    = nadir::lnr_gam,
                     glmnet = nadir::lnr_glmnet)

## ---- delta-independent cross-fits (vaccine arm) ------------------------------
## All models are fit on the FULL vaccinee data with HT weights that are zero
## off the relevant subset, so every row receives an out-of-fold prediction
## while only the intended subset (with its 1/pw reweighting) influences the
## base learners and the ensemble weights.
fit_vaccinee_nuisances <- function(df_v, n_folds = 5, inner_n_folds = 5,
                                   learners = sve_learners) {
  w_QC <- df_v$R * df_v$C / df_v$pw   # phase-two & uncensored, HT-weighted
  w_R  <- df_v$R / df_v$pw            # phase-two, HT-weighted

  ## outcome regression Qbar(S, L) = E[CY | S, L, C = 1]
  cf_Q <- crossfit_super_learner(
    data = df_v,
    formulas = list(.default = CY ~ S + L1 + L2,
                    gam      = CY ~ s(S) + s(L1) + L2),
    learners = learners, outcome_type = "binary",
    weights = w_QC, n_folds = n_folds, inner_n_folds = inner_n_folds,
    cv_schema = cv_strat_CY_schema, inner_cv_schema = cv_strat_CY_schema)

  ## censoring mechanism pC(S, L) = P(C = 1 | S, L)
  cf_C <- crossfit_super_learner(
    data = df_v,
    formulas = list(.default = C ~ S + L1 + L2,
                    gam      = C ~ s(S) + s(L1) + L2),
    learners = learners, outcome_type = "binary",
    weights = w_R, n_folds = n_folds, inner_n_folds = inner_n_folds,
    cv_schema = cv_strat_CY_schema, inner_cv_schema = cv_strat_CY_schema)

  ## Gaussian working model for the marker density g(s | L): cross-fitted
  ## conditional mean, residual scale from HT-weighted out-of-fold residuals
  cf_S <- crossfit_super_learner(
    data = df_v,
    formulas = list(.default = S ~ L1 + L2,
                    gam      = S ~ s(L1) + L2),
    learners = learners, outcome_type = "continuous",
    weights = w_R, n_folds = n_folds, inner_n_folds = inner_n_folds,
    cv_schema = cv_strat_CY_schema, inner_cv_schema = cv_strat_CY_schema)
  mu  <- cf_S$oof_predictions
  sig <- sqrt(weighted.mean((df_v$S - mu)^2, w_R))

  list(cf_Q = cf_Q, cf_C = cf_C, mu = mu, sig = sig,
       n_folds = n_folds, inner_n_folds = inner_n_folds, learners = learners)
}

## ---- per-delta nuisance evaluation (out-of-fold, no refitting) ---------------
shift_S <- function(delta) function(d) { d$S <- d$S + delta; d }

## stages 1 (censoring) and 2 (MTP): out-of-fold evaluations of the
## delta-independent cross-fits. Stage 2's regression Q_2 = E[h_1 | S, L] is
## closed-form (see the header), so its only nuisances are the shifted
## outcome-regression evaluation and the density ratio at observed/shifted
## exposure.
build_censoring_nuis <- function(nf) {
  list(m  = nf$cf_Q$oof_predictions,     # Qbar at observed (S, L)
       pC = nf$cf_C$oof_predictions)     # P(C = 1 | S, L)
}

build_mtp_nuis <- function(nf, df_v, delta) {
  list(
    m_shift    = nf$cf_Q$oof_predict_modified(shift_S(delta)),
    ## alpha at observed data: g(S - delta | L) / g(S | L)
    ratio      = dnorm(df_v$S - delta, nf$mu, nf$sig) /
                 dnorm(df_v$S,         nf$mu, nf$sig),
    ## alpha under the intervention S -> S + delta
    ratio_star = dnorm(df_v$S,         nf$mu, nf$sig) /
                 dnorm(df_v$S + delta, nf$mu, nf$sig)
  )
}

## stage 3: G(V) = E[h_2 | R = 1, V], V = (CY, C, L). Indicator weights mean
## only sampled rows train the model (fit "among R = 1"), yet every row --
## sampled or not -- receives an out-of-fold prediction, exactly what the
## plug-in E[G(V)] requires. h2 here is the MTP stage's counterfactual
## prediction (the same vector ComposedRieszCurve exposes to stage 3 as
## `h[j=2]`).
build_sampling_nuis <- function(nf, df_v, h2) {
  d_G <- df_v
  d_G$h2 <- h2
  cf_G <- crossfit_super_learner(
    data = d_G,
    formulas = list(.default = h2 ~ CY + C + L1 + L2,
                    gam      = h2 ~ CY + C + s(L1) + L2),
    learners = nf$learners, outcome_type = "continuous",
    weights = as.numeric(df_v$R == 1),
    n_folds = nf$n_folds, inner_n_folds = nf$inner_n_folds,
    cv_schema = cv_strat_CY_schema, inner_cv_schema = cv_strat_CY_schema)
  list(G = cf_G$oof_predictions)
}

## ---- denominator nuisances (full sample) -------------------------------------
fit_denominator_nuisances <- function(df, n_folds = 5, inner_n_folds = 5,
                                      learners = sve_learners) {
  d <- df[, c("CY", "C", "A", "L1", "L2")]
  cf_m <- crossfit_super_learner(
    data = d,
    formulas = list(.default = CY ~ A * (L1 + L2),
                    gam      = CY ~ A * (L1 + L2) + s(L1)),
    learners = learners, outcome_type = "binary",
    weights = as.numeric(df$C == 1),
    n_folds = n_folds, inner_n_folds = inner_n_folds,
    cv_schema = cv_strat_CY_schema, inner_cv_schema = cv_strat_CY_schema)
  cf_pC <- crossfit_super_learner(
    data = d,
    formulas = list(.default = C ~ A * (L1 + L2),
                    gam      = C ~ A * (L1 + L2) + s(L1)),
    learners = learners, outcome_type = "binary",
    n_folds = n_folds, inner_n_folds = inner_n_folds,
    cv_schema = cv_strat_CY_schema, inner_cv_schema = cv_strat_CY_schema)
  set_A0 <- function(d) { d$A <- 0; d }
  list(m0  = cf_m$oof_predict_modified(set_A0),
       pC0 = cf_pC$oof_predict_modified(set_A0),
       gA  = mean(df$A))
}

## =============================================================================
## 3. Riesz curves (precomputed nuisance lists shared by one-step and TMLE)
## =============================================================================

## numerator: composed THREE-stage curve, one literature primitive per stage.
## Stage order innermost-first; the composed EIF (Theorem 2) is
##   G(V) + (R/pi) * (h_2 - G)
##        + (R/pi) * ratio * (h_1 - Q_2)                 [ = 0: Q_2 = h_1 ]
##        + (R/pi) * ratio * (C/pC) * (CY - Qbar),
## identical to the EIF one would derive for the combined problem -- the
## middle residual vanishes because Q_2 is closed-form -- but obtained here
## without deriving any combined representer.
make_rc_numerator <- function(nuis_cens, nuis_mtp, nuis_samp) {
  rc_cens <- RieszCurve$new(                    # 1. censoring primitive
    nuis       = nuis_cens,
    alpha      = ~ C / pC,
    alpha_star = ~ 1 / pC,
    f          = ~ m,
    h          = ~ m,                           # evaluation under C := 1
    ic_expr    = ~ h + alpha * (CY - f)
  )
  rc_mtp <- RieszCurve$new(                     # 2. MTP primitive
    nuis       = nuis_mtp,
    alpha      = ~ ratio,
    alpha_star = ~ ratio_star,
    f          = ~ `h[j=1]`,                    # Q_2 = h_1, in closed form
    h          = ~ m_shift,                     # Q_2 evaluated at S + delta
    ic_expr    = ~ h + alpha * (`h[j=1]` - f)   # residual identically zero
  )
  rc_samp <- RieszCurve$new(                    # 3. two-phase primitive
    nuis       = nuis_samp,
    alpha      = ~ R / pw,                      # pw read off the data
    alpha_star = ~ 1 / pw,                      # all units get the update
    f          = ~ G,
    h          = ~ G,
    ic_expr    = ~ h + alpha * (`h[j=2]` - f)
  )
  ComposedRieszCurve$new(rc_list = list(rc_cens, rc_mtp, rc_samp),
                         outcome_col = "CY")
}

## denominator: placebo-arm counterfactual risk with IPCW; a single curve
## with a product representer within the stage.
make_rc_denominator <- function(nuis_den) {
  RieszCurve$new(
    nuis       = nuis_den,
    alpha      = ~ (1 - A) / (1 - gA) * C / pC0,
    alpha_star = ~ 1 / ((1 - gA) * pC0),        # at (A = 0, C = 1)
    f          = ~ m0,
    h          = ~ m0,
    ic_expr    = ~ h + alpha * (CY - f),
    targeting_steps = list(m0 = list(alpha_star = ~ 1 / ((1 - gA) * pC0)))
  )
}

## ---- expand a vaccinee-subset RieszFit to a full-sample RieszFit -------------
## psi1 is a subgroup (A = 1) parameter; its full-sample EIF is
## (A / P(A=1)) * (phi_v - psi1), so riesz_delta() can combine it with the
## full-sample denominator fit.
pad_fit_to_full <- function(fit_v, A_full, idx_v, significance_alpha = 0.05) {
  n  <- length(A_full)
  pA <- mean(A_full)
  ic_v_centered <- fit_v$ic - mean(fit_v$ic)
  ic_c <- numeric(n)
  ic_c[idx_v] <- ic_v_centered / pA
  v  <- stats::var(ic_c) / n
  se <- sqrt(v)
  RieszFit$new(
    estimate = fit_v$estimate, var = v, se = se,
    ci_low  = fit_v$estimate + qnorm(significance_alpha / 2) * se,
    ci_high = fit_v$estimate + qnorm(1 - significance_alpha / 2) * se,
    significance_alpha = significance_alpha,
    estimator = fit_v$estimator, n = n,
    ic = fit_v$estimate + ic_c
  )
}

## ---- lightweight phase logging ------------------------------------------------
## Prints wall time and this process's current/peak resident memory (VmRSS /
## VmHWM from /proc, Linux) at each phase, so an out-of-memory failure can be
## localized from the .out log. Note: this sees the PARENT process only --
## multisession workers are separate processes -- so it is fully informative
## under a sequential run (SVE_WORKERS=1), which is the recommended way to
## run the one-off primary analysis.
mem_note <- function(label) {
  st <- tryCatch(readLines("/proc/self/status"), error = function(e) character())
  pick <- function(k) {
    v <- grep(paste0("^", k, ":"), st, value = TRUE)
    if (length(v)) gsub("\\s+", " ", sub(paste0("^", k, ":\\s+"), "", v[1])) else "NA"
  }
  cat(sprintf("[mem %s] %-32s VmRSS=%s  VmHWM=%s\n",
              format(Sys.time(), "%H:%M:%S"), label,
              pick("VmRSS"), pick("VmHWM")))
}

## =============================================================================
## 4. Full analysis over a delta grid
## =============================================================================
## Fits all delta-independent nuisances ONCE; per delta only cheap out-of-fold
## re-evaluations (and the small stage-2 cross-fit) are computed, and the
## one-step and TMLE share the same fitted curves.
sve_analysis <- function(df, delta_grid,
                         estimators = c("onestep", "tmle"),
                         n_folds = 5, inner_n_folds = 5,
                         learners = sve_learners) {
  df_v  <- df[df$A == 1, ]
  idx_v <- which(df$A == 1)

  ## Fill the unmeasured marker with an arbitrary finite value so that
  ## length-n stage-1 predictions exist for every row; these rows are inert
  ## because every term in which they appear -- the composed influence curve,
  ## the stage-1 fluctuation score (omega_1), and the stage-2 regression (fit
  ## with weights R) -- carries the factor R = 0.
  df_v$S[df_v$R == 0] <- 0

  mem_note("sve_analysis: start")
  nf_v   <- fit_vaccinee_nuisances(df_v, n_folds, inner_n_folds, learners)
  mem_note("vaccinee nuisances fit")
  nu_den <- fit_denominator_nuisances(df, n_folds, inner_n_folds, learners)
  mem_note("denominator nuisances fit")
  rc_den <- make_rc_denominator(nu_den)

  fit_den <- list(
    onestep = riesz_estimate(df, rc_den),
    tmle    = riesz_tmle(df, rc_den, fluctuation_type = "logistic",
                         fluctuation_method = "weighted",
                         bounds = c(0, 1), outcome_col = "CY"))

  out <- list()
  for (delta in delta_grid) {
    nuis_cens <- build_censoring_nuis(nf_v)
    nuis_mtp  <- build_mtp_nuis(nf_v, df_v, delta)
    nuis_samp <- build_sampling_nuis(nf_v, df_v, h2 = nuis_mtp$m_shift)
    rc_num <- make_rc_numerator(nuis_cens, nuis_mtp, nuis_samp)

    for (est in estimators) {
      fit_num <- if (est == "onestep") {
        riesz_estimate(df_v, rc_num)
      } else {
        ## Sequential (Algorithm 1) TMLE on the composed curve; no
        ## `fluctuation_weights` needed -- the HT reweighting rides inside
        ## the suffix product omega_1 = ratio * (C/pC) * (R/pw).
        riesz_tmle(df_v, rc_num, fluctuation_type = "logistic",
                   fluctuation_method = "weighted",
                   bounds = c(0, 1), outcome_col = "CY")
      }
      fit_num_full <- pad_fit_to_full(fit_num, df$A, idx_v)
      fit_sve <- riesz_delta(
        fit_num_full, fit_den[[est]],
        g    = function(t) 1 - t[1] / t[2],
        grad = function(t) c(-1 / t[2], t[1] / t[2]^2),
        parameter = sprintf("SVE(delta = %.2f)", delta),
        estimator = est)
      mem_note(sprintf("delta % .2f (%s) done", delta, est))
      out[[length(out) + 1L]] <- data.frame(
        delta = delta, estimator = est,
        sve = fit_sve$estimate, se = fit_sve$se,
        lo = fit_sve$ci_low, hi = fit_sve$ci_high,
        psi1 = fit_num$estimate, psi1_se = fit_num_full$se,
        psi0 = fit_den[[est]]$estimate, psi0_se = fit_den[[est]]$se)
    }
  }
  out <- do.call(rbind, out)
  ## marker-support diagnostics over a fine grid spanning delta_grid, for the
  ## figure's support strip (see sve_support_data / sve_figure)
  attr(out, "support") <- sve_support_data(
    nf_v, df_v, seq(min(delta_grid), max(delta_grid), by = 0.05))
  out
}

## =============================================================================
## 5. Presentation: figure and table
## =============================================================================

## ---- support diagnostics ------------------------------------------------------
## Where does the shifted counterfactual marker S + delta leave the observed
## data? For each delta, `extrap_lo` is the share of counterfactual evaluation
## points S + delta falling BELOW the q-th percentile of the observed
## phase-two marker distribution (computed with the design sampling weights,
## so the share refers to the full vaccine arm); the "rug" attribute holds,
## per sampled subject, the shift q_lo - S_i at which that subject's
## evaluation point exits the observed support.
sve_support_data <- function(nf, df_v, delta_grid, q = 0.01) {
  obs  <- df_v$R == 1
  w    <- (1 / df_v$pw)[obs]
  S2   <- df_v$S[obs]
  q_lo <- quantile(S2, q, names = FALSE)
  q_hi <- quantile(S2, 1 - q, names = FALSE)
  out <- data.frame(
    delta = delta_grid,
    extrap_lo = sapply(delta_grid, function(d) weighted.mean(S2 + d < q_lo, w)),
    extrap_hi = sapply(delta_grid, function(d) weighted.mean(S2 + d > q_hi, w)),
    ess = sapply(delta_grid, function(d) {
      ratio <- dnorm(df_v$S - d, nf$mu, nf$sig) / dnorm(df_v$S, nf$mu, nf$sig)
      ww <- (df_v$R / df_v$pw) * (df_v$C / nf$cf_C$oof_predictions) * ratio
      sum(ww)^2 / sum(ww^2)
    }))
  attr(out, "rug") <- q_lo - S2
  attr(out, "q")   <- q
  out
}

sve_support_panel <- function(support) {
  rug <- data.frame(x = attr(support, "rug"))
  dr  <- range(support$delta)
  rug <- rug[rug$x >= dr[1] & rug$x <= dr[2], , drop = FALSE]
  ggplot(support, aes(x = delta)) +
    geom_area(aes(y = extrap_lo), fill = "#b30000", alpha = 0.25) +
    geom_line(aes(y = extrap_lo), color = "#b30000", linewidth = 0.6) +
    geom_rug(data = rug, aes(x = x), inherit.aes = FALSE,
             sides = "b", alpha = 0.35, length = grid::unit(0.14, "npc"),
             color = "#b30000") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = expression("Shift in standardized immune marker, " * delta *
                          " (SD units)"),
         y = "Below observed\nmarker support") +
    theme_minimal(base_size = 11.5) +
    theme(panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 9))
}

sve_figure <- function(res, truth = NULL, n_label = NULL, support = NULL) {
  res$Estimator <- factor(ifelse(res$estimator == "tmle", "TMLE", "One-step"),
                          levels = c("TMLE", "One-step"))
  subtitle <- paste0(
    "COVE-calibrated simulated trial",
    if (!is.null(n_label)) paste0(" (", n_label, ")"), ": two-phase marker ",
    "sampling,\nright-censored outcome, cross-fitted super-learner ",
    "nuisances ({nadir}).\n",
    "Shaded band: pointwise 95% Wald confidence intervals (TMLE).")
  p <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dotted",
               color = "grey60") +
    geom_ribbon(data = subset(res, Estimator == "TMLE"),
                aes(x = delta, ymin = lo, ymax = hi),
                fill = "#2c7fb8", alpha = 0.18) +
    geom_line(data = subset(res, Estimator == "TMLE"),
              aes(x = delta, y = sve), color = "#2c7fb8", linewidth = 0.7) +
    geom_point(data = res,
               aes(x = delta, y = sve, shape = Estimator, color = Estimator),
               size = 2.1, stroke = 0.8,
               position = position_dodge(width = 0.06)) +
    scale_color_manual(values = c("TMLE" = "#2c7fb8", "One-step" = "#d95f02")) +
    scale_shape_manual(values = c("TMLE" = 16, "One-step" = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = expression("Shift in standardized immune marker, " * delta *
                       " (SD units)"),
      y = expression("Stochastic-interventional vaccine efficacy, SVE(" *
                       delta * ")"),
      title = "SVE under shifts of a vaccine-induced immune marker",
      subtitle = subtitle,
      caption = paste0("Estimation: {RieszCML} composed Riesz curves; ",
                       "delta method for 1 - \u03c8\u2081(\u03b4)/\u03c8\u2080.")) +
    theme_minimal(base_size = 11.5) +
    theme(
      legend.position = c(0.985, 0.03),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.75), color = NA),
      legend.title = element_blank(),
      legend.margin = margin(0, 4, 2, 4),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30", size = 9.5),
      plot.caption = element_text(color = "grey45", size = 8.5),
      panel.grid.minor = element_blank())
  if (!is.null(truth)) {
    p <- p +
      geom_line(data = truth, aes(x = delta, y = sve,
                                  linetype = "Truth (DGP)"),
                linewidth = 0.5, color = "grey25") +
      scale_linetype_manual(NULL, values = c("Truth (DGP)" = "22"))
  }
  p <- p + coord_cartesian(ylim = c(min(res$lo) - 0.02, 1))
  if (is.null(support)) return(p)

  ## two-panel layout: SVE curve above, marker-support strip below
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("The support panel requires the {patchwork} package.")
  }
  p_main <- p +
    labs(x = NULL, caption = NULL,
         y = expression("SVE(" * delta * ")")) +
    theme(axis.text.x = element_blank())
  p_sup <- sve_support_panel(support) +
    labs(caption = paste0(
      "Bottom: share of counterfactual markers S + \u03b4 below ",
      "the 1st percentile of the observed\nphase-two marker distribution; ",
      "rug ticks mark the shift at which each sampled subject's evaluation\n",
      "point exits the observed support. Estimation: {RieszCML} composed ",
      "Riesz curves + delta method.")) +
    theme(plot.caption = element_text(color = "grey45", size = 8, hjust = 0))
  patchwork::wrap_plots(p_main, p_sup, ncol = 1, heights = c(3.1, 1))
}

## kable table of the analysis results (LaTeX booktabs by default; use
## format = "html" when knitting to HTML).
sve_table <- function(res, truth = NULL, format = "latex") {
  fmt <- function(x) sprintf("%.1f", 100 * x)
  wide <- reshape(res[, c("delta", "estimator", "sve", "lo", "hi")],
                  idvar = "delta", timevar = "estimator", direction = "wide")
  wide <- wide[order(wide$delta), ]
  tab <- data.frame(
    delta   = sprintf("%+.2f", wide$delta),
    tmle    = sprintf("%s (%s, %s)", fmt(wide$sve.tmle),
                      fmt(wide$lo.tmle), fmt(wide$hi.tmle)),
    onestep = sprintf("%s (%s, %s)", fmt(wide$sve.onestep),
                      fmt(wide$lo.onestep), fmt(wide$hi.onestep)))
  cn <- c("$\\delta$ (SD)", "TMLE", "One-step")
  if (!is.null(truth)) {
    tab$truth <- fmt(truth$sve[match(wide$delta, truth$delta)])
    cn <- c(cn, "Truth")
  }
  kableExtra::kbl(
    tab, format = format, booktabs = TRUE, escape = FALSE,
    align = c("r", rep("c", ncol(tab) - 1L)), col.names = cn,
    caption = paste0(
      "Stochastic-interventional vaccine efficacy, SVE($\\delta$), percent, ",
      "with pointwise 95\\% Wald confidence intervals; COVE-calibrated ",
      "simulated trial, cross-fitted super-learner nuisances."))
}

## =============================================================================
## 6. Monte Carlo study (optional; run on a cluster / multicore machine)
## =============================================================================
## Replicates are parallelized with future.apply; nadir's own fold-level
## future_lapply runs sequentially inside each worker under a single-level
## multisession plan, so no nested-future configuration is needed.
## Summarize a stack of per-replicate results (from run_simulation() or from
## job-array output collected by `sve_cove_mc_collect.R`).
summarize_sve_mc <- function(results, scenario = "cove") {
  truth <- get_truth(sort(unique(results$delta)), scenario)
  summ <- do.call(rbind, lapply(
    split(results, interaction(results$delta, results$estimator)),
    function(g) {
      tr <- truth[truth$delta == g$delta[1], ]
      data.frame(delta = g$delta[1], estimator = g$estimator[1],
                 truth = tr$sve, n_reps = nrow(g),
                 mean_est = mean(g$sve),
                 bias = mean(g$sve) - tr$sve,
                 emp_sd = sd(g$sve), mean_se = mean(g$se),
                 rmse = sqrt(mean((g$sve - tr$sve)^2)),
                 coverage = mean(g$lo <= tr$sve & tr$sve <= g$hi))
    }))
  rownames(summ) <- NULL
  summ[order(summ$delta, summ$estimator), ]
}

run_simulation <- function(n = 28000, M = 200, delta_grid = c(-1, 0, 1),
                           scenario = "cove", n_folds = 3, inner_n_folds = 3,
                           seed = 3000) {
  reps <- future.apply::future_lapply(seq_len(M), function(m) {
    set.seed(seed + m)
    df <- sim_data(n, scenario)
    res <- tryCatch(
      sve_analysis(df, delta_grid, n_folds = n_folds,
                   inner_n_folds = inner_n_folds),
      error = function(e) NULL)
    if (is.null(res)) return(NULL)
    attr(res, "support") <- NULL
    res$rep <- m
    res
  }, future.seed = TRUE)
  results <- do.call(rbind, reps)
  list(results = results,
       summary = summarize_sve_mc(results, scenario),
       truth = get_truth(delta_grid, scenario))
}

mc_table <- function(summ, format = "latex") {
  tab <- data.frame(
    delta = sprintf("%+.2f", summ$delta),
    est   = ifelse(summ$estimator == "tmle", "TMLE", "One-step"),
    truth = sprintf("%.3f", summ$truth),
    bias  = sprintf("%.4f", summ$bias),
    emp_sd = sprintf("%.4f", summ$emp_sd),
    mean_se = sprintf("%.4f", summ$mean_se),
    rmse  = sprintf("%.4f", summ$rmse),
    cover = sprintf("%.1f", 100 * summ$coverage))
  # tab <- tab[order(tab$delta, tab$est), ]
  kableExtra::kbl(
    tab, format = format, booktabs = TRUE, escape = FALSE, row.names = FALSE,
    align = c("r", "l", rep("c", 6)),
    col.names = c("$\\delta$", "Estimator", "Truth", "Bias", "Emp.\\ SD",
                  "Mean SE", "RMSE", "Cover (\\%)"),
    caption = paste0(
      "Monte Carlo operating characteristics of the SVE($\\delta$) ",
      "estimators; COVE-calibrated scenario."))
}

## =============================================================================
## 7. Run
## =============================================================================
if (sys.nframe() == 0L) {
  ## Parallel cross-fitting over the outer folds. Worker count: SVE_WORKERS
  ## (if set) overrides SLURM_CPUS_PER_TASK, so a memory-tight allocation can
  ## run fewer workers than CPUs; capped at the number of outer folds --
  ## extra workers only add memory, not speed. Each concurrent multisession
  ## worker is a separate R session costing ~1-1.5 GB at the sizes used here.
  ## SVE_WORKERS=1 skips {future} entirely (recommended for one-off runs).
  ## Do NOT set options(future.globals.maxSize = Inf): the default 500 MB
  ## guard is a tripwire -- if it fires, an abnormally large object is being
  ## shipped to every worker, and raising the limit converts that diagnostic
  ## into a multi-worker out-of-memory kill.
  n_workers <- suppressWarnings(as.integer(
    Sys.getenv("SVE_WORKERS", unset = Sys.getenv("SLURM_CPUS_PER_TASK", "1"))))
  if (isTRUE(n_workers > 1L)) {
    future::plan(future::multisession, workers = min(n_workers, 5L))
  }

  ## --- primary analysis: one COVE-calibrated trial ---------------------------
  ## Grid restricted to [-1, 1.5]: see the positivity note in the header.
  delta_grid <- seq(-1, 1.5, by = 0.25)
  set.seed(505)
  df <- sim_data(28000, "cove")
  cat(sprintf("simulated trial: n = %d, vaccine cases = %d, phase-two = %d\n",
              nrow(df), sum(df$CY[df$A == 1]), sum(df$R)))

  res     <- sve_analysis(df, delta_grid, n_folds = 5, inner_n_folds = 5)
  truth   <- get_truth(delta_grid, "cove")
  support <- attr(res, "support")

  library(here)
  res_dir <- here("sve_cove_cluster", "results")
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

  write.csv(res, file.path(res_dir, "sve_cove_results.csv"), row.names = FALSE)
  saveRDS(support, file.path(res_dir, "sve_cove_support.rds"))   # keeps the rug attribute
  print(merge(res, truth[, c("delta", "sve")], by = "delta",
              suffixes = c("", "_true")), digits = 3)

  p <- sve_figure(res, truth, n_label = "n = 28,000", support = support)
  ggsave(file.path(res_dir, "sve_cove_figure.pdf"), p, width = 7.6, height = 6.3,
         device = cairo_pdf)
  ggsave(file.path(res_dir, "sve_cove_figure.png"), p, width = 7.6, height = 6.3, dpi = 300)
  writeLines(as.character(sve_table(res, truth)),
             file.path(res_dir, "sve_cove_table.tex"))
}
