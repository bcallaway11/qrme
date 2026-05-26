# Smoke tests based on the simulation DGP in simulations/simulations.R.
# These only check that the functions run without error and return
# well-formed output — not that estimates are close to truth.
# Keep n, tau, and MCMC settings small so tests run in a few seconds.

# --- Shared DGP setup ---------------------------------------------------------
local({
  set.seed(42)
  n   <- 100
  tau <- c(0.25, 0.5, 0.75)

  b0Y <- function(u) 1 + 3 * u - u^2
  b1Y <- function(u) exp(u)
  b0T <- function(u) 1 + sqrt(u)
  b1T <- function(u) 0.1 * exp(u)

  X   <- runif(n)
  cop <- copula::claytonCopula(param = 1.5)
  draws <- copula::rCopula(n, cop)
  eY  <- draws[, 1]
  eT  <- draws[, 2]

  Ystar <- b0Y(eY) + b1Y(eY) * X
  Tstar <- b0T(eT) + b1T(eT) * X
  Y <- Ystar + rnorm(n, sd = 0.5)
  T <- Tstar + rnorm(n, sd = 0.5)

  dd <<- data.frame(Y = Y, T = T, X = X)
  tau <<- tau
})

ctrl <- list(
  se            = FALSE,
  n_mix          = 1,
  mcmc_draws    = 60L,
  mcmc_burn_in   = 30L,
  conv_criterion     = "params",
  max_em_iters         = 20L,
  verbose       = FALSE
)

# --- Test 1: qrme -------------------------------------------------------------
test_that("qrme runs and returns a well-formed merr object", {
  fit <- suppressWarnings(do.call(
    qrme,
    c(list(formula = Y ~ X, data = dd, tau = tau), ctrl)
  ))

  expect_s3_class(fit, "merr")
  expect_true(is.matrix(fit$bet))
  expect_equal(nrow(fit$bet), length(tau))
  expect_true(is.numeric(fit$n_iter) && fit$n_iter >= 1)
  expect_equal(fit$conv_criterion, ctrl$conv_criterion)
  expect_equal(fit$tol, 1e-2)
  expect_true(is.numeric(fit$conv_criteria))
  expect_true(is.logical(fit$conv_converged))
  expect_true(all(fit$sig > 0))
  expect_true(all(fit$pi  > 0) && abs(sum(fit$pi) - 1) < 1e-8)
})

test_that("qrme controls finite-mixture output and warns on nonconvergence", {
  fit <- NULL
  expect_output(
    suppressWarnings(fit <- qrme(
      formula = Y ~ X,
      data = dd,
      tau = tau[c(1, 3)],
      n_mix = 2L,
      tol = Inf,
      max_em_iters = 1L,
      mcmc_draws = 40L,
      mcmc_burn_in = 20L,
      verbose = FALSE
    )),
    NA
  )
  expect_s3_class(fit, "merr")
  expect_true(is.numeric(fit$mix_n_iter))
  expect_true(is.logical(fit$mix_converged))

  expect_warning(
    qrme(
      formula = Y ~ X,
      data = dd,
      tau = tau[c(1, 3)],
      n_mix = 1L,
      tol = 0,
      max_em_iters = 1L,
      mcmc_draws = 40L,
      mcmc_burn_in = 20L,
      verbose = FALSE
    ),
    "EM algorithm failed to converge"
  )
})

# --- Test 2: tsme -------------------------------------------------------------
test_that("tsme runs and returns a well-formed result list", {
  t_values   <- mean(dd$T)
  pov_line <- quantile(dd$Y, 0.2)

  res <- suppressWarnings(tsme(
    data      = dd,
    y_formula   = Y ~ X,
    t_formula   = T ~ X,
    tau       = tau,
    t_values     = t_values,
    pov_line   = pov_line,
    y_n_mix     = 1L,
    t_n_mix     = 1L,
    mcmc_draws  = ctrl$mcmc_draws,
    mcmc_burn_in = ctrl$mcmc_burn_in,
    conv_criterion = ctrl$conv_criterion,
    tol       = Inf,
    max_em_iters     = 1L,
    se        = FALSE,
    verbose   = FALSE
  ))

  expect_type(res, "list")
  expect_s3_class(res$me_qyx, "merr")
  expect_s3_class(res$me_qtx, "merr")
  expect_equal(res$tol[["Y"]], res$me_qyx$tol)
  expect_equal(res$tol[["T"]], res$me_qtx$tol)
  expect_equal(res$conv_criterion[["Y"]], ctrl$conv_criterion)
  expect_true(is.numeric(res$conv_criteria))
  expect_true(is.logical(res$conv_converged))
  expect_true(is.numeric(res$mix_n_iter))
  expect_true(is.logical(res$mix_converged))
  expect_true(is.numeric(res$me_cop_param))
  expect_true(is.matrix(res$me_tmat))
  expect_equal(dim(res$me_tmat), c(4L, 4L))
})

# --- Shared tsme result for plot tests ----------------------------------------
# Created once at file scope with minimal settings; shared by autoplot/plot tests.
local({
  tsme_res <<- suppressWarnings(tsme(
    data           = dd,
    y_formula      = Y ~ X,
    t_formula      = T ~ X,
    tau            = tau,
    t_values       = mean(dd$T),
    pov_line       = quantile(dd$Y, 0.2),
    y_n_mix        = 1L,
    t_n_mix        = 1L,
    mcmc_draws     = ctrl$mcmc_draws,
    mcmc_burn_in   = ctrl$mcmc_burn_in,
    conv_criterion = ctrl$conv_criterion,
    tol            = Inf,
    max_em_iters   = 1L,
    se             = FALSE,
    verbose        = FALSE
  ))
})

# --- Test 3: qrme_start_search ------------------------------------------------
test_that("qrme_start_search returns best_fit and a ll-sorted table", {
  out <- suppressWarnings(qrme_start_search(
    formula        = Y ~ X,
    data           = dd,
    tau            = tau,
    n_starts       = 3L,
    mcmc_draws     = ctrl$mcmc_draws,
    mcmc_burn_in   = ctrl$mcmc_burn_in,
    max_em_iters   = ctrl$max_em_iters,
    conv_criterion = ctrl$conv_criterion
  ))

  expect_type(out, "list")
  expect_s3_class(out$best_fit, "merr")
  expect_true(is.data.frame(out$table))
  expect_true(all(c("start", "ll") %in% names(out$table)))
  expect_gte(nrow(out$table), 1L)
  expect_lte(nrow(out$table), 3L)
  if (nrow(out$table) > 1L)
    expect_gte(out$table$ll[1L], out$table$ll[nrow(out$table)])
})

# --- Test 4: qrme_nmix_select -------------------------------------------------
test_that("qrme_nmix_select returns a BIC-sorted comparison table", {
  out <- suppressWarnings(qrme_nmix_select(
    formula        = Y ~ X,
    data           = dd,
    tau            = tau,
    n_mix_vals     = 0:2,
    mcmc_draws     = ctrl$mcmc_draws,
    mcmc_burn_in   = ctrl$mcmc_burn_in,
    max_em_iters   = ctrl$max_em_iters,
    conv_criterion = ctrl$conv_criterion
  ))

  expect_type(out, "list")
  expect_true(is.data.frame(out$table))
  expect_true(all(c("n_mix", "k_me", "ll", "AIC", "BIC") %in% names(out$table)))
  expect_equal(nrow(out$table), 3L)
  expect_lte(out$table$BIC[1L], out$table$BIC[nrow(out$table)])
})

# --- Test 5: tsme_model_select ------------------------------------------------
test_that("tsme_model_select returns a bic-sorted comparison table", {
  out <- suppressWarnings(tsme_model_select(
    data             = dd,
    y_formula        = Y ~ X,
    t_formula        = T ~ X,
    tau              = tau,
    t_values         = mean(dd$T),
    pov_line         = quantile(dd$Y, 0.2),
    copulas          = "gaussian",
    me_distributions = "gaussian",
    mcmc_draws       = ctrl$mcmc_draws,
    mcmc_burn_in     = ctrl$mcmc_burn_in,
    tol              = Inf,
    max_em_iters     = 1L,
    verbose          = FALSE
  ))

  expect_type(out, "list")
  expect_true(is.data.frame(out$table))
  expect_true(all(c("copula", "me_distribution", "y_n_mix", "t_n_mix",
                    "k_params", "ll_y", "ll_t", "ll_cop", "ll",
                    "aic", "bic") %in% names(out$table)))
  expect_equal(nrow(out$table), 1L)
  expect_equal(out$table$copula[1L], "gaussian")
})

# --- Test 6: tmat and upMob ---------------------------------------------------
test_that("tmat returns a numeric transition matrix with correct structure", {
  m <- tmat(dd$Y, dd$T)

  expect_true(is.matrix(m))
  expect_true(is.numeric(m))
  expect_equal(dim(m), c(4L, 4L))
  expect_true(all(m >= -1e-10))
  expect_true(all(abs(colSums(m) - 1) < 1e-8))
})

test_that("upMob returns a numeric probability vector of correct length", {
  mob <- upMob(dd$Y, dd$T)

  expect_true(is.numeric(mob))
  expect_equal(length(mob), 4L)
  expect_true(all(mob >= 0 & mob <= 1))
})

# --- Test 7: autoplot.tsme and plot.tsme -------------------------------------
test_that("autoplot.tsme returns a ggplot object", {
  p <- autoplot(tsme_res)
  expect_s3_class(p, "ggplot")

  p2 <- autoplot(tsme_res, type = "cond_quant", which = "nome")
  expect_s3_class(p2, "ggplot")
})

test_that("plot.tsme runs without error", {
  expect_no_error(plot(tsme_res))
})
