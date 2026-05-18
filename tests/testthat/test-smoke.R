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
