# Smoke tests based on the simulation DGP in simulations/simulations.R.
# These only check that the functions run without error and return
# well-formed output — not that estimates are close to truth.
# Keep n, tau, and MCMC settings small so tests run in a few seconds.

# --- Shared DGP setup ---------------------------------------------------------
local({
  set.seed(42)
  n   <- 150
  tau <- seq(0.1, 0.9, length.out = 5)

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
  nmix          = 1,
  mcmc_draws    = 100L,
  mcmc_burnin   = 50L,
  conv_crit     = "params",
  maxit         = 20L,
  verbose       = FALSE
)

# --- Test 1: qrme -------------------------------------------------------------
test_that("qrme runs and returns a well-formed merr object", {
  fit <- do.call(qrme, c(list(formula = Y ~ X, data = dd, tau = tau), ctrl))

  expect_s3_class(fit, "merr")
  expect_true(is.matrix(fit$bet))
  expect_equal(nrow(fit$bet), length(tau))
  expect_equal(length(fit$mu),  ctrl$nmix)
  expect_equal(length(fit$sig), ctrl$nmix)
  expect_equal(length(fit$pi),  ctrl$nmix)
  expect_true(is.numeric(fit$n_iter) && fit$n_iter >= 1)
  expect_equal(fit$conv_crit, ctrl$conv_crit)
  expect_equal(fit$tol, sqrt(length(tau) * ncol(model.matrix(Y ~ X, dd)) + 3 * ctrl$nmix) * 0.1)
  expect_true(is.numeric(fit$conv_criteria))
  expect_true(is.logical(fit$conv_converged))
  expect_true(is.numeric(fit$mix_n_iter))
  expect_true(is.logical(fit$mix_converged))
  expect_true(all(fit$sig > 0))
  expect_true(all(fit$pi  > 0) && abs(sum(fit$pi) - 1) < 1e-8)
})

test_that("qrme controls finite-mixture output and reports progress", {
  fit <- NULL
  expect_output(
    fit <- qrme(
      formula = Y ~ X,
      data = dd,
      tau = tau[c(2, 4)],
      nmix = 2L,
      tol = Inf,
      maxit = 1L,
      mcmc_draws = 40L,
      mcmc_burnin = 20L,
      verbose = FALSE
    ),
    NA
  )
  expect_s3_class(fit, "merr")
  expect_true(is.numeric(fit$mix_n_iter))
  expect_true(is.logical(fit$mix_converged))

  expect_message(
    qrme(
      formula = Y ~ X,
      data = dd,
      tau = tau[c(2, 4)],
      nmix = 1L,
      tol = Inf,
      maxit = 1L,
      mcmc_draws = 40L,
      mcmc_burnin = 20L,
      verbose = TRUE
    ),
    "QRME EM algorithm"
  )

  expect_warning(
    qrme(
      formula = Y ~ X,
      data = dd,
      tau = tau[c(2, 4)],
      nmix = 1L,
      tol = 0,
      maxit = 1L,
      mcmc_draws = 40L,
      mcmc_burnin = 20L,
      verbose = FALSE
    ),
    "EM algorithm failed to converge"
  )
})

# --- Test 2: tsme -------------------------------------------------------------
test_that("tsme runs and returns a well-formed result list", {
  tvals   <- mean(dd$T)
  povline <- quantile(dd$Y, 0.2)

  res <- tsme(
    data      = dd,
    Yformla   = Y ~ X,
    Tformla   = T ~ X,
    tau       = tau,
    tvals     = tvals,
    povline   = povline,
    Ynmix     = 1L,
    Tnmix     = 1L,
    mcmc_draws  = ctrl$mcmc_draws,
    mcmc_burnin = ctrl$mcmc_burnin,
    conv_crit = ctrl$conv_crit,
    se        = FALSE,
    verbose   = FALSE
  )

  expect_type(res, "list")
  expect_s3_class(res$meQyx, "merr")
  expect_s3_class(res$meQtx, "merr")
  expect_equal(res$tol[["Y"]], res$meQyx$tol)
  expect_equal(res$tol[["T"]], res$meQtx$tol)
  expect_equal(res$conv_crit[["Y"]], ctrl$conv_crit)
  expect_true(is.numeric(res$conv_criteria))
  expect_true(is.logical(res$conv_converged))
  expect_true(is.numeric(res$mix_n_iter))
  expect_true(is.logical(res$mix_converged))
  expect_true(is.numeric(res$meCopParam))
  expect_true(res$meCopParam >= 0 && res$meCopParam <= 1)
  expect_true(is.matrix(res$meTmat))
  expect_equal(dim(res$meTmat), c(4L, 4L))
})
