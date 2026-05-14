# Slow tests: accuracy check for qrme and bootstrap smoke test.
# Skipped by default; opt in two ways:
#   1. devtools::test(filter = "slow")          — always runs slow tests
#   2. QRME_RUN_SLOW=1 devtools::test()         — runs everything including slow

skip_if_not(
  nzchar(Sys.getenv("QRME_RUN_SLOW")),
  "Set QRME_RUN_SLOW=1 to include slow tests in a full devtools::test() run"
)

# --- Shared DGP ---------------------------------------------------------------
# Same DGP as simulations/simulations.R: Clayton copula, Gaussian ME (nmix=1).
# No measurement error on Y so qrme should recover coefficients close to truth.
local({
  set.seed(2025)
  n   <- 500
  tau <- seq(0.1, 0.9, length.out = 9)

  b0Y <<- function(u) 1 + 3 * u - u^2
  b1Y <<- function(u) exp(u)
  b0T <<- function(u) 1 + sqrt(u)
  b1T <<- function(u) 0.1 * exp(u)

  X   <- runif(n)
  cop <- copula::claytonCopula(param = 1.5)
  draws <- copula::rCopula(n, cop)
  eY  <- draws[, 1]
  eT  <- draws[, 2]

  Ystar <- b0Y(eY) + b1Y(eY) * X
  Tstar <- b0T(eT) + b1T(eT) * X

  # Small ME on Y only
  me_sd <- 0.3
  Y <- Ystar + rnorm(n, sd = me_sd)
  T <- Tstar + rnorm(n, sd = me_sd)

  dd  <<- data.frame(Y = Y, T = T, X = X)
  tau <<- tau
})

# --- Test 1: qrme accuracy ----------------------------------------------------
test_that("qrme coefficients are close to true values", {
  fit <- qrme(
    formula     = Y ~ X,
    data        = dd,
    tau         = tau,
    nmix        = 1,
    mcmc_draws  = 400L,
    mcmc_burnin = 200L,
    conv_crit   = "params",
    maxit       = 100L,
    se          = FALSE,
    messages    = FALSE
  )

  est   <- fit$bet                        # tau x 2 matrix: intercept, slope
  true  <- cbind(b0Y(tau), b1Y(tau))

  # Loose tolerance: mean absolute error < 0.3 per coefficient across quantiles
  mae_intercept <- mean(abs(est[, 1] - true[, 1]))
  mae_slope     <- mean(abs(est[, 2] - true[, 2]))

  expect_lt(mae_intercept, 0.3)
  expect_lt(mae_slope,     0.3)
})

# --- Test 2: bootstrap smoke test ---------------------------------------------
# Only checks that qrme with se=TRUE runs and attaches SE fields — not accuracy.
# n_boot=3 keeps this as fast as possible while still exercising the code path.
test_that("qrme bootstrap runs and attaches SE fields", {
  fit <- qrme(
    formula     = Y ~ X,
    data        = dd[1:100, ],     # small subset to keep it fast
    tau         = tau[c(1, 5, 9)], # three quantiles only
    nmix        = 1,
    mcmc_draws  = 100L,
    mcmc_burnin = 50L,
    conv_crit   = "params",
    maxit       = 20L,
    se          = TRUE,
    n_boot      = 3L,
    messages    = FALSE
  )

  expect_true(!is.null(fit$sig.se))
  expect_true(!is.null(fit$mu.se))
  expect_true(!is.null(fit$pi.se))
  expect_true(!is.null(fit$bet.se))
  expect_true(all(is.numeric(fit$sig.se)))
  expect_true(all(is.numeric(fit$bet.se)))
})
