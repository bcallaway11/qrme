# =============================================================================
# Title: EM Algorithm for QRME
# Description: Pseudo EM algorithm and helpers for quantile regression with
#              mixture-of-normals measurement error. em.algo() is the outer
#              loop; em.algo.inner() does one E+M step.
# Author: Brant Callaway
# Last update: 2026-05-15
# Date created: 2026-05-07
# =============================================================================

qrme_verbose_level <- function(verbose) {
  if (isTRUE(verbose)) {
    return(1L)
  }
  if (isFALSE(verbose) || is.null(verbose)) {
    return(0L)
  }
  level <- suppressWarnings(as.integer(verbose[1]))
  if (is.na(level) || level < 0L) {
    stop("verbose must be TRUE, FALSE, or a nonnegative integer")
  }
  level
}

qrme_progress <- function(verbose, ..., level = 1L) {
  if (qrme_verbose_level(verbose) >= level) {
    message(..., appendLF = TRUE)
  }
}

fit_normal_mixture <- function(x, m, epsilon = 1e-03, verbose = FALSE) {
  if (m == 1) {
    return(list(
      fit = list(m = 1, lambda = 1, mu = 0, sigma = sd(x)),
      n_iter = 0L,
      loglik = NA_real_,
      converged = TRUE
    ))
  }

  # normalmixEM can fail with "Too many tries!" on awkward data; catch it and
  # return a data-informed fallback so the outer EM can continue.
  result <- tryCatch({
    nm_output <- capture.output({
      nm <- mixtools::normalmixEM(
        x,
        k = m,
        epsilon = epsilon,
        verb = qrme_verbose_level(verbose) >= 3L
      )
    })
    list(nm = nm, nm_output = nm_output, failed = FALSE)
  }, error = function(e) {
    list(nm = NULL, nm_output = character(0), failed = TRUE,
         msg = conditionMessage(e))
  })

  if (result$failed) {
    warning("normalmixEM failed (", result$msg, "); ",
            "using empirical fallback for mixture parameters.", call. = FALSE)
    probs  <- seq(1 / (m + 1), m / (m + 1), length.out = m)
    mu_fb  <- as.numeric(quantile(x, probs))
    mu_fb  <- mu_fb - mean(mu_fb)
    return(list(
      fit = list(m = m, lambda = rep(1 / m, m), mu = mu_fb,
                 sigma = rep(sd(x), m)),
      n_iter = NA_integer_,
      loglik = NA_real_,
      converged = FALSE
    ))
  }

  nm       <- result$nm
  nm_output <- result$nm_output

  if (qrme_verbose_level(verbose) >= 3L && length(nm_output) > 0L) {
    message(paste(nm_output, collapse = "\n"))
  }

  not_converged <- any(grepl("NOT CONVERGENT", nm_output))
  mix_iter <- max(length(nm$all.loglik) - 1L, 0L)
  qrme_progress(
    verbose,
    "  finite mixture iterations: ", mix_iter,
    "; log-likelihood: ", round(nm$loglik, 4),
    "; converged: ", !not_converged,
    level = 2L
  )

  list(
    fit = nm,
    n_iter = mix_iter,
    loglik = nm$loglik,
    converged = !not_converged
  )
}


#' @title em.algo
#'
#' @description A pseudo EM algorithm for quantile regression with measurement error.  The measurement error here follows a mixture of normals.
#' @param betmatguess Initial values for the beta parameters.  This should be
#'  an LxK matrix where L is the number of quantiles and K is the dimension
#'  of the covariates
#' @param tau An L-vector indicating which quantiles have been estimated
#' @param me_distribution The distribution of the measurement error.  "gaussian" is the
#'  default and supports a mixture of normals.  "laplace" is also supported.
#' @param m The number of components of the mixture distribution for the
#'  measurement error
#' @param piguess Starting value for the probabilities of each mixture distribution (should have length equal to k)
#' @param muguess Starting value for the mean of each mixture component (should
#'  have length equal to k)
#' @param sigguess Starting value for the standard deviation of each mixture
#'  component (should have length equal to k)
#' @param mcmc_method Simulation step to use in the EM algorithm. Currently
#'  only \code{"MH"} is supported.
#' @param max_em_iters Maximum number of EM outer iterations. If convergence is not
#'  reached, the estimates from the final iteration are returned (default is 100)
#' @param tol Convergence tolerance.  When \code{NULL} (default), \code{1e-2}
#'  is used for both criteria.  For \code{"params"}, \code{tol} is compared
#'  against \code{max(|theta_new - theta_old| / (1 + |theta_old|))}, the
#'  maximum scaled parameter change across all parameters; this is
#'  dimensionless and application-invariant.  For \code{"loglik"}, it is
#'  compared against the relative log-likelihood change
#'  \code{|ll_new - ll_old| / |ll_old|}.  Tighten to \code{1e-3} or smaller
#'  for higher accuracy at the cost of more iterations.
#' @param conv_criterion Convergence criterion.  \code{"params"} (default) stops
#'  when the maximum scaled parameter change falls below \code{tol} (see
#'  \code{tol} above).  \code{"loglik"} stops when the relative change in the
#'  observed-data log-likelihood, \code{|ll_new - ll_old| / |ll_old|}, falls
#'  below \code{tol}; the log-likelihood is evaluated using Monte Carlo
#'  integration with \code{loglik_draws} draws per outer iteration.
#' @param loglik_draws Number of Monte Carlo draws used to evaluate the
#'  log-likelihood at each outer EM iteration when \code{conv_criterion = "loglik"}.
#'  More draws reduce noise in the convergence check at the cost of speed.
#'  Default is 100.
#' @param conv_patience Integer. Number of consecutive iterations that must
#'  satisfy the convergence criterion before the algorithm stops. Default is
#'  1 (current behaviour). Setting this to 2 requires two back-to-back
#'  iterations below \code{tol}, guarding against false convergence caused by
#'  Monte Carlo noise in the log-likelihood or parameter updates.
#' @param proposal_sd Standard deviation of the Metropolis-Hastings proposal
#'  (random-walk step size). When \code{NULL} (default), initialised to the
#'  SD of the ME mixture from the starting parameters and updated after each
#'  M-step to track the current ME scale:
#'  \code{sqrt(sum(pi_k * (sig_k^2 + mu_k^2)))}. Pass a positive numeric to
#'  fix the proposal at that value for all iterations.
#' @inheritParams qrme
#' @return QRME object
#'
#' @export
em.algo <- function(formula, data,
                    betmatguess, tau,
                    me_distribution = "gaussian",
                    m = 1, piguess = 1, muguess = 0,
                    sigguess = 1, mcmc_method = "MH", tol = NULL,
                    conv_criterion = "params", loglik_draws = 100L,
                    conv_patience = 1L,
                    mcmc_draws = 200, mcmc_burn_in = 100, proposal_sd = NULL, n_cores = 1,
                    max_em_iters = 100, verbose = FALSE) {
    # some checks
    if (me_distribution == "laplace" & mcmc_method != "MH") {
        stop("Laplace distribution only supported with MH algorithm")
    }
    if (me_distribution == "laplace" & m > 1) {
        stop("Laplace distribution only supported with m=1")
    }
    if (mcmc_method == "ImpSamp") {
        stop("mcmc_method = 'ImpSamp' is not currently supported", call. = FALSE)
    }
    if (mcmc_method != "MH") {
        stop('mcmc_method must be "MH"', call. = FALSE)
    }
    if (!conv_criterion %in% c("params", "loglik")) {
        stop('conv_criterion must be "params" or "loglik"')
    }

    # Extract y and x once; needed for default tol computation and/or loglik criterion
    xformula_rhs <- formula
    xformula_rhs[[2]] <- NULL
    x_for_ll <- model.matrix(xformula_rhs, data)
    y_name <- as.character(formula[[2]])
    y_for_ll <- data[, y_name]

    # When proposal_sd is not supplied, initialise it to the SD of the ME
    # mixture implied by the starting parameters, then update it after each
    # M-step to track the current ME scale. Each EM iteration starts a fresh
    # chain, so the proposal should match the current ME distribution rather
    # than a fixed value tied to the total outcome variance.
    # Var[mixture] = sum_k pi_k*(sig_k^2 + mu_k^2) when sum(pi*mu)=0.
    adapt_proposal_sd <- is.null(proposal_sd)
    if (adapt_proposal_sd) {
        proposal_sd <- sqrt(sum(piguess * (sigguess^2 + muguess^2)))
    }

    # Set default tol if not provided by the user
    if (is.null(tol)) {
        if (conv_criterion == "loglik") {
            # Relative log-likelihood change; 1e-2 is deliberately loose as a starting point
            tol <- 1e-2
        } else {
            # Maximum scaled parameter change; 1e-2 matches the loglik default
            # and is application-invariant (dimensionless by construction)
            tol <- 1e-2
        }
    }

    counter <- 1
    criteria <- NA_real_
    ll_old <- NULL        # used only when conv_criterion == "loglik"
    consec_count <- 0L    # consecutive iterations satisfying the convergence criterion
    last_accept_rate <- NA_real_  # updated each iteration; checked at exit

    add_convergence_info <- function(obj, n_iter, converged) {
        if (isFALSE(obj$mix_converged)) {
            warning(
                "Finite mixture model did not converge in the final EM iteration",
                call. = FALSE
            )
        }
        # Warn once if the final-iteration acceptance rate is extreme.
        # For 1-D random-walk MH, ~20-70% is healthy; outside 10-90% suggests
        # proposal_sd is badly mis-scaled and effective sample size is very low.
        if (!is.na(last_accept_rate)) {
            if (last_accept_rate < 0.10) {
                warning(sprintf(
                    "MH acceptance rate was %.0f%% in the final EM iteration. ",
                    last_accept_rate * 100
                ), "proposal_sd may be too large; consider passing a smaller value.",
                call. = FALSE)
            } else if (last_accept_rate > 0.90) {
                warning(sprintf(
                    "MH acceptance rate was %.0f%% in the final EM iteration. ",
                    last_accept_rate * 100
                ), "proposal_sd may be too small; consider passing a larger value.",
                call. = FALSE)
            }
        }
        obj$n_iter <- n_iter
        obj$tol <- tol
        obj$conv_criterion <- conv_criterion
        obj$conv_criteria <- criteria
        obj$conv_converged <- converged
        obj
    }

    qrme_progress(verbose, "QRME EM algorithm")
    qrme_progress(verbose, "  convergence criterion: ", conv_criterion)
    qrme_progress(verbose, "  tolerance: ", signif(tol, 4), level = 2L)
    qrme_progress(verbose, "  max iterations: ", max_em_iters, level = 2L)
    qrme_progress(verbose, "  MCMC draws: ", mcmc_draws, "; burnin: ", mcmc_burn_in, level = 2L)

    # run em algorithm
    while (counter <= max_em_iters) {
        qrme_progress(verbose, "EM iteration ", counter)
        newone <- em.algo.inner(formula, data,
            betmatguess, tau,
            me_distribution,
            m, piguess, muguess,
            sigguess, mcmc_method,
            mcmc_draws = mcmc_draws, mcmc_burn_in = mcmc_burn_in,
            proposal_sd = proposal_sd,
            n_cores = n_cores,
            verbose = verbose
        )
        newbet <- newone$bet
        newpi <- newone$pi
        newmu <- newone$mu
        newsig <- newone$sig
        last_accept_rate <- newone$accept_rate

        qrme_progress(verbose, "  MH acceptance rate: ",
            sprintf("%.1f%%", last_accept_rate * 100), level = 2L)
        qrme_progress(verbose, "  pi: ", paste(signif(newpi, 4), collapse = ", "), level = 2L)
        qrme_progress(verbose, "  mu: ", paste(signif(newmu, 4), collapse = ", "), level = 2L)
        qrme_progress(verbose, "  sig: ", paste(signif(newsig, 4), collapse = ", "), level = 2L)

        if (conv_criterion == "loglik") {
            ll_new <- loglik_raw(y_for_ll, x_for_ll, newbet, tau,
                                 newpi, newmu, newsig, me_distribution, ndraws = loglik_draws)
            qrme_progress(verbose, "  log-likelihood: ", round(ll_new, 4), level = 2L)
            # Skip convergence check on first iteration (no ll_old yet)
            if (!is.null(ll_old)) {
                criteria <- abs(ll_new - ll_old) / abs(ll_old)
                qrme_progress(verbose, "  relative log-likelihood change: ", signif(criteria, 4))
                if (criteria <= tol) {
                    consec_count <- consec_count + 1L
                    if (consec_count >= conv_patience) {
                        qrme_progress(verbose, "EM algorithm converged")
                        return(add_convergence_info(newone, counter, TRUE))
                    }
                } else {
                    consec_count <- 0L
                }
            }
            ll_old <- ll_new
        } else {
            old_params <- c(betmatguess, sigguess, piguess, muguess)
            new_params <- c(newbet,      newsig,   newpi,   newmu)
            criteria <- max(abs(new_params - old_params) / (1 + abs(old_params)))
            qrme_progress(verbose, "  max scaled parameter change: ", signif(criteria, 4))
            if (criteria <= tol) {
                consec_count <- consec_count + 1L
                if (consec_count >= conv_patience) {
                    qrme_progress(verbose, "EM algorithm converged")
                    return(add_convergence_info(newone, counter, TRUE))
                }
            } else {
                consec_count <- 0L
            }
        }

        counter <- counter + 1
        betmatguess <- newbet
        sigguess <- newsig
        muguess <- newmu
        piguess <- newpi
        if (adapt_proposal_sd) {
            proposal_sd <- sqrt(sum(newpi * (newsig^2 + newmu^2)))
            qrme_progress(verbose, "  proposal_sd (next iter): ",
                signif(proposal_sd, 4), level = 2L)
        }
    }
    warning("EM algorithm failed to converge after ", max_em_iters, " iterations", call. = FALSE)
    return(add_convergence_info(newone, max_em_iters, FALSE))
}

#' @title Inner part of EM-algorithm for QR with measurement error
#'
#' @description Does the heavy-lifting of the EM-algorithm for QR with
#'  measurment error
#' @inheritParams em.algo
#' @inheritParams fv.yx
#'
#' @keywords internal
#'
#' @return A list of QR parameters and parameters for mixture of normals for
#'  the measurement error term
#' @export
em.algo.inner <- function(formula, data,
                          betmat, tau,
                          me_distribution = "gaussian",
                          m = 1, pi = 1, mu = 0, sig = 1,
                          mcmc_method = "MH",
                          mcmc_draws = 200, mcmc_burn_in = 100, proposal_sd = NULL, n_cores = 1, verbose = FALSE) {
    xformula <- formula
    xformula[[2]] <- NULL # drop y variable
    X <- model.matrix(xformula, data)
    y_name <- as.character(formula[[2]])
    Y <- data[, y_name]
    k <- ncol(X) # number of x variables
    n <- length(Y)
    qrme_progress(verbose, "  simulating measurement error...")

    if (mcmc_method == "MH") {
        startval <- 0

        edraws <- mh_mcmcC(Y, X,
            startval = startval, mcmc_draws = mcmc_draws, mcmc_burn_in = mcmc_burn_in,
            proposal_sd = proposal_sd, betmat = betmat, me_distribution = me_distribution, m = m,
            pi = pi, mu = mu, sig = sig, tau = tau
        )

        # Average per-observation MH acceptance rate from post-burnin draws.
        # Consecutive equal values indicate a rejection.
        n_kept <- mcmc_draws - mcmc_burn_in
        accept_rate <- mean(vapply(seq_len(n), function(i) {
            draws_i <- edraws[((i - 1L) * n_kept + 1L):(i * n_kept)]
            mean(diff(draws_i) != 0)
        }, numeric(1L)))

        newids <- unlist(lapply(1:n, function(i) rep(i, n_kept)))

        newdta1 <- as.data.frame(cbind(Y = (Y[newids] - edraws), X = X[newids, ], e = edraws))

        # Note: this currently just works for one X; will need to update
        qrme_progress(verbose, "  fitting QR including simulated measurement error...")
        # newdta1 <- do.call(rbind.data.frame, newdta)
        colnames(newdta1) <- c(y_name, colnames(X), "e")

        # thetime <- Sys.time()
        # out  <- quantreg::rq(formla, tau=tau, data=newdta1, method="fn")
        # Sys.time() - thetime
        # out <- quantreg::rq(formla, tau=tau, data=newdta1, method="pfn", weights=rep(1,nrow(newdta1)))
        # out <- rq.fit.pfnb(x=model.matrix(formla, data=data),
        #              y=model.response(model.frame(formla, data=data)),
        #              tau=tau
        newdta1$w <- 1
        out <- quantreg::rq(
            formula = formula,
            tau = tau,
            weights = w,  # w is a column of newdta1, found via model.frame NSE
            data = newdta1,
            method = "pfn"
        )
        # this is part I am not sure about, once you have a new beta then estimate a new sigma??
        # also should probably restrict overall mean of measurement error term to be equal to 0
        if (me_distribution == "laplace") {
            # Laplace M-step: MLE of scale parameter with mu fixed to 0 is
            # mean(|epsilon|). Fitting a normal mixture here would be wrong
            # because it recovers the SD (~b*sqrt(2)) rather than the Laplace
            # scale (b), causing the EM to inflate sigma each iteration.
            qrme_progress(verbose, "  fitting Laplace scale (MLE)...")
            pi_ord  <- 1
            mu_ord  <- 0
            sig_ord <- mean(abs(newdta1$e))
            qrme_progress(verbose, "  sig: ", signif(sig_ord, 4), level = 2L)
            bet_out <- t(coef(out))
            return(list(
                bet = bet_out, m = m, pi = pi_ord, mu = mu_ord, sig = sig_ord,
                mix_n_iter = NA_integer_,
                mix_loglik = NA_real_,
                mix_converged = TRUE,
                accept_rate = accept_rate
            ))
        }

        qrme_progress(verbose, "  fitting finite mixture model...")
        nm_fit <- fit_normal_mixture(newdta1$e, m = m, epsilon = 1e-03, verbose = verbose)
        nm <- nm_fit$fit
        nmorder <- order(nm$mu) # reorder results by mean of each component
        pi_ord  <- nm$lambda[nmorder]
        mu_ord  <- nm$mu[nmorder]
        sig_ord <- nm$sigma[nmorder]
        qrme_progress(verbose, "  mixture pi: ", paste(signif(pi_ord, 4), collapse = ", "), level = 2L)
        qrme_progress(verbose, "  mixture mu: ", paste(signif(mu_ord, 4), collapse = ", "), level = 2L)
        qrme_progress(verbose, "  mixture sig: ", paste(signif(sig_ord, 4), collapse = ", "), level = 2L)

        # Enforce mean-zero constraint: sum(pi_k * mu_k) = 0.
        # Without this, the intercept of beta and the ME mean drift together
        # along a flat ridge in the observed-data likelihood (they are not
        # separately identified). A uniform shift by delta is exact in the
        # observed-data likelihood but is an approximation to the true
        # constrained M-step optimum (which requires component-specific shifts
        # weighted by pi_k * sigma_k^2 / n_k). The approximation is exact when
        # all components share the same variance; in practice it is close and
        # the algorithm still converges to the correct constrained solution.
        delta        <- sum(pi_ord * mu_ord)
        mu_ord       <- mu_ord - delta
        bet_out      <- t(coef(out))
        bet_out[, 1] <- bet_out[, 1] + delta  # column 1 is always the intercept

        return(list(
            bet = bet_out, m = m, pi = pi_ord, mu = mu_ord, sig = sig_ord,
            mix_n_iter = nm_fit$n_iter,
            mix_loglik = nm_fit$loglik,
            mix_converged = nm_fit$converged,
            accept_rate = accept_rate
        ))
    } else if (mcmc_method == "ImpSamp") {
        # importance sampling
        edraws <- rnorm((mcmc_draws * n), 0, proposal_sd)

        newids <- unlist(lapply(1:n, function(i) rep(i, mcmc_draws))) # just replicates Y and X over and over

        newdta1 <- as.data.frame(cbind(Y = (Y[newids] - edraws), X = X[newids, ], e = edraws)) # prepopulate some fields in dataset

        # compute weights using importance sampling
        newdta1$w <- imp_sampC(
            Y = Y[newids], X = X[newids, ], V = edraws, mcmc_draws = mcmc_draws, proposal_sd = proposal_sd,
            betmat = betmat, m = m, pi = pi, mu = mu, sig = sig, tau = tau
        ) # but use original versions of the data (not adjusted for measurement errors) to compute weights
        newdta1$w <- sapply(1:length(edraws), function(i) max(1e-05, newdta1$w[i])) # drop negative weights (not many of these...)
        # run weighted quantile regression
        out <- quantreg::rq(formula, tau = tau, data = newdta1, method = "fn", weights = newdta1$w)
        qrme_progress(verbose, "  fitting finite mixture model...")
        ##
        # need to make actual draws here
        ##
        # set up newids again
        # this is going to have a different length from
        # newids above because we are going to use the length of tau rather than number of iterations
        newids <- unlist(lapply(1:n, function(i) rep(i, length(tau))))
        # draws of Y^* using X'\beta(U) and having U take all possible values of tau
        Ystar <- c(t(X %*% coef(out)))
        # finally, recover draws of measurement error
        U <- Y[newids] - Ystar
        nm_fit <- fit_normal_mixture(U, m = m, epsilon = 1e-03, verbose = verbose)
        nm <- nm_fit$fit
        nmorder <- order(nm$mu) # reorder results by mean of each component
        pi_ord  <- nm$lambda[nmorder]
        mu_ord  <- nm$mu[nmorder]
        sig_ord <- nm$sigma[nmorder]
        qrme_progress(verbose, "  mixture pi: ", paste(signif(pi_ord, 4), collapse = ", "), level = 2L)
        qrme_progress(verbose, "  mixture mu: ", paste(signif(mu_ord, 4), collapse = ", "), level = 2L)
        qrme_progress(verbose, "  mixture sig: ", paste(signif(sig_ord, 4), collapse = ", "), level = 2L)

        # Enforce mean-zero constraint (same as MH branch; see comment there)
        delta        <- sum(pi_ord * mu_ord)
        mu_ord       <- mu_ord - delta
        bet_out      <- t(coef(out))
        bet_out[, 1] <- bet_out[, 1] + delta

        return(list(
            bet = bet_out, m = m, pi = pi_ord, mu = mu_ord, sig = sig_ord,
            mix_n_iter = nm_fit$n_iter,
            mix_loglik = nm_fit$loglik,
            mix_converged = nm_fit$converged
        ))
    } else {
        stop("provided mcmc_method not supported")
    }
}


##############################################################
##
## Code for computing densities used in MH algorithm
##
## Note: most of this has been replaced by C++ code to
## improve performance, but keeping it here to have R
## versions of these.
##
##############################################################
## X should be n x k
## the density of y conditional on x
#' @title fy.x
#'
#' @description return the density of Y (for particular value y) conditional
#'  on X (which can include n observations) when Q(Y|X) has been estimated
#'  using QR.  This is used in our simulation approach.
#' @param y particular value of y to estimate f(y|x)
#' @param betmat LxK matrix of parameter values with L the number of quantiles
#'  and K the dimension of the covariates
#' @param X An nxK matrix with n the number of observations of X
#' @param tau an L-vector containing the quantile at which Q(Y|X) was estimated
#'
#' @return An nx1 vector that contains f(y|X)
#'
#' @export
fy.x <- function(y, betmat, X, tau) {
    X <- as.matrix(X)

    fout <- apply(X, 1, FUN = function(x) {
        ## take a particular row of X
        x <- as.matrix(x)

        ## the index for the "X" part
        xtb <- t(x) %*% t(betmat)

        ## figure out if we are in an "inner case" (i.e. standard case)
        ## or in one of the tails
        ul.pos <- tail(which(xtb - y <= 0), 1) ## position of ul (in text)
        uu.pos <- head(which(xtb - y > 0), 1) ## position of uu (in text)

        ## code to handle tails uniquely
        lam1 <- 1 - min(tau)
        lam2 <- max(tau)
        if (length(uu.pos) > 0) { ## add extra check for some missing cases; don't
            ## need for other side because of the ordering
            if (uu.pos == 1) { ## we are way on left tail
                return(min(tau) * lam1 * exp(lam1 * (y - t(x) %*% betmat[1, ])))
            }
        }
        if (ul.pos == length(tau)) { ## we are way on the right tail
            return((1 - max(tau)) * lam2 * exp(-lam2 * (y - t(x) %*% betmat[length(tau), ])))
        }
        ## standard case ("inner case")
        (tau[uu.pos] - tau[ul.pos]) / t(x) %*% (betmat[uu.pos, ] - betmat[ul.pos, ])
    })
    fout
}


#' @title fv.yx
#'
#' @description Return the density of the measurement error conditional
#'  on y and x; this takes as given some QR parameters from Y* (the true
#'  outcome) conditional on X.  Here, we also presume that the distribution
#'  of the measurement error is a mixture of normal distributions
#' @param v A particular value of the measurement error to estimate f(v|y,x)
#' @inheritParams fy.x
#' @param m The dimension of the measurement error
#' @param pi The probability of each mixture component (should have length
#'  equal to m)
#' @param mu The mean of each mixture component (should have length equal
#'  to m)
#' @param sig The standard deviation of each mixture component (should have
#'  length equal to m)
#' @param Y An nx1 vector of outcomes
#' @param X An nxK matrix of covariates
#' @param tau an L-vector of all the quantiles where betas were estimated
#' @return n x 1 vector of f(v|Y,X)
#'
#' @keywords internal
#' @export
fv.yx <- function(v, betmat, m, pi, mu, sig, Y, X, tau) {
    fy.xvals <- sapply(1:length(Y), function(i) {
        fy.x(y = (Y[i] - v), X = t(X[i, ]), betmat = betmat, tau = tau)
    })
    fy.xvals * fv(v, m, pi, mu, sig)
}


#' @title fv
#'
#' @description The distribution of the measurement error using a mixture of
#'  normal distributions
#'
#' @inheritParams fv.yx
#' @return scalar f(v)
#' @keywords internal
#' @export
fv <- function(v, m = 1, pi = 1, mu = 0, sig = 1) {
    ## mixture of normals
    sum(sapply(1:m, function(i) {
        pi[i] / sig[i] * dnorm((v - mu[i]) / sig[i])
    }))
}


#' @title mh_mcmc
#' @description A Metropolis-Hastings algorithm for drawing measurment errors.
#'
#' @param startval The first value in the markov chain
#' @param mcmc_draws The total number of measurement error draws to make
#' @param mcmc_burn_in The number of draws to drop
#' @param proposal_sd Standard deviation of the random-walk MH proposal. Passed
#'  from \code{em.algo} after automatic scaling; see \code{\link{em.algo}}.
#' @inheritParams fv.yx
#' @param y particular value of y
#' @param x particular value of x
#' @return vector of draws of measurement error
#' @export
mh_mcmc <- function(startval = 0, mcmc_draws = 200, mcmc_burn_in = 100, proposal_sd = NULL, betmat, m, pi, mu, sig, y, x, tau) {
    x <- t(x)
    out <- rep(NA, mcmc_draws)
    out[1] <- startval
    for (i in 2:mcmc_draws) {
        trialval <- out[i - 1] + rnorm(1, sd = proposal_sd)
        fvold <- fv.yx(out[i - 1], betmat, m, pi, mu, sig, y, x, tau)
        fvnew <- fv.yx(trialval, betmat, m, pi, mu, sig, y, x, tau)
        if (fvnew > fvold) {
            out[i] <- trialval
        } else {
            if ((fvnew / fvold) >= runif(1)) {
                out[i] <- trialval
            } else {
                out[i] <- out[i - 1]
            }
        }
    }
    return(tail(out, mcmc_draws - mcmc_burn_in))
}
