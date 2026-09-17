#' Estimate the MSE of Benchmarked Small Area Estimates
#'
#' @description Computes the mean squared error (MSE) of small area estimates benchmarked
#' with [sae_benchmarking()], for the ratio, difference, or optimum method,
#' when the original estimates are EBLUPs under the Fay-Herriot model or
#' hierarchical Bayes (HB) estimates.
#'
#' @param object An object of class `"sae_benchmarking"` returned by
#'   [sae_benchmarking()]. `vardir` must have been supplied in that call when
#'   `estimator = "eblup"`.
#' @param estimator Character. `"eblup"` or `"hb"`.
#' @param z Design matrix (areas in rows, including the intercept column if
#'   any) used to fit the Fay-Herriot model. Required when
#'   `estimator = "eblup"`.
#' @param sigma2_v Estimated random effect variance from the original
#'   Fay-Herriot fit. Required when `estimator = "eblup"`.
#' @param fitting_method Character. `"REML"` or `"PR"` (Prasad-Rao moment
#'   method), matching how `sigma2_v` was estimated. Used only for
#'   `method = "difference"` under `estimator = "eblup"`.
#' @param theta_hb Numeric vector of HB estimates (posterior means) in the
#'   same order as the areas of `object`. Used when `estimator = "hb"`.
#' @param V_hb Numeric vector of HB posterior variances in the same order as
#'   the areas of `object`. Used when `estimator = "hb"`.
#' @param posterior_draws Matrix of MCMC draws of the unbenchmarked small
#'   area means (areas in rows, draws in columns). Used when
#'   `estimator = "hb"` and `theta_hb` or `V_hb` is not supplied.
#' @param B Positive integer. Number of bootstrap replicates, used for the
#'   ratio and optimum methods under `estimator = "eblup"`.
#' @param seed Optional integer seed for the bootstrap. The random number
#'   generator state of the session is restored on exit.
#'
#' @return A list of class `"mse_benchmarking"` with elements:
#' * `method`: the benchmarking method of `object`.
#' * `estimator`: `"eblup"` or `"hb"`.
#' * `estimate`: named numeric vector of benchmarked estimates.
#' * `mse`: named numeric vector of estimated MSEs.
#' * `B`: number of bootstrap replicates requested, or `NA` when no
#'   bootstrap was used.
#' * `n_failed`: number of discarded bootstrap replicates, or `NA` when no
#'   bootstrap was used.
#'
#' @seealso [sae_benchmarking()]
#'
#' @examples
#' # EBLUP, difference benchmarking (closed form)
#' bm_diff <- sae_benchmarking(method = "difference", direct = "direct",
#'                             weight = "weight", estimate = "eblup",
#'                             mse = "mse", vardir = "vardir",
#'                             data = data_eblup)
#' Z <- cbind(1, data_eblup$z)
#' s2v <- attr(data_eblup, "sigma2_v")
#' res_diff <- mse_benchmarking(bm_diff, estimator = "eblup", z = Z,
#'                              sigma2_v = s2v, fitting_method = "REML")
#' head(res_diff$mse)
#'
#' @export
mse_benchmarking <- function(object, estimator = c("eblup", "hb"),
                             z = NULL, sigma2_v = NULL,
                             fitting_method = c("REML", "PR"),
                             theta_hb = NULL, V_hb = NULL,
                             posterior_draws = NULL,
                             B = 1000, seed = NULL) {

  estimator      <- match.arg(estimator)
  fitting_method <- match.arg(fitting_method)
  .validate_mse_inputs(object, estimator, z, sigma2_v, posterior_draws,
                       theta_hb, V_hb, B)

  method   <- object$method
  B_out    <- NA_integer_
  n_failed <- NA_integer_

  if (estimator == "eblup") {
    z <- as.matrix(z)
    if (method == "difference") {
      mse_vec <- .mse_eblup_difference(object, z, sigma2_v, fitting_method)
    } else {
      boot     <- .mse_eblup_bootstrap(object, z, sigma2_v, B, seed)
      mse_vec  <- boot$mse
      B_out    <- as.integer(B)
      n_failed <- boot$n_failed
    }
    estimate_out <- object$estimates$Estimate
  } else {
    hb           <- .mse_benchmarking_hb(object, theta_hb = theta_hb,
                                         V_hb = V_hb,
                                         posterior_draws = posterior_draws)
    mse_vec      <- hb$mse_hb
    estimate_out <- hb$estimate_hb
  }

  mse_vec      <- as.numeric(mse_vec)
  estimate_out <- as.numeric(estimate_out)
  names(mse_vec)      <- object$inputs$area_names
  names(estimate_out) <- object$inputs$area_names

  result <- list(method = method, estimator = estimator,
                 estimate = estimate_out, mse = mse_vec,
                 B = B_out, n_failed = n_failed)
  class(result) <- "mse_benchmarking"
  result
}

#' @export
print.mse_benchmarking <- function(x, ...) {
  cat("Benchmarked MSE estimates (method = '", x$method, "', estimator = '",
      x$estimator, "')\n", sep = "")
  if (!is.na(x$B)) {
    cat("Bootstrap replicates:", x$B, "(failed:", x$n_failed, ")\n")
  }
  cat("\n")
  print(data.frame(Estimate = x$estimate, MSE = x$mse,
                   row.names = names(x$mse)), ...)
  invisible(x)
}

#' @noRd
#' @noRd
.validate_mse_inputs <- function(object, estimator, z, sigma2_v,
                                 posterior_draws, theta_hb, V_hb, B) {

  if (!inherits(object, "sae_benchmarking")) {
    stop("'object' must be an object of class 'sae_benchmarking'.",
         call. = FALSE)
  }

  D <- nrow(object$estimates)

  if (estimator == "eblup") {

    if (is.null(object$inputs$vardir) || anyNA(object$inputs$vardir)) {
      stop("'vardir' must be supplied in sae_benchmarking() when ",
           "estimator = 'eblup'.", call. = FALSE)
    }
    if (is.null(z) || is.null(sigma2_v)) {
      stop("For estimator = 'eblup', arguments 'z' and 'sigma2_v' are ",
           "required.", call. = FALSE)
    }
    z <- as.matrix(z)
    if (!is.numeric(z) || nrow(z) != D) {
      stop("'z' must be a numeric matrix with ", D, " rows.", call. = FALSE)
    }
    if (!is.numeric(sigma2_v) || length(sigma2_v) != 1L || is.na(sigma2_v) ||
        sigma2_v < 0) {
      stop("'sigma2_v' must be a single non-negative number.", call. = FALSE)
    }
    if (object$method != "difference" &&
        (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1)) {
      stop("'B' must be a positive integer.", call. = FALSE)
    }

  } else {

    has_summary <- !is.null(theta_hb) && !is.null(V_hb)

    if (!has_summary && is.null(posterior_draws)) {
      stop("For estimator = 'hb', supply 'theta_hb' and 'V_hb' ",
           "(posterior means and variances), or 'posterior_draws'.",
           call. = FALSE)
    }

    if (has_summary) {
      if (length(theta_hb) != D || length(V_hb) != D) {
        stop("'theta_hb' and 'V_hb' must both have length ", D, ".",
             call. = FALSE)
      }
      if (anyNA(V_hb) || any(V_hb < 0)) {
        stop("'V_hb' must contain non-negative values.", call. = FALSE)
      }
      if (!isTRUE(all.equal(as.numeric(theta_hb),
                            as.numeric(object$inputs$estimate)))) {
        warning("'theta_hb' differs from the 'estimate' used in ",
                "sae_benchmarking().", call. = FALSE)
      }
    } else {
      posterior_draws <- as.matrix(posterior_draws)
      if (nrow(posterior_draws) != D || ncol(posterior_draws) < 2L) {
        stop("'posterior_draws' must have ", D, " rows and at least two ",
             "columns.", call. = FALSE)
      }
    }
  }

  invisible(TRUE)
}

#' @noRd
.fh_gamma <- function(sigma2_v, vardir) sigma2_v / (sigma2_v + vardir)

#' @noRd
.fh_g1 <- function(sigma2_v, vardir) .fh_gamma(sigma2_v, vardir) * vardir

#' @noRd
.fh_sigma_beta <- function(z, sigma2_v, vardir) {
  w <- 1 / (sigma2_v + vardir)
  solve(t(z * w) %*% z)
}

#' @noRd
.fh_g2 <- function(z, sigma2_v, vardir) {
  gamma_i     <- .fh_gamma(sigma2_v, vardir)
  Sigma_beta  <- .fh_sigma_beta(z, sigma2_v, vardir)
  h_ii        <- rowSums((z %*% Sigma_beta) * z)
  (1 - gamma_i)^2 * h_ii
}

#' @noRd
.fh_g3 <- function(sigma2_v, vardir, fitting_method = c("REML", "PR")) {
  fitting_method <- match.arg(fitting_method)
  D     <- length(vardir)
  denom <- (sigma2_v + vardir)^3
  # Asymptotic variance of sigma2_v_hat: differs between REML and the PR moment method
  v_sigma2 <- if (fitting_method == "REML") {
    2 / sum((sigma2_v + vardir)^-2)
  } else {
    (2 / D^2) * sum((sigma2_v + vardir)^2)
  }
  (vardir^2 / denom) * v_sigma2
}

#' @noRd
.fh_mse_h <- function(z, sigma2_v, vardir, fitting_method) {
  g1 <- .fh_g1(sigma2_v, vardir)
  g2 <- .fh_g2(z, sigma2_v, vardir)
  g3 <- .fh_g3(sigma2_v, vardir, fitting_method)
  g1 + g2 + 2 * g3
}

#' @noRd
.mse_eblup_difference <- function(object, z, sigma2_v, fitting_method) {

  vardir  <- object$inputs$vardir
  weight  <- object$inputs$weight
  gamma_i <- .fh_gamma(sigma2_v, vardir)

  mse_H      <- .fh_mse_h(z, sigma2_v, vardir, fitting_method)
  Sigma_beta <- .fh_sigma_beta(z, sigma2_v, vardir)
  h_mat      <- z %*% Sigma_beta %*% t(z)

  term1 <- sum((weight^2) * (1 - gamma_i) * vardir)
  wg    <- weight * (1 - gamma_i)
  term2 <- as.numeric(t(wg) %*% h_mat %*% wg)
  g4    <- term1 - term2

  mse_H + g4
}

#' @noRd
.mse_eblup_bootstrap <- function(object, z, sigma2_v, B, seed = NULL) {

  if (!is.null(seed)) {
    return(withr::with_seed(
      seed,
      .mse_eblup_bootstrap(object, z, sigma2_v, B, seed = NULL)
    ))
  }

  method     <- object$method
  vardir     <- object$inputs$vardir
  weight     <- object$inputs$weight
  direct     <- object$inputs$direct
  phi_source <- object$inputs$phi_source
  D          <- length(vardir)

  need_mse <- method == "optimum" && identical(phi_source, "mse")

  w        <- 1 / (sigma2_v + vardir)
  Sigma_b  <- solve(t(z * w) %*% z)
  beta_hat <- as.numeric(Sigma_b %*% t(z * w) %*% direct)
  mu_hat   <- as.numeric(z %*% beta_hat)

  sq_err   <- matrix(NA_real_, nrow = B, ncol = D)
  n_failed <- 0L

  for (b in seq_len(B)) {

    theta_star  <- mu_hat + stats::rnorm(D, 0, sqrt(sigma2_v))
    direct_star <- theta_star + stats::rnorm(D, 0, sqrt(vardir))

    fit_star <- tryCatch(
      if (need_mse) {
        sae::mseFH(direct_star ~ z - 1, vardir, method = "REML")
      } else {
        sae::eblupFH(direct_star ~ z - 1, vardir, method = "REML")
      },
      error = function(e) NULL
    )

    if (is.null(fit_star)) {
      n_failed <- n_failed + 1L
      next
    }

    if (need_mse) {
      converged  <- isTRUE(fit_star$est$fit$convergence)
      eblup_star <- as.numeric(fit_star$est$eblup)
      phi_star   <- as.numeric(fit_star$mse)
    } else {
      converged  <- isTRUE(fit_star$fit$convergence)
      eblup_star <- as.numeric(fit_star$eblup)
      phi_star   <- vardir
    }

    if (!converged) {
      n_failed <- n_failed + 1L
      next
    }

    target_star <- sum(weight * direct_star)
    agg_star    <- sum(weight * eblup_star)
    gap_star    <- target_star - agg_star

    bench_star <- switch(
      method,
      ratio   = eblup_star * (target_star / agg_star),
      optimum = {
        lambda_star <- (phi_star * weight) / sum(phi_star * weight^2)
        eblup_star + lambda_star * gap_star
      }
    )

    sq_err[b, ] <- (bench_star - theta_star)^2
  }

  if (n_failed == B) {
    stop("All bootstrap replicates failed to fit the Fay-Herriot model.",
         call. = FALSE)
  }
  if (n_failed > 0L) {
    warning(n_failed, " of ", B, " bootstrap replicates failed and were ",
            "discarded.", call. = FALSE)
  }

  list(mse = colMeans(sq_err, na.rm = TRUE), n_failed = n_failed)
}

#' @noRd
.mse_benchmarking_hb <- function(object, theta_hb = NULL, V_hb = NULL,
                                 posterior_draws = NULL) {

  method <- object$method
  weight <- object$inputs$weight
  target <- object$target

  if (is.null(theta_hb) || is.null(V_hb)) {
    posterior_draws <- as.matrix(posterior_draws)
    if (is.null(theta_hb)) theta_hb <- rowMeans(posterior_draws)
    if (is.null(V_hb))     V_hb     <- apply(posterior_draws, 1, stats::var)
  }
  theta_hb <- as.numeric(theta_hb)
  V_hb     <- as.numeric(V_hb)

  agg_hb <- sum(weight * theta_hb)
  gap_hb <- target - agg_hb

  theta_bhb <- switch(
    method,
    ratio      = theta_hb * (target / agg_hb),
    difference = theta_hb + gap_hb,
    optimum    = {
      phi <- switch(object$inputs$phi_source,
                    mse    = object$inputs$mse,
                    vardir = object$inputs$vardir)
      lambda <- (phi * weight) / sum(phi * weight^2)
      theta_hb + lambda * gap_hb
    }
  )

  list(estimate_hb = theta_bhb, mse_hb = V_hb + (theta_hb - theta_bhb)^2)
}
