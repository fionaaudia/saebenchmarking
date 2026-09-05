#' Estimate the mean squared error of benchmarked small area estimates
#'
#' Estimates the MSE of the benchmarked estimates produced by
#' [sae_benchmarking()]. Two approaches are available depending on the
#' information supplied: a residual bootstrap for model-based
#' estimators whose variance parameters are estimated (e.g. EBLUP/EB),
#' following Butar & Lahiri (2003) as extended by Steorts & Ghosh
#' (2013); or direct simulation from posterior draws for Hierarchical
#' Bayes estimators, following the general decomposition of
#' Datta, Ghosh, Steorts & Maples (2011) and Bell, Datta & Ghosh
#' (2013).
#'
#' @param bench_result A list returned by [sae_benchmarking()].
#' @param posterior_samples Optional numeric matrix with `G` rows
#'   (posterior draws) and one column per area, containing draws of
#'   the small area parameter *before* benchmarking. If supplied, the
#'   posterior-simulation approach is used; otherwise the bootstrap
#'   approach is used.
#' @param B Integer. Number of bootstrap replications. Only used when
#'   `posterior_samples` is `NULL`.
#' @param seed Optional integer seed for the bootstrap.
#'
#' @return A list with elements `method`, `approach` (`"bootstrap"` or
#'   `"posterior"`), `MSE` (a data frame with one row per area), and
#'   additional diagnostic elements specific to the approach used.
#'
#' @seealso [sae_benchmarking()] for computing the benchmarked point
#'   estimates required by this function.
#'
#' @examples
#' est <- c(10.2, 8.7, 12.1)
#' mse_val <- c(0.5, 0.4, 0.6)
#' w <- c(0.3, 0.4, 0.3)
#' direct <- c(10.5, 8.2, 12.8)
#' vardir <- c(0.6, 0.5, 0.7)
#'
#' hasil <- sae_benchmarking(method = "difference", estimate = est,
#'                         mse = mse_val, weight = w, direct = direct,
#'                         vardir = vardir)
#' mse_benchmarking(hasil, B = 200, seed = 1)
#'
#' @export
mse_benchmarking <- function(bench_result, posterior_samples = NULL,
                              B = 1000, seed = NULL) {

  if (!is.list(bench_result) || is.null(bench_result$method))
    stop("bench_result must be the output of sae_benchmarking().")

  if (!is.null(posterior_samples)) {
    result <- .mse_via_posterior(bench_result, posterior_samples)
  } else {
    if (!is.numeric(B) || length(B) != 1 || B < 1)
      stop("B must be a single positive integer.")
    result <- .mse_via_bootstrap(bench_result, B = B, seed = seed)
  }

  result
}

#' @keywords internal
#' @noRd
.mse_via_bootstrap <- function(bench_result, B = 1000, seed = NULL) {

  method     <- bench_result$method
  inputs     <- bench_result$inputs
  estimate   <- inputs$estimate
  mse        <- inputs$mse
  weight     <- inputs$weight
  direct     <- inputs$direct
  vardir     <- inputs$vardir
  area_names <- inputs$area_names
  target     <- bench_result$benchmark$target
  D          <- length(estimate)

  if (!is.null(seed)) set.seed(seed)

  have_resid_pool <- !is.null(direct) && !is.null(vardir) &&
    any(!is.na(direct) & !is.na(vardir))

  if (have_resid_pool) {
    ok <- !is.na(direct) & !is.na(vardir)
    std_resid <- (direct[ok] - estimate[ok]) / sqrt(vardir[ok])
    std_resid <- std_resid - mean(std_resid)
    draw_error <- function() sample(std_resid, size = D, replace = TRUE) * sqrt(mse)
  } else {
    warning("No usable 'direct'/'vardir' residual pool found; falling back ",
            "to a parametric bootstrap using N(0, mse) pseudo-errors.")
    draw_error <- function() stats::rnorm(D, mean = 0, sd = sqrt(mse))
  }

  boot_estimation <- matrix(NA_real_, nrow = B, ncol = D)
  for (b in seq_len(B)) {
    estimate_star <- estimate + draw_error()
    point_star    <- .compute_benchmark_point(method, estimate_star, mse, weight, target)
    boot_estimation[b, ] <- point_star$SAE
  }

  var_boot  <- apply(boot_estimation, 2, stats::var, na.rm = TRUE)
  MSE_bench <- mse + var_boot

  list(
    method   = method,
    approach = "bootstrap",
    MSE      = data.frame(MSE = MSE_bench, row.names = area_names),
    B        = B,
    n_failed = sum(is.na(boot_estimation[, 1]))
  )
}

#' @keywords internal
#' @noRd
.mse_via_posterior <- function(bench_result, posterior_samples) {

  method     <- bench_result$method
  inputs     <- bench_result$inputs
  weight     <- inputs$weight
  mse        <- inputs$mse
  area_names <- inputs$area_names
  target     <- bench_result$benchmark$target
  D          <- length(inputs$estimate)

  posterior_samples <- as.matrix(posterior_samples)
  if (ncol(posterior_samples) != D)
    stop("Number of columns in posterior_samples must match the number ",
         "of areas (", D, ").")

  G <- nrow(posterior_samples)
  bench_samples <- matrix(NA_real_, nrow = G, ncol = D)

  for (g in seq_len(G)) {
    theta_g <- posterior_samples[g, ]
    point_g <- .compute_benchmark_point(method, theta_g, mse, weight, target)
    bench_samples[g, ] <- point_g$SAE
  }

  Estimation_HB <- colMeans(bench_samples, na.rm = TRUE)
  MSE_HB        <- apply(bench_samples, 2, stats::var, na.rm = TRUE)

  list(
    method        = method,
    approach      = "posterior",
    Estimation    = data.frame(Estimation = Estimation_HB, row.names = area_names),
    MSE           = data.frame(MSE = MSE_HB, row.names = area_names),
    G             = G,
    bench_samples = bench_samples
  )
}
