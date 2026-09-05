#' Hitung titik estimasi hasil benchmarking (internal, tidak diekspor)
#' @keywords internal
#' @noRd
.compute_benchmark_point <- function(method, estimate, mse, weight, target) {

  if (method == "difference") {
    shift  <- target - sum(weight * estimate)
    SAE    <- estimate + shift
    params <- list(shift = shift)

  } else if (method == "ratio") {
    denom <- sum(weight * estimate)
    if (denom == 0) {
      SAE    <- rep(NA_real_, length(estimate))
      params <- list(phi = NA_real_)
    } else {
      phi    <- target / denom
      SAE    <- estimate * phi
      params <- list(phi = phi)
    }

  } else if (method == "optimum") {
    phi_inv <- mse
    lambda  <- (phi_inv * weight) / sum(phi_inv * weight^2)
    alpha   <- target - sum(weight * estimate)
    SAE     <- estimate + lambda * alpha
    params  <- list(lambda = lambda, alpha = alpha)
  }

  list(SAE = SAE, params = params)
}
