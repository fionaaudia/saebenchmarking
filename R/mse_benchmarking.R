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
#' @details
#' For **external benchmarking** (`target` supplied directly to
#' [sae_benchmarking()]), the bootstrap MSE is fully model-agnostic:
#' only `estimate` is perturbed across replications, since `target` is
#' treated as a fixed, known constant.
#'
#' For **internal benchmarking** (`target` computed from `direct`),
#' `target` is itself an estimate with its own sampling uncertainty
#' (`vardir`), so it is re-simulated at every bootstrap replication
#' alongside `estimate`. Because the correlation between the sampling
#' errors of `estimate` and `direct` depends on the (unknown, model-
#' specific) shrinkage structure of `estimate`, a default correlation
#' `rho_i = sqrt(mse_i / vardir_i)` is assumed, motivated by the
#' shrinkage relationship of Fay-Herriot-type EBLUP estimators (Rao &
#' Molina, 2015, Ch. 6). **This is an approximation.** If `estimate`
#' comes from a materially different model, supply `rho` explicitly.
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
#' @param rho Optional numeric scalar or vector (length equal to the
#'   number of areas) giving the assumed correlation between the
#'   sampling errors of `estimate` and `direct`. Only relevant for
#'   **internal** benchmarking. If `NULL` (default), it is set to
#'   `sqrt(mse / vardir)` (capped at 1). Ignored for external
#'   benchmarking.
#'
#' @return A list with elements `method`, `approach` (`"bootstrap"` or
#'   `"posterior"`), `estimator_type` (`"EBLUP/EB"` or `"HB"`), `MSE`
#'   (a data frame with one row per area), and additional diagnostic
#'   elements specific to the approach used. The bootstrap approach
#'   additionally reports `n_failed` and, for internal benchmarking,
#'   `rho_used`. The posterior approach additionally reports
#'   `Estimation` (posterior mean of the benchmarked draws) and `G`.
#'
#' @references
#' Butar, F.B. & Lahiri, P. (2003). On measures of uncertainty of
#' empirical Bayes small-area estimators. \emph{Journal of Statistical
#' Planning and Inference}, 112, 63-76.
#'
#' Steorts, R.C. & Ghosh, M. (2013). On estimation of mean squared
#' errors of benchmarked empirical Bayes estimators. \emph{Statistica
#' Sinica}, 23(2), 749-767.
#'
#' Datta, G.S., Ghosh, M., Steorts, R. & Maples, J. (2011). Bayesian
#' benchmarking with applications to small area estimation.
#' \emph{TEST}, 20, 574-588.
#'
#' Bell, W.R., Datta, G.S. & Ghosh, M. (2013). Benchmarking small area
#' estimators. \emph{Biometrika}, 100(1), 189-202.
#'
#' Rao, J.N.K. & Molina, I. (2015). \emph{Small Area Estimation}, 2nd
#' Edition. Wiley, Chapter 6.
#'
#' @seealso [sae_benchmarking()] for computing the benchmarked point
#' estimates required by this function.
#'
#' @examples
#' # --- External benchmarking, bootstrap MSE ---
#' est <- c(10.2, 8.7, 12.1)
#' mse_val <- c(0.5, 0.4, 0.6)
#' w <- c(0.3, 0.4, 0.3)
#' hasil <- sae_benchmarking(method = "difference", estimate = est,
#'                            mse = mse_val, weight = w, target = 10)
#' mse_benchmarking(hasil, B = 200, seed = 1)
#'
#' # --- Internal benchmarking, bootstrap MSE ---
#' direct <- c(10.5, 8.2, 12.8)
#' vardir <- c(0.6, 0.5, 0.7)
#' hasil_internal <- sae_benchmarking(method = "difference", estimate = est,
#'                                     mse = mse_val, weight = w,
#'                                     direct = direct, vardir = vardir)
#' mse_benchmarking(hasil_internal, B = 200, seed = 1)
#'
#' # --- Hierarchical Bayes, posterior simulation ---
#' set.seed(1)
#' post <- matrix(rnorm(300 * 3, rep(est, each = 300), sqrt(mse_val)), ncol = 3)
#' mse_benchmarking(hasil, posterior_samples = post)
#'
#' @export
mse_benchmarking <- function(bench_result, posterior_samples = NULL,
                             B = 1000, seed = NULL, rho = NULL) {

  if (!is.list(bench_result) || is.null(bench_result$method))
    stop("bench_result must be the output of sae_benchmarking().")

  if (!is.null(posterior_samples)) {
    result <- .mse_via_posterior(bench_result, posterior_samples)
  } else {
    if (!is.numeric(B) || length(B) != 1 || B < 1)
      stop("B must be a single positive integer.")
    result <- .mse_via_bootstrap(bench_result, B = B, seed = seed, rho = rho)
  }

  result
}


#' Estimasi MSE hasil benchmarking via bootstrap residual (internal)
#'
#' Berlaku untuk estimator EBLUP maupun EB. Untuk external
#' benchmarking, hanya `estimate` yang disimulasikan ulang di tiap
#' replikasi (target dianggap konstanta pasti). Untuk internal
#' benchmarking, `direct` (dan karenanya `target`) ikut disimulasikan
#' ulang, dengan korelasi `rho` terhadap error `estimate` (lihat
#' `@details` pada [mse_benchmarking()]).
#'
#' PERBAIKAN vs versi sebelumnya: `MSE_bench` dihitung sebagai
#' `Var(boot_estimation)` secara LANGSUNG, TANPA ditambah `mse` lagi.
#' `boot_estimation` sudah berisi estimator PENUH hasil benchmark
#' (bukan hanya komponen koreksinya), sehingga variansnya di seluruh
#' replikasi SUDAH merupakan estimasi bootstrap dari
#' MSE(theta_hat_bench) secara langsung. Bukti kasus batas: jika
#' weight_i = 0 (area tidak menyerap koreksi apa pun), maka
#' Var_boot_i harus persis sama dengan mse_i asli -- menambahkan mse_i
#' lagi (seperti versi lama) akan menghitungnya dua kali.
#'
#' @keywords internal
#' @noRd
.mse_via_bootstrap <- function(bench_result, B = 1000, seed = NULL,
                               rho = NULL) {

  method         <- bench_result$method
  inputs         <- bench_result$inputs
  estimate       <- inputs$estimate
  mse            <- inputs$mse
  weight         <- inputs$weight
  direct         <- inputs$direct
  vardir         <- inputs$vardir
  area_names     <- inputs$area_names
  target         <- bench_result$benchmark$target
  benchmark_type <- bench_result$benchmark$benchmark_type
  D              <- length(estimate)

  if (!is.null(seed)) set.seed(seed)

  is_internal <- identical(benchmark_type, "internal") &&
    !is.null(direct) && !is.null(vardir)

  if (is_internal) {
    if (is.null(rho)) {
      rho <- sqrt(pmin(mse / vardir, 1))
    } else {
      if (length(rho) == 1) rho <- rep(rho, D)
      if (length(rho) != D)
        stop("'rho' harus berupa skalar atau vektor sepanjang jumlah area (", D, ").")
      if (any(rho < -1 | rho > 1))
        stop("'rho' harus berada pada rentang [-1, 1].")
    }
  }

  have_resid_pool <- !is.null(direct) && !is.null(vardir) &&
    any(!is.na(direct) & !is.na(vardir))

  if (have_resid_pool) {
    ok <- !is.na(direct) & !is.na(vardir)
    std_resid <- (direct[ok] - estimate[ok]) / sqrt(vardir[ok])
    std_resid <- std_resid - mean(std_resid)
    draw_z <- function() sample(std_resid, size = D, replace = TRUE)
  } else {
    warning("No usable 'direct'/'vardir' residual pool found; falling back ",
            "to a parametric bootstrap using N(0, 1) pseudo-errors.")
    draw_z <- function() stats::rnorm(D)
  }

  boot_estimation <- matrix(NA_real_, nrow = B, ncol = D)

  for (b in seq_len(B)) {

    z_est     <- draw_z()
    error_est <- z_est * sqrt(mse)
    estimate_star <- estimate + error_est

    if (is_internal) {
      # Bangkitkan noise 'direct' berkorelasi rho dengan noise
      # 'estimate' (per area), lewat dekomposisi 2-komponen standar:
      # error_dir_i = rho_i * z_est_i * sqrt(vardir_i)
      #               + sqrt(1 - rho_i^2) * z_indep_i * sqrt(vardir_i)
      z_indep   <- stats::rnorm(D)
      error_dir <- rho * z_est * sqrt(vardir) +
        sqrt(1 - rho^2) * z_indep * sqrt(vardir)
      target_star <- sum(weight * (direct + error_dir))
    } else {
      target_star <- target
    }

    point_star <- .compute_benchmark_point(method, estimate_star, mse,
                                           weight, target_star)
    boot_estimation[b, ] <- point_star$SAE
  }

  # --- MSE_bench = Var_boot langsung (TIDAK ditambah 'mse' lagi) ---
  MSE_bench <- apply(boot_estimation, 2, stats::var, na.rm = TRUE)

  list(
    method         = method,
    approach       = "bootstrap",
    estimator_type = "EBLUP/EB",
    MSE            = data.frame(MSE = MSE_bench, row.names = area_names),
    B              = B,
    n_failed       = sum(is.na(boot_estimation[, 1])),
    rho_used       = if (is_internal) rho else NA
  )
}


#' Estimasi titik dan MSE hasil benchmarking dari sampel posterior (internal)
#'
#' Berlaku untuk estimator Hierarchical Bayes. Varians posterior
#' setelah benchmarking dihitung langsung dari sebaran posterior
#' penuh, sehingga secara definisi merupakan risiko Bayes (MSE) di
#' bawah squared loss -- tidak terpengaruh isu double-counting pada
#' jalur bootstrap.
#'
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

  if (method == "optimum") {
    post_var <- apply(posterior_samples, 2, stats::var)
    rel_diff <- abs(post_var - mse) / pmax(mse, .Machine$double.eps)
    if (any(rel_diff > 0.5, na.rm = TRUE)) {
      warning(
        "Untuk method = 'optimum', bobot lambda_i dihitung dari 'mse' ",
        "yang disuplai ke sae_benchmarking(), BUKAN dari varians ",
        "posterior_samples. Ditemukan selisih relatif > 50% pada ",
        "sejumlah area antara 'mse' dan var(posterior_samples). ",
        "Pastikan 'mse' berasal dari sumber yang konsisten dengan ",
        "posterior_samples.",
        call. = FALSE
      )
    }
  }

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
    method         = method,
    approach       = "posterior",
    estimator_type = "HB",
    Estimation     = data.frame(Estimation = Estimation_HB, row.names = area_names),
    MSE            = data.frame(MSE = MSE_HB, row.names = area_names),
    G              = G,
    bench_samples  = bench_samples
  )
}
