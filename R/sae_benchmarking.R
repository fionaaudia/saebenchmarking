#' Benchmark small area estimates
#'
#' Applies difference, ratio, or optimum benchmarking to small area
#' estimates so that their weighted aggregate matches a known target,
#' following You & Rao (2002), Wang, Fuller & Qu (2008), and
#' Rao & Molina (2015, Section 6.4.6).
#'
#' @param method Character. One of `"difference"`, `"ratio"`, or
#'   `"optimum"`.
#' @param estimate Numeric vector of indirect small area estimates
#'   (e.g. EBLUP, EB, or HB posterior means), one value per area.
#' @param mse Numeric vector of mean squared errors corresponding to
#'   `estimate`, one value per area.
#' @param weight Numeric vector of benchmarking weights (e.g.
#'   population proportions), one per area. Automatically normalised
#'   to sum to 1 with a warning if it does not.
#' @param target Optional single numeric value used as the external
#'   benchmark target. If `NULL`, `direct` must be supplied.
#' @param direct Optional numeric vector of direct estimates, used to
#'   compute an internal target when `target` is `NULL`.
#' @param vardir Optional numeric vector of sampling variances of
#'   `direct`. Only used later by [mse_benchmarking()] for bootstrap
#'   residuals; unused here.
#' @param area_names Optional character vector of area labels.
#'
#' @return A list with elements `method`, `benchmark` (containing
#'   `estimation`, `params`, `target`, `benchmark_type`),
#'   `agregation`, and `inputs` (raw inputs retained for
#'   [mse_benchmarking()]).
#'
#' @seealso [mse_benchmarking()] for estimating the mean squared
#'   error of the benchmarked estimates returned by this function.
#'
#' @examples
#' est <- c(10.2, 8.7, 12.1)
#' mse_val <- c(0.5, 0.4, 0.6)
#' w <- c(0.3, 0.4, 0.3)
#' sae_benchmarking(method = "difference", estimate = est, mse = mse_val,
#'               weight = w, target = 10)
#'
#' @export
sae_benchmarking <- function(method = c("difference", "ratio", "optimum"),
                          estimate, mse, weight,
                          target = NULL, direct = NULL, vardir = NULL,
                          area_names = NULL) {

  method <- match.arg(method)
  D <- length(estimate)

  if (length(mse) != D)
    stop("Length of mse is not appropriate, the length must be ", D)
  if (length(weight) != D)
    stop("Length of weight is not appropriate, the length must be ", D)
  if (any(is.na(estimate)))
    stop("Object estimate contains NA values.")
  if (any(is.na(mse)))
    stop("Object mse contains NA values.")
  if (any(is.na(weight)))
    stop("Object weight contains NA values.")
  if (any(mse < 0))
    stop("Object mse is not appropriate, all values must be non-negative")
  if (any(weight < 0))
    stop("Object weight is not appropriate, all values must be non-negative")

  if (!is.null(direct) && length(direct) != D)
    stop("Length of direct is not appropriate, the length must be ", D)
  if (!is.null(vardir) && length(vardir) != D)
    stop("Length of vardir is not appropriate, the length must be ", D)
  if (!is.null(vardir) && any(vardir[!is.na(vardir)] <= 0))
    stop("Object vardir is not appropriate, all non-NA values must be strictly positive (> 0)")

  if (is.null(area_names)) {
    area_names <- paste0("Area_", seq_len(D))
  } else {
    if (length(area_names) != D)
      stop("Length of area_names is not appropriate, the length must be ", D)
    area_names <- as.character(area_names)
  }

  wsum <- sum(weight)
  if (abs(wsum - 1) > 1e-6) {
    warning("Sum of weight is not 1 (", round(wsum, 6),
            "); weights are normalised automatically.")
    weight <- weight / wsum
  }

  if (!is.null(target)) {
    if (length(target) != 1 || is.na(target))
      stop("target must be a single, non-NA numeric value.")
    benchmark_type <- "external"
  } else {
    if (is.null(direct))
      stop("Either 'target' (external benchmarking) or a complete 'direct' ",
           "vector (internal benchmarking) must be supplied.")
    if (any(is.na(direct)))
      stop("'direct' contains NA values (e.g. out-of-sample areas), so an ",
           "internal target cannot be computed reliably.")
    target <- sum(weight * direct)
    benchmark_type <- "internal"
  }

  if (method == "ratio" && sum(weight * estimate) == 0)
    stop("sum(weight * estimate) is 0; ratio benchmarking is not defined ",
         "(division by zero). Consider using method = 'difference' instead.")

  point  <- .compute_benchmark_point(method, estimate, mse, weight, target)
  SAE    <- data.frame(Estimation = point$SAE, row.names = area_names)
  params <- point$params

  Aggregation_Bench    <- sum(SAE$Estimation * weight)
  Aggregation_Estimate <- sum(estimate * weight)
  Aggregation_Direct   <- if (!is.null(direct)) sum(direct * weight, na.rm = TRUE) else NA

  Aggregation <- matrix(c(Aggregation_Direct, Aggregation_Bench, Aggregation_Estimate, target),
                        4, 1)
  rownames(Aggregation) <- c("Aggregation_Direct", "Aggregation_Bench",
                             "Aggregation_Estimate", "Target")
  colnames(Aggregation) <- "Y"

  result <- list()
  result$method                   <- method
  result$benchmark$estimation     <- SAE
  result$benchmark$params         <- params
  result$benchmark$target         <- target
  result$benchmark$benchmark_type <- benchmark_type
  result$aggregation               <- Aggregation
  result$inputs <- list(estimate = estimate, mse = mse, weight = weight,
                        direct = direct, vardir = vardir, area_names = area_names)

  return(result)
}
