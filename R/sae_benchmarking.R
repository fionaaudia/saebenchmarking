#' @title Benchmark small area estimates
#'
#' @description Applies difference, ratio, or optimum benchmarking to small area estimates
#' so that their weighted aggregate equals the weighted aggregate of the
#' direct estimates.
#'
#' @param method Character. One of `"ratio"`, `"difference"`, or `"optimum"`.
#' @param direct Numeric vector of direct estimates, or the name of a column
#'   of `data`. Used to compute the target `sum(weight * direct)`.
#' @param weight Numeric vector of non-negative benchmarking weights (e.g.
#'   population proportions), or the name of a column of `data`. Rescaled to
#'   sum to one.
#' @param estimate Numeric vector of model-based small area estimates (e.g.
#'   EBLUP or HB), or the name of a column of `data`.
#' @param mse Optional numeric vector of mean squared errors of `estimate`,
#'   or the name of a column of `data`. Required when `method = "optimum"`
#'   and `phi_source = "mse"`.
#' @param vardir Optional numeric vector of sampling variances of `direct`,
#'   or the name of a column of `data`. Required when `method = "optimum"`
#'   and `phi_source = "vardir"`, and by [mse_benchmarking()] with
#'   `estimator = "eblup"`.
#' @param phi_source Character. One of `"mse"` or `"vardir"`; the quantity
#'   used to build the optimum benchmarking factor. Ignored for other methods.
#' @param area_names Optional character vector of area labels. Defaults to
#'   `"Area_1"`, `"Area_2"`, and so on.
#' @param data Optional data frame. When supplied, `direct`, `weight`,
#'   `estimate`, `mse`, and `vardir` may be given as column names.
#'
#' @return This function returns a list of class `"sae_benchmarking"` with elements:
#' * `method`: the benchmarking method used.
#' * `target`: numeric value of the target.
#' * `estimates`: a one-column data frame (`Estimate`) of benchmarked
#'   estimates, with area names as row names.
#' * `params`: a list with the adjustment parameter (`adjustment`, `ratio`,
#'   or `lambda`, depending on `method`).
#' * `aggregation`: a one-row matrix with columns `Direct`, `Estimate`,
#'   `Bench`, and `Target` giving the weighted aggregates.
#' * `inputs`: a list of the inputs used, retained for [mse_benchmarking()].
#'
#' @seealso [mse_benchmarking()] to estimate the mean squared error of the
#'   benchmarked estimates.
#'
#' @examples
#' est    <- c(10.2, 8.7, 12.1)
#' direct <- c(10.5, 8.2, 12.8)
#' mse    <- c(0.5, 0.4, 0.6)
#' w      <- c(0.3, 0.4, 0.3)
#'
#' res <- sae_benchmarking(method = "difference", direct = direct,
#'                         weight = w, estimate = est, mse = mse)
#' res
#' coef(res)
#' summary(res)
#'
#' @export
sae_benchmarking <- function(method = c("ratio", "difference", "optimum"),
                             direct, weight, estimate,
                             mse = NULL, vardir = NULL,
                             phi_source = c("mse", "vardir"),
                             area_names = NULL, data = NULL) {

  method     <- match.arg(method)
  phi_source <- match.arg(phi_source)

  if (missing(direct) || missing(weight) || missing(estimate)) {
    stop("Arguments 'direct', 'weight', and 'estimate' are required.",
         call. = FALSE)
  }

  if (!is.null(data) && !is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  get_col <- function(x, arg) {
    if (is.null(x)) return(NULL)
    if (is.character(x) && length(x) == 1L) {
      if (is.null(data)) {
        stop(sprintf("'%s' is a column name but 'data' was not supplied.",
                     arg), call. = FALSE)
      }
      if (!x %in% names(data)) {
        stop(sprintf("Column '%s' not found in 'data'.", x), call. = FALSE)
      }
      return(as.numeric(data[[x]]))
    }
    if (is.matrix(x) || is.data.frame(x)) x <- as.numeric(unlist(x))
    x
  }

  direct   <- get_col(direct, "direct")
  weight   <- get_col(weight, "weight")
  estimate <- get_col(estimate, "estimate")
  mse      <- get_col(mse, "mse")
  vardir   <- get_col(vardir, "vardir")

  D <- length(estimate)

  area_names <- if (is.null(area_names)) paste0("Area_", seq_len(D)) else as.character(area_names)
  if (length(area_names) != D) {
    stop("'area_names' must have length ", D, ".", call. = FALSE)
  }

  phi_inv <- .validate_sae_inputs(method = method, direct = direct, weight = weight,
                                  estimate = estimate, mse = mse, vardir = vardir,
                                  phi_source = phi_source, D = D)

  if (sum(weight) <= 0) {
    stop("The sum of 'weight' must be positive.", call. = FALSE)
  }

  weight <- weight / sum(weight)

  target <- sum(weight * direct)

  aggregate_est <- sum(weight * estimate)
  gap           <- target - aggregate_est

  if (method == "ratio") {
    if (aggregate_est == 0) {
      stop("Weighted sum of 'estimate' is zero; ratio method is undefined.",
           call. = FALSE)
    }
    ratio_factor <- target / aggregate_est
    bench  <- estimate * ratio_factor
    params <- list(ratio = ratio_factor)

  } else if (method == "difference") {
    bench  <- estimate + gap
    params <- list(adjustment = gap)

  } else {
    denom  <- sum(phi_inv * weight^2)
    lambda <- (phi_inv * weight) / denom
    bench  <- estimate + lambda * gap
    params <- list(lambda = lambda)
  }

  if (is.null(mse)) mse <- rep(NA_real_, D)

  aggregation <- matrix(
    c(target, aggregate_est, sum(weight * bench), target),
    nrow = 1,
    dimnames = list("Weighted aggregate",
                    c("Direct", "Estimate", "Bench", "Target"))
  )

  result <- list(
    method      = method,
    target      = target,
    estimates   = data.frame(Estimate = bench, row.names = area_names),
    params      = params,
    aggregation = aggregation,
    inputs      = list(estimate = estimate, mse = mse, weight = weight,
                       direct = direct, vardir = vardir,
                       phi_source = if (method == "optimum") phi_source else NA,
                       area_names = area_names)
  )
  class(result) <- "sae_benchmarking"
  result
}

#' Validate the raw inputs to sae_benchmarking()
#'
#' sae_benchmarking()'s own, self-contained validation function -- it does
#' not call any shared/external helper. A small local closure
#' (`check_size_nonneg`) avoids writing out the same length + non-negativity
#' check twice for `mse` and `vardir`.
#'
#' @noRd
.validate_sae_inputs <- function(method, direct, weight, estimate, mse, vardir,
                                 phi_source, D) {

  if (!is.numeric(estimate)) stop("'estimate' must be numeric.", call. = FALSE)
  if (!is.numeric(weight))   stop("'weight' must be numeric.", call. = FALSE)
  if (!is.numeric(direct))   stop("'direct' must be numeric.", call. = FALSE)

  if (length(weight) != D || length(direct) != D) {
    stop("'weight' and 'direct' must have the same length as 'estimate' (",
         D, ").", call. = FALSE)
  }
  if (any(weight < 0, na.rm = TRUE)) {
    stop("'weight' must be non-negative.", call. = FALSE)
  }

  check_size_nonneg <- function(x, arg) {
    if (is.null(x)) return(invisible(NULL))
    if (length(x) != D) {
      stop("'", arg, "' must have length ", D, ".", call. = FALSE)
    }
    if (any(x < 0, na.rm = TRUE)) {
      stop("'", arg, "' must be non-negative.", call. = FALSE)
    }
    invisible(NULL)
  }
  check_size_nonneg(mse, "mse")
  check_size_nonneg(vardir, "vardir")

  phi_inv <- NULL
  if (method == "optimum") {
    phi_inv <- if (phi_source == "mse") mse else vardir
    if (is.null(phi_inv)) {
      stop("method = 'optimum' with phi_source = '", phi_source,
           "' requires '", phi_source, "' to be supplied.", call. = FALSE)
    }
    if (any(phi_inv <= 0, na.rm = TRUE)) {
      stop("'", phi_source, "' must be strictly positive for method = 'optimum'.",
           call. = FALSE)
    }
  }

  na_check <- list(estimate = estimate, weight = weight, direct = direct)
  if (!is.null(phi_inv)) na_check[[phi_source]] <- phi_inv
  has_na <- vapply(na_check, anyNA, logical(1))
  if (any(has_na)) {
    stop("Missing values found in: ",
         paste(names(na_check)[has_na], collapse = ", "), ".", call. = FALSE)
  }

  invisible(phi_inv)
}

#' @export
print.sae_benchmarking <- function(x, ...) {
  cat("Small area benchmarking (method = '", x$method, "')\n", sep = "")
  cat("Target:", format(x$target, digits = 6), "\n\n")
  print(x$estimates)
  invisible(x)
}

#' @export
summary.sae_benchmarking <- function(object, ...) {
  out <- list(method = object$method,
              target = object$target,
              aggregation = object$aggregation,
              n_areas = nrow(object$estimates))
  class(out) <- "summary.sae_benchmarking"
  out
}

#' @export
print.summary.sae_benchmarking <- function(x, ...) {
  cat("Small area benchmarking summary\n")
  cat("Method          :", x$method, "\n")
  cat("Number of areas :", x$n_areas, "\n\n")
  print(x$aggregation)
  invisible(x)
}

#' @export
coef.sae_benchmarking <- function(object, ...) {
  stats::setNames(object$estimates$Estimate, rownames(object$estimates))
}

#' @export
fitted.sae_benchmarking <- function(object, ...)
  coef.sae_benchmarking(object)

#' @export
as.data.frame.sae_benchmarking <- function(x, row.names = NULL,
                                           optional = FALSE, ...) {
  x$estimates
}
