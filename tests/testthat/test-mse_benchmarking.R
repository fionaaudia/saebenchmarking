Z   <- cbind(1, data_eblup$z)
s2v <- attr(data_eblup, "sigma2_v")

bench <- function(method, ...) {
  sae_benchmarking(method = method, direct = "direct", weight = "weight",
                   estimate = "eblup", mse = "mse", vardir = "vardir",
                   data = data_eblup, ...)
}

test_that("difference method returns finite MSEs for every area", {
  res <- mse_benchmarking(bench("difference"), estimator = "eblup", z = Z,
                          sigma2_v = s2v)
  expect_s3_class(res, "mse_benchmarking")
  expect_length(res$mse, nrow(data_eblup))
  expect_true(all(is.finite(res$mse)))
  expect_true(is.na(res$B))
})

test_that("bootstrap is reproducible and leaves the RNG state untouched", {
  skip_on_cran()
  bm <- bench("ratio")
  set.seed(42)
  before <- .Random.seed
  r1 <- mse_benchmarking(bm, estimator = "eblup", z = Z, sigma2_v = s2v,
                         B = 20, seed = 1)
  expect_identical(.Random.seed, before)
  r2 <- mse_benchmarking(bm, estimator = "eblup", z = Z, sigma2_v = s2v,
                         B = 20, seed = 1)
  expect_equal(r1$mse, r2$mse)
  expect_true(all(r1$mse > 0))
})

test_that("HB posterior MSE is at least the posterior variance", {
  for (m in c("difference", "ratio", "optimum")) {
    bm <- sae_benchmarking(method = m, direct = "direct", weight = "weight",
                           estimate = "theta_hb", vardir = "vardir",
                           phi_source = "vardir", data = data_hb)
    res <- mse_benchmarking(bm, estimator = "hb",
                            theta_hb = data_hb$theta_hb,
                            V_hb = data_hb$var_hb)
    expect_true(all(res$mse >= data_hb$var_hb))
    expect_equal(sum(bm$inputs$weight * res$estimate), bm$target)
  }
})

test_that("invalid inputs give informative errors", {
  bm <- bench("difference")
  expect_error(mse_benchmarking(list()), "class 'sae_benchmarking'")
  expect_error(mse_benchmarking(bm, estimator = "eblup", sigma2_v = s2v),
               "'z' and 'sigma2_v'")
  expect_error(mse_benchmarking(bm, estimator = "eblup", z = Z[1:5, ],
                                sigma2_v = s2v), "rows")
  expect_error(mse_benchmarking(bm, estimator = "eblup", z = Z,
                                sigma2_v = -1), "non-negative")
  expect_error(mse_benchmarking(bm, estimator = "hb"), "theta_hb")

  bm_nov <- sae_benchmarking("difference", direct = "direct",
                             weight = "weight", estimate = "eblup",
                             data = data_eblup)
  expect_error(mse_benchmarking(bm_nov, estimator = "eblup", z = Z,
                                sigma2_v = s2v), "'vardir' must be supplied")
})
