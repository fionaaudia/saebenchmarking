test_that("MSE bootstrap selalu >= mse asli (benchmarking menambah ketidakpastian)", {
  est <- c(10, 8, 12); mse_val <- c(0.5, 0.4, 0.6); w <- c(0.3, 0.4, 0.3)
  direct <- c(10.5, 8.2, 12.8); vardir <- c(0.6, 0.5, 0.7)

  hasil <- sae_benchmarking(method = "difference", estimate = est, mse = mse_val,
                         weight = w, direct = direct, vardir = vardir)
  mse_hasil <- mse_benchmarking(hasil, B = 500, seed = 1)

  expect_true(all(mse_hasil$MSE$MSE >= mse_val - 1e-8))
})

test_that("jalur posterior menghasilkan MSE dan estimation dengan panjang sesuai", {
  est <- c(10, 8, 12); mse_val <- c(0.5, 0.4, 0.6); w <- c(0.3, 0.4, 0.3)
  hasil <- sae_benchmarking(method = "difference", estimate = est,
                         mse = mse_val, weight = w, target = 10)

  set.seed(1)
  post <- matrix(rnorm(300 * 3, rep(est, each = 300), sqrt(mse_val)), ncol = 3)
  mse_hb <- mse_benchmarking(hasil, posterior_samples = post)

  expect_equal(nrow(mse_hb$MSE), 3)
  expect_equal(mse_hb$approach, "posterior")
})

test_that("bench_result yang tidak valid memicu error", {
  expect_error(mse_benchmarking(list(a = 1)))
})
