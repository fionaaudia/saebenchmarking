test_that("difference benchmarking hasilkan agregat sama dengan target", {
  est <- c(10, 8, 12); mse_val <- c(0.5, 0.4, 0.6); w <- c(0.3, 0.4, 0.3)
  hasil <- sae_benchmarking(method = "difference", estimate = est,
                         mse = mse_val, weight = w, target = 10)
  agg <- sum(hasil$benchmark$estimation$Estimation * w)
  expect_equal(agg, 10, tolerance = 1e-8)
})

test_that("ratio benchmarking hasilkan agregat sama dengan target", {
  est <- c(10, 8, 12); mse_val <- c(0.5, 0.4, 0.6); w <- c(0.3, 0.4, 0.3)
  hasil <- sae_benchmarking(method = "ratio", estimate = est,
                         mse = mse_val, weight = w, target = 10)
  agg <- sum(hasil$benchmark$estimation$Estimation * w)
  expect_equal(agg, 10, tolerance = 1e-8)
})

test_that("optimum benchmarking hasilkan agregat sama dengan target", {
  est <- c(10, 8, 12); mse_val <- c(0.5, 0.4, 0.6); w <- c(0.3, 0.4, 0.3)
  hasil <- sae_benchmarking(method = "optimum", estimate = est,
                         mse = mse_val, weight = w, target = 10)
  agg <- sum(hasil$benchmark$estimation$Estimation * w)
  expect_equal(agg, 10, tolerance = 1e-8)
})

test_that("target dan direct NULL sekaligus memicu error", {
  expect_error(
    sae_benchmarking(method = "difference", estimate = c(1, 2),
                  mse = c(0.1, 0.1), weight = c(0.5, 0.5))
  )
})

test_that("ratio benchmarking dengan denominator nol memicu error", {
  expect_error(
    sae_benchmarking(method = "ratio", estimate = c(1, -1),
                  mse = c(0.1, 0.1), weight = c(0.5, 0.5), target = 5)
  )
})
