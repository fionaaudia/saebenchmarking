est    <- c(10.2, 8.7, 12.1)
direct <- c(10.5, 8.2, 12.8)
mse    <- c(0.5, 0.4, 0.6)
vardir <- c(0.6, 0.5, 0.7)
w      <- c(0.3, 0.4, 0.3)

test_that("every method satisfies the benchmarking constraint", {
  for (m in c("difference", "ratio", "optimum")) {
    res <- sae_benchmarking(method = m, direct = direct, weight = w,
                            estimate = est, mse = mse, vardir = vardir)
    expect_s3_class(res, "sae_benchmarking")
    expect_equal(sum(w * res$estimates$Estimate), sum(w * direct))
    expect_equal(res$target, sum(w * direct))
  }
})

test_that("difference and ratio use a common adjustment", {
  d <- sae_benchmarking("difference", direct = direct, weight = w,
                        estimate = est)
  expect_equal(d$estimates$Estimate - est, rep(d$params$adjustment, 3))
  r <- sae_benchmarking("ratio", direct = direct, weight = w, estimate = est)
  expect_equal(r$estimates$Estimate / est, rep(r$params$ratio, 3))
})

test_that("weights are rescaled to sum to one", {
  a <- sae_benchmarking("ratio", direct = direct, weight = w, estimate = est)
  b <- sae_benchmarking("ratio", direct = direct, weight = w * 100,
                        estimate = est)
  expect_equal(a$estimates, b$estimates)
})

test_that("column names of 'data' can be used", {
  df <- data.frame(d = direct, w = w, e = est, m = mse)
  a <- sae_benchmarking("optimum", direct = "d", weight = "w",
                        estimate = "e", mse = "m", data = df)
  b <- sae_benchmarking("optimum", direct = direct, weight = w,
                        estimate = est, mse = mse)
  expect_equal(a$estimates, b$estimates)
})

test_that("invalid inputs give informative errors", {
  expect_error(sae_benchmarking("ratio", direct = direct, weight = w[-1],
                                estimate = est), "same length")
  expect_error(sae_benchmarking("ratio", direct = direct, weight = -w,
                                estimate = est), "non-negative")
  expect_error(sae_benchmarking("ratio", direct = direct, weight = c(0, 0, 0),
                                estimate = est), "must be positive")
  expect_error(sae_benchmarking("ratio", direct = c(NA, 1, 2), weight = w,
                                estimate = est), "Missing values")
  expect_error(sae_benchmarking("optimum", direct = direct, weight = w,
                                estimate = est), "requires 'mse'")
  expect_error(sae_benchmarking("ratio", direct = "x", weight = w,
                                estimate = est,
                                data = data.frame(y = 1:3)), "not found")
  expect_error(sae_benchmarking("ratio", weight = w, estimate = est),
               "required")
})

test_that("S3 methods work", {
  res <- sae_benchmarking("difference", direct = direct, weight = w,
                          estimate = est, area_names = c("A", "B", "C"))
  expect_named(coef(res), c("A", "B", "C"))
  expect_equal(fitted(res), coef(res))
  expect_s3_class(as.data.frame(res), "data.frame")
  expect_s3_class(summary(res), "summary.sae_benchmarking")
  expect_output(print(res), "Small area benchmarking")
  expect_output(print(summary(res)), "Number of areas")
})
