
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saebenchmarking

<!-- badges: start -->

[![R-CMD-check](https://github.com/fionaaudia/saebenchmarking/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fionaaudia/saebenchmarking/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

saebenchmarking provides **difference**, **ratio**, and **optimum**
benchmarking for small area estimates, so that their weighted aggregate
agrees with the weighted aggregate of the direct estimates, together
with mean squared error (MSE) estimation of the benchmarked estimates.

- Point benchmarking follows Wang, Fuller and Qu (2008) and Rao and
  Molina (2015).
- The MSE of benchmarked EBLUPs under the Fay-Herriot model uses the
  second-order approximation (difference method) or a parametric
  bootstrap (ratio and optimum methods), following Steorts and Ghosh
  (2013).
- The posterior MSE of benchmarked hierarchical Bayes (HB) estimates
  follows Datta, Ghosh, Steorts and Maples (2011).

## Installation

You can install the development version of saebenchmarking from
[GitHub](https://github.com/fionaaudia/saebenchmarking) with:

``` r
# install.packages("pak")
pak::pak("fionaaudia/saebenchmarking")
```

## Example

The package ships with `data_eblup`, a simulated Fay-Herriot dataset
with EBLUP estimates for 50 areas.

``` r
library(saebenchmarking)

bm <- sae_benchmarking(
  method   = "difference",
  direct   = "direct",
  weight   = "weight",
  estimate = "eblup",
  mse      = "mse",
  vardir   = "vardir",
  data     = data_eblup
)

summary(bm)
#> Small area benchmarking summary
#> Method          : difference 
#> Number of areas : 50 
#> 
#>                      Direct Estimate    Bench   Target
#> Weighted aggregate 9.307983  9.33162 9.307983 9.307983
head(coef(bm))
#>    Area_1    Area_2    Area_3    Area_4    Area_5    Area_6 
#>  6.964284  8.636035 13.415885 10.476780  9.714784 15.000472
```

### Ratio and optimum benchmarking

``` r
bm_ratio <- sae_benchmarking("ratio", direct = "direct", weight = "weight",
                             estimate = "eblup", vardir = "vardir",
                             data = data_eblup)

bm_opt <- sae_benchmarking("optimum", direct = "direct", weight = "weight",
                           estimate = "eblup", mse = "mse",
                           vardir = "vardir", phi_source = "mse",
                           data = data_eblup)

bm_opt$aggregation
#>                      Direct Estimate    Bench   Target
#> Weighted aggregate 9.307983  9.33162 9.307983 9.307983
```

### MSE of the benchmarked estimates

``` r
Z   <- cbind(1, data_eblup$z)
s2v <- attr(data_eblup, "sigma2_v")

mse_diff <- mse_benchmarking(bm, estimator = "eblup", z = Z,
                             sigma2_v = s2v, fitting_method = "REML")
head(mse_diff$mse)
#>    Area_1    Area_2    Area_3    Area_4    Area_5    Area_6 
#> 0.2157345 0.2536077 0.2714558 0.2781556 0.3027576 0.3194359

mse_ratio <- mse_benchmarking(bm_ratio, estimator = "eblup", z = Z,
                              sigma2_v = s2v, B = 100, seed = 1)
head(mse_ratio$mse)
#>    Area_1    Area_2    Area_3    Area_4    Area_5    Area_6 
#> 0.2569607 0.2312479 0.2697205 0.3118135 0.2699223 0.3808344
```

For hierarchical Bayes estimates, supply the posterior means and
variances:

``` r
bm_hb <- sae_benchmarking("ratio", direct = "direct", weight = "weight",
                          estimate = "theta_hb", data = data_hb)

mse_hb <- mse_benchmarking(bm_hb, estimator = "hb",
                           theta_hb = data_hb$theta_hb,
                           V_hb = data_hb$var_hb)
head(mse_hb$mse)
#>    Area_1    Area_2    Area_3    Area_4    Area_5    Area_6 
#> 0.2145888 0.2378206 0.2616146 0.3014484 0.2893667 0.3220340
```

## Learn more

See `vignette("saebenchmarking", package = "saebenchmarking")`.
