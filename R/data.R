#' @title Sample Data for Fay-Herriot Model with EBLUP Estimation
#'
#' @description Dataset to simulate benchmarking of small area estimates
#'   obtained by the empirical best linear unbiased predictor (EBLUP) under
#'   the Fay-Herriot model.
#'
#' This data is generated based on the Fay-Herriot model by these following
#' steps:
#' \enumerate{
#'   \item Take the population size \code{pop} (in thousands) and the sample
#'     size \code{n} of 50 areas from Table 1 of Wang, Fuller and Qu (2008).
#'   \item Calculate the known sampling variance
#'     \eqn{\psi_i = 16 / n_i}{psi_i = 16 / n_i} and the benchmarking weight
#'     \eqn{w_i = pop_i / \sum pop_i}{w_i = pop_i / sum(pop)}.
#'   \item Generate the auxiliary variable \code{z ~ N(1, 1)} once and keep
#'     it fixed.
#'   \item Generate the random effect \code{v ~ N(0, 1)}. Set
#'     \eqn{\beta_0}{beta0} = 6 and \eqn{\beta_1}{beta1} = 3, then calculate
#'     the true area mean
#'     \eqn{\theta_i = \beta_0 + \beta_1 z_i + v_i}{
#'     theta_i = beta0 + beta1 z_i + v_i}.
#'   \item For each area, generate \eqn{n_i}{n_i} unit observations
#'     \eqn{\theta_i + \varepsilon_{ij}}{theta_i + eps_ij} with
#'     \eqn{\varepsilon_{ij} \sim N(0, 16)}{eps_ij ~ N(0, 16)}. Calculate
#'     the direct estimation \code{direct} as their mean and the sampling
#'     error as \code{e = direct - theta}.
#'   \item Estimate the random effect variance by restricted maximum
#'     likelihood (REML), then calculate the EBLUP \code{eblup} and its mean
#'     squared error \code{mse} using [sae::mseFH()].
#'   \item Combine all variables in a data frame named \code{data_eblup}.
#'     The REML estimate of the random effect variance is stored as the
#'     attribute \code{sigma2_v}.
#' }
#'
#' @docType data
#'
#' @usage data(data_eblup)
#'
#' @format A data frame with 50 rows and 12 variables:
#' \describe{
#'   \item{area}{Area identifier, 1 to 50}
#'   \item{pop}{Population size of each area (in thousands)}
#'   \item{n}{Sample size of each area}
#'   \item{weight}{Benchmarking weight, \code{pop / sum(pop)}}
#'   \item{z}{Auxiliary variable}
#'   \item{v}{True random effect}
#'   \item{e}{True sampling error}
#'   \item{theta}{True small area mean}
#'   \item{direct}{Direct estimation}
#'   \item{vardir}{Known sampling variance, \code{16 / n}}
#'   \item{eblup}{EBLUP estimation}
#'   \item{mse}{Mean squared error of the EBLUP estimation}
#' }
#'
#' @details The attribute \code{sigma2_v} is dropped when the data frame is
#'   subset. Read it first with \code{attr(data_eblup, "sigma2_v")}.
#'
#' @references
#' Molina, I. and Marhuenda, Y. (2015). sae: An R package for small area
#' estimation. \emph{The R Journal}, 7(1), 81-98.
#'
#' Wang, J., Fuller, W. A. and Qu, Y. (2008). Small area estimation under a
#' restriction. \emph{Survey Methodology}, 34(1), 29-36.
#'
#' You, Y., Rao, J. N. K. and Hidiroglou, M. (2013). On the performance of
#' self benchmarked small area estimators under the Fay-Herriot area level
#' model. \emph{Survey Methodology}, 39(1), 217-229.
#'
#' @keywords datasets
"data_eblup"


#' @title Sample Data for Fay-Herriot Model with Hierarchical Bayes Estimation
#'
#' @description Dataset to simulate benchmarking of small area estimates
#'   obtained by the hierarchical Bayes (HB) method under the Fay-Herriot
#'   model.
#'
#' This data is generated based on the Fay-Herriot model by these following
#' steps:
#' \enumerate{
#'   \item Generate \code{pop}, \code{n}, \code{weight}, \code{z}, \code{v},
#'     \code{e}, \code{theta}, \code{direct}, and \code{vardir} by steps 1
#'     to 5 of \code{\link{data_eblup}}. Both datasets share the same values
#'     of these variables.
#'   \item Set the uniform prior
#'     \eqn{\pi(\beta, \tau^2) \propto 1}{pi(beta, tau2) proportional to 1}
#'     on the regression coefficients and the random effect variance
#'     \eqn{\tau^2}{tau2}, as in Sugasawa, Tamae and Kubokawa (2017), and
#'     treat the sampling variance \code{vardir} as known.
#'   \item Run the Gibbs sampler by updating in turn
#'     \eqn{\tau^2}{tau2} from the inverse gamma distribution,
#'     \eqn{\beta}{beta} from the bivariate normal distribution, and
#'     \eqn{\theta_i}{theta_i} from the normal distribution, using their
#'     full conditional distributions.
#'   \item Discard the first 1,000 iterations as burn-in and keep the next
#'     5,000 draws.
#'   \item Calculate the HB estimation \code{theta_hb} as the mean of the
#'     kept draws and the posterior variance \code{var_hb} as their variance.
#'   \item Combine all variables in a data frame named \code{data_hb}.
#' }
#'
#' @docType data
#'
#' @usage data(data_hb)
#'
#' @format A data frame with 50 rows and 12 variables:
#' \describe{
#'   \item{area}{Area identifier, 1 to 50}
#'   \item{pop}{Population size of each area (in thousands)}
#'   \item{n}{Sample size of each area}
#'   \item{weight}{Benchmarking weight, \code{pop / sum(pop)}}
#'   \item{z}{Auxiliary variable}
#'   \item{v}{True random effect}
#'   \item{e}{True sampling error}
#'   \item{theta}{True small area mean}
#'   \item{direct}{Direct estimation}
#'   \item{vardir}{Known sampling variance, \code{16 / n}}
#'   \item{theta_hb}{HB estimation (posterior mean)}
#'   \item{var_hb}{Posterior variance of the HB estimation}
#' }
#'
#' @details Convergence was checked by comparing three chains with different
#'   starting values and update orders.
#'
#' @references
#' Sugasawa, S., Tamae, H. and Kubokawa, T. (2017). Bayesian estimators for
#' small area models shrinking both means and variances. \emph{Scandinavian
#' Journal of Statistics}, 44(1), 150-167. \doi{10.1111/sjos.12246}
#'
#' Wang, J., Fuller, W. A. and Qu, Y. (2008). Small area estimation under a
#' restriction. \emph{Survey Methodology}, 34(1), 29-36.
#'
#' You, Y., Rao, J. N. K. and Hidiroglou, M. (2013). On the performance of
#' self benchmarked small area estimators under the Fay-Herriot area level
#' model. \emph{Survey Methodology}, 39(1), 217-229.
#'
#' @keywords datasets
"data_hb"
