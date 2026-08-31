#' prescco: Semiparametric Prediction Intervals Under a Censored Covariate
#'
#' Constructs prediction intervals for an outcome \code{Y} when a
#' time-to-event covariate \code{X} is subject to right-censoring, so that only
#' \code{W = min(X, C)} and \code{delta = I(X <= C)} are observed. Fully
#' observed covariates \code{Z} (continuous and discrete) are optional.
#'
#' The package provides:
#' \itemize{
#'   \item maximum-likelihood estimation of the working \code{X | Z} and
#'     \code{C | Z} models (\code{\link{find_alpha1_MLE}},
#'     \code{\link{find_alpha2_MLE}});
#'   \item two coefficient estimators: complete-case least squares
#'     (\code{\link{find_beta_cc}}) and the semiparametrically efficient SPARCC
#'     estimator (\code{\link{find_beta_sparcc}});
#'   \item two families of prediction intervals: PRESCCO
#'     (\code{\link{PRESCCO_prediction_interval}}) and conformal --- split, full, and
#'     jackknife+ (\code{\link{split_conformal_prediction_interval}},
#'     \code{\link{full_conformal_prediction_interval}}, \code{\link{jackknife_plus_prediction_interval}}).
#' }
#'
#' A separate, self-contained, no-covariate simulation study reproducing the
#' paper's simulation studies ships under \code{sim/}. It is run by sourcing
#' \code{sim/run.R} and is not part of the installed namespace.
#'
#' @section Name:
#' \strong{PRE}diction with \strong{S}emiparametric Efficiency under a
#' right-\strong{C}ensored \strong{CO}variate.
#'
#' @section Imports:
#' The tags below are the single source of truth for the package NAMESPACE.
#' Keep them in sync when adding dependencies, then run
#' \code{devtools::document()}.
#'
#' @importFrom stats IQR approx dnorm lm.fit optim optimize pnorm qnorm
#'   quantile rnorm sd setNames var
#' @importFrom utils write.csv
#' @importFrom MASS ginv
#' @importFrom nleqslv nleqslv
#' @importFrom truncnorm dtruncnorm ptruncnorm qtruncnorm rtruncnorm
#'
#' @keywords internal
"_PACKAGE"
