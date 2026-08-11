# =============================================================================
#  Prediction helpers: covariate basis, interval centers, and residuals
#  ---------------------------------------------------------------------------
#  General (data-in) building blocks used by both the semiparametric and
#  conformal prediction methods in this package. All functions accept optional
#  fully-observed covariates: a continuous block `z_c` and a discrete block
#  `z_d`. When both are NULL the model reduces to the no-covariate case
#  Y | X ~ N((1, X) beta, sigma^2), X ~ TN(mu_X; [w_min, w_max]).
#
#  The outcome model is linear in the basis
#     phi_xz(x, z_c, z_d) = (1, x, z_c, z_d, x * z_c, x * z_d),
#  so `beta` has length 2 in the no-covariate case and grows with the number
#  of covariates and their interactions with X.
#
#  Package imports are declared centrally in R/censcovpred-package.R.
# =============================================================================

#' Covariate basis for the outcome model
#'
#' @param x Numeric event-time covariate value (scalar).
#' @param z_c Optional continuous covariate(s): a scalar/vector for one point.
#' @param z_d Optional discrete covariate(s): a scalar/vector for one point.
#' @return The numeric basis vector (1, x, z_c, z_d, x*z_c, x*z_d).
#' @examples
#' # No covariates: the basis is (1, x)
#' phi_xz(0.5)
#'
#' # One continuous and one discrete covariate:
#' # (1, x, z_c, z_d, x * z_c, x * z_d)
#' phi_xz(0.5, z_c = 2, z_d = 1)
#'
#' # The basis fixes the layout of beta, so a fitted coefficient vector
#' # is read in this order.
#' length(phi_xz(0.5, z_c = c(2, 3), z_d = 1))
#' @export
phi_xz <- function(x, z_c = NULL, z_d = NULL) {
  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)
  
  c(
    1,
    x,
    z_c_vec,
    z_d_vec,
    x * z_c_vec,
    x * z_d_vec
  )
}
#' Interval center m0: the outcome mean phi_xz(x, z) %*% beta
#'
#' @param x Event-time covariate value.
#' @param beta Outcome-model coefficient vector matching phi_xz.
#' @param z_c,z_d Optional covariates.
#' @return The scalar mean outcome at (x, z).
#' @noRd
m0 = function(x, beta, z_c = NULL, z_d = NULL) {
  sum(phi_xz(x, z_c, z_d) * beta)
}
#' Interval center m1: mean outcome given the observed (W, Delta)
#'
#' For a censored observation (delta == 0) with w below the upper bound, the
#' covariate X is integrated over its truncated-normal tail (w, w_max] under
#' the working model N(mu_X(Z), tau1^2); otherwise m1 reduces to m0 at W.
#' @param w,delta Observed time and event indicator (scalars).
#' @param beta Outcome coefficient vector.
#' @param alpha1_star Coefficient vector for the eta1 mean mu_X(Z).
#' @param tau1 Standard deviation of X.
#' @param z_c,z_d Optional covariates.
#' @param w_min,w_max Truncation bounds for X (default -1, 1).
#' @return The scalar center m1.
#' @noRd
m1 = function(w, delta, beta,
              alpha1_star, tau1,
              z_c = NULL, z_d = NULL,
              w_min, w_max) {
  
  ## alpha1_star is the coefficient vector for mu_X(Z):
  ## mu_X(Z) = (1, z_c, z_d)^T alpha1_star
  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)
  
  mu_x <- sum(c(1, z_c_vec, z_d_vec) * alpha1_star)
  v_x  <- tau1^2
  
  if (delta == 0) {
    ## Censored case: need E[X | X > w, X in [w_min, w_max], Z]
    if (w >= w_max) {
      ## Tail essentially empty: any value is fine, use m0 at W
      return(m0(w, beta, z_c_vec, z_d_vec))
    } else {
      lower <- max(w, w_min)
      a <- (lower  - mu_x) / tau1
      b <- (w_max - mu_x) / tau1
      
      denom <- pnorm(b) - pnorm(a)
      if (denom <= 0) {
        ## Fallback if numerically degenerate
        return(m0(lower, beta, z_c_vec, z_d_vec))
      }
      
      ## Truncated normal mean: E[X | a < (X-mu)/tau1 < b]
      x_exp <- mu_x + tau1 * (dnorm(a) - dnorm(b)) / denom
      
      return(m0(x_exp, beta, z_c_vec, z_d_vec))
    }
  } else {
    ## Uncensored: X = W
    return(m0(w, beta, z_c_vec, z_d_vec))
  }
}
#' Residual r1 = |Y - m1|
#'
#' @param y Observed outcome (scalar).
#' @param w,delta Observed time and event indicator (scalars).
#' @param beta Outcome coefficient vector matching \code{phi_xz}.
#' @param alpha1_star Coefficient vector for the eta1 mean mu_X(Z). Supplying a
#'   deliberately wrong value here is what produces the \code{r1*} residual.
#' @param tau1 Standard deviation of X.
#' @param z_c,z_d Optional covariates.
#' @param w_min,w_max Truncation bounds for X.
#' @return The absolute residual |y - m1(...)|.
#' @examples
#' set.seed(1)
#' n  <- 200
#' x  <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' cc <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w  <- pmin(x, cc)
#' delta <- as.integer(x <= cc)
#' y  <- 3 * x + rnorm(n)          # beta = c(0, 3)
#'
#' beta <- c(0, 3)
#'
#' # Uncensored observation: r1 reduces to |y - m0(w)|
#' r1(y[delta == 1][1], w[delta == 1][1], 1L, beta,
#'    alpha1_star = 0, tau1 = 1, w_min = -1, w_max = 1)
#'
#' # Censored observation: X is integrated over the range above w
#' i <- which(delta == 0)[1]
#' r1(y[i], w[i], 0L, beta, alpha1_star = 0, tau1 = 1, w_min = -1, w_max = 1)
#'
#' # r1* is the same call under a working model for X | Z that may be wrong
#' r1(y[i], w[i], 0L, beta, alpha1_star = -2, tau1 = 1, w_min = -1, w_max = 1)
#' @export
r1 = function(y, w, delta, beta,
              alpha1_star, tau1,
              z_c = NULL, z_d = NULL,
              w_min, w_max) {
  
  abs(y - m1(w, delta, beta,
             alpha1_star, tau1,
             z_c, z_d,
             w_min, w_max))
}
#' Residual r2 = |Y - m0(W)|
#'
#' @param y,w,delta Observed outcome, time, and event indicator.
#' @param beta Outcome coefficient vector.
#' @param alpha1_star,tau1 Present for a common signature with r1; unused.
#' @param z_c,z_d Optional covariates.
#' @param w_min,w_max Present for a common signature with r1; unused.
#' @return The absolute residual |y - m0(w)|.
#' @examples
#' set.seed(1)
#' n  <- 200
#' x  <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' cc <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w  <- pmin(x, cc)
#' delta <- as.integer(x <= cc)
#' y  <- 3 * x + rnorm(n)          # beta = c(0, 3)
#'
#' beta <- c(0, 3)
#' i <- which(delta == 0)[1]
#'
#' # r2 plugs w in for x even when censored, so it never touches the
#' # X | Z model -- unlike r1, which integrates over X when delta == 0.
#' r2(y[i], w[i], 0L, beta, alpha1_star = 0, tau1 = 1, w_min = -1, w_max = 1)
#' r1(y[i], w[i], 0L, beta, alpha1_star = 0, tau1 = 1, w_min = -1, w_max = 1)
#'
#' # When the event is observed the two agree exactly
#' j <- which(delta == 1)[1]
#' r2(y[j], w[j], 1L, beta, alpha1_star = 0, tau1 = 1, w_min = -1, w_max = 1)
#' r1(y[j], w[j], 1L, beta, alpha1_star = 0, tau1 = 1, w_min = -1, w_max = 1)
#' @export
r2 = function(y, w, delta, beta,
              alpha1_star, tau1,
              z_c = NULL, z_d = NULL,
              w_min, w_max) {
  abs(y - m0(w, beta, z_c, z_d))
}
