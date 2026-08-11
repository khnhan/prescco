# =============================================================================
#  Conformal prediction intervals (covariate-aware, data-in)
#  ---------------------------------------------------------------------------
#  Three conformal comparison methods operating on given data, all optionally
#  covariate-aware (continuous z_c, discrete z_d):
#
#    * split_conformal_prediction_interval()  split conformal, with beta fit by
#                                        SPARCC on the calibration split.
#    * full_conformal_prediction_interval()                full conformal.
#    * jackknife_plus_prediction_interval()          jackknife+.
#
#  The vectorized centers/residuals m0_vec/m1_vec/r1_vec/r2_vec are the
#  fast counterparts of m0/m1/r1/r2, evaluated over whole test sets. Inference
#  and variance outputs from the original scripts are omitted: split conformal
#  returns the half-length; full conformal and jackknife+ return the empirical
#  coverage rate (and half-length) directly.
#
#  Package imports are declared centrally in R/censcovpred-package.R.
# =============================================================================

# ---- Vectorized centers and residuals ----

#' Vectorized interval center m0
#'
#' @param w Numeric vector of event-time covariate values.
#' @param beta Outcome coefficient vector.
#' @param z_c,z_d Optional covariate matrices aligned with w.
#' @return Numeric vector of centers.
#' @noRd
m0_vec <- function(w, beta, z_c = NULL, z_d = NULL) {
  w <- as.numeric(w)
  n <- length(w)
  
  ## coerce z_c, z_d to matrices with n rows
  if (is.null(z_c)) {
    z_c_mat <- matrix(0, nrow = n, ncol = 0); p_c <- 0
  } else if (is.null(dim(z_c))) {
    if (length(z_c) != n) stop("z_c must have length n")
    z_c_mat <- matrix(as.numeric(z_c), ncol = 1); p_c <- 1
  } else {
    if (nrow(z_c) != n) stop("z_c must have n rows")
    z_c_mat <- as.matrix(z_c); p_c <- ncol(z_c_mat)
  }
  
  if (is.null(z_d)) {
    z_d_mat <- matrix(0, nrow = n, ncol = 0); p_d <- 0
  } else if (is.null(dim(z_d))) {
    if (length(z_d) != n) stop("z_d must have length n")
    z_d_mat <- matrix(as.numeric(z_d), ncol = 1); p_d <- 1
  } else {
    if (nrow(z_d) != n) stop("z_d must have n rows")
    z_d_mat <- as.matrix(z_d); p_d <- ncol(z_d_mat)
  }
  
  ## indices in beta: (1, X, Z_c, Z_d, X Z_c, X Z_d)
  idx_x   <- 2L
  idx_zc  <- if (p_c > 0) 2L + seq_len(p_c) else integer(0)
  idx_zd  <- if (p_d > 0) 2L + p_c + seq_len(p_d) else integer(0)
  idx_xzc <- if (p_c > 0) 2L + p_c + p_d + seq_len(p_c) else integer(0)
  idx_xzd <- if (p_d > 0) 2L + p_c + p_d + p_c + seq_len(p_d) else integer(0)
  
  ## intercept and slope in X
  # rep_len keeps these length n even with no covariates, so that logical
  # subsetting below (intercept[A], slope[C], ...) stays well defined.
  intercept <- rep_len(beta[1] +
    (if (p_c > 0) as.vector(z_c_mat %*% beta[idx_zc]) else 0) +
    (if (p_d > 0) as.vector(z_d_mat %*% beta[idx_zd]) else 0), n)
  
  slope <- rep_len(beta[idx_x] +
    (if (p_c > 0) as.vector(z_c_mat %*% beta[idx_xzc]) else 0) +
    (if (p_d > 0) as.vector(z_d_mat %*% beta[idx_xzd]) else 0), n)
  
  intercept + slope * w
}
#' Vectorized interval center m1
#'
#' @param w,delta Numeric vectors of times and event indicators.
#' @param beta Outcome coefficient vector.
#' @param alpha1_star,tau1 eta1 mean coefficients and SD.
#' @param z_c,z_d Optional covariate matrices.
#' @param w_max Upper truncation bound.
#' @return Numeric vector of centers.
#' @noRd
m1_vec <- function(w, delta, beta, alpha1_star, tau1 = 1,
                   z_c = NULL, z_d = NULL, w_max = 1) {
  w     <- as.numeric(w)
  delta <- as.integer(delta)
  n     <- length(w)
  
  ## coerce z_c, z_d same way as in m0_vec
  if (is.null(z_c)) {
    z_c_mat <- matrix(0, nrow = n, ncol = 0); p_c <- 0
  } else if (is.null(dim(z_c))) {
    if (length(z_c) != n) stop("z_c must have length n")
    z_c_mat <- matrix(as.numeric(z_c), ncol = 1); p_c <- 1
  } else {
    if (nrow(z_c) != n) stop("z_c must have n rows")
    z_c_mat <- as.matrix(z_c); p_c <- ncol(z_c_mat)
  }
  
  if (is.null(z_d)) {
    z_d_mat <- matrix(0, nrow = n, ncol = 0); p_d <- 0
  } else if (is.null(dim(z_d))) {
    if (length(z_d) != n) stop("z_d must have length n")
    z_d_mat <- matrix(as.numeric(z_d), ncol = 1); p_d <- 1
  } else {
    if (nrow(z_d) != n) stop("z_d must have n rows")
    z_d_mat <- as.matrix(z_d); p_d <- ncol(z_d_mat)
  }
  
  ## indices for outcome model (same as m0_vec)
  idx_x   <- 2L
  idx_zc  <- if (p_c > 0) 2L + seq_len(p_c) else integer(0)
  idx_zd  <- if (p_d > 0) 2L + p_c + seq_len(p_d) else integer(0)
  idx_xzc <- if (p_c > 0) 2L + p_c + p_d + seq_len(p_c) else integer(0)
  idx_xzd <- if (p_d > 0) 2L + p_c + p_d + p_c + seq_len(p_d) else integer(0)
  
  # rep_len keeps these length n even with no covariates, so that logical
  # subsetting below (intercept[A], slope[C], ...) stays well defined.
  intercept <- rep_len(beta[1] +
    (if (p_c > 0) as.vector(z_c_mat %*% beta[idx_zc]) else 0) +
    (if (p_d > 0) as.vector(z_d_mat %*% beta[idx_zd]) else 0), n)
  
  slope <- rep_len(beta[idx_x] +
    (if (p_c > 0) as.vector(z_c_mat %*% beta[idx_xzc]) else 0) +
    (if (p_d > 0) as.vector(z_d_mat %*% beta[idx_xzd]) else 0), n)
  
  ## X|Z ~ N(mu_x_prior, tau1^2), with mu_x_prior = (1, z_c, z_d) %*% alpha1_star
  alpha1_star <- as.numeric(alpha1_star)
  if (length(alpha1_star) != (1 + p_c + p_d))
    stop("length(alpha1_star) must be 1 + p_c + p_d")
  
  idx_a_zc <- 1L + seq_len(p_c)
  idx_a_zd <- 1L + p_c + seq_len(p_d)
  
  mu_x_prior <- rep_len(alpha1_star[1] +
    (if (p_c > 0) as.vector(z_c_mat %*% alpha1_star[idx_a_zc]) else 0) +
    (if (p_d > 0) as.vector(z_d_mat %*% alpha1_star[idx_a_zd]) else 0), n)
  
  out <- numeric(n)
  
  # Masks
  A <- (delta == 1L)
  B <- (delta == 0L) & (w >= w_max)
  C <- (delta == 0L) & (w <  w_max)
  
  # A, B: just plug in w into m0
  if (any(A)) out[A] <- intercept[A] + slope[A] * w[A]
  if (any(B)) out[B] <- intercept[B] + slope[B] * w[B]
  
  # C: truncated mean E[X | w < X < w_max], X ~ N(mu_x_prior, tau1^2)
  if (any(C)) {
    s  <- tau1
    aC <- (w[C]        - mu_x_prior[C]) / s
    bC <- (w_max       - mu_x_prior[C]) / s
    
    denom <- pnorm(bC) - pnorm(aC)          # P(w < X < w_max)
    numer <- dnorm(aC) - dnorm(bC)          # φ(aC) - φ(bC)
    
    tiny <- 1e-12
    good <- denom > tiny
    
    x_exp <- numeric(sum(C))
    x_exp[good]  <- mu_x_prior[C][good] + s * numer[good] / denom[good]
    x_exp[!good] <- w[C][!good]   # fallback if interval prob ~ 0
    
    out[C] <- intercept[C] + slope[C] * x_exp
  }
  
  out
}
#' Vectorized residual r1
#'
#' @param y Numeric vector of outcomes.
#' @param w,delta Numeric vectors of observed times and event indicators.
#' @param beta Outcome coefficient vector.
#' @param alpha1_star Coefficient vector for the eta1 mean.
#' @param tau1 Standard deviation of X.
#' @param z_c,z_d Optional covariate matrices.
#' @param w_max Upper truncation bound.
#' @return Numeric vector |y - m1|.
#' @noRd
r1_vec <- function(y, w, delta, beta, alpha1_star, tau1 = 1,
                   z_c = NULL, z_d = NULL, w_max = 1) {
  y <- as.numeric(y)
  if (length(y) != length(w) || length(w) != length(delta))
    stop("y, w, delta must have same length.")
  
  abs(y - m1_vec(w, delta, beta, alpha1_star, tau1,
                 z_c = z_c, z_d = z_d, w_max = w_max))
}
#' Vectorized residual r2
#'
#' @param y,w,delta Numeric vectors.
#' @param beta Outcome coefficient vector.
#' @param alpha1_star,tau1 Unused; kept for a common signature.
#' @param z_c,z_d Optional covariate matrices.
#' @param w_max Unused; kept for a common signature.
#' @return Numeric vector |y - m0(w)|.
#' @noRd
r2_vec <- function(y, w, delta, beta, alpha1_star, tau1 = 1,
                   z_c = NULL, z_d = NULL, w_max = 1) {
  # alpha1_star, tau1, w_max kept for signature compatibility
  y <- as.numeric(y)
  if (length(y) != length(w) || length(w) != length(delta))
    stop("y, w, delta must have same length.")
  
  abs(y - m0_vec(w, beta, z_c = z_c, z_d = z_d))
}


# ---- Shared coverage-rate helper ----

#' Empirical coverage rate of fixed half-lengths on a test set
#'
#' Used by the methods that produce a single half-length from training data
#' (semiparametric, split conformal) to attach a coverage rate when the caller
#' supplies test data. Returns a vector of NAs when no test set is given.
#'
#' @param zeta Named half-length vector; names select the residual
#'   (\code{"r1"}, \code{"r2"}, \code{"r1star"}).
#' @param beta Fitted outcome coefficients.
#' @param test_y_data,test_w_data,test_delta_data Test set; \code{NULL} to skip.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param alpha1,tau1 Nuisance mean coefficients and SD for r1 / r2.
#' @param alpha1_star_r,tau1_r Residual center and SD for r1*.
#' @param w_min,w_max Support bounds.
#' @return Named numeric vector matching \code{names(zeta)}.
#' @noRd
.coverage_on_test <- function(zeta, beta,
                              test_y_data, test_w_data, test_delta_data,
                              test_z_c_data = NULL, test_z_d_data = NULL,
                              alpha1, tau1, alpha1_star_r, tau1_r,
                              w_min, w_max) {
  if (is.null(test_y_data) || is.null(test_w_data) || is.null(test_delta_data)) {
    return(stats::setNames(rep(NA_real_, length(zeta)), names(zeta)))
  }

  vapply(names(zeta), function(nm) {
    res_fun <- if (nm == "r2") r2 else r1
    center  <- if (nm == "r1star") alpha1_star_r else alpha1
    scale   <- if (nm == "r1star") tau1_r        else tau1
    prediction_coverage_rate(
      zeta = zeta[[nm]], r = res_fun,
      test_y_data = test_y_data, test_w_data = test_w_data,
      test_delta_data = test_delta_data,
      test_z_c_data = test_z_c_data, test_z_d_data = test_z_d_data,
      beta = beta, alpha1_star_r = center, tau1 = scale,
      w_min = w_min, w_max = w_max
    )
  }, numeric(1))
}


# ---- Split conformal ----

#' Split conformal prediction half-length
#'
#' Fits beta by SPARCC on the first split, then calibrates the half-length on
#' the held-out split using empirical residual quantiles. Covariate-aware.
#' @param y_data,w_data,delta_data Observed outcome, time, event indicator.
#' @param z_c_data,z_d_data Optional covariates.
#' @param alpha1,alpha2 Nuisance mean coefficient vectors.
#' @param alpha1_star_r Center used inside the residual.
#' @param alpha Miscoverage level.
#' @param split_rate Fraction used for fitting (rest calibrates).
#' @param tau1,tau2,tau1_r Standard deviations (tau1_r for the residual center).
#' @param test_y_data,test_w_data,test_delta_data Optional test set. When
#'   supplied, \code{coverage_rate} is evaluated on it; otherwise it is
#'   \code{NA}.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param tt,m Integration grid size and covariate-grid resolution.
#' @param w_min,w_max Truncation bounds.
#' @return A list with \code{method}, \code{alpha}, \code{residual},
#'   \code{zeta} (named half-lengths for r1, r2, r1*), \code{coverage_rate}
#'   (\code{NA} unless test data is supplied), \code{beta}, \code{sigma},
#'   \code{n_fit}, and \code{n_calibrate}.
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 200
#' x  <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' cc <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w  <- pmin(x, cc)
#' delta <- as.integer(x <= cc)
#' y  <- 3 * x + rnorm(n)          # beta = c(0, 3)
#'
#' a1 <- unname(find_alpha1_MLE(w, delta, w_min = -1, w_max = 1))
#' a2 <- unname(find_alpha2_MLE(w, delta, w_min = -1, w_max = 1))
#'
#' fit <- split_conformal_prediction_interval(
#'   y, w, delta,
#'   alpha1 = a1[1], alpha2 = a2[1], alpha1_star_r = a1[1],
#'   alpha = 0.1, tau1 = a1[2], tau2 = a2[2], tau1_r = a1[2],
#'   w_min = -1, w_max = 1
#' )
#' fit$zeta            # half-lengths for r1, r2 and r1*
#' fit$coverage_rate   # NA unless test data is supplied
#'
#' # Supply a test set and the coverage rate is filled in
#' set.seed(2)
#' xt <- truncnorm::rtruncnorm(300, a = -1, b = 1, mean = 0, sd = 1)
#' ct <- truncnorm::rtruncnorm(300, a = -1, b = 1, mean = 1, sd = 1)
#' fit2 <- split_conformal_prediction_interval(
#'   y, w, delta,
#'   alpha1 = a1[1], alpha2 = a2[1], alpha1_star_r = a1[1],
#'   alpha = 0.1, tau1 = a1[2], tau2 = a2[2], tau1_r = a1[2],
#'   test_y_data = 3 * xt + rnorm(300),
#'   test_w_data = pmin(xt, ct),
#'   test_delta_data = as.integer(xt <= ct),
#'   w_min = -1, w_max = 1
#' )
#' fit2$coverage_rate
#' }
#' @export
split_conformal_prediction_interval <- function(y_data, w_data, delta_data,
                                           z_c_data = NULL, z_d_data = NULL,
                                           alpha1, alpha2, alpha1_star_r,
                                           alpha,
                                           split_rate = 0.5,
                                           tau1, tau2, tau1_r,
                                           test_y_data = NULL,
                                           test_w_data = NULL,
                                           test_delta_data = NULL,
                                           test_z_c_data = NULL,
                                           test_z_d_data = NULL,
                                           tt = 20, m = 20,
                                           w_min, w_max) {
  # basic checks
  n <- length(y_data)
  if (length(w_data) != n || length(delta_data) != n) {
    stop("y_data, w_data, delta_data must have same length.")
  }
  
  # split index
  n1 = floor(n * split_rate)
  set.seed(1)
  train_idx <- sample(1:n, size = n1)
  if (n1 < length(beta)) {
    warning("Very small estimation sample relative to beta dimension.")
  }
  
  # X-grid
  x_a <- seq(w_min, w_max, length.out = m)
  
  # slice z_c, z_d for first half (pass through as-is, robust to vector/matrix/NULL)
  z_c_1 <-  if (is.null(z_c_data)) {
    NULL
  } else if (is.null(dim(z_c_data))) {
    z_c_data[train_idx]
  } else {
    z_c_data[train_idx, ]
  }
  
  z_d_1 <-     if (is.null(z_d_data)) {
    NULL
  } else if (is.null(dim(z_d_data))) {
    z_d_data[train_idx]
  } else {
    z_d_data[train_idx, ]
  }
  
  result_sparcc = find_beta_sparcc(y_data[train_idx], w_data[train_idx], delta_data[train_idx],
                                            z_c_1,
                                            z_d_1,
                                            alpha1,
                                            alpha2,
                                            tau1,
                                            tau2,
                                            tt, m,
                                            w_min, w_max)
  beta_hat       = result_sparcc$beta_hat
  sigma_hat      = result_sparcc$sigma_hat
  ## helpers to extract z_c, z_d for obs i
  get_zc_i <- function(i) {
    if (is.null(z_c_data)) {
      NULL
    } else if (is.null(dim(z_c_data))) {
      z_c_data[i]
    } else {
      z_c_data[i, ]
    }
  }
  get_zd_i <- function(i) {
    if (is.null(z_d_data)) {
      NULL
    } else if (is.null(dim(z_d_data))) {
      z_d_data[i]
    } else {
      z_d_data[i, ]
    }
  }
  
  ## 2. Residuals on second split
  r1_vec     <- numeric(n - n1)
  r2_vec     <- numeric(n - n1)
  r1star_vec <- numeric(n - n1)
  
  test_idx = (1:n)[-train_idx]
  for (idx in 1:(n-n1)) {
    i = test_idx[idx]
    z_c_i <- get_zc_i(i)
    z_d_i <- get_zd_i(i)
    
    r1_vec[idx] <- r1(
      y     = y_data[i],
      w     = w_data[i],
      delta = delta_data[i],
      beta  = beta_hat,
      alpha1_star = alpha1,
      tau1  = tau1,
      z_c   = z_c_i,
      z_d   = z_d_i,
      w_min = w_min,
      w_max = w_max
    )
    
    r2_vec[idx] <- r2(
      y     = y_data[i],
      w     = w_data[i],
      delta = delta_data[i],
      beta  = beta_hat,
      alpha1_star = alpha1,
      tau1  = tau1,
      z_c   = z_c_i,
      z_d   = z_d_i,
      w_min = w_min,
      w_max = w_max
    )
    
    r1star_vec[idx] <- r1(
      y     = y_data[i],
      w     = w_data[i],
      delta = delta_data[i],
      beta  = beta_hat,
      alpha1_star = alpha1_star_r,
      tau1  = tau1_r,
      z_c   = z_c_i,
      z_d   = z_d_i,
      w_min = w_min,
      w_max = w_max
    )
  }
  
  ## 3. Split-conformal quantiles
  zeta_r1     <- stats::quantile(c(r1_vec, Inf),     1 - alpha)
  zeta_r2     <- stats::quantile(c(r2_vec, Inf),     1 - alpha)
  zeta_r1star <- stats::quantile(c(r1star_vec, Inf), 1 - alpha)
  
  zeta <- c(r1 = unname(zeta_r1), r2 = unname(zeta_r2),
            r1star = unname(zeta_r1star))

  ## Optional empirical coverage rate on a supplied test set.
  coverage_rate <- .coverage_on_test(
    zeta = zeta, beta = beta_hat,
    test_y_data = test_y_data, test_w_data = test_w_data,
    test_delta_data = test_delta_data,
    test_z_c_data = test_z_c_data, test_z_d_data = test_z_d_data,
    alpha1 = alpha1, tau1 = tau1,
    alpha1_star_r = alpha1_star_r, tau1_r = tau1_r,
    w_min = w_min, w_max = w_max
  )

  list(
    method        = "split_conformal",
    alpha         = alpha,
    residual      = c("r1", "r2", "r1star"),
    zeta          = zeta,
    coverage_rate = coverage_rate,
    beta          = beta_hat,
    sigma         = sigma_hat,
    n_fit         = length(train_idx),
    n_calibrate   = n - n1
  )
}


# ---- Full conformal ----

#' Full conformal prediction coverage
#'
#' For each test point, augments the training data and computes conformal
#' ranks to decide inclusion; returns empirical coverage. Covariate-aware.
#' @param alpha Miscoverage level.
#' @param y_data,w_data,delta_data Training outcome, time, indicator.
#' @param z_c_data,z_d_data Optional training covariates.
#' @param test_y_data,test_w_data,test_delta_data Test outcome, time, indicator.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param alpha1,tau1 eta1 mean coefficients and SD.
#' @param alpha1_star_r,tau1_r Residual center coefficients and SD.
#' @param w_min,w_max Truncation bounds.
#' @return A list with \code{method}, \code{alpha}, \code{residual},
#'   \code{zeta} (half-lengths averaged over the test set), and
#'   \code{coverage_rate}, plus \code{n_train} and \code{n_test}.
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 200
#' x  <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' cc <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w  <- pmin(x, cc)
#' delta <- as.integer(x <= cc)
#' y  <- 3 * x + rnorm(n)          # beta = c(0, 3)
#'
#' a1 <- unname(find_alpha1_MLE(w, delta, w_min = -1, w_max = 1))
#'
#' # Full conformal refits at every test point, so a test set is required
#' # and the cost grows quickly -- keep it small.
#' set.seed(2)
#' xt <- truncnorm::rtruncnorm(20, a = -1, b = 1, mean = 0, sd = 1)
#' ct <- truncnorm::rtruncnorm(20, a = -1, b = 1, mean = 1, sd = 1)
#'
#' fit <- full_conformal_prediction_interval(
#'   alpha = 0.1,
#'   y, w, delta,
#'   test_y_data = 3 * xt + rnorm(20),
#'   test_w_data = pmin(xt, ct),
#'   test_delta_data = as.integer(xt <= ct),
#'   alpha1 = a1[1], tau1 = a1[2],
#'   alpha1_star_r = a1[1], tau1_r = a1[2],
#'   w_min = -1, w_max = 1
#' )
#' fit$zeta            # averaged over the test set
#' fit$coverage_rate
#' }
#' @export
full_conformal_prediction_interval <- function(alpha,
                             y_data, w_data, delta_data,
                             z_c_data = NULL, z_d_data = NULL,
                             test_y_data, test_w_data, test_delta_data,
                             test_z_c_data = NULL, test_z_d_data = NULL,
                             alpha1,
                             tau1,
                             alpha1_star_r,
                             tau1_r,
                             w_min, w_max) {
  n      <- length(y_data)
  n_test <- length(test_y_data)
  
  if (length(w_data) != n || length(delta_data) != n)
    stop("y_data, w_data, delta_data must have same length.")
  if (length(test_w_data) != n_test || length(test_delta_data) != n_test)
    stop("test_y_data, test_w_data, test_delta_data must have same length.")
  
  ## ---------- Coerce training z_c, z_d ----------
  if (is.null(z_c_data)) {
    z_c_mat <- matrix(0, nrow = n, ncol = 0)
    p_c     <- 0
  } else if (is.null(dim(z_c_data))) {
    if (length(z_c_data) != n) stop("z_c_data must have length n.")
    z_c_mat <- matrix(as.numeric(z_c_data), ncol = 1)
    p_c     <- 1
  } else {
    if (nrow(z_c_data) != n) stop("z_c_data must have n rows.")
    z_c_mat <- as.matrix(z_c_data)
    p_c     <- ncol(z_c_mat)
  }
  
  if (is.null(z_d_data)) {
    z_d_mat <- matrix(0, nrow = n, ncol = 0)
    p_d     <- 0
  } else if (is.null(dim(z_d_data))) {
    if (length(z_d_data) != n) stop("z_d_data must have length n.")
    z_d_mat <- matrix(as.numeric(z_d_data), ncol = 1)
    p_d     <- 1
  } else {
    if (nrow(z_d_data) != n) stop("z_d_data must have n rows.")
    z_d_mat <- as.matrix(z_d_data)
    p_d     <- ncol(z_d_mat)
  }
  
  ## ---------- Coerce test z_c, z_d ----------
  if (is.null(test_z_c_data)) {
    test_z_c_mat <- matrix(0, nrow = n_test, ncol = 0)
    if (p_c != 0) stop("Training has z_c but test_z_c_data is NULL.")
  } else if (is.null(dim(test_z_c_data))) {
    if (length(test_z_c_data) != n_test) stop("test_z_c_data must have length n_test.")
    if (p_c > 1) stop("test_z_c_data must be matrix with p_c columns.")
    test_z_c_mat <- matrix(as.numeric(test_z_c_data), ncol = p_c)
  } else {
    if (nrow(test_z_c_data) != n_test) stop("test_z_c_data must have n_test rows.")
    test_z_c_mat <- as.matrix(test_z_c_data)
    if (ncol(test_z_c_mat) != p_c) stop("test_z_c_data has wrong number of columns.")
  }
  
  if (is.null(test_z_d_data)) {
    test_z_d_mat <- matrix(0, nrow = n_test, ncol = 0)
    if (p_d != 0) stop("Training has z_d but test_z_d_data is NULL.")
  } else if (is.null(dim(test_z_d_data))) {
    if (length(test_z_d_data) != n_test) stop("test_z_d_data must have length n_test.")
    if (p_d > 1) stop("test_z_d_data must be matrix with p_d columns.")
    test_z_d_mat <- matrix(as.numeric(test_z_d_data), ncol = p_d)
  } else {
    if (nrow(test_z_d_data) != n_test) stop("test_z_d_data must have n_test rows.")
    test_z_d_mat <- as.matrix(test_z_d_data)
    if (ncol(test_z_d_mat) != p_d) stop("test_z_d_data has wrong number of columns.")
  }
  
  ## ---------- Loop over test points ----------
  pred_cvg_list <- lapply(seq_len(n_test), function(j) {
    ## Augment with one test point j
    y_aug     <- c(test_y_data[j],     y_data)
    w_aug     <- c(test_w_data[j],     w_data)
    delta_aug <- c(test_delta_data[j], delta_data)
    
    z_c_aug <- if (p_c > 0) {
      rbind(test_z_c_mat[j, , drop = FALSE], z_c_mat)
    } else {
      matrix(0, nrow = length(y_aug), ncol = 0)
    }
    z_d_aug <- if (p_d > 0) {
      rbind(test_z_d_mat[j, , drop = FALSE], z_d_mat)
    } else {
      matrix(0, nrow = length(y_aug), ncol = 0)
    }
    
    ## Complete cases (Δ = 1)
    cc_idx <- (delta_aug == 1)
    if (!any(cc_idx)) stop("No complete cases in augmented sample.")
    
    y_cc   <- y_aug[cc_idx]
    x_cc   <- w_aug[cc_idx]
    z_c_cc <- z_c_aug[cc_idx, , drop = FALSE]
    z_d_cc <- z_d_aug[cc_idx, , drop = FALSE]
    
    ## Design: 1 + X + Z_c + Z_d + X:Z_c + X:Z_d
    X <- cbind(
      1,
      x_cc,
      if (p_c > 0) z_c_cc else NULL,
      if (p_d > 0) z_d_cc else NULL,
      if (p_c > 0) sweep(z_c_cc, 1, x_cc, `*`) else NULL,
      if (p_d > 0) sweep(z_d_cc, 1, x_cc, `*`) else NULL
    )
    Y <- y_cc
    
    beta_cc <- as.numeric(solve(crossprod(X), crossprod(X, Y)))
    
    ## Compute residual radii for all augmented points
    r1_vector <- r1_vec(
      y     = y_aug,
      w     = w_aug,
      delta = delta_aug,
      beta  = beta_cc,
      alpha1_star = alpha1,
      tau1  = tau1,
      z_c   = z_c_aug,
      z_d   = z_d_aug,
      w_max = w_max
    )
    
    r2_vector <- r2_vec(
      y     = y_aug,
      w     = w_aug,
      delta = delta_aug,
      beta  = beta_cc,
      alpha1_star = alpha1,
      tau1  = tau1,
      z_c   = z_c_aug,
      z_d   = z_d_aug,
      w_max = w_max
    )
    
    r1star_vector <- r1_vec(
      y     = y_aug,
      w     = w_aug,
      delta = delta_aug,
      beta  = beta_cc,
      alpha1_star = alpha1_star_r,
      tau1  = tau1_r,
      z_c   = z_c_aug,
      z_d   = z_d_aug,
      w_max = w_max
    )
    
    # unname() on the quantiles matters: quantile() returns a named value
    # ("90%"), and c(zeta_r1 = <named>) composes the two into "zeta_r1.90%".
    c(
      cvg_r1     = as.numeric(r1_vector[1]     <= quantile(r1_vector,     1 - alpha)),
      cvg_r2     = as.numeric(r2_vector[1]     <= quantile(r2_vector,     1 - alpha)),
      cvg_r1star = as.numeric(r1star_vector[1] <= quantile(r1star_vector, 1 - alpha)),
      zeta_r1     = unname(quantile(r1_vector,     1 - alpha)),
      zeta_r2     = unname(quantile(r2_vector,     1 - alpha)),
      zeta_r1star = unname(quantile(r1star_vector, 1 - alpha))
    )
  })
  
  pred_cvg_mat <- do.call(rbind, pred_cvg_list)
  avg <- apply(pred_cvg_mat, 2, mean)

  list(
    method        = "full_conformal",
    alpha         = alpha,
    residual      = c("r1", "r2", "r1star"),
    zeta          = c(r1     = unname(avg[["zeta_r1"]]),
                      r2     = unname(avg[["zeta_r2"]]),
                      r1star = unname(avg[["zeta_r1star"]])),
    coverage_rate = c(r1     = unname(avg[["cvg_r1"]]),
                      r2     = unname(avg[["cvg_r2"]]),
                      r1star = unname(avg[["cvg_r1star"]])),
    n_train       = length(y_data),
    n_test        = length(test_y_data)
  )
}


# ---- Jackknife+ ----

#' Jackknife+ prediction coverage
#'
#' Builds leave-one-out complete-case fits and forms jackknife+ prediction
#' intervals; returns empirical coverage and half-length. Covariate-aware.
#' @param alpha Miscoverage level.
#' @param y_data,w_data,delta_data Training outcome, time, indicator.
#' @param z_c_data,z_d_data Optional training covariates.
#' @param test_y_data,test_w_data,test_delta_data Test outcome, time, indicator.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param alpha1,tau1 eta1 mean coefficients and SD.
#' @param alpha1_star_r,tau1_r Residual center coefficients and SD.
#' @param w_min,w_max Truncation bounds.
#' @return A list with \code{method}, \code{alpha}, \code{residual},
#'   \code{zeta} (half-lengths averaged over the test set), and
#'   \code{coverage_rate}, plus \code{n_train} and \code{n_test}.
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 200
#' x  <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' cc <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w  <- pmin(x, cc)
#' delta <- as.integer(x <= cc)
#' y  <- 3 * x + rnorm(n)          # beta = c(0, 3)
#'
#' a1 <- unname(find_alpha1_MLE(w, delta, w_min = -1, w_max = 1))
#'
#' set.seed(2)
#' xt <- truncnorm::rtruncnorm(20, a = -1, b = 1, mean = 0, sd = 1)
#' ct <- truncnorm::rtruncnorm(20, a = -1, b = 1, mean = 1, sd = 1)
#'
#' # Same arguments as full_conformal_prediction_interval(); builds the
#' # interval from leave-one-out fits instead of refitting per test point.
#' fit <- jackknife_plus_prediction_interval(
#'   alpha = 0.1,
#'   y, w, delta,
#'   test_y_data = 3 * xt + rnorm(20),
#'   test_w_data = pmin(xt, ct),
#'   test_delta_data = as.integer(xt <= ct),
#'   alpha1 = a1[1], tau1 = a1[2],
#'   alpha1_star_r = a1[1], tau1_r = a1[2],
#'   w_min = -1, w_max = 1
#' )
#' fit$zeta
#' fit$coverage_rate
#' }
#' @export
jackknife_plus_prediction_interval = function(alpha,
                                  y_data, w_data, delta_data,
                                  z_c_data = NULL, z_d_data = NULL,
                                  test_y_data, test_w_data, test_delta_data,
                                  test_z_c_data = NULL, test_z_d_data = NULL,
                                  alpha1, tau1,
                                  alpha1_star_r, tau1_r,
                                  w_min, w_max) {
  n      <- length(y_data)
  n_test <- length(test_y_data)
  
  if (length(w_data) != n || length(delta_data) != n)
    stop("y_data, w_data, delta_data must have same length.")
  if (length(test_w_data) != n_test || length(test_delta_data) != n_test)
    stop("test_y_data, test_w_data, test_delta_data must have same length.")
  
  ## ---- coerce training z_c, z_d to matrices ----
  if (is.null(z_c_data)) {
    z_c_mat <- matrix(0, nrow = n, ncol = 0)
    p_c     <- 0
  } else if (is.null(dim(z_c_data))) {
    if (length(z_c_data) != n) stop("z_c_data must have length n.")
    z_c_mat <- matrix(as.numeric(z_c_data), ncol = 1)
    p_c     <- 1
  } else {
    if (nrow(z_c_data) != n) stop("z_c_data must have n rows.")
    z_c_mat <- as.matrix(z_c_data)
    p_c     <- ncol(z_c_mat)
  }
  
  if (is.null(z_d_data)) {
    z_d_mat <- matrix(0, nrow = n, ncol = 0)
    p_d     <- 0
  } else if (is.null(dim(z_d_data))) {
    if (length(z_d_data) != n) stop("z_d_data must have length n.")
    z_d_mat <- matrix(as.numeric(z_d_data), ncol = 1)
    p_d     <- 1
  } else {
    if (nrow(z_d_data) != n) stop("z_d_data must have n rows.")
    z_d_mat <- as.matrix(z_d_data)
    p_d     <- ncol(z_d_mat)
  }
  
  ## ---- coerce test z_c, z_d to matrices ----
  if (is.null(test_z_c_data)) {
    test_z_c_mat <- matrix(0, nrow = n_test, ncol = 0)
    if (p_c != 0) stop("Training has z_c but test_z_c_data is NULL.")
  } else if (is.null(dim(test_z_c_data))) {
    if (length(test_z_c_data) != n_test) stop("test_z_c_data must have length n_test.")
    if (p_c > 1) stop("test_z_c_data must be matrix with p_c columns.")
    test_z_c_mat <- matrix(as.numeric(test_z_c_data), ncol = p_c)
  } else {
    if (nrow(test_z_c_data) != n_test) stop("test_z_c_data must have n_test rows.")
    test_z_c_mat <- as.matrix(test_z_c_data)
    if (ncol(test_z_c_mat) != p_c) stop("test_z_c_data has wrong number of columns.")
  }
  
  if (is.null(test_z_d_data)) {
    test_z_d_mat <- matrix(0, nrow = n_test, ncol = 0)
    if (p_d != 0) stop("Training has z_d but test_z_d_data is NULL.")
  } else if (is.null(dim(test_z_d_data))) {
    if (length(test_z_d_data) != n_test) stop("test_z_d_data must have length n_test.")
    if (p_d > 1) stop("test_z_d_data must be matrix with p_d columns.")
    test_z_d_mat <- matrix(as.numeric(test_z_d_data), ncol = p_d)
  } else {
    if (nrow(test_z_d_data) != n_test) stop("test_z_d_data must have n_test rows.")
    test_z_d_mat <- as.matrix(test_z_d_data)
    if (ncol(test_z_d_mat) != p_d) stop("test_z_d_data has wrong number of columns.")
  }
  
  ## ---- jackknife+ loop over training i ----
  m_list <- lapply(seq_len(n), function(i) {
    cc_idx <- delta_data[-i] == 1
    if (!any(cc_idx)) stop("No complete cases in training set (minus i).")
    
    y_cc   <- y_data[-i][cc_idx]
    x_cc   <- w_data[-i][cc_idx]
    z_c_cc <- z_c_mat[-i, , drop = FALSE][cc_idx, , drop = FALSE]
    z_d_cc <- z_d_mat[-i, , drop = FALSE][cc_idx, , drop = FALSE]
    
    X <- cbind(
      1,
      x_cc,
      if (p_c > 0) z_c_cc else NULL,
      if (p_d > 0) z_d_cc else NULL,
      if (p_c > 0) sweep(z_c_cc, 1, x_cc, `*`) else NULL,
      if (p_d > 0) sweep(z_d_cc, 1, x_cc, `*`) else NULL
    )
    Y <- y_cc
    
    beta_cc <- solve(crossprod(X), crossprod(X, Y))
    
    ## predictions for all test points
    m1_vec_test <- m1_vec(
      w      = test_w_data,
      delta  = test_delta_data,
      z_c    = test_z_c_mat,
      z_d    = test_z_d_mat,
      beta   = beta_cc,
      alpha1_star = alpha1,
      tau1   = tau1,
      w_max  = w_max
    )
    
    m2_vec_test <- m0_vec(
      w    = test_w_data,
      z_c  = test_z_c_mat,
      z_d  = test_z_d_mat,
      beta = beta_cc
    )
    
    m1star_vec_test <- m1_vec(
      w      = test_w_data,
      delta  = test_delta_data,
      z_c    = test_z_c_mat,
      z_d    = test_z_d_mat,
      beta   = beta_cc,
      alpha1_star = alpha1_star_r,
      tau1   = tau1_r,
      w_max  = w_max
    )
    
    ## residuals for left-out training observation i
    z_c_i <- if (p_c > 0) z_c_mat[i, ] else NULL
    z_d_i <- if (p_d > 0) z_d_mat[i, ] else NULL
    
    r1_i <- r1(
      y     = y_data[i],
      w     = w_data[i],
      delta = delta_data[i],
      z_c   = z_c_i,
      z_d   = z_d_i,
      beta  = beta_cc,
      alpha1_star = alpha1,
      tau1  = tau1,
      w_min = w_min,
      w_max = w_max
    )
    
    r2_i <- r2(
      y     = y_data[i],
      w     = w_data[i],
      delta = delta_data[i],
      z_c   = z_c_i,
      z_d   = z_d_i,
      beta  = beta_cc,
      alpha1_star = alpha1,
      tau1  = tau1,
      w_min = w_min,
      w_max = w_max
    )
    
    r1star_i <- r1(
      y     = y_data[i],
      w     = w_data[i],
      delta = delta_data[i],
      z_c   = z_c_i,
      z_d   = z_d_i,
      beta  = beta_cc,
      alpha1_star = alpha1_star_r,
      tau1  = tau1_r,
      w_min = w_min,
      w_max = w_max
    )
    
    cbind(
      m1_vec_test - r1_i,
      m1_vec_test + r1_i,
      m2_vec_test - r2_i,
      m2_vec_test + r2_i,
      m1star_vec_test - r1star_i,
      m1star_vec_test + r1star_i
    )
  })
  
  ## m_vals: dim = (n_test, 6, n)
  m_vals <- simplify2array(m_list)
  
  lower <- apply(m_vals[, c(1, 3, 5), , drop = FALSE],
                 c(1, 2),
                 function(x) quantile(x, probs = alpha))
  upper <- apply(m_vals[, c(2, 4, 6), , drop = FALSE],
                 c(1, 2),
                 function(x) quantile(x, probs = 1 - alpha))
  
  ## coverage for each of the three radii
  cvg <- sapply(1:3, function(j) {
    mean(test_y_data >= lower[, j] & test_y_data <= upper[, j])
  })
  names(cvg) <- c("r1", "r2", "r1star")
  
  ## average half-lengths
  zeta_sol <- colMeans((upper - lower) / 2)
  names(zeta_sol) <- c("r1", "r2", "r1star")
  
  list(
    method        = "jackknife_plus",
    alpha         = alpha,
    residual      = c("r1", "r2", "r1star"),
    zeta          = zeta_sol,
    coverage_rate = cvg,
    n_train       = length(y_data),
    n_test        = length(test_y_data)
  )
}
