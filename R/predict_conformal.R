# =============================================================================
# Conformal prediction intervals
#
# Split conformal prediction, full conformal prediction, and jackknife+ for a
# right-censored covariate, with optional continuous and discrete covariates.
# The center and residual functions in predict_helpers.R are vectorized and are
# used directly throughout this file.
# =============================================================================


# Empirical coverage on an optional test set.

.coverage_on_test <- function(zeta, beta,
                              test_y_data, test_w_data, test_delta_data,
                              test_z_c_data = NULL, test_z_d_data = NULL,
                              alpha1, tau1,
                              alpha1_star_r, tau1_r,
                              w_min, w_max) {

  if (
    is.null(test_y_data) ||
    is.null(test_w_data) ||
    is.null(test_delta_data)
  ) {
    return(
      stats::setNames(
        rep(NA_real_, length(zeta)),
        names(zeta)
      )
    )
  }

  if (
    length(test_y_data) != length(test_w_data) ||
    length(test_w_data) != length(test_delta_data)
  ) {
    stop(
      "test_y_data, test_w_data, and test_delta_data must have the same length."
    )
  }

  vapply(
    names(zeta),
    function(nm) {

      if (nm == "r2") {

        r_vec <- r2(
          y = test_y_data,
          w = test_w_data,
          delta = test_delta_data,
          beta = beta,
          alpha1_star = alpha1,
          tau1 = tau1,
          z_c = test_z_c_data,
          z_d = test_z_d_data,
          w_min = w_min,
          w_max = w_max
        )

      } else {

        alpha1_use <- if (nm == "r1star") alpha1_star_r else alpha1
        tau1_use <- if (nm == "r1star") tau1_r else tau1

        r_vec <- r1(
          y = test_y_data,
          w = test_w_data,
          delta = test_delta_data,
          beta = beta,
          alpha1_star = alpha1_use,
          tau1 = tau1_use,
          z_c = test_z_c_data,
          z_d = test_z_d_data,
          w_min = w_min,
          w_max = w_max
        )
      }

      mean(r_vec <= zeta[[nm]])
    },
    numeric(1)
  )
}


# =============================================================================
# Split conformal prediction
# =============================================================================

#' Split conformal prediction interval
#'
#' Fits the outcome-model coefficients by SPARCC on one split and estimates
#' the prediction half-length from residuals on the calibration split.
#'
#' @param y_data,w_data,delta_data Outcome, observed time, and event indicator.
#' @param z_c_data,z_d_data Optional continuous and discrete covariates.
#' @param alpha Miscoverage level.
#' @param alpha1,tau1 Coefficients and standard deviation for the working model
#'   for \code{X | Z}.
#' @param alpha1_star_r,tau1_r Coefficients and SD for \code{r1*}.
#' @param alpha2,tau2 Coefficients and standard deviation for the working model
#'   for \code{C | Z}.
#' @param test_y_data,test_w_data,test_delta_data Optional test data for
#'   evaluating empirical coverage.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param split_rate Fraction of observations used to estimate the outcome
#'   model.
#' @param tt,m Integration grid size and covariate-grid resolution.
#' @param w_min,w_max Truncation bounds.
#'
#' @return A list containing the estimated half-lengths, optional empirical
#'   coverage rates, fitted outcome-model coefficients and scale, and split
#'   sizes.
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' c_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w_vec <- pmin(x_vec, c_vec)
#' delta_vec <- as.integer(x_vec <= c_vec)
#' y_vec <- 3 * x_vec + rnorm(n)
#'
#' alpha1_fit <- unname(find_alpha1_MLE(
#'   w_vec, delta_vec, w_min = -1, w_max = 1
#' ))
#' alpha2_fit <- unname(find_alpha2_MLE(
#'   w_vec, delta_vec, w_min = -1, w_max = 1
#' ))
#'
#' fit <- split_conformal_prediction_interval(
#'   y_data = y_vec,
#'   w_data = w_vec,
#'   delta_data = delta_vec,
#'   alpha = 0.1,
#'   alpha1 = alpha1_fit[1],
#'   tau1 = alpha1_fit[2],
#'   alpha1_star_r = alpha1_fit[1],
#'   tau1_r = alpha1_fit[2],
#'   alpha2 = alpha2_fit[1],
#'   tau2 = alpha2_fit[2],
#'   tt = 5,
#'   m = 5,
#'   w_min = -1,
#'   w_max = 1
#' )
#' fit$zeta
#' fit$coverage_rate
#'
#' @export
split_conformal_prediction_interval <- function(
    y_data, w_data, delta_data,
    z_c_data = NULL, z_d_data = NULL,
    alpha,
    alpha1, tau1,
    alpha1_star_r, tau1_r,
    alpha2, tau2,
    test_y_data = NULL,
    test_w_data = NULL,
    test_delta_data = NULL,
    test_z_c_data = NULL,
    test_z_d_data = NULL,
    split_rate = 0.5,
    tt = 20, m = 20,
    w_min, w_max) {

  n <- length(y_data)

  if (length(w_data) != n || length(delta_data) != n) {
    stop("y_data, w_data, and delta_data must have the same length.")
  }

  n_fit <- floor(n * split_rate)

  if (n_fit <= 0 || n_fit >= n) {
    stop("split_rate must leave observations in both splits.")
  }

  set.seed(1)
  fit_idx <- sample(seq_len(n), size = n_fit)
  calibration_idx <- setdiff(seq_len(n), fit_idx)

  z_c_fit <- if (is.null(z_c_data)) {
    NULL
  } else if (is.null(dim(z_c_data))) {
    z_c_data[fit_idx]
  } else {
    z_c_data[fit_idx, , drop = FALSE]
  }

  z_d_fit <- if (is.null(z_d_data)) {
    NULL
  } else if (is.null(dim(z_d_data))) {
    z_d_data[fit_idx]
  } else {
    z_d_data[fit_idx, , drop = FALSE]
  }

  fit_sparcc <- find_beta_sparcc(
    y_data = y_data[fit_idx],
    w_data = w_data[fit_idx],
    delta_data = delta_data[fit_idx],
    z_c_data = z_c_fit,
    z_d_data = z_d_fit,
    alpha1_star = alpha1,
    alpha2_star = alpha2,
    tau1 = tau1,
    tau2 = tau2,
    tt = tt,
    m = m,
    w_min = w_min,
    w_max = w_max
  )

  beta_hat <- fit_sparcc$beta_hat
  sigma_hat <- fit_sparcc$sigma_hat

  z_c_calibration <- if (is.null(z_c_data)) {
    NULL
  } else if (is.null(dim(z_c_data))) {
    z_c_data[calibration_idx]
  } else {
    z_c_data[calibration_idx, , drop = FALSE]
  }

  z_d_calibration <- if (is.null(z_d_data)) {
    NULL
  } else if (is.null(dim(z_d_data))) {
    z_d_data[calibration_idx]
  } else {
    z_d_data[calibration_idx, , drop = FALSE]
  }

  res_r1 <- r1(
    y = y_data[calibration_idx],
    w = w_data[calibration_idx],
    delta = delta_data[calibration_idx],
    beta = beta_hat,
    alpha1_star = alpha1,
    tau1 = tau1,
    z_c = z_c_calibration,
    z_d = z_d_calibration,
    w_min = w_min,
    w_max = w_max
  )

  res_r2 <- r2(
    y = y_data[calibration_idx],
    w = w_data[calibration_idx],
    delta = delta_data[calibration_idx],
    beta = beta_hat,
    alpha1_star = alpha1,
    tau1 = tau1,
    z_c = z_c_calibration,
    z_d = z_d_calibration,
    w_min = w_min,
    w_max = w_max
  )

  res_r1star <- r1(
    y = y_data[calibration_idx],
    w = w_data[calibration_idx],
    delta = delta_data[calibration_idx],
    beta = beta_hat,
    alpha1_star = alpha1_star_r,
    tau1 = tau1_r,
    z_c = z_c_calibration,
    z_d = z_d_calibration,
    w_min = w_min,
    w_max = w_max
  )

  zeta <- c(
    r1 = unname(stats::quantile(c(res_r1, Inf), 1 - alpha)),
    r2 = unname(stats::quantile(c(res_r2, Inf), 1 - alpha)),
    r1star = unname(stats::quantile(c(res_r1star, Inf), 1 - alpha))
  )

  coverage_rate <- .coverage_on_test(
    zeta = zeta,
    beta = beta_hat,
    test_y_data = test_y_data,
    test_w_data = test_w_data,
    test_delta_data = test_delta_data,
    test_z_c_data = test_z_c_data,
    test_z_d_data = test_z_d_data,
    alpha1 = alpha1,
    tau1 = tau1,
    alpha1_star_r = alpha1_star_r,
    tau1_r = tau1_r,
    w_min = w_min,
    w_max = w_max
  )

  list(
    method = "split_conformal",
    alpha = alpha,
    residual = c("r1", "r2", "r1star"),
    zeta = zeta,
    coverage_rate = coverage_rate,
    beta = beta_hat,
    sigma = sigma_hat,
    n_fit = length(fit_idx),
    n_calibrate = length(calibration_idx)
  )
}


# =============================================================================
# Full conformal prediction
# =============================================================================

#' Full conformal prediction interval
#'
#' Refits the outcome model after augmenting the training sample with each test
#' observation and evaluates conformal inclusion using the augmented residuals.
#'
#' @param y_data,w_data,delta_data Training outcome, observed time, and event
#'   indicator.
#' @param z_c_data,z_d_data Optional training covariates.
#' @param alpha Miscoverage level.
#' @param alpha1,tau1 Coefficients and standard deviation for the working model
#'   for \code{X | Z}.
#' @param alpha1_star_r,tau1_r Coefficients and SD for \code{r1*}.
#' @param test_y_data,test_w_data,test_delta_data Test outcome, observed time,
#'   and event indicator.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param w_min,w_max Truncation bounds.
#'
#' @return A list containing half-lengths and empirical coverage rates for
#'   \code{r1}, \code{r2}, and \code{r1*}, together with training and test
#'   sample sizes.
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' c_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w_vec <- pmin(x_vec, c_vec)
#' delta_vec <- as.integer(x_vec <= c_vec)
#' y_vec <- 3 * x_vec + rnorm(n)
#'
#' alpha1_fit <- unname(find_alpha1_MLE(
#'   w_vec, delta_vec, w_min = -1, w_max = 1
#' ))
#'
#' set.seed(2)
#' n_test <- 5
#' x_test_vec <- truncnorm::rtruncnorm(
#'   n_test, a = -1, b = 1, mean = 0, sd = 1
#' )
#' c_test_vec <- truncnorm::rtruncnorm(
#'   n_test, a = -1, b = 1, mean = 1, sd = 1
#' )
#' y_test_vec <- 3 * x_test_vec + rnorm(n_test)
#' w_test_vec <- pmin(x_test_vec, c_test_vec)
#' delta_test_vec <- as.integer(x_test_vec <= c_test_vec)
#'
#' fit <- full_conformal_prediction_interval(
#'   y_data = y_vec,
#'   w_data = w_vec,
#'   delta_data = delta_vec,
#'   alpha = 0.1,
#'   alpha1 = alpha1_fit[1],
#'   tau1 = alpha1_fit[2],
#'   alpha1_star_r = alpha1_fit[1],
#'   tau1_r = alpha1_fit[2],
#'   test_y_data = y_test_vec,
#'   test_w_data = w_test_vec,
#'   test_delta_data = delta_test_vec,
#'   w_min = -1,
#'   w_max = 1
#' )
#' fit$zeta
#' fit$coverage_rate
#'
#' @export
full_conformal_prediction_interval <- function(
    y_data, w_data, delta_data,
    z_c_data = NULL, z_d_data = NULL,
    alpha,
    alpha1, tau1,
    alpha1_star_r, tau1_r,
    test_y_data, test_w_data, test_delta_data,
    test_z_c_data = NULL, test_z_d_data = NULL,
    w_min, w_max) {

  n <- length(y_data)
  n_test <- length(test_y_data)

  if (length(w_data) != n || length(delta_data) != n) {
    stop("y_data, w_data, and delta_data must have the same length.")
  }

  if (
    length(test_w_data) != n_test ||
    length(test_delta_data) != n_test
  ) {
    stop(
      "test_y_data, test_w_data, and test_delta_data must have the same length."
    )
  }

  z_c_mat <- .as_covariate_matrix(
    z_c_data,
    n,
    "z_c_data"
  )

  z_d_mat <- .as_covariate_matrix(
    z_d_data,
    n,
    "z_d_data"
  )

  test_z_c_mat <- .as_covariate_matrix(
    test_z_c_data,
    n_test,
    "test_z_c_data"
  )

  test_z_d_mat <- .as_covariate_matrix(
    test_z_d_data,
    n_test,
    "test_z_d_data"
  )

  if (ncol(test_z_c_mat) != ncol(z_c_mat)) {
    stop("Training and test continuous covariates must have the same columns.")
  }

  if (ncol(test_z_d_mat) != ncol(z_d_mat)) {
    stop("Training and test discrete covariates must have the same columns.")
  }

  pred_cvg_list <- lapply(
    seq_len(n_test),
    function(j) {

      y_aug <- c(
        test_y_data[j],
        y_data
      )

      w_aug <- c(
        test_w_data[j],
        w_data
      )

      delta_aug <- c(
        test_delta_data[j],
        delta_data
      )

      z_c_aug <- rbind(
        test_z_c_mat[j, , drop = FALSE],
        z_c_mat
      )

      z_d_aug <- rbind(
        test_z_d_mat[j, , drop = FALSE],
        z_d_mat
      )

      cc_idx <- delta_aug == 1

      if (!any(cc_idx)) {
        stop("No complete cases in the augmented sample.")
      }

      x_cc <- w_aug[cc_idx]

      X <- phi_xz(
        x_cc,
        z_c = z_c_aug[cc_idx, , drop = FALSE],
        z_d = z_d_aug[cc_idx, , drop = FALSE]
      )

      beta_cc <- as.numeric(
        solve(
          crossprod(X),
          crossprod(X, y_aug[cc_idx])
        )
      )

      res_r1 <- r1(
        y = y_aug,
        w = w_aug,
        delta = delta_aug,
        beta = beta_cc,
        alpha1_star = alpha1,
        tau1 = tau1,
        z_c = z_c_aug,
        z_d = z_d_aug,
        w_min = w_min,
        w_max = w_max
      )

      res_r2 <- r2(
        y = y_aug,
        w = w_aug,
        delta = delta_aug,
        beta = beta_cc,
        alpha1_star = alpha1,
        tau1 = tau1,
        z_c = z_c_aug,
        z_d = z_d_aug,
        w_min = w_min,
        w_max = w_max
      )

      res_r1star <- r1(
        y = y_aug,
        w = w_aug,
        delta = delta_aug,
        beta = beta_cc,
        alpha1_star = alpha1_star_r,
        tau1 = tau1_r,
        z_c = z_c_aug,
        z_d = z_d_aug,
        w_min = w_min,
        w_max = w_max
      )

      c(
        cvg_r1 = as.numeric(
          res_r1[1] <= stats::quantile(res_r1, 1 - alpha)
        ),
        cvg_r2 = as.numeric(
          res_r2[1] <= stats::quantile(res_r2, 1 - alpha)
        ),
        cvg_r1star = as.numeric(
          res_r1star[1] <= stats::quantile(res_r1star, 1 - alpha)
        ),
        zeta_r1 = unname(
          stats::quantile(res_r1, 1 - alpha)
        ),
        zeta_r2 = unname(
          stats::quantile(res_r2, 1 - alpha)
        ),
        zeta_r1star = unname(
          stats::quantile(res_r1star, 1 - alpha)
        )
      )
    }
  )

  pred_cvg_mat <- do.call(
    rbind,
    pred_cvg_list
  )

  avg <- colMeans(
    pred_cvg_mat
  )

  list(
    method = "full_conformal",
    alpha = alpha,
    residual = c("r1", "r2", "r1star"),
    zeta = c(
      r1 = unname(avg[["zeta_r1"]]),
      r2 = unname(avg[["zeta_r2"]]),
      r1star = unname(avg[["zeta_r1star"]])
    ),
    coverage_rate = c(
      r1 = unname(avg[["cvg_r1"]]),
      r2 = unname(avg[["cvg_r2"]]),
      r1star = unname(avg[["cvg_r1star"]])
    ),
    n_train = n,
    n_test = n_test
  )
}


# =============================================================================
# Jackknife+
# =============================================================================

#' Jackknife+ prediction interval
#'
#' Fits the outcome model after leaving out each training observation and forms
#' jackknife+ prediction intervals from the resulting residuals and predictions.
#'
#' @param y_data,w_data,delta_data Training outcome, observed time, and event
#'   indicator.
#' @param z_c_data,z_d_data Optional training covariates.
#' @param alpha Miscoverage level.
#' @param alpha1,tau1 Coefficients and standard deviation for the working model
#'   for \code{X | Z}.
#' @param alpha1_star_r,tau1_r Coefficients and SD for \code{r1*}.
#' @param test_y_data,test_w_data,test_delta_data Test outcome, observed time,
#'   and event indicator.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param w_min,w_max Truncation bounds.
#'
#' @return A list containing half-lengths and empirical coverage rates for
#'   \code{r1}, \code{r2}, and \code{r1*}, together with training and test
#'   sample sizes.
#'
#' @examples
#' set.seed(1)
#' n <- 40
#' x_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' c_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w_vec <- pmin(x_vec, c_vec)
#' delta_vec <- as.integer(x_vec <= c_vec)
#' y_vec <- 3 * x_vec + rnorm(n)
#'
#' alpha1_fit <- unname(find_alpha1_MLE(
#'   w_vec, delta_vec, w_min = -1, w_max = 1
#' ))
#'
#' set.seed(2)
#' n_test <- 5
#' x_test_vec <- truncnorm::rtruncnorm(
#'   n_test, a = -1, b = 1, mean = 0, sd = 1
#' )
#' c_test_vec <- truncnorm::rtruncnorm(
#'   n_test, a = -1, b = 1, mean = 1, sd = 1
#' )
#' y_test_vec <- 3 * x_test_vec + rnorm(n_test)
#' w_test_vec <- pmin(x_test_vec, c_test_vec)
#' delta_test_vec <- as.integer(x_test_vec <= c_test_vec)
#'
#' fit <- jackknife_plus_prediction_interval(
#'   y_data = y_vec,
#'   w_data = w_vec,
#'   delta_data = delta_vec,
#'   alpha = 0.1,
#'   alpha1 = alpha1_fit[1],
#'   tau1 = alpha1_fit[2],
#'   alpha1_star_r = alpha1_fit[1],
#'   tau1_r = alpha1_fit[2],
#'   test_y_data = y_test_vec,
#'   test_w_data = w_test_vec,
#'   test_delta_data = delta_test_vec,
#'   w_min = -1,
#'   w_max = 1
#' )
#' fit$zeta
#' fit$coverage_rate
#'
#' @export
jackknife_plus_prediction_interval <- function(
    y_data, w_data, delta_data,
    z_c_data = NULL, z_d_data = NULL,
    alpha,
    alpha1, tau1,
    alpha1_star_r, tau1_r,
    test_y_data, test_w_data, test_delta_data,
    test_z_c_data = NULL, test_z_d_data = NULL,
    w_min, w_max) {

  n <- length(y_data)
  n_test <- length(test_y_data)

  if (length(w_data) != n || length(delta_data) != n) {
    stop("y_data, w_data, and delta_data must have the same length.")
  }

  if (
    length(test_w_data) != n_test ||
    length(test_delta_data) != n_test
  ) {
    stop(
      "test_y_data, test_w_data, and test_delta_data must have the same length."
    )
  }

  z_c_mat <- .as_covariate_matrix(
    z_c_data,
    n,
    "z_c_data"
  )

  z_d_mat <- .as_covariate_matrix(
    z_d_data,
    n,
    "z_d_data"
  )

  test_z_c_mat <- .as_covariate_matrix(
    test_z_c_data,
    n_test,
    "test_z_c_data"
  )

  test_z_d_mat <- .as_covariate_matrix(
    test_z_d_data,
    n_test,
    "test_z_d_data"
  )

  if (ncol(test_z_c_mat) != ncol(z_c_mat)) {
    stop("Training and test continuous covariates must have the same columns.")
  }

  if (ncol(test_z_d_mat) != ncol(z_d_mat)) {
    stop("Training and test discrete covariates must have the same columns.")
  }

  m_list <- lapply(
    seq_len(n),
    function(i) {

      cc_idx <- delta_data[-i] == 1

      if (!any(cc_idx)) {
        stop("No complete cases after leaving out a training observation.")
      }

      x_cc <- w_data[-i][cc_idx]

      X <- phi_xz(
        x_cc,
        z_c = z_c_mat[-i, , drop = FALSE][cc_idx, , drop = FALSE],
        z_d = z_d_mat[-i, , drop = FALSE][cc_idx, , drop = FALSE]
      )

      beta_cc <- as.numeric(
        solve(
          crossprod(X),
          crossprod(X, y_data[-i][cc_idx])
        )
      )

      m1_test <- m1(
        w = test_w_data,
        delta = test_delta_data,
        beta = beta_cc,
        alpha1_star = alpha1,
        tau1 = tau1,
        z_c = test_z_c_mat,
        z_d = test_z_d_mat,
        w_min = w_min,
        w_max = w_max
      )

      m2_test <- m0(
        x = test_w_data,
        beta = beta_cc,
        z_c = test_z_c_mat,
        z_d = test_z_d_mat
      )

      m1star_test <- m1(
        w = test_w_data,
        delta = test_delta_data,
        beta = beta_cc,
        alpha1_star = alpha1_star_r,
        tau1 = tau1_r,
        z_c = test_z_c_mat,
        z_d = test_z_d_mat,
        w_min = w_min,
        w_max = w_max
      )

      z_c_i <- if (ncol(z_c_mat) > 0) {
        z_c_mat[i, ]
      } else {
        NULL
      }

      z_d_i <- if (ncol(z_d_mat) > 0) {
        z_d_mat[i, ]
      } else {
        NULL
      }

      r1_i <- r1(
        y = y_data[i],
        w = w_data[i],
        delta = delta_data[i],
        beta = beta_cc,
        alpha1_star = alpha1,
        tau1 = tau1,
        z_c = z_c_i,
        z_d = z_d_i,
        w_min = w_min,
        w_max = w_max
      )

      r2_i <- r2(
        y = y_data[i],
        w = w_data[i],
        delta = delta_data[i],
        beta = beta_cc,
        alpha1_star = alpha1,
        tau1 = tau1,
        z_c = z_c_i,
        z_d = z_d_i,
        w_min = w_min,
        w_max = w_max
      )

      r1star_i <- r1(
        y = y_data[i],
        w = w_data[i],
        delta = delta_data[i],
        beta = beta_cc,
        alpha1_star = alpha1_star_r,
        tau1 = tau1_r,
        z_c = z_c_i,
        z_d = z_d_i,
        w_min = w_min,
        w_max = w_max
      )

      cbind(
        m1_test - r1_i,
        m1_test + r1_i,
        m2_test - r2_i,
        m2_test + r2_i,
        m1star_test - r1star_i,
        m1star_test + r1star_i
      )
    }
  )

  m_vals <- simplify2array(
    m_list
  )

  lower <- apply(
    m_vals[, c(1, 3, 5), , drop = FALSE],
    c(1, 2),
    function(x) {
      stats::quantile(
        x,
        probs = alpha
      )
    }
  )

  upper <- apply(
    m_vals[, c(2, 4, 6), , drop = FALSE],
    c(1, 2),
    function(x) {
      stats::quantile(
        x,
        probs = 1 - alpha
      )
    }
  )

  coverage_rate <- vapply(
    seq_len(3),
    function(j) {
      mean(
        test_y_data >= lower[, j] &
          test_y_data <= upper[, j]
      )
    },
    numeric(1)
  )

  names(coverage_rate) <- c(
    "r1",
    "r2",
    "r1star"
  )

  zeta <- colMeans(
    (upper - lower) / 2
  )

  names(zeta) <- c(
    "r1",
    "r2",
    "r1star"
  )

  list(
    method = "jackknife_plus",
    alpha = alpha,
    residual = c("r1", "r2", "r1star"),
    zeta = zeta,
    coverage_rate = coverage_rate,
    n_train = n,
    n_test = n_test
  )
}
