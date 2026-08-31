# =============================================================================
#  Coefficient estimation: complete-case and SPARCC
#  ---------------------------------------------------------------------------
#  Two estimators of the outcome-model coefficient vector beta (and residual
#  scale sigma) from right-censored-covariate data, both covariate-aware:
#
#    * find_beta_cc()               complete-case estimtaor on the
#                                   uncensored (delta == 1) observations.
#    * find_beta_sparcc()           SPARCC estimator, solving the efficient estimating
#                                   equation under working nuisance models
#                                   (alpha1_star, alpha2_star, tau1, tau2).
#
#  Both accept optional continuous (z_c) and discrete (z_d) covariates.
#
#  Package imports are declared centrally in R/prescco-package.R.
# =============================================================================

#' Complete-case fit of the outcome model
#'
#' Fits beta by ordinary least squares on the uncensored observations
#' (delta == 1), where X = W is observed, using the phi_xz design.
#'
#' @param y_data,w_data,delta_data Outcome, observed time, and event indicator
#'   (length n).
#' @param z_c_data,z_d_data Optional continuous / discrete covariates.
#' @return A list with `beta_cc` (coefficient vector), `sigma_cc` (residual
#'   standard deviation), `idx_cc` (indices of the complete cases used), and
#'   `p_c` / `p_d` (number of continuous / discrete covariates).
#' @examples
#' set.seed(1)
#' n  <- 200
#' x  <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' cc <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w  <- pmin(x, cc)
#' delta <- as.integer(x <= cc)
#' y  <- 3 * x + rnorm(n)          # beta = c(0, 3)
#'
#' fit <- find_beta_cc(y, w, delta)
#' fit$beta_cc          # compare with the true c(0, 3)
#' fit$sigma_cc
#' length(fit$idx_cc)   # how many observations were actually used
#'
#' # Equivalent to least squares on the uncensored rows only
#' coef(lm(y[delta == 1] ~ w[delta == 1]))
#' @export
find_beta_cc <- function(y_data, w_data, delta_data,
                         z_c_data = NULL, z_d_data = NULL) {
  n <- length(y_data)
  if (length(w_data) != n || length(delta_data) != n) {
    stop("y_data, w_data, and delta_data must have the same length.")
  }

  ## Handle z_c_data and z_d_data as matrices (possibly 0 cols)
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
    z_d_mat <- matrix(z_d_data, ncol = 1)
    p_d     <- 1
  } else {
    if (nrow(z_d_data) != n) stop("z_d_data must have n rows.")
    z_d_mat <- as.matrix(z_d_data)
    p_d     <- ncol(z_d_mat)
  }

  ## Complete cases: X = W when delta = 1
  idx_cc <- which(delta_data == 1)
  if (length(idx_cc) == 0) {
    stop("No complete cases (delta_data == 1); cannot compute CC estimator.")
  }

  y_cc <- y_data[idx_cc]
  x_cc <- w_data[idx_cc]

  z_c_cc <- z_c_mat[idx_cc, , drop = FALSE]
  z_d_cc <- z_d_mat[idx_cc, , drop = FALSE]

  ## Design matrix for:
  ## Y ~ 1 + X + Z_c + Z_d + X:Z_c + X:Z_d
  design_cc <- cbind(
    1,
    x_cc,
    if (p_c > 0) z_c_cc else NULL,
    if (p_d > 0) z_d_cc else NULL,
    if (p_c > 0) sweep(z_c_cc, 1, x_cc, `*`) else NULL,
    if (p_d > 0) sweep(z_d_cc, 1, x_cc, `*`) else NULL
  )

  d_beta <- ncol(design_cc)

  ## OLS on complete cases
  fit_cc   <- stats::lm.fit(design_cc, y_cc)
  beta_cc  <- as.numeric(fit_cc$coefficients)

  ## Replace NA coefficients (if any dropped) with 0 and ensure length d_beta
  if (length(beta_cc) < d_beta) {
    beta_cc <- c(beta_cc, rep(NA_real_, d_beta - length(beta_cc)))
  }
  beta_cc[is.na(beta_cc)] <- 0

  ## Residual-based sigma
  resid_cc  <- y_cc - as.vector(design_cc %*% beta_cc)
  sigma_cc  <- sd(resid_cc)
  if (!is.finite(sigma_cc) || sigma_cc <= 0) {
    warning("CC residual sd is not positive/finite; using IQR/1.349 as fallback.")
    sigma_cc <- IQR(resid_cc) / 1.349
  }

  list(
    beta_cc    = beta_cc,
    sigma_cc   = sigma_cc,
    design_cc  = design_cc,
    idx_cc     = idx_cc,
    p_c        = p_c,
    p_d        = p_d
  )
}

#' SPARCC fit of the outcome model
#'
#' Solves the efficient estimating equation for beta under the working nuisance
#' models for X | Z and C | Z. Covariate-aware via the phi_xz basis.
#'
#' @param y_data,w_data,delta_data Outcome, observed time, event indicator.
#' @param z_c_data,z_d_data Optional continuous / discrete covariates.
#' @param alpha1_star,alpha2_star Coefficient vectors for the eta1 / eta2 means
#'   (typically from find_alpha1_MLE / find_alpha2_MLE).
#' @param tau1,tau2 Standard deviations of X and C.
#' @param tt,m Integration grid size and covariate-grid resolution.
#' @param w_min,w_max Truncation bounds for X and C.
#' @param verbose If \code{TRUE}, report the \code{nleqslv} convergence
#'   message. Defaults to \code{FALSE}.
#' @return A list with `beta_hat` and `sigma_hat`.
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
#' # Fit the two working models first
#' a1 <- unname(find_alpha1_MLE(w, delta, w_min = -1, w_max = 1))
#' a2 <- unname(find_alpha2_MLE(w, delta, w_min = -1, w_max = 1))
#'
#' fit <- find_beta_sparcc(y, w, delta, NULL, NULL,
#'                         alpha1_star = a1[1], alpha2_star = a2[1],
#'                         tau1 = a1[2], tau2 = a2[2],
#'                         w_min = -1, w_max = 1)
#' fit$beta_hat     # uses every observation, unlike find_beta_cc
#' fit$sigma_hat
#' }
#' @export
find_beta_sparcc = function(y_data, w_data, delta_data,
                                     z_c_data, z_d_data,
                                     alpha1_star, alpha2_star,
                                     tau1, tau2,
                                     tt = 20, m = 20,
                                     w_min, w_max,
                                     verbose = FALSE) {

  n <- length(y_data)
  if (length(w_data) != n || length(delta_data) != n) {
    stop("y_data, w_data, and delta_data must have the same length.")
  }

  ## ---------- Handle z_c_data and z_d_data as matrices ----------
  if (is.null(z_c_data)) {
    z_c_mat <- matrix(0, nrow = n, ncol = 0)
    p_c     <- 0
  } else if (is.null(dim(z_c_data))) {
    if (length(z_c_data) != n) stop("z_c_data must have length n")
    z_c_mat <- matrix(as.numeric(z_c_data), ncol = 1)
    p_c     <- 1
  } else {
    if (nrow(z_c_data) != n) stop("z_c_data must have n rows")
    z_c_mat <- as.matrix(z_c_data)
    p_c     <- ncol(z_c_mat)
  }

  if (is.null(z_d_data)) {
    z_d_mat <- matrix(0, nrow = n, ncol = 0)
    p_d     <- 0
  } else if (is.null(dim(z_d_data))) {
    if (length(z_d_data) != n) stop("z_d_data must have length n")
    z_d_mat <- matrix(z_d_data, ncol = 1)
    p_d     <- 1
  } else {
    if (nrow(z_d_data) != n) stop("z_d_data must have n rows")
    z_d_mat <- as.matrix(z_d_data)
    p_d     <- ncol(z_d_mat)
  }

  ## ---------- Beta dimension implied by design ----------
  ## Model: Y ~ 1 + X + Z_c + Z_d + X:Z_c + X:Z_d
  d_beta <- 2 + 2 * p_c + 2 * p_d   # (1, X, Z_c, Z_d, XZ_c, XZ_d)

  ## ---------- Complete cases: X = W when delta = 1 ----------
  idx_cc <- which(delta_data == 1)

  if (length(idx_cc) == 0) {
    warning("No complete cases (delta = 1); using zeros for beta start.")
    start_beta  <- rep(0, d_beta)
    start_sigma <- (w_max - w_min) / 4
  } else {
    if (length(idx_cc) < d_beta) {
      warning("Very few complete cases relative to beta dimension; CC start may be unstable.")
    }

    y_cc <- y_data[idx_cc]
    x_cc <- w_data[idx_cc]

    z_c_cc <- z_c_mat[idx_cc, , drop = FALSE]
    z_d_cc <- z_d_mat[idx_cc, , drop = FALSE]

    ## Design for: Y ~ 1 + X + Z_c + Z_d + X:Z_c + X:Z_d
    design_cc <- cbind(
      1,
      x_cc,
      if (p_c > 0) z_c_cc else NULL,
      if (p_d > 0) z_d_cc else NULL,
      if (p_c > 0) sweep(z_c_cc, 1, x_cc, `*`) else NULL,
      if (p_d > 0) sweep(z_d_cc, 1, x_cc, `*`) else NULL
    )

    ## OLS on complete cases
    cc_fit     <- stats::lm.fit(design_cc, y_cc)
    start_beta <- as.numeric(cc_fit$coefficients)

    ## Ensure correct length and replace NAs with 0
    if (length(start_beta) < ncol(design_cc)) {
      start_beta <- c(start_beta, rep(NA_real_, ncol(design_cc) - length(start_beta)))
    }
    start_beta[is.na(start_beta)] <- 0

    if (length(start_beta) < d_beta) {
      start_beta <- c(start_beta, rep(0, d_beta - length(start_beta)))
    } else if (length(start_beta) > d_beta) {
      start_beta <- start_beta[seq_len(d_beta)]
    }

    ## Residual-based sigma from CC fit
    resid_cc    <- y_cc - as.vector(design_cc %*% start_beta)
    start_sigma <- sd(resid_cc)
    if (!is.finite(start_sigma) || start_sigma <= 0) {
      start_sigma <- (w_max - w_min) / 4
    }
  }

  ## ---------- X-grid ----------
  x_a <- seq(w_min, w_max, length.out = m)

  ## ---------- Equation solve for (beta, log(sigma)) via pe_gauss_param12 ----------
  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("Package 'nleqslv' is required.")
  }

  # Parameter for nleqslv: betasigma = (beta, log_sigma)
  start_betasigma <- c(start_beta, log(start_sigma))

  result <- nleqslv::nleqslv(
    start_betasigma,
    pe_gauss_param12,
    y_data      = y_data,
    w_data      = w_data,
    delta_data  = delta_data,
    x_a         = x_a,
    z_c_data    = z_c_data,
    z_d_data    = z_d_data,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    tau1        = tau1,
    tau2        = tau2,
    w_min       = w_min,
    w_max       = w_max,
    tt          = tt
  )

  if (verbose) message("nleqslv: ", result$message)

  ## ---------- Extract estimates ----------
  betasigma_hat <- result$x
  beta_hat      <- betasigma_hat[1:d_beta]
  sigma_hat     <- exp(betasigma_hat[d_beta + 1])

  list(
    beta_hat       = beta_hat,
    sigma_hat      = sigma_hat,
    beta_start     = start_beta,
    sigma_start    = start_sigma,
    betasigma_hat  = betasigma_hat,
    betasigma_start = start_betasigma,
    nleqslv_result = result
  )
}
