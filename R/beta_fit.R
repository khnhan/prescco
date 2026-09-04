# =============================================================================
# Coefficient estimation: complete-case and SPARCC
# =============================================================================

#' Complete-case fit of the outcome model
#'
#' Estimates the outcome-model coefficient vector using only uncensored
#' observations (\code{delta == 1}), for which \code{X = W} is observed.
#' For the linear outcome model implemented here, the complete-case estimator
#' is obtained by least squares
#'
#' @param y_data,w_data,delta_data Outcome, observed time, and event indicator
#'   (length \code{n}).
#' @param z_c_data,z_d_data Optional continuous and discrete covariates.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{beta_cc}: estimated coefficient vector;
#'   \item \code{sigma_cc}: residual standard deviation;
#'   \item \code{design_cc}: complete-case design matrix;
#'   \item \code{idx_cc}: indices of the complete cases;
#'   \item \code{p_c}, \code{p_d}: numbers of continuous and discrete
#'     covariates.
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' c_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w_vec <- pmin(x_vec, c_vec)
#' delta_vec <- as.integer(x_vec <= c_vec)
#' y_vec <- 3 * x_vec + rnorm(n)
#'
#' fit_cc <- find_beta_cc(y_vec, w_vec, delta_vec)
#' fit_cc$beta_cc
#' fit_cc$sigma_cc
#' length(fit_cc$idx_cc)
#'
#' # Under the linear outcome model, this agrees with a least-squares fit
#' # using the uncensored observations.
#' stats::coef(stats::lm(
#'   y_vec[delta_vec == 1] ~ w_vec[delta_vec == 1]
#' ))
#'
#' @export
find_beta_cc <- function(y_data, w_data, delta_data,
                         z_c_data = NULL, z_d_data = NULL) {

  n <- length(y_data)

  if (length(w_data) != n || length(delta_data) != n) {
    stop("y_data, w_data, and delta_data must have the same length.")
  }

  # Continuous covariates.
  if (is.null(z_c_data)) {

    z_c_mat <- matrix(0, nrow = n, ncol = 0)
    p_c <- 0

  } else if (is.null(dim(z_c_data))) {

    if (length(z_c_data) != n) {
      stop("z_c_data must have length n.")
    }

    z_c_mat <- matrix(as.numeric(z_c_data), ncol = 1)
    p_c <- 1

  } else {

    if (nrow(z_c_data) != n) {
      stop("z_c_data must have n rows.")
    }

    z_c_mat <- as.matrix(z_c_data)
    p_c <- ncol(z_c_mat)
  }

  # Discrete covariates.
  if (is.null(z_d_data)) {

    z_d_mat <- matrix(0, nrow = n, ncol = 0)
    p_d <- 0

  } else if (is.null(dim(z_d_data))) {

    if (length(z_d_data) != n) {
      stop("z_d_data must have length n.")
    }

    z_d_mat <- matrix(z_d_data, ncol = 1)
    p_d <- 1

  } else {

    if (nrow(z_d_data) != n) {
      stop("z_d_data must have n rows.")
    }

    z_d_mat <- as.matrix(z_d_data)
    p_d <- ncol(z_d_mat)
  }

  # Complete cases: X = W when delta = 1.
  idx_cc <- which(delta_data == 1)

  if (length(idx_cc) == 0) {
    stop("No complete cases (delta_data == 1); cannot compute CC estimator.")
  }

  y_cc <- y_data[idx_cc]
  x_cc <- w_data[idx_cc]

  z_c_cc <- z_c_mat[idx_cc, , drop = FALSE]
  z_d_cc <- z_d_mat[idx_cc, , drop = FALSE]

  # Design matrix for Y ~ 1 + X + Z_c + Z_d + X:Z_c + X:Z_d.
  design_cc <- cbind(
    1,
    x_cc,
    if (p_c > 0) z_c_cc else NULL,
    if (p_d > 0) z_d_cc else NULL,
    if (p_c > 0) sweep(z_c_cc, 1, x_cc, `*`) else NULL,
    if (p_d > 0) sweep(z_d_cc, 1, x_cc, `*`) else NULL
  )

  d_beta <- ncol(design_cc)

  fit_cc <- stats::lm.fit(design_cc, y_cc)
  beta_cc <- as.numeric(fit_cc$coefficients)

  if (length(beta_cc) < d_beta) {
    beta_cc <- c(beta_cc, rep(NA_real_, d_beta - length(beta_cc)))
  }

  beta_cc[is.na(beta_cc)] <- 0

  resid_cc <- y_cc - as.vector(design_cc %*% beta_cc)
  sigma_cc <- sd(resid_cc)

  if (!is.finite(sigma_cc) || sigma_cc <= 0) {
    warning("CC residual sd is not positive/finite; using IQR/1.349 as fallback.")
    sigma_cc <- IQR(resid_cc) / 1.349
  }

  list(
    beta_cc = beta_cc,
    sigma_cc = sigma_cc,
    design_cc = design_cc,
    idx_cc = idx_cc,
    p_c = p_c,
    p_d = p_d
  )
}


#' SPARCC fit of the outcome model
#'
#' Solves the efficient estimating equation for beta under working nuisance
#' models for \code{X | Z} and \code{C | Z}.
#'
#' @param y_data,w_data,delta_data Outcome, observed time, and event indicator.
#' @param z_c_data,z_d_data Optional continuous and discrete covariates.
#' @param alpha1_star,alpha2_star Coefficient vectors for the working
#'   \code{X | Z} and \code{C | Z} models, typically obtained from
#'   \code{\link{find_alpha1_MLE}} and \code{\link{find_alpha2_MLE}}.
#' @param tau1,tau2 Standard deviations of \code{X} and \code{C}.
#' @param tt,m Integration grid size and covariate-grid resolution.
#' @param w_min,w_max Truncation bounds for \code{X} and \code{C}.
#' @param verbose If \code{TRUE}, report the \code{nleqslv} convergence
#'   message. Defaults to \code{FALSE}.
#'
#' @return A list containing the SPARCC estimates \code{beta_hat} and
#'   \code{sigma_hat}, their starting values, the combined parameter vectors,
#'   and the \code{nleqslv} result.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' c_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w_vec <- pmin(x_vec, c_vec)
#' delta_vec <- as.integer(x_vec <= c_vec)
#' y_vec <- 3 * x_vec + rnorm(n)
#'
#' # Fit the two working nuisance models.
#' a1 <- unname(find_alpha1_MLE(
#'   w_vec, delta_vec, w_min = -1, w_max = 1
#' ))
#' a2 <- unname(find_alpha2_MLE(
#'   w_vec, delta_vec, w_min = -1, w_max = 1
#' ))
#'
#' fit_sparcc <- find_beta_sparcc(
#'   y_vec, w_vec, delta_vec,
#'   z_c_data = NULL, z_d_data = NULL,
#'   alpha1_star = a1[1], alpha2_star = a2[1],
#'   tau1 = a1[2], tau2 = a2[2],
#'   w_min = -1, w_max = 1
#' )
#' fit_sparcc$beta_hat
#' fit_sparcc$sigma_hat
#'
#' @export
find_beta_sparcc <- function(y_data, w_data, delta_data,
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

  # Continuous covariates.
  if (is.null(z_c_data)) {

    z_c_mat <- matrix(0, nrow = n, ncol = 0)
    p_c <- 0

  } else if (is.null(dim(z_c_data))) {

    if (length(z_c_data) != n) {
      stop("z_c_data must have length n")
    }

    z_c_mat <- matrix(as.numeric(z_c_data), ncol = 1)
    p_c <- 1

  } else {

    if (nrow(z_c_data) != n) {
      stop("z_c_data must have n rows")
    }

    z_c_mat <- as.matrix(z_c_data)
    p_c <- ncol(z_c_mat)
  }

  # Discrete covariates.
  if (is.null(z_d_data)) {

    z_d_mat <- matrix(0, nrow = n, ncol = 0)
    p_d <- 0

  } else if (is.null(dim(z_d_data))) {

    if (length(z_d_data) != n) {
      stop("z_d_data must have length n")
    }

    z_d_mat <- matrix(z_d_data, ncol = 1)
    p_d <- 1

  } else {

    if (nrow(z_d_data) != n) {
      stop("z_d_data must have n rows")
    }

    z_d_mat <- as.matrix(z_d_data)
    p_d <- ncol(z_d_mat)
  }

  # Model: Y ~ 1 + X + Z_c + Z_d + X:Z_c + X:Z_d.
  d_beta <- 2 + 2 * p_c + 2 * p_d

  # Complete-case starting values.
  idx_cc <- which(delta_data == 1)

  if (length(idx_cc) == 0) {

    warning("No complete cases (delta = 1); using zeros for beta start.")
    start_beta <- rep(0, d_beta)
    start_sigma <- (w_max - w_min) / 4

  } else {

    if (length(idx_cc) < d_beta) {
      warning("Very few complete cases relative to beta dimension; CC start may be unstable.")
    }

    y_cc <- y_data[idx_cc]
    x_cc <- w_data[idx_cc]

    z_c_cc <- z_c_mat[idx_cc, , drop = FALSE]
    z_d_cc <- z_d_mat[idx_cc, , drop = FALSE]

    design_cc <- cbind(
      1,
      x_cc,
      if (p_c > 0) z_c_cc else NULL,
      if (p_d > 0) z_d_cc else NULL,
      if (p_c > 0) sweep(z_c_cc, 1, x_cc, `*`) else NULL,
      if (p_d > 0) sweep(z_d_cc, 1, x_cc, `*`) else NULL
    )

    cc_fit <- stats::lm.fit(design_cc, y_cc)
    start_beta <- as.numeric(cc_fit$coefficients)

    if (length(start_beta) < ncol(design_cc)) {
      start_beta <- c(
        start_beta,
        rep(NA_real_, ncol(design_cc) - length(start_beta))
      )
    }

    start_beta[is.na(start_beta)] <- 0

    if (length(start_beta) < d_beta) {
      start_beta <- c(start_beta, rep(0, d_beta - length(start_beta)))
    } else if (length(start_beta) > d_beta) {
      start_beta <- start_beta[seq_len(d_beta)]
    }

    resid_cc <- y_cc - as.vector(design_cc %*% start_beta)
    start_sigma <- sd(resid_cc)

    if (!is.finite(start_sigma) || start_sigma <= 0) {
      start_sigma <- (w_max - w_min) / 4
    }
  }

  # Grid for X.
  x_a <- seq(w_min, w_max, length.out = m)

  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("Package 'nleqslv' is required.")
  }

  # Solve for beta and log(sigma).
  start_betasigma <- c(start_beta, log(start_sigma))

  result <- nleqslv::nleqslv(
    start_betasigma,
    pe_gauss_param12,
    y_data = y_data,
    w_data = w_data,
    delta_data = delta_data,
    x_a = x_a,
    z_c_data = z_c_data,
    z_d_data = z_d_data,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    tau1 = tau1,
    tau2 = tau2,
    w_min = w_min,
    w_max = w_max,
    tt = tt
  )

  if (verbose) {
    message("nleqslv: ", result$message)
  }

  betasigma_hat <- result$x
  beta_hat <- betasigma_hat[seq_len(d_beta)]
  sigma_hat <- exp(betasigma_hat[d_beta + 1])

  list(
    beta_hat = beta_hat,
    sigma_hat = sigma_hat,
    beta_start = start_beta,
    sigma_start = start_sigma,
    betasigma_hat = betasigma_hat,
    betasigma_start = start_betasigma,
    nleqslv_result = result
  )
}
