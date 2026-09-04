# =============================================================================
# Prediction helpers: covariate basis, interval centers, and residuals
#
# The outcome model uses the basis
#
#   phi_xz(x, z_c, z_d) = (1, x, z_c, z_d, x * z_c, x * z_d).
#
# Continuous and discrete covariates are optional. The center and residual
# functions accept either one observation or vectors of observations.
# =============================================================================


# Convert optional covariates to an n-row matrix.

.as_covariate_matrix <- function(z, n, name) {

  if (is.null(z)) {
    return(
      matrix(0, nrow = n, ncol = 0)
    )
  }

  if (is.null(dim(z))) {

    if (n == 1) {
      return(
        matrix(as.numeric(z), nrow = 1)
      )
    }

    if (length(z) != n) {
      stop(name, " must have length n or n rows.")
    }

    return(
      matrix(as.numeric(z), ncol = 1)
    )
  }

  if (nrow(z) != n) {
    stop(name, " must have n rows.")
  }

  as.matrix(z)
}


#' Covariate basis for the outcome model
#'
#' Constructs the basis
#' \code{(1, X, Z_c, Z_d, X Z_c, X Z_d)} used by the outcome model.
#'
#' @param x Event-time covariate value or numeric vector.
#' @param z_c Optional continuous covariates. For multiple observations, supply
#'   a vector for one covariate or a matrix with one row per observation.
#' @param z_d Optional discrete covariates, supplied in the same form as
#'   \code{z_c}.
#'
#' @return A numeric basis vector for one observation or a design matrix for
#'   multiple observations.
#'
#' @examples
#' phi_xz(0.5)
#' phi_xz(0.5, z_c = 2, z_d = 1)
#'
#' x_vec <- c(0.2, 0.5, 0.8)
#' z_c_vec <- c(1, 2, 3)
#' phi_xz(x_vec, z_c = z_c_vec)
#'
#' @export
phi_xz <- function(x, z_c = NULL, z_d = NULL) {

  x <- as.numeric(x)
  n <- length(x)

  z_c_mat <- .as_covariate_matrix(
    z_c,
    n,
    "z_c"
  )

  z_d_mat <- .as_covariate_matrix(
    z_d,
    n,
    "z_d"
  )

  design <- cbind(
    1,
    x,
    z_c_mat,
    z_d_mat,
    if (ncol(z_c_mat) > 0) sweep(z_c_mat, 1, x, `*`) else NULL,
    if (ncol(z_d_mat) > 0) sweep(z_d_mat, 1, x, `*`) else NULL
  )

  if (n == 1) {
    return(
      as.numeric(design[1, ])
    )
  }

  design
}


# Interval center m0.

m0 <- function(x, beta, z_c = NULL, z_d = NULL) {

  design <- phi_xz(
    x,
    z_c = z_c,
    z_d = z_d
  )

  if (is.null(dim(design))) {

    if (length(beta) != length(design)) {
      stop("beta has the wrong length for the supplied covariates.")
    }

    return(
      sum(design * beta)
    )
  }

  if (length(beta) != ncol(design)) {
    stop("beta has the wrong length for the supplied covariates.")
  }

  as.vector(
    design %*% beta
  )
}


# Interval center m1.

m1 <- function(w, delta, beta,
               alpha1_star, tau1,
               z_c = NULL, z_d = NULL,
               w_min, w_max) {

  w <- as.numeric(w)
  delta <- as.integer(delta)

  n <- length(w)

  if (length(delta) != n) {
    stop("w and delta must have the same length.")
  }

  z_c_mat <- .as_covariate_matrix(
    z_c,
    n,
    "z_c"
  )

  z_d_mat <- .as_covariate_matrix(
    z_d,
    n,
    "z_d"
  )

  alpha1_star <- as.numeric(alpha1_star)

  if (
    length(alpha1_star) !=
    1 + ncol(z_c_mat) + ncol(z_d_mat)
  ) {
    stop(
      "alpha1_star has the wrong length for the supplied covariates."
    )
  }

  mu_x <- alpha1_star[1] +
    (if (ncol(z_c_mat) > 0) {
      as.vector(
        z_c_mat %*%
          alpha1_star[
            1 + seq_len(ncol(z_c_mat))
          ]
      )
    } else {
      0
    }) +
    (if (ncol(z_d_mat) > 0) {
      offset <- 1 + ncol(z_c_mat)

      as.vector(
        z_d_mat %*%
          alpha1_star[
            offset + seq_len(ncol(z_d_mat))
          ]
      )
    } else {
      0
    })

  mu_x <- rep_len(
    mu_x,
    n
  )

  x_use <- w

  censored <- delta == 0L & w < w_max

  if (any(censored)) {

    lower <- pmax(
      w[censored],
      w_min
    )

    a <- (
      lower -
        mu_x[censored]
    ) / tau1

    b <- (
      w_max -
        mu_x[censored]
    ) / tau1

    denom <- pnorm(b) -
      pnorm(a)

    numer <- dnorm(a) -
      dnorm(b)

    good <- denom > 1e-12

    x_exp <- lower

    x_exp[good] <-
      mu_x[censored][good] +
      tau1 *
      numer[good] /
      denom[good]

    x_use[censored] <- x_exp
  }

  m0(
    x_use,
    beta,
    z_c = z_c_mat,
    z_d = z_d_mat
  )
}


#' Residual r1
#'
#' Computes \code{|Y - m1|}. For censored observations, \code{m1} averages
#' over the unobserved values of \code{X} under the working model for
#' \code{X | Z}. Supplying alternative values of \code{alpha1_star} and
#' \code{tau1} gives \code{r1*}.
#'
#' @param y Outcome value or numeric vector.
#' @param w,delta Observed time and event indicator.
#' @param beta Outcome-model coefficient vector.
#' @param alpha1_star Coefficient vector for the working model for
#'   \code{X | Z}.
#' @param tau1 Standard deviation for the working model for \code{X | Z}.
#' @param z_c,z_d Optional continuous and discrete covariates.
#' @param w_min,w_max Truncation bounds for \code{X}.
#'
#' @return A numeric residual or vector of residuals.
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
#' beta <- c(0, 3)
#'
#' res_r1 <- r1(
#'   y_vec, w_vec, delta_vec, beta,
#'   alpha1_star = 0, tau1 = 1,
#'   w_min = -1, w_max = 1
#' )
#' head(res_r1)
#'
#' res_r1star <- r1(
#'   y_vec, w_vec, delta_vec, beta,
#'   alpha1_star = -2, tau1 = 1,
#'   w_min = -1, w_max = 1
#' )
#' head(res_r1star)
#'
#' @export
r1 <- function(y, w, delta, beta,
               alpha1_star, tau1,
               z_c = NULL, z_d = NULL,
               w_min, w_max) {

  y <- as.numeric(y)

  if (
    length(y) != length(w) ||
    length(w) != length(delta)
  ) {
    stop("y, w, and delta must have the same length.")
  }

  abs(
    y -
      m1(
        w,
        delta,
        beta,
        alpha1_star,
        tau1,
        z_c = z_c,
        z_d = z_d,
        w_min = w_min,
        w_max = w_max
      )
  )
}


#' Residual r2
#'
#' Computes \code{|Y - m0(W)|}.
#'
#' @param y Outcome value or numeric vector.
#' @param w,delta Observed time and event indicator.
#' @param beta Outcome-model coefficient vector.
#' @param alpha1_star,tau1 Included for a common interface with \code{r1};
#'   unused by \code{r2}.
#' @param z_c,z_d Optional continuous and discrete covariates.
#' @param w_min,w_max Included for a common interface with \code{r1};
#'   unused by \code{r2}.
#'
#' @return A numeric residual or vector of residuals.
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
#' beta <- c(0, 3)
#'
#' res_r2 <- r2(
#'   y_vec, w_vec, delta_vec, beta,
#'   alpha1_star = 0, tau1 = 1,
#'   w_min = -1, w_max = 1
#' )
#' head(res_r2)
#'
#' @export
r2 <- function(y, w, delta, beta,
               alpha1_star, tau1,
               z_c = NULL, z_d = NULL,
               w_min, w_max) {

  y <- as.numeric(y)

  if (
    length(y) != length(w) ||
    length(w) != length(delta)
  ) {
    stop("y, w, and delta must have the same length.")
  }

  abs(
    y -
      m0(
        w,
        beta,
        z_c = z_c,
        z_d = z_d
      )
  )
}
