# Shared fixtures for the test suite.
#
# `make_data()` mirrors the no-covariate data-generating process used in the
# simulation study, optionally with fully observed covariates attached, so that
# the tests exercise the same code paths the paper does.

W_MIN <- -1
W_MAX <-  1

#' Generate a small test dataset
#'
#' @param n Sample size.
#' @param seed Seed (set internally, so results are reproducible).
#' @param p_c,p_d Number of continuous / discrete covariates.
#' @param alpha1,alpha2 True means of X and C.
#' @param censor If FALSE, no observation is censored (delta == 1 throughout).
make_data <- function(n = 200, seed = 1, p_c = 0, p_d = 0,
                      alpha1 = 0, alpha2 = 1, sigma = 1,
                      tau1 = 1, tau2 = 1, censor = TRUE) {
  set.seed(seed)

  z_c <- if (p_c > 0) matrix(runif(n * p_c, -1, 1), n, p_c) else NULL
  z_d <- if (p_d > 0) matrix(rbinom(n * p_d, 1, 0.5), n, p_d) else NULL

  # mu_X(Z) = (1, z_c, z_d) %*% alpha1_vec
  alpha1_vec <- c(alpha1, rep(0.3, p_c), rep(-0.2, p_d))
  alpha2_vec <- c(alpha2, rep(0.1, p_c), rep(0.1, p_d))

  zc0 <- if (p_c > 0) z_c else matrix(0, n, 0)
  zd0 <- if (p_d > 0) z_d else matrix(0, n, 0)
  mu_x <- as.vector(cbind(1, zc0, zd0) %*% alpha1_vec)
  mu_c <- as.vector(cbind(1, zc0, zd0) %*% alpha2_vec)

  x <- truncnorm::rtruncnorm(n, a = W_MIN, b = W_MAX, mean = mu_x, sd = tau1)
  cc <- if (censor)
    truncnorm::rtruncnorm(n, a = W_MIN, b = W_MAX, mean = mu_c, sd = tau2)
  else
    rep(W_MAX + 10, n)

  w     <- pmin(x, cc)
  delta <- as.integer(x <= cc)

  # beta matches phi_xz = (1, x, z_c, z_d, x*z_c, x*z_d)
  beta <- c(0, 3, rep(0.5, p_c), rep(-0.7, p_d), rep(0.4, p_c), rep(0.6, p_d))

  design <- cbind(1, x, zc0, zd0,
                  if (p_c > 0) sweep(zc0, 1, x, `*`) else NULL,
                  if (p_d > 0) sweep(zd0, 1, x, `*`) else NULL)
  y <- as.vector(design %*% beta) + rnorm(n, 0, sigma)

  list(y = y, w = w, delta = delta, x = x,
       z_c = z_c, z_d = z_d,
       beta = beta, sigma = sigma,
       alpha1 = alpha1_vec, alpha2 = alpha2_vec,
       tau1 = tau1, tau2 = tau2,
       w_min = W_MIN, w_max = W_MAX,
       p_c = p_c, p_d = p_d, n = n)
}

#' Row i of a covariate block, or NULL when the block is absent
row_or_null <- function(z, i) if (is.null(z)) NULL else as.numeric(z[i, ])
