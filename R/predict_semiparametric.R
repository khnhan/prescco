# =============================================================================
#  Semiparametric prediction intervals (covariate-aware, data-in)
#  ---------------------------------------------------------------------------
#  Given fitted coefficients and working nuisance models, these routines build
#  the SPARCC influence-function arrays, solve for the efficient half-length
#  zeta on the observed data, and evaluate coverage on a test set. All are
#  covariate-aware (optional z_c, z_d). The user-facing entry point is
#  semiparametric_prediction_interval(); the c*/M*/b* helpers are internal.
#
#  Package imports are declared centrally in R/censcovpred-package.R.
# =============================================================================

#' Internal semiparametric-interval helper (covariate-aware).
#' @noRd
c1_xz_gauss_param12 = function(zeta, alpha, r,
                               beta, x_a,
                               z_c, z_d,
                               alpha1_star, alpha2_star, alpha1_star_r,
                               sigma, tau1, tau2, tau1_r,
                               w_min, w_max,
                               tt = 20, tt2 = 100) {
  
  cc <- gauss(tt2)  # Gaussian quadrature in Y
  
  ## Local basis for mu_Y(X,Z): 1, X, Z_c, Z_d, X Z_c, X Z_d
  phi_xz_local <- function(x, z_c = NULL, z_d = NULL) {
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
  
  ## temp2: E_C[ I(X < C) I{ r(Y, X, 1) <= zeta } + I(X >= C) I{ r(Y, C, 0) <= zeta } | Y, X, Z ]
  ## where C | Z ~ TN( (1,Z)^T alpha2_star, tau2^2; [w_min, w_max] )
  temp2 <- function(x, y) {
    z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
    z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)
    
    if (x < w_min) {
      ## Below support; treat as uncensored at X = x
      return(
        as.numeric(
          r(y, x, 1,
            beta,
            alpha1_star_r, tau1,
            z_c_vec, z_d_vec,
            w_min, w_max) <= zeta
        )
      )
    } else {
      ## Integrate over censoring time C in [w_min, w_max]
      c_grid <- seq(w_min, w_max, length.out = tt)
      
      mean_C <- sum(c(1, z_c_vec, z_d_vec) * alpha2_star)
      
      dens_c <- truncnorm::dtruncnorm(c_grid,
                           a    = w_min,
                           b    = w_max,
                           mean = mean_C,
                           sd   = tau2)
      
      ## Normalize in case of numerical drift
      dens_c <- dens_c / sum(dens_c)
      
      ## For each C, define observed (W, delta)
      w_grid     <- pmin(x, c_grid)
      delta_grid <- as.numeric(x <= c_grid)
      
      RHS <- sapply(1:tt, function(i) {
        as.numeric(
          r(y,
            w_grid[i],
            delta_grid[i],
            beta,
            alpha1_star_r, tau1,
            z_c_vec, z_d_vec,
            w_min, w_max) <= zeta
        )
      })
      
      return(as.numeric(RHS %*% dens_c))
    }
  }
  
  ## temp3: E[ temp2(Y, X, Z) | X, Z ]
  temp3 <- function(x) {
    z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
    z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)
    
    ## Mean of Y | X = x, Z = (z_c, z_d)
    mu_xz <- sum(phi_xz_local(x, z_c_vec, z_d_vec) * beta)
    
    ## Quadrature over Y = mu_xz + sigma * Z, Z ~ N(0,1)
    y_grid <- mu_xz + sigma * cc$x
    
    vals <- sapply(y_grid, function(y_norm) temp2(x, y_norm))
    
    as.numeric(vals %*% cc$w)
  }
  
  temp3 <- Vectorize(temp3, vectorize.args = "x")
  
  ## Return equation evaluated on x_a minus (1 - alpha)
  temp3(x_a) - (1 - alpha)
}
#' Internal semiparametric-interval helper (covariate-aware).
#' @noRd
M1_xz_gauss_param12 = function(beta, x_a,
                               z_c, z_d,
                               alpha1_star, alpha2_star,
                               sigma, tau1, tau2,
                               w_min, w_max,
                               tt = 20, tt2 = 20) {
  
  cc <- gauss(tt2)  # Gaussian quadrature for Y
  
  ## Ensure z_c, z_d are numeric vectors
  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)
  
  p_c <- length(z_c_vec)
  p_d <- length(z_d_vec)
  
  ## Indexing for beta under basis:
  ## (1, X, Z_c, Z_d, X Z_c, X Z_d)
  idx_x   <- 2L
  idx_zc  <- if (p_c > 0) 2L + seq_len(p_c) else integer(0)
  idx_zd  <- if (p_d > 0) 2L + p_c + seq_len(p_d) else integer(0)
  idx_xzc <- if (p_c > 0) 2L + p_c + p_d + seq_len(p_c) else integer(0)
  idx_xzd <- if (p_d > 0) 2L + p_c + p_d + p_c + seq_len(p_d) else integer(0)
  
  ## Outcome mean Y|X,Z: intercept + slope * X
  intercept_x <- beta[1] +
    (if (p_c > 0) sum(beta[idx_zc]  * z_c_vec) else 0) +
    (if (p_d > 0) sum(beta[idx_zd]  * z_d_vec) else 0)
  
  slope_x <- beta[idx_x] +
    (if (p_c > 0) sum(beta[idx_xzc] * z_c_vec) else 0) +
    (if (p_d > 0) sum(beta[idx_xzd] * z_d_vec) else 0)
  
  ## Prior mean of X|Z
  mu_x_prior <- sum(c(1, z_c_vec, z_d_vec) * alpha1_star)
  
  ## Mean of C|Z
  mu_c <- sum(c(1, z_c_vec, z_d_vec) * alpha2_star)
  
  ## temp[a(X,Z)] = E[ I(X > C) a(X,Z) | C, Y, Z ] / E[ I(X > C) | C, Y, Z ]
  temp <- function(x_grid, c, y) {
    ## Posterior X|Y,Z: N(eta_x, v_x), then truncated to [w_min, w_max]
    v_x  <- 1 / (slope_x^2 / sigma^2 + 1 / tau1^2)
    eta_x <- v_x * (slope_x * (y - intercept_x) / sigma^2 +
                      mu_x_prior / tau1^2)
    
    p <- truncnorm::dtruncnorm(x_grid,
                    a    = w_min,
                    b    = w_max,
                    mean = eta_x,
                    sd   = sqrt(v_x))
    
    num   <- (x_grid > c) * p
    denom <- sum(num)
    
    ifelse(is.nan(num / denom), 0, num / denom)
  }
  
  ## temp2[a(X,Z)] = E[ I(X > C) temp(x_a, C, Y, Z) | X, Y, Z ]
  temp2 <- function(x_grid, x, y) {
    c_grid <- seq(w_min, w_max, length.out = tt)
    
    dens <- truncnorm::dtruncnorm(c_grid,
                       a    = w_min,
                       b    = w_max,
                       mean = mu_c,
                       sd   = tau2)
    dens <- dens / sum(dens)
    
    sapply(c_grid, function(c_norm) {
      (x > c_norm) * temp(x_grid, c_norm, y)
    }) %*% dens
  }
  
  ## temp3[a(X,Z)] = E[ temp2(X, Y, Z) | X, Z ]
  temp3 <- function(x_grid, x) {
    ## Mean of Y | X = x, Z = (z_c, z_d)
    mu_xz <- intercept_x + slope_x * x
    
    y_grid <- mu_xz + sigma * cc$x
    
    sapply(y_grid, function(y_norm) {
      temp2(x_grid, x, y_norm)
    }) %*% cc$w
  }
  temp3 <- Vectorize(temp3, vectorize.args = "x")
  
  ## temp4[a(X,Z)] = E[ I(X < C) | Y, X, Z ] a(X,Z)  (here I(X < C) only depends on C|Z)
  temp4 <- function(x, y) {
    c_grid <- seq(w_min, w_max, length.out = tt)
    
    dens <- truncnorm::dtruncnorm(c_grid,
                       a    = w_min,
                       b    = w_max,
                       mean = mu_c,
                       sd   = tau2)
    dens <- dens / sum(dens)
    
    sum((x <= c_grid) * dens)
  }
  temp4 <- Vectorize(temp4, vectorize.args = "y")
  
  ## temp5[a(X,Z)] = E[ temp4(X, Y, Z) | X, Z ] a(X,Z)
  temp5 <- function(x) {
    mu_xz  <- intercept_x + slope_x * x
    y_grid <- mu_xz + sigma * cc$x
    
    sum(temp4(x, y_grid) * cc$w)
  }
  temp5 <- Vectorize(temp5, vectorize.args = "x")
  
  ## Final matrix: diag(temp5(x_a)) + t(temp3(x_a, x_a))
  res <- diag(temp5(x_a)) + t(temp3(x_a, x_a))
  res
}
#' Internal semiparametric-interval helper (covariate-aware).
#' @noRd
b1_gauss_param12 = function(zeta, alpha, r,
                            beta, x_a,
                            z_c, z_d,
                            alpha1_star, alpha2_star, alpha1_star_r,
                            sigma, tau1, tau2,
                            w_min, w_max,
                            tt = 20, tt2 = 100) {
  
  M1_mat <- M1_xz_gauss_param12(
    beta        = beta,
    x_a         = x_a,
    z_c         = z_c,
    z_d         = z_d,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    sigma       = sigma,
    tau1        = tau1,
    tau2        = tau2,
    w_min       = w_min,
    w_max       = w_max,
    tt          = tt,
    tt2         = tt2
  )
  
  c1_vec <- c1_xz_gauss_param12(
    zeta        = zeta,
    alpha       = alpha,
    r           = r,
    beta        = beta,
    x_a         = x_a,
    z_c         = z_c,
    z_d         = z_d,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    alpha1_star_r = alpha1_star_r,
    sigma       = sigma,
    tau1        = tau1,
    tau2        = tau2,
    w_min       = w_min,
    w_max       = w_max,
    tt          = tt,
    tt2         = tt2
  )
  
  MASS::ginv(M1_mat) %*% c1_vec
}
#' Internal semiparametric-interval helper (covariate-aware).
#' @noRd
c2_xz_gauss_param12 = function(zeta, alpha, r,
                               beta, x_a,
                               z_c, z_d,
                               alpha1_star, alpha2_star, alpha1_star_r,
                               sigma, tau1, tau2,
                               w_min, w_max,
                               tt = 20, tt2 = 100) {
  
  cc <- gauss(tt2)
  
  ## Ensure z_c, z_d are numeric vectors
  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)
  
  p_c <- length(z_c_vec)
  p_d <- length(z_d_vec)
  
  ## Basis indexing for beta under:
  ## (1, X, Z_c, Z_d, X Z_c, X Z_d)
  idx_x   <- 2L
  idx_zc  <- if (p_c > 0) 2L + seq_len(p_c) else integer(0)
  idx_zd  <- if (p_d > 0) 2L + p_c + seq_len(p_d) else integer(0)
  idx_xzc <- if (p_c > 0) 2L + p_c + p_d + seq_len(p_c) else integer(0)
  idx_xzd <- if (p_d > 0) 2L + p_c + p_d + p_c + seq_len(p_d) else integer(0)
  
  ## Y | X,Z has mean: intercept_x + slope_x * X
  intercept_x <- beta[1] +
    (if (p_c > 0) sum(beta[idx_zc]  * z_c_vec) else 0) +
    (if (p_d > 0) sum(beta[idx_zd]  * z_d_vec) else 0)
  
  slope_x <- beta[idx_x] +
    (if (p_c > 0) sum(beta[idx_xzc] * z_c_vec) else 0) +
    (if (p_d > 0) sum(beta[idx_xzd] * z_d_vec) else 0)
  
  ## X | Z ~ N(mu_x_prior, tau1^2), truncated to [w_min, w_max]
  mu_x_prior <- sum(c(1, z_c_vec, z_d_vec) * alpha1_star)
  
  ## temp2(c, x) = E_Y[ I{ r(Y, W, delta) <= zeta } | C = c, X = x, Z ]
  ## where W = min(X, C), delta = 1{ X <= C }
  temp2 <- function(c, x) {
    # mean of Y | X = x, Z
    mu_xz  <- intercept_x + slope_x * x
    y_grid <- mu_xz + sigma * cc$x
    
    w_obs    <- min(x, c)
    delta_obs <- as.numeric(x <= c)
    
    vals <- sapply(y_grid, function(y_norm) {
      as.numeric(
        r(y_norm,
          w_obs,
          delta_obs,
          beta,
          alpha1_star_r, tau1,
          z_c_vec, z_d_vec,
          w_min, w_max) <= zeta
      )
    })
    
    as.numeric(vals %*% cc$w)
  }
  
  ## temp3(c) = E_X[ temp2(c, X) | Z ] over X|Z
  temp3 <- function(c) {
    x_grid <- seq(w_min, w_max, length.out = tt)
    
    dens_x <- truncnorm::dtruncnorm(x_grid,
                         a    = w_min,
                         b    = w_max,
                         mean = mu_x_prior,
                         sd   = tau1)
    dens_x <- dens_x / sum(dens_x)
    
    RHS <- vapply(seq_len(tt), function(i) {
      temp2(c, x_grid[i])
    }, numeric(1))
    
    as.numeric(RHS %*% dens_x)
  }
  
  temp3 <- Vectorize(temp3, vectorize.args = "c")
  
  ## Equation evaluated on c-grid x_a
  temp3(x_a) - (1 - alpha)
}
#' Internal semiparametric-interval helper (covariate-aware).
#' @noRd
M2_xz_gauss_param12 = function(beta, x_a,
                               z_c, z_d,
                               alpha1_star, alpha2_star,
                               sigma, tau1, tau2,
                               w_min, w_max,
                               tt = 50, tt2 = 50) {
  
  ## z_c, z_d are one covariate point (can be vectors)
  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)
  
  ## Means for X|Z and C|Z
  mu_x <- sum(c(1, z_c_vec, z_d_vec) * alpha1_star)
  mu_c <- sum(c(1, z_c_vec, z_d_vec) * alpha2_star)
  
  ## temp[a(X,Z)] = E[ I(C >= x) a(C,Z) | "C on grid x_a", Z ] /
  ##                E[ I(C >= x) | Z ]  (implemented via discretization on x_a)
  temp <- function(x_grid, x) {
    p <- truncnorm::dtruncnorm(x_grid,
                    a    = w_min,
                    b    = w_max,
                    mean = mu_c,
                    sd   = tau2)
    
    num   <- (x_grid >= x) * p
    denom <- sum(num)
    
    ifelse(is.nan(num / denom), 0, num / denom)
  }
  
  ## temp2[a(X,Z)] = E[ I(X <= C) temp(x_a, X) | C, Z ]
  ## integrate over X|Z
  temp2 <- function(x_grid, c) {
    x_grid_int <- seq(w_min, w_max, length.out = tt)
    
    dens_x <- truncnorm::dtruncnorm(x_grid_int,
                         a    = w_min,
                         b    = w_max,
                         mean = mu_x,
                         sd   = tau1)
    dens_x <- dens_x / sum(dens_x)
    
    sapply(x_grid_int, function(x_norm) {
      (x_norm <= c) * temp(x_grid, x_norm)
    }) %*% dens_x
  }
  temp2 <- Vectorize(temp2, vectorize.args = "c")
  
  ## temp3(c) = E[ I(X > C) | C, Z ]
  temp3 <- function(c) {
    x_grid_int <- seq(w_min, w_max, length.out = tt)
    
    dens_x <- truncnorm::dtruncnorm(x_grid_int,
                         a    = w_min,
                         b    = w_max,
                         mean = mu_x,
                         sd   = tau1)
    dens_x <- dens_x / sum(dens_x)
    
    sum((x_grid_int > c) * dens_x)
  }
  temp3 <- Vectorize(temp3, vectorize.args = "c")
  
  ## Final operator: diag(temp3(x_a)) + t(temp2(x_a, x_a))
  res <- diag(temp3(x_a)) + t(temp2(x_a, x_a))
  res
}
#' Internal semiparametric-interval helper (covariate-aware).
#' @noRd
b2_gauss_param12 = function(zeta, alpha, r,
                            beta, x_a,
                            z_c, z_d,
                            alpha1_star, alpha2_star, alpha1_star_r,
                            sigma, tau1, tau2,
                            w_min, w_max,
                            tt = 50, tt2 = 50) {
  
  M2_mat <- M2_xz_gauss_param12(
    beta        = beta,
    x_a         = x_a,
    z_c         = z_c,
    z_d         = z_d,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    sigma       = sigma,
    tau1        = tau1,
    tau2        = tau2,
    w_min       = w_min,
    w_max       = w_max,
    tt          = tt,
    tt2         = tt2
  )
  
  c2_vec <- c2_xz_gauss_param12(
    zeta        = zeta,
    alpha       = alpha,
    r           = r,
    beta        = beta,
    x_a         = x_a,
    z_c         = z_c,
    z_d         = z_d,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    alpha1_star_r = alpha1_star_r,
    sigma       = sigma,
    tau1        = tau1,
    tau2        = tau2,
    w_min       = w_min,
    w_max       = w_max,
    tt          = tt,
    tt2         = tt2
  )
  
  MASS::ginv(M2_mat) %*% c2_vec
}
#' Internal semiparametric-interval helper (covariate-aware).
#' @noRd
b3_gauss_param12 = function(zeta, alpha, r,
                            beta, x_a,
                            z_c, z_d,
                            alpha1_star, alpha2_star, alpha1_star_r,
                            sigma, tau1, tau2,
                            w_min, w_max,
                            tt = 50, tt2 = 50, tt3 = 50) {
  
  cc <- gauss(tt2)  # quadrature for Y
  
  ## z_c, z_d are one covariate point
  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)
  
  p_c <- length(z_c_vec)
  p_d <- length(z_d_vec)
  
  ## Basis indexing for beta under:
  ## (1, X, Z_c, Z_d, X Z_c, X Z_d)
  idx_x   <- 2L
  idx_zc  <- if (p_c > 0) 2L + seq_len(p_c) else integer(0)
  idx_zd  <- if (p_d > 0) 2L + p_c + seq_len(p_d) else integer(0)
  idx_xzc <- if (p_c > 0) 2L + p_c + p_d + seq_len(p_c) else integer(0)
  idx_xzd <- if (p_d > 0) 2L + p_c + p_d + p_c + seq_len(p_d) else integer(0)
  
  ## Y | X,Z has mean: intercept_x + slope_x * X
  intercept_x <- beta[1] +
    (if (p_c > 0) sum(beta[idx_zc]  * z_c_vec) else 0) +
    (if (p_d > 0) sum(beta[idx_zd]  * z_d_vec) else 0)
  
  slope_x <- beta[idx_x] +
    (if (p_c > 0) sum(beta[idx_xzc] * z_c_vec) else 0) +
    (if (p_d > 0) sum(beta[idx_xzd] * z_d_vec) else 0)
  
  ## X | Z ~ TN(mu_x, tau1^2; [w_min, w_max])
  mu_x <- sum(c(1, z_c_vec, z_d_vec) * alpha1_star)
  
  ## C | Z ~ TN(mu_c, tau2^2; [w_min, w_max])
  mu_c <- sum(c(1, z_c_vec, z_d_vec) * alpha2_star)
  
  ## temp2(c, x) = E_Y[ I{ r(Y, W, delta) <= zeta } | C=c, X=x, Z ]
  ## where W = min(X, C), delta = 1{ X <= C }
  temp2 <- function(c, x) {
    mu_xz  <- intercept_x + slope_x * x
    y_grid <- mu_xz + sigma * cc$x
    
    w_obs     <- min(x, c)
    delta_obs <- as.numeric(x <= c)
    
    vals <- sapply(y_grid, function(y_norm) {
      as.numeric(
        r(y_norm,
          w_obs,
          delta_obs,
          beta,
          alpha1_star_r, tau1,
          z_c_vec, z_d_vec,
          w_min, w_max) <= zeta
      )
    })
    
    as.numeric(vals %*% cc$w)
  }
  
  ## temp3(c) = E_X[ temp2(c, X) | Z ]
  temp3 <- function(c) {
    x_grid <- seq(w_min, w_max, length.out = tt)
    
    dens_x <- truncnorm::dtruncnorm(x_grid,
                         a    = w_min,
                         b    = w_max,
                         mean = mu_x,
                         sd   = tau1)
    dens_x <- dens_x / sum(dens_x)
    
    RHS <- vapply(seq_len(tt), function(i) {
      temp2(c, x_grid[i])
    }, numeric(1))
    
    as.numeric(RHS %*% dens_x)
  }
  
  ## temp4 = E_C[ temp3(C) | Z ]
  temp4 <- function() {
    c_grid <- seq(w_min, w_max, length.out = tt3)
    
    dens_c <- truncnorm::dtruncnorm(c_grid,
                         a    = w_min,
                         b    = w_max,
                         mean = mu_c,
                         sd   = tau2)
    dens_c <- dens_c / sum(dens_c)
    
    vals <- vapply(c_grid, temp3, numeric(1))
    
    as.numeric(vals %*% dens_c)
  }
  
  temp4() - (1 - alpha)
}
#' Internal semiparametric-interval helper (covariate-aware).
#' @noRd
trans_phi_eff_gauss_param12 = function(beta, y, w, delta,
                                       x_a, b1, b2, b3,
                                       z_c, z_d,
                                       alpha1_star, alpha2_star,
                                       sigma, tau1, tau2,
                                       w_min, w_max) {
  
  ## z_c, z_d: one covariate point (can be vectors)
  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)
  
  p_c <- length(z_c_vec)
  p_d <- length(z_d_vec)
  
  ## ---- Outcome model indexing: basis (1, X, Z_c, Z_d, X Z_c, X Z_d) ----
  idx_x   <- 2L
  idx_zc  <- if (p_c > 0) 2L + seq_len(p_c) else integer(0)
  idx_zd  <- if (p_d > 0) 2L + p_c + seq_len(p_d) else integer(0)
  idx_xzc <- if (p_c > 0) 2L + p_c + p_d + seq_len(p_c) else integer(0)
  idx_xzd <- if (p_d > 0) 2L + p_c + p_d + p_c + seq_len(p_d) else integer(0)
  
  ## Intercept and slope in X for Y|X,Z
  intercept_x <- beta[1] +
    (if (p_c > 0) sum(beta[idx_zc]  * z_c_vec) else 0) +
    (if (p_d > 0) sum(beta[idx_zd]  * z_d_vec) else 0)
  
  slope_x <- beta[idx_x] +
    (if (p_c > 0) sum(beta[idx_xzc] * z_c_vec) else 0) +
    (if (p_d > 0) sum(beta[idx_xzd] * z_d_vec) else 0)
  
  ## Prior means: X|Z and C|Z
  mu_x_prior <- sum(c(1, z_c_vec, z_d_vec) * alpha1_star)
  mu_c       <- sum(c(1, z_c_vec, z_d_vec) * alpha2_star)
  
  if (delta == 1) {
    ## -------- Uncensored: X = W --------
    ## temp[a(C,Z)] = E[ I(C >= X) a(C,Z) | X,Z ] / E[ I(C >= X) | X,Z ]
    temp <- function(x_grid, x, y) {
      ## C|Z ~ TN(mu_c, tau2^2; [w_min, w_max]) evaluated on x_grid (as C grid)
      p <- truncnorm::dtruncnorm(x_grid,
                      a    = w_min,
                      b    = w_max,
                      mean = mu_c,
                      sd   = tau2)
      
      num   <- (x_grid >= x) * p
      denom <- sum(num)
      
      ifelse(is.nan(num / denom), 0, num / denom)
    }
    
    ## Interpolate b1 at w
    b1_w <- approx(x_a, b1, w, rule = 2)$y
    
    ## Return: b1(w) + b2^T temp - b3
    return(b1_w + as.numeric(t(b2) %*% temp(x_a, w, y)) - b3)
    
  } else {
    ## -------- Censored: X > W, only W and Y observed --------
    ## temp[a(X,Z)] = E[ I(X > C) a(X,Z) | C, Y, Z ] / E[ I(X > C) | C, Y, Z ]
    temp <- function(x_grid, c, y) {
      ## Posterior X|Y,Z: N(eta_x, v_x), truncated to [w_min, w_max]
      v_x  <- 1 / (slope_x^2 / sigma^2 + 1 / tau1^2)
      eta_x <- v_x * (slope_x * (y - intercept_x) / sigma^2 +
                        mu_x_prior / tau1^2)
      
      p <- truncnorm::dtruncnorm(x_grid,
                      a    = w_min,
                      b    = w_max,
                      mean = eta_x,
                      sd   = sqrt(v_x))
      
      num   <- (x_grid > c) * p
      denom <- sum(num)
      
      ifelse(is.nan(num / denom), 0, num / denom)
    }
    
    ## Interpolate b2 at w
    b2_w <- approx(x_a, b2, w, rule = 2)$y
    
    ## Return: b1^T temp + b2(w) - b3
    return(as.numeric(t(b1) %*% temp(x_a, w, y)) + b2_w - b3)
  }
}
#' Build the b1/b2/b3 estimating-function arrays
#'
#' Constructs the three influence-function arrays over the half-length grid
#' `zeta_seq` and a covariate grid, for one residual. Covariate-aware.
#' @param zeta_seq Numeric grid of candidate half-lengths.
#' @param alpha Miscoverage level.
#' @param r Residual function (r1 or r2).
#' @param beta Fitted outcome coefficients.
#' @param x_a Numeric grid over the event-time covariate.
#' @param z_c_data,z_d_data Optional covariates from the fitting sample.
#' @param alpha1_star,alpha2_star Nuisance mean coefficient vectors.
#' @param alpha1_star_r Center used inside the residual.
#' @param sigma,tau1,tau2 Model standard deviations.
#' @param w_min,w_max Truncation bounds.
#' @param n_c_grid Covariate-grid resolution for continuous covariates.
#' @param b1_tt,b1_tt2,b2_tt,b2_tt2,b3_tt,b3_tt2,b3_tt3 Integration grid sizes.
#'   Defaults match those declared by \code{b1_gauss_param12},
#'   \code{b2_gauss_param12}, and \code{b3_gauss_param12}, so the function can
#'   be called without specifying them. Larger values are more accurate and
#'   considerably slower; the simulation study uses
#'   \code{b1_tt = b2_tt = 100}, \code{b1_tt2 = b2_tt2 = b3_tt2 = 500},
#'   \code{b3_tt = 20}, \code{b3_tt3 = 200}.
#' @param verbose If \code{TRUE}, report progress through the zeta and
#'   covariate grids. Defaults to \code{FALSE}.
#' @return A list with b1_array, b2_array, b3_array.
#' @noRd
build_b123_arrays <- function(zeta_seq, alpha, r,
                              beta, x_a,
                              z_c_data, z_d_data,
                              alpha1_star, alpha2_star, alpha1_star_r,
                              sigma, tau1, tau2,
                              w_min, w_max,
                              n_c_grid = 10,
                              b1_tt = 20, b1_tt2 = 100,
                              b2_tt = 50, b2_tt2 = 50,
                              b3_tt = 50, b3_tt2 = 50, b3_tt3 = 50,
                              verbose = FALSE) {
  ## --- coerce z_c, z_d to matrices ---
  n <- length(x_a)  # just to have a length for NULL case
  
  if (is.null(z_c_data)) {
    z_c_mat <- matrix(0, nrow = if (is.null(z_d_data)) 1 else length(z_d_data), ncol = 0)
    p_c     <- 0
  } else if (is.null(dim(z_c_data))) {
    z_c_mat <- matrix(as.numeric(z_c_data), ncol = 1)
    p_c     <- 1
  } else {
    z_c_mat <- as.matrix(z_c_data)
    p_c     <- ncol(z_c_mat)
  }
  
  if (is.null(z_d_data)) {
    z_d_mat <- matrix(0, nrow = nrow(z_c_mat), ncol = 0)
    p_d     <- 0
  } else if (is.null(dim(z_d_data))) {
    z_d_mat <- matrix(z_d_data, ncol = 1)
    p_d     <- 1
  } else {
    z_d_mat <- as.matrix(z_d_data)
    p_d     <- ncol(z_d_mat)
  }
  
  ## --- continuous grid for z_c: tensor product of n_c_grid knots per dim ---
  if (p_c > 0) {
    z_c_knots_list <- lapply(seq_len(p_c), function(j) {
      rng <- range(z_c_mat[, j], na.rm = TRUE)
      seq(rng[1], rng[2], length.out = n_c_grid)
    })
    z_c_grid <- as.matrix(do.call(expand.grid, z_c_knots_list))
  } else {
    z_c_grid <- matrix(0, nrow = 1, ncol = 0)
  }
  n_zc_grid <- nrow(z_c_grid)
  
  ## --- all combinations of discrete z_d ---
  if (p_d > 0) {
    z_d_levels_list <- lapply(seq_len(p_d), function(j) sort(unique(z_d_mat[, j])))
    z_d_grid <- as.matrix(do.call(expand.grid, z_d_levels_list))
  } else {
    z_d_grid <- matrix(0, nrow = 1, ncol = 0)
  }
  n_zd_comb <- nrow(z_d_grid)
  
  ## --- allocate arrays ---
  n_zeta <- length(zeta_seq)
  m      <- length(x_a)
  
  b1_array <- array(NA_real_, dim = c(n_zeta, m, n_zc_grid, n_zd_comb))
  b2_array <- array(NA_real_, dim = c(n_zeta, m, n_zc_grid, n_zd_comb))
  b3_array <- array(NA_real_, dim = c(n_zeta,     n_zc_grid, n_zd_comb))
  
  for (k in seq_along(zeta_seq)) {
    zeta <- zeta_seq[k]
    if (verbose) message("  zeta ", k, "/", length(zeta_seq))
    for (g in seq_len(n_zc_grid)) {
      z_c_pt <- if (p_c > 0) z_c_grid[g, ] else NULL
      if (verbose) message("    z_c grid point ", g, "/", n_zc_grid)
      for (h in seq_len(n_zd_comb)) {
        z_d_pt <- if (p_d > 0) z_d_grid[h, ] else NULL
        if (verbose) message("      z_d combination ", h, "/", n_zd_comb)
        b1_array[k, , g, h] <- b1_gauss_param12(
          zeta, alpha, r,
          beta, x_a,
          z_c = z_c_pt, z_d = z_d_pt,
          alpha1_star   = alpha1_star,
          alpha2_star   = alpha2_star,
          alpha1_star_r = alpha1_star_r,
          sigma = sigma, tau1 = tau1, tau2 = tau2,
          w_min = w_min, w_max = w_max,
          tt   = b1_tt,  tt2  = b1_tt2
        )
        b2_array[k, , g, h] <- b2_gauss_param12(
          zeta, alpha, r,
          beta, x_a,
          z_c = z_c_pt, z_d = z_d_pt,
          alpha1_star   = alpha1_star,
          alpha2_star   = alpha2_star,
          alpha1_star_r = alpha1_star_r,
          sigma = sigma, tau1 = tau1, tau2 = tau2,
          w_min = w_min, w_max = w_max,
          tt   = b2_tt,  tt2  = b2_tt2
        )
        b3_array[k, g, h] <- b3_gauss_param12(
          zeta, alpha, r,
          beta, x_a,
          z_c = z_c_pt, z_d = z_d_pt,
          alpha1_star   = alpha1_star,
          alpha2_star   = alpha2_star,
          alpha1_star_r = alpha1_star_r,
          sigma = sigma, tau1 = tau1, tau2 = tau2,
          w_min = w_min, w_max = w_max,
          tt   = b3_tt,  tt2  = b3_tt2, tt3 = b3_tt3
        )
      }
    }
  }
  
  list(
    b1_array = b1_array,   # (n_zeta, m, n_zc_grid, n_zd_comb)
    b2_array = b2_array,   # (n_zeta, m, n_zc_grid, n_zd_comb)
    b3_array = b3_array,   # (n_zeta,     n_zc_grid, n_zd_comb)
    z_c_grid = z_c_grid,   # for later interpolation
    z_d_grid = z_d_grid
  )
}
#' Internal semiparametric-interval helper (covariate-aware).
#' @noRd
find_sol <- function(x_vec, y_vec) {
  sgn  <- sign(y_vec)
  idxc <- which(diff(sgn) != 0)
  
  if (!length(idxc)) {
    if (all(y_vec > 0)) {
      idx <- which.min(y_vec)
      x1 <- x_vec[idx];   x2 <- x_vec[idx + 1]
      y1 <- y_vec[idx];   y2 <- y_vec[idx + 1]
      return(x1 - y1 * (x2 - x1) / (y2 - y1))
    } else {
      idx <- which.max(y_vec)
      x1 <- x_vec[idx - 1]; x2 <- x_vec[idx]
      y1 <- y_vec[idx - 1]; y2 <- y_vec[idx]
      return(x2 - y2 * (x2 - x1) / (y2 - y1))
    }
  } else {
    idx <- idxc[which.min(abs(y_vec[idxc]))]
    x1 <- x_vec[idx];   x2 <- x_vec[idx + 1]
    y1 <- y_vec[idx];   y2 <- y_vec[idx + 1]
    return(x1 - y1 * (x2 - x1) / (y2 - y1))
  }
}
#' Solve for the estimated half-length on given data
#'
#' Interpolates the b1/b2/b3 arrays at the fitted beta and locates the zeta
#' root of the efficient estimating function. Covariate-aware.
#' @param y_data,w_data,delta_data Fitting-sample outcome, time, indicator.
#' @param z_c_data,z_d_data Optional covariates.
#' @param b1_array,b2_array,b3_array Arrays from build_b123_arrays.
#' @param beta_temp Fitted coefficient vector.
#' @param zeta_seq Half-length grid.
#' @param alpha1_star,alpha2_star Nuisance mean coefficient vectors.
#' @param sigma,tau1,tau2 Model standard deviations.
#' @param m Covariate-grid resolution.
#' @param w_min,w_max Truncation bounds.
#' @return A list with `vec` (grid and objective values) and `sol` (the root).
#' @noRd
find_zeta_param_int <- function(y_data, w_data, delta_data,
                                z_c_data, z_d_data,
                                b1_array, b2_array, b3_array,
                                beta_temp,
                                zeta_seq,
                                alpha1_star,
                                alpha2_star,
                                sigma, tau1, tau2,
                                m = 20,
                                w_min, w_max) {
  n <- length(y_data)
  if (length(w_data) != n || length(delta_data) != n)
    stop("y_data, w_data, delta_data must have same length.")
  
  ## --- coerce z_c, z_d into matrices ---
  if (is.null(z_c_data)) {
    z_c_mat <- matrix(0, nrow = n, ncol = 0); p_c <- 0
  } else if (is.null(dim(z_c_data))) {
    if (length(z_c_data) != n) stop("z_c_data must have length n.")
    z_c_mat <- matrix(as.numeric(z_c_data), ncol = 1); p_c <- 1
  } else {
    if (nrow(z_c_data) != n) stop("z_c_data must have n rows.")
    z_c_mat <- as.matrix(z_c_data); p_c <- ncol(z_c_mat)
  }
  
  if (is.null(z_d_data)) {
    z_d_mat <- matrix(0, nrow = n, ncol = 0); p_d <- 0
  } else if (is.null(dim(z_d_data))) {
    if (length(z_d_data) != n) stop("z_d_data must have length n.")
    z_d_mat <- matrix(z_d_data, ncol = 1); p_d <- 1
  } else {
    if (nrow(z_d_data) != n) stop("z_d_data must have n rows.")
    z_d_mat <- as.matrix(z_d_data); p_d <- ncol(z_d_mat)
  }
  
  ## --- reconstruct grids from array dims + data ranges ---
  dims_b1    <- dim(b1_array)   # (n_zeta, m, n_zc_grid, n_zd_comb)
  n_zc_grid  <- dims_b1[3]
  n_zd_comb  <- dims_b1[4]
  
  if (p_c > 0) {
    n_c_grid <- round(n_zc_grid^(1 / p_c))
    z_c_knots_list <- lapply(seq_len(p_c), function(j) {
      rng <- range(z_c_mat[, j], na.rm = TRUE)
      seq(rng[1], rng[2], length.out = n_c_grid)
    })
    z_c_grid <- as.matrix(do.call(expand.grid, z_c_knots_list))
    
    knots_list <- lapply(seq_len(p_c), function(j) sort(unique(z_c_grid[, j])))
    L_vec      <- vapply(knots_list, length, integer(1))
    cumprod_L  <- if (p_c == 1) 1L else c(1L, cumprod(L_vec[-length(L_vec)]))
    corner_bits <- as.matrix(expand.grid(rep(list(c(0, 1)), p_c)))
  } else {
    z_c_grid <- matrix(0, nrow = 1, ncol = 0)
    knots_list <- list(); cumprod_L <- 1L
    corner_bits <- matrix(0, nrow = 1, ncol = 0)
  }
  
  if (p_d > 0) {
    z_d_levels_list <- lapply(seq_len(p_d), function(j) sort(unique(z_d_mat[, j])))
    z_d_grid <- as.matrix(do.call(expand.grid, z_d_levels_list))
    if (nrow(z_d_grid) != n_zd_comb)
      warning("nrow(z_d_grid) != dim(b1_array)[4]; check discrete grid.")
  } else {
    z_d_grid <- matrix(0, nrow = 1, ncol = 0)
  }
  
  ## precompute z_d index for each obs
  idx_zd_vec <- integer(n)
  if (p_d == 0) {
    idx_zd_vec[] <- 1L
  } else {
    for (i in seq_len(n)) {
      row_match <- which(apply(z_d_grid, 1, function(r) all(r == z_d_mat[i, ])))
      if (!length(row_match)) stop("z_d pattern not found for obs ", i)
      idx_zd_vec[i] <- row_match[1]
    }
  }
  
  ## precompute low/high indices & t for z_c for each obs
  if (p_c > 0) {
    idx_low_mat  <- matrix(0L, n, p_c)
    idx_high_mat <- matrix(0L, n, p_c)
    t_mat        <- matrix(0, n, p_c)
    for (i in seq_len(n)) {
      z_i <- z_c_mat[i, ]
      for (j in seq_len(p_c)) {
        kz  <- knots_list[[j]]
        val <- z_i[j]
        if (val <= kz[1]) {
          idx_low_mat[i, j]  <- 1L
          idx_high_mat[i, j] <- 1L
          t_mat[i, j]        <- 0
        } else if (val >= kz[length(kz)]) {
          idx_low_mat[i, j]  <- length(kz)
          idx_high_mat[i, j] <- length(kz)
          t_mat[i, j]        <- 0
        } else {
          hi <- which(kz >= val)[1]
          lo <- hi - 1L
          idx_low_mat[i, j]  <- lo
          idx_high_mat[i, j] <- hi
          t_mat[i, j]        <- (val - kz[lo]) / (kz[hi] - kz[lo])
        }
      }
    }
  }
  
  ## local helper: for a given zeta_idx, get (b1,b2,b3) for all obs
  get_b123_for_zeta <- function(zeta_idx, x_a) {
    m <- length(x_a)
    b1_obs <- matrix(NA_real_, n, m)
    b2_obs <- matrix(NA_real_, n, m)
    b3_obs <- numeric(n)
    
    for (i in seq_len(n)) {
      idx_zd <- idx_zd_vec[i]
      
      if (p_c == 0) {
        idx_zc <- 1L
        b1_obs[i, ] <- b1_array[zeta_idx, , idx_zc, idx_zd]
        b2_obs[i, ] <- b2_array[zeta_idx, , idx_zc, idx_zd]
        b3_obs[i]   <- b3_array[zeta_idx, idx_zc, idx_zd]
      } else {
        # multilinear interpolation in z_c
        idx_low  <- idx_low_mat[i, ]
        idx_high <- idx_high_mat[i, ]
        t_vec    <- t_mat[i, ]
        
        # 2^p_c corners
        v3    <- 0; w3 <- 0
        for (r in seq_len(nrow(corner_bits))) {
          idx_j <- integer(p_c)
          w     <- 1
          for (j in seq_len(p_c)) {
            if (corner_bits[r, j] == 0) {
              idx_j[j] <- idx_low[j];  w <- w * (1 - t_vec[j])
            } else {
              idx_j[j] <- idx_high[j]; w <- w * t_vec[j]
            }
          }
          row_idx <- 1L + sum((idx_j - 1L) * cumprod_L)
          
          # b1,b2: need over x_a; b3: scalar
          for (xj in seq_len(m)) {
            b1_obs[i, xj] <- if (is.na(b1_obs[i, xj])) 0 else b1_obs[i, xj]
            b2_obs[i, xj] <- if (is.na(b2_obs[i, xj])) 0 else b2_obs[i, xj]
            b1_obs[i, xj] <- b1_obs[i, xj] + w * b1_array[zeta_idx, xj, row_idx, idx_zd]
            b2_obs[i, xj] <- b2_obs[i, xj] + w * b2_array[zeta_idx, xj, row_idx, idx_zd]
          }
          v3 <- v3 + w * b3_array[zeta_idx, row_idx, idx_zd]
          w3 <- w3 + w
        }
        b1_obs[i, ] <- b1_obs[i, ] / w3
        b2_obs[i, ] <- b2_obs[i, ] / w3
        b3_obs[i]   <- v3 / w3
      }
    }
    list(b1 = b1_obs, b2 = b2_obs, b3 = b3_obs)
  }
  
  x_a <- seq(w_min, w_max, length.out = m)
  func_vals <- numeric(length(zeta_seq))
  
  for (zeta_idx in seq_along(zeta_seq)) {
    zeta <- zeta_seq[zeta_idx]
    
    b123 <- get_b123_for_zeta(zeta_idx, x_a)
    b1_obs <- b123$b1
    b2_obs <- b123$b2
    b3_obs <- b123$b3
    
    val <- 0
    for (i in seq_len(n)) {
      z_c_i <- if (p_c > 0) z_c_mat[i, ] else NULL
      z_d_i <- if (p_d > 0) z_d_mat[i, ] else NULL
      
      val <- val + trans_phi_eff_gauss_param12(
        beta        = beta_temp,
        y           = y_data[i],
        w           = w_data[i],
        delta       = delta_data[i],
        x_a         = x_a,
        b1          = b1_obs[i, ],
        b2          = b2_obs[i, ],
        b3          = b3_obs[i],
        z_c         = z_c_i,
        z_d         = z_d_i,
        alpha1_star = alpha1_star,
        alpha2_star = alpha2_star,
        sigma       = sigma,
        tau1        = tau1,
        tau2        = tau2,
        w_min       = w_min,
        w_max       = w_max
      )
    }
    func_vals[zeta_idx] <- val
  }
  
  list(
    vec = rbind(zeta_seq, func_vals),
    sol = find_sol(zeta_seq, func_vals)
  )
}
#' Empirical coverage of a half-length on a test set
#'
#' @param zeta Half-length to evaluate.
#' @param r Residual function (r1 or r2).
#' @param test_y_data,test_w_data,test_delta_data Test outcome, time, indicator.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param beta Fitted coefficients.
#' @param alpha1_star_r Center used inside the residual.
#' @param tau1 Standard deviation of X.
#' @param w_min,w_max Truncation bounds.
#' @return The empirical coverage rate (fraction of residuals <= zeta).
#' @examples
#' set.seed(1)
#' n  <- 200
#' x  <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' cc <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w  <- pmin(x, cc)
#' delta <- as.integer(x <= cc)
#' y  <- 3 * x + rnorm(n)          # beta = c(0, 3)
#'
#' beta <- find_beta_cc(y, w, delta)$beta_cc
#'
#' # Fresh test data
#' set.seed(2)
#' xt  <- truncnorm::rtruncnorm(500, a = -1, b = 1, mean = 0, sd = 1)
#' ct  <- truncnorm::rtruncnorm(500, a = -1, b = 1, mean = 1, sd = 1)
#' wt  <- pmin(xt, ct); dt <- as.integer(xt <= ct)
#' yt  <- 3 * xt + rnorm(500)
#'
#' # Coverage rate rises with zeta, from 0 to 1
#' prediction_coverage_rate(0.5, r1, yt, wt, dt, NULL, NULL,
#'                          beta, alpha1_star_r = 0, tau1 = 1,
#'                          w_min = -1, w_max = 1)
#' prediction_coverage_rate(3.0, r1, yt, wt, dt, NULL, NULL,
#'                          beta, alpha1_star_r = 0, tau1 = 1,
#'                          w_min = -1, w_max = 1)
#' @export
prediction_coverage_rate = function(zeta, r,
                    test_y_data, test_w_data, test_delta_data,
                    test_z_c_data, test_z_d_data,
                    beta,
                    alpha1_star_r,
                    tau1,
                    w_min, w_max){
  n = length(test_y_data)
  ## ----- Handle z_c_data and z_d_data as matrices -----
  if (is.null(test_z_c_data)) {
    z_c_mat <- matrix(0, nrow = n, ncol = 0)
    p_c     <- 0
  } else if (is.null(dim(test_z_c_data))) {
    if (length(test_z_c_data) != n) stop("z_c_data must have length n")
    z_c_mat <- matrix(as.numeric(test_z_c_data), ncol = 1)
    p_c     <- 1
  } else {
    if (nrow(test_z_c_data) != n) stop("z_c_data must have n rows")
    z_c_mat <- as.matrix(test_z_c_data)
    p_c     <- ncol(z_c_mat)
  }
  
  if (is.null(test_z_d_data)) {
    z_d_mat <- matrix(0, nrow = n, ncol = 0)
    p_d     <- 0
  } else if (is.null(dim(test_z_d_data))) {
    if (length(test_z_d_data) != n) stop("z_d_data must have length n")
    z_d_mat <- matrix(test_z_d_data, ncol = 1)
    p_d     <- 1
  } else {
    if (nrow(test_z_d_data) != n) stop("z_d_data must have n rows")
    z_d_mat <- as.matrix(test_z_d_data)
    p_d     <- ncol(z_d_mat)
  }
  mean(sapply(1:n, function(i){
    r(test_y_data[i], test_w_data[i], test_delta_data[i],
      beta,
      alpha1_star_r, tau1,
      z_c_mat[i,], z_d_mat[i,],
      w_min, w_max)}) <= zeta)
}


#' Build a half-length search grid from observed residuals
#'
#' Data-in replacement for the simulation grid builder: forms the search grid
#' for zeta from the empirical (1 - alpha) residual quantile of the observed
#' data, widened by a relative margin. No data generation is performed.
#'
#' @param alpha Miscoverage level.
#' @param r Residual function (r1 or r2).
#' @param y_data,w_data,delta_data Observed outcome, time, event indicator.
#' @param z_c_data,z_d_data Optional covariates.
#' @param beta Fitted coefficients.
#' @param alpha1_star_r Center used inside the residual.
#' @param tau1 Standard deviation of X.
#' @param w_min,w_max Truncation bounds.
#' @param margin Relative half-width of the grid around the empirical quantile.
#' @param length Number of grid points.
#' @return A numeric vector of candidate half-lengths.
#' @noRd
get_zeta_seq_r <- function(alpha, r,
                           y_data, w_data, delta_data,
                           z_c_data = NULL, z_d_data = NULL,
                           beta, alpha1_star_r, tau1,
                           w_min = -1, w_max = 1,
                           margin = 0.1, length = 5) {
  n <- length(y_data)
  z_c_mat <- if (is.null(z_c_data)) matrix(0, n, 0) else as.matrix(z_c_data)
  z_d_mat <- if (is.null(z_d_data)) matrix(0, n, 0) else as.matrix(z_d_data)
  res <- vapply(seq_len(n), function(i) {
    zc <- if (ncol(z_c_mat)) z_c_mat[i, ] else NULL
    zd <- if (ncol(z_d_mat)) z_d_mat[i, ] else NULL
    r(y_data[i], w_data[i], delta_data[i], beta, alpha1_star_r, tau1,
      zc, zd, w_min, w_max)
  }, numeric(1))
  q <- stats::quantile(res, 1 - alpha, names = FALSE)
  seq(q * (1 - margin), q * (1 + margin), length.out = length)
}

#' Semiparametric (SPARCC) prediction half-length, end to end
#'
#' Convenience wrapper that fits the nuisance models for X | Z and C | Z,
#' estimates beta by SPARCC, builds the b-arrays, and solves for the efficient
#' half-length zeta -- all from a single training data set with optional
#' covariates. Any pre-fitted piece can be supplied to skip its step.
#'
#' @param y_data,w_data,delta_data Training outcome, observed time, indicator.
#' @param z_c_data,z_d_data Optional continuous / discrete covariates.
#' @param residual Residual function to use (default \code{r1}). Pass
#'   \code{r2} for the r2 residual; the \code{r1*} residual is \code{r1} with a
#'   deliberately wrong \code{alpha1_star_r}.
#' @param residual_name Label recorded in the result and used to name
#'   \code{zeta}. Defaults to the deparsed \code{residual} argument; set it to
#'   \code{"r1star"} when using a misspecified residual center.
#' @param alpha Miscoverage level.
#' @param x_a Grid over the event-time covariate (default seq(w_min, w_max)).
#' @param alpha1_star,tau1 Optional pre-fitted eta1 mean / SD; fitted by
#'   \code{find_alpha1_MLE} when NULL.
#' @param alpha2_star,tau2 Optional pre-fitted eta2 mean / SD; fitted by
#'   \code{find_alpha2_MLE} when NULL.
#' @param beta,sigma Optional pre-fitted outcome coefficients / SD; fitted by
#'   \code{find_beta_sparcc} when NULL.
#' @param alpha1_star_r Center used inside the residual (defaults to
#'   \code{alpha1_star}).
#' @param tau1_r SD used inside the residual (defaults to \code{tau1}).
#' @param test_y_data,test_w_data,test_delta_data Optional test set. When
#'   supplied, \code{coverage_rate} is evaluated on it; otherwise it is
#'   \code{NA}.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param w_min,w_max Truncation bounds.
#' @param tt,m Integration grid size and covariate-grid resolution.
#' @param seq_length Number of half-length grid points.
#' @param b_args Optional named list overriding build_b123_arrays grid sizes.
#' @return A list with \code{method}, \code{alpha}, \code{residual},
#'   \code{zeta} (the half-length, named by the residual),
#'   \code{coverage_rate}, the fitted \code{beta}, \code{sigma} and nuisance
#'   parameters, plus \code{zeta_seq} and \code{zeta_list} from the solve.
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
#' # Runs end to end: fits both working models, fits beta by SPARCC,
#' # then solves for the half-length.
#' fit <- semiparametric_prediction_interval(
#'   y, w, delta,
#'   residual = r1, alpha = 0.1,
#'   w_min = -1, w_max = 1, seq_length = 3,
#'   b_args = list(b1_tt = 5, b1_tt2 = 10, b2_tt = 5, b2_tt2 = 10,
#'                 b3_tt = 5, b3_tt2 = 10, b3_tt3 = 5)
#' )
#' fit$zeta             # half-length, named by the residual
#' fit$coverage_rate    # NA: no test set was supplied
#'
#' # r1* is the same call under a working model for X | Z; label it so the
#' # result records which residual was used.
#' fit_star <- semiparametric_prediction_interval(
#'   y, w, delta,
#'   residual = r1, residual_name = "r1star",
#'   alpha1_star_r = -2, alpha = 0.1,
#'   w_min = -1, w_max = 1, seq_length = 3,
#'   b_args = list(b1_tt = 5, b1_tt2 = 10, b2_tt = 5, b2_tt2 = 10,
#'                 b3_tt = 5, b3_tt2 = 10, b3_tt3 = 5)
#' )
#' fit_star$zeta
#' }
#' @export
semiparametric_prediction_interval <- function(y_data, w_data, delta_data,
                            z_c_data = NULL, z_d_data = NULL,
                            residual = r1, residual_name = NULL,
                            alpha = 0.1,
                            x_a = NULL,
                            alpha1_star = NULL, tau1 = NULL,
                            alpha2_star = NULL, tau2 = NULL,
                            beta = NULL, sigma = NULL,
                            alpha1_star_r = NULL, tau1_r = NULL,
                            test_y_data = NULL,
                            test_w_data = NULL,
                            test_delta_data = NULL,
                            test_z_c_data = NULL,
                            test_z_d_data = NULL,
                            w_min = -1, w_max = 1,
                            tt = 20, m = 20, seq_length = 5,
                            b_args = list()) {
  if (is.null(residual_name)) residual_name <- deparse(substitute(residual))
  if (is.null(x_a)) x_a <- seq(w_min, w_max, length.out = m)

  if (is.null(alpha1_star) || is.null(tau1)) {
    fit1 <- find_alpha1_MLE(w_data, delta_data, z_c_data, z_d_data, w_min, w_max)
    alpha1_star <- fit1[-length(fit1)]; tau1 <- fit1[length(fit1)]
  }
  if (is.null(alpha2_star) || is.null(tau2)) {
    fit2 <- find_alpha2_MLE(w_data, delta_data, z_c_data, z_d_data, w_min, w_max)
    alpha2_star <- fit2[-length(fit2)]; tau2 <- fit2[length(fit2)]
  }
  if (is.null(beta) || is.null(sigma)) {
    fitb <- find_beta_sparcc(y_data, w_data, delta_data,
                                      z_c_data, z_d_data,
                                      alpha1_star, alpha2_star, tau1, tau2,
                                      tt = tt, m = m, w_min = w_min, w_max = w_max)
    beta <- fitb$beta_hat; sigma <- fitb$sigma_hat
  }
  if (is.null(alpha1_star_r)) alpha1_star_r <- alpha1_star
  if (is.null(tau1_r)) tau1_r <- tau1

  zeta_seq <- get_zeta_seq_r(alpha, residual, y_data, w_data, delta_data,
                             z_c_data, z_d_data, beta, alpha1_star_r, tau1,
                             w_min, w_max, length = seq_length)

  b <- do.call(build_b123_arrays, c(list(
    zeta_seq = zeta_seq, alpha = alpha, r = residual, beta = beta, x_a = x_a,
    z_c_data = z_c_data, z_d_data = z_d_data,
    alpha1_star = alpha1_star, alpha2_star = alpha2_star,
    alpha1_star_r = alpha1_star_r, sigma = sigma, tau1 = tau1, tau2 = tau2,
    w_min = w_min, w_max = w_max), b_args))

  zl <- find_zeta_param_int(y_data, w_data, delta_data, z_c_data, z_d_data,
                            b$b1_array, b$b2_array, b$b3_array,
                            beta_temp = beta, zeta_seq = zeta_seq,
                            alpha1_star = alpha1_star, alpha2_star = alpha2_star,
                            sigma = sigma, tau1 = tau1, tau2 = tau2,
                            m = m, w_min = w_min, w_max = w_max)

  zeta <- stats::setNames(zl$sol, residual_name)

  coverage_rate <- .coverage_on_test(
    zeta = zeta, beta = beta,
    test_y_data = test_y_data, test_w_data = test_w_data,
    test_delta_data = test_delta_data,
    test_z_c_data = test_z_c_data, test_z_d_data = test_z_d_data,
    alpha1 = alpha1_star_r, tau1 = tau1,
    alpha1_star_r = alpha1_star_r, tau1_r = tau1_r,
    w_min = w_min, w_max = w_max
  )

  list(method        = "semiparametric",
       alpha         = alpha,
       residual      = residual_name,
       zeta          = zeta,
       coverage_rate = coverage_rate,
       beta          = beta, sigma = sigma,
       alpha1_star   = alpha1_star, tau1 = tau1,
       alpha2_star   = alpha2_star, tau2 = tau2,
       alpha1_star_r = alpha1_star_r,
       zeta_seq      = zeta_seq, zeta_list = zl)
}
