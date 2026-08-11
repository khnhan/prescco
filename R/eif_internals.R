# =============================================================================
#  Efficient influence-function internals (covariate-aware, data-in)
#  ---------------------------------------------------------------------------
#  Low-level building blocks of the SPARCC efficient score used by the
#  semiparametric estimator and interval. These operate on given data and
#  support optional covariates via the phi_xz basis. They are internal helpers
#  (not exported); see find_beta_sparcc() and
#  semiparametric_prediction_interval() for the user-facing entry points.
#
#  Package imports are declared centrally in R/censcovpred-package.R.
# =============================================================================

#' Internal SPARCC efficient-score helper (covariate-aware).
#' @noRd
gauss = function(tt, len = 3){
  grid = seq(-len, len, length.out = tt)
  d = dnorm(grid)
  list(x = grid,
       w = d / sum(d))
}
#' Internal SPARCC efficient-score helper (covariate-aware).
#' @noRd
S_beta_f = function(beta, y, x, z_c, z_d, sigma){
  # Ensure these are numeric vectors
  z_c = as.numeric(z_c)
  z_d = as.numeric(z_d)
  
  # Design vector for regression part:
  # (1, x, z_c, z_d, x*z_c, x*z_d)
  design = c(1,
             x,
             z_c,
             z_d,
             x * z_c,
             x * z_d)
  
  # beta must have length equal to length(design)
  # (sigma^2 is *not* part of beta; it comes as `sigma`)
  mu  = sum(beta * design)
  
  res = y - mu
  
  # score wrt (all beta's, then sigma^2)
  c(
    res * design,
    res^2 / (2 * sigma^2) - 1 / (2)
  )
}
#' Internal SPARCC efficient-score helper (covariate-aware).
#' @noRd
S_beta = function(beta, y, w, delta, alpha1_star,
                  z_c, z_d,
                  sigma, tau1,
                  w_min, w_max,
                  n_grid = 20){
  
  # Make sure z_c, z_d are numeric vectors
  z_c = as.numeric(z_c)
  z_d = as.numeric(z_d)
  p_c = length(z_c)
  p_d = length(z_d)
  
  ## Extract coefficients according to design:
  ## design = c(1, x, z_c, z_d, x*z_c, x*z_d)
  idx_beta0   = 1
  idx_betax   = 2
  idx_betazc  = if (p_c > 0) 3:(2 + p_c) else integer(0)
  idx_betazd  = if (p_d > 0) (3 + p_c):(2 + p_c + p_d) else integer(0)
  idx_betaxzc = if (p_c > 0) (3 + p_c + p_d):(2 + 2*p_c + p_d) else integer(0)
  idx_betaxzd = if (p_d > 0) (3 + 2*p_c + p_d):(2 + 2*p_c + 2*p_d) else integer(0)
  
  beta0   = beta[idx_beta0]
  beta_x  = beta[idx_betax]
  beta_zc  = if (length(idx_betazc)  > 0) beta[idx_betazc]  else numeric(0)
  beta_zd  = if (length(idx_betazd)  > 0) beta[idx_betazd]  else numeric(0)
  beta_xzc = if (length(idx_betaxzc) > 0) beta[idx_betaxzc] else numeric(0)
  beta_xzd = if (length(idx_betaxzd) > 0) beta[idx_betaxzd] else numeric(0)
  
  ## Effective intercept and slope in X given (z_c, z_d)
  intercept_x = beta0 +
    (if (p_c > 0) sum(beta_zc  * z_c) else 0) +
    (if (p_d > 0) sum(beta_zd  * z_d) else 0)
  
  slope_x = beta_x +
    (if (p_c > 0) sum(beta_xzc * z_c) else 0) +
    (if (p_d > 0) sum(beta_xzd * z_d) else 0)
  
  ## Posterior var and mean of X | Y, Z (Gaussian–Gaussian conjugacy)
  v_x  = 1 / (slope_x^2 / sigma^2 + 1 / tau1^2)
  eta_x = v_x * (slope_x * (y - intercept_x) / sigma^2 +
                   sum(c(1, z_c, z_d) * alpha1_star) / tau1^2)
  
  if (delta == 0) {
    ## Censored: integrate over X > max(W, w_min) up to w_max
    lower = max(w, w_min)
    
    if (lower >= w_max) {
      # Region {X > W} is essentially empty; any value is fine
      return(as.vector(S_beta_f(beta, y, w, z_c, z_d, sigma)))
    } else {
      x_grid = seq(lower, w_max, length.out = n_grid)
      x_std  = (x_grid - eta_x) / sqrt(v_x)
      dens   = dnorm(x_std)
      
      # S_beta_f is scalar in x: evaluate it at each grid point to build a
      # (len_param x n_grid) matrix, then integrate against dens.
      S_mat = vapply(x_grid,
                     function(xx) S_beta_f(beta, y, xx, z_c, z_d, sigma),
                     numeric(length(beta) + 1L))
      num   = S_mat %*% dens
      denom = sum(dens)
      
      return(as.vector(num / denom))
    }
  } else {
    ## Uncensored
    return(as.vector(S_beta_f(beta, y, w, z_c, z_d, sigma)))
  }
}
#' Internal SPARCC efficient-score helper (covariate-aware).
#' @noRd
c0_xz_gauss_param12 = function(beta, x_a, z_c, z_d,
                               alpha1_star, alpha2_star,
                               sigma, tau1, tau2,
                               w_min, w_max,
                               tt = 20) {
  cc  = gauss(tt)   # Gaussian quadrature nodes/weights (standard normal)
  len = 20          # grid size for integrating over C
  
  # Make sure z_c, z_d are numeric vectors
  z_c = as.numeric(z_c)
  z_d = as.numeric(z_d)
  
  ## temp2: E[ I(X < C) S_beta(Y, X, Z) + I(X >= C) S_beta(Y, C, Z) | Y, X, Z ]
  temp2 = function(x, y) {
    if (x < w_min) {
      ## If X is below support, just evaluate score at (X = x)
      return(S_beta_f(beta, y, x, z_c, z_d, sigma))
    } else {
      ## Integrate over C in [w_min, min(x, w_max)]
      c_grid = seq(w_min, min(x, w_max), length.out = len)
      
      ## C | Z ~ TruncNormal( mean(1, z_c, z_d), tau2^2 ) on [w_min, w_max]
      mean_C = sum(c(1, z_c, z_d) * alpha2_star)
      dens_c = truncnorm::dtruncnorm(c_grid,
                          a    = w_min,
                          b    = w_max,
                          mean = mean_C,
                          sd   = tau2)
      
      ## S_beta at each C-grid point (delta = 0 case, censored)
      Sbeta = sapply(c_grid, function(w)
        S_beta(beta, y, w, 0, alpha1_star,
               z_c, z_d,
               sigma, tau1,
               w_min, w_max))
      
      ## Trapezoidal correction at endpoints
      Sbeta[, c(1, len)] = Sbeta[, c(1, len)] / 2
      by = c_grid[2] - c_grid[1]
      
      ## P(C > x | Z) for the uncensored part
      tail_prob = 1 - truncnorm::ptruncnorm(x,
                                 a    = w_min,
                                 b    = w_max,
                                 mean = mean_C,
                                 sd   = tau2)
      
      ## E[ I(X < C) S_beta(Y, X, Z) + I(X >= C) S_beta(Y, C, Z) | Y, X, Z ]
      res = S_beta_f(beta, y, x, z_c, z_d, sigma) * tail_prob +
        Sbeta %*% (dens_c * by)
      
      as.vector(res)
    }
  }
  
  ## temp3: E[ temp2(Y, X, Z) | X, Z ]
  temp3 = function(x) {
    ## Mean of Y | X = x, Z = (z_c, z_d)
    design_x = c(1,
                 x,
                 z_c,
                 z_d,
                 x * z_c,
                 x * z_d)
    mu_xz = sum(beta * design_x)
    
    ## Quadrature nodes for Y: Y = mu_xz + sigma * Z, Z ~ N(0,1)
    y_grid = mu_xz + sigma * cc$x
    
    ## Integrand: temp2(x, y) at each y-grid
    vals = sapply(y_grid, function(y_norm) temp2(x, y_norm))
    
    ## Integrate over Y with quadrature weights
    as.vector(vals %*% cc$w)
  }
  
  temp3 = Vectorize(temp3, vectorize.args = "x")
  temp3(x_a)
}
#' Internal SPARCC efficient-score helper (covariate-aware).
#' @noRd
L_xz_gauss_param12 = function(beta, x_a, z_c, z_d,
                              alpha1_star, alpha2_star,
                              sigma, tau1, tau2,
                              w_min, w_max,
                              tt = 20){
  cc = gauss(tt)
  
  # Ensure covariates are numeric vectors
  z_c = as.numeric(z_c)
  z_d = as.numeric(z_d)
  p_c = length(z_c)
  p_d = length(z_d)
  
  ## Precompute intercept and slope in X for the linear model
  ## mu(x, z) = beta' * design(x, z)
  ## design(x, z) = (1, x, z_c, z_d, x*z_c, x*z_d)
  design_x0   = c(1, 0,              z_c,          z_d,          0 * z_c, 0 * z_d)
  d_design_dx = c(0, 1, rep(0, p_c + p_d),        z_c,          z_d)
  
  intercept_x = sum(beta * design_x0)
  slope_x     = sum(beta * d_design_dx)
  
  ## temp[a(X,Z)] = E[ I(X > C) a(X,Z) | C, Y, Z ] / E[ I(X > C) | C, Y, Z ]
  temp = function(x_a, c, y){
    # Posterior X | Y, Z parameters
    v_x  = 1 / (slope_x^2 / sigma^2 + 1 / tau1^2)
    eta_x = v_x * (slope_x * (y - intercept_x) / sigma^2 +
                     sum(c(1, z_c, z_d) * alpha1_star) / tau1^2)
    
    # f_{X|Y,Z} truncated to [w_min, w_max] evaluated on grid x_a
    p = truncnorm::dtruncnorm(x_a,
                   a    = w_min,
                   b    = w_max,
                   mean = eta_x,
                   sd   = sqrt(v_x))
    
    num   = (x_a > c) * p
    denom = sum(num)
    
    ifelse(is.nan(num / denom),
           0, num / denom)
  }
  
  ## temp2[a(X,Z)] = E[ I(X > C) temp(x_a, C, Y, Z) | X, Y, Z ]
  temp2 = function(x_a, x, y){
    c_grid = seq(w_min, w_max, length = 20)
    
    mean_C = sum(c(1, z_c, z_d) * alpha2_star)
    dens = truncnorm::dtruncnorm(c_grid,
                      a    = w_min,
                      b    = w_max,
                      mean = mean_C,
                      sd   = tau2)
    dens = dens / sum(dens)
    
    sapply(c_grid, function(c_norm){
      (x > c_norm) * temp(x_a, c_norm, y)
    }) %*% dens
  }
  
  ## temp3[a(X,Z)] = E[ temp2(X, Y, Z) | X, Z ]
  temp3 = function(x_a, x){
    # mean of Y | X = x, Z = (z_c, z_d)
    design_x = c(1,
                 x,
                 z_c,
                 z_d,
                 x * z_c,
                 x * z_d)
    mu_xz = sum(beta * design_x)
    
    y_grid = mu_xz + sigma * cc$x
    
    sapply(y_grid, function(y_norm){
      temp2(x_a, x, y_norm)
    }) %*% cc$w
  }
  temp3 = Vectorize(temp3, vectorize.args = "x")
  
  ## temp4[a(X,Z)] = E[ I(X < C) | Y, X, Z ] a(X,Z)
  temp4 = function(x, y){
    c_grid = seq(w_min, w_max, length = 20)
    
    mean_C = sum(c(1, z_c, z_d) * alpha2_star)
    dens = truncnorm::dtruncnorm(c_grid,
                      a    = w_min,
                      b    = w_max,
                      mean = mean_C,
                      sd   = tau2)
    dens = dens / sum(dens)
    
    # P(C >= x | Z) via discrete approximation
    sum((x <= c_grid) * dens)
  }
  temp4 = Vectorize(temp4, vectorize.args = "y")
  
  ## temp5[a(X,Z)] = E[ temp4(X, Y, Z) | X, Z ] a(X,Z)
  temp5 = function(x){
    design_x = c(1,
                 x,
                 z_c,
                 z_d,
                 x * z_c,
                 x * z_d)
    mu_xz = sum(beta * design_x)
    
    y_grid = mu_xz + sigma * cc$x
    
    sum(temp4(x, y_grid) * cc$w)
  }
  temp5 = Vectorize(temp5, vectorize.args = "x")
  
  res = diag(temp5(x_a)) + t(temp3(x_a, x_a))
  res
}
#' Internal SPARCC efficient-score helper (covariate-aware).
#' @noRd
a_gauss_param12 = function(beta, x_a,
                           z_c, z_d,
                           alpha1_star, alpha2_star,
                           sigma, tau1, tau2,
                           tt = 20,
                           w_min = 0, w_max = 12) {
  
  # L(x_a, z_c, z_d): m x m matrix, m = length(x_a)
  L_mat = L_xz_gauss_param12(
    beta      = beta,
    x_a       = x_a,
    z_c       = z_c,
    z_d       = z_d,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    sigma     = sigma,
    tau1      = tau1,
    tau2      = tau2,
    w_min     = w_min,
    w_max     = w_max,
    tt        = tt
  )
  
  # b(x_a, z_c, z_d): row vector of length m
  b_row = c0_xz_gauss_param12(
    beta      = beta,
    x_a       = x_a,
    z_c       = z_c,
    z_d       = z_d,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    sigma     = sigma,
    tau1      = tau1,
    tau2      = tau2,
    w_min     = w_min,
    w_max     = w_max,
    tt        = tt
  )
  
  # a(x_a, z_c, z_d) = L^{-1} b^T  → length(x_a) vector
  a_vec = MASS::ginv(L_mat) %*% t(b_row)
  a_vec
}
#' Internal SPARCC efficient-score helper (covariate-aware).
#' @noRd
S_eff_gauss_param12 = function(beta, y, w, delta,
                               x_a, a0,
                               z_c, z_d,
                               alpha1_star,
                               sigma, tau1, tau2,
                               w_min, w_max,
                               tt = 20) {
  
  ## Number of parameters: all betas + sigma^2 score
  len_beta = length(beta) + 1
  
  # Make sure covariates are numeric vectors
  z_c = as.numeric(z_c)
  z_d = as.numeric(z_d)
  p_c = length(z_c)
  p_d = length(z_d)
  
  ## Precompute intercept and slope in X for the linear model
  ## mu(x, z) = beta' * design(x, z)
  ## design(x, z) = (1, x, z_c, z_d, x*z_c, x*z_d)
  design_x0   = c(1, 0,              z_c,          z_d,          0 * z_c, 0 * z_d)
  d_design_dx = c(0, 1, rep(0, p_c + p_d),        z_c,          z_d)
  
  intercept_x = sum(beta * design_x0)
  slope_x     = sum(beta * d_design_dx)
  
  if (delta == 0) {
    # temp[a(X,Z)] = E[ I(X > C) a(X,Z) | C, Y, Z ] / E[ I(X > C) | C, Y, Z ]
    temp = function(x_a, c, y) {
      # Posterior X | Y, Z parameters under full interaction model
      v_x  = 1 / (slope_x^2 / sigma^2 + 1 / tau1^2)
      eta_x = v_x * (slope_x * (y - intercept_x) / sigma^2 +
                       sum(c(1, z_c, z_d) * alpha1_star) / tau1^2)
      
      # Discretization of f(X | Y, Z) on x_a, truncated to [w_min, w_max]
      p = truncnorm::dtruncnorm(x_a,
                     a    = w_min,
                     b    = w_max,
                     mean = eta_x,
                     sd   = sqrt(v_x))
      
      num   = (x_a > c) * p
      denom = sum(num)
      
      ifelse(is.nan(num / denom), 0, num / denom)
    }
    
    # Full score S_beta (including sigma^2 component)
    sbeta = S_beta(beta, y, w, 0,
                   alpha1_star,
                   z_c, z_d,
                   sigma, tau1,
                   w_min, w_max)
    
    # S_eff = S_beta - A0 * temp
    return(sbeta - as.vector(t(a0) %*% temp(x_a, w, y)))
    
  } else {  # delta == 1, uncensored: interpolate a0 at w
    m    = length(x_a)
    a0_w = numeric(len_beta)
    
    for (j in 1:len_beta) {
      a0_w[j] = approx(x_a, a0[, j], w, rule = 2)$y  # linear interpolation
    }
    
    # Full score at observed X = w
    sbeta = S_beta_f(beta, y, w, z_c, z_d, sigma)
    
    # S_eff = S_beta - a0(w)
    return(sbeta - a0_w)
  }
}
#' Internal SPARCC efficient-score helper (covariate-aware).
#' @noRd
pe_gauss_param12 = function(betasigma,
                            y_data, w_data, delta_data,
                            x_a,
                            z_c_data, z_d_data,
                            alpha1_star, alpha2_star,
                            tau1, tau2,
                            w_min, w_max,
                            tt = 20) {
  
  n        <- length(y_data)
  p_theta  <- length(betasigma)
  p_beta   <- p_theta - 1L          # all but last are beta
  beta     <- betasigma[1:p_beta]
  log_sig  <- betasigma[p_theta]
  sigma    <- exp(log_sig)          # enforce sigma > 0
  
  len_beta <- p_beta + 1L           # beta components + sigma^2 score
  x_len    <- length(x_a)
  
  ## ----- Handle z_c_data and z_d_data as matrices -----
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
  
  ## ======================================================
  ## Precompute a0_array over (z_c grid, z_d levels)
  ## a0 is an x_len x len_beta matrix for each (z_c, z_d)
  ## ------------------------------------------------------
  ## You already have this logic; here is the structure:
  
  n_c_grid <- 10  # or your choice
  
  # Continuous grids
  if (p_c > 0) {
    z_c_grid_list <- lapply(1:p_c, function(j) {
      range_j <- range(z_c_mat[, j])
      seq(range_j[1], range_j[2], length.out = n_c_grid)
    })
  } else {
    z_c_grid_list <- list()
  }
  
  # Discrete levels per z_d dimension
  if (p_d > 0) {
    levels_d_list <- lapply(1:p_d, function(j) sort(unique(z_d_mat[, j])))
    n_d_levels    <- vapply(levels_d_list, length, integer(1))
  } else {
    levels_d_list <- list()
    n_d_levels    <- integer(0)
  }
  
  n_zc_comb <- if (p_c > 0) n_c_grid^p_c else 1L
  n_zd_comb <- if (p_d > 0) prod(n_d_levels) else 1L
  
  # a0_array: (x, parameter, z_c-comb, z_d-comb)
  a0_array <- array(
    NA_real_,
    dim = c(x_len, len_beta, n_zc_comb, n_zd_comb)
  )
  
  ## Helpers to map indices <-> discrete z_d patterns
  z_d_from_index <- function(l) {
    if (p_d == 0) return(numeric(0))
    tmp   <- l - 1L
    idx_d <- integer(p_d)
    for (j in 1:p_d) {
      base_j   <- n_d_levels[j]
      idx_d[j] <- (tmp %% base_j) + 1L
      tmp      <- tmp %/% base_j
    }
    vapply(1:p_d, function(j) levels_d_list[[j]][idx_d[j]], numeric(1))
  }
  
  index_from_zd <- function(z_d_i) {
    if (p_d == 0) return(1L)
    if (length(z_d_i) != p_d) stop("z_d_i has wrong length")
    
    idx_d <- integer(p_d)
    for (j in 1:p_d) {
      idx_j <- match(z_d_i[j], levels_d_list[[j]])
      if (is.na(idx_j)) stop("z_d_i[j] not among observed levels")
      idx_d[j] <- idx_j
    }
    
    l    <- 1L
    mult <- 1L
    for (j in 1:p_d) {
      l    <- l + (idx_d[j] - 1L) * mult
      mult <- mult * n_d_levels[j]
    }
    l
  }
  
  # Linear index for z_c grid multi-index
  lin_index_zc <- function(k_vec) {
    if (p_c == 0) return(1L)
    idx  <- 1L
    mult <- 1L
    for (j in 1:p_c) {
      idx  <- idx + (k_vec[j] - 1L) * mult
      mult <- mult * n_c_grid
    }
    idx
  }
  
  # Precompute a0_array
  if (p_c > 0) {
    k_c_vec <- rep(1L, p_c)
  } else {
    k_c_vec <- integer(0)
  }
  
  for (r in 1:n_zc_comb) {
    if (p_c > 0) {
      z_c_point <- vapply(
        1:p_c,
        function(j) z_c_grid_list[[j]][k_c_vec[j]],
        numeric(1)
      )
    } else {
      z_c_point <- numeric(0)
    }
    
    for (l in 1:n_zd_comb) {
      z_d_point <- z_d_from_index(l)
      
      a0_mat <- a_gauss_param12(
        beta        = beta,
        x_a         = x_a,
        z_c         = z_c_point,
        z_d         = z_d_point,
        alpha1_star = alpha1_star,
        alpha2_star = alpha2_star,
        sigma       = sigma,
        tau1        = tau1,
        tau2        = tau2,
        tt          = tt,
        w_min       = w_min,
        w_max       = w_max
      )
      
      if (!is.matrix(a0_mat) ||
          !all(dim(a0_mat) == c(x_len, len_beta))) {
        stop("a_gauss_param12 must return an x_len x len_beta matrix.")
      }
      
      a0_array[, , r, l] <- a0_mat
    }
    
    if (p_c > 0) {
      for (j in 1:p_c) {
        k_c_vec[j] <- k_c_vec[j] + 1L
        if (k_c_vec[j] <= n_c_grid) break
        k_c_vec[j] <- 1L
      }
    }
  }
  
  ## Interpolator: returns a0(x,·) matrix for given (z_c_i,z_d_i)
  get_a0_for_obs <- function(z_c_i, z_d_i) {
    l <- index_from_zd(z_d_i)
    
    if (p_c == 0) {
      return(a0_array[, , 1, l])
    }
    
    k_low  <- integer(p_c)
    k_high <- integer(p_c)
    t_vec  <- numeric(p_c)
    
    for (j in 1:p_c) {
      grid_j <- z_c_grid_list[[j]]
      z_ij   <- z_c_i[j]
      
      if (z_ij <= grid_j[1]) {
        k_low[j]  <- 1L
        k_high[j] <- 1L
        t_vec[j]  <- 0
      } else if (z_ij >= grid_j[n_c_grid]) {
        k_low[j]  <- n_c_grid
        k_high[j] <- n_c_grid
        t_vec[j]  <- 0
      } else {
        k_low[j]  <- max(which(grid_j <= z_ij))
        k_high[j] <- k_low[j] + 1L
        t_vec[j]  <- (z_ij - grid_j[k_low[j]]) /
          (grid_j[k_high[j]] - grid_j[k_low[j]])
      }
    }
    
    res    <- matrix(0, nrow = x_len, ncol = len_beta)
    n_vert <- 2^p_c
    
    for (v in 0:(n_vert - 1L)) {
      k_vec_v <- integer(p_c)
      weight  <- 1
      
      for (j in 1:p_c) {
        bit_j <- (v %/% 2^(j - 1L)) %% 2L
        if (bit_j == 0L) {
          k_vec_v[j] <- k_low[j]
          wj         <- 1 - t_vec[j]
        } else {
          k_vec_v[j] <- k_high[j]
          wj         <- t_vec[j]
        }
        weight <- weight * wj
      }
      
      r_idx <- lin_index_zc(k_vec_v)
      res   <- res + weight * a0_array[, , r_idx, l]
    }
    
    res
  }
  
  ## ======================================================
  ## Main estimating equation: sum_i S_eff(beta,sigma^2; data_i) = 0
  ## ------------------------------------------------------
  val_full <- rep(0, len_beta)
  
  for (i in 1:n) {
    z_c_i <- if (p_c > 0) z_c_mat[i, ] else numeric(0)
    z_d_i <- if (p_d > 0) z_d_mat[i, ] else numeric(0)
    
    a0_i <- get_a0_for_obs(z_c_i, z_d_i)
    
    S_i <- S_eff_gauss_param12(
      beta       = beta,
      y          = y_data[i],
      w          = w_data[i],
      delta      = delta_data[i],
      x_a        = x_a,
      a0         = a0_i,
      z_c        = z_c_i,
      z_d        = z_d_i,
      alpha1_star = alpha1_star,
      sigma      = sigma,
      tau1       = tau1,
      tau2       = tau2,
      w_min      = w_min,
      w_max      = w_max,
      tt         = tt
    )
    
    val_full <- val_full + S_i  # length len_beta
  }
  
  ## IMPORTANT: length(val_full) == length(betasigma),
  ## so nleqslv sees matching dimensions.
  return(val_full)
}
