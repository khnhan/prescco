# =============================================================================
# Efficient influence function calculations for SPARCC
#
# Internal numerical functions used by find_beta_sparcc().
# These functions support optional continuous and discrete covariates.
# =============================================================================


# Gaussian quadrature grid.

gauss <- function(tt, len = 3) {
  grid <- seq(-len, len, length.out = tt)
  d <- dnorm(grid)

  list(
    x = grid,
    w = d / sum(d)
  )
}


# Complete-data score for beta and sigma^2.

S_beta_f <- function(beta, y, x, z_c, z_d, sigma) {

  z_c <- as.numeric(z_c)
  z_d <- as.numeric(z_d)

  design <- c(
    1,
    x,
    z_c,
    z_d,
    x * z_c,
    x * z_d
  )

  mu <- sum(beta * design)
  res <- y - mu

  c(
    res * design,
    res^2 / (2 * sigma^2) - 1 / 2
  )
}


# Observed-data score.

S_beta <- function(beta, y, w, delta, alpha1_star,
                   z_c, z_d,
                   sigma, tau1,
                   w_min, w_max,
                   n_grid = 20) {

  z_c <- as.numeric(z_c)
  z_d <- as.numeric(z_d)

  p_c <- length(z_c)
  p_d <- length(z_d)

  idx_beta0 <- 1
  idx_betax <- 2

  idx_betazc <- if (p_c > 0) {
    3:(2 + p_c)
  } else {
    integer(0)
  }

  idx_betazd <- if (p_d > 0) {
    (3 + p_c):(2 + p_c + p_d)
  } else {
    integer(0)
  }

  idx_betaxzc <- if (p_c > 0) {
    (3 + p_c + p_d):(2 + 2 * p_c + p_d)
  } else {
    integer(0)
  }

  idx_betaxzd <- if (p_d > 0) {
    (3 + 2 * p_c + p_d):(2 + 2 * p_c + 2 * p_d)
  } else {
    integer(0)
  }

  beta0 <- beta[idx_beta0]
  beta_x <- beta[idx_betax]

  beta_zc <- if (length(idx_betazc) > 0) {
    beta[idx_betazc]
  } else {
    numeric(0)
  }

  beta_zd <- if (length(idx_betazd) > 0) {
    beta[idx_betazd]
  } else {
    numeric(0)
  }

  beta_xzc <- if (length(idx_betaxzc) > 0) {
    beta[idx_betaxzc]
  } else {
    numeric(0)
  }

  beta_xzd <- if (length(idx_betaxzd) > 0) {
    beta[idx_betaxzd]
  } else {
    numeric(0)
  }

  intercept_x <- beta0 +
    (if (p_c > 0) sum(beta_zc * z_c) else 0) +
    (if (p_d > 0) sum(beta_zd * z_d) else 0)

  slope_x <- beta_x +
    (if (p_c > 0) sum(beta_xzc * z_c) else 0) +
    (if (p_d > 0) sum(beta_xzd * z_d) else 0)

  v_x <- 1 / (
    slope_x^2 / sigma^2 +
      1 / tau1^2
  )

  eta_x <- v_x * (
    slope_x * (y - intercept_x) / sigma^2 +
      sum(c(1, z_c, z_d) * alpha1_star) / tau1^2
  )

  if (delta == 0) {

    lower <- max(w, w_min)

    if (lower >= w_max) {

      return(
        as.vector(
          S_beta_f(
            beta,
            y,
            w,
            z_c,
            z_d,
            sigma
          )
        )
      )
    }

    x_grid <- seq(
      lower,
      w_max,
      length.out = n_grid
    )

    x_std <- (x_grid - eta_x) / sqrt(v_x)
    dens <- dnorm(x_std)

    S_mat <- vapply(
      x_grid,
      function(xx) {
        S_beta_f(
          beta,
          y,
          xx,
          z_c,
          z_d,
          sigma
        )
      },
      numeric(length(beta) + 1L)
    )

    num <- S_mat %*% dens
    denom <- sum(dens)

    return(
      as.vector(
        num / denom
      )
    )
  }

  as.vector(
    S_beta_f(
      beta,
      y,
      w,
      z_c,
      z_d,
      sigma
    )
  )
}


# Right-hand side of the integral equation.

c0_xz_gauss_param12 <- function(beta, x_a, z_c, z_d,
                                alpha1_star, alpha2_star,
                                sigma, tau1, tau2,
                                w_min, w_max,
                                tt = 20) {

  cc <- gauss(tt)
  len <- 20

  z_c <- as.numeric(z_c)
  z_d <- as.numeric(z_d)

  temp2 <- function(x, y) {

    if (x < w_min) {

      return(
        S_beta_f(
          beta,
          y,
          x,
          z_c,
          z_d,
          sigma
        )
      )
    }

    c_grid <- seq(
      w_min,
      min(x, w_max),
      length.out = len
    )

    mean_C <- sum(
      c(1, z_c, z_d) *
        alpha2_star
    )

    dens_c <- truncnorm::dtruncnorm(
      c_grid,
      a = w_min,
      b = w_max,
      mean = mean_C,
      sd = tau2
    )

    Sbeta <- sapply(
      c_grid,
      function(w) {
        S_beta(
          beta,
          y,
          w,
          0,
          alpha1_star,
          z_c,
          z_d,
          sigma,
          tau1,
          w_min,
          w_max
        )
      }
    )

    Sbeta[, c(1, len)] <- Sbeta[, c(1, len)] / 2
    by <- c_grid[2] - c_grid[1]

    tail_prob <- 1 -
      truncnorm::ptruncnorm(
        x,
        a = w_min,
        b = w_max,
        mean = mean_C,
        sd = tau2
      )

    res <- S_beta_f(
      beta,
      y,
      x,
      z_c,
      z_d,
      sigma
    ) * tail_prob +
      Sbeta %*% (dens_c * by)

    as.vector(res)
  }

  temp3 <- function(x) {

    design_x <- c(
      1,
      x,
      z_c,
      z_d,
      x * z_c,
      x * z_d
    )

    mu_xz <- sum(beta * design_x)
    y_grid <- mu_xz + sigma * cc$x

    vals <- sapply(
      y_grid,
      function(y_norm) {
        temp2(
          x,
          y_norm
        )
      }
    )

    as.vector(
      vals %*% cc$w
    )
  }

  temp3 <- Vectorize(
    temp3,
    vectorize.args = "x"
  )

  temp3(x_a)
}


# Linear operator in the integral equation.

L_xz_gauss_param12 <- function(beta, x_a, z_c, z_d,
                               alpha1_star, alpha2_star,
                               sigma, tau1, tau2,
                               w_min, w_max,
                               tt = 20) {

  cc <- gauss(tt)

  z_c <- as.numeric(z_c)
  z_d <- as.numeric(z_d)

  p_c <- length(z_c)
  p_d <- length(z_d)

  design_x0 <- c(
    1,
    0,
    z_c,
    z_d,
    0 * z_c,
    0 * z_d
  )

  d_design_dx <- c(
    0,
    1,
    rep(0, p_c + p_d),
    z_c,
    z_d
  )

  intercept_x <- sum(beta * design_x0)
  slope_x <- sum(beta * d_design_dx)

  temp <- function(x_a, c, y) {

    v_x <- 1 / (
      slope_x^2 / sigma^2 +
        1 / tau1^2
    )

    eta_x <- v_x * (
      slope_x * (y - intercept_x) / sigma^2 +
        sum(c(1, z_c, z_d) * alpha1_star) / tau1^2
    )

    p <- truncnorm::dtruncnorm(
      x_a,
      a = w_min,
      b = w_max,
      mean = eta_x,
      sd = sqrt(v_x)
    )

    num <- (x_a > c) * p
    denom <- sum(num)

    ifelse(
      is.nan(num / denom),
      0,
      num / denom
    )
  }

  temp2 <- function(x_a, x, y) {

    c_grid <- seq(
      w_min,
      w_max,
      length = 20
    )

    mean_C <- sum(
      c(1, z_c, z_d) *
        alpha2_star
    )

    dens <- truncnorm::dtruncnorm(
      c_grid,
      a = w_min,
      b = w_max,
      mean = mean_C,
      sd = tau2
    )

    dens <- dens / sum(dens)

    sapply(
      c_grid,
      function(c_norm) {
        (x > c_norm) *
          temp(
            x_a,
            c_norm,
            y
          )
      }
    ) %*% dens
  }

  temp3 <- function(x_a, x) {

    design_x <- c(
      1,
      x,
      z_c,
      z_d,
      x * z_c,
      x * z_d
    )

    mu_xz <- sum(beta * design_x)
    y_grid <- mu_xz + sigma * cc$x

    sapply(
      y_grid,
      function(y_norm) {
        temp2(
          x_a,
          x,
          y_norm
        )
      }
    ) %*% cc$w
  }

  temp3 <- Vectorize(
    temp3,
    vectorize.args = "x"
  )

  temp4 <- function(x, y) {

    c_grid <- seq(
      w_min,
      w_max,
      length = 20
    )

    mean_C <- sum(
      c(1, z_c, z_d) *
        alpha2_star
    )

    dens <- truncnorm::dtruncnorm(
      c_grid,
      a = w_min,
      b = w_max,
      mean = mean_C,
      sd = tau2
    )

    dens <- dens / sum(dens)

    sum(
      (x <= c_grid) *
        dens
    )
  }

  temp4 <- Vectorize(
    temp4,
    vectorize.args = "y"
  )

  temp5 <- function(x) {

    design_x <- c(
      1,
      x,
      z_c,
      z_d,
      x * z_c,
      x * z_d
    )

    mu_xz <- sum(beta * design_x)
    y_grid <- mu_xz + sigma * cc$x

    sum(
      temp4(
        x,
        y_grid
      ) *
        cc$w
    )
  }

  temp5 <- Vectorize(
    temp5,
    vectorize.args = "x"
  )

  diag(
    temp5(x_a)
  ) +
    t(
      temp3(
        x_a,
        x_a
      )
    )
}


# Solution of the integral equation.

a_gauss_param12 <- function(beta, x_a,
                            z_c, z_d,
                            alpha1_star, alpha2_star,
                            sigma, tau1, tau2,
                            tt = 20,
                            w_min = 0, w_max = 12) {

  L_mat <- L_xz_gauss_param12(
    beta = beta,
    x_a = x_a,
    z_c = z_c,
    z_d = z_d,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2,
    w_min = w_min,
    w_max = w_max,
    tt = tt
  )

  b_row <- c0_xz_gauss_param12(
    beta = beta,
    x_a = x_a,
    z_c = z_c,
    z_d = z_d,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2,
    w_min = w_min,
    w_max = w_max,
    tt = tt
  )

  MASS::ginv(L_mat) %*%
    t(b_row)
}


# Efficient observed-data score.

S_eff_gauss_param12 <- function(beta, y, w, delta,
                                x_a, a0,
                                z_c, z_d,
                                alpha1_star,
                                sigma, tau1, tau2,
                                w_min, w_max,
                                tt = 20) {

  len_beta <- length(beta) + 1

  z_c <- as.numeric(z_c)
  z_d <- as.numeric(z_d)

  p_c <- length(z_c)
  p_d <- length(z_d)

  design_x0 <- c(
    1,
    0,
    z_c,
    z_d,
    0 * z_c,
    0 * z_d
  )

  d_design_dx <- c(
    0,
    1,
    rep(0, p_c + p_d),
    z_c,
    z_d
  )

  intercept_x <- sum(beta * design_x0)
  slope_x <- sum(beta * d_design_dx)

  if (delta == 0) {

    temp <- function(x_a, c, y) {

      v_x <- 1 / (
        slope_x^2 / sigma^2 +
          1 / tau1^2
      )

      eta_x <- v_x * (
        slope_x * (y - intercept_x) / sigma^2 +
          sum(c(1, z_c, z_d) * alpha1_star) / tau1^2
      )

      p <- truncnorm::dtruncnorm(
        x_a,
        a = w_min,
        b = w_max,
        mean = eta_x,
        sd = sqrt(v_x)
      )

      num <- (x_a > c) * p
      denom <- sum(num)

      ifelse(
        is.nan(num / denom),
        0,
        num / denom
      )
    }

    sbeta <- S_beta(
      beta,
      y,
      w,
      0,
      alpha1_star,
      z_c,
      z_d,
      sigma,
      tau1,
      w_min,
      w_max
    )

    return(
      sbeta -
        as.vector(
          t(a0) %*%
            temp(
              x_a,
              w,
              y
            )
        )
    )
  }

  a0_w <- numeric(len_beta)

  for (j in seq_len(len_beta)) {
    a0_w[j] <- approx(
      x_a,
      a0[, j],
      w,
      rule = 2
    )$y
  }

  sbeta <- S_beta_f(
    beta,
    y,
    w,
    z_c,
    z_d,
    sigma
  )

  sbeta - a0_w
}


# Efficient estimating equation for beta and sigma.

pe_gauss_param12 <- function(betasigma,
                             y_data, w_data, delta_data,
                             x_a,
                             z_c_data, z_d_data,
                             alpha1_star, alpha2_star,
                             tau1, tau2,
                             w_min, w_max,
                             tt = 20) {

  n <- length(y_data)

  p_theta <- length(betasigma)
  p_beta <- p_theta - 1L

  beta <- betasigma[seq_len(p_beta)]
  log_sig <- betasigma[p_theta]
  sigma <- exp(log_sig)

  len_beta <- p_beta + 1L
  x_len <- length(x_a)

  # Continuous covariates.
  if (is.null(z_c_data)) {

    z_c_mat <- matrix(0, nrow = n, ncol = 0)
    p_c <- 0

  } else if (is.null(dim(z_c_data))) {

    if (length(z_c_data) != n) {
      stop("z_c_data must have length n")
    }

    z_c_mat <- matrix(
      as.numeric(z_c_data),
      ncol = 1
    )

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

    z_d_mat <- matrix(
      z_d_data,
      ncol = 1
    )

    p_d <- 1

  } else {

    if (nrow(z_d_data) != n) {
      stop("z_d_data must have n rows")
    }

    z_d_mat <- as.matrix(z_d_data)
    p_d <- ncol(z_d_mat)
  }

  # Precompute a0 on a grid over the observed covariates.
  n_c_grid <- 10

  if (p_c > 0) {

    z_c_grid_list <- lapply(
      seq_len(p_c),
      function(j) {
        range_j <- range(
          z_c_mat[, j]
        )

        seq(
          range_j[1],
          range_j[2],
          length.out = n_c_grid
        )
      }
    )

  } else {

    z_c_grid_list <- list()
  }

  if (p_d > 0) {

    levels_d_list <- lapply(
      seq_len(p_d),
      function(j) {
        sort(
          unique(
            z_d_mat[, j]
          )
        )
      }
    )

    n_d_levels <- vapply(
      levels_d_list,
      length,
      integer(1)
    )

  } else {

    levels_d_list <- list()
    n_d_levels <- integer(0)
  }

  n_zc_comb <- if (p_c > 0) {
    n_c_grid^p_c
  } else {
    1L
  }

  n_zd_comb <- if (p_d > 0) {
    prod(n_d_levels)
  } else {
    1L
  }

  a0_array <- array(
    NA_real_,
    dim = c(
      x_len,
      len_beta,
      n_zc_comb,
      n_zd_comb
    )
  )

  z_d_from_index <- function(l) {

    if (p_d == 0) {
      return(numeric(0))
    }

    tmp <- l - 1L
    idx_d <- integer(p_d)

    for (j in seq_len(p_d)) {
      base_j <- n_d_levels[j]
      idx_d[j] <- (tmp %% base_j) + 1L
      tmp <- tmp %/% base_j
    }

    vapply(
      seq_len(p_d),
      function(j) {
        levels_d_list[[j]][idx_d[j]]
      },
      numeric(1)
    )
  }

  index_from_zd <- function(z_d_i) {

    if (p_d == 0) {
      return(1L)
    }

    if (length(z_d_i) != p_d) {
      stop("z_d_i has wrong length")
    }

    idx_d <- integer(p_d)

    for (j in seq_len(p_d)) {

      idx_j <- match(
        z_d_i[j],
        levels_d_list[[j]]
      )

      if (is.na(idx_j)) {
        stop("z_d_i[j] not among observed levels")
      }

      idx_d[j] <- idx_j
    }

    l <- 1L
    mult <- 1L

    for (j in seq_len(p_d)) {
      l <- l + (idx_d[j] - 1L) * mult
      mult <- mult * n_d_levels[j]
    }

    l
  }

  lin_index_zc <- function(k_vec) {

    if (p_c == 0) {
      return(1L)
    }

    idx <- 1L
    mult <- 1L

    for (j in seq_len(p_c)) {
      idx <- idx + (k_vec[j] - 1L) * mult
      mult <- mult * n_c_grid
    }

    idx
  }

  if (p_c > 0) {
    k_c_vec <- rep(1L, p_c)
  } else {
    k_c_vec <- integer(0)
  }

  for (r in seq_len(n_zc_comb)) {

    if (p_c > 0) {

      z_c_point <- vapply(
        seq_len(p_c),
        function(j) {
          z_c_grid_list[[j]][k_c_vec[j]]
        },
        numeric(1)
      )

    } else {

      z_c_point <- numeric(0)
    }

    for (l in seq_len(n_zd_comb)) {

      z_d_point <- z_d_from_index(l)

      a0_mat <- a_gauss_param12(
        beta = beta,
        x_a = x_a,
        z_c = z_c_point,
        z_d = z_d_point,
        alpha1_star = alpha1_star,
        alpha2_star = alpha2_star,
        sigma = sigma,
        tau1 = tau1,
        tau2 = tau2,
        tt = tt,
        w_min = w_min,
        w_max = w_max
      )

      if (
        !is.matrix(a0_mat) ||
        !all(
          dim(a0_mat) ==
          c(
            x_len,
            len_beta
          )
        )
      ) {
        stop("a_gauss_param12 must return an x_len x len_beta matrix.")
      }

      a0_array[, , r, l] <- a0_mat
    }

    if (p_c > 0) {

      for (j in seq_len(p_c)) {

        k_c_vec[j] <- k_c_vec[j] + 1L

        if (k_c_vec[j] <= n_c_grid) {
          break
        }

        k_c_vec[j] <- 1L
      }
    }
  }

  get_a0_for_obs <- function(z_c_i, z_d_i) {

    l <- index_from_zd(z_d_i)

    if (p_c == 0) {
      return(
        a0_array[, , 1, l]
      )
    }

    k_low <- integer(p_c)
    k_high <- integer(p_c)
    t_vec <- numeric(p_c)

    for (j in seq_len(p_c)) {

      grid_j <- z_c_grid_list[[j]]
      z_ij <- z_c_i[j]

      if (z_ij <= grid_j[1]) {

        k_low[j] <- 1L
        k_high[j] <- 1L
        t_vec[j] <- 0

      } else if (z_ij >= grid_j[n_c_grid]) {

        k_low[j] <- n_c_grid
        k_high[j] <- n_c_grid
        t_vec[j] <- 0

      } else {

        k_low[j] <- max(
          which(
            grid_j <= z_ij
          )
        )

        k_high[j] <- k_low[j] + 1L

        t_vec[j] <- (
          z_ij -
            grid_j[k_low[j]]
        ) /
          (
            grid_j[k_high[j]] -
              grid_j[k_low[j]]
          )
      }
    }

    res <- matrix(
      0,
      nrow = x_len,
      ncol = len_beta
    )

    n_vert <- 2^p_c

    for (v in 0:(n_vert - 1L)) {

      k_vec_v <- integer(p_c)
      weight <- 1

      for (j in seq_len(p_c)) {

        bit_j <- (
          v %/%
            2^(j - 1L)
        ) %% 2L

        if (bit_j == 0L) {

          k_vec_v[j] <- k_low[j]
          wj <- 1 - t_vec[j]

        } else {

          k_vec_v[j] <- k_high[j]
          wj <- t_vec[j]
        }

        weight <- weight * wj
      }

      r_idx <- lin_index_zc(
        k_vec_v
      )

      res <- res +
        weight *
        a0_array[, , r_idx, l]
    }

    res
  }

  val_full <- rep(
    0,
    len_beta
  )

  for (i in seq_len(n)) {

    z_c_i <- if (p_c > 0) {
      z_c_mat[i, ]
    } else {
      numeric(0)
    }

    z_d_i <- if (p_d > 0) {
      z_d_mat[i, ]
    } else {
      numeric(0)
    }

    a0_i <- get_a0_for_obs(
      z_c_i,
      z_d_i
    )

    S_i <- S_eff_gauss_param12(
      beta = beta,
      y = y_data[i],
      w = w_data[i],
      delta = delta_data[i],
      x_a = x_a,
      a0 = a0_i,
      z_c = z_c_i,
      z_d = z_d_i,
      alpha1_star = alpha1_star,
      sigma = sigma,
      tau1 = tau1,
      tau2 = tau2,
      w_min = w_min,
      w_max = w_max,
      tt = tt
    )

    val_full <- val_full +
      S_i
  }

  val_full
}
