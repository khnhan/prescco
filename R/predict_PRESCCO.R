# =============================================================================
# PRESCCO prediction intervals
#
# Internal numerical functions construct the PRESCCO estimating-function
# arrays and solve for the prediction half-length. The user-facing function is
# PRESCCO_prediction_interval().
# =============================================================================


# Outcome-model intercept and slope in X at one covariate value.

.outcome_intercept_slope <- function(beta, z_c = NULL, z_d = NULL) {

  intercept <- m0(
    0,
    beta,
    z_c = z_c,
    z_d = z_d
  )

  slope <- m0(
    1,
    beta,
    z_c = z_c,
    z_d = z_d
  ) - intercept

  c(
    intercept = intercept,
    slope = slope
  )
}


# =============================================================================
# Estimating-function components
# =============================================================================

c1_xz_gauss_param12 <- function(
    zeta, alpha, r,
    beta, x_a,
    z_c, z_d,
    alpha1_star, alpha2_star,
    alpha1_star_r, tau1_r,
    sigma, tau1, tau2,
    w_min, w_max,
    tt = 20, tt2 = 100) {

  cc <- gauss(tt2)

  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)

  mean_c <- sum(
    c(1, z_c_vec, z_d_vec) *
      alpha2_star
  )

  temp2 <- function(x, y) {

    if (x < w_min) {

      return(
        as.numeric(
          r(
            y = y,
            w = x,
            delta = 1,
            beta = beta,
            alpha1_star = alpha1_star_r,
            tau1 = tau1_r,
            z_c = z_c_vec,
            z_d = z_d_vec,
            w_min = w_min,
            w_max = w_max
          ) <= zeta
        )
      )
    }

    c_grid <- seq(
      w_min,
      w_max,
      length.out = tt
    )

    dens_c <- truncnorm::dtruncnorm(
      c_grid,
      a = w_min,
      b = w_max,
      mean = mean_c,
      sd = tau2
    )

    dens_c <- dens_c / sum(dens_c)

    w_grid <- pmin(
      x,
      c_grid
    )

    delta_grid <- as.numeric(
      x <= c_grid
    )

    residuals <- r(
      y = rep(y, tt),
      w = w_grid,
      delta = delta_grid,
      beta = beta,
      alpha1_star = alpha1_star_r,
      tau1 = tau1_r,
      z_c = if (length(z_c_vec)) {
        matrix(
          rep(z_c_vec, each = tt),
          nrow = tt
        )
      } else {
        NULL
      },
      z_d = if (length(z_d_vec)) {
        matrix(
          rep(z_d_vec, each = tt),
          nrow = tt
        )
      } else {
        NULL
      },
      w_min = w_min,
      w_max = w_max
    )

    sum(
      (residuals <= zeta) *
        dens_c
    )
  }

  temp3 <- function(x) {

    mu_xz <- m0(
      x,
      beta,
      z_c = z_c_vec,
      z_d = z_d_vec
    )

    y_grid <- mu_xz +
      sigma * cc$x

    vals <- vapply(
      y_grid,
      function(y_norm) {
        temp2(
          x,
          y_norm
        )
      },
      numeric(1)
    )

    sum(
      vals *
        cc$w
    )
  }

  temp3 <- Vectorize(
    temp3,
    vectorize.args = "x"
  )

  temp3(x_a) -
    (1 - alpha)
}


M1_xz_gauss_param12 <- function(
    beta, x_a,
    z_c, z_d,
    alpha1_star, alpha2_star,
    sigma, tau1, tau2,
    w_min, w_max,
    tt = 20, tt2 = 20) {

  cc <- gauss(tt2)

  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)

  outcome <- .outcome_intercept_slope(
    beta,
    z_c = z_c_vec,
    z_d = z_d_vec
  )

  intercept_x <- outcome[["intercept"]]
  slope_x <- outcome[["slope"]]

  mu_x <- sum(
    c(1, z_c_vec, z_d_vec) *
      alpha1_star
  )

  mu_c <- sum(
    c(1, z_c_vec, z_d_vec) *
      alpha2_star
  )

  temp <- function(x_grid, c, y) {

    v_x <- 1 / (
      slope_x^2 / sigma^2 +
        1 / tau1^2
    )

    eta_x <- v_x * (
      slope_x * (y - intercept_x) / sigma^2 +
        mu_x / tau1^2
    )

    p <- truncnorm::dtruncnorm(
      x_grid,
      a = w_min,
      b = w_max,
      mean = eta_x,
      sd = sqrt(v_x)
    )

    num <- (x_grid > c) * p
    denom <- sum(num)

    ifelse(
      is.nan(num / denom),
      0,
      num / denom
    )
  }

  temp2 <- function(x_grid, x, y) {

    c_grid <- seq(
      w_min,
      w_max,
      length.out = tt
    )

    dens <- truncnorm::dtruncnorm(
      c_grid,
      a = w_min,
      b = w_max,
      mean = mu_c,
      sd = tau2
    )

    dens <- dens / sum(dens)

    sapply(
      c_grid,
      function(c_norm) {
        (x > c_norm) *
          temp(
            x_grid,
            c_norm,
            y
          )
      }
    ) %*%
      dens
  }

  temp3 <- function(x_grid, x) {

    y_grid <- (
      intercept_x +
        slope_x * x
    ) +
      sigma * cc$x

    sapply(
      y_grid,
      function(y_norm) {
        temp2(
          x_grid,
          x,
          y_norm
        )
      }
    ) %*%
      cc$w
  }

  temp3 <- Vectorize(
    temp3,
    vectorize.args = "x"
  )

  temp4 <- function(x, y) {

    c_grid <- seq(
      w_min,
      w_max,
      length.out = tt
    )

    dens <- truncnorm::dtruncnorm(
      c_grid,
      a = w_min,
      b = w_max,
      mean = mu_c,
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

    y_grid <- (
      intercept_x +
        slope_x * x
    ) +
      sigma * cc$x

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


b1_gauss_param12 <- function(
    zeta, alpha, r,
    beta, x_a,
    z_c, z_d,
    alpha1_star, alpha2_star,
    alpha1_star_r, tau1_r,
    sigma, tau1, tau2,
    w_min, w_max,
    tt = 20, tt2 = 100) {

  M1_mat <- M1_xz_gauss_param12(
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
    tt = tt,
    tt2 = tt2
  )

  c1_vec <- c1_xz_gauss_param12(
    zeta = zeta,
    alpha = alpha,
    r = r,
    beta = beta,
    x_a = x_a,
    z_c = z_c,
    z_d = z_d,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    alpha1_star_r = alpha1_star_r,
    tau1_r = tau1_r,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2,
    w_min = w_min,
    w_max = w_max,
    tt = tt,
    tt2 = tt2
  )

  MASS::ginv(M1_mat) %*%
    c1_vec
}


c2_xz_gauss_param12 <- function(
    zeta, alpha, r,
    beta, x_a,
    z_c, z_d,
    alpha1_star, alpha2_star,
    alpha1_star_r, tau1_r,
    sigma, tau1, tau2,
    w_min, w_max,
    tt = 20, tt2 = 100) {

  cc <- gauss(tt2)

  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)

  outcome <- .outcome_intercept_slope(
    beta,
    z_c = z_c_vec,
    z_d = z_d_vec
  )

  intercept_x <- outcome[["intercept"]]
  slope_x <- outcome[["slope"]]

  mu_x <- sum(
    c(1, z_c_vec, z_d_vec) *
      alpha1_star
  )

  temp2 <- function(c, x) {

    y_grid <- (
      intercept_x +
        slope_x * x
    ) +
      sigma * cc$x

    w_obs <- min(x, c)
    delta_obs <- as.numeric(x <= c)

    residuals <- r(
      y = y_grid,
      w = rep(w_obs, length(y_grid)),
      delta = rep(delta_obs, length(y_grid)),
      beta = beta,
      alpha1_star = alpha1_star_r,
      tau1 = tau1_r,
      z_c = if (length(z_c_vec)) {
        matrix(
          rep(z_c_vec, each = length(y_grid)),
          nrow = length(y_grid)
        )
      } else {
        NULL
      },
      z_d = if (length(z_d_vec)) {
        matrix(
          rep(z_d_vec, each = length(y_grid)),
          nrow = length(y_grid)
        )
      } else {
        NULL
      },
      w_min = w_min,
      w_max = w_max
    )

    sum(
      (residuals <= zeta) *
        cc$w
    )
  }

  temp3 <- function(c) {

    x_grid <- seq(
      w_min,
      w_max,
      length.out = tt
    )

    dens_x <- truncnorm::dtruncnorm(
      x_grid,
      a = w_min,
      b = w_max,
      mean = mu_x,
      sd = tau1
    )

    dens_x <- dens_x / sum(dens_x)

    vals <- vapply(
      x_grid,
      function(x) {
        temp2(
          c,
          x
        )
      },
      numeric(1)
    )

    sum(
      vals *
        dens_x
    )
  }

  temp3 <- Vectorize(
    temp3,
    vectorize.args = "c"
  )

  temp3(x_a) -
    (1 - alpha)
}


M2_xz_gauss_param12 <- function(
    beta, x_a,
    z_c, z_d,
    alpha1_star, alpha2_star,
    sigma, tau1, tau2,
    w_min, w_max,
    tt = 50, tt2 = 50) {

  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)

  mu_x <- sum(
    c(1, z_c_vec, z_d_vec) *
      alpha1_star
  )

  mu_c <- sum(
    c(1, z_c_vec, z_d_vec) *
      alpha2_star
  )

  temp <- function(x_grid, x) {

    p <- truncnorm::dtruncnorm(
      x_grid,
      a = w_min,
      b = w_max,
      mean = mu_c,
      sd = tau2
    )

    num <- (x_grid >= x) * p
    denom <- sum(num)

    ifelse(
      is.nan(num / denom),
      0,
      num / denom
    )
  }

  temp2 <- function(x_grid, c) {

    x_grid_int <- seq(
      w_min,
      w_max,
      length.out = tt
    )

    dens_x <- truncnorm::dtruncnorm(
      x_grid_int,
      a = w_min,
      b = w_max,
      mean = mu_x,
      sd = tau1
    )

    dens_x <- dens_x / sum(dens_x)

    sapply(
      x_grid_int,
      function(x_norm) {
        (x_norm <= c) *
          temp(
            x_grid,
            x_norm
          )
      }
    ) %*%
      dens_x
  }

  temp2 <- Vectorize(
    temp2,
    vectorize.args = "c"
  )

  temp3 <- function(c) {

    x_grid_int <- seq(
      w_min,
      w_max,
      length.out = tt
    )

    dens_x <- truncnorm::dtruncnorm(
      x_grid_int,
      a = w_min,
      b = w_max,
      mean = mu_x,
      sd = tau1
    )

    dens_x <- dens_x / sum(dens_x)

    sum(
      (x_grid_int > c) *
        dens_x
    )
  }

  temp3 <- Vectorize(
    temp3,
    vectorize.args = "c"
  )

  diag(
    temp3(x_a)
  ) +
    t(
      temp2(
        x_a,
        x_a
      )
    )
}


b2_gauss_param12 <- function(
    zeta, alpha, r,
    beta, x_a,
    z_c, z_d,
    alpha1_star, alpha2_star,
    alpha1_star_r, tau1_r,
    sigma, tau1, tau2,
    w_min, w_max,
    tt = 50, tt2 = 50) {

  M2_mat <- M2_xz_gauss_param12(
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
    tt = tt,
    tt2 = tt2
  )

  c2_vec <- c2_xz_gauss_param12(
    zeta = zeta,
    alpha = alpha,
    r = r,
    beta = beta,
    x_a = x_a,
    z_c = z_c,
    z_d = z_d,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    alpha1_star_r = alpha1_star_r,
    tau1_r = tau1_r,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2,
    w_min = w_min,
    w_max = w_max,
    tt = tt,
    tt2 = tt2
  )

  MASS::ginv(M2_mat) %*%
    c2_vec
}


b3_gauss_param12 <- function(
    zeta, alpha, r,
    beta, x_a,
    z_c, z_d,
    alpha1_star, alpha2_star,
    alpha1_star_r, tau1_r,
    sigma, tau1, tau2,
    w_min, w_max,
    tt = 50, tt2 = 50, tt3 = 50) {

  cc <- gauss(tt2)

  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)

  outcome <- .outcome_intercept_slope(
    beta,
    z_c = z_c_vec,
    z_d = z_d_vec
  )

  intercept_x <- outcome[["intercept"]]
  slope_x <- outcome[["slope"]]

  mu_x <- sum(
    c(1, z_c_vec, z_d_vec) *
      alpha1_star
  )

  mu_c <- sum(
    c(1, z_c_vec, z_d_vec) *
      alpha2_star
  )

  temp2 <- function(c, x) {

    y_grid <- (
      intercept_x +
        slope_x * x
    ) +
      sigma * cc$x

    w_obs <- min(x, c)
    delta_obs <- as.numeric(x <= c)

    residuals <- r(
      y = y_grid,
      w = rep(w_obs, length(y_grid)),
      delta = rep(delta_obs, length(y_grid)),
      beta = beta,
      alpha1_star = alpha1_star_r,
      tau1 = tau1_r,
      z_c = if (length(z_c_vec)) {
        matrix(
          rep(z_c_vec, each = length(y_grid)),
          nrow = length(y_grid)
        )
      } else {
        NULL
      },
      z_d = if (length(z_d_vec)) {
        matrix(
          rep(z_d_vec, each = length(y_grid)),
          nrow = length(y_grid)
        )
      } else {
        NULL
      },
      w_min = w_min,
      w_max = w_max
    )

    sum(
      (residuals <= zeta) *
        cc$w
    )
  }

  temp3 <- function(c) {

    x_grid <- seq(
      w_min,
      w_max,
      length.out = tt
    )

    dens_x <- truncnorm::dtruncnorm(
      x_grid,
      a = w_min,
      b = w_max,
      mean = mu_x,
      sd = tau1
    )

    dens_x <- dens_x / sum(dens_x)

    vals <- vapply(
      x_grid,
      function(x) {
        temp2(
          c,
          x
        )
      },
      numeric(1)
    )

    sum(
      vals *
        dens_x
    )
  }

  c_grid <- seq(
    w_min,
    w_max,
    length.out = tt3
  )

  dens_c <- truncnorm::dtruncnorm(
    c_grid,
    a = w_min,
    b = w_max,
    mean = mu_c,
    sd = tau2
  )

  dens_c <- dens_c / sum(dens_c)

  vals <- vapply(
    c_grid,
    temp3,
    numeric(1)
  )

  sum(
    vals *
      dens_c
  ) -
    (1 - alpha)
}


trans_phi_eff_gauss_param12 <- function(
    beta, y, w, delta,
    x_a, b1, b2, b3,
    z_c, z_d,
    alpha1_star, alpha2_star,
    sigma, tau1, tau2,
    w_min, w_max) {

  z_c_vec <- if (is.null(z_c)) numeric(0) else as.numeric(z_c)
  z_d_vec <- if (is.null(z_d)) numeric(0) else as.numeric(z_d)

  outcome <- .outcome_intercept_slope(
    beta,
    z_c = z_c_vec,
    z_d = z_d_vec
  )

  intercept_x <- outcome[["intercept"]]
  slope_x <- outcome[["slope"]]

  mu_x <- sum(
    c(1, z_c_vec, z_d_vec) *
      alpha1_star
  )

  mu_c <- sum(
    c(1, z_c_vec, z_d_vec) *
      alpha2_star
  )

  if (delta == 1) {

    temp <- function(x_grid, x) {

      p <- truncnorm::dtruncnorm(
        x_grid,
        a = w_min,
        b = w_max,
        mean = mu_c,
        sd = tau2
      )

      num <- (x_grid >= x) * p
      denom <- sum(num)

      ifelse(
        is.nan(num / denom),
        0,
        num / denom
      )
    }

    b1_w <- approx(
      x_a,
      b1,
      w,
      rule = 2
    )$y

    return(
      b1_w +
        as.numeric(
          t(b2) %*%
            temp(
              x_a,
              w
            )
        ) -
        b3
    )
  }

  temp <- function(x_grid, c, y) {

    v_x <- 1 / (
      slope_x^2 / sigma^2 +
        1 / tau1^2
    )

    eta_x <- v_x * (
      slope_x * (y - intercept_x) / sigma^2 +
        mu_x / tau1^2
    )

    p <- truncnorm::dtruncnorm(
      x_grid,
      a = w_min,
      b = w_max,
      mean = eta_x,
      sd = sqrt(v_x)
    )

    num <- (x_grid > c) * p
    denom <- sum(num)

    ifelse(
      is.nan(num / denom),
      0,
      num / denom
    )
  }

  b2_w <- approx(
    x_a,
    b2,
    w,
    rule = 2
  )$y

  as.numeric(
    t(b1) %*%
      temp(
        x_a,
        w,
        y
      )
  ) +
    b2_w -
    b3
}


# =============================================================================
# Build and interpolate b1, b2, and b3
# =============================================================================

build_b123_arrays <- function(
    zeta_seq, alpha, r,
    beta, x_a,
    z_c_data, z_d_data,
    alpha1_star, alpha2_star,
    alpha1_star_r, tau1_r,
    sigma, tau1, tau2,
    w_min, w_max,
    n_c_grid = 10,
    b1_tt = 20, b1_tt2 = 100,
    b2_tt = 50, b2_tt2 = 50,
    b3_tt = 50, b3_tt2 = 50, b3_tt3 = 50,
    verbose = FALSE) {

  n_obs <- if (!is.null(z_c_data)) {
    if (is.null(dim(z_c_data))) length(z_c_data) else nrow(z_c_data)
  } else if (!is.null(z_d_data)) {
    if (is.null(dim(z_d_data))) length(z_d_data) else nrow(z_d_data)
  } else {
    1L
  }

  z_c_mat <- .as_covariate_matrix(
    z_c_data,
    n_obs,
    "z_c_data"
  )

  z_d_mat <- .as_covariate_matrix(
    z_d_data,
    n_obs,
    "z_d_data"
  )

  p_c <- ncol(z_c_mat)
  p_d <- ncol(z_d_mat)

  if (p_c > 0) {

    z_c_knots <- lapply(
      seq_len(p_c),
      function(j) {

        rng <- range(
          z_c_mat[, j],
          na.rm = TRUE
        )

        seq(
          rng[1],
          rng[2],
          length.out = n_c_grid
        )
      }
    )

    z_c_grid <- as.matrix(
      do.call(
        expand.grid,
        z_c_knots
      )
    )

  } else {

    z_c_grid <- matrix(
      0,
      nrow = 1,
      ncol = 0
    )
  }

  if (p_d > 0) {

    z_d_levels <- lapply(
      seq_len(p_d),
      function(j) {
        sort(
          unique(
            z_d_mat[, j]
          )
        )
      }
    )

    z_d_grid <- as.matrix(
      do.call(
        expand.grid,
        z_d_levels
      )
    )

  } else {

    z_d_grid <- matrix(
      0,
      nrow = 1,
      ncol = 0
    )
  }

  n_zeta <- length(zeta_seq)
  m <- length(x_a)
  n_zc <- nrow(z_c_grid)
  n_zd <- nrow(z_d_grid)

  b1_array <- array(
    NA_real_,
    dim = c(
      n_zeta,
      m,
      n_zc,
      n_zd
    )
  )

  b2_array <- array(
    NA_real_,
    dim = c(
      n_zeta,
      m,
      n_zc,
      n_zd
    )
  )

  b3_array <- array(
    NA_real_,
    dim = c(
      n_zeta,
      n_zc,
      n_zd
    )
  )

  for (k in seq_along(zeta_seq)) {

    zeta <- zeta_seq[k]

    if (verbose) {
      message(
        "zeta ",
        k,
        "/",
        n_zeta
      )
    }

    for (g in seq_len(n_zc)) {

      z_c_pt <- if (p_c > 0) {
        z_c_grid[g, ]
      } else {
        NULL
      }

      for (h in seq_len(n_zd)) {

        z_d_pt <- if (p_d > 0) {
          z_d_grid[h, ]
        } else {
          NULL
        }

        b1_array[k, , g, h] <- b1_gauss_param12(
          zeta = zeta,
          alpha = alpha,
          r = r,
          beta = beta,
          x_a = x_a,
          z_c = z_c_pt,
          z_d = z_d_pt,
          alpha1_star = alpha1_star,
          alpha2_star = alpha2_star,
          alpha1_star_r = alpha1_star_r,
          tau1_r = tau1_r,
          sigma = sigma,
          tau1 = tau1,
          tau2 = tau2,
          w_min = w_min,
          w_max = w_max,
          tt = b1_tt,
          tt2 = b1_tt2
        )

        b2_array[k, , g, h] <- b2_gauss_param12(
          zeta = zeta,
          alpha = alpha,
          r = r,
          beta = beta,
          x_a = x_a,
          z_c = z_c_pt,
          z_d = z_d_pt,
          alpha1_star = alpha1_star,
          alpha2_star = alpha2_star,
          alpha1_star_r = alpha1_star_r,
          tau1_r = tau1_r,
          sigma = sigma,
          tau1 = tau1,
          tau2 = tau2,
          w_min = w_min,
          w_max = w_max,
          tt = b2_tt,
          tt2 = b2_tt2
        )

        b3_array[k, g, h] <- b3_gauss_param12(
          zeta = zeta,
          alpha = alpha,
          r = r,
          beta = beta,
          x_a = x_a,
          z_c = z_c_pt,
          z_d = z_d_pt,
          alpha1_star = alpha1_star,
          alpha2_star = alpha2_star,
          alpha1_star_r = alpha1_star_r,
          tau1_r = tau1_r,
          sigma = sigma,
          tau1 = tau1,
          tau2 = tau2,
          w_min = w_min,
          w_max = w_max,
          tt = b3_tt,
          tt2 = b3_tt2,
          tt3 = b3_tt3
        )
      }
    }
  }

  list(
    b1_array = b1_array,
    b2_array = b2_array,
    b3_array = b3_array,
    z_c_grid = z_c_grid,
    z_d_grid = z_d_grid
  )
}


find_sol <- function(x_vec, y_vec) {

  sign_change <- which(
    diff(
      sign(y_vec)
    ) != 0
  )

  if (!length(sign_change)) {

    if (all(y_vec > 0)) {

      idx <- which.min(y_vec)

      x1 <- x_vec[idx]
      x2 <- x_vec[idx + 1]

      y1 <- y_vec[idx]
      y2 <- y_vec[idx + 1]

      return(
        x1 -
          y1 *
          (x2 - x1) /
          (y2 - y1)
      )
    }

    idx <- which.max(y_vec)

    x1 <- x_vec[idx - 1]
    x2 <- x_vec[idx]

    y1 <- y_vec[idx - 1]
    y2 <- y_vec[idx]

    return(
      x2 -
        y2 *
        (x2 - x1) /
        (y2 - y1)
    )
  }

  idx <- sign_change[
    which.min(
      abs(
        y_vec[sign_change]
      )
    )
  ]

  x1 <- x_vec[idx]
  x2 <- x_vec[idx + 1]

  y1 <- y_vec[idx]
  y2 <- y_vec[idx + 1]

  x1 -
    y1 *
    (x2 - x1) /
    (y2 - y1)
}


find_zeta_param_int <- function(
    y_data, w_data, delta_data,
    z_c_data, z_d_data,
    b1_array, b2_array, b3_array,
    z_c_grid, z_d_grid,
    beta_temp,
    zeta_seq,
    alpha1_star, alpha2_star,
    sigma, tau1, tau2,
    m = 20,
    w_min, w_max) {

  n <- length(y_data)

  if (
    length(w_data) != n ||
    length(delta_data) != n
  ) {
    stop(
      "y_data, w_data, and delta_data must have the same length."
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

  p_c <- ncol(z_c_mat)
  p_d <- ncol(z_d_mat)

  if (p_c > 0) {

    knots_list <- lapply(
      seq_len(p_c),
      function(j) {
        sort(
          unique(
            z_c_grid[, j]
          )
        )
      }
    )

    L_vec <- vapply(
      knots_list,
      length,
      integer(1)
    )

    cumprod_L <- if (p_c == 1) {
      1L
    } else {
      c(
        1L,
        cumprod(
          L_vec[-length(L_vec)]
        )
      )
    }

    corner_bits <- as.matrix(
      expand.grid(
        rep(
          list(
            c(0, 1)
          ),
          p_c
        )
      )
    )

  } else {

    knots_list <- list()
    cumprod_L <- 1L

    corner_bits <- matrix(
      0,
      nrow = 1,
      ncol = 0
    )
  }

  idx_zd <- integer(n)

  if (p_d == 0) {

    idx_zd[] <- 1L

  } else {

    for (i in seq_len(n)) {

      row_match <- which(
        apply(
          z_d_grid,
          1,
          function(x) {
            all(
              x ==
                z_d_mat[i, ]
            )
          }
        )
      )

      if (!length(row_match)) {
        stop(
          "Discrete covariate pattern not found for observation ",
          i,
          "."
        )
      }

      idx_zd[i] <- row_match[1]
    }
  }

  if (p_c > 0) {

    idx_low <- matrix(
      0L,
      n,
      p_c
    )

    idx_high <- matrix(
      0L,
      n,
      p_c
    )

    weight_high <- matrix(
      0,
      n,
      p_c
    )

    for (i in seq_len(n)) {

      for (j in seq_len(p_c)) {

        knots <- knots_list[[j]]
        value <- z_c_mat[i, j]

        if (value <= knots[1]) {

          idx_low[i, j] <- 1L
          idx_high[i, j] <- 1L
          weight_high[i, j] <- 0

        } else if (value >= knots[length(knots)]) {

          idx_low[i, j] <- length(knots)
          idx_high[i, j] <- length(knots)
          weight_high[i, j] <- 0

        } else {

          high <- which(
            knots >= value
          )[1]

          low <- high - 1L

          idx_low[i, j] <- low
          idx_high[i, j] <- high

          weight_high[i, j] <- (
            value -
              knots[low]
          ) /
            (
              knots[high] -
                knots[low]
            )
        }
      }
    }
  }

  x_a <- seq(
    w_min,
    w_max,
    length.out = m
  )

  func_vals <- numeric(
    length(zeta_seq)
  )

  for (zeta_idx in seq_along(zeta_seq)) {

    b1_obs <- matrix(
      0,
      nrow = n,
      ncol = m
    )

    b2_obs <- matrix(
      0,
      nrow = n,
      ncol = m
    )

    b3_obs <- numeric(n)

    for (i in seq_len(n)) {

      if (p_c == 0) {

        b1_obs[i, ] <- b1_array[
          zeta_idx,
          ,
          1,
          idx_zd[i]
        ]

        b2_obs[i, ] <- b2_array[
          zeta_idx,
          ,
          1,
          idx_zd[i]
        ]

        b3_obs[i] <- b3_array[
          zeta_idx,
          1,
          idx_zd[i]
        ]

      } else {

        weight_sum <- 0

        for (corner in seq_len(nrow(corner_bits))) {

          idx_corner <- integer(p_c)
          weight <- 1

          for (j in seq_len(p_c)) {

            if (corner_bits[corner, j] == 0) {

              idx_corner[j] <- idx_low[i, j]
              weight <- weight *
                (
                  1 -
                    weight_high[i, j]
                )

            } else {

              idx_corner[j] <- idx_high[i, j]
              weight <- weight *
                weight_high[i, j]
            }
          }

          grid_idx <- 1L +
            sum(
              (idx_corner - 1L) *
                cumprod_L
            )

          b1_obs[i, ] <- b1_obs[i, ] +
            weight *
            b1_array[
              zeta_idx,
              ,
              grid_idx,
              idx_zd[i]
            ]

          b2_obs[i, ] <- b2_obs[i, ] +
            weight *
            b2_array[
              zeta_idx,
              ,
              grid_idx,
              idx_zd[i]
            ]

          b3_obs[i] <- b3_obs[i] +
            weight *
            b3_array[
              zeta_idx,
              grid_idx,
              idx_zd[i]
            ]

          weight_sum <- weight_sum +
            weight
        }

        b1_obs[i, ] <- b1_obs[i, ] /
          weight_sum

        b2_obs[i, ] <- b2_obs[i, ] /
          weight_sum

        b3_obs[i] <- b3_obs[i] /
          weight_sum
      }
    }

    value <- 0

    for (i in seq_len(n)) {

      z_c_i <- if (p_c > 0) {
        z_c_mat[i, ]
      } else {
        NULL
      }

      z_d_i <- if (p_d > 0) {
        z_d_mat[i, ]
      } else {
        NULL
      }

      value <- value +
        trans_phi_eff_gauss_param12(
          beta = beta_temp,
          y = y_data[i],
          w = w_data[i],
          delta = delta_data[i],
          x_a = x_a,
          b1 = b1_obs[i, ],
          b2 = b2_obs[i, ],
          b3 = b3_obs[i],
          z_c = z_c_i,
          z_d = z_d_i,
          alpha1_star = alpha1_star,
          alpha2_star = alpha2_star,
          sigma = sigma,
          tau1 = tau1,
          tau2 = tau2,
          w_min = w_min,
          w_max = w_max
        )
    }

    func_vals[zeta_idx] <- value
  }

  list(
    vec = rbind(
      zeta_seq,
      func_vals
    ),
    sol = find_sol(
      zeta_seq,
      func_vals
    )
  )
}


# =============================================================================
# Residual grid and empirical coverage
# =============================================================================

#' Empirical prediction coverage
#'
#' Evaluates the empirical coverage of a fixed half-length on a test set.
#'
#' @param zeta Half-length.
#' @param r Residual function, such as \code{r1} or \code{r2}.
#' @param test_y_data,test_w_data,test_delta_data Test outcome, observed time,
#'   and event indicator.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param beta Fitted outcome-model coefficient vector.
#' @param alpha1_star_r Coefficients used by the residual.
#' @param tau1 Standard deviation used by the residual.
#' @param w_min,w_max Truncation bounds.
#'
#' @return The empirical coverage rate.
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
#' beta_fit <- find_beta_cc(
#'   y_vec, w_vec, delta_vec
#' )$beta_cc
#'
#' prediction_coverage_rate(
#'   zeta = 2,
#'   r = r1,
#'   test_y_data = y_vec,
#'   test_w_data = w_vec,
#'   test_delta_data = delta_vec,
#'   beta = beta_fit,
#'   alpha1_star_r = 0,
#'   tau1 = 1,
#'   w_min = -1,
#'   w_max = 1
#' )
#'
#' @export
prediction_coverage_rate <- function(
    zeta, r,
    test_y_data, test_w_data, test_delta_data,
    test_z_c_data = NULL, test_z_d_data = NULL,
    beta,
    alpha1_star_r,
    tau1,
    w_min, w_max) {

  n <- length(test_y_data)

  if (
    length(test_w_data) != n ||
    length(test_delta_data) != n
  ) {
    stop(
      "test_y_data, test_w_data, and test_delta_data must have the same length."
    )
  }

  residuals <- r(
    y = test_y_data,
    w = test_w_data,
    delta = test_delta_data,
    beta = beta,
    alpha1_star = alpha1_star_r,
    tau1 = tau1,
    z_c = test_z_c_data,
    z_d = test_z_d_data,
    w_min = w_min,
    w_max = w_max
  )

  mean(
    residuals <= zeta
  )
}


get_zeta_seq_r <- function(
    alpha, r,
    y_data, w_data, delta_data,
    z_c_data = NULL, z_d_data = NULL,
    beta,
    alpha1_star_r, tau1_r,
    w_min = -1, w_max = 1,
    margin = 0.1,
    length = 5) {

  residuals <- r(
    y = y_data,
    w = w_data,
    delta = delta_data,
    beta = beta,
    alpha1_star = alpha1_star_r,
    tau1 = tau1_r,
    z_c = z_c_data,
    z_d = z_d_data,
    w_min = w_min,
    w_max = w_max
  )

  q <- stats::quantile(
    residuals,
    1 - alpha,
    names = FALSE
  )

  seq(
    q * (1 - margin),
    q * (1 + margin),
    length.out = length
  )
}


# =============================================================================
# PRESCCO
# =============================================================================

#' PRESCCO prediction interval
#'
#' Estimates a prediction half-length under a right-censored covariate. Nuisance
#' models and the outcome model are fitted when their parameters are not
#' supplied.
#'
#' @param y_data,w_data,delta_data Outcome, observed time, and event indicator.
#' @param z_c_data,z_d_data Optional continuous and discrete covariates.
#' @param alpha Miscoverage level.
#' @param alpha1_star,tau1 Optional coefficients and standard deviation for the
#'   working model for \code{X | Z}. Estimated when not supplied.
#' @param alpha1_star_r,tau1_r Optional coefficients and SD used for the
#'   residual. Use alternative values for \code{r1*}. By default,
#'   \code{alpha1_star_r = alpha1_star} and \code{tau1_r = tau1}.
#' @param alpha2_star,tau2 Optional coefficients and standard deviation for the
#'   working model for \code{C | Z}. Estimated when not supplied.
#' @param beta,sigma Optional outcome-model coefficients and standard deviation.
#'   Estimated by SPARCC when not supplied.
#' @param test_y_data,test_w_data,test_delta_data Optional test data for
#'   evaluating empirical coverage.
#' @param test_z_c_data,test_z_d_data Optional test covariates.
#' @param residual Residual function. Use \code{r1} or \code{r2}.
#' @param residual_name Label for the residual in the returned object. Set to
#'   \code{"r1star"} when using \code{r1*}.
#' @param x_a Grid over the censored covariate. By default, an equally spaced
#'   grid between \code{w_min} and \code{w_max}.
#' @param w_min,w_max Truncation bounds.
#' @param tt,m Integration grid size and covariate-grid resolution used when
#'   fitting the outcome model.
#' @param seq_length Number of candidate half-lengths.
#' @param b_args Optional named list overriding numerical settings in
#'   \code{build_b123_arrays}.
#'
#' @return A list containing the estimated half-length, optional empirical
#'   coverage rate, fitted model parameters, and numerical quantities used in
#'   the half-length solve.
#'
#' @examples
#' set.seed(1)
#' n <- 30
#' x_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' c_vec <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w_vec <- pmin(x_vec, c_vec)
#' delta_vec <- as.integer(x_vec <= c_vec)
#' y_vec <- 3 * x_vec + rnorm(n)
#'
#' fit <- PRESCCO_prediction_interval(
#'   y_data = y_vec,
#'   w_data = w_vec,
#'   delta_data = delta_vec,
#'   alpha = 0.1,
#'   alpha1_star = 0,
#'   tau1 = 1,
#'   alpha1_star_r = 0,
#'   tau1_r = 1,
#'   alpha2_star = 1,
#'   tau2 = 1,
#'   beta = c(0, 3),
#'   sigma = 1,
#'   residual = r1,
#'   w_min = -1,
#'   w_max = 1,
#'   m = 5,
#'   seq_length = 3,
#'   b_args = list(
#'     n_c_grid = 3,
#'     b1_tt = 3,
#'     b1_tt2 = 5,
#'     b2_tt = 3,
#'     b2_tt2 = 5,
#'     b3_tt = 3,
#'     b3_tt2 = 5,
#'     b3_tt3 = 3
#'   )
#' )
#' fit$zeta
#'
#' @export
PRESCCO_prediction_interval <- function(
    y_data, w_data, delta_data,
    z_c_data = NULL, z_d_data = NULL,
    alpha = 0.1,
    alpha1_star = NULL, tau1 = NULL,
    alpha1_star_r = NULL, tau1_r = NULL,
    alpha2_star = NULL, tau2 = NULL,
    beta = NULL, sigma = NULL,
    test_y_data = NULL,
    test_w_data = NULL,
    test_delta_data = NULL,
    test_z_c_data = NULL,
    test_z_d_data = NULL,
    residual = r1,
    residual_name = NULL,
    x_a = NULL,
    w_min = -1, w_max = 1,
    tt = 20, m = 20,
    seq_length = 5,
    b_args = list()) {

  n <- length(y_data)

  if (
    length(w_data) != n ||
    length(delta_data) != n
  ) {
    stop(
      "y_data, w_data, and delta_data must have the same length."
    )
  }

  if (is.null(residual_name)) {
    residual_name <- deparse(
      substitute(residual)
    )
  }

  if (is.null(x_a)) {
    x_a <- seq(
      w_min,
      w_max,
      length.out = m
    )
  }

  if (
    is.null(alpha1_star) ||
    is.null(tau1)
  ) {

    fit1 <- find_alpha1_MLE(
      w_data = w_data,
      delta_data = delta_data,
      z_c_data = z_c_data,
      z_d_data = z_d_data,
      w_min = w_min,
      w_max = w_max
    )

    alpha1_star <- fit1[
      -length(fit1)
    ]

    tau1 <- fit1[
      length(fit1)
    ]
  }

  if (
    is.null(alpha2_star) ||
    is.null(tau2)
  ) {

    fit2 <- find_alpha2_MLE(
      w_data = w_data,
      delta_data = delta_data,
      z_c_data = z_c_data,
      z_d_data = z_d_data,
      w_min = w_min,
      w_max = w_max
    )

    alpha2_star <- fit2[
      -length(fit2)
    ]

    tau2 <- fit2[
      length(fit2)
    ]
  }

  if (
    is.null(beta) ||
    is.null(sigma)
  ) {

    fit_beta <- find_beta_sparcc(
      y_data = y_data,
      w_data = w_data,
      delta_data = delta_data,
      z_c_data = z_c_data,
      z_d_data = z_d_data,
      alpha1_star = alpha1_star,
      alpha2_star = alpha2_star,
      tau1 = tau1,
      tau2 = tau2,
      tt = tt,
      m = m,
      w_min = w_min,
      w_max = w_max
    )

    beta <- fit_beta$beta_hat
    sigma <- fit_beta$sigma_hat
  }

  if (is.null(alpha1_star_r)) {
    alpha1_star_r <- alpha1_star
  }

  if (is.null(tau1_r)) {
    tau1_r <- tau1
  }

  zeta_seq <- get_zeta_seq_r(
    alpha = alpha,
    r = residual,
    y_data = y_data,
    w_data = w_data,
    delta_data = delta_data,
    z_c_data = z_c_data,
    z_d_data = z_d_data,
    beta = beta,
    alpha1_star_r = alpha1_star_r,
    tau1_r = tau1_r,
    w_min = w_min,
    w_max = w_max,
    length = seq_length
  )

  b <- do.call(
    build_b123_arrays,
    c(
      list(
        zeta_seq = zeta_seq,
        alpha = alpha,
        r = residual,
        beta = beta,
        x_a = x_a,
        z_c_data = z_c_data,
        z_d_data = z_d_data,
        alpha1_star = alpha1_star,
        alpha2_star = alpha2_star,
        alpha1_star_r = alpha1_star_r,
        tau1_r = tau1_r,
        sigma = sigma,
        tau1 = tau1,
        tau2 = tau2,
        w_min = w_min,
        w_max = w_max
      ),
      b_args
    )
  )

  zeta_fit <- find_zeta_param_int(
    y_data = y_data,
    w_data = w_data,
    delta_data = delta_data,
    z_c_data = z_c_data,
    z_d_data = z_d_data,
    b1_array = b$b1_array,
    b2_array = b$b2_array,
    b3_array = b$b3_array,
    z_c_grid = b$z_c_grid,
    z_d_grid = b$z_d_grid,
    beta_temp = beta,
    zeta_seq = zeta_seq,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2,
    m = m,
    w_min = w_min,
    w_max = w_max
  )

  zeta <- stats::setNames(
    zeta_fit$sol,
    residual_name
  )

  coverage_rate <- if (
    is.null(test_y_data) ||
    is.null(test_w_data) ||
    is.null(test_delta_data)
  ) {

    stats::setNames(
      NA_real_,
      residual_name
    )

  } else {

    stats::setNames(
      prediction_coverage_rate(
        zeta = zeta_fit$sol,
        r = residual,
        test_y_data = test_y_data,
        test_w_data = test_w_data,
        test_delta_data = test_delta_data,
        test_z_c_data = test_z_c_data,
        test_z_d_data = test_z_d_data,
        beta = beta,
        alpha1_star_r = alpha1_star_r,
        tau1 = tau1_r,
        w_min = w_min,
        w_max = w_max
      ),
      residual_name
    )
  }

  list(
    method = "PRESCCO",
    alpha = alpha,
    residual = residual_name,
    zeta = zeta,
    coverage_rate = coverage_rate,
    beta = beta,
    sigma = sigma,
    alpha1_star = alpha1_star,
    tau1 = tau1,
    alpha2_star = alpha2_star,
    tau2 = tau2,
    alpha1_star_r = alpha1_star_r,
    tau1_r = tau1_r,
    zeta_seq = zeta_seq,
    zeta_list = zeta_fit
  )
}
