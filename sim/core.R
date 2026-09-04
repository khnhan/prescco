# =============================================================================
# Core functions for the PRESCCO simulation study
#
# This file contains:
#
#   1. Data generation
#   2. SPARCC estimating equations
#   3. Prediction centers and residuals
#   4. Numerical integral-equation solvers for PRESCCO
#   5. Half-length estimation
#   6. Shared simulation settings
#
# This file is sourced by the other simulation scripts and does not run a
# simulation by itself.
# =============================================================================


# =============================================================================
# Data generation
# =============================================================================

data_generating <- function(k, n, d, beta, alpha1, alpha2,
                            sigma = 4, tau1 = 1, tau2 = 1) {

  set.seed(k)

  # Event-time covariate X.
  x_data <- truncnorm::rtruncnorm(
    n,
    a = -1,
    b = 1,
    mean = alpha1,
    sd = tau1
  )

  # Outcome Y.
  trunc_norm <- truncnorm::rtruncnorm(
    n,
    a = -3,
    b = 3,
    mean = 0,
    sd = 1
  )

  y_data <- cbind(1, x_data) %*% beta + sigma * trunc_norm

  # Censoring time C.
  c_data <- truncnorm::rtruncnorm(
    n,
    a = -1,
    b = 1,
    mean = alpha2,
    sd = tau2
  )

  # Observed covariate and censoring indicator.
  w_data <- pmin(c_data, x_data)
  delta_data <- as.numeric(x_data <= c_data)

  cens_rate <- sum(c_data <= x_data) / n

  list(
    x_data = x_data,
    y_data = y_data,
    c_data = c_data,
    w_data = w_data,
    delta_data = delta_data,
    cens_rate = cens_rate
  )
}


# =============================================================================
# SPARCC beta estimation
# =============================================================================

# Score for beta when X is observed.

sim_S_beta_f <- function(beta, y, x, sigma = 4) {

  res <- y - (beta[1] + beta[2] * x)

  res * c(1, x) / sigma^2
}

sim_S_beta_f <- Vectorize(
  sim_S_beta_f,
  vectorize.args = "x"
)


# Observed-data score for beta.

sim_S_beta <- function(beta, y, w, delta, alpha1_star,
                       sigma = 4, tau1 = 1) {

  v_x <- 1 / (
    beta[2]^2 / sigma^2 +
      1 / tau1^2
  )

  eta_x <- v_x * (
    beta[2] * (y - beta[1]) / sigma^2 +
      alpha1_star / tau1^2
  )

  if (delta == 0) {

    if (w > 1) {

      sim_S_beta_f(
        beta,
        y,
        w,
        sigma
      )

    } else {

      x_norm <- (
        seq(
          w,
          1,
          length.out = 20
        ) -
          eta_x
      ) / sqrt(v_x)

      d <- dnorm(x_norm)

      num <- sim_S_beta_f(
        beta,
        y,
        x_norm * sqrt(v_x) + eta_x,
        sigma
      ) %*% d

      denom <- sum(d)

      num / denom
    }

  } else {

    sim_S_beta_f(
      beta,
      y,
      w,
      sigma
    )
  }
}


# Gaussian quadrature grid.

sim_gauss <- function(tt, len = 3) {

  grid <- seq(
    -len,
    len,
    length.out = tt
  )

  d <- dnorm(grid)

  list(
    x = grid,
    w = d / sum(d)
  )
}


# Right-hand side of the integral equation for beta.

sim_c0_xz_gauss_param12 <- function(
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20) {

  cc <- sim_gauss(tt)
  len <- 20

  temp2 <- function(x, y) {

    v_x <- 1 / (
      beta[2]^2 / sigma^2 +
        1 / tau1^2
    )

    eta_x <- v_x * (
      beta[2] * (y - beta[1]) / sigma^2 +
        alpha1_star / tau1^2
    )

    if (x < -1) {

      sim_S_beta_f(
        beta,
        y,
        x,
        sigma
      )

    } else {

      c_grid <- seq(
        -1,
        min(x, 1),
        length.out = len
      )

      dens_c <- truncnorm::dtruncnorm(
        c_grid,
        a = -1,
        b = 1,
        mean = alpha2_star,
        sd = tau2
      )

      Sbeta <- sapply(
        c_grid,
        function(w) {
          sim_S_beta(
            beta,
            y,
            w,
            0,
            alpha1_star,
            sigma,
            tau1
          )
        }
      )

      Sbeta[, c(1, len)] <- Sbeta[, c(1, len)] / 2

      by <- c_grid[2] - c_grid[1]

      sim_S_beta_f(
        beta,
        y,
        x,
        sigma
      ) *
        (
          1 -
            truncnorm::ptruncnorm(
              x,
              a = -1,
              b = 1,
              mean = alpha2_star,
              sd = tau2
            )
        ) +
        Sbeta %*% dens_c * by
    }
  }


  temp3 <- function(x) {

    sapply(
      sigma * cc$x + sum(c(1, x) * beta),
      function(y_norm) {
        temp2(
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

  temp3(x_a)
}


# Linear operator in the integral equation for beta.

sim_L_xz_gauss_param12 <- function(
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20) {

  cc <- sim_gauss(tt)


  temp <- function(x_a, c, y) {

    v_x <- 1 / (
      beta[2]^2 / sigma^2 +
        1 / tau1^2
    )

    eta_x <- v_x * (
      beta[2] * (y - beta[1]) / sigma^2 +
        alpha1_star / tau1^2
    )

    p <- truncnorm::dtruncnorm(
      x_a,
      a = -1,
      b = 1,
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
      -1,
      1,
      length = 20
    )

    dens <- truncnorm::dtruncnorm(
      c_grid,
      a = -1,
      b = 1,
      mean = alpha2_star,
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

    sapply(
      sigma * cc$x + sum(beta * c(1, x)),
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
      -1,
      1,
      length = 20
    )

    dens <- truncnorm::dtruncnorm(
      c_grid,
      a = -1,
      b = 1,
      mean = alpha2_star,
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

    sum(
      temp4(
        x,
        sigma * cc$x +
          sum(beta * c(1, x))
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


# Solution of the beta integral equation.

sim_a_gauss_param12 <- function(
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20) {

  MASS::ginv(
    sim_L_xz_gauss_param12(
      beta,
      x_a,
      alpha1_star,
      alpha2_star,
      sigma,
      tau1,
      tau2,
      tt
    )
  ) %*%
    t(
      sim_c0_xz_gauss_param12(
        beta,
        x_a,
        alpha1_star,
        alpha2_star,
        sigma,
        tau1,
        tau2,
        tt
      )
    )
}


# Efficient score for beta.

sim_S_eff_gauss_param12 <- function(
    beta,
    y,
    w,
    delta,
    x_a,
    a0,
    alpha1_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20) {

  len_beta <- length(beta)

  if (delta == 0) {

    temp <- function(x_a, c, y) {

      v_x <- 1 / (
        beta[2]^2 / sigma^2 +
          1 / tau1^2
      )

      eta_x <- v_x * (
        beta[2] * (y - beta[1]) / sigma^2 +
          alpha1_star / tau1^2
      )

      p <- truncnorm::dtruncnorm(
        x_a,
        a = -1,
        b = 1,
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


    sbeta <- sim_S_beta(
      beta,
      y,
      w,
      0,
      alpha1_star,
      sigma,
      tau1
    )

    sbeta -
      t(a0) %*%
      temp(
        x_a,
        w,
        y
      )

  } else {

    a0_w <- vector(
      length = len_beta
    )

    for (j in seq_len(len_beta)) {

      a0_w[j] <- approx(
        x_a,
        a0[, j],
        w,
        rule = 2
      )$y
    }

    sbeta <- sim_S_beta_f(
      beta,
      y,
      w,
      sigma
    )

    sbeta - a0_w
  }
}


# Estimating equation for beta.

sim_pe_gauss_param12 <- function(
    beta,
    y_data,
    w_data,
    delta_data,
    x_a,
    alpha1_star,
    alpha2_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20) {

  len_beta <- length(beta)

  a0 <- sim_a_gauss_param12(
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    sigma,
    tau1,
    tau2,
    tt
  )

  val <- rep(
    0,
    len_beta
  )

  n <- length(y_data)

  for (i in seq_len(n)) {

    val <- val +
      sim_S_eff_gauss_param12(
        beta,
        y_data[i],
        w_data[i],
        delta_data[i],
        x_a,
        a0,
        alpha1_star,
        sigma,
        tau1,
        tau2,
        tt
      )
  }

  val
}


# Estimate alpha1 for a fixed beta.

sim_find_alpha1_MLE <- function(
    beta,
    y_data,
    w_data,
    delta_data,
    sigma = 4,
    tau1 = 1) {

  log_llhd <- function(
    alpha1_star,
    y,
    w,
    delta) {

    v_x <- 1 / (
      beta[2]^2 / sigma^2 +
        1 / tau1^2
    )

    eta_x <- v_x * (
      beta[2] * (y - beta[1]) / sigma^2 +
        alpha1_star / tau1^2
    )

    if (delta == 1) {

      log(
        truncnorm::dtruncnorm(
          w,
          a = -1,
          b = 1,
          mean = eta_x,
          sd = sqrt(v_x)
        )
      )

    } else {

      log(
        pnorm(
          1,
          mean = eta_x,
          sd = sqrt(v_x)
        ) -
          pnorm(
            w,
            mean = eta_x,
            sd = sqrt(v_x)
          )
      ) +
        dnorm(
          y,
          mean = sum(
            c(
              1,
              alpha1_star
            ) *
              beta
          ),
          sd = sqrt(
            sigma^2 +
              beta[2]^2 * tau1^2
          ),
          log = TRUE
        ) -
        log(
          pnorm(
            1,
            mean = alpha1_star,
            sd = tau1
          ) -
            pnorm(
              -1,
              mean = alpha1_star,
              sd = tau1
            )
        )
    }
  }


  log_llhd_sum <- function(alpha1) {

    n <- length(y_data)
    llhd <- 0

    for (i in seq_len(n)) {

      llhd <- llhd +
        log_llhd(
          alpha1,
          y_data[i],
          w_data[i],
          delta_data[i]
        )
    }

    -llhd
  }


  optimize(
    log_llhd_sum,
    interval = c(-5, 5)
  )$minimum
}


# Estimate alpha2.

sim_find_alpha2_MLE <- function(
    y_data,
    w_data,
    delta_data,
    sigma = 4,
    tau2 = 1) {

  log_llhd <- function(
    alpha2_star,
    y,
    w,
    delta) {

    if (delta == 1) {

      log(
        1 -
          truncnorm::ptruncnorm(
            w,
            a = -1,
            b = 1,
            mean = sum(
              c(1, y) *
                alpha2_star
            ),
            sd = tau2
          )
      )

    } else {

      log(
        truncnorm::dtruncnorm(
          w,
          a = -1,
          b = 1,
          mean = sum(
            c(1, y) *
              alpha2_star
          ),
          sd = tau2
        )
      )
    }
  }


  log_llhd_sum <- function(alpha2) {

    n <- length(y_data)
    llhd <- 0

    for (i in seq_len(n)) {

      llhd <- llhd +
        log_llhd(
          alpha2,
          y_data[i],
          w_data[i],
          delta_data[i]
        )
    }

    -llhd
  }


  optim(
    rep(0, d + 1),
    log_llhd_sum
  )$par
}


# Estimating equation with alpha1 estimated by maximum likelihood.

pe_gauss_alpha1_MLE <- function(
    beta,
    y_data,
    w_data,
    delta_data,
    x_a,
    alpha2_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20) {

  len_beta <- length(beta)

  alpha1_star <- sim_find_alpha1_MLE(
    beta,
    y_data,
    w_data,
    delta_data,
    sigma,
    tau1
  )

  a0 <- sim_a_gauss_param12(
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    sigma,
    tau1,
    tau2,
    tt
  )

  val <- rep(
    0,
    len_beta
  )

  n <- length(y_data)

  for (i in seq_len(n)) {

    val <- val +
      sim_S_eff_gauss_param12(
        beta,
        y_data[i],
        w_data[i],
        delta_data[i],
        x_a,
        a0,
        alpha1_star,
        sigma,
        tau1,
        tau2,
        tt
      )
  }

  val
}


# SPARCC estimate of beta for one simulation replicate.

get_beta_param_param12 <- function(
    k,
    n,
    d,
    beta,
    alpha1,
    alpha2,
    alpha1_star = alpha1,
    alpha2_star = alpha2,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20,
    m = 20) {

  data_k <- data_generating(
    k,
    n,
    d,
    beta,
    alpha1,
    alpha2,
    sigma,
    tau1,
    tau2
  )

  y_data <- data_k$y_data
  w_data <- data_k$w_data
  delta_data <- data_k$delta_data

  x_a <- seq(
    -1,
    1,
    length.out = m
  )

  set.seed(k)

  result <- nleqslv::nleqslv(
    beta + rnorm(d + 1) * 0.1,
    sim_pe_gauss_param12,
    y_data = y_data,
    w_data = w_data,
    delta_data = delta_data,
    x_a = x_a,
    alpha1_star = alpha1_star,
    alpha2_star = alpha2_star,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2,
    tt = tt
  )

  print(result$message)

  result$x
}


# =============================================================================
# Prediction centers and residuals
# =============================================================================

# Center m0.

sim_m0 <- function(x, beta) {

  beta[1] +
    beta[2] * x
}


# Center m1.

sim_m1 <- function(
    w,
    delta,
    beta,
    alpha1_star,
    tau1 = 1) {

  v_x <- tau1^2
  eta_x <- alpha1_star

  if (delta == 0) {

    if (w >= 1) {

      sim_m0(
        w,
        beta
      )

    } else {

      x_exp <- alpha1_star +
        tau1 *
        (
          dnorm(
            (w - alpha1_star) /
              tau1
          ) -
            dnorm(
              (1 - alpha1_star) /
                tau1
            )
        ) /
        (
          pnorm(
            1,
            mean = eta_x,
            sd = sqrt(v_x)
          ) -
            pnorm(
              w,
              mean = eta_x,
              sd = sqrt(v_x)
            )
        )

      beta[1] +
        beta[2] * x_exp
    }

  } else {

    sim_m0(
      w,
      beta
    )
  }
}


# Residual r1.

sim_r1 <- function(
    y,
    w,
    delta,
    beta,
    alpha1_star,
    tau1 = 1) {

  abs(
    y -
      sim_m1(
        w,
        delta,
        beta,
        alpha1_star,
        tau1
      )
  )
}

sim_r1 <- Vectorize(
  sim_r1,
  vectorize.args = "w"
)


# Residual r2.

sim_r2 <- function(
    y,
    w,
    delta,
    beta,
    alpha1_star,
    tau1 = 1) {

  abs(
    y -
      sim_m0(
        w,
        beta
      )
  )
}

sim_r2 <- Vectorize(
  sim_r2,
  vectorize.args = "w"
)


# Vectorized center m0.

sim_m0_vec <- function(w, beta) {

  beta[1] +
    beta[2] * w
}


# Vectorized center m1.

sim_m1_vec <- function(
    w,
    delta,
    beta,
    alpha1_star,
    tau1 = 1) {

  w <- as.numeric(w)
  delta <- as.integer(delta)

  n <- length(w)
  out <- numeric(n)

  observed <- delta == 1L
  censored_boundary <- delta == 0L & w >= 1
  censored_inside <- delta == 0L & w < 1

  if (any(observed)) {

    out[observed] <- sim_m0_vec(
      w[observed],
      beta
    )
  }

  if (any(censored_boundary)) {

    out[censored_boundary] <- sim_m0_vec(
      w[censored_boundary],
      beta
    )
  }

  if (any(censored_inside)) {

    a <- alpha1_star
    s <- tau1

    z1 <- (1 - a) / s
    zw <- (
      w[censored_inside] -
        a
    ) / s

    denom <- pnorm(z1) -
      pnorm(zw)

    numer <- dnorm(zw) -
      dnorm(z1)

    tiny <- 1e-12
    good <- denom > tiny

    x_exp <- numeric(
      sum(censored_inside)
    )

    x_exp[good] <- a +
      s *
      (
        numer[good] /
          denom[good]
      )

    x_exp[!good] <- w[censored_inside][!good]

    out[censored_inside] <- beta[1] +
      beta[2] * x_exp
  }

  out
}


# Vectorized residual r1.

sim_r1_vec <- function(
    y,
    w,
    delta,
    beta,
    alpha1_star,
    tau1 = 1) {

  abs(
    y -
      sim_m1_vec(
        w,
        delta,
        beta,
        alpha1_star,
        tau1
      )
  )
}


# Vectorized residual r2.

sim_r2_vec <- function(
    y,
    w,
    delta,
    beta,
    alpha1_star,
    tau1 = 1) {

  abs(
    y -
      sim_m0_vec(
        w,
        beta
      )
  )
}


# =============================================================================
# Monte Carlo true half-lengths
# =============================================================================

get_zeta_r1 <- function(
    alpha,
    alpha1,
    alpha2,
    beta = c(0, 3),
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    d = 1,
    N = 1000000) {

  data_mc <- data_generating(
    9999999,
    N,
    d,
    beta,
    alpha1,
    alpha2,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2
  )

  y_data <- data_mc$y_data
  w_data <- data_mc$w_data
  delta_data <- data_mc$delta_data

  quantile(
    sapply(
      seq_len(N),
      function(i) {
        sim_r1(
          y_data[i],
          w_data[i],
          delta_data[i],
          beta,
          alpha1,
          tau1
        )
      }
    ),
    1 - alpha
  )
}


get_zeta_r2 <- function(
    alpha,
    alpha1,
    alpha2,
    beta = c(0, 3),
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    d = 1,
    N = 1000000) {

  data_mc <- data_generating(
    9999999,
    N,
    d,
    beta,
    alpha1,
    alpha2,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2
  )

  y_data <- data_mc$y_data
  w_data <- data_mc$w_data
  delta_data <- data_mc$delta_data

  quantile(
    sapply(
      seq_len(N),
      function(i) {
        sim_r2(
          y_data[i],
          w_data[i],
          delta_data[i],
          beta,
          alpha1,
          tau1
        )
      }
    ),
    1 - alpha
  )
}


# =============================================================================
# Integral-equation quantities for PRESCCO
# =============================================================================

sim_c1_xz_gauss_param12 <- function(
    zeta,
    alpha,
    r,
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    alpha1_star_r,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20,
    tt2 = 100) {

  cc <- sim_gauss(tt2)


  temp2 <- function(x, y) {

    if (x < -1) {

      r(
        y,
        x,
        1,
        beta,
        alpha1_star_r,
        tau1
      ) <= zeta

    } else {

      c_grid <- seq(
        -1,
        1,
        length.out = tt
      )

      dens_c <- truncnorm::dtruncnorm(
        c_grid,
        a = -1,
        b = 1,
        mean = alpha2_star,
        sd = tau2
      )

      w_grid <- pmin(
        x,
        c_grid
      )

      delta_grid <- as.numeric(
        x <= c_grid
      )

      RHS <- sapply(
        seq_len(tt),
        function(i) {
          r(
            y,
            w_grid[i],
            delta_grid[i],
            beta,
            alpha1_star_r,
            tau1
          ) <= zeta
        }
      )

      RHS %*%
        dens_c /
        sum(dens_c)
    }
  }


  temp3 <- function(x) {

    sapply(
      sigma * cc$x +
        sum(c(1, x) * beta),
      function(y_norm) {
        temp2(
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

  temp3(x_a) -
    (1 - alpha)
}


make_seq <- function(n) {

  seq_len(n) *
    2 /
    (2 * n + 1)
}


c1_xz_gauss_param12_2 <- function(
    zeta,
    alpha,
    r,
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    alpha1_star_r,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20,
    tt2 = 20) {

  c_vec <- truncnorm::qtruncnorm(
    make_seq(tt),
    a = -1,
    b = 1,
    mean = alpha2_star,
    sd = tau2
  )

  sapply(
    x_a,
    function(x) {

      y_vec <- sigma *
        truncnorm::qtruncnorm(
          make_seq(tt2),
          a = -3,
          b = 3,
          mean = 0,
          sd = 1
        ) +
        sum(
          c(1, x) *
            beta
        )

      mean(
        sapply(
          seq_len(tt),
          function(i) {

            mean(
              sapply(
                seq_len(tt2),
                function(j) {

                  r(
                    y_vec[j],
                    min(
                      x,
                      c_vec[i]
                    ),
                    x <= c_vec[i],
                    beta,
                    alpha1_star,
                    tau1
                  ) <= zeta
                }
              )
            )
          }
        )
      ) -
        (1 - alpha)
    }
  )
}


sim_M1_xz_gauss_param12 <- function(
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20,
    tt2 = 20) {

  cc <- sim_gauss(tt2)


  temp <- function(x_a, c, y) {

    v_x <- 1 / (
      beta[2]^2 / sigma^2 +
        1 / tau1^2
    )

    eta_x <- v_x * (
      beta[2] * (y - beta[1]) / sigma^2 +
        alpha1_star / tau1^2
    )

    p <- truncnorm::dtruncnorm(
      x_a,
      a = -1,
      b = 1,
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
      -1,
      1,
      length = tt
    )

    dens <- truncnorm::dtruncnorm(
      c_grid,
      a = -1,
      b = 1,
      mean = alpha2_star,
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
    ) %*%
      dens
  }


  temp3 <- function(x_a, x) {

    sapply(
      sigma * cc$x +
        sum(
          beta *
            c(1, x)
        ),
      function(y_norm) {

        temp2(
          x_a,
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
      -1,
      1,
      length = tt
    )

    dens <- truncnorm::dtruncnorm(
      c_grid,
      a = -1,
      b = 1,
      mean = alpha2_star,
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

    sum(
      temp4(
        x,
        sigma * cc$x +
          sum(
            beta *
              c(1, x)
          )
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


sim_b1_gauss_param12 <- function(
    zeta,
    alpha,
    r,
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    alpha1_star_r,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20,
    tt2 = 100) {

  MASS::ginv(
    sim_M1_xz_gauss_param12(
      beta,
      x_a,
      alpha1_star,
      alpha2_star,
      sigma,
      tau1,
      tau2,
      tt,
      tt2
    )
  ) %*%
    sim_c1_xz_gauss_param12(
      zeta,
      alpha,
      r,
      beta,
      x_a,
      alpha1_star,
      alpha2_star,
      alpha1_star_r,
      sigma,
      tau1,
      tau2,
      tt,
      tt2
    )
}


sim_c2_xz_gauss_param12 <- function(
    zeta,
    alpha,
    r,
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    alpha1_star_r,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20,
    tt2 = 100) {

  cc <- sim_gauss(tt2)


  temp2 <- function(c, x) {

    sapply(
      sigma * cc$x +
        sum(
          c(1, x) *
            beta
        ),
      function(y) {

        r(
          y,
          min(
            x,
            c
          ),
          x <= c,
          beta,
          alpha1_star_r,
          tau1
        ) <= zeta
      }
    ) %*%
      cc$w
  }


  temp3 <- function(c) {

    x_grid <- seq(
      -1,
      1,
      length.out = tt
    )

    dens_x <- truncnorm::dtruncnorm(
      x_grid,
      a = -1,
      b = 1,
      mean = alpha1_star,
      sd = tau1
    )

    RHS <- sapply(
      seq_len(tt),
      function(i) {
        temp2(
          c,
          x_grid[i]
        )
      }
    )

    RHS %*%
      dens_x /
      sum(dens_x)
  }

  temp3 <- Vectorize(
    temp3,
    vectorize.args = "c"
  )

  temp3(x_a) -
    (1 - alpha)
}


sim_M2_xz_gauss_param12 <- function(
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 50,
    tt2 = 50) {

  cc <- sim_gauss(tt2)


  temp <- function(x_a, x) {

    p <- truncnorm::dtruncnorm(
      x_a,
      a = -1,
      b = 1,
      mean = alpha2_star,
      sd = tau2
    )

    num <- (x_a >= x) * p
    denom <- sum(num)

    ifelse(
      is.nan(num / denom),
      0,
      num / denom
    )
  }


  temp2 <- function(x_a, c) {

    x_grid <- seq(
      -1,
      1,
      length = tt
    )

    dens <- truncnorm::dtruncnorm(
      x_grid,
      a = -1,
      b = 1,
      mean = alpha1_star,
      sd = tau1
    )

    dens <- dens / sum(dens)

    sapply(
      x_grid,
      function(x_norm) {

        (x_norm <= c) *
          temp(
            x_a,
            x_norm
          )
      }
    ) %*%
      dens
  }

  temp2 <- Vectorize(
    temp2,
    vectorize.args = "c"
  )


  temp3 <- function(c) {

    x_grid <- seq(
      -1,
      1,
      length = tt
    )

    dens <- truncnorm::dtruncnorm(
      x_grid,
      a = -1,
      b = 1,
      mean = alpha1_star,
      sd = tau1
    )

    dens <- dens / sum(dens)

    sum(
      (x_grid > c) *
        dens
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


sim_b2_gauss_param12 <- function(
    zeta,
    alpha,
    r,
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    alpha1_star_r,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 50,
    tt2 = 50) {

  MASS::ginv(
    sim_M2_xz_gauss_param12(
      beta,
      x_a,
      alpha1_star,
      alpha2_star,
      sigma,
      tau1,
      tau2,
      tt,
      tt2
    )
  ) %*%
    sim_c2_xz_gauss_param12(
      zeta,
      alpha,
      r,
      beta,
      x_a,
      alpha1_star,
      alpha2_star,
      alpha1_star_r,
      sigma,
      tau1,
      tau2,
      tt,
      tt2
    )
}


sim_b3_gauss_param12 <- function(
    zeta,
    alpha,
    r,
    beta,
    x_a,
    alpha1_star,
    alpha2_star,
    alpha1_star_r,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 50,
    tt2 = 50,
    tt3 = 50) {

  cc <- sim_gauss(tt2)


  temp2 <- function(c, x) {

    sapply(
      sigma * cc$x +
        sum(
          c(1, x) *
            beta
        ),
      function(y) {

        r(
          y,
          min(
            x,
            c
          ),
          x <= c,
          beta,
          alpha1_star_r,
          tau1
        ) <= zeta
      }
    ) %*%
      cc$w
  }


  temp3 <- function(c) {

    x_grid <- seq(
      -1,
      1,
      length.out = tt
    )

    dens_x <- truncnorm::dtruncnorm(
      x_grid,
      a = -1,
      b = 1,
      mean = alpha1_star,
      sd = tau1
    )

    RHS <- sapply(
      seq_len(tt),
      function(i) {
        temp2(
          c,
          x_grid[i]
        )
      }
    )

    RHS %*%
      dens_x /
      sum(dens_x)
  }


  temp4 <- function() {

    c_grid <- seq(
      -1,
      1,
      length = tt3
    )

    dens <- truncnorm::dtruncnorm(
      c_grid,
      a = -1,
      b = 1,
      mean = alpha2_star,
      sd = tau2
    )

    dens <- dens / sum(dens)

    sapply(
      c_grid,
      function(c) {
        temp3(c)
      }
    ) %*%
      dens
  }


  temp4() -
    (1 - alpha)
}


# Observed-data estimating function for zeta.

sim_trans_phi_eff_gauss_param12 <- function(
    zeta,
    alpha,
    r,
    beta,
    y,
    w,
    delta,
    x_a,
    b1,
    b2,
    b3,
    alpha1_star,
    alpha2_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1) {

  if (delta == 1) {

    temp <- function(x_a, x, y) {

      p <- truncnorm::dtruncnorm(
        x_a,
        a = -1,
        b = 1,
        mean = alpha2_star,
        sd = tau2
      )

      num <- (x_a >= x) * p
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


    b1_w +
      t(b2) %*%
      temp(
        x_a,
        w,
        y
      ) -
      b3

  } else {

    temp <- function(x_a, c, y) {

      v_x <- 1 / (
        beta[2]^2 / sigma^2 +
          1 / tau1^2
      )

      eta_x <- v_x * (
        beta[2] *
          (y - beta[1]) /
          sigma^2 +
          alpha1_star /
          tau1^2
      )

      p <- truncnorm::dtruncnorm(
        x_a,
        a = -1,
        b = 1,
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


    b2_w <- approx(
      x_a,
      b2,
      w,
      rule = 2
    )$y


    t(b1) %*%
      temp(
        x_a,
        w,
        y
      ) +
      b2_w -
      b3
  }
}


# =============================================================================
# Half-length estimation
# =============================================================================

# Find the root by linear interpolation or extrapolation.

sim_find_sol <- function(
    x_vec,
    y_vec) {

  sign_vec <- sign(y_vec)

  sign_change <- which(
    diff(sign_vec) != 0
  )


  if (length(sign_change) == 0) {

    if (all(y_vec > 0)) {

      idx <- which.min(y_vec)

      x1 <- x_vec[idx]
      x2 <- x_vec[idx + 1]

      y1 <- y_vec[idx]
      y2 <- y_vec[idx + 1]

      x1 -
        y1 *
        (x2 - x1) /
        (y2 - y1)

    } else {

      idx <- which.max(y_vec)

      x1 <- x_vec[idx - 1]
      x2 <- x_vec[idx]

      y1 <- y_vec[idx - 1]
      y2 <- y_vec[idx]

      x2 -
        y2 *
        (x2 - x1) /
        (y2 - y1)
    }

  } else {

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
}


# Bilinear interpolation over the beta grid.

interpolate <- function(
    x_array,
    y_array,
    x) {

  x_vec1 <- x_array[1, ]
  x_vec2 <- x_array[2, ]

  x1 <- x[1]
  x2 <- x[2]

  n1 <- length(x_vec1)
  n2 <- length(x_vec2)

  idx1 <- max(
    which(
      x_vec1 <= x1
    )
  )

  idx2 <- max(
    which(
      x_vec2 <= x2
    )
  )


  if (x1 < min(x_vec1)) {
    idx1 <- 1
  }

  if (x2 < min(x_vec2)) {
    idx2 <- 1
  }

  if (idx1 == n1) {
    idx1 <- n1 - 1
  }

  if (idx2 == n2) {
    idx2 <- n2 - 1
  }


  x11 <- x_vec1[idx1]
  x12 <- x_vec1[idx1 + 1]

  x21 <- x_vec2[idx2]
  x22 <- x_vec2[idx2 + 1]

  y11 <- y_array[idx1, idx2]
  y12 <- y_array[idx1, idx2 + 1]

  y21 <- y_array[idx1 + 1, idx2]
  y22 <- y_array[idx1 + 1, idx2 + 1]


  y1 <- y11 +
    (
      y12 -
        y11
    ) *
    (
      x2 -
        x21
    ) /
    (
      x22 -
        x21
    )

  y2 <- y21 +
    (
      y22 -
        y21
    ) *
    (
      x2 -
        x21
    ) /
    (
      x22 -
        x21
    )


  y1 +
    (
      y2 -
        y1
    ) *
    (
      x1 -
        x11
    ) /
    (
      x12 -
        x11
    )
}


# Estimate zeta using interpolation over the beta grid.

get_zeta_param_param12_int <- function(
    k,
    n,
    d,
    alpha,
    r,
    beta,
    alpha1,
    alpha2,
    b1_array,
    b2_array,
    b3_array,
    beta_array,
    beta_temp,
    zeta_seq,
    alpha1_star,
    alpha2_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    m = 20) {

  data_k <- data_generating(
    k,
    n,
    d,
    beta,
    alpha1,
    alpha2,
    sigma,
    tau1,
    tau2
  )

  y_data <- data_k$y_data
  w_data <- data_k$w_data
  delta_data <- data_k$delta_data

  x_a <- seq(
    -1,
    1,
    length.out = m
  )

  func_vals <- rep(
    0,
    length(zeta_seq)
  )


  for (zeta_idx in seq_along(zeta_seq)) {

    zeta <- zeta_seq[zeta_idx]

    b1 <- rep(
      0,
      m
    )

    b2 <- rep(
      0,
      m
    )


    for (i in seq_len(
      length(
        b1_array[
          zeta_idx,
          ,
          1,
          1
        ]
      )
    )) {

      b1[i] <- interpolate(
        beta_array,
        b1_array[
          zeta_idx,
          i,
          ,
        ],
        beta_temp
      )

      b2[i] <- interpolate(
        beta_array,
        b2_array[
          zeta_idx,
          i,
          ,
        ],
        beta_temp
      )
    }


    b3 <- interpolate(
      beta_array,
      b3_array[
        zeta_idx,
        ,
      ],
      beta_temp
    )


    val <- 0

    n_temp <- length(y_data)

    for (i in seq_len(n_temp)) {

      val <- val +
        sim_trans_phi_eff_gauss_param12(
          zeta,
          alpha,
          r,
          beta_temp,
          y_data[i],
          w_data[i],
          delta_data[i],
          x_a,
          b1,
          b2,
          b3,
          alpha1_star,
          alpha2_star,
          sigma,
          tau1,
          tau2
        )
    }

    func_vals[zeta_idx] <- val
  }


  list(
    vec = rbind(
      zeta_seq,
      func_vals
    ),
    sol = sim_find_sol(
      zeta_seq,
      func_vals
    )
  )
}


# Construct the zeta search grid.

sim_get_zeta_seq_r <- function(
    alpha,
    alpha1,
    alpha2,
    r,
    alpha1_star_r = alpha1,
    beta = c(0, 3),
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    d = 1,
    length = 5,
    N = 1000,
    reps = 10) {

  zeta_r_vals <- rep(
    0,
    reps
  )


  for (k in seq_len(reps)) {

    data_mc <- data_generating(
      k,
      N,
      d,
      beta,
      alpha1,
      alpha2,
      sigma = sigma,
      tau1 = tau1,
      tau2 = tau2
    )

    y_data <- data_mc$y_data
    w_data <- data_mc$w_data
    delta_data <- data_mc$delta_data


    zeta_r_vals[k] <- quantile(
      sapply(
        seq_len(N),
        function(i) {

          r(
            y_data[i],
            w_data[i],
            delta_data[i],
            beta,
            alpha1_star_r,
            tau1
          )
        }
      ),
      1 - alpha
    )
  }


  seq(
    mean(zeta_r_vals) -
      sd(zeta_r_vals),
    mean(zeta_r_vals) +
      sd(zeta_r_vals),
    length.out = length
  )
}


# Estimate zeta at the true beta.

get_zeta_param_param12_int_true <- function(
    k,
    n,
    d,
    alpha,
    r,
    beta,
    alpha1,
    alpha2,
    b1_array,
    b2_array,
    b3_array,
    zeta_seq,
    alpha1_star,
    alpha2_star,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    m = 20) {

  data_k <- data_generating(
    k,
    n,
    d,
    beta,
    alpha1,
    alpha2,
    sigma,
    tau1,
    tau2
  )

  y_data <- data_k$y_data
  w_data <- data_k$w_data
  delta_data <- data_k$delta_data

  x_a <- seq(
    -1,
    1,
    length.out = m
  )

  func_vals <- rep(
    0,
    length(zeta_seq)
  )


  for (zeta_idx in seq_along(zeta_seq)) {

    zeta <- zeta_seq[zeta_idx]

    b1 <- rep(
      0,
      m
    )

    b2 <- rep(
      0,
      m
    )


    for (i in seq_len(
      length(
        b1_array[
          zeta_idx,
          ,
          1,
          1
        ]
      )
    )) {

      b1[i] <- b1_array[
        zeta_idx,
        i,
        ,
      ]

      b2[i] <- b2_array[
        zeta_idx,
        i,
        ,
      ]
    }


    b3 <- b3_array[
      zeta_idx,
      ,
    ]


    val <- 0

    n_temp <- length(y_data)

    for (i in seq_len(n_temp)) {

      val <- val +
        sim_trans_phi_eff_gauss_param12(
          zeta,
          alpha,
          r,
          beta,
          y_data[i],
          w_data[i],
          delta_data[i],
          x_a,
          b1,
          b2,
          b3,
          alpha1_star,
          alpha2_star,
          sigma,
          tau1,
          tau2
        )
    }

    func_vals[zeta_idx] <- val
  }


  list(
    vec = rbind(
      zeta_seq,
      func_vals
    ),
    sol = sim_find_sol(
      zeta_seq,
      func_vals
    )
  )
}


# =============================================================================
# Shared simulation settings
# =============================================================================

# True alpha2 values corresponding to the five censoring levels.

.censor_alpha2 <- c(
  low = 2,
  lowmid = 1,
  mid = 0,
  highmid = -1,
  high = -2
)


# Misspecified eta2 mean.

.gamma_star <- function(gamma) {

  if (gamma != 0) {
    0
  } else {
    2
  }
}


# Misspecified eta1 mean.

.ALPHA1_MIS <- -2
