# =============================================================================
# Conformal prediction simulations
#
# Run from the package root:
#
#   Rscript sim/conformal.R
#
# Outputs are saved in:
#
#   sim/Rdata_files/
#
# The three comparison methods are:
#   1. Split conformal prediction
#   2. Full conformal prediction
#   3. Jackknife+
#
# Each method is evaluated under five censoring levels and three residual
# definitions: r1, r2, and r1*.
# =============================================================================


# ---- setup -------------------------------------------------------------------

source("sim/core.R")

RDATA_DIR <- "sim/Rdata_files"
dir.create(RDATA_DIR, showWarnings = FALSE, recursive = TRUE)


# ---- simulation settings -----------------------------------------------------

M <- 1000
n <- 1000
d <- 1

beta <- c(0, 3)

alpha1 <- 0
alpha <- 0.1

sigma <- 4
tau1 <- 1
tau2 <- 1

tt <- 20
m <- 20

N <- 10000
seed_offset <- 98765

split_rate <- 0.5

# Center used for the r1* residual.
alpha1_star_r <- -2


# =============================================================================
# Split conformal prediction
# =============================================================================

# k: replicate index.
# n, d: sample size and covariate dimension.
# beta, alpha1, alpha2: data-generating parameters.
# alpha: miscoverage level.
# split_rate: fraction of observations used to estimate beta.
# alpha1_star_r: center used for the r1* residual.
# N: size of the independent test set used to estimate prediction coverage.
# seed_offset: offset used to generate the common test set.
# sigma, tau1, tau2, tt, m: model and numerical settings.

sim_get_zeta_beta_param_param12_cp <- function(
    k, n, d,
    beta, alpha1, alpha2, alpha,
    split_rate = 0.5,
    alpha1_star_r = -2,
    N = 10000,
    seed_offset = 98765,
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    tt = 20,
    m = 20) {

  data_k <- data_generating(
    k, n, d, beta, alpha1, alpha2,
    sigma, tau1, tau2
  )

  y_data <- data_k$y_data
  w_data <- data_k$w_data
  delta_data <- data_k$delta_data

  n1 <- floor(n * split_rate)
  x_a <- seq(-1, 1, length.out = m)

  set.seed(k)

  result <- nleqslv::nleqslv(
    beta + rnorm(d + 1) * 0.1,
    sim_pe_gauss_param12,
    y_data = y_data[1:n1],
    w_data = w_data[1:n1],
    delta_data = delta_data[1:n1],
    x_a = x_a,
    alpha1_star = alpha1,
    alpha2_star = alpha2,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2,
    tt = tt
  )

  beta_hat <- result$x


  # Calibration residuals.

  idx <- (n1 + 1):n

  res_r1 <- sim_r1(
    y_data[idx],
    w_data[idx],
    delta_data[idx],
    beta_hat,
    alpha1,
    tau1
  )

  res_r2 <- sim_r2(
    y_data[idx],
    w_data[idx],
    delta_data[idx],
    beta_hat,
    alpha1,
    tau1
  )

  res_r1star <- sim_r1(
    y_data[idx],
    w_data[idx],
    delta_data[idx],
    beta_hat,
    alpha1_star_r,
    tau1
  )


  # Estimated half-lengths.

  zeta_r1 <- unname(
    quantile(c(res_r1, Inf), 1 - alpha)
  )

  zeta_r2 <- unname(
    quantile(c(res_r2, Inf), 1 - alpha)
  )

  zeta_r1star <- unname(
    quantile(c(res_r1star, Inf), 1 - alpha)
  )


  # Prediction coverage on an independent test set.

  test <- data_generating(
    k + seed_offset,
    N,
    d,
    beta,
    alpha1,
    alpha2,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2
  )

  t_r1 <- sim_r1_vec(
    test$y_data,
    test$w_data,
    test$delta_data,
    beta_hat,
    alpha1,
    tau1
  )

  t_r2 <- sim_r2_vec(
    test$y_data,
    test$w_data,
    test$delta_data,
    beta_hat,
    alpha1,
    tau1
  )

  t_r1star <- sim_r1_vec(
    test$y_data,
    test$w_data,
    test$delta_data,
    beta_hat,
    alpha1_star_r,
    tau1
  )


  c(
    beta_hat = beta_hat,
    cvg_r1 = mean(t_r1 <= zeta_r1),
    cvg_r2 = mean(t_r2 <= zeta_r2),
    cvg_r1star = mean(t_r1star <= zeta_r1star),
    zeta_r1 = zeta_r1,
    zeta_r2 = zeta_r2,
    zeta_r1star = zeta_r1star
  )
}


# ---- run split conformal prediction -----------------------------------------

for (censor in names(.censor_alpha2)) {

  alpha2 <- .censor_alpha2[[censor]]

  message(
    "Split conformal: ", censor,
    "  (alpha2 = ", alpha2, ")"
  )

  res <- vector("list", M)

  for (k in seq_len(M)) {

    set.seed(k)

    res[[k]] <- sim_get_zeta_beta_param_param12_cp(
      k,
      n,
      d,
      beta,
      alpha1,
      alpha2,
      alpha,
      split_rate = split_rate,
      alpha1_star_r = alpha1_star_r,
      N = N,
      seed_offset = seed_offset,
      sigma = sigma,
      tau1 = tau1,
      tau2 = tau2,
      tt = tt,
      m = m
    )

    if (k %% 50 == 0) {
      message("  replicate ", k, "/", M)
    }
  }

  res <- do.call(rbind, res)


  # Beta estimates.

  result_beta <- t(
    res[, 1:length(beta), drop = FALSE]
  )

  save(
    result_beta,
    file = file.path(
      RDATA_DIR,
      paste0(
        "beta_param12_",
        censor,
        "_cp.Rdata"
      )
    )
  )


  # r1.

  zeta_sol <- res[, "zeta_r1"]

  save(
    zeta_sol,
    file = file.path(
      RDATA_DIR,
      paste0(
        "zeta_param12_r1_",
        alpha, "_",
        censor,
        "_cp.Rdata"
      )
    )
  )

  pred_cvg <- res[, "cvg_r1"]

  save(
    pred_cvg,
    file = file.path(
      RDATA_DIR,
      paste0(
        "pred_cvg_r1_",
        alpha, "_",
        censor,
        "_cp.Rdata"
      )
    )
  )


  # r2.

  zeta_sol <- res[, "zeta_r2"]

  save(
    zeta_sol,
    file = file.path(
      RDATA_DIR,
      paste0(
        "zeta_param12_r2_",
        alpha, "_",
        censor,
        "_cp.Rdata"
      )
    )
  )

  pred_cvg <- res[, "cvg_r2"]

  save(
    pred_cvg,
    file = file.path(
      RDATA_DIR,
      paste0(
        "pred_cvg_r2_",
        alpha, "_",
        censor,
        "_cp.Rdata"
      )
    )
  )


  # r1*.

  zeta_sol <- res[, "zeta_r1star"]

  save(
    zeta_sol,
    file = file.path(
      RDATA_DIR,
      paste0(
        "zeta_param12_r1_",
        alpha, "_",
        censor,
        "_r1star_cp.Rdata"
      )
    )
  )

  pred_cvg <- res[, "cvg_r1star"]

  save(
    pred_cvg,
    file = file.path(
      RDATA_DIR,
      paste0(
        "pred_cvg_r1_",
        alpha, "_",
        censor,
        "_r1star_cp.Rdata"
      )
    )
  )
}


# =============================================================================
# Full conformal prediction
# =============================================================================

# k: replicate index.
# n, d: sample size and covariate dimension.
# beta, alpha1, alpha2: data-generating parameters.
# alpha: miscoverage level.
# N: size of the independent test set.
# alpha1_star_r: center used for the r1* residual.
# seed_offset: offset used to generate the common test set.
# sigma, tau1, tau2: model settings.

sim_get_pred_cvg_fcp <- function(
    k, n, d,
    beta, alpha1, alpha2, alpha,
    N = 10000,
    alpha1_star_r = -2,
    seed_offset = 98765,
    sigma = 4,
    tau1 = 1,
    tau2 = 1) {

  data_k <- data_generating(
    k, n, d, beta, alpha1, alpha2,
    sigma, tau1, tau2
  )

  y_data <- data_k$y_data
  w_data <- data_k$w_data
  delta_data <- data_k$delta_data


  newdata <- data_generating(
    k + seed_offset,
    N,
    d,
    beta,
    alpha1,
    alpha2,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2
  )

  y_new <- newdata$y_data
  w_new <- newdata$w_data
  delta_new <- newdata$delta_data


  pred_cvg_list <- lapply(
    seq_len(length(y_new)),
    function(j) {

      y_aug <- c(
        y_new[j],
        y_data
      )

      w_aug <- c(
        w_new[j],
        w_data
      )

      delta_aug <- c(
        delta_new[j],
        delta_data
      )


      # Complete-case estimate using the augmented sample.

      cc_idx <- delta_aug == 1

      X <- cbind(
        1,
        w_aug[cc_idx]
      )

      Y <- y_aug[cc_idx]

      beta_cc <- solve(
        crossprod(X),
        crossprod(X, Y)
      )


      # Residuals in the augmented sample.

      r1_vector <- sim_r1_vec(
        y_aug,
        w_aug,
        delta_aug,
        beta_cc,
        alpha1,
        tau1
      )

      r2_vector <- sim_r2_vec(
        y_aug,
        w_aug,
        delta_aug,
        beta_cc,
        alpha1,
        tau1
      )

      r1star_vector <- sim_r1_vec(
        y_aug,
        w_aug,
        delta_aug,
        beta_cc,
        alpha1_star_r,
        tau1
      )


      c(
        cvg_r1 = unname(
          r1_vector[1] <=
            quantile(r1_vector, 1 - alpha)
        ),

        cvg_r2 = unname(
          r2_vector[1] <=
            quantile(r2_vector, 1 - alpha)
        ),

        cvg_r1star = unname(
          r1star_vector[1] <=
            quantile(r1star_vector, 1 - alpha)
        ),

        zeta_r1 = unname(
          quantile(r1_vector, 1 - alpha)
        ),

        zeta_r2 = unname(
          quantile(r2_vector, 1 - alpha)
        ),

        zeta_r1star = unname(
          quantile(r1star_vector, 1 - alpha)
        )
      )
    }
  )

  pred_cvg_list <- do.call(
    "rbind",
    pred_cvg_list
  )

  apply(
    pred_cvg_list,
    2,
    mean
  )
}


# ---- run full conformal prediction ------------------------------------------

for (censor in names(.censor_alpha2)) {

  alpha2 <- .censor_alpha2[[censor]]

  message(
    "Full conformal: ", censor,
    "  (alpha2 = ", alpha2, ")"
  )

  res <- vector("list", M)

  for (k in seq_len(M)) {

    set.seed(k)

    res[[k]] <- sim_get_pred_cvg_fcp(
      k,
      n,
      d,
      beta,
      alpha1,
      alpha2,
      alpha,
      N = N,
      alpha1_star_r = alpha1_star_r,
      seed_offset = seed_offset,
      sigma = sigma,
      tau1 = tau1,
      tau2 = tau2
    )

    if (k %% 50 == 0) {
      message("  replicate ", k, "/", M)
    }
  }

  res <- do.call(rbind, res)


  # r1.

  pred_cvg <- res[, "cvg_r1"]

  save(
    pred_cvg,
    file = file.path(
      RDATA_DIR,
      paste0(
        "pred_cvg_r1_",
        alpha, "_",
        censor,
        "_fcp.Rdata"
      )
    )
  )

  zeta_sol <- res[, "zeta_r1"]

  save(
    zeta_sol,
    file = file.path(
      RDATA_DIR,
      paste0(
        "zeta_r1_",
        alpha, "_",
        censor,
        "_fcp.Rdata"
      )
    )
  )


  # r2.

  pred_cvg <- res[, "cvg_r2"]

  save(
    pred_cvg,
    file = file.path(
      RDATA_DIR,
      paste0(
        "pred_cvg_r2_",
        alpha, "_",
        censor,
        "_fcp.Rdata"
      )
    )
  )

  zeta_sol <- res[, "zeta_r2"]

  save(
    zeta_sol,
    file = file.path(
      RDATA_DIR,
      paste0(
        "zeta_r2_",
        alpha, "_",
        censor,
        "_fcp.Rdata"
      )
    )
  )


  # r1*.

  pred_cvg <- res[, "cvg_r1star"]

  save(
    pred_cvg,
    file = file.path(
      RDATA_DIR,
      paste0(
        "pred_cvg_r1_",
        alpha, "_",
        censor,
        "_r1star_fcp.Rdata"
      )
    )
  )

  zeta_sol <- res[, "zeta_r1star"]

  save(
    zeta_sol,
    file = file.path(
      RDATA_DIR,
      paste0(
        "zeta_r1_",
        alpha, "_",
        censor,
        "_r1star_fcp.Rdata"
      )
    )
  )
}


# =============================================================================
# Jackknife+
# =============================================================================

# k: replicate index.
# n, d: sample size and covariate dimension.
# beta, alpha1, alpha2: data-generating parameters.
# alpha: miscoverage level.
# N: size of the independent test set.
# alpha1_star_r: center used for the r1* residual.
# seed_offset: offset used to generate the common test set.
# sigma, tau1, tau2: model settings.

sim_get_pred_cvg_jackknife <- function(
    k, n, d,
    beta, alpha1, alpha2, alpha,
    N = 10000,
    alpha1_star_r = -2,
    seed_offset = 98765,
    sigma = 4,
    tau1 = 1,
    tau2 = 1) {

  data_k <- data_generating(
    k, n, d, beta, alpha1, alpha2,
    sigma, tau1, tau2
  )

  y_data <- data_k$y_data
  w_data <- data_k$w_data
  delta_data <- data_k$delta_data


  newdata <- data_generating(
    k + seed_offset,
    N,
    d,
    beta,
    alpha1,
    alpha2,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2
  )

  y_new <- newdata$y_data
  w_new <- newdata$w_data
  delta_new <- newdata$delta_data


  m_vals <- lapply(
    seq_len(length(y_data)),
    function(i) {

      cc_idx <- delta_data[-i] == 1

      X <- cbind(
        1,
        w_data[-i][cc_idx]
      )

      Y <- y_data[-i][cc_idx]

      beta_cc <- solve(
        crossprod(X),
        crossprod(X, Y)
      )


      m1_vector <- sim_m1_vec(
        w_new,
        delta_new,
        beta_cc,
        alpha1,
        tau1
      )

      m2_vector <- sim_m0_vec(
        w_new,
        beta_cc
      )

      m1star_vector <- sim_m1_vec(
        w_new,
        delta_new,
        beta_cc,
        alpha1_star_r,
        tau1
      )


      r1_val <- sim_r1(
        y_data[i],
        w_data[i],
        delta_data[i],
        beta_cc,
        alpha1,
        tau1
      )

      r2_val <- sim_r2(
        y_data[i],
        w_data[i],
        delta_data[i],
        beta_cc,
        alpha1,
        tau1
      )

      r1star_val <- sim_r1(
        y_data[i],
        w_data[i],
        delta_data[i],
        beta_cc,
        alpha1_star_r,
        tau1
      )


      cbind(
        m1_vector - r1_val,
        m1_vector + r1_val,
        m2_vector - r2_val,
        m2_vector + r2_val,
        m1star_vector - r1star_val,
        m1star_vector + r1star_val
      )
    }
  )

  m_vals <- simplify2array(m_vals)


  lower <- apply(
    m_vals[, c(1, 3, 5), ],
    c(1, 2),
    function(x) {
      quantile(
        x,
        probs = alpha
      )
    }
  )

  upper <- apply(
    m_vals[, c(2, 4, 6), ],
    c(1, 2),
    function(x) {
      quantile(
        x,
        probs = 1 - alpha
      )
    }
  )


  cvg <- sapply(
    1:3,
    function(j) {
      mean(
        (y_new >= lower[, j]) &
          (y_new <= upper[, j])
      )
    }
  )

  names(cvg) <- c(
    "cvg_r1",
    "cvg_r2",
    "cvg_r1star"
  )


  zeta_sol <- apply(
    (upper - lower) / 2,
    2,
    mean
  )

  names(zeta_sol) <- c(
    "zeta_r1",
    "zeta_r2",
    "zeta_r1star"
  )


  c(
    cvg,
    zeta_sol
  )
}


# ---- run Jackknife+ ----------------------------------------------------------

for (censor in names(.censor_alpha2)) {

  alpha2 <- .censor_alpha2[[censor]]

  message(
    "Jackknife+: ", censor,
    "  (alpha2 = ", alpha2, ")"
  )

  res <- vector("list", M)

  for (k in seq_len(M)) {

    set.seed(k)

    res[[k]] <- sim_get_pred_cvg_jackknife(
      k,
      n,
      d,
      beta,
      alpha1,
      alpha2,
      alpha,
      N = N,
      alpha1_star_r = alpha1_star_r,
      seed_offset = seed_offset,
      sigma = sigma,
      tau1 = tau1,
      tau2 = tau2
    )

    if (k %% 50 == 0) {
      message("  replicate ", k, "/", M)
    }
  }

  res <- do.call(rbind, res)


  # r1.

  pred_cvg <- res[, "cvg_r1"]

  save(
    pred_cvg,
    file = file.path(
      RDATA_DIR,
      paste0(
        "pred_cvg_r1_",
        alpha, "_",
        censor,
        "_jackknife.Rdata"
      )
    )
  )

  zeta_sol <- res[, "zeta_r1"]

  save(
    zeta_sol,
    file = file.path(
      RDATA_DIR,
      paste0(
        "zeta_r1_",
        alpha, "_",
        censor,
        "_jackknife.Rdata"
      )
    )
  )


  # r2.

  pred_cvg <- res[, "cvg_r2"]

  save(
    pred_cvg,
    file = file.path(
      RDATA_DIR,
      paste0(
        "pred_cvg_r2_",
        alpha, "_",
        censor,
        "_jackknife.Rdata"
      )
    )
  )

  zeta_sol <- res[, "zeta_r2"]

  save(
    zeta_sol,
    file = file.path(
      RDATA_DIR,
      paste0(
        "zeta_r2_",
        alpha, "_",
        censor,
        "_jackknife.Rdata"
      )
    )
  )


  # r1*.

  pred_cvg <- res[, "cvg_r1star"]

  save(
    pred_cvg,
    file = file.path(
      RDATA_DIR,
      paste0(
        "pred_cvg_r1_",
        alpha, "_",
        censor,
        "_r1star_jackknife.Rdata"
      )
    )
  )

  zeta_sol <- res[, "zeta_r1star"]

  save(
    zeta_sol,
    file = file.path(
      RDATA_DIR,
      paste0(
        "zeta_r1_",
        alpha, "_",
        censor,
        "_r1star_jackknife.Rdata"
      )
    )
  )
}
