# =============================================================================
#  Conformal-baseline simulation driver for censcovpred
#  ---------------------------------------------------------------------------
#  Runs the three conformal comparison methods (split conformal, full conformal,
#  jackknife+) across censoring levels and writes each method's half-length and
#  empirical coverage rate outputs. Sequential execution with
#  progress messages, matching the semiparametric stage; the original
#  foreach/doSNOW parallelism and all variance/inference outputs are omitted.
# =============================================================================

# ---- Split conformal prediction ---------------------------------------------

#' Split conformal prediction half-lengths
#'
#' Fits beta by SPARCC on the first \code{split_rate} fraction of the data and
#' takes the \code{1 - alpha} empirical quantile of the residuals on the held-out
#' fraction (with the usual \code{+Inf} finite-sample correction) for each of
#' \code{r1}, \code{r2}, and the misspecified-center \code{r1*}.
#'
#' @param k Integer seed / replicate index.
#' @param n,d Sample size and covariate dimension.
#' @param beta,alpha1,alpha2 Data-generating parameters.
#' @param alpha Miscoverage level.
#' @param split_rate Fraction used for fitting (default 0.5).
#' @param alpha1_star_r Misspecified center for the r1* residual (default -2).
#' @param sigma,tau1,tau2,tt,m Model and numerical settings.
#' @return Named numeric vector: \code{beta_hat} (length 2) and
#'   \code{zeta_r1}, \code{zeta_r2}, \code{zeta_r1star}.
sim_get_zeta_beta_param_param12_cp <- function(k, n, d,
                                           beta, alpha1, alpha2, alpha,
                                           split_rate = 0.5,
                                           alpha1_star_r = -2,
                                           N = 10000, seed_offset = 98765,
                                           sigma = 4, tau1 = 1, tau2 = 1,
                                           tt = 20, m = 20) {
  data_k <- data_generating(k, n, d, beta, alpha1, alpha2, sigma, tau1, tau2)
  y_data <- data_k$y_data
  w_data <- data_k$w_data
  delta_data <- data_k$delta_data
  n1 <- floor(n * split_rate)
  x_a <- seq(-1, 1, length.out = m)

  set.seed(k)
  result <- nleqslv::nleqslv(beta + rnorm(d + 1) * 0.1, sim_pe_gauss_param12,
                             y_data = y_data[1:n1],
                             w_data = w_data[1:n1],
                             delta_data = delta_data[1:n1],
                             x_a = x_a,
                             alpha1_star = alpha1, alpha2_star = alpha2,
                             sigma = sigma, tau1 = tau1, tau2 = tau2, tt = tt)
  beta_hat <- result$x

  idx <- (n1 + 1):n
  res_r1     <- sim_r1(y_data[idx], w_data[idx], delta_data[idx], beta_hat, alpha1, tau1)
  res_r2     <- sim_r2(y_data[idx], w_data[idx], delta_data[idx], beta_hat, alpha1, tau1)
  res_r1star <- sim_r1(y_data[idx], w_data[idx], delta_data[idx], beta_hat, alpha1_star_r, tau1)

  zeta_r1     <- unname(quantile(c(res_r1,     Inf), 1 - alpha))
  zeta_r2     <- unname(quantile(c(res_r2,     Inf), 1 - alpha))
  zeta_r1star <- unname(quantile(c(res_r1star, Inf), 1 - alpha))

  ## --- empirical coverage rate on the shared test set ----------------------
  ## Same test set every other method sees: seed k + seed_offset, size N,
  ## same data-generating parameters. Generated after the fit, so it cannot
  ## disturb the RNG stream the fit depends on.
  test <- data_generating(k + seed_offset, N, d, beta, alpha1, alpha2,
                          sigma = sigma, tau1 = tau1, tau2 = tau2)
  t_r1     <- sim_r1_vec(test$y_data, test$w_data, test$delta_data, beta_hat, alpha1, tau1)
  t_r2     <- sim_r2_vec(test$y_data, test$w_data, test$delta_data, beta_hat, alpha1, tau1)
  t_r1star <- sim_r1_vec(test$y_data, test$w_data, test$delta_data, beta_hat, alpha1_star_r, tau1)

  c(beta_hat    = beta_hat,
    cvg_r1      = mean(t_r1     <= zeta_r1),
    cvg_r2      = mean(t_r2     <= zeta_r2),
    cvg_r1star  = mean(t_r1star <= zeta_r1star),
    zeta_r1     = zeta_r1,
    zeta_r2     = zeta_r2,
    zeta_r1star = zeta_r1star)
}


# ---- Full conformal prediction ----------------------------------------------

#' Full conformal prediction coverage rate and half-length
#'
#' For each of \code{N} test points, augments the training data with the test
#' point, refits beta by complete-case least squares, and forms the conformal
#' interval from the augmented residuals. Returns the empirical coverage rate
#' and mean half-length, averaged over the test set, for \code{r1}, \code{r2},
#' and \code{r1*}.
#'
#' @param k Integer seed / replicate index.
#' @param n,d Sample size and covariate dimension.
#' @param beta,alpha1,alpha2 Data-generating parameters.
#' @param alpha Miscoverage level.
#' @param N Test-set size.
#' @param alpha1_star_r Misspecified center for r1* (default -2).
#' @param seed_offset Offset added to the replicate index for the test-set seed.
#' @param sigma,tau1,tau2,tt,m Model and numerical settings.
#' @return Named numeric vector of means: \code{cvg_r1}, \code{cvg_r2},
#'   \code{cvg_r1star}, \code{zeta_r1}, \code{zeta_r2}, \code{zeta_r1star}.
sim_get_pred_cvg_fcp <- function(k, n, d,
                             beta, alpha1, alpha2, alpha,
                             N = 10000, alpha1_star_r = -2, seed_offset = 98765,
                             sigma = 4, tau1 = 1, tau2 = 1, tt = 20, m = 20) {
  data_k <- data_generating(k, n, d, beta, alpha1, alpha2, sigma, tau1, tau2)
  y_data <- data_k$y_data
  w_data <- data_k$w_data
  delta_data <- data_k$delta_data

  newdata <- data_generating(k + seed_offset, N, d, beta, alpha1, alpha2,
                             sigma = sigma, tau1 = tau1, tau2 = tau2)
  y_new <- newdata$y_data
  w_new <- newdata$w_data
  delta_new <- newdata$delta_data

  pred_cvg_list <- lapply(seq_len(length(y_new)), function(j) {
    y_aug     <- c(y_new[j], y_data)
    w_aug     <- c(w_new[j], w_data)
    delta_aug <- c(delta_new[j], delta_data)

    cc_idx  <- (delta_aug == 1)
    X       <- cbind(1, w_aug[cc_idx])
    Y       <- y_aug[cc_idx]
    beta_cc <- solve(crossprod(X), crossprod(X, Y))

    r1_vector     <- sim_r1_vec(y_aug, w_aug, delta_aug, beta_cc, alpha1, tau1)
    r2_vector     <- sim_r2_vec(y_aug, w_aug, delta_aug, beta_cc, alpha1, tau1)
    r1star_vector <- sim_r1_vec(y_aug, w_aug, delta_aug, beta_cc, alpha1_star_r, tau1)

    # unname() matters: quantile() returns a named value ("90%"), and
    # c(zeta_r1 = <named>) composes the two into "zeta_r1.90%", which then
    # breaks the res[, "zeta_r1"] lookup in run_conformal_scenario().
    c(cvg_r1     = unname(r1_vector[1]     <= quantile(r1_vector,     1 - alpha)),
      cvg_r2     = unname(r2_vector[1]     <= quantile(r2_vector,     1 - alpha)),
      cvg_r1star = unname(r1star_vector[1] <= quantile(r1star_vector, 1 - alpha)),
      zeta_r1     = unname(quantile(r1_vector,     1 - alpha)),
      zeta_r2     = unname(quantile(r2_vector,     1 - alpha)),
      zeta_r1star = unname(quantile(r1star_vector, 1 - alpha)))
  })
  pred_cvg_list <- do.call("rbind", pred_cvg_list)
  apply(pred_cvg_list, 2, mean)
}


# ---- Jackknife+ -------------------------------------------------------------

#' Jackknife+ coverage rate and half-length
#'
#' Builds leave-one-out complete-case fits, forms jackknife+ prediction
#' intervals on a test set of size \code{N}, and returns the empirical coverage rate
#' rate and mean half-length for \code{r1}, \code{r2}, and \code{r1*}.
#'
#' @param k Integer seed / replicate index.
#' @param n,d Sample size and covariate dimension.
#' @param beta,alpha1,alpha2 Data-generating parameters.
#' @param alpha Miscoverage level.
#' @param N Test-set size.
#' @param alpha1_star_r Misspecified center for r1* (default -2).
#' @param seed_offset Offset added to the replicate index for the test-set seed.
#' @param sigma,tau1,tau2,tt,m Model and numerical settings.
#' @return Named numeric vector: \code{cvg_r1}, \code{cvg_r2}, \code{cvg_r1star}
#'   and the corresponding mean half-lengths.
sim_get_pred_cvg_jackknife <- function(k, n, d,
                                   beta, alpha1, alpha2, alpha,
                                   N = 10000, alpha1_star_r = -2, seed_offset = 98765,
                                   sigma = 4, tau1 = 1, tau2 = 1, tt = 20, m = 20) {
  data_k <- data_generating(k, n, d, beta, alpha1, alpha2, sigma, tau1, tau2)
  y_data <- data_k$y_data
  w_data <- data_k$w_data
  delta_data <- data_k$delta_data

  newdata <- data_generating(k + seed_offset, N, d, beta, alpha1, alpha2,
                             sigma = sigma, tau1 = tau1, tau2 = tau2)
  y_new <- newdata$y_data
  w_new <- newdata$w_data
  delta_new <- newdata$delta_data

  m_vals <- lapply(seq_len(length(y_data)), function(i) {
    cc_idx  <- (delta_data[-i] == 1)
    X       <- cbind(1, w_data[-i][cc_idx])
    Y       <- y_data[-i][cc_idx]
    beta_cc <- solve(crossprod(X), crossprod(X, Y))

    m1_vector     <- sim_m1_vec(w_new, delta_new, beta_cc, alpha1, tau1)
    m2_vector     <- sim_m0_vec(w_new, beta_cc)
    m1star_vector <- sim_m1_vec(w_new, delta_new, beta_cc, alpha1_star_r, tau1)

    r1_val     <- sim_r1(y_data[i], w_data[i], delta_data[i], beta_cc, alpha1, tau1)
    r2_val     <- sim_r2(y_data[i], w_data[i], delta_data[i], beta_cc, alpha1, tau1)
    r1star_val <- sim_r1(y_data[i], w_data[i], delta_data[i], beta_cc, alpha1_star_r, tau1)

    cbind(m1_vector - r1_val,     m1_vector + r1_val,
          m2_vector - r2_val,     m2_vector + r2_val,
          m1star_vector - r1star_val, m1star_vector + r1star_val)
  })
  m_vals <- simplify2array(m_vals)

  lower <- apply(m_vals[, c(1, 3, 5), ], c(1, 2), function(x) quantile(x, probs = alpha))
  upper <- apply(m_vals[, c(2, 4, 6), ], c(1, 2), function(x) quantile(x, probs = 1 - alpha))

  cvg <- sapply(1:3, function(j) mean((y_new >= lower[, j]) & (y_new <= upper[, j])))
  names(cvg) <- c("cvg_r1", "cvg_r2", "cvg_r1star")
  zeta_sol <- apply((upper - lower) / 2, 2, mean)
  # upper/lower carry no column names, so name these explicitly for the
  # res[, "zeta_<r>"] lookup in run_conformal_scenario().
  names(zeta_sol) <- c("zeta_r1", "zeta_r2", "zeta_r1star")
  c(cvg, zeta_sol)
}


# ---- Simulation driver ------------------------------------------------------

#' Run one conformal-method scenario across replicates
#'
#' @param method One of \code{"cp"} (split conformal), \code{"fcp"} (full
#'   conformal), or \code{"jackknife"}.
#' @param censor_tag Censoring level: \code{"low"}, \code{"lowmid"},
#'   \code{"mid"}, \code{"highmid"}, or \code{"high"}.
#' @param alpha2 True censoring-time mean (gamma) for this level.
#' @param M Number of Monte Carlo replicates.
#' @param alpha Miscoverage level.
#' @param alpha1,beta,sigma,tau1,tau2 Model parameters.
#' @param n,d,m,tt Sample size, dimension, grid resolution, integration size.
#' @param N Test-set size. All three methods evaluate their coverage rate on
#'   this same test set.
#' @param seed_offset Offset added to the replicate index for the test-set
#'   seed, so every method sees identical test data for a given replicate.
#' @param split_rate Fitting fraction for split conformal.
#' @param alpha1_star_r Misspecified center for the r1* variant.
#' @param out_dir Output directory.
#' @param verbose If \code{TRUE}, print progress every 50 replicates.
#' @return Invisibly, a list of the saved objects.
#'
#' @details Every method writes both its half-length and its empirical coverage rate
#'   rate, for each of \code{r1}, \code{r2}, and \code{r1*}. Outputs written to
#'   \code{out_dir}:
#'   \itemize{
#'     \item split conformal: \code{beta_param12_<censor>_cp.Rdata},
#'       \code{zeta_param12_<r>_<alpha>_<censor>[_r1star]_cp.Rdata}, and
#'       \code{pred_cvg_<r>_<alpha>_<censor>[_r1star]_cp.Rdata}.
#'     \item full conformal / jackknife+:
#'       \code{zeta_<r>_<alpha>_<censor>[_r1star]_<method>.Rdata} and
#'       \code{pred_cvg_<r>_<alpha>_<censor>[_r1star]_<method>.Rdata}.
#'   }
#' @export
run_conformal_scenario <- function(method, censor_tag, alpha2,
                                   M = 1000, alpha = 0.1,
                                   alpha1 = 0, beta = c(0, 3),
                                   sigma = 4, tau1 = 1, tau2 = 1,
                                   n = 1000, d = 1, m = 20, tt = 20,
                                   N = 10000, seed_offset = 98765,
                                   split_rate = 0.5,
                                   alpha1_star_r = -2,
                                   out_dir = ".", verbose = TRUE) {
  method <- match.arg(method, c("cp", "fcp", "jackknife"))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  runner <- switch(
    method,
    cp        = function(k) sim_get_zeta_beta_param_param12_cp(
                  k, n, d, beta, alpha1, alpha2, alpha,
                  split_rate = split_rate, alpha1_star_r = alpha1_star_r,
                  N = N, seed_offset = seed_offset,
                  sigma = sigma, tau1 = tau1, tau2 = tau2, tt = tt, m = m),
    fcp       = function(k) sim_get_pred_cvg_fcp(
                  k, n, d, beta, alpha1, alpha2, alpha,
                  N = N, alpha1_star_r = alpha1_star_r, seed_offset = seed_offset,
                  sigma = sigma, tau1 = tau1, tau2 = tau2, tt = tt, m = m),
    jackknife = function(k) sim_get_pred_cvg_jackknife(
                  k, n, d, beta, alpha1, alpha2, alpha,
                  N = N, alpha1_star_r = alpha1_star_r, seed_offset = seed_offset,
                  sigma = sigma, tau1 = tau1, tau2 = tau2, tt = tt, m = m))

  res <- vector("list", M)
  for (k in seq_len(M)) {
    set.seed(k)
    res[[k]] <- runner(k)
    if (verbose && k %% 50 == 0) message("  ", method, " replicate ", k, "/", M)
  }
  res <- do.call(rbind, res)

  a <- alpha
  saved <- list()

  if (method == "cp") {
    # Columns: beta_hat1, beta_hat2, cvg_r1/r2/r1star, zeta_r1/r2/r1star
    result_beta <- t(res[, 1:length(beta), drop = FALSE])
    save(result_beta,
         file = file.path(out_dir, paste0("beta_param12_", censor_tag, "_cp.Rdata")))
    for (rr in c("r1", "r2", "r1star")) {
      infix <- if (rr == "r1star") "_r1star" else ""
      rlab  <- if (rr == "r1star") "r1" else rr

      zeta_sol <- res[, paste0("zeta_", rr)]
      z_name   <- paste0("zeta_param12_", rlab, "_", a, "_", censor_tag, infix, "_cp.Rdata")
      save(zeta_sol, file = file.path(out_dir, z_name))

      pred_cvg <- res[, paste0("cvg_", rr)]
      pc_name  <- paste0("pred_cvg_", rlab, "_", a, "_", censor_tag, infix, "_cp.Rdata")
      save(pred_cvg, file = file.path(out_dir, pc_name))

      saved[[z_name]]  <- zeta_sol
      saved[[pc_name]] <- pred_cvg
    }
  } else {
    # full conformal / jackknife+: columns cvg_r1,cvg_r2,cvg_r1star, zeta_r1,zeta_r2,zeta_r1star
    for (rr in c("r1", "r2", "r1star")) {
      infix <- if (rr == "r1star") "_r1star" else ""
      rlab  <- if (rr == "r1star") "r1" else rr

      pred_cvg <- res[, paste0("cvg_", rr)]
      pc_name  <- paste0("pred_cvg_", rlab, "_", a, "_", censor_tag, infix, "_", method, ".Rdata")
      save(pred_cvg, file = file.path(out_dir, pc_name))

      zeta_sol <- res[, paste0("zeta_", rr)]
      z_name   <- paste0("zeta_", rlab, "_", a, "_", censor_tag, infix, "_", method, ".Rdata")
      save(zeta_sol, file = file.path(out_dir, z_name))

      saved[[pc_name]] <- pred_cvg
      saved[[z_name]]  <- zeta_sol
    }
  }

  invisible(saved)
}

#' Build the conformal-method configuration table
#'
#' One row per (method, censoring level). Split conformal, full conformal, and
#' jackknife+ are each run at all five censoring levels; each run internally
#' produces the r1, r2, and r1* outputs.
#'
#' @return A data frame with columns \code{method}, \code{censor_tag},
#'   \code{alpha2}.
#' @export
build_conformal_scenarios <- function() {
  rows <- list()
  for (method in c("cp", "fcp", "jackknife")) {
    for (censor in c("low", "lowmid", "mid", "highmid", "high")) {
      rows[[length(rows) + 1]] <- data.frame(
        method = method, censor_tag = censor,
        alpha2 = .censor_alpha2[[censor]], stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

#' Run all (or a subset of) conformal-method scenarios
#'
#' @param out_dir Output directory.
#' @param which Optional integer vector selecting rows of
#'   \code{build_conformal_scenarios} (default: all).
#' @param ... Extra arguments forwarded to \code{run_conformal_scenario}
#'   (e.g. \code{M}, \code{N}, \code{verbose}).
#' @return Invisibly, the conformal-scenario table that was run.
#' @export
run_all_conformal_scenarios <- function(out_dir = ".", which = NULL, ...) {
  scenarios <- build_conformal_scenarios()
  if (!is.null(which)) scenarios <- scenarios[which, , drop = FALSE]
  for (i in seq_len(nrow(scenarios))) {
    s <- scenarios[i, ]
    message(sprintf("[conformal %d/%d] %s  censor=%s",
                    i, nrow(scenarios), s$method, s$censor_tag))
    run_conformal_scenario(method = s$method, censor_tag = s$censor_tag,
                           alpha2 = s$alpha2, out_dir = out_dir, ...)
  }
  invisible(scenarios)
}
