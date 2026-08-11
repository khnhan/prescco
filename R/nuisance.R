# =============================================================================
#  Nuisance-model estimation: X | Z and C | Z
#  ---------------------------------------------------------------------------
#  Maximum-likelihood estimation of the working truncated-normal nuisance
#  models for the event-time covariate X and the censoring time C given the
#  fully-observed covariates Z = (z_c, z_d):
#     X | Z ~ TN( (1, z_c, z_d) alpha1, tau1^2 ; [w_min, w_max] )
#     C | Z ~ TN( (1, z_c, z_d) alpha2, tau2^2 ; [w_min, w_max] )
#  Censored observations (delta == 0) contribute survival terms P(X > w | Z);
#  uncensored observations contribute the truncated-normal density. Both z_c
#  and z_d are optional; with neither, the mean reduces to an intercept.
#
#  Package imports are declared centrally in R/censcovpred-package.R.
# =============================================================================

#' Maximum-likelihood fit of the eta1 model X | Z
#'
#' @param w_data,delta_data Observed times and event indicators (length n).
#' @param z_c_data Optional continuous covariates (length-n vector or n-by-p matrix).
#' @param z_d_data Optional discrete covariates (length-n vector or n-by-q matrix).
#' @param w_min,w_max Truncation bounds for X.
#' @return A named numeric vector: the coefficient estimates alpha1_hat
#'   (intercept, then z_c then z_d) followed by tau1_hat.
#' @examples
#' set.seed(1)
#' n  <- 200
#' x  <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' cc <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w  <- pmin(x, cc)
#' delta <- as.integer(x <= cc)
#' y  <- 3 * x + rnorm(n)          # beta = c(0, 3)
#'
#' # Recovers the mean and SD of X from censored data; here the truth
#' # is mean 0, sd 1.
#' find_alpha1_MLE(w, delta, w_min = -1, w_max = 1)
#' @export
find_alpha1_MLE = function(w_data, delta_data,
                           z_c_data = NULL, z_d_data = NULL,
                           w_min, w_max) {
  
  n = length(w_data)
  if (length(delta_data) != n) stop("delta_data must have length n")
  
  ## ---------- Handle z_c_data (continuous) ----------
  if (is.null(z_c_data)) {
    z_c_mat = matrix(0, nrow = n, ncol = 0)
    p_c     = 0
  } else if (is.null(dim(z_c_data))) {
    if (length(z_c_data) != n) stop("z_c_data must have length n")
    z_c_mat = matrix(as.numeric(z_c_data), ncol = 1)
    p_c     = 1
  } else {
    if (nrow(z_c_data) != n) stop("z_c_data must have n rows")
    z_c_mat = as.matrix(z_c_data)
    p_c     = ncol(z_c_mat)
  }
  
  ## ---------- Handle z_d_data (discrete / arbitrary levels) ----------
  if (is.null(z_d_data)) {
    z_d_mat = matrix(0, nrow = n, ncol = 0)
    p_d     = 0
  } else if (is.null(dim(z_d_data))) {
    if (length(z_d_data) != n) stop("z_d_data must have length n")
    z_d_mat = matrix(z_d_data, ncol = 1)
    p_d     = 1
  } else {
    if (nrow(z_d_data) != n) stop("z_d_data must have n rows")
    z_d_mat = as.matrix(z_d_data)
    p_d     = ncol(z_d_mat)
  }
  
  ## Dimension of alpha1: intercept + all z_c + all z_d
  p_alpha1 = 1 + p_c + p_d
  
  ## ---------- Single-observation log-likelihood ----------
  log_llhd_one = function(alpha1, tau1, w, delta, z_c_i, z_d_i) {
    z_c_i = as.numeric(z_c_i)
    z_d_i = as.numeric(z_d_i)
    
    mu_i = sum(c(1, z_c_i, z_d_i) * alpha1)
    
    if (delta == 1) {
      # log f_trunc(w | Z) for X|Z ~ TN(mu_i, tau1^2; [w_min, w_max])
      dens = truncnorm::dtruncnorm(w,
                        a    = w_min,
                        b    = w_max,
                        mean = mu_i,
                        sd   = tau1)
      if (dens <= 0 || !is.finite(dens)) return(-1e12)
      return(log(dens))
    } else {
      # log P(X > w | Z) under truncated normal
      F_trunc_w = truncnorm::ptruncnorm(w,
                             a    = w_min,
                             b    = w_max,
                             mean = mu_i,
                             sd   = tau1)
      surv = 1 - F_trunc_w
      if (surv <= 0 || !is.finite(surv)) return(-1e12)
      return(log(surv))
    }
  }
  
  ## ---------- Total negative log-likelihood ----------
  neg_loglik = function(theta) {
    alpha1   = theta[1:p_alpha1]
    log_tau1 = theta[p_alpha1 + 1]
    tau1     = exp(log_tau1)
    if (tau1 <= 0) return(1e12)
    
    ll = 0
    for (i in 1:n) {
      z_c_i = if (p_c > 0) z_c_mat[i, ] else numeric(0)
      z_d_i = if (p_d > 0) z_d_mat[i, ] else numeric(0)
      
      ll = ll + log_llhd_one(alpha1, tau1,
                             w      = w_data[i],
                             delta  = delta_data[i],
                             z_c_i  = z_c_i,
                             z_d_i  = z_d_i)
    }
    -ll
  }
  
  ## ---------- Initial values ----------
  # Use only delta == 1 (where X is observed exactly) for OLS init
  if (any(delta_data == 1)) {
    idx1       = which(delta_data == 1)
    w_obs      = w_data[idx1]
    z_c_obs    = z_c_mat[idx1, , drop = FALSE]
    z_d_obs    = z_d_mat[idx1, , drop = FALSE]
    design_mat = cbind(
      rep(1, length(idx1)),
      if (p_c > 0) z_c_obs else NULL,
      if (p_d > 0) z_d_obs else NULL
    )
    alpha_init = as.numeric(stats::lm.fit(design_mat, w_obs)$coefficients)
  } else {
    # If no uncensored data, fall back to crude init
    design_mat = cbind(
      1,
      if (p_c > 0) z_c_mat else NULL,
      if (p_d > 0) z_d_mat else NULL
    )
    alpha_init = rep(0, ncol(design_mat))
  }
  
  if (any(is.na(alpha_init)) || length(alpha_init) == 0) {
    alpha_init = rep(0, p_alpha1)
  } else if (length(alpha_init) < p_alpha1) {
    alpha_init = c(alpha_init, rep(0, p_alpha1 - length(alpha_init)))
  }
  
  # Residual sd as starting tau1
  resid_init = w_data - as.vector(
    cbind(
      1,
      if (p_c > 0) z_c_mat else NULL,
      if (p_d > 0) z_d_mat else NULL
    ) %*% alpha_init
  )
  tau1_init  = sd(resid_init)
  if (!is.finite(tau1_init) || tau1_init <= 0) {
    tau1_init = (w_max - w_min) / 4
  }
  
  theta_init = c(alpha_init, log(tau1_init))
  
  ## ---------- Optimize ----------
  opt = optim(theta_init, neg_loglik, 
              method = "L-BFGS-B",
              lower = c(rep(-Inf, p_alpha1), log(1e-3)),
              upper = c(rep( Inf, p_alpha1), log(100)))
  
  alpha1_hat = opt$par[1:p_alpha1]
  tau1_hat   = exp(opt$par[p_alpha1 + 1])
  c(alpha1_hat  = alpha1_hat,
    tau1_hat    = tau1_hat)
}

#' Maximum-likelihood fit of the eta2 model C | Z
#'
#' @param w_data,delta_data Observed times and event indicators (length n).
#' @param z_c_data Optional continuous covariates (length-n vector or n-by-p matrix).
#' @param z_d_data Optional discrete covariates (length-n vector or n-by-q matrix).
#' @param w_min,w_max Truncation bounds for C.
#' @param verbose If \code{TRUE}, report the \code{optim} convergence
#'   message. Defaults to \code{FALSE}.
#' @return A named numeric vector: the coefficient estimates alpha2_hat
#'   (intercept, then z_c then z_d) followed by tau2_hat.
#' @examples
#' set.seed(1)
#' n  <- 200
#' x  <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' cc <- truncnorm::rtruncnorm(n, a = -1, b = 1, mean = 1, sd = 1)
#' w  <- pmin(x, cc)
#' delta <- as.integer(x <= cc)
#' y  <- 3 * x + rnorm(n)          # beta = c(0, 3)
#'
#' # Recovers the mean and SD of the censoring time C; here mean 1, sd 1.
#' find_alpha2_MLE(w, delta, w_min = -1, w_max = 1)
#' @export
find_alpha2_MLE = function(w_data, delta_data,
                           z_c_data = NULL, z_d_data = NULL,
                           w_min = 0, w_max = 12,
                           verbose = FALSE) {
  n = length(w_data)
  if (length(delta_data) != n) stop("delta_data must have length n")
  
  ## ---------- Handle z_c_data (continuous) ----------
  if (is.null(z_c_data)) {
    z_c_mat = matrix(0, nrow = n, ncol = 0)
    p_c     = 0
  } else if (is.null(dim(z_c_data))) {
    if (length(z_c_data) != n) stop("z_c_data must have length n")
    z_c_mat = matrix(as.numeric(z_c_data), ncol = 1)
    p_c     = 1
  } else {
    if (nrow(z_c_data) != n) stop("z_c_data must have n rows")
    z_c_mat = as.matrix(z_c_data)
    p_c     = ncol(z_c_mat)
  }
  
  ## ---------- Handle z_d_data (discrete / arbitrary levels) ----------
  if (is.null(z_d_data)) {
    z_d_mat = matrix(0, nrow = n, ncol = 0)
    p_d     = 0
  } else if (is.null(dim(z_d_data))) {
    if (length(z_d_data) != n) stop("z_d_data must have length n")
    z_d_mat = matrix(z_d_data, ncol = 1)
    p_d     = 1
  } else {
    if (nrow(z_d_data) != n) stop("z_d_data must have n rows")
    z_d_mat = as.matrix(z_d_data)
    p_d     = ncol(z_d_mat)
  }
  
  ## Dimension of alpha2: intercept + all z_c + all z_d
  p_alpha2 = 1 + p_c + p_d
  
  ## ---------- Single-observation log-likelihood ----------
  log_llhd_one = function(alpha2, tau2, w, delta, z_c_i, z_d_i) {
    z_c_i = as.numeric(z_c_i)
    z_d_i = as.numeric(z_d_i)
    
    mu_i = sum(c(1, z_c_i, z_d_i) * alpha2)
    
    # Truncated normal C|Z ~ TN(mu_i, tau2^2; [w_min, w_max])
    
    if (delta == 1) {
      # log P(C > w | Z) under truncated normal
      F_trunc_w = truncnorm::ptruncnorm(w,
                             a    = w_min,
                             b    = w_max,
                             mean = mu_i,
                             sd   = tau2)
      surv = 1 - F_trunc_w
      if (surv <= 0 || !is.finite(surv)) return(-1e12)
      return(log(surv))
    } else {
      # log f_trunc(w | Z)
      dens = truncnorm::dtruncnorm(w,
                        a    = w_min,
                        b    = w_max,
                        mean = mu_i,
                        sd   = tau2)
      if (dens <= 0 || !is.finite(dens)) return(-1e12)
      return(log(dens))
    }
  }
  
  ## ---------- Total negative log-likelihood ----------
  neg_loglik = function(theta) {
    alpha2   = theta[1:p_alpha2]
    log_tau2 = theta[p_alpha2 + 1]
    tau2     = exp(log_tau2)
    if (tau2 <= 0) return(1e12)
    
    ll = 0
    for (i in 1:n) {
      z_c_i = if (p_c > 0) z_c_mat[i, ] else numeric(0)
      z_d_i = if (p_d > 0) z_d_mat[i, ] else numeric(0)
      
      ll = ll + log_llhd_one(alpha2, tau2,
                             w      = w_data[i],
                             delta  = delta_data[i],
                             z_c_i  = z_c_i,
                             z_d_i  = z_d_i)
    }
    -ll
  }
  
  ## ---------- Initial values ----------
  # Use only delta == 0 (where we see exact C) for a crude regression init
  if (any(delta_data == 0)) {
    idx0       = which(delta_data == 0)
    w_obs      = w_data[idx0]
    z_c_obs    = z_c_mat[idx0, , drop = FALSE]
    z_d_obs    = z_d_mat[idx0, , drop = FALSE]
    design_mat = cbind(
      rep(1, length(idx0)),
      if (p_c > 0) z_c_obs else NULL,
      if (p_d > 0) z_d_obs else NULL
    )
    alpha_init = as.numeric(stats::lm.fit(design_mat, w_obs)$coefficients)
  } else {
    # If no uncensored data, fall back to crude init
    design_mat = cbind(
      1,
      if (p_c > 0) z_c_mat else NULL,
      if (p_d > 0) z_d_mat else NULL
    )
    alpha_init = rep(0, ncol(design_mat))
  }
  
  if (any(is.na(alpha_init)) || length(alpha_init) == 0) {
    alpha_init = rep(0, p_alpha2)
  } else if (length(alpha_init) < p_alpha2) {
    alpha_init = c(alpha_init, rep(0, p_alpha2 - length(alpha_init)))
  }
  
  # Residual sd as starting tau2
  # Residual sd as starting tau1
  resid_init = w_data - as.vector(
    cbind(
      1,
      if (p_c > 0) z_c_mat else NULL,
      if (p_d > 0) z_d_mat else NULL
    ) %*% alpha_init
  )
  tau1_init  = sd(resid_init)
  tau2_init  = sd(resid_init)
  if (!is.finite(tau2_init) || tau2_init <= 0) {
    tau2_init = (w_max - w_min) / 4
  }
  
  theta_init = c(alpha_init, log(tau2_init))
  
  ## ---------- Optimize ----------
  opt = optim(theta_init, neg_loglik, 
              method = "L-BFGS-B",
              lower = c(rep(-Inf, p_alpha2), log(1e-3)),
              upper = c(rep( Inf, p_alpha2), log(100)))
  if (verbose) message("optim (alpha2): ", opt$message)
  alpha2_hat = opt$par[1:p_alpha2]
  tau2_hat   = exp(opt$par[p_alpha2 + 1])
  
  c(alpha2_hat  = alpha2_hat,
    tau2_hat    = tau2_hat)
}
