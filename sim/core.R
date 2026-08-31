# =============================================================================
#  prescco: semiparametric prediction intervals under a censored covariate
#  ---------------------------------------------------------------------------
#  Core numerical routines for the simulation study (Section: Simulation
#  Studies). This file implements the data-generating process, the efficient
#  score / influence-function machinery solved by SPARCC, the three interval
#  centers (m0/m1/m2) and residuals (r1/r2), and the routines that solve the
#  integral equations for the estimating function b1/b2/b3 and invert them for
#  the half-length zeta.
#
#  Notation follows the paper: X is the time-to-event covariate, C the
#  censoring time, W = min(X, C), Delta = I(X <= C), Y the outcome, and
#  (eta1, eta2) the working nuisance models for X|Z and C|Z. alpha1_star and
#  alpha2_star index the (possibly misspecified) nuisance means; alpha1_star_r
#  indexes the (possibly misspecified) center used inside the residual.
#
#  All heavy numerical defaults (grid sizes tt, tt2, tt3; Gauss node counts)
#  are exposed as arguments and kept at the values used to produce the paper's
#  .Rdata outputs.
# =============================================================================

#' @importFrom stats approx dnorm optim optimize pnorm qnorm quantile rnorm sd setNames var
#' @importFrom utils write.csv
#' @importFrom truncnorm rtruncnorm dtruncnorm ptruncnorm qtruncnorm
#' @importFrom MASS ginv
#' @importFrom nleqslv nleqslv
NULL

#' Generate one simulated data set for the censored-covariate model.
#'
#' @param k Integer seed for reproducibility.
#' @param n Sample size.
#' @param d Covariate dimension (number of fully observed covariates Z); 1 in the simulations with no Z.
#' @param beta Outcome-model coefficients (beta_1, beta_2).
#' @param alpha1 Mean of the truncated-normal event-time covariate X.
#' @param alpha2 Mean of the truncated-normal censoring time C (governs the censoring rate).
#' @param sigma Outcome standard deviation.
#' @param tau1,tau2 Standard deviations of X and C.
#' @return A list with x_data, y_data, c_data, w_data, delta_data, and the realized cens_rate.
#' @export
data_generating = function(k, n, d, beta, alpha1, alpha2,
                           sigma = 4, tau1 = 1, tau2 = 1){
  set.seed(k)
  # X|Z ~ N(alpha1, 1) truncated in [-1,1]
  x_data = truncnorm::rtruncnorm(n,
                      a = -1,
                      b = 1,
                      mean = alpha1,
                      sd = tau1)
  # Y|X ~ N(beta_0 + beta_1 X, 4)
  trunc_norm = truncnorm::rtruncnorm(n,
                          a = -3,
                          b = 3,
                          mean = 0,
                          sd = 1)
  y_data = cbind(1,x_data) %*% beta + sigma * trunc_norm
  # C|Y ~ N(alpha2_0 + alpha2_1 Y,1)
  c_data = truncnorm::rtruncnorm(n,
                      a = -1,
                      b = 1,
                      mean = alpha2,
                      sd = tau2)
  # generate W and Delta
  w_data = pmin(c_data, x_data)
  delta_data = as.numeric(x_data<=c_data)
  # calculate the censoring rate
  cens_rate = sum(c_data<=x_data)/n
  list(x_data = x_data,
       y_data = y_data,
       c_data = c_data,
       w_data = w_data,
       delta_data = delta_data,
       cens_rate = cens_rate)
}

# Define score function
#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_S_beta_f = function(beta, y, x, sigma = 4){
  res = y - (beta[1] + beta[2]*x)
  res * c(1,x) / sigma^2
}
sim_S_beta_f = Vectorize(sim_S_beta_f, vectorize.args = 'x')

#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_S_beta = function(beta, y, w, delta, alpha1_star,
                  sigma = 4, tau1 = 1){
  #E[I(X>W)S_beta_f(Y,X,Z) | W,Y,Z] / E[I(X>W) | W,Y,Z]
  v_x = 1/(beta[2]^2 / sigma^2 + 1/tau1^2)
  eta_x = v_x * (beta[2] * (y - beta[1]) / sigma^2
                 + alpha1_star / tau1^2)
  if(delta==0){
    if(w>1){
      sim_S_beta_f(beta, y, w, sigma) # Any value possible
    } else{ #Middle part: integration
      x_norm = (seq(w, 1, length.out = 20) - eta_x) / sqrt(v_x)
      by = x_norm[2] - x_norm[1]
      d = dnorm(x_norm)
      num = sim_S_beta_f(beta, y, x_norm*sqrt(v_x)+eta_x, sigma) %*% d
      denom = sum(d)
      num / denom
    }
  }else{
    sim_S_beta_f(beta, y, w, sigma)
  }
}

#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_gauss = function(tt, len = 3){
  grid = seq(-len, len, length.out = tt)
  d = dnorm(grid)
  list(x = grid,
       w = d / sum(d))
}

# b_xz function will be evaluated on grid x_a
#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_c0_xz_gauss_param12 = function(beta, x_a,
                               alpha1_star, alpha2_star,
                               sigma = 4, tau1 = 1, tau2 = 1,
                               tt = 20){
  cc = sim_gauss(tt)
  len = 20
  # temp2: E2[I(X<C) | Y,X,Z]S_beta(Y,X,Z) + E2[I(X>C)S_beta(Y,C,Z)| Y,X,Z]
  temp2 = function(x, y){
    v_x = 1/(beta[2]^2 / sigma^2 + 1/tau1^2)
    eta_x = v_x * (beta[2] * (y - beta[1]) / sigma^2
                   + alpha1_star / tau1^2)
    if(x < -1){
      sim_S_beta_f(beta, y, x, sigma)
    }else{
      c_grid = seq(-1, min(x, 1),
                   length.out = len)
      dens_c = truncnorm::dtruncnorm(c_grid,
                          a = -1,
                          b = 1,
                          mean = alpha2_star,
                          sd = tau2)
      Sbeta = sapply(c_grid, function(w)
        sim_S_beta(beta, y, w, 0, alpha1_star,
               sigma, tau1))
      Sbeta[,c(1,len)] = Sbeta[,c(1,len)] / 2
      by = c_grid[2] - c_grid[1]
      return(sim_S_beta_f(beta, y, x, sigma) *
               (1 - truncnorm::ptruncnorm(x,
                               a = -1,
                               b = 1,
                               mean = alpha2_star,
                               sd = tau2)) +
               Sbeta %*% dens_c * by)
    }
  }

  # # temp3: E[temp2(Y,X,Z)| X,Z]
  temp3 = function(x){
    sapply(sigma * cc$x + sum(c(1, x) * beta),
           function(y_norm)
             temp2(x, y_norm)) %*% cc$w
  }
  temp3 = Vectorize(temp3, vectorize.args = 'x')
  return(temp3(x_a))
}

#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_L_xz_gauss_param12 = function(beta, x_a,
                              alpha1_star, alpha2_star,
                              sigma = 4, tau1 = 1, tau2 = 1,
                              tt = 20){
  cc = sim_gauss(tt)
  #temp[a(X,Z)] =  E[I(X>C)a(X,Z)|C,Y,Z] / E[I(X>C)|C,Y,Z]
  temp = function(x_a, c, y){
    v_x = 1/(beta[2]^2 / sigma^2 + 1/tau1^2)
    eta_x = v_x * (beta[2] * (y - beta[1]) / sigma^2
                   + alpha1_star / tau1^2)
    # p: discretization of f(x|y,z) on x_a
    p = truncnorm::dtruncnorm(x_a,
                   a = -1,
                   b = 1,
                   mean = eta_x,
                   sd = sqrt(v_x))
    num = (x_a > c) * p
    denom = sum(num)
    ifelse(is.nan(num/denom),
           0, num/denom)
  }
  # temp2[a(X,Z)] = E[I(X>C)temp(x_a,C,Y,Z)| X,Y,Z]
  temp2 = function(x_a, x, y){
    c_grid = seq(-1, 1, length = 20)
    dens = truncnorm::dtruncnorm(c_grid,
                      a = -1,
                      b = 1,
                      mean = alpha2_star,
                      sd = tau2)
    dens = dens / sum(dens)
    sapply(c_grid, function(c_norm){
      (x > c_norm) * temp(x_a, c_norm, y)}) %*% dens
  }
  # temp3[a(X,Z)] = E[temp2(X,Y,Z)| X,Z]
  temp3 = function(x_a, x){
    sapply(sigma * cc$x + sum(beta * c(1,x)), function(y_norm){
      temp2(x_a, x, y_norm)}) %*% cc$w
  }
  temp3 = Vectorize(temp3, vectorize.args = 'x')
  # temp4 = E[I(X<C) | Y,X,Z]a(X,Z)
  temp4 = function(x, y){
    c_grid = seq(-1, 1, length = 20)
    dens = truncnorm::dtruncnorm(c_grid,
                      a = -1,
                      b = 1,
                      mean = alpha2_star,
                      sd = tau2)
    dens = dens / sum(dens)
    sum((x <= c_grid) * dens)
    # 1-truncnorm::ptruncnorm(x,
    #              a = -1,
    #              b = 1,
    #              mean = sum(c(1,y) * alpha2_star),
    #              sd = tau2)
  }
  temp4 = Vectorize(temp4, vectorize.args = 'y')
  # temp5 = E[temp4(X,Y,Z) | X,Z]a(X,Z)
  temp5 = function(x){
    sum(temp4(x, sigma * cc$x + sum(beta * c(1,x))) * cc$w)
  }
  temp5 = Vectorize(temp5, vectorize.args = 'x')
  res = diag(temp5(x_a)) + t(temp3(x_a, x_a))
  res
}

# Solution of integral equation
#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_a_gauss_param12 = function(beta, x_a,
                           alpha1_star, alpha2_star,
                           sigma = 4, tau1 = 1, tau2 = 1,
                           tt = 20){
  MASS::ginv(sim_L_xz_gauss_param12(beta, x_a,
                                alpha1_star, alpha2_star,
                                sigma, tau1, tau2,
                                tt)) %*%
    t(sim_c0_xz_gauss_param12(beta, x_a,
                          alpha1_star, alpha2_star,
                          sigma, tau1, tau2,
                          tt))
}


#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_S_eff_gauss_param12 = function(beta, y, w, delta,
                               x_a, a0,
                               alpha1_star,
                               sigma = 4, tau1 = 1, tau2 = 1,
                               tt = 20){
  len_beta = length(beta)
  if(delta==0){
    #temp[a(X,Z)] =  E[I(X>C)a(X,Z)|C,Y,Z] / E[I(X>C)|C,Y,Z]
    temp = function(x_a, c, y){
      v_x = 1/(beta[2]^2 / sigma^2 + 1/tau1^2)
      eta_x = v_x * (beta[2] * (y - beta[1]) / sigma^2
                     + alpha1_star / tau1^2)
      # p: discretization of f(x|y,z) on x_a
      p = truncnorm::dtruncnorm(x_a,
                     a = -1,
                     b = 1,
                     mean = eta_x,
                     sd = sqrt(v_x))
      num = (x_a > c) * p
      denom = sum(num)
      ifelse(is.nan(num/denom),
             0, num/denom)
    }
    #S_beta
    # sbeta = S_beta(beta, y, w, 0, z, alpha1_star, sigma, tau1)
    sbeta = sim_S_beta(beta, y, w, 0,
                   alpha1_star,
                   sigma, tau1)
    #S_eff = S_beta - a0
    return(sbeta - t(a0) %*% temp(x_a, w, y))
  } else{ # If delta=1, linearize from the grid
    m = length(x_a)
    a0_w = vector(length = len_beta)
    for(j in 1:len_beta){
      a0_w[j] = approx(x_a, a0[,j], w, rule = 2)$y #linear interpolation
    }
    sbeta = sim_S_beta_f(beta, y, w, sigma)
    #S_eff = S_beta - a0
    return(sbeta - a0_w)
  }
}

#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_pe_gauss_param12 = function(beta, y_data, w_data, delta_data,
                            x_a,
                            alpha1_star, alpha2_star,
                            sigma = 4, tau1 = 1, tau2 = 1,
                            tt = 20){
  len_beta = length(beta)
  # Define a0
  a0 = sim_a_gauss_param12(beta, x_a,
                       alpha1_star, alpha2_star,
                       sigma, tau1, tau2, tt)
  val = rep(0, len_beta)
  n = length(y_data)
  for(i in 1:n){
    val = val + sim_S_eff_gauss_param12(beta, y_data[i], w_data[i], delta_data[i],
                                    x_a, a0,
                                    alpha1_star,
                                    sigma, tau1, tau2, tt)
  }
  return(val)
}


#' Profile MLE of the eta1 mean given the current beta.
#'
#' @keywords internal
sim_find_alpha1_MLE = function(beta, y_data, w_data, delta_data,
                           sigma = 4, tau1 = 1){
  log_llhd = function(alpha1_star, y, w, delta){
    v_x = 1/(beta[2]^2 / sigma^2 + 1/tau1^2)
    eta_x = v_x * (beta[2] * (y - beta[1]) / sigma^2
                   + alpha1_star / tau1^2)
    if(delta == 1){ #log(eta1(w|z))
      return(log(truncnorm::dtruncnorm(w,
                            a = -1,
                            b = 1,
                            mean = eta_x,
                            sd = sqrt(v_x))))
    }else{ #log(P(X > w|y,z) * f(y|z))
      v_x = 1/(beta[2]^2 / sigma^2 + 1/tau1^2)
      eta_x = v_x * (beta[2] * (y - beta[1]) / sigma^2
                     + alpha1_star / tau1^2)
      # Y|Z ~ N(beta_0 + beta_1 E(X|Z) + beta_2^T Z, 1 + beta_1^2)
      return(log(pnorm(1,
                       mean = eta_x,
                       sd = sqrt(v_x))
                 - pnorm(w,
                         mean = eta_x,
                         sd = sqrt(v_x))) +
               dnorm(y,
                     mean = sum(c(1, alpha1_star) * beta),
                     sd = sqrt(sigma^2 + beta[2]^2 * tau1^2),
                     log = T) -
               log(pnorm(1,
                         mean = alpha1_star,
                         sd = tau1) -
                     pnorm(-1,
                           mean = alpha1_star,
                           sd = tau1))
      )
    }
  }
  log_llhd_sum = function(alpha1){
    n = length(y_data)
    llhd = 0
    for(i in 1:n){
      llhd = llhd +
        log_llhd(alpha1, y_data[i], w_data[i], delta_data[i])
    }
    return(-llhd)
  }
  optimize(log_llhd_sum, interval = c(-5,5))$minimum
}


# Find MLE of alpha2
#' MLE of the eta2 (censoring) mean.
#'
#' @keywords internal
sim_find_alpha2_MLE = function(y_data, w_data, delta_data,
                           sigma = 4, tau2 = 1){
  log_llhd = function(alpha2_star, y, w, delta){
    if(delta == 1){ #log(P(C > w | y,z))
      return(log(1-truncnorm::ptruncnorm(w,
                              a = -1,
                              b = 1,
                              mean = sum(c(1, y) * alpha2_star),
                              sd = tau2)))
    }else{ #log(eta2(w|y,z))
      # Y|Z ~ N(beta_0 + beta_1 E(X|Z) + beta_2^T Z, 1 + beta_1^2)
      return(log(truncnorm::dtruncnorm(w,
                            a = -1,
                            b = 1,
                            mean = sum(c(1, y) * alpha2_star),
                            sd = tau2)))
    }
  }
  log_llhd_sum = function(alpha2){
    n = length(y_data)
    llhd = 0
    for(i in 1:n){
      llhd = llhd + log_llhd(alpha2, y_data[i], w_data[i], delta_data[i])
    }
    return(-llhd)
  }
  optim(rep(0,d+1), log_llhd_sum)$par
}


#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
pe_gauss_alpha1_MLE = function(beta, y_data, w_data, delta_data,
                               x_a,
                               alpha2_star,
                               sigma = 4, tau1 = 1, tau2 = 1,
                               tt = 20){
  len_beta = length(beta)
  alpha1_star = sim_find_alpha1_MLE(beta, y_data, w_data, delta_data,
                                sigma, tau1)
  # Define a0
  a0 = sim_a_gauss_param12(beta, x_a,
                       alpha1_star, alpha2_star,
                       sigma, tau1, tau2, tt)
  val = rep(0, len_beta)
  n = length(y_data)
  for(i in 1:n){
    val = val + sim_S_eff_gauss_param12(beta, y_data[i], w_data[i], delta_data[i],
                                    x_a, a0,
                                    alpha1_star,
                                    sigma, tau1, tau2, tt)
  }
  return(val)
}


#' Estimate beta by solving the efficient estimating equation (SPARCC).
#'
#' @param k Integer seed / replicate index.
#' @param n,d Sample size and covariate dimension.
#' @param beta True coefficient vector (used to generate data and as the
#'   starting value for the solver).
#' @param alpha1,alpha2 True nuisance means used to \emph{generate} the data.
#' @param alpha1_star,alpha2_star Nuisance means used in the \emph{estimating}
#'   equation. Set them equal to \code{alpha1}/\code{alpha2} for correctly
#'   specified nuisance models, or to perturbed values to obtain the misspecified
#'   (mis1/mis2) beta estimates.
#' @param sigma,tau1,tau2 Model standard deviations.
#' @param tt Integration grid size forwarded to the influence-function solver.
#' @param m Resolution of the covariate grid \code{x_a}.
#' @return The estimated coefficient vector for replicate \code{k}.
#' @export
get_beta_param_param12 = function(k, n, d,
                                  beta, alpha1, alpha2,
                                  alpha1_star = alpha1, alpha2_star = alpha2,
                                  sigma = 4, tau1 = 1, tau2 = 1, tt = 20, m = 20){
  data_k = data_generating(k, n, d, beta, alpha1, alpha2, sigma, tau1, tau2)
  y_data = data_k$y_data
  w_data = data_k$w_data
  delta_data = data_k$delta_data

  x_a = seq(-1, 1, length.out = m)
  #equation solve
  set.seed(k)
  result = nleqslv::nleqslv(beta + rnorm(d+1) * 0.1, sim_pe_gauss_param12,
                            y_data = y_data,
                            w_data = w_data,
                            delta_data = delta_data,
                            x_a=x_a,
                            alpha1_star = alpha1_star,
                            alpha2_star = alpha2_star,
                            sigma = sigma,
                            tau1 = tau1,
                            tau2 = tau2,
                            tt = tt)
  print(result$message)
  return(result$x)
}
# get_beta_param_param12(1, n, d,
#                        beta, alpha1, alpha2,
#                        sigma = 4, tau1 = 1, tau2 = 1, tt = 20, m = 20)


#' Interval center m0: the outcome mean beta_1 + beta_2 x.
#'
sim_m0 = function(x, beta){
  beta[1] + beta[2]*x
}
#' Interval center m1: mean outcome given the observed W and Delta under eta1.
#'
sim_m1 = function(w, delta, beta,
              alpha1_star, tau1 = 1){
  #E[I(X>W)m0(X,Z) | W,Z] / E[I(X>W) | W,Z]
  v_x = tau1^2 #1/(beta[2]^2 / sigma^2 + 1/tau1^2)
  eta_x = alpha1_star
  if(delta==0){
    if(w>=1){
      sim_m0(w, beta) # Any value possible
    } else{ #Middle part: integration
      # x_norm = (seq(w, 1, length.out = 50) - alpha1_star) / tau1
      # d = dnorm(x_norm)
      # num = m0(x_norm*tau1 + alpha1_star, beta) %*% d
      # denom = sum(d)
      # num / denom

      #x_exp is the truncated mean
      x_exp = alpha1_star +
        tau1 * (dnorm((w - alpha1_star) / tau1) - dnorm((1 - alpha1_star) / tau1)) /
        (pnorm(1, mean = eta_x, sd = sqrt(v_x)) -
           pnorm(w, mean = eta_x, sd = sqrt(v_x)))
      beta[1] + beta[2] * x_exp
    }
  }else{
    sim_m0(w, beta)
  }
}
#' Residual r1 = |Y - m1|.
#'
sim_r1 = function(y, w, delta, beta,
              alpha1_star, tau1 = 1){
  return(abs(y - sim_m1(w, delta, beta,
                    alpha1_star, tau1)))
}
sim_r1 = Vectorize(sim_r1, vectorize.args = 'w')

#' Residual r2 = |Y - m0(W)|.
#'
sim_r2 = function(y, w, delta, beta,
              alpha1_star, tau1 = 1){
  return(abs(y - sim_m0(w, beta)))
}
sim_r2 = Vectorize(sim_r2, vectorize.args = 'w')


# ---- Vectorized centers and residuals ---------------------------------------
# Faster, vector-in/vector-out equivalents of m0/m1/r1/r2, used by the
# conformal routines that evaluate residuals over large test sets.

#' Vectorized outcome-mean center m0
#'
#' @param w Numeric vector of observed times.
#' @param beta Length-2 coefficient vector (intercept, slope).
#' @return Numeric vector \code{beta[1] + beta[2] * w}.
sim_m0_vec <- function(w, beta) {
  beta[1] + beta[2] * w
}

#' Vectorized center m1
#'
#' Vectorized form of \code{m1}: for \code{delta == 1} or \code{w >= 1} it
#' returns \code{m0(w)}; for \code{delta == 0} with \code{w < 1} it returns
#' \code{m0()} applied to the truncated-normal mean of X on \code{(w, 1)} under
#' \code{N(alpha1_star, tau1^2)}.
#'
#' @param w,delta Numeric/integer vectors of observed times and event flags.
#' @param beta Length-2 coefficient vector.
#' @param alpha1_star Center (eta1 mean) used inside the residual.
#' @param tau1 Standard deviation of X.
#' @return Numeric vector of centers.
sim_m1_vec <- function(w, delta, beta, alpha1_star, tau1 = 1) {
  w     <- as.numeric(w)
  delta <- as.integer(delta)
  n     <- length(w)
  out   <- numeric(n)

  # Masks: A = observed; B = censored past the boundary; C = censored inside.
  A <- (delta == 1L)
  B <- (delta == 0L) & (w >= 1)
  C <- (delta == 0L) & (w < 1)

  if (any(A)) out[A] <- sim_m0_vec(w[A], beta)
  if (any(B)) out[B] <- sim_m0_vec(w[B], beta)

  # Case C: truncated-normal mean on (w, 1).
  if (any(C)) {
    a  <- alpha1_star
    s  <- tau1
    z1 <- (1 - a) / s

    zw    <- (w[C] - a) / s
    denom <- pnorm(z1) - pnorm(zw)
    numer <- dnorm(zw) - dnorm(z1)

    tiny  <- 1e-12
    good  <- denom > tiny
    x_exp <- numeric(sum(C))
    x_exp[good]  <- a + s * (numer[good] / denom[good])
    x_exp[!good] <- w[C][!good]   # fallback when the interval has ~zero mass

    out[C] <- beta[1] + beta[2] * x_exp
  }

  out
}

#' Vectorized residual r1 = |Y - m1|
#'
#' @param y,w,delta Numeric vectors.
#' @param beta Length-2 coefficient vector.
#' @param alpha1_star Center used inside the residual.
#' @param tau1 Standard deviation of X.
#' @return Numeric vector of residuals.
sim_r1_vec <- function(y, w, delta, beta, alpha1_star, tau1 = 1) {
  abs(y - sim_m1_vec(w, delta, beta, alpha1_star, tau1))
}

#' Vectorized residual r2 = |Y - m0(W)|
#'
#' @inheritParams r1_vec
#' @return Numeric vector of residuals.
sim_r2_vec <- function(y, w, delta, beta, alpha1_star, tau1 = 1) {
  # alpha1_star, tau1 unused; kept for a common signature with r1_vec.
  abs(y - sim_m0_vec(w, beta))
}


#' Monte Carlo approximation of the true half-length for residual r1.
#'
#' @export
get_zeta_r1 = function(alpha, alpha1, alpha2,
                       beta = c(0, 3), sigma = 4,
                       tau1 = 1, tau2 = 1, d = 1,
                       N = 1000000){
  data_mc = data_generating(9999999, N, d, beta, alpha1, alpha2,
                            sigma = sigma, tau1 = tau1, tau2 = tau2)
  y_data = data_mc$y_data
  w_data = data_mc$w_data
  delta_data = data_mc$delta_data
  return(quantile(sapply(1:N, function(i){
    sim_r1(y_data[i], w_data[i], delta_data[i],
       beta, alpha1, tau1)
  }), 1-alpha))
}

#' Monte Carlo approximation of the true half-length for residual r2.
#'
#' @export
get_zeta_r2 = function(alpha, alpha1, alpha2,
                       beta = c(0, 3), sigma = 4,
                       tau1 = 1, tau2 = 1, d = 1,
                       N = 1000000){
  data_mc = data_generating(9999999, N, d, beta, alpha1, alpha2,
                            sigma = sigma, tau1 = tau1, tau2 = tau2)
  y_data = data_mc$y_data
  w_data = data_mc$w_data
  delta_data = data_mc$delta_data
  return(quantile(sapply(1:N, function(i){
    sim_r2(y_data[i], w_data[i], delta_data[i],
       beta, alpha1, tau1)
  }), 1-alpha))
}
# zeta_true_r1 = get_zeta_r1(0.2, alpha1, alpha2 = 2)
# zeta_true_r2 = get_zeta_r2(0.2, alpha1, alpha2 = 2)



# c1_xz function will be evaluated on grid x_a
#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_c1_xz_gauss_param12 = function(zeta, alpha, r,
                               beta, x_a,
                               alpha1_star, alpha2_star, alpha1_star_r,
                               sigma = 4, tau1 = 1, tau2 = 1,
                               tt = 20, tt2 = 100){
  cc = sim_gauss(tt2)
  # temp2: E2[I(X<C) | Y,X,Z]I(Y,X,Z) + E2[I(X>C)I(Y,C,Z)| Y,X,Z]
  temp2 = function(x, y){
    if(x < -1){
      (r(y, x, 1, beta, alpha1_star_r, tau1) <= zeta)
    }else{
      c_grid = seq(-1, 1,
                   length.out = tt)
      dens_c = truncnorm::dtruncnorm(c_grid,
                          a = -1,
                          b = 1,
                          mean = alpha2_star,
                          sd = tau2)
      w_grid = pmin(x, c_grid)
      delta_grid = as.numeric(x <= c_grid)
      RHS = sapply(1:tt, function(i){
        (r(y, w_grid[i], delta_grid[i], beta, alpha1_star_r, tau1) <= zeta)})
      return(RHS %*% dens_c / sum(dens_c))
    }
  }
  # # temp3: E[temp2(Y,X,Z)| X,Z]
  temp3 = function(x){
    sapply(sigma * cc$x + sum(c(1, x) * beta),
           function(y_norm)
             temp2(x, y_norm)) %*% cc$w
  }
  temp3 = Vectorize(temp3, vectorize.args = 'x')
  return(temp3(x_a) - (1-alpha))
}


#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
make_seq = function(n){
  seq(1:n) * 2 / (2*n+1)
}
# c1_xz function will be evaluated on grid x_a
#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
c1_xz_gauss_param12_2 = function(zeta, alpha, r,
                                 beta, x_a,
                                 alpha1_star, alpha2_star, alpha1_star_r,
                                 sigma = 4, tau1 = 1, tau2 = 1,
                                 tt = 20, tt2 = 20){

  # x_vec = truncnorm::qtruncnorm(make_seq(tt), a = -1,
  #                    b = 1,
  #                    mean = alpha1_star,
  #                    sd = tau1)
  c_vec = truncnorm::qtruncnorm(make_seq(tt), a = -1,
                     b = 1,
                     mean = alpha2_star,
                     sd = tau2)
  return(sapply(x_a, function(x){
    y_vec = sigma * truncnorm::qtruncnorm(make_seq(tt2), a = -3,
                               b = 3,
                               mean = 0,
                               sd = 1) + sum(c(1, x) * beta)
    mean(sapply(1:tt, function(i){
      mean(sapply(1:tt2, function(j){
        (r(y_vec[j], min(x,c_vec[i]), (x<=c_vec[i]),
           beta, alpha1_star, tau1)
         <= zeta)}))})) - (1-alpha)}))
}


#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_M1_xz_gauss_param12 = function(beta, x_a,
                               alpha1_star, alpha2_star,
                               sigma = 4, tau1 = 1, tau2 = 1,
                               tt = 20, tt2 = 20){
  cc = sim_gauss(tt2)
  #temp[a(X,Z)] =  E[I(X>C)a(X,Z)|C,Y,Z] / E[I(X>C)|C,Y,Z]
  temp = function(x_a, c, y){
    v_x = 1/(beta[2]^2 / sigma^2 + 1/tau1^2)
    eta_x = v_x * (beta[2] * (y - beta[1]) / sigma^2
                   + alpha1_star / tau1^2)
    # p: discretization of f(x|y,z) on x_a
    p = truncnorm::dtruncnorm(x_a,
                   a = -1,
                   b = 1,
                   mean = eta_x,
                   sd = sqrt(v_x))
    num = (x_a > c) * p
    denom = sum(num)
    ifelse(is.nan(num/denom),
           0, num/denom)
  }
  # temp2[a(X,Z)] = E[I(X>C)temp(x_a,C,Y,Z)| X,Y,Z]
  temp2 = function(x_a, x, y){
    c_grid = seq(-1, 1, length = tt)
    dens = truncnorm::dtruncnorm(c_grid,
                      a = -1,
                      b = 1,
                      mean = alpha2_star,
                      sd = tau2)
    dens = dens / sum(dens)
    sapply(c_grid, function(c_norm){
      (x > c_norm) * temp(x_a, c_norm, y)}) %*% dens
  }
  # temp3[a(X,Z)] = E[temp2(X,Y,Z)| X,Z]
  temp3 = function(x_a, x){
    sapply(sigma * cc$x + sum(beta * c(1,x)), function(y_norm){
      temp2(x_a, x, y_norm)}) %*% cc$w
  }
  temp3 = Vectorize(temp3, vectorize.args = 'x')
  # temp4 = E[I(X<C) | Y,X,Z]a(X,Z)
  temp4 = function(x, y){
    c_grid = seq(-1, 1, length = tt)
    dens = truncnorm::dtruncnorm(c_grid,
                      a = -1,
                      b = 1,
                      mean = alpha2_star,
                      sd = tau2)
    dens = dens / sum(dens)
    sum((x <= c_grid) * dens)
    # 1-truncnorm::ptruncnorm(x,
    #              a = -1,
    #              b = 1,
    #              mean = sum(c(1,y) * alpha2_star),
    #              sd = tau2)
  }
  temp4 = Vectorize(temp4, vectorize.args = 'y')
  # temp5 = E[temp4(X,Y,Z) | X,Z]a(X,Z)
  temp5 = function(x){
    sum(temp4(x, sigma * cc$x + sum(beta * c(1,x))) * cc$w)
  }
  temp5 = Vectorize(temp5, vectorize.args = 'x')
  res = diag(temp5(x_a)) + t(temp3(x_a, x_a))
  res
}


# Solution of integral equation
#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_b1_gauss_param12 = function(zeta, alpha, r,
                            beta, x_a,
                            alpha1_star, alpha2_star, alpha1_star_r,
                            sigma = 4, tau1 = 1, tau2 = 1,
                            tt = 20, tt2 = 100){
  MASS::ginv(sim_M1_xz_gauss_param12(beta, x_a,
                                 alpha1_star, alpha2_star,
                                 sigma, tau1, tau2,
                                 tt, tt2)) %*%
    (sim_c1_xz_gauss_param12(zeta, alpha, r,
                         beta, x_a,
                         alpha1_star, alpha2_star, alpha1_star_r,
                         sigma, tau1, tau2,
                         tt, tt2))
}


# c2_xz function will be evaluated on grid x_a
#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_c2_xz_gauss_param12 = function(zeta, alpha, r,
                               beta, x_a,
                               alpha1_star, alpha2_star, alpha1_star_r,
                               sigma = 4, tau1 = 1, tau2 = 1,
                               tt = 20, tt2 = 100){

  cc = sim_gauss(tt2)
  # temp2: E2[I(X<C) | Y,X,Z]S_beta(Y,X,Z) + E2[I(X>C)S_beta(Y,C,Z)| Y,X,Z]
  temp2 = function(c, x){
    return(sapply(sigma * cc$x + sum(c(1, x) * beta),
                  function(y){
                    (r(y, min(x,c), (x<=c),
                       beta, alpha1_star_r, tau1) <= zeta)}) %*% cc$w)

  }
  # # temp3: E[temp2(C,X,Z)| C,Z]
  temp3 = function(c){
    x_grid = seq(-1, 1,
                 length.out = tt)
    dens_x = truncnorm::dtruncnorm(x_grid,
                        a = -1,
                        b = 1,
                        mean = alpha1_star,
                        sd = tau1)
    RHS = sapply(1:tt, function(i){
      temp2(c, x_grid[i])})
    return(RHS %*% dens_x / sum(dens_x))

  }
  temp3 = Vectorize(temp3, vectorize.args = 'c')
  return(temp3(x_a) - (1-alpha))
}
# c2_xz_gauss_param12(zeta, alpha, r,
#                                beta, x_a,
#                                alpha1_star, alpha2_star,
#                                sigma = 4, tau1 = 1, tau2 = 1,
#                                tt = 20)


#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_M2_xz_gauss_param12 = function(beta, x_a,
                               alpha1_star, alpha2_star,
                               sigma = 4, tau1 = 1, tau2 = 1,
                               tt = 50, tt2 = 50){
  cc = sim_gauss(tt2)
  #temp[a(X,Z)] =  E[I(X>C)a(X,Z)|C,Y,Z] / E[I(X>C)|C,Y,Z]
  temp = function(x_a, x){
    # p: discretization of f(x|y,z) on x_a
    p = truncnorm::dtruncnorm(x_a,
                   a = -1,
                   b = 1,
                   mean = alpha2_star,
                   sd = tau2)
    num = (x_a >= x) * p
    denom = sum(num)
    ifelse(is.nan(num/denom),
           0, num/denom)
  }
  # temp2[a(X,Z)] = E[I(X>C)temp(x_a,C,Y,Z)| X,Y,Z]
  temp2 = function(x_a, c){
    x_grid = seq(-1, 1, length = tt)
    dens = truncnorm::dtruncnorm(x_grid,
                      a = -1,
                      b = 1,
                      mean = alpha1_star,
                      sd = tau1)
    dens = dens / sum(dens)
    sapply(x_grid, function(x_norm){
      (x_norm <= c) * temp(x_a, x_norm)}) %*% dens
  }
  temp2 = Vectorize(temp2, vectorize.args = 'c')
  # temp4 = E[I(X>C) | C,Z]a(X,Z)
  temp3 = function(c){
    x_grid = seq(-1, 1, length = tt)
    dens = truncnorm::dtruncnorm(x_grid,
                      a = -1,
                      b = 1,
                      mean = alpha1_star,
                      sd = tau1)
    dens = dens / sum(dens)
    sum((x_grid > c) * dens)
    # 1-truncnorm::ptruncnorm(x,
    #              a = -1,
    #              b = 1,
    #              mean = sum(c(1,y) * alpha2_star),
    #              sd = tau2)
  }
  temp3 = Vectorize(temp3, vectorize.args = 'c')
  res = diag(temp3(x_a)) + t(temp2(x_a, x_a))
  res
}


# Solution of integral equation
#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_b2_gauss_param12 = function(zeta, alpha, r,
                            beta, x_a,
                            alpha1_star, alpha2_star, alpha1_star_r,
                            sigma = 4, tau1 = 1, tau2 = 1,
                            tt = 50, tt2 = 50){
  MASS::ginv(sim_M2_xz_gauss_param12(beta, x_a,
                                 alpha1_star, alpha2_star,
                                 sigma, tau1, tau2,
                                 tt, tt2)) %*%
    (sim_c2_xz_gauss_param12(zeta, alpha, r,
                         beta, x_a,
                         alpha1_star, alpha2_star, alpha1_star_r,
                         sigma, tau1, tau2,
                         tt, tt2))
}


# c2_xz function will be evaluated on grid x_a
#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_b3_gauss_param12 = function(zeta, alpha, r,
                            beta, x_a,
                            alpha1_star, alpha2_star, alpha1_star_r,
                            sigma = 4, tau1 = 1, tau2 = 1,
                            tt = 50, tt2 = 50, tt3 = 50){
  cc = sim_gauss(tt2)
  # temp2: E2[I(X<C) | Y,X,Z]S_beta(Y,X,Z) + E2[I(X>C)S_beta(Y,C,Z)| Y,X,Z]
  temp2 = function(c, x){
    return(sapply(sigma * cc$x + sum(c(1, x) * beta),
                  function(y){
                    (r(y, min(x,c), (x<=c),
                       beta, alpha1_star_r, tau1) <= zeta)}) %*% cc$w)

  }
  # # temp3: E[temp2(C,X,Z)| C,Z]
  temp3 = function(c){
    x_grid = seq(-1, 1,
                 length.out = tt)
    dens_x = truncnorm::dtruncnorm(x_grid,
                        a = -1,
                        b = 1,
                        mean = alpha1_star,
                        sd = tau1)
    RHS = sapply(1:tt, function(i) {temp2(c, x_grid[i])})
    return(RHS %*% dens_x / sum(dens_x))

  }
  temp4 = function(){
    c_grid = seq(-1, 1, length = tt3)
    dens = truncnorm::dtruncnorm(c_grid,
                      a = -1,
                      b = 1,
                      mean = alpha2_star,
                      sd = tau2)
    dens = dens / sum(dens)
    sapply(c_grid, function(c) {temp3(c)}) %*% dens
    # 1-truncnorm::ptruncnorm(x,
    #              a = -1,
    #              b = 1,
    #              mean = sum(c(1,y) * alpha2_star),
    #              sd = tau2)
  }
  return(temp4() - (1-alpha))
}


#' Internal numerical helper for the efficient influence-function / integral-equation solvers.
#'
#' @keywords internal
sim_trans_phi_eff_gauss_param12 = function(zeta, alpha, r,
                                       beta, y, w, delta,
                                       x_a, b1, b2, b3,
                                       alpha1_star, alpha2_star,
                                       sigma = 4, tau1 = 1, tau2 = 1){
  if(delta == 1){
    #temp[a(X,Z)] =  E[I(C>=X)a(C,Z)|X,Z] / E[I(C>=X)|X,Z]
    temp = function(x_a, x, y){
      # p: discretization of f(c|z) on x_a
      p = truncnorm::dtruncnorm(x_a,
                     a = -1,
                     b = 1,
                     mean = alpha2_star,
                     sd = tau2)
      num = (x_a >= x) * p
      denom = sum(num)
      ifelse(is.nan(num/denom),
             0, num/denom)
    }
    b1_w = approx(x_a, b1, w, rule = 2)$y #linear interpolation
    return(b1_w + t(b2) %*% temp(x_a, w, y) - b3)
  }else{
    #temp[a(X,Z)] =  E[I(X>C)a(X,Z)|C,Y,Z] / E[I(X>C)|C,Y,Z]
    temp = function(x_a, c, y){
      v_x = 1/(beta[2]^2 / sigma^2 + 1/tau1^2)
      eta_x = v_x * (beta[2] * (y - beta[1]) / sigma^2
                     + alpha1_star / tau1^2)
      # p: discretization of f(x|y,z) on x_a
      p = truncnorm::dtruncnorm(x_a,
                     a = -1,
                     b = 1,
                     mean = eta_x,
                     sd = sqrt(v_x))
      num = (x_a > c) * p
      denom = sum(num)
      ifelse(is.nan(num/denom),
             0, num/denom)
    }
    b2_w = approx(x_a, b2, w, rule = 2)$y #linear interpolation
    return(t(b1) %*% temp(x_a, w, y) + b2_w - b3)
  }
}


#' Locate the zeta root of the estimating function by sign-change interpolation/extrapolation.
#'
#' @keywords internal
sim_find_sol = function(x_vec, y_vec){
  #find the root
  sign_vec = sign(y_vec)
  sign_change = which(diff(sign_vec) != 0)
  if(length(sign_change) == 0){
    #all positive or all negative: extrapolate
    if(all(y_vec > 0)){
      idx = which.min(y_vec)
      x1 = x_vec[idx]; x2 = x_vec[idx+1]
      y1 = y_vec[idx]; y2 = y_vec[idx+1]
      #linear extrapolation
      return(x1 - y1 * (x2 - x1) / (y2 - y1))
    }else{
      idx = which.max(y_vec)
      x1 = x_vec[idx-1]; x2 = x_vec[idx]
      y1 = y_vec[idx-1]; y2 = y_vec[idx]
      #linear extrapolation
      return(x2 - y2 * (x2 - x1) / (y2 - y1))
    }
  }else{
    idx = sign_change[which.min(abs(y_vec[sign_change]))]
    x1 = x_vec[idx]; x2 = x_vec[idx+1]
    y1 = y_vec[idx]; y2 = y_vec[idx+1]
    #linear interpolation
    return(x1 - y1 * (x2 - x1) / (y2 - y1))
  }
}
#' Bilinear interpolation of a value array over the two-dimensional beta grid.
#'
#' @keywords internal
interpolate = function(x_array, y_array, x){
  x_vec1 = x_array[1,]; x_vec2 = x_array[2,]
  x1 = x[1]; x2 = x[2]
  #interpolate y_array on (x_vec1, x_vec2)
  n1 = length(x_vec1)
  n2 = length(x_vec2)
  idx1 = max(which(x_vec1 <= x1))
  idx2 = max(which(x_vec2 <= x2))
  if(x1 < min(x_vec1)){
    idx1 = 1
  }
  if(x2 < min(x_vec2)){
    idx2 = 1
  }
  if(idx1 == n1){
    idx1 = n1 - 1
  }
  if(idx2 == n2){
    idx2 = n2 - 1
  }
  x11 = x_vec1[idx1]; x12 = x_vec1[idx1+1]
  x21 = x_vec2[idx2]; x22 = x_vec2[idx2+1]
  y11 = y_array[idx1, idx2]; y12 = y_array[idx1, idx2+1]
  y21 = y_array[idx1+1, idx2]; y22 = y_array[idx1+1, idx2+1]
  y1 = y11 + (y12 - y11) * (x2 - x21) / (x22 - x21)
  y2 = y21 + (y22 - y21) * (x2 - x21) / (x22 - x21)
  return(y1 + (y2 - y1) * (x1 - x11) / (x12 - x11))
}

#' Solve for the estimated half-length zeta by interpolating b1/b2/b3 over the beta grid.
#'
#' @return A list with the (zeta, objective) grid and the interpolated root `sol`.
#' @export
get_zeta_param_param12_int = function(k, n, d,
                                      alpha, r,
                                      beta, alpha1, alpha2,
                                      b1_array, b2_array, b3_array,
                                      beta_array,
                                      beta_temp,
                                      zeta_seq,
                                      alpha1_star,
                                      alpha2_star, 
                                      sigma = 4, tau1 = 1, tau2 = 1, m = 20){
  data_k = data_generating(k, n, d, beta, alpha1, alpha2, sigma, tau1, tau2)
  y_data = data_k$y_data
  w_data = data_k$w_data
  delta_data = data_k$delta_data
  x_a = seq(-1, 1, length.out = m)
  #interpolation
  func_vals = rep(0, length(zeta_seq))
  for(zeta_idx in 1:length(zeta_seq)){
    zeta = zeta_seq[zeta_idx]
    b1 = rep(0, m)
    b2 = rep(0, m)
    for(i in 1:length(b1_array[zeta_idx, , 1, 1])){
      b1[i] = interpolate(beta_array, b1_array[zeta_idx,i,,], beta_temp)
      b2[i] = interpolate(beta_array, b2_array[zeta_idx,i,,], beta_temp)
    }
    b3 = interpolate(beta_array, b3_array[zeta_idx,,], beta_temp)
    #Find the sum
    val = 0
    n = length(y_data)
    for(i in 1:n){
      val = val + sim_trans_phi_eff_gauss_param12(zeta, alpha, r,
                                              beta_temp, y_data[i], w_data[i], delta_data[i],
                                              x_a, b1, b2, b3,
                                              alpha1_star,
                                              alpha2_star,
                                              sigma, tau1, tau2)
    }
    func_vals[zeta_idx] = val
  }
  return(list(vec = rbind(zeta_seq, func_vals),
              sol = sim_find_sol(zeta_seq, func_vals)))
}

#' Build the search grid of zeta values for a given residual.
#'
#' @param alpha1_star_r Center used inside the residual; set to a misspecified value to reproduce r1*.
sim_get_zeta_seq_r = function(alpha, alpha1, alpha2, r,
                          alpha1_star_r = alpha1,
                          beta = c(0, 3), sigma = 4,
                          tau1 = 1, tau2 = 1, d = 1,
                          length = 5, N = 1000, reps = 10){
  zeta_r_vals = rep(0, reps)
  for(k in 1:reps){
    data_mc = data_generating(k, N, d, beta, alpha1, alpha2,
                              sigma = sigma, tau1 = tau1, tau2 = tau2)
    y_data = data_mc$y_data
    w_data = data_mc$w_data
    delta_data = data_mc$delta_data
    zeta_r_vals[k] = quantile(sapply(1:N, function(i){
      r(y_data[i], w_data[i], delta_data[i],
        beta, alpha1_star_r, tau1)
    }), 1-alpha)
  }
  return(seq(mean(zeta_r_vals) - sd(zeta_r_vals),
             mean(zeta_r_vals) + sd(zeta_r_vals),
             length.out = length))
}



################################
######### True beta #######################
################################

#' As get_zeta_param_param12_int but evaluated at the true beta (no interpolation).
#'
#' @export
get_zeta_param_param12_int_true = function(k, n, d,
                                           alpha, r,
                                           beta, alpha1, alpha2,
                                           b1_array, b2_array, b3_array,
                                           zeta_seq,
                                           alpha1_star,
                                           alpha2_star,
                                           sigma = 4, tau1 = 1, tau2 = 1, m = 20){
  data_k = data_generating(k, n, d, beta, alpha1, alpha2, sigma, tau1, tau2)
  y_data = data_k$y_data
  w_data = data_k$w_data
  delta_data = data_k$delta_data
  x_a = seq(-1, 1, length.out = m)
  #interpolation
  func_vals = rep(0, length(zeta_seq))
  for(zeta_idx in 1:length(zeta_seq)){
    zeta = zeta_seq[zeta_idx]
    b1 = rep(0, m)
    b2 = rep(0, m)
    for(i in 1:length(b1_array[zeta_idx, , 1, 1])){
      b1[i] = b1_array[zeta_idx,i,,]
      b2[i] = b2_array[zeta_idx,i,,]
    }
    b3 = b3_array[zeta_idx,,]
    #Find the sum
    val = 0
    n = length(y_data)
    for(i in 1:n){
      val = val + sim_trans_phi_eff_gauss_param12(zeta, alpha, r,
                                              beta, y_data[i], w_data[i], delta_data[i],
                                              x_a, b1, b2, b3,
                                              alpha1_star,
                                              alpha2_star,
                                              sigma, tau1, tau2)
    }
    func_vals[zeta_idx] = val
  }
  return(list(vec = rbind(zeta_seq, func_vals),
              sol = sim_find_sol(zeta_seq, func_vals)))
}



# =============================================================================
#  Shared scenario constants
#  ---------------------------------------------------------------------------
#  Used across the beta, semiparametric, and conformal simulation code.
# =============================================================================

# Censoring level -> true alpha2 (gamma). Governs the censoring rate:
# low ~ 20-30%, lowmid ~ 30-40%, mid ~ 45-55%, highmid ~ 60-70%, high ~ 70-80%.
.censor_alpha2 <- c(low = 2, lowmid = 1, mid = 0, highmid = -1, high = -2)

# Misspecified eta2 (mis2) mean: gamma* = 0 if gamma != 0, else gamma* = 2.
.gamma_star <- function(gamma) if (gamma != 0) 0 else 2

# Misspecified eta1 mean used by mis1 and by the r1* residual center.
.ALPHA1_MIS <- -2
