# =============================================================================
#  Beta estimation for prescco: SPARCC and complete-case
#  ---------------------------------------------------------------------------
#  Two estimators of the outcome coefficients:
#
#    * SPARCC  -- SPARCC estimator, one file per
#                 (censoring level, nuisance specification). Prerequisite for
#                 the PRESCCO prediction.
#    * CC      -- complete-case estimator, one file per censoring level.
#                 A comparison estimator; nothing downstream consumes it.
#
#  ---------------------------------------------------------------------------
#  SPARCC. Prerequisite for the PRESCCO method: produces the
#  beta_param12_<tag>.Rdata files consumed by the PRESCCO prediction stage
#  (run_scenario). Each file contains a matrix `result_beta` with
#  rows = length(beta) and columns = M replicates, column k being the SPARCC
#  estimate for replicate k under the given nuisance specification.
# =============================================================================

#' Estimate beta over M replicates for one nuisance specification
#'
#' Loops \code{get_beta_param_param12} over \code{M} replicates and writes
#' \code{beta_param12_<tag>.Rdata} containing the matrix \code{result_beta}.
#'
#' @param tag Output tag, e.g. \code{"highmid"}, \code{"highmid_mis1"}.
#' @param alpha2 True censoring-time mean (gamma) used to generate the data.
#' @param alpha1_star,alpha2_star Nuisance means used in the estimating
#'   equation; set equal to \code{alpha1}/\code{alpha2} for the correct
#'   specification, or perturbed for mis1/mis2.
#' @param M Number of Monte Carlo replicates.
#' @param alpha1 True mean of X.
#' @param beta True coefficient vector.
#' @param n,d,sigma,tau1,tau2,tt,m Model and numerical settings.
#' @param out_dir Directory for the output file.
#' @param verbose If \code{TRUE}, print progress every 100 replicates.
#' @return Invisibly, the \code{result_beta} matrix (also saved to disk).
#' @export
run_beta_scenario <- function(tag, alpha2,
                              alpha1_star, alpha2_star,
                              M = 1000, alpha1 = 0, beta = c(0, 3),
                              n = 1000, d = 1, sigma = 4, tau1 = 1, tau2 = 1,
                              tt = 20, m = 20,
                              out_dir = ".", verbose = TRUE) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  result_beta <- matrix(0, nrow = length(beta), ncol = M)
  for (k in seq_len(M)) {
    set.seed(k)
    result_beta[, k] <- get_beta_param_param12(
      k, n, d, beta, alpha1, alpha2,
      alpha1_star = alpha1_star, alpha2_star = alpha2_star,
      sigma = sigma, tau1 = tau1, tau2 = tau2, tt = tt, m = m)
    if (verbose && k %% 100 == 0) message("  beta replicate ", k, "/", M)
  }
  out_file <- file.path(out_dir, paste0("beta_param12_", tag, ".Rdata"))
  save(result_beta, file = out_file)
  invisible(result_beta)
}

#' Build the beta-estimation configuration table
#'
#' Every censoring level x nuisance specification the PRESCCO method
#' consumes: 5 levels x \{correct, mis1, mis2\} = 15 files. Beta does not
#' depend on the residual, so the \code{r1}, \code{r2}, and \code{r1*} runs at
#' a given (level, spec) all read the same file. The misspecified-eta2 mean
#' follows \code{gamma* = 0} when \code{gamma != 0} and \code{gamma* = 2} when
#' \code{gamma == 0}.
#'
#' @param alpha1 True mean of X.
#' @return A data frame; each row is a set of arguments for
#'   \code{run_beta_scenario}.
#' @export
build_beta_scenarios <- function(alpha1 = 0) {
  rows <- list()
  add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

  for (censor in names(.censor_alpha2)) {
    gamma <- .censor_alpha2[[censor]]
    for (spec in c("correct", "mis1", "mis2")) {
      if (spec == "correct") {
        a1s <- alpha1;       a2s <- gamma;              tag <- censor
      } else if (spec == "mis1") {
        a1s <- .ALPHA1_MIS;  a2s <- gamma;              tag <- paste0(censor, "_mis1")
      } else {                                        # mis2
        a1s <- alpha1;       a2s <- .gamma_star(gamma); tag <- paste0(censor, "_mis2")
      }
      add(tag = tag, alpha2 = gamma, alpha1_star = a1s, alpha2_star = a2s)
    }
  }
  do.call(rbind, rows)
}

#' Run all (or a subset of) beta-estimation scenarios
#'
#' @param out_dir Directory for the \code{beta_param12_*.Rdata} outputs.
#' @param which Optional integer vector selecting rows of
#'   \code{build_beta_scenarios} to run (default: all).
#' @param alpha1 True mean of X.
#' @param ... Extra arguments forwarded to \code{run_beta_scenario}
#'   (e.g. \code{M}, \code{verbose}).
#' @return Invisibly, the beta-scenario table that was run.
#' @export
run_all_beta_scenarios <- function(out_dir = ".", which = NULL, alpha1 = 0, ...) {
  scenarios <- build_beta_scenarios(alpha1 = alpha1)
  if (!is.null(which)) scenarios <- scenarios[which, , drop = FALSE]
  for (i in seq_len(nrow(scenarios))) {
    s <- scenarios[i, ]
    message(sprintf("[beta %d/%d] %s  (alpha2=%s, alpha1*=%s, alpha2*=%s)",
                    i, nrow(scenarios), s$tag, s$alpha2, s$alpha1_star, s$alpha2_star))
    run_beta_scenario(tag = s$tag, alpha2 = s$alpha2,
                      alpha1_star = s$alpha1_star, alpha2_star = s$alpha2_star,
                      alpha1 = alpha1, out_dir = out_dir, ...)
  }
  invisible(scenarios)
}


# =============================================================================
#  Complete-case (CC) beta estimation
#  ---------------------------------------------------------------------------
#  Ordinary least squares of Y on (1, W) over the uncensored (delta == 1)
#  observations. This is the same estimator the full-conformal and jackknife+
#  methods fit internally. It does not depend on the nuisance specification, so
#  there is one file per censoring level rather than one per (level, spec).
#
#  CC consumes no random numbers: replicate k's data comes from
#  data_generating(k, ...), which sets its own seed, so these results are
#  independent of loop order and of anything else in the session.
# =============================================================================

#' Complete-case beta for one replicate
#'
#' @param k Integer seed / replicate index.
#' @param n,d Sample size and covariate dimension.
#' @param beta,alpha1,alpha2 Data-generating parameters.
#' @param sigma,tau1,tau2 Model standard deviations.
#' @return The complete-case coefficient vector.
#' @export
sim_get_beta_cc <- function(k, n, d, beta, alpha1, alpha2,
                            sigma = 4, tau1 = 1, tau2 = 1) {
  data_k <- data_generating(k, n, d, beta, alpha1, alpha2, sigma, tau1, tau2)
  cc_idx <- (data_k$delta_data == 1)
  X <- cbind(1, data_k$w_data[cc_idx])
  Y <- data_k$y_data[cc_idx]
  as.vector(solve(crossprod(X), crossprod(X, Y)))
}

#' Estimate CC beta over M replicates for one censoring level
#'
#' Writes \code{beta_cc_<tag>.Rdata} containing the matrix \code{result_beta}
#' (rows = \code{length(beta)}, columns = \code{M} replicates), matching the
#' layout of the SPARCC beta files.
#'
#' @param tag Output tag, i.e. the censoring level.
#' @param alpha2 True censoring-time mean (gamma) used to generate the data.
#' @param M Number of Monte Carlo replicates.
#' @param alpha1 True mean of X.
#' @param beta True coefficient vector.
#' @param n,d,sigma,tau1,tau2 Model settings.
#' @param out_dir Directory for the output file.
#' @param verbose If \code{TRUE}, print progress every 100 replicates.
#' @return Invisibly, the \code{result_beta} matrix (also saved to disk).
#' @export
run_cc_scenario <- function(tag, alpha2,
                            M = 1000, alpha1 = 0, beta = c(0, 3),
                            n = 1000, d = 1, sigma = 4, tau1 = 1, tau2 = 1,
                            out_dir = ".", verbose = TRUE) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  result_beta <- matrix(0, nrow = length(beta), ncol = M)
  for (k in seq_len(M)) {
    result_beta[, k] <- sim_get_beta_cc(k, n, d, beta, alpha1, alpha2,
                                        sigma = sigma, tau1 = tau1, tau2 = tau2)
    if (verbose && k %% 100 == 0) message("  cc replicate ", k, "/", M)
  }
  out_file <- file.path(out_dir, paste0("beta_cc_", tag, ".Rdata"))
  save(result_beta, file = out_file)
  invisible(result_beta)
}

#' Run CC beta estimation at every censoring level
#'
#' @param out_dir Directory for the \code{beta_cc_*.Rdata} outputs.
#' @param which Optional integer vector selecting censoring levels (1-5).
#' @param alpha1 True mean of X.
#' @param ... Extra arguments forwarded to \code{run_cc_scenario}.
#' @return Invisibly, the censoring levels that were run.
#' @export
run_all_cc_scenarios <- function(out_dir = ".", which = NULL, alpha1 = 0, ...) {
  levels <- names(.censor_alpha2)
  idx <- if (is.null(which)) seq_along(levels) else which
  for (i in idx) {
    censor <- levels[[i]]
    message(sprintf("[cc %d/%d] %s  (alpha2=%s)", i, length(levels),
                    censor, .censor_alpha2[[censor]]))
    run_cc_scenario(tag = censor, alpha2 = .censor_alpha2[[censor]],
                    alpha1 = alpha1, out_dir = out_dir, ...)
  }
  invisible(levels[idx])
}
