#' Run one PRESCCO method simulation scenario
#'
#' Reproduces a single driver from the original simulation scripts: it builds
#' the \eqn{b_1, b_2, b_3} arrays on a grid of half-length values \eqn{\zeta}
#' and \eqn{\beta}, then solves the estimating equation for \eqn{\widehat\zeta}
#' across all Monte Carlo replicates. A scenario is one combination of
#' \emph{residual} (\code{r1}, \code{r2}, or the misspecified-center
#' \code{r1*}), \emph{nuisance specification} (correct, \code{mis1}, or
#' \code{mis2}), and \emph{censoring level} (via \code{alpha2}).
#'
#' The function expects an upstream file \code{beta_param12_<beta_tag>.Rdata}
#' (produced separately) that contains a matrix \code{result_beta} of dimension
#' \code{length(beta) x M}, whose column \code{k} is the estimated coefficient
#' vector for replicate \code{k}. The \eqn{\beta} interpolation grid is built
#' from its per-coordinate mean +/- 3 SD, matching the original scripts.
#'
#' Two \code{.Rdata} files are written to \code{out_dir}: an intermediate with
#' the \code{b*} arrays and \code{zeta_seq}, and a final one additionally
#' containing \code{zeta_sol} (the estimated half-length per replicate) and
#' \code{zeta_vals} (the full objective grid per replicate).
#'
#' @param residual The residual function, one of \code{r1} or \code{r2}.
#' @param residual_name Short label used in the output filename
#'   (\code{"r1"}, \code{"r2"}, or \code{"r1"} for the r1* runs, matching the
#'   original file naming).
#' @param censor_tag Censoring-level tag used to locate the input beta file and
#'   name outputs: one of \code{"low"}, \code{"lowmid"}, \code{"mid"},
#'   \code{"highmid"}, \code{"high"}.
#' @param spec_tag Output tag: \code{""} (correct), \code{"mis1"},
#'   \code{"mis2"}, or the \code{r1*} variants \code{"r1star"},
#'   \code{"r1star_mis1"}, \code{"r1star_mis2"}.
#' @param beta_spec_tag Nuisance-specification tag used to locate the upstream
#'   beta file: \code{""}, \code{"mis1"}, or \code{"mis2"}. Defaults to
#'   \code{spec_tag} with any \code{"r1star"} prefix stripped, since beta does
#'   not depend on the residual.
#' @param alpha2 Censoring-time mean used to generate the evaluation data
#'   (the true \eqn{\gamma}).
#' @param alpha1_star,alpha2_star Nuisance means passed to the \eqn{b} solvers
#'   and to \code{get_zeta_param_param12_int} (set to misspecified values to
#'   reproduce mis1/mis2).
#' @param alpha1_star_r Center used inside the residual; equals \code{alpha1}
#'   for r1/r2 and the misspecified value (e.g. \code{-2}) for r1*.
#' @param alpha Miscoverage level (\code{0.1} in the paper).
#' @param alpha1 True mean of X (0 in the simulations).
#' @param beta,sigma,tau1,tau2 Model parameters (paper defaults).
#' @param n,d,m Sample size, covariate dimension, and grid resolution for
#'   \code{x_a}.
#' @param n_rep Number of Monte Carlo replicates to solve (1000 in the paper).
#' @param tt,tt2,tt3 Numerical grid sizes forwarded to the \eqn{b} solvers.
#' @param seq_length Number of \eqn{\zeta} grid points.
#' @param in_dir,out_dir Directories for the input beta file and the outputs.
#' @param verbose If \code{TRUE}, print progress.
#'
#' @return Invisibly, a list with \code{beta_array}, \code{zeta_seq},
#'   \code{zeta_sol}, and \code{zeta_vals}. Called mainly for its side effect of
#'   writing the two \code{.Rdata} files.
#' @export
run_scenario <- function(residual, residual_name,
                         censor_tag, spec_tag = "",
                         beta_spec_tag = NULL,
                         alpha2,
                         alpha1_star, alpha2_star, alpha1_star_r,
                         alpha = 0.1, alpha1 = 0,
                         beta = c(0, 3), sigma = 4, tau1 = 1, tau2 = 1,
                         n = 1000, d = 1, m = 20, n_rep = 1000,
                         tt = 100, tt2 = 500, tt3 = 200,
                         seq_length = 5,
                         in_dir = ".", out_dir = ".",
                         verbose = TRUE) {

  x_a <- seq(-1, 1, length.out = m)

  ## --- locate and load the upstream beta estimates -------------------------
  # The residual label ("r1star") is part of the OUTPUT tag only. Beta is
  # estimated per (censoring, nuisance spec) and does not depend on which
  # residual is used, so it must be stripped from the INPUT tag:
  #   spec_tag ""             -> beta_param12_<censor>.Rdata
  #   spec_tag "mis1"         -> beta_param12_<censor>_mis1.Rdata
  #   spec_tag "r1star"       -> beta_param12_<censor>.Rdata
  #   spec_tag "r1star_mis1"  -> beta_param12_<censor>_mis1.Rdata
  beta_spec <- if (!is.null(beta_spec_tag)) beta_spec_tag
               else sub("^r1star_?", "", spec_tag)
  beta_tag <- if (nzchar(beta_spec)) paste0(censor_tag, "_", beta_spec) else censor_tag
  in_file  <- file.path(in_dir, paste0("beta_param12_", beta_tag, ".Rdata"))
  if (!file.exists(in_file)) {
    stop("Upstream beta file not found: ", in_file,
         "\nRun the beta-estimation step first to create it.")
  }
  env <- new.env()
  load(in_file, envir = env)
  if (!exists("result_beta", envir = env)) {
    stop("`result_beta` not found in ", in_file)
  }
  result_beta <- get("result_beta", envir = env)

  ## --- beta interpolation grid: per-coordinate mean +/- 3 SD ----------------
  beta_array <- apply(result_beta, 1, mean) +
    apply(result_beta, 1, sd) %*% t(-3:3)

  ## --- zeta search grid -----------------------------------------------------
  zeta_seq <- sim_get_zeta_seq_r(alpha, alpha1, alpha2, residual,
                             alpha1_star_r = alpha1_star_r,
                             beta = beta, sigma = sigma,
                             tau1 = tau1, tau2 = tau2, d = d,
                             length = seq_length)

  n_zeta <- length(zeta_seq)
  n_grid <- length(beta_array[1, ])

  b1_array <- array(0, dim = c(n_zeta, length(x_a), n_grid, n_grid))
  b2_array <- array(0, dim = c(n_zeta, length(x_a), n_grid, n_grid))
  b3_array <- array(0, dim = c(n_zeta, n_grid, n_grid))

  ## --- build b1/b2/b3 over the (zeta, beta_2) grid --------------------------
  for (l in seq_len(n_zeta)) {
    zeta <- zeta_seq[l]
    if (verbose) message("zeta grid point ", l, "/", n_zeta)
    for (k in seq_len(n_grid)) {
      beta_grid <- c(beta_array[1, 1], beta_array[2, k])

      b1 <- sim_b1_gauss_param12(zeta, alpha, residual, beta_grid, x_a,
                             alpha1_star = alpha1_star,
                             alpha2_star = alpha2_star,
                             alpha1_star_r = alpha1_star_r,
                             sigma = sigma, tau1 = tau1, tau2 = tau2,
                             tt = tt, tt2 = tt2)
      b2 <- sim_b2_gauss_param12(zeta, alpha, residual, beta_grid, x_a,
                             alpha1_star = alpha1_star,
                             alpha2_star = alpha2_star,
                             alpha1_star_r = alpha1_star_r,
                             sigma = sigma, tau1 = tau1, tau2 = tau2,
                             tt = tt, tt2 = tt2)
      b3 <- sim_b3_gauss_param12(zeta, alpha, residual, beta_grid, x_a,
                             alpha1_star = alpha1_star,
                             alpha2_star = alpha2_star,
                             alpha1_star_r = alpha1_star_r,
                             sigma = sigma, tau1 = tau1, tau2 = tau2,
                             tt = 20, tt2 = tt2, tt3 = tt3)

      for (j in seq_len(n_grid)) {
        b1_array[l, , j, k] <- b1
        b2_array[l, , j, k] <- b2
        b3_array[l, j, k]   <- b3
      }
    }
  }

  ## --- output filename ------------------------------------------------------
  out_tag  <- if (nzchar(spec_tag)) paste0(censor_tag, "_", spec_tag) else censor_tag
  out_file <- file.path(
    out_dir,
    paste0("zeta_param12_", residual_name, "_", alpha, "_", out_tag, ".Rdata"))

  save(beta_array, b1_array, b2_array, b3_array, zeta_seq, file = out_file)

  ## --- solve for zeta at each replicate -------------------------------------
  zeta_vals <- vector("list", n_rep)
  zeta_sol  <- rep(0, n_rep)
  for (k in seq_len(n_rep)) {
    beta_temp <- result_beta[, k]
    zeta_vals[[k]] <- get_zeta_param_param12_int(
      k, n, d, alpha, residual, beta, alpha1, alpha2,
      b1_array, b2_array, b3_array, beta_array, beta_temp, zeta_seq,
      alpha1_star = alpha1_star, alpha2_star = alpha2_star,
      sigma = sigma, tau1 = tau1, tau2 = tau2, m = m)
    zeta_sol[k] <- zeta_vals[[k]]$sol
    if (verbose && k %% 100 == 0) message("solved replicate ", k, "/", n_rep)
  }

  save(beta_array, b1_array, b2_array, b3_array, zeta_seq,
       zeta_sol, zeta_vals, file = out_file)

  invisible(list(beta_array = beta_array, zeta_seq = zeta_seq,
                 zeta_sol = zeta_sol, zeta_vals = zeta_vals))
}


# ---- Scenario configuration and drivers -------------------------------------

#' Build the full simulation configuration table
#'
#' The full factorial: 5 censoring levels (\code{low}, \code{lowmid},
#' \code{mid}, \code{highmid}, \code{high}) x 3 residuals (\code{r1},
#' \code{r2}, \code{r1*}) x 3 nuisance specifications (correct, \code{mis1},
#' \code{mis2}) = 45 runs. \code{mis1} perturbs the \code{eta1} mean to
#' \code{.ALPHA1_MIS}; \code{mis2} perturbs the \code{eta2} mean to
#' \code{.gamma_star(gamma)}. The \code{r1*} runs additionally use a
#' misspecified residual center (\code{.ALPHA1_MIS}) under all three
#' specifications.
#'
#' @param alpha1 True mean of X (0 in the paper).
#' @return A data frame; each row is a set of arguments for \code{run_scenario}.
#'   The \code{residual} column holds the residual-function name as a string
#'   (resolve with \code{get()} at call time).
#' @export
build_scenarios <- function(alpha1 = 0) {
  rows <- list()
  add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

  ## ---- 5 censoring x 3 residuals x 3 nuisance specifications = 45 runs -----
  for (censor in names(.censor_alpha2)) {
    gamma <- .censor_alpha2[[censor]]
    for (res in c("r1", "r2", "r1star")) {
      # r1* uses residual r1 with a misspecified center; r1/r2 use their own.
      res_fun     <- if (res == "r1star") "r1" else res
      center_r    <- if (res == "r1star") .ALPHA1_MIS else alpha1
      # spec_tag prefix in the output filename ("" for r1/r2, "r1star" for r1*).
      tag_prefix  <- if (res == "r1star") "r1star" else ""

      specs <- list(
        correct = list(a1s = alpha1,       a2s = gamma),
        mis1    = list(a1s = .ALPHA1_MIS,  a2s = gamma),
        mis2    = list(a1s = alpha1,       a2s = .gamma_star(gamma))
      )
      for (spec in names(specs)) {
        s <- specs[[spec]]
        spec_tag <- if (spec == "correct") tag_prefix
                    else if (nzchar(tag_prefix)) paste0(tag_prefix, "_", spec)
                    else spec
        # Beta file depends only on (censoring, nuisance spec), not the residual.
        beta_spec_tag <- if (spec == "correct") "" else spec
        add(residual = res_fun, residual_name = res_fun,
            censor_tag = censor, spec_tag = spec_tag,
            beta_spec_tag = beta_spec_tag,
            alpha2 = gamma,
            alpha1_star = s$a1s, alpha2_star = s$a2s, alpha1_star_r = center_r)
      }
    }
  }

  do.call(rbind, rows)
}

#' Run all (or a subset of) simulation scenarios
#'
#' @param in_dir Directory holding the upstream \code{beta_param12_*.Rdata}.
#' @param out_dir Directory for the \code{zeta_param12_*.Rdata} outputs.
#' @param which Optional integer vector selecting rows of \code{build_scenarios}
#'   to run (default: all).
#' @param alpha1 True mean of X.
#' @param ... Extra arguments forwarded to \code{run_scenario} (e.g. \code{n_rep},
#'   \code{tt}, \code{verbose}).
#' @return Invisibly, the scenario table that was run.
#' @export
run_all_scenarios <- function(in_dir = ".", out_dir = ".",
                              which = NULL, alpha1 = 0, ...) {
  scenarios <- build_scenarios(alpha1 = alpha1)
  if (!is.null(which)) scenarios <- scenarios[which, , drop = FALSE]
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  for (i in seq_len(nrow(scenarios))) {
    s <- scenarios[i, ]
    message(sprintf("[%d/%d] %s  censor=%s  spec=%s",
                    i, nrow(scenarios), s$residual_name,
                    s$censor_tag, ifelse(nzchar(s$spec_tag), s$spec_tag, "correct")))
    run_scenario(
      residual      = get(paste0("sim_", s$residual)),
      residual_name = s$residual_name,
      censor_tag    = s$censor_tag,
      spec_tag      = s$spec_tag,
      beta_spec_tag = s$beta_spec_tag,
      alpha2        = s$alpha2,
      alpha1_star   = s$alpha1_star,
      alpha2_star   = s$alpha2_star,
      alpha1_star_r = s$alpha1_star_r,
      alpha1        = alpha1,
      in_dir        = in_dir,
      out_dir       = out_dir,
      ...)
  }
  invisible(scenarios)
}
