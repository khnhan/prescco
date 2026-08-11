# =============================================================================
#  Prediction coverage rate for the semiparametric method
#  ---------------------------------------------------------------------------
#  Computes the empirical prediction coverage rate of a fitted half-length
#  sequence on a fresh test set. This stage exists only for the semiparametric
#  method (correct / mis1 / mis2), which stores a `zeta_sol` and a
#  `result_beta` but does not evaluate a coverage rate of its own.
#
#  All three conformal methods -- split conformal, full conformal, and
#  jackknife+ -- compute and save their own coverage rate in conformal.R, on
#  the same test set used here: seed `k + seed_offset`, size `N`, same
#  data-generating parameters.
# =============================================================================

#' Empirical prediction coverage rate for a fitted half-length sequence
#'
#' For each replicate \code{k}, generates an independent test set of size
#' \code{N} (seed \code{k + seed_offset}), evaluates the residual with the
#' replicate's fitted coefficients \code{result_beta[, k]}, and records the
#' fraction of test residuals not exceeding \code{zeta_sol[k]}.
#'
#' @param zeta_sol Numeric vector of per-replicate half-lengths.
#' @param result_beta Coefficient matrix, columns indexed by replicate.
#' @param r_char Residual type: \code{"r1"}, \code{"r2"}, or \code{"r1star"}.
#' @param alpha2 True censoring-time mean (gamma) for the test data.
#' @param alpha1 True mean of X.
#' @param beta True coefficient vector (used only to generate test data).
#' @param sigma,tau1,tau2 Model standard deviations.
#' @param n,d Sample size and covariate dimension for the test generator.
#' @param N Test-set size.
#' @param alpha1_star_r Misspecified center for the r1* residual (default -2).
#' @param seed_offset Offset added to the replicate index for the test seed.
#' @param verbose If \code{TRUE}, print progress every 100 replicates.
#' @return Numeric vector of empirical coverage rates, one per replicate.
#' @export
get_pred_cvg <- function(zeta_sol, result_beta, r_char, alpha2,
                         alpha1 = 0, beta = c(0, 3),
                         sigma = 4, tau1 = 1, tau2 = 1,
                         n = 1000, d = 1, N = 10000,
                         alpha1_star_r = -2, seed_offset = 98765,
                         verbose = FALSE) {
  r_char <- match.arg(r_char, c("r1", "r2", "r1star"))
  center <- if (r_char == "r1star") alpha1_star_r else alpha1
  res_fun <- if (r_char == "r2") sim_r2_vec else sim_r1_vec

  M <- length(zeta_sol)
  pred_cvg <- rep(0, M)
  for (k in seq_len(M)) {
    test <- data_generating(k + seed_offset, N, d, beta, alpha1, alpha2,
                            sigma = sigma, tau1 = tau1, tau2 = tau2)
    beta_hat <- result_beta[, k]
    r_vec <- res_fun(test$y_data, test$w_data, test$delta_data,
                     beta_hat, center, tau1)
    pred_cvg[k] <- mean(r_vec <= zeta_sol[k])
    if (verbose && k %% 100 == 0) message("  coverage replicate ", k, "/", M)
  }
  pred_cvg
}

#' Evaluate and save prediction coverage rate for one censoring level
#'
#' Loads the fitted \code{zeta_sol} and \code{result_beta} for the
#' semiparametric method (true / mis1 / mis2) and split conformal (cp), computes
#' coverage rate with \code{get_pred_cvg}, and writes
#' \code{pred_cvg_<r>_0.1_<censor>[...].Rdata} for each residual and method.
#'
#' @param censor_tag Censoring level: \code{"low"}, \code{"lowmid"},
#'   \code{"mid"}, \code{"highmid"}, or \code{"high"}.
#' @param alpha2 True censoring-time mean (gamma) for this level.
#' @param alpha Miscoverage level (used only in the filename, default 0.1).
#' @param in_dir Directory holding the zeta \code{.Rdata} inputs (and the split
#'   conformal beta file, which is written alongside them).
#' @param out_dir Directory for the \code{pred_cvg_*} outputs.
#' @param beta_dir Directory holding the SPARCC \code{beta_param12_*.Rdata}
#'   files; defaults to \code{in_dir}.
#' @param methods Which nuisance specifications to evaluate; subset of
#'   \code{c("true", "mis1", "mis2")}.
#' @param residuals Which residuals to evaluate; subset of
#'   \code{c("r1", "r2", "r1star")}.
#' @param verbose Passed to \code{get_pred_cvg}.
#' @param ... Extra arguments forwarded to \code{get_pred_cvg}
#'   (e.g. \code{N}, \code{beta}).
#' @return Invisibly, a named list of the coverage rate vectors that were saved.
#' @export
run_pred_cvg_scenario <- function(censor_tag, alpha2, alpha = 0.1,
                                  in_dir = ".", out_dir = ".",
                                  beta_dir = in_dir,
                                  methods = c("true", "mis1", "mis2"),
                                  residuals = c("r1", "r2", "r1star"),
                                  verbose = FALSE, ...) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Method -> (beta file suffix, zeta file infix before ".Rdata")
  beta_suffix <- c(true = "", mis1 = "_mis1", mis2 = "_mis2")
  saved <- list()

  load_one <- function(path, want) {
    env <- new.env(); load(path, envir = env)
    if (!exists(want, envir = env)) stop("`", want, "` not found in ", path)
    get(want, envir = env)
  }

  for (r_char in residuals) {
    rlab   <- if (r_char == "r1star") "r1" else r_char
    infix  <- if (r_char == "r1star") "_r1star" else ""
    for (mth in methods) {
      beta_file <- file.path(beta_dir,
        paste0("beta_param12_", censor_tag, beta_suffix[[mth]], ".Rdata"))
      # zeta filename: zeta_param12_<r>_<alpha>_<censor>[_r1star][_mis1|_mis2]
      zeta_file <- file.path(in_dir, paste0(
        "zeta_param12_", rlab, "_", alpha, "_", censor_tag, infix,
        if (mth == "true") "" else paste0("_", sub("^_", "", beta_suffix[[mth]])),
        ".Rdata"))

      result_beta <- load_one(beta_file, "result_beta")
      zeta_sol    <- load_one(zeta_file, "zeta_sol")

      pred_cvg <- get_pred_cvg(zeta_sol, result_beta, r_char = r_char,
                               alpha2 = alpha2, verbose = verbose, ...)

      out_name <- paste0("pred_cvg_", rlab, "_", alpha, "_", censor_tag, infix,
                         if (mth == "true") "" else paste0("_", sub("^_", "", beta_suffix[[mth]])),
                         ".Rdata")
      save(pred_cvg, file = file.path(out_dir, out_name))
      saved[[out_name]] <- pred_cvg
      if (verbose) message("saved ", out_name)
    }
  }
  invisible(saved)
}

#' Run prediction coverage rate evaluation across all censoring levels
#'
#' @param in_dir Directory with beta/zeta inputs.
#' @param out_dir Directory for the \code{pred_cvg_*} outputs.
#' @param which Optional integer vector selecting censoring levels (1-5) to run.
#' @param ... Forwarded to \code{run_pred_cvg_scenario}.
#' @return Invisibly, the censoring levels that were run.
#' @export
run_all_pred_cvg <- function(in_dir = ".", out_dir = ".", beta_dir = in_dir,
                             which = NULL, ...) {
  levels <- names(.censor_alpha2)
  idx <- if (is.null(which)) seq_along(levels) else which
  for (i in idx) {
    censor <- levels[[i]]
    message(sprintf("[coverage %d/%d] %s", i, length(levels), censor))
    run_pred_cvg_scenario(censor_tag = censor,
                          alpha2 = .censor_alpha2[[censor]],
                          in_dir = in_dir, out_dir = out_dir,
                          beta_dir = beta_dir, ...)
  }
  invisible(levels[idx])
}
