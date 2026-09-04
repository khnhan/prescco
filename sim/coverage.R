# =============================================================================
# Prediction coverage for PRESCCO
#
# Run from the package root:
#
#   Rscript sim/coverage.R
#
# Inputs and outputs are stored in:
#
#   sim/Rdata_files/
#
# For each censoring level and nuisance specification, this script loads the
# fitted beta coefficients and PRESCCO half-lengths, evaluates prediction
# coverage on an independent test set, and saves the resulting coverage rates.
# =============================================================================


# ---- setup -------------------------------------------------------------------

source("sim/core.R")

RDATA_DIR <- "sim/Rdata_files"
dir.create(RDATA_DIR, showWarnings = FALSE, recursive = TRUE)


# ---- simulation settings -----------------------------------------------------

alpha <- 0.1

alpha1 <- 0
beta <- c(0, 3)

sigma <- 4
tau1 <- 1
tau2 <- 1

n <- 1000
d <- 1
N <- 10000

alpha1_star_r <- -2
seed_offset <- 98765


# =============================================================================
# Prediction coverage
# =============================================================================

# zeta_sol: estimated half-lengths across simulation replicates.
# result_beta: estimated beta coefficients, with one column per replicate.
# r_char: residual definition, one of "r1", "r2", or "r1star".
# alpha2: true censoring-model mean for the current censoring level.
# N: size of the independent test set.
# alpha1_star_r: center used for the r1* residual.
# seed_offset: offset used to generate the independent test set.

get_pred_cvg <- function(
    zeta_sol,
    result_beta,
    r_char,
    alpha2,
    alpha1 = 0,
    beta = c(0, 3),
    sigma = 4,
    tau1 = 1,
    tau2 = 1,
    n = 1000,
    d = 1,
    N = 10000,
    alpha1_star_r = -2,
    seed_offset = 98765) {

  center <- if (r_char == "r1star") {
    alpha1_star_r
  } else {
    alpha1
  }

  res_fun <- if (r_char == "r2") {
    sim_r2_vec
  } else {
    sim_r1_vec
  }

  M <- length(zeta_sol)
  pred_cvg <- rep(0, M)

  for (k in seq_len(M)) {

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

    beta_hat <- result_beta[, k]

    r_vec <- res_fun(
      test$y_data,
      test$w_data,
      test$delta_data,
      beta_hat,
      center,
      tau1
    )

    pred_cvg[k] <- mean(
      r_vec <= zeta_sol[k]
    )

    if (k %% 100 == 0) {
      message("  replicate ", k, "/", M)
    }
  }

  pred_cvg
}


# =============================================================================
# Run all censoring levels and nuisance specifications
# =============================================================================

for (censor in names(.censor_alpha2)) {

  alpha2 <- .censor_alpha2[[censor]]

  for (method in c("true", "mis1", "mis2")) {

    suffix <- if (method == "true") {
      ""
    } else if (method == "mis1") {
      "_mis1"
    } else {
      "_mis2"
    }

    message(
      "Coverage: ", censor,
      if (method == "true") "" else paste0(" ", method)
    )


    # Load beta estimates.

    beta_file <- file.path(
      RDATA_DIR,
      paste0(
        "beta_param12_",
        censor,
        suffix,
        ".Rdata"
      )
    )

    load(beta_file)


    for (r_char in c("r1", "r2", "r1star")) {

      rlab <- if (r_char == "r1star") {
        "r1"
      } else {
        r_char
      }

      infix <- if (r_char == "r1star") {
        "_r1star"
      } else {
        ""
      }


      # Load estimated half-lengths.

      zeta_file <- file.path(
        RDATA_DIR,
        paste0(
          "zeta_param12_",
          rlab, "_",
          alpha, "_",
          censor,
          infix,
          suffix,
          ".Rdata"
        )
      )

      load(zeta_file)


      # Evaluate prediction coverage.

      pred_cvg <- get_pred_cvg(
        zeta_sol = zeta_sol,
        result_beta = result_beta,
        r_char = r_char,
        alpha2 = alpha2,
        alpha1 = alpha1,
        beta = beta,
        sigma = sigma,
        tau1 = tau1,
        tau2 = tau2,
        n = n,
        d = d,
        N = N,
        alpha1_star_r = alpha1_star_r,
        seed_offset = seed_offset
      )


      # Save prediction coverage.

      save(
        pred_cvg,
        file = file.path(
          RDATA_DIR,
          paste0(
            "pred_cvg_",
            rlab, "_",
            alpha, "_",
            censor,
            infix,
            suffix,
            ".Rdata"
          )
        )
      )
    }
  }
}
