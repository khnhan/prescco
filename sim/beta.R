# =============================================================================
# Beta estimation for PRESCCO: SPARCC and complete-case
#
# Run from the package root:
#
#   Rscript sim/beta.R
#
# Outputs are saved in:
#
#   sim/Rdata_files/
#
# SPARCC:
#   beta_param12_<censor>.Rdata
#   beta_param12_<censor>_mis1.Rdata
#   beta_param12_<censor>_mis2.Rdata
#
# Complete case:
#   beta_cc_<censor>.Rdata
#
# Each file contains `result_beta`, with rows corresponding to coefficients
# and columns corresponding to Monte Carlo replicates.
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

sigma <- 4
tau1 <- 1
tau2 <- 1

tt <- 20
m <- 20


# =============================================================================
# SPARCC beta estimation
# =============================================================================
#
# For each censoring level, three nuisance specifications are considered:
#
#   correct:  alpha1* = alpha1,       alpha2* = alpha2
#   mis1:     alpha1* = .ALPHA1_MIS,  alpha2* = alpha2
#   mis2:     alpha1* = alpha1,       alpha2* = .gamma_star(alpha2)
#
# The same beta file is used for r1, r2, and r1* because beta estimation
# does not depend on the residual definition.
# =============================================================================

for (censor in names(.censor_alpha2)) {

  alpha2 <- .censor_alpha2[[censor]]

  for (spec in c("correct", "mis1", "mis2")) {

    if (spec == "correct") {

      alpha1_star <- alpha1
      alpha2_star <- alpha2
      tag <- censor

    } else if (spec == "mis1") {

      alpha1_star <- .ALPHA1_MIS
      alpha2_star <- alpha2
      tag <- paste0(censor, "_mis1")

    } else if (spec == "mis2") {

      alpha1_star <- alpha1
      alpha2_star <- .gamma_star(alpha2)
      tag <- paste0(censor, "_mis2")
    }


    message(
      "SPARCC beta: ", tag,
      "  (alpha2 = ", alpha2,
      ", alpha1* = ", alpha1_star,
      ", alpha2* = ", alpha2_star, ")"
    )


    result_beta <- matrix(
      0,
      nrow = length(beta),
      ncol = M
    )


    for (k in seq_len(M)) {

      set.seed(k)

      result_beta[, k] <- get_beta_param_param12(
        k,
        n,
        d,
        beta,
        alpha1,
        alpha2,
        alpha1_star = alpha1_star,
        alpha2_star = alpha2_star,
        sigma = sigma,
        tau1 = tau1,
        tau2 = tau2,
        tt = tt,
        m = m
      )

      if (k %% 100 == 0) {
        message("  replicate ", k, "/", M)
      }
    }


    save(
      result_beta,
      file = file.path(
        RDATA_DIR,
        paste0("beta_param12_", tag, ".Rdata")
      )
    )
  }
}


# =============================================================================
# Complete-case beta estimation
# =============================================================================
#
# Ordinary least squares of Y on (1, W) among observations with delta = 1.
#
# There is one file for each censoring level because the complete-case
# estimator does not involve either nuisance specification.
#
# Replicate k uses data_generating(k, ...), exactly as in the original code.
# =============================================================================

for (censor in names(.censor_alpha2)) {

  alpha2 <- .censor_alpha2[[censor]]

  message(
    "Complete-case beta: ", censor,
    "  (alpha2 = ", alpha2, ")"
  )


  result_beta <- matrix(
    0,
    nrow = length(beta),
    ncol = M
  )


  for (k in seq_len(M)) {

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

    cc_idx <- data_k$delta_data == 1

    X <- cbind(
      1,
      data_k$w_data[cc_idx]
    )

    Y <- data_k$y_data[cc_idx]

    result_beta[, k] <- as.vector(
      solve(
        crossprod(X),
        crossprod(X, Y)
      )
    )

    if (k %% 100 == 0) {
      message("  replicate ", k, "/", M)
    }
  }


  save(
    result_beta,
    file = file.path(
      RDATA_DIR,
      paste0("beta_cc_", censor, ".Rdata")
    )
  )
}
