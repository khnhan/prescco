# =============================================================================
# PRESCCO simulation
#
# Run from the package root:
#
#   Rscript sim/PRESCCO.R
#
# Inputs and outputs are stored in:
#
#   sim/Rdata_files/
#
# The simulation considers:
#   - five censoring levels,
#   - three residual definitions: r1, r2, and r1*,
#   - three nuisance specifications: correct, mis1, and mis2.
#
# For each scenario, the script constructs b1, b2, and b3 on the zeta and beta
# grids, solves the estimating equation for each Monte Carlo replicate, and
# saves the resulting half-length estimates.
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
m <- 20

n_rep <- 1000

tt <- 100
tt2 <- 500
tt3 <- 200

seq_length <- 5

x_a <- seq(-1, 1, length.out = m)


# =============================================================================
# Run all PRESCCO scenarios
# =============================================================================

for (censor in names(.censor_alpha2)) {

  alpha2 <- .censor_alpha2[[censor]]

  for (r_char in c("r1", "r2", "r1star")) {

    if (r_char == "r1") {

      residual <- sim_r1
      residual_name <- "r1"
      alpha1_star_r <- alpha1
      residual_suffix <- ""

    } else if (r_char == "r2") {

      residual <- sim_r2
      residual_name <- "r2"
      alpha1_star_r <- alpha1
      residual_suffix <- ""

    } else {

      residual <- sim_r1
      residual_name <- "r1"
      alpha1_star_r <- .ALPHA1_MIS
      residual_suffix <- "_r1star"
    }


    for (spec in c("correct", "mis1", "mis2")) {

      if (spec == "correct") {

        alpha1_star <- alpha1
        alpha2_star <- alpha2
        spec_suffix <- ""

      } else if (spec == "mis1") {

        alpha1_star <- .ALPHA1_MIS
        alpha2_star <- alpha2
        spec_suffix <- "_mis1"

      } else {

        alpha1_star <- alpha1
        alpha2_star <- .gamma_star(alpha2)
        spec_suffix <- "_mis2"
      }


      message(
        "PRESCCO: ",
        censor,
        "  ",
        r_char,
        "  ",
        spec
      )


      # -----------------------------------------------------------------------
      # Load beta estimates
      # -----------------------------------------------------------------------

      beta_file <- file.path(
        RDATA_DIR,
        paste0(
          "beta_param12_",
          censor,
          spec_suffix,
          ".Rdata"
        )
      )

      load(beta_file)


      # -----------------------------------------------------------------------
      # Beta interpolation grid
      # -----------------------------------------------------------------------

      beta_array <- apply(
        result_beta,
        1,
        mean
      ) +
        apply(
          result_beta,
          1,
          sd
        ) %*%
        t(-3:3)


      # -----------------------------------------------------------------------
      # Zeta search grid
      # -----------------------------------------------------------------------

      zeta_seq <- sim_get_zeta_seq_r(
        alpha,
        alpha1,
        alpha2,
        residual,
        alpha1_star_r = alpha1_star_r,
        beta = beta,
        sigma = sigma,
        tau1 = tau1,
        tau2 = tau2,
        d = d,
        length = seq_length
      )

      n_zeta <- length(zeta_seq)
      n_grid <- length(beta_array[1, ])


      # -----------------------------------------------------------------------
      # Construct b1, b2, and b3
      # -----------------------------------------------------------------------

      b1_array <- array(
        0,
        dim = c(
          n_zeta,
          length(x_a),
          n_grid,
          n_grid
        )
      )

      b2_array <- array(
        0,
        dim = c(
          n_zeta,
          length(x_a),
          n_grid,
          n_grid
        )
      )

      b3_array <- array(
        0,
        dim = c(
          n_zeta,
          n_grid,
          n_grid
        )
      )


      for (l in seq_len(n_zeta)) {

        zeta <- zeta_seq[l]

        message(
          "  zeta grid point ",
          l,
          "/",
          n_zeta
        )


        for (k in seq_len(n_grid)) {

          beta_grid <- c(
            beta_array[1, 1],
            beta_array[2, k]
          )


          b1 <- sim_b1_gauss_param12(
            zeta,
            alpha,
            residual,
            beta_grid,
            x_a,
            alpha1_star = alpha1_star,
            alpha2_star = alpha2_star,
            alpha1_star_r = alpha1_star_r,
            sigma = sigma,
            tau1 = tau1,
            tau2 = tau2,
            tt = tt,
            tt2 = tt2
          )


          b2 <- sim_b2_gauss_param12(
            zeta,
            alpha,
            residual,
            beta_grid,
            x_a,
            alpha1_star = alpha1_star,
            alpha2_star = alpha2_star,
            alpha1_star_r = alpha1_star_r,
            sigma = sigma,
            tau1 = tau1,
            tau2 = tau2,
            tt = tt,
            tt2 = tt2
          )


          b3 <- sim_b3_gauss_param12(
            zeta,
            alpha,
            residual,
            beta_grid,
            x_a,
            alpha1_star = alpha1_star,
            alpha2_star = alpha2_star,
            alpha1_star_r = alpha1_star_r,
            sigma = sigma,
            tau1 = tau1,
            tau2 = tau2,
            tt = 20,
            tt2 = tt2,
            tt3 = tt3
          )


          for (j in seq_len(n_grid)) {

            b1_array[l, , j, k] <- b1
            b2_array[l, , j, k] <- b2
            b3_array[l, j, k] <- b3
          }
        }
      }


      # -----------------------------------------------------------------------
      # Solve for zeta in each Monte Carlo replicate
      # -----------------------------------------------------------------------

      zeta_vals <- vector(
        "list",
        n_rep
      )

      zeta_sol <- rep(
        0,
        n_rep
      )


      for (k in seq_len(n_rep)) {

        beta_temp <- result_beta[, k]

        zeta_vals[[k]] <- get_zeta_param_param12_int(
          k,
          n,
          d,
          alpha,
          residual,
          beta,
          alpha1,
          alpha2,
          b1_array,
          b2_array,
          b3_array,
          beta_array,
          beta_temp,
          zeta_seq,
          alpha1_star = alpha1_star,
          alpha2_star = alpha2_star,
          sigma = sigma,
          tau1 = tau1,
          tau2 = tau2,
          m = m
        )

        zeta_sol[k] <- zeta_vals[[k]]$sol

        if (k %% 100 == 0) {
          message(
            "  replicate ",
            k,
            "/",
            n_rep
          )
        }
      }


      # -----------------------------------------------------------------------
      # Save results
      # -----------------------------------------------------------------------

      out_file <- file.path(
        RDATA_DIR,
        paste0(
          "zeta_param12_",
          residual_name, "_",
          alpha, "_",
          censor,
          residual_suffix,
          spec_suffix,
          ".Rdata"
        )
      )


      save(
        beta_array,
        b1_array,
        b2_array,
        b3_array,
        zeta_seq,
        zeta_sol,
        zeta_vals,
        file = out_file
      )
    }
  }
}
