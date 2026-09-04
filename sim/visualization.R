# =============================================================================
# Visualization / post-processing for PRESCCO
#
# Run from the package root:
#
#   Rscript sim/visualization.R
#
# Input:
#   sim/Rdata_files/
#
# Output:
#   sim/results/
# =============================================================================

library(ggplot2)
library(patchwork)

RDATA_DIR <- "sim/Rdata_files"
RESULTS_DIR <- "sim/results"

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.1

cens_rate_list <- c(
  "low",
  "lowmid",
  "mid",
  "highmid",
  "high"
)

r_char_list <- c(
  "r1",
  "r2",
  "r1star"
)

method_list <- c(
  "true",
  "mis1",
  "mis2",
  "scp",
  "fcp",
  "jackknife"
)


# =============================================================================
# True half-lengths
# =============================================================================

zeta_true_mat <- NULL

true_zeta_file <- file.path(
  RDATA_DIR,
  "true_zeta.Rdata"
)

if (file.exists(true_zeta_file)) {

  load(true_zeta_file)

  zeta_true_mat <- matrix(
    apply(true_zeta[, -(1:2)], 1, median),
    ncol = length(cens_rate_list)
  )

  rownames(zeta_true_mat) <- r_char_list
  colnames(zeta_true_mat) <- cens_rate_list
}


# =============================================================================
# Read all simulation results
# =============================================================================

zeta_all <- NULL
pred_cvg_all <- NULL

for (censor in cens_rate_list) {

  message("Loading: ", censor)

  for (r_char in r_char_list) {

    rlab <- if (r_char == "r1star") "r1" else r_char
    infix <- if (r_char == "r1star") "_r1star" else ""

    for (method in method_list) {

      # -----------------------------------------------------------------------
      # zeta filename
      # -----------------------------------------------------------------------

      if (method == "true") {

        zeta_file <- paste0(
          "zeta_param12_",
          rlab, "_", alpha, "_", censor,
          infix,
          ".Rdata"
        )

      } else if (method == "mis1") {

        zeta_file <- paste0(
          "zeta_param12_",
          rlab, "_", alpha, "_", censor,
          infix,
          "_mis1.Rdata"
        )

      } else if (method == "mis2") {

        zeta_file <- paste0(
          "zeta_param12_",
          rlab, "_", alpha, "_", censor,
          infix,
          "_mis2.Rdata"
        )

      } else if (method == "scp") {

        zeta_file <- paste0(
          "zeta_param12_",
          rlab, "_", alpha, "_", censor,
          infix,
          "_cp.Rdata"
        )

      } else if (method == "fcp") {

        zeta_file <- paste0(
          "zeta_",
          rlab, "_", alpha, "_", censor,
          infix,
          "_fcp.Rdata"
        )

      } else if (method == "jackknife") {

        zeta_file <- paste0(
          "zeta_",
          rlab, "_", alpha, "_", censor,
          infix,
          "_jackknife.Rdata"
        )
      }


      # -----------------------------------------------------------------------
      # prediction coverage filename
      # -----------------------------------------------------------------------

      if (method == "true") {

        pred_file <- paste0(
          "pred_cvg_",
          rlab, "_", alpha, "_", censor,
          infix,
          ".Rdata"
        )

      } else if (method == "mis1") {

        pred_file <- paste0(
          "pred_cvg_",
          rlab, "_", alpha, "_", censor,
          infix,
          "_mis1.Rdata"
        )

      } else if (method == "mis2") {

        pred_file <- paste0(
          "pred_cvg_",
          rlab, "_", alpha, "_", censor,
          infix,
          "_mis2.Rdata"
        )

      } else if (method == "scp") {

        pred_file <- paste0(
          "pred_cvg_",
          rlab, "_", alpha, "_", censor,
          infix,
          "_cp.Rdata"
        )

      } else if (method == "fcp") {

        pred_file <- paste0(
          "pred_cvg_",
          rlab, "_", alpha, "_", censor,
          infix,
          "_fcp.Rdata"
        )

      } else if (method == "jackknife") {

        pred_file <- paste0(
          "pred_cvg_",
          rlab, "_", alpha, "_", censor,
          infix,
          "_jackknife.Rdata"
        )
      }


      # -----------------------------------------------------------------------
      # Load zeta_sol
      # -----------------------------------------------------------------------

      zeta_path <- file.path(
        RDATA_DIR,
        zeta_file
      )

      if (!file.exists(zeta_path)) {
        stop("Missing file: ", zeta_path)
      }

      load(zeta_path)

      zeta_all <- rbind(
        zeta_all,
        data.frame(
          censoring = censor,
          r_char = r_char,
          method = method,
          zeta_sol = zeta_sol
        )
      )


      # -----------------------------------------------------------------------
      # Load pred_cvg
      # -----------------------------------------------------------------------

      pred_path <- file.path(
        RDATA_DIR,
        pred_file
      )

      if (!file.exists(pred_path)) {
        stop("Missing file: ", pred_path)
      }

      load(pred_path)

      pred_cvg_all <- rbind(
        pred_cvg_all,
        data.frame(
          censoring = censor,
          r_char = r_char,
          method = method,
          pred_cvg = pred_cvg
        )
      )
    }
  }
}


# =============================================================================
# Summary table
# =============================================================================

result <- NULL

for (censor in cens_rate_list) {

  for (r_char in r_char_list) {

    for (method in method_list) {

      z <- zeta_all$zeta_sol[
        zeta_all$censoring == censor &
          zeta_all$r_char == r_char &
          zeta_all$method == method
      ]

      p <- pred_cvg_all$pred_cvg[
        pred_cvg_all$censoring == censor &
          pred_cvg_all$r_char == r_char &
          pred_cvg_all$method == method
      ]

      zeta_mean <- mean(z)

      if (!is.null(zeta_true_mat)) {

        zeta_true <- zeta_true_mat[
          r_char,
          censor
        ]

        bias <- zeta_mean - zeta_true

      } else {

        bias <- NA
      }

      result <- rbind(
        result,
        data.frame(
          censoring = censor,
          r_char = r_char,
          method = method,
          zeta_mean = zeta_mean,
          bias = bias,
          emp_sd = sd(z),
          mean_pred_cvg = mean(p),
          sd_pred_cvg = sd(p)
        )
      )
    }
  }
}


result[, c(
  "zeta_mean",
  "bias",
  "emp_sd",
  "mean_pred_cvg",
  "sd_pred_cvg"
)] <- round(
  result[, c(
    "zeta_mean",
    "bias",
    "emp_sd",
    "mean_pred_cvg",
    "sd_pred_cvg"
  )],
  3
)

write.csv(
  result,
  file = file.path(
    RESULTS_DIR,
    "simulation_result.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Figures
# =============================================================================

method_order <- c(
  "true",
  "mis1",
  "mis2",
  "scp",
  "fcp",
  "jackknife"
)

method_labels <- expression(
  "PRESCCO (" * eta[1] * "," * eta[2] * ")",
  "PRESCCO (" * eta[1] * "*" * "," * eta[2] * ")",
  "PRESCCO (" * eta[1] * "," * eta[2]^"\u2605" * ")",
  "Split CP",
  "Full CP",
  "Jackknife+"
)


for (censor in cens_rate_list) {

  message("Plotting: ", censor)


  # ===========================================================================
  # Half-length: r1
  # ===========================================================================

  d1 <- zeta_all[
    zeta_all$censoring == censor &
      zeta_all$r_char == "r1",
  ]

  p1 <- ggplot(
    d1,
    aes(
      y = method,
      x = zeta_sol
    )
  ) +
    geom_boxplot(
      fill = "magenta",
      linewidth = 0.5
    ) +
    scale_y_discrete(
      limits = rev(method_order),
      labels = rev(method_labels)
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_minimal(
      base_size = 12
    )


  # ===========================================================================
  # Half-length: r2
  # ===========================================================================

  d2 <- zeta_all[
    zeta_all$censoring == censor &
      zeta_all$r_char == "r2",
  ]

  p2 <- ggplot(
    d2,
    aes(
      y = method,
      x = zeta_sol
    )
  ) +
    geom_boxplot(
      fill = "cyan",
      linewidth = 0.5
    ) +
    scale_y_discrete(
      limits = rev(method_order),
      labels = rev(method_labels)
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_minimal(
      base_size = 12
    )


  # ===========================================================================
  # Half-length: r1*
  # ===========================================================================

  d3 <- zeta_all[
    zeta_all$censoring == censor &
      zeta_all$r_char == "r1star",
  ]

  p3 <- ggplot(
    d3,
    aes(
      y = method,
      x = zeta_sol
    )
  ) +
    geom_boxplot(
      fill = "purple",
      linewidth = 0.5
    ) +
    scale_y_discrete(
      limits = rev(method_order),
      labels = rev(method_labels)
    ) +
    labs(
      x = "estimated half-length",
      y = NULL
    ) +
    theme_minimal(
      base_size = 12
    )


  left_panel <- p1 / p2 / p3


  # ===========================================================================
  # Coverage panel
  # ===========================================================================

  dcvg <- pred_cvg_all[
    pred_cvg_all$censoring == censor,
  ]

  levels_vert <- c(
    paste0(
      rev(method_order),
      "_r1star"
    ),
    "gap1",
    paste0(
      rev(method_order),
      "_r2"
    ),
    "gap2",
    paste0(
      rev(method_order),
      "_r1"
    )
  )

  dcvg$cat <- paste0(
    dcvg$method,
    "_",
    dcvg$r_char
  )

  dcvg$cat <- factor(
    dcvg$cat,
    levels = levels_vert
  )

  dcvg$r_char <- factor(
    dcvg$r_char,
    levels = c(
      "r1",
      "r2",
      "r1star"
    )
  )

  right_panel <- ggplot(
    dcvg,
    aes(
      x = pred_cvg,
      y = cat,
      fill = r_char
    )
  ) +
    geom_boxplot(
      linewidth = 0.5
    ) +
    geom_vline(
      xintercept = 1 - alpha,
      color = "red"
    ) +
    scale_fill_manual(
      values = c(
        r1 = "magenta",
        r2 = "cyan",
        r1star = "purple"
      ),
      breaks = c(
        "r1",
        "r2",
        "r1star"
      ),
      labels = c(
        expression(r[1]),
        expression(r[2]),
        expression(r[1] * "*")
      ),
      name = NULL
    ) +
    scale_y_discrete(
      breaks = levels_vert,
      labels = c(
        rev(method_labels),
        "",
        rev(method_labels),
        "",
        rev(method_labels)
      ),
      drop = FALSE
    ) +
    labs(
      x = "coverage rate",
      y = NULL
    ) +
    theme_minimal(
      base_size = 12
    ) +
    theme(
      legend.position = "right"
    )


  # ===========================================================================
  # Combined figure
  # ===========================================================================

  final_plot <- left_panel | right_panel

  ggsave(
    filename = file.path(
      RESULTS_DIR,
      paste0(
        "boxplot_",
        censor,
        ".pdf"
      )
    ),
    plot = final_plot,
    width = 8,
    height = 6
  )
}
