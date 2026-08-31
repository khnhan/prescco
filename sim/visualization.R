# =============================================================================
#  Visualization / post-processing for prescco
#  ---------------------------------------------------------------------------
#  Aggregates the per-method .Rdata outputs (half-length zeta_sol and
#  coverage rate pred_cvg) into:
#    1. a summary table (mean / bias / empirical SD of zeta, mean / SD of the
#       empirical coverage rate), written to CSV, and
#    2. the combined "transposed" boxplot figure per censoring level: a left
#       column of three stacked half-length panels (r1, r2, r1*) and a right
#       panel of empirical coverage rates for all method x residual combinations.
#
#  Inference/variance summaries (estimated SD, zeta coverage rate) from the original
#  script are omitted, matching the rest of the package.
#
#  Requires ggplot2 and patchwork. The six methods are, in
#  order: PRESCCO (eta1,eta2) / (eta1*,eta2) / (eta1,eta2*), split
#  conformal, full conformal, jackknife+.
# =============================================================================

# Fixed method order and display labels shared by table and plots.
.method_keys   <- c("true", "mis1", "mis2", "scp", "fcp", "jackknife")
.residual_keys <- c("r1", "r2", "r1star")

# Plain-text method labels (used in the CSV / non-expression contexts).
.method_labels_plain <- c(
  true = "(eta1,eta2)", mis1 = "(eta1*,eta2)", mis2 = "(eta1,eta2*)",
  scp = "Split CP", fcp = "Full CP", jackknife = "Jackknife+")

# ---- Filename helpers --------------------------------------------------------
# Resolve the on-disk file names produced by the earlier stages.

.zeta_file <- function(r_char, method, censor, alpha = 0.1) {
  rlab  <- if (r_char == "r1star") "r1" else r_char
  infix <- if (r_char == "r1star") "_r1star" else ""
  switch(method,
    true      = paste0("zeta_param12_", rlab, "_", alpha, "_", censor, infix, ".Rdata"),
    mis1      = paste0("zeta_param12_", rlab, "_", alpha, "_", censor, infix, "_mis1.Rdata"),
    mis2      = paste0("zeta_param12_", rlab, "_", alpha, "_", censor, infix, "_mis2.Rdata"),
    scp       = paste0("zeta_param12_", rlab, "_", alpha, "_", censor, infix, "_cp.Rdata"),
    fcp       = paste0("zeta_", rlab, "_", alpha, "_", censor, infix, "_fcp.Rdata"),
    jackknife = paste0("zeta_", rlab, "_", alpha, "_", censor, infix, "_jackknife.Rdata"))
}

.pred_cvg_file <- function(r_char, method, censor, alpha = 0.1) {
  rlab  <- if (r_char == "r1star") "r1" else r_char
  infix <- if (r_char == "r1star") "_r1star" else ""
  switch(method,
    true      = paste0("pred_cvg_", rlab, "_", alpha, "_", censor, infix, ".Rdata"),
    mis1      = paste0("pred_cvg_", rlab, "_", alpha, "_", censor, infix, "_mis1.Rdata"),
    mis2      = paste0("pred_cvg_", rlab, "_", alpha, "_", censor, infix, "_mis2.Rdata"),
    scp       = paste0("pred_cvg_", rlab, "_", alpha, "_", censor, infix, "_cp.Rdata"),
    fcp       = paste0("pred_cvg_", rlab, "_", alpha, "_", censor, infix, "_fcp.Rdata"),
    jackknife = paste0("pred_cvg_", rlab, "_", alpha, "_", censor, infix, "_jackknife.Rdata"))
}

.load_var <- function(path, want) {
  if (!file.exists(path)) stop("missing input file: ", path)
  env <- new.env(); load(path, envir = env)
  if (!exists(want, envir = env)) stop("`", want, "` not found in ", path)
  get(want, envir = env)
}

#' Collect half-length and coverage rate draws for one censoring level
#'
#' Reads every method x residual pair for one censoring level and returns tidy
#' long data frames of the half-length and coverage rate draws.
#'
#' @param censor Censoring level tag.
#' @param in_dir Directory holding the \code{.Rdata} inputs.
#' @param alpha Miscoverage level used in the filenames.
#' @return A list with \code{zeta} and \code{cvg} data frames, each with columns
#'   \code{r_char}, \code{method}, and the drawn value.
#' @export
collect_sim_draws <- function(censor, in_dir = ".", alpha = 0.1) {
  zeta <- list(); cvg <- list()
  for (r_char in .residual_keys) {
    for (mth in .method_keys) {
      zs <- .load_var(file.path(in_dir, .zeta_file(r_char, mth, censor, alpha)), "zeta_sol")
      pc <- .load_var(file.path(in_dir, .pred_cvg_file(r_char, mth, censor, alpha)), "pred_cvg")
      zeta[[length(zeta) + 1]] <- data.frame(r_char = r_char, method = mth, zeta_sol = zs)
      cvg[[length(cvg) + 1]]   <- data.frame(r_char = r_char, method = mth, pred_cvg = pc)
    }
  }
  list(zeta = do.call(rbind, zeta), cvg = do.call(rbind, cvg))
}

#' Summarize one censoring level into a table of means and SDs
#'
#' @param draws Output of \code{collect_sim_draws}.
#' @param censor Censoring level tag (added as a column).
#' @param zeta_true Optional named vector of true half-lengths by residual
#'   (\code{r1}, \code{r2}, \code{r1star}); when supplied, a \code{bias} column
#'   is added.
#' @return A data frame with one row per residual x method.
#' @export
summarize_sim_level <- function(draws, censor, zeta_true = NULL) {
  rows <- list()
  for (r_char in .residual_keys) {
    for (mth in .method_keys) {
      z <- draws$zeta$zeta_sol[draws$zeta$r_char == r_char & draws$zeta$method == mth]
      p <- draws$cvg$pred_cvg[draws$cvg$r_char == r_char & draws$cvg$method == mth]
      zmean <- mean(z)
      row <- data.frame(
        censoring     = censor,
        r_char        = r_char,
        method        = mth,
        zeta_mean     = zmean,
        emp_sd        = sd(z),
        mean_pred_cvg = mean(p),
        sd_pred_cvg   = sd(p),
        stringsAsFactors = FALSE)
      if (!is.null(zeta_true)) row$bias <- zmean - zeta_true[[r_char]]
      rows[[length(rows) + 1]] <- row
    }
  }
  do.call(rbind, rows)
}

# ---- Combined transposed plot ------------------------------------------------

# Expression labels for the six methods (used on plot axes).
.method_labels_expr <- function() {
  expression(
    "(" * eta[1] * "," * eta[2] * " )",
    "(" * eta[1] * "*" * "," * eta[2] * ")",
    "(" * eta[1] * "," * eta[2]^"\u2605" * ")",
    "Split CP", "Full CP", "Jackknife+")
}

#' Build the combined transposed boxplot for one censoring level
#'
#' Left column: three stacked half-length panels (r1, r2, r1*). Right column:
#' one tall panel of empirical coverage rates for all 18 method x residual
#' combinations (with small gaps between residual blocks) and a vertical
#' reference line at \code{1 - alpha}.
#'
#' @param draws Output of \code{collect_sim_draws}.
#' @param alpha Miscoverage level (reference line at \code{1 - alpha}).
#' @param fills Named fill colors for \code{r1}, \code{r2}, \code{r1star}.
#' @return A patchwork object combining the panels.
#' @export
plot_combined_transposed <- function(draws, alpha = 0.1,
                                     fills = c(r1 = "magenta", r2 = "cyan",
                                               r1star = "purple")) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("plot_combined_transposed requires the 'ggplot2' and 'patchwork' packages.")
  }
  labels <- .method_labels_expr()
  rev_methods <- rev(.method_keys)

  # --- left: three stacked half-length panels ---
  half_panel <- function(rc, fill, xlab = NULL) {
    d <- draws$zeta[draws$zeta$r_char == rc, ]
    ggplot2::ggplot(d, ggplot2::aes(y = ggplot2::.data$method, x = ggplot2::.data$zeta_sol)) +
      ggplot2::geom_boxplot(fill = fill, linewidth = 0.5) +
      ggplot2::scale_y_discrete(limits = rev_methods, labels = rev(labels)) +
      ggplot2::labs(x = xlab, y = NULL) +
      ggplot2::theme_minimal(base_size = 12)
  }
  p1 <- half_panel("r1",     fills[["r1"]])
  p2 <- half_panel("r2",     fills[["r2"]])
  p3 <- half_panel("r1star", fills[["r1star"]],
                   xlab = "estimated half-length")
  left_panel <- patchwork::wrap_plots(p1, p2, p3, ncol = 1)

  # --- right: tall coverage rate panel with gaps between residual blocks ---
  levels_vert <- c(
    paste0(rev_methods, "_r1star"), "gap1",
    paste0(rev_methods, "_r2"),     "gap2",
    paste0(rev_methods, "_r1"))
  d456 <- draws$cvg
  d456$cat <- factor(paste0(d456$method, "_", d456$r_char), levels = levels_vert)
  d456$r_char <- factor(d456$r_char, levels = .residual_keys)

  right_panel <- ggplot2::ggplot(
      d456, ggplot2::aes(x = ggplot2::.data$pred_cvg, y = ggplot2::.data$cat, fill = ggplot2::.data$r_char)) +
    ggplot2::geom_boxplot(linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 1 - alpha, color = "red") +
    ggplot2::scale_fill_manual(
      values = fills, breaks = .residual_keys,
      labels = c(expression(r[1]), expression(r[2]), expression(r[1] * "*")),
      name = NULL) +
    ggplot2::scale_y_discrete(
      breaks = levels_vert,
      labels = c(rev(labels), "", rev(labels), "", rev(labels)),
      drop = FALSE) +
    ggplot2::labs(x = "coverage rate", y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right")

  patchwork::wrap_plots(left_panel, right_panel, ncol = 2)
}

# ---- Top-level driver --------------------------------------------------------

#' Reproduce the simulation summary table and combined figures
#'
#' Iterates over all censoring levels, writes the combined transposed boxplot
#' \code{boxplot_combined_transposed_<censor>.png} for each, and writes the
#' rounded summary table to CSV.
#'
#' @param in_dir Directory holding the \code{.Rdata} inputs.
#' @param out_dir Directory for the figures and CSV.
#' @param alpha Miscoverage level.
#' @param zeta_true_file Optional path to a \code{true_zeta.Rdata} (matrix
#'   \code{true_zeta} with residuals in rows and censoring levels in columns,
#'   after two leading id columns) for the \code{bias} column; if \code{NULL},
#'   the bias column is omitted.
#' @param width,height Figure dimensions in inches.
#' @return Invisibly, the combined summary table (also written to CSV).
#' @export
run_visualization <- function(in_dir = ".", out_dir = ".", alpha = 0.1,
                              zeta_true_file = NULL, width = 8, height = 6) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  zeta_true_mat <- NULL
  if (!is.null(zeta_true_file) && file.exists(zeta_true_file)) {
    true_zeta <- .load_var(zeta_true_file, "true_zeta")
    zeta_true_mat <- matrix(apply(true_zeta[, -(1:2)], 1, median),
                            ncol = length(.censor_alpha2))
    rownames(zeta_true_mat) <- .residual_keys
    colnames(zeta_true_mat) <- names(.censor_alpha2)
  }

  result <- list()
  levels <- names(.censor_alpha2)
  for (i in seq_along(levels)) {
    censor <- levels[[i]]
    message(sprintf("[viz %d/%d] %s", i, length(levels), censor))
    draws <- collect_sim_draws(censor, in_dir = in_dir, alpha = alpha)

    zt <- if (!is.null(zeta_true_mat))
      stats::setNames(zeta_true_mat[, censor], .residual_keys) else NULL
    result[[length(result) + 1]] <- summarize_sim_level(draws, censor, zeta_true = zt)

    p <- plot_combined_transposed(draws, alpha = alpha)
    ggplot2::ggsave(
      filename = file.path(out_dir,
                           paste0("boxplot_combined_transposed_", censor, ".png")),
      plot = p, width = width, height = height)
  }

  result <- do.call(rbind, result)
  num_cols <- vapply(result, is.numeric, logical(1))
  result[num_cols] <- round(result[num_cols], 3)
  utils::write.csv(result, file = file.path(out_dir, "simulation_result.csv"),
                   row.names = FALSE)
  invisible(result)
}
