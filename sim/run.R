# =============================================================================
#  prescco simulation study: single entry point
#  ---------------------------------------------------------------------------
#  Sourcing this file defines every simulation function (in dependency order)
#  in the current session and does nothing else:
#
#      source("sim/run.R")
#      run_all_beta_scenarios(out_dir = "beta")
#      run_all_scenarios(in_dir = "beta", out_dir = "zeta")
#
#  Running it as a script also dispatches one stage:
#
#      Rscript run.R beta       <beta_out>
#      Rscript run.R cc         <beta_out>
#      Rscript run.R PRESCCO       <beta_in>  <zeta_out>
#      Rscript run.R conformal  <zeta_out>
#      Rscript run.R coverage   <beta_in>  <zeta_in_out>
#      Rscript run.R viz        <zeta_in>  <results_out>
#      Rscript run.R all        <workdir>
#
#  The simulation code is self-contained: it does NOT require installing the
#  prescco package. Only `truncnorm`, `MASS`, and `nleqslv` are needed
#  (plus `ggplot2` and `patchwork` for the `viz` stage).
# =============================================================================

## ---- locate this file, so it can be sourced from anywhere ------------------
## Preference order matters: `ofile` (set by source()) is checked first so that
## sourcing this file from *another* Rscript resolves to this file, not to the
## outer script named in `--file=`.

.sim_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)

.sim_locate_self <- function() {
  for (i in seq_len(sys.nframes())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, mustWork = FALSE))
  }
  if (length(.sim_file_arg))
    return(normalizePath(sub("^--file=", "", .sim_file_arg[1]), mustWork = FALSE))
  NA_character_
}

if (!exists("SIM_DIR")) {            # honour a user-set SIM_DIR if present
  .self <- .sim_locate_self()
  SIM_DIR <- if (!is.na(.self)) dirname(.self) else getwd()
}

## ---- external dependencies --------------------------------------------------
for (pkg in c("truncnorm", "MASS", "nleqslv")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' is required for the simulation study.", pkg))
}

## ---- source the simulation files in dependency order -----------------------
.source_order <- c(
  "core.R",           # DGP, EIF/SPARCC machinery, centers/residuals, b-solvers, constants
  "beta.R",           # SPARCC and complete-case beta, configs, drivers
  "PRESCCO.R",        # PRESCCO zeta solver, 45-run config, drivers
  "conformal.R",      # split CP / full CP / jackknife+, config, drivers
  "coverage.R",       # coverage rate for the PRESCCO method
  "visualization.R"   # summary table + combined boxplot figures
)

for (f in .source_order) {
  fp <- file.path(SIM_DIR, f)
  if (!file.exists(fp)) stop("Simulation file not found: ", fp)
  source(fp, local = FALSE)
}

## ---- command-line dispatch --------------------------------------------------
## Script mode only: true when R was launched on THIS file. When this file is
## sourced (interactively or from another script) the functions are defined and
## nothing runs.

.sim_is_script <- local({
  if (!length(.sim_file_arg)) return(FALSE)
  self <- .sim_locate_self()
  invoked <- normalizePath(sub("^--file=", "", .sim_file_arg[1]), mustWork = FALSE)
  !is.na(self) && identical(self, invoked)
})

if (.sim_is_script) {
  .args  <- commandArgs(trailingOnly = TRUE)
  .stage <- if (length(.args) >= 1L) .args[[1]] else ""
  .a <- function(i, default) if (length(.args) >= i) .args[[i]] else default

  .usage <- function() {
    cat("Usage: Rscript run.R <stage> [dirs...]\n\n",
        "  beta       <beta_out>              SPARCC beta   (15 files)\n",
        "  cc         <beta_out>              CC beta       (5 files)\n",
        "  PRESCCO    <beta_in> <zeta_out>    PRESCCO intervals (45 runs)\n",
        "  conformal  <zeta_out>              split CP, full CP, jackknife+ (15 runs)\n",
        "  coverage   <beta_in> <zeta_dir>    coverage rate for the PRESCCO method\n",
        "  viz        <zeta_in> <results_out> summary table + figures\n",
        "  all        <workdir>               every stage above, in order\n",
        sep = "")
  }

  switch(.stage,
    beta      = run_all_beta_scenarios(out_dir = .a(2, "beta")),
    cc        = run_all_cc_scenarios(out_dir = .a(2, "beta")),
    PRESCCO   = run_all_scenarios(in_dir = .a(2, "beta"), out_dir = .a(3, "zeta")),
    conformal = run_all_conformal_scenarios(out_dir = .a(2, "zeta")),
    coverage  = run_all_pred_cvg(beta_dir = .a(2, "beta"), in_dir = .a(3, "zeta"),
                                 out_dir = .a(3, "zeta")),
    viz       = run_visualization(in_dir = .a(2, "zeta"), out_dir = .a(3, "results")),
    all       = {
      wd   <- .a(2, ".")
      beta <- file.path(wd, "beta")
      zeta <- file.path(wd, "zeta")
      res  <- file.path(wd, "results")
      run_all_beta_scenarios(out_dir = beta)
      run_all_cc_scenarios(out_dir = beta)
      run_all_scenarios(in_dir = beta, out_dir = zeta)
      run_all_conformal_scenarios(out_dir = zeta)
      run_all_pred_cvg(in_dir = zeta, beta_dir = beta, out_dir = zeta)
      run_visualization(in_dir = zeta, out_dir = res)
    },
    .usage())
}
