# Simulation study

This directory contains the code used for the simulation study accompanying **prescco**. The simulations are self-contained and do not require installing the `prescco` package.

## Files

```text
sim/
├── core.R             # data generation and shared numerical functions
├── beta.R             # SPARCC and complete-case beta estimation
├── PRESCCO.R          # PRESCCO half-length estimation
├── conformal.R        # split conformal, full conformal, and jackknife+
├── coverage.R         # prediction coverage for PRESCCO
├── visualization.R    # summary table and figures
├── Rdata_files/       # simulation outputs used by the scripts
└── results/           # summary table and figures
```

The simulation considers five censoring levels (`low`, `lowmid`, `mid`, `highmid`, and `high`), three residual definitions (`r1`, `r2`, and `r1*`), and three nuisance specifications for PRESCCO (correct, `mis1`, and `mis2`).

The censoring levels correspond to the following values of the censoring-model mean:

| Level   | `alpha2` | Approximate censoring rate |
| ------- | -------: | -------------------------: |
| low     |        2 |                     20–30% |
| lowmid  |        1 |                     30–40% |
| mid     |        0 |                     45–55% |
| highmid |       -1 |                     60–70% |
| high    |       -2 |                     70–80% |

For `mis1`, the mean of the working model for the censored covariate is set to -2. For `mis2`, the censoring-model mean is set to 0 when the true value is nonzero and to 2 when the true value is 0. The `r1*` residual also uses -2 in place of the true center.

## Running the simulations

Run the scripts from the package root.

```sh
Rscript sim/beta.R
Rscript sim/PRESCCO.R
Rscript sim/conformal.R
Rscript sim/coverage.R
Rscript sim/visualization.R
```

`beta.R` should be run before `PRESCCO.R`, and `PRESCCO.R` should be run before `coverage.R`. The conformal simulations are computed independently in `conformal.R`.

All `.Rdata` outputs are written to `sim/Rdata_files/`. The supplied files in this directory can also be used directly to reproduce the simulation summary and figures without rerunning the simulations:

```sh
Rscript sim/visualization.R
```

The resulting summary table and figures are written to `sim/results/`.

## Reproducibility

Simulation replicate `k` is generated using seed `k`. Prediction coverage is evaluated on an independent test set generated using seed `k + 98765`, with test-set size 10,000. The same test-set construction is used for PRESCCO and all three conformal methods.

The main simulation settings are `n = 1000`, 1,000 Monte Carlo replicates, and miscoverage level `alpha = 0.1`. Numerical settings used by PRESCCO are specified directly in `PRESCCO.R`.

Required R packages are `truncnorm`, `MASS`, and `nleqslv`, together with `ggplot2` and `patchwork` for visualization.
