# Simulation study (companion to censcovpred)

The simulation study behind the package. **Self-contained** and
**no-covariate**: nothing here needs `library(censcovpred)`. Shared helpers are
prefixed `sim_` so they never collide with the general library in `../R/`.

## Layout

```
sim/
├── run.R              # entry point: sources the files below, dispatches a stage
├── core.R             # DGP, EIF/SPARCC machinery, centers and residuals
│                      #   (scalar + vectorized), b-solvers, shared constants
├── beta.R             # SPARCC and complete-case beta: estimators, configs, drivers
├── semiparametric.R   # zeta solver, 45-run config, drivers
├── conformal.R        # split CP / full CP / jackknife+, config, drivers
├── coverage.R         # coverage rate for the semiparametric method
└── visualization.R    # summary table + figures, all six methods (cross-cutting)
```

Each of the three estimation stages is one file holding its estimator, its
configuration table, and its drivers. All three conformal methods evaluate
their own coverage rate, so `conformal.R` is fully self-contained.
`coverage.R` exists only for the semiparametric method, which fits a
half-length without evaluating one. `visualization.R` is separate because it
spans all six methods.

`run.R` sources the other files in dependency order (`core.R` first). Sourcing
it defines every function and runs nothing; running it as a script also dispatches
one stage.

## The design

Full factorial, everything crossed:

* **5 censoring levels** — `low`, `lowmid`, `mid`, `highmid`, `high`
* **3 residuals** — `r1`, `r2`, `r1*`
* **3 nuisance specifications** (semiparametric) — correct, `mis1`, `mis2`
* **3 conformal methods** — split CP, full CP, jackknife+

| Stage | Driver | Output | Count |
|-------|--------|--------|-------|
| SPARCC beta | `run_all_beta_scenarios()` | `beta_param12_<censor>[_mis1\|_mis2].Rdata` | 15 |
| CC beta | `run_all_cc_scenarios()` | `beta_cc_<censor>.Rdata` | 5 |
| Semiparametric | `run_all_scenarios()` | `zeta_param12_*.Rdata` | 45 |
| Conformal | `run_all_conformal_scenarios()` | `zeta_*` + `pred_cvg_*` `.Rdata` | 15 |
| Coverage rate | `run_all_pred_cvg()` | `pred_cvg_*.Rdata` | 45 |
| Visualization | `run_visualization()` | `simulation_result.csv`, 5 figures | — |

Only one dependency between stages: the semiparametric method needs the SPARCC
beta estimates. The conformal methods are self-contained.

Beta is estimated per (censoring level, nuisance specification) and does not
depend on the residual, so the `r1`, `r2`, and `r1*` runs at a given cell all
read the same beta file. CC does not depend on the specification either, hence
one file per censoring level.

Censoring level maps to the true censoring-time mean `alpha2` (gamma):

| Level | `alpha2` | Approx. censoring | `mis2` uses `gamma*` |
|-------|----------|-------------------|----------------------|
| low | 2 | 20–30% | 0 |
| lowmid | 1 | 30–40% | 0 |
| mid | 0 | 45–55% | 2 |
| highmid | -1 | 60–70% | 0 |
| high | -2 | 70–80% | 0 |

`mis1` sets the `eta1` mean to -2, which is also the misspecified center used
by the `r1*` residual.

## Running

From the command line, inside `sim/`:

```sh
Rscript run.R beta       beta          # SPARCC beta   (prerequisite for semi)
Rscript run.R cc         beta          # CC beta
Rscript run.R semi       beta  zeta    # semiparametric intervals
Rscript run.R conformal  zeta          # split CP, full CP, jackknife+ (+ their coverage rates)
Rscript run.R coverage   beta  zeta    # coverage rate for the semiparametric method
Rscript run.R viz        zeta  results # summary table + figures
```

or everything in order:

```sh
Rscript run.R all .
```

From an R session:

```r
source("sim/run.R")                                    # defines the functions only
run_all_beta_scenarios(out_dir = "beta")
run_all_cc_scenarios(out_dir = "beta")
run_all_scenarios(in_dir = "beta", out_dir = "zeta")
run_all_conformal_scenarios(out_dir = "zeta")
run_all_pred_cvg(beta_dir = "beta", in_dir = "zeta", out_dir = "zeta")
run_visualization(in_dir = "zeta", out_dir = "results")
```

`run_all_pred_cvg()` reads the SPARCC beta files from `beta_dir` and the fitted
half-lengths from `in_dir`.

Dependencies: `truncnorm`, `MASS`, `nleqslv`, plus `ggplot2` and `patchwork`
for the `viz` stage.

## Running a subset or a smoke test

```r
scen <- build_scenarios()

# the moderate-to-high, r1, correct-nuisance run, with few replicates
idx <- which(scen$censor_tag    == "highmid" &
             scen$residual_name == "r1" &
             scen$spec_tag      == "")

run_all_scenarios(in_dir = "beta", out_dir = "zeta",
                  which = idx, n_rep = 20, verbose = TRUE)
```

`run_scenario()` exposes the numerical controls (`tt`, `tt2`, `tt3`,
`seq_length`, `m`, `n_rep`); the defaults are the ones the reference
figures were produced with.
Changing `M`, `n`, or any of the grid sizes changes every downstream number,
since the semiparametric interpolation grid is built from the mean and SD of
the whole beta file.

## Reproducibility

Every replicate is seed-self-contained: `data_generating(k, ...)` calls
`set.seed(k)` internally, and `get_beta_param_param12()` reseeds before drawing
its `nleqslv` starting perturbation. Results therefore do not depend on loop
order, on how many scenarios are run, or on anything earlier in the session.

All six methods are evaluated on the **same test set** per replicate: seed
`k + seed_offset` (default 98765), size `N` (default 10000), generated from the
same DGP parameters. The three conformal methods compute their coverage rate
inside `conformal.R`; the semiparametric method's is computed in `coverage.R`
from the same generator with the same offset and size.

Reference figures from this study are in `../inst/results/figures/`.
