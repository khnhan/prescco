# prescco

<u>PRE</u>diction with <u>S</u>emiparametric Efficiency under a
right-<u>C</u>ensored <u>CO</u>variate.

Prediction intervals for an outcome when one of the covariates is a
right-censored time-to-event variable.

## The problem

You want to predict an outcome $Y$ from a covariate $X$ that is a time to some
event — time to progression, time to failure, age at onset. For some subjects
that event has not happened by the end of follow-up, so $X$ is never observed.
What you observe instead is

$$
W = \min(X, C), \qquad \Delta = I(X \le C)
$$

the earlier of the event time and the censoring time, together with an
indicator of which one it was. When $\Delta = 1$ the event happened and
$W = X$. When $\Delta = 0$ all you know is that $X > W$.

Two obvious approaches both fail:

* **Regress $Y$ on $W$.** For censored subjects $W$ understates $X$, so the
  fit is attenuated and the intervals end up centered in the wrong place.
* **Drop the censored subjects.** Unbiased under independent censoring, but it
  throws away every observation with $\Delta = 0$ — often most of the sample
  under heavy censoring — and the resulting intervals are needlessly wide.

This package keeps the censored observations and uses what they do tell you:
that $X$ exceeded $W$.

A prediction interval here has two pieces: a **center** $m$, which predicts $Y$
from the observed data, and a **half-length** $\zeta$, chosen so that

$$
P\left( \lvert Y - m \rvert \le \zeta \right) = 1 - \alpha
$$

The interval is $[m - \zeta,\ m + \zeta]$. Estimating $m$ needs the outcome
coefficients $\beta$; estimating $\zeta$ needs the distribution of the residual
$\lvert Y - m \rvert$. Both are complicated by $X$ being only partly observed.

## The model

The outcome is linear in a basis over $X$ and any fully observed covariates
$Z$, split into a continuous block $Z_c$ and a discrete block $Z_d$:

$$
\phi_{XZ}(x, z_c, z_d) = (1,\ x,\ z_c,\ z_d,\ x z_c,\ x z_d)
$$

$$
Y \mid X, Z \ \sim\ N\left( \phi_{XZ}(X, Z)^\top \beta,\ \sigma^2 \right)
$$

With no covariates this reduces to
$Y \mid X \sim N(\beta_1 + \beta_2 X, \sigma^2)$ and $\beta$ has length 2.
`phi_xz()` is exported precisely so you can read off which coefficient is which
when $\beta$ is longer.

Two **working nuisance models** describe the event and censoring times. Both
are truncated normal on $[w_{\min}, w_{\max}]$ with means linear in $Z$:

$$
\eta_1: \quad X \mid Z \ \sim\ \mathrm{TruncNormal}\left( (1, z_c, z_d)^\top \alpha_1,\ \tau_1^2 \right)
$$

$$
\eta_2: \quad C \mid Z \ \sim\ \mathrm{TruncNormal}\left( (1, z_c, z_d)^\top \alpha_2,\ \tau_2^2 \right)
$$

They are *working* models: $f_{X \mid Z}$ and $f_{C \mid Z}$ are generally
unknown, so what the methods actually use is a specification of them that may
or may not be right. $\eta_1$ is written $f^{*}_{X \mid Z}$ below when the
distinction matters. The simulation study evaluates each method under a correct
specification and under two deliberate misspecifications.

## Estimating $\beta$

Two estimators, deliberately named alike:

* **`find_beta_cc()`** — complete-case estimator. Regresses $Y$ on the
  design built from $(W, Z)$ using only the $\Delta = 1$ rows, where $W = X$.
  Simple and consistent under independent censoring, but discards data.

* **`find_beta_sparcc()`** — SPARCC estimator. It solves
  an estimating equation built from the model's efficient influence function,
  which uses *every* observation: for a censored subject it integrates the
  score over the conditional distribution of $X$ given $X > W$ and $Z$ implied
  by $\eta_1$. This recovers information the complete-case fit discards, at the
  cost of needing the working models and a numerical solve.

## Interval centers and residuals

Choosing the center is not straightforward here. The natural choice is the
conditional mean of the outcome given the covariates,

$$
m_0(X, Z, \beta) = E\left( Y \mid X, Z, \beta \right) = \phi_{XZ}(X, Z)^\top \beta
$$

but $X$ is not always observed, so $m_0$ is unusable as a center. The center
has to be a function of the observed data $(W, \Delta, Z)$. Three such centers
are available.

**$m_1(W, \Delta, Z, \beta)$** — the conditional mean given what was actually
observed:

$$
m_1 = E\left( Y \mid W, \Delta, Z, \beta, f_{X \mid Z} \right)
$$

Under noninformative censoring, $C$ independent of $(X, Y)$ given $Z$, this
splits into the two cases:

$$
m_1 =
\begin{cases}
m_0(W, Z, \beta), & \Delta = 1 \\
\dfrac{E\left[\, I(X > W)\, m_0(X, Z, \beta) \mid W, Z, f_{X \mid Z} \,\right]}
      {E\left[\, I(X > W) \mid W, Z, f_{X \mid Z} \,\right]}, & \Delta = 0
\end{cases}
$$

That is, when the event was seen ($X$ was observed, $X = W$) you use $m_0$ at
the observed $X$; when it was censored you average $m_0(X, Z, \beta)$ over the
part of the distribution of $X$ that is still possible, namely $X > W$.

**$m_1^{*}(W, \Delta, Z, \beta)$** — the same thing computed under a *working*
model. $m_1$ requires the true $f_{X \mid Z}$, which is generally unknown.
Replacing it with a possibly misspecified working model $f^{*}_{X \mid Z}$
gives

$$
m_1^* = E\left( Y \mid W, \Delta, Z, \beta, f^*_{X \mid Z} \right)
$$

This is the center you can actually compute in practice, and it is why the
package treats the $X \mid Z$ model as a *working* model throughout.

**$m_2(W, \Delta, Z, \beta)$** — sidesteps the integration entirely:

$$
m_2 = m_0(W, Z, \beta)
$$

It uses $W$ in place of $X$ even when $\Delta = 0$. It is biased when censoring
is present, but it is cheap and it does not depend on the $X \mid Z$ model at
all, so nothing about it can be misspecified.

All three depend on the unknown $\beta$, so the center is written generically
as $m(W, \Delta, Z, \beta)$.

Each center gives a residual, and the package labels them:

| Center | Residual | Depends on $f_{X \mid Z}$ |
|---|---|---|
| $m_1(W, \Delta, Z, \beta)$ | $r_1 = \lvert Y - m_1 \rvert$ | yes, the true model |
| $m_1^{*}(W, \Delta, Z, \beta)$ | $r_1^{*} = \lvert Y - m_1^{*} \rvert$ | yes, a working model |
| $m_2(W, \Delta, Z, \beta)$ | $r_2 = \lvert Y - m_2 \rvert$ | no |

Other centers are possible; these three are the natural ones in this setting.

Throughout, math names the estimands and backticks name the code: $r_1$ is the
residual, `r1()` is the function that computes it.

In code, $m_0$, $m_1$, and $m_2$ are internal — you reach them through `r1()`
and `r2()` and through the interval functions. $m_2$ has no separate function,
since it is just $m_0$ evaluated at $W$; that is what `r2()` computes.
$r_1^{*}$ has no separate function either: it is `r1()` evaluated under a
working model, which you select through the `alpha1_star_r` (and optionally
`tau1_r`) argument. Passing the true $X \mid Z$ parameters there gives $r_1$;
passing working-model values gives $r_1^{*}$.

## Estimating the half-length

Two families, with different guarantees.

**PRESCCO.** `prescco_prediction_interval()` treats
$P(\lvert Y - m \rvert \le \zeta) = 1 - \alpha$ as an estimating equation in
$\zeta$ and solves it using the efficient influence function for that
probability, again integrating over the unobserved $X$ for censored subjects.
Efficient when the working models are right; its accuracy degrades as they get
worse, which is what the misspecification settings are designed to probe.

**Conformal.** Three methods that calibrate $\zeta$ from observed residual
quantiles rather than from a model, and so stay valid in finite samples whether
or not the working models are right:

* `split_conformal_prediction_interval()` — fits $\beta$ on one half of the
  data and takes the $1 - \alpha$ empirical quantile of the residuals on the
  other half. One fit, so it is fast, but it pays for the split in width.
* `full_conformal_prediction_interval()` — for each candidate outcome value,
  refits with that point appended and asks whether its residual is extreme
  among the augmented sample. No data splitting, far more computation, and the
  interval differs at every test point.
* `jackknife_plus_prediction_interval()` — builds the interval from
  leave-one-out fits and their residuals; a middle ground between the two.

**Coverage rate.** Given a fitted $\zeta$, `prediction_coverage_rate()` reports
the fraction of a test set whose residual falls inside it — the empirical check
on whether the nominal $1 - \alpha$ was achieved.

## Installation

```r
# install.packages("remotes")
remotes::install_github("<your-account>/prescco")
```

Dependencies: `MASS`, `truncnorm`, `nleqslv`, `stats`.

## Usage

Every function takes data vectors plus optional covariates. Pass
`z_c_data` / `z_d_data` for the covariate case, or omit them for the
no-covariate case — nothing else changes.

### Fit the working models

```r
a1 <- find_alpha1_MLE(w_data, delta_data, z_c_data, z_d_data, w_min = 0, w_max = 8)
a2 <- find_alpha2_MLE(w_data, delta_data, z_c_data, z_d_data, w_min = 0, w_max = 8)

alpha1_star <- a1[-length(a1)];  tau1 <- a1[length(a1)]
alpha2_star <- a2[-length(a2)];  tau2 <- a2[length(a2)]
```

Each returns the mean coefficients followed by the standard deviation.

### Estimate $\beta$

```r
cc     <- find_beta_cc(y_data, w_data, delta_data, z_c_data, z_d_data)
sparcc <- find_beta_sparcc(y_data, w_data, delta_data,
                           z_c_data, z_d_data,
                           alpha1_star, alpha2_star, tau1, tau2,
                           w_min = 0, w_max = 8)

cc$beta_cc
sparcc$beta_hat
```

### PRESCCO interval

Runs end to end — fits the nuisance models, fits $\beta$ by SPARCC, builds the
influence-function arrays, solves for $\zeta$. Supply any piece to skip its
step.

```r
fit <- PRESCCO_prediction_interval(
  y_data, w_data, delta_data, z_c_data, z_d_data,
  residual = r1, alpha = 0.1, w_min = 0, w_max = 8
)
fit$zeta            # half-length, named by the residual
fit$coverage_rate   # NA unless test data was supplied
```

$r_1^{*}$ is the same call with working-model values in `alpha1_star_r` (and
`tau1_r`) instead of the true $X \mid Z$ parameters — label it so the result
records which residual it is:

```r
fit_star <- PRESCCO_prediction_interval(
  y_data, w_data, delta_data,
  residual = r1, residual_name = "r1star",
  alpha1_star_r = wrong_alpha1, alpha = 0.1, w_min = 0, w_max = 8
)
```

### Conformal intervals

```r
scp <- split_conformal_prediction_interval(
  y_data, w_data, delta_data, z_c_data, z_d_data,
  alpha1 = alpha1_star, alpha2 = alpha2_star, alpha1_star_r = alpha1_star,
  alpha = 0.1, tau1 = tau1, tau2 = tau2, tau1_r = tau1,
  w_min = 0, w_max = 8
)

fcp <- full_conformal_prediction_interval(
  alpha = 0.1,
  y_data, w_data, delta_data, z_c_data, z_d_data,
  test_y_data, test_w_data, test_delta_data, test_z_c_data, test_z_d_data,
  alpha1 = alpha1_star, tau1 = tau1,
  alpha1_star_r = alpha1_star, tau1_r = tau1,
  w_min = 0, w_max = 8
)
```

`jackknife_plus_prediction_interval()` takes the same arguments as
`full_conformal_prediction_interval()`.

### One return shape

All four interval methods return the same list:

```r
list(
  method        = "split_conformal",          # which method produced this
  alpha         = 0.1,                        # miscoverage level
  residual      = c("r1", "r2", "r1star"),    # residuals reported
  zeta          = c(r1 = ..., r2 = ..., r1star = ...),   # half-length
  coverage_rate = c(r1 = ..., r2 = ..., r1star = ...),   # empirical rate
  ...                                         # method-specific extras
)
```

| Function | Test set | Residuals per call |
|---|---|---|
| `PRESCCO_prediction_interval()` | optional | one (the `residual` argument) |
| `split_conformal_prediction_interval()` | optional | all three |
| `full_conformal_prediction_interval()` | **required** | all three |
| `jackknife_plus_prediction_interval()` | **required** | all three |

The PRESCCO and split conformal prediction
methods produce a single half-length from training data alone, so a test set is
optional: supply one and `coverage_rate` is filled in, omit it and the entries
are `NA`. Full conformal and jackknife+ construct a *different* interval at
every test point, so there is no half-length to report without test data; what
they return is the average over the test set. To score a fitted half-length
later, call `prediction_coverage_rate()` directly.

## Exported functions

Twelve. Everything else is internal.

| Function | Purpose |
|---|---|
| `find_alpha1_MLE()` | MLE for the working $X \mid Z$ model |
| `find_alpha2_MLE()` | MLE for the working $C \mid Z$ model |
| `find_beta_cc()` | Complete-case beta estimator |
| `find_beta_sparcc()` | SPARCC beta estimator |
| `PRESCCO_prediction_interval()` | PRESCCO prediction interval |
| `split_conformal_prediction_interval()` | Split conformal prediction interval |
| `full_conformal_prediction_interval()` | Full conformal prediction interval |
| `jackknife_plus_prediction_interval()` | Jackknife+ prediction interval |
| `prediction_coverage_rate()` | Coverage of a fitted half-length on a test set |
| `phi_xz()` | Covariate basis fixing the layout of $\beta$ |
| `r1()`, `r2()` | Residuals, passed to the two functions that take one |

## Package layout

```
prescco/
├── DESCRIPTION                     # metadata and dependencies
├── NAMESPACE                       # exports / imports (regenerate with devtools::document())
├── LICENSE, LICENSE.md             # MIT license
├── README.md
├── .Rbuildignore                   # keeps sim/ out of the installed package
├── R/                              # the installed library
│   ├── prescco-package.R               # package doc; central @importFrom declarations
│   ├── predict_helpers.R               # phi_xz basis, m0/m1, r1/r2 (Z-optional)
│   ├── nuisance.R                      # find_alpha1_MLE, find_alpha2_MLE
│   ├── eif_internals.R                 # efficient influence function helpers
│   ├── beta_fit.R                      # find_beta_cc, find_beta_sparcc
│   ├── predict_PRESCCO.R               # PRESCCO_prediction_interval
│   └── predict_conformal.R             # the three conformal methods
├── tests/                          # testthat suite
├── sim/                            # simulation study (not installed)
└── inst/results/figures/           # reference figures, one per censoring level
```

## Testing

```r
devtools::test()          # Ctrl/Cmd+Shift+T in RStudio
devtools::check()         # runs the suite as part of R CMD check
```

Tests across five files. Most run in seconds; the end-to-end paths fit nuisance
models on every call and are gated behind an environment variable:

```r
Sys.setenv(PRESCCO_SLOW_TESTS = "true")
devtools::test()
```

The suite checks properties that should hold regardless of implementation:
$m_1$ collapses to $m_0(W, Z, \beta)$ when $\Delta = 1$; the vectorized
internals agree with their scalar counterparts elementwise across every
combination of covariate blocks; `find_beta_cc` reproduces `lm()` on the
complete cases; `prediction_coverage_rate` is monotone in $\zeta$;
split conformal prediction half-lengths widen as $\alpha$ shrinks; and all four methods
return the same shape.

## Simulation study

The simulation study lives in `sim/`. It is excluded from the
package build, has no covariates, prefixes its helpers `sim_`, and is run by
sourcing `sim/run.R` rather than through `library(prescco)` — so it stays
insulated from the library API and can run without installing the package.

It crosses five censoring levels with three residuals ($r_1$, $r_2$,
$r_1^{*}$), running the PRESCCO method under three nuisance
specifications (correct, `mis1`, `mis2`) alongside all three conformal methods,
then compares half-lengths and coverage rates.

```sh
cd sim
Rscript run.R all .
```

See `sim/README.md` for the stage-by-stage breakdown, the censoring-level
parameters, how to run a subset, and the reproducibility guarantees. Reference
figures from a full run are in `inst/results/figures/`.

## Author

Kihyun Han (khnhan@psu.edu)

## License

MIT (c) 2026 Kihyun Han. See `LICENSE.md`.
