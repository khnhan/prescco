# prescco

<u>PRE</u>diction with <u>S</u>emiparametric Efficiency under a
right-<u>C</u>ensored <u>CO</u>variate.

`prescco` constructs prediction intervals for an outcome when one covariate is
a right-censored time-to-event variable.

## The problem

Let \(X\) be a time-to-event covariate and \(C\) its censoring time. Instead of
observing \(X\) for everyone, we observe

\[
W = \min(X,C), \qquad \Delta = I(X \le C).
\]

When \(\Delta=1\), \(X=W\) is observed. When \(\Delta=0\), we only know that
\(X>W\).

A prediction interval has a center \(m\) and half-length \(\zeta\),

\[
[m-\zeta,\;m+\zeta],
\]

with target coverage

\[
P\{|Y-m|\le \zeta\}=1-\alpha.
\]

The package estimates the outcome model, working nuisance models for \(X\mid Z\)
and \(C\mid Z\), and prediction half-lengths while retaining the censored
observations.

## Model

The outcome model is linear in the basis

\[
\phi_{XZ}(x,z_c,z_d)
=
(1,\;x,\;z_c,\;z_d,\;xz_c,\;xz_d),
\]

so that

\[
Y\mid X,Z
\sim
N\{\phi_{XZ}(X,Z)^\top\beta,\sigma^2\}.
\]

The fully observed covariates \(Z\) may contain a continuous block \(Z_c\), a
discrete block \(Z_d\), both, or neither. With no additional covariates,

\[
Y\mid X \sim N(\beta_1+\beta_2X,\sigma^2).
\]

The working nuisance models are truncated normal on
\([w_{\min},w_{\max}]\):

\[
X\mid Z
\sim
\operatorname{TruncNormal}
\{(1,z_c,z_d)^\top\alpha_1,\tau_1^2\},
\]

\[
C\mid Z
\sim
\operatorname{TruncNormal}
\{(1,z_c,z_d)^\top\alpha_2,\tau_2^2\}.
\]

`phi_xz()` returns the outcome-model basis and fixes the coefficient ordering
used throughout the package.

## Outcome-model coefficient estimation

Two estimators are available.

- `find_beta_cc()` fits the outcome model using only complete cases
  (\(\Delta=1\)), for which \(X=W\) is observed. Under the linear outcome model
  implemented here, this fit is obtained by least squares.

- `find_beta_sparcc()` estimates the outcome-model coefficients with SPARCC,
  using all observations and the working nuisance models.

## Centers and residuals

The complete-data outcome mean is

\[
m_0(X,Z,\beta)
=
E(Y\mid X,Z,\beta)
=
\phi_{XZ}(X,Z)^\top\beta.
\]

Because \(X\) may be censored, the package uses centers based on the observed
data \((W,\Delta,Z)\).

For \(r_1\), the center is

\[
m_1(W,\Delta,Z,\beta)
=
E\{m_0(X,Z,\beta)\mid W,\Delta,Z\},
\]

where censored observations average over the values \(X>W\) allowed by the
model for \(X\mid Z\). The residual is

\[
r_1 = |Y-m_1|.
\]

A working version \(r_1^*\) is obtained by evaluating the same residual with
alternative coefficients and, if desired, an alternative standard deviation
for the working model for \(X\mid Z\). In code this is done by calling `r1()`
with the corresponding `alpha1_star` and `tau1`.

The second residual uses

\[
m_2(W,\Delta,Z,\beta)=m_0(W,Z,\beta),
\qquad
r_2=|Y-m_2|.
\]

`r1()` and `r2()` are vectorized and accept either single observations or
vectors of observations.

## Prediction intervals

### PRESCCO

`PRESCCO_prediction_interval()` treats the target half-length as a parameter of
the observed-data model and solves its efficient estimating equation.

PRESCCO is doubly robust: the half-length estimator is consistent when either
working nuisance model is correctly specified, under the conditions of the
method. When both nuisance models are correctly specified, it attains the
semiparametric efficiency bound.

### Conformal methods

Three comparison methods are provided:

- `split_conformal_prediction_interval()`
- `full_conformal_prediction_interval()`
- `jackknife_plus_prediction_interval()`

The split method estimates the outcome model by SPARCC on a fitting split and
calibrates the half-length on the remaining observations. Full conformal and
jackknife+ use complete-case outcome-model fits internally.

`prediction_coverage_rate()` evaluates the empirical coverage of a fitted
half-length on a supplied test set.

## Installation

```r
# install.packages("remotes")
remotes::install_github("khnhan/prescco")
```

Required packages are `MASS`, `truncnorm`, and `nleqslv`.

## Usage

### Fit the working nuisance models

```r
alpha1_fit <- find_alpha1_MLE(
  w_data, delta_data,
  z_c_data = z_c_data,
  z_d_data = z_d_data,
  w_min = 0,
  w_max = 8
)

alpha2_fit <- find_alpha2_MLE(
  w_data, delta_data,
  z_c_data = z_c_data,
  z_d_data = z_d_data,
  w_min = 0,
  w_max = 8
)

alpha1_star <- alpha1_fit[-length(alpha1_fit)]
tau1 <- alpha1_fit[length(alpha1_fit)]

alpha2_star <- alpha2_fit[-length(alpha2_fit)]
tau2 <- alpha2_fit[length(alpha2_fit)]
```

Each fit returns the mean-model coefficients followed by the estimated standard
deviation.

### Estimate the outcome-model coefficients

```r
fit_cc <- find_beta_cc(
  y_data, w_data, delta_data,
  z_c_data = z_c_data,
  z_d_data = z_d_data
)

fit_sparcc <- find_beta_sparcc(
  y_data, w_data, delta_data,
  z_c_data = z_c_data,
  z_d_data = z_d_data,
  alpha1_star = alpha1_star,
  alpha2_star = alpha2_star,
  tau1 = tau1,
  tau2 = tau2,
  w_min = 0,
  w_max = 8
)

fit_cc$beta_cc
fit_sparcc$beta_hat
```

### PRESCCO

The function can fit the nuisance and outcome models internally, or use
pre-fitted quantities supplied by the user.

```r
fit <- PRESCCO_prediction_interval(
  y_data = y_data,
  w_data = w_data,
  delta_data = delta_data,
  z_c_data = z_c_data,
  z_d_data = z_d_data,
  alpha = 0.1,
  residual = r1,
  w_min = 0,
  w_max = 8
)

fit$zeta
fit$coverage_rate
```

To evaluate \(r_1^*\), supply the coefficients and SD to be used for that
residual and label the result accordingly:

```r
fit_star <- PRESCCO_prediction_interval(
  y_data = y_data,
  w_data = w_data,
  delta_data = delta_data,
  alpha = 0.1,
  alpha1_star_r = alpha1_star_r,
  tau1_r = tau1_r,
  residual = r1,
  residual_name = "r1star",
  w_min = 0,
  w_max = 8
)
```

### Split conformal prediction

```r
fit_split <- split_conformal_prediction_interval(
  y_data = y_data,
  w_data = w_data,
  delta_data = delta_data,
  z_c_data = z_c_data,
  z_d_data = z_d_data,
  alpha = 0.1,
  alpha1 = alpha1_star,
  tau1 = tau1,
  alpha1_star_r = alpha1_star,
  tau1_r = tau1,
  alpha2 = alpha2_star,
  tau2 = tau2,
  w_min = 0,
  w_max = 8
)

fit_split$zeta
```

### Full conformal prediction and jackknife+

Both methods require test observations.

```r
fit_full <- full_conformal_prediction_interval(
  y_data = y_data,
  w_data = w_data,
  delta_data = delta_data,
  z_c_data = z_c_data,
  z_d_data = z_d_data,
  alpha = 0.1,
  alpha1 = alpha1_star,
  tau1 = tau1,
  alpha1_star_r = alpha1_star,
  tau1_r = tau1,
  test_y_data = test_y_data,
  test_w_data = test_w_data,
  test_delta_data = test_delta_data,
  test_z_c_data = test_z_c_data,
  test_z_d_data = test_z_d_data,
  w_min = 0,
  w_max = 8
)
```

`jackknife_plus_prediction_interval()` uses the same arguments as
`full_conformal_prediction_interval()`.

## Return values

All four prediction methods return a list containing at least

```r
method
alpha
residual
zeta
coverage_rate
```

PRESCCO returns one residual per call. The three conformal functions report
\(r_1\), \(r_2\), and \(r_1^*\) together.

For PRESCCO and split conformal prediction, test data are optional. If test data
are omitted, `coverage_rate` is `NA`. Full conformal prediction and jackknife+
require test data and report average half-lengths and empirical coverage over
that test set.

## Exported functions

| Function | Purpose |
|---|---|
| `find_alpha1_MLE()` | Fit the working model for \(X\mid Z\) |
| `find_alpha2_MLE()` | Fit the working model for \(C\mid Z\) |
| `find_beta_cc()` | Complete-case outcome-model fit |
| `find_beta_sparcc()` | SPARCC outcome-model fit |
| `PRESCCO_prediction_interval()` | PRESCCO prediction interval |
| `split_conformal_prediction_interval()` | Split conformal prediction interval |
| `full_conformal_prediction_interval()` | Full conformal prediction interval |
| `jackknife_plus_prediction_interval()` | Jackknife+ prediction interval |
| `prediction_coverage_rate()` | Empirical coverage of a fitted half-length |
| `phi_xz()` | Outcome-model basis |
| `r1()` | Residual \(r_1\), or \(r_1^*\) under alternative working-model values |
| `r2()` | Residual \(r_2\) |

## Package structure

```text
prescco/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── LICENSE.md
├── README.md
├── .Rbuildignore
├── R/
│   ├── prescco-package.R
│   ├── predict_helpers.R
│   ├── nuisance.R
│   ├── eif_internals.R
│   ├── beta_fit.R
│   ├── predict_PRESCCO.R
│   └── predict_conformal.R
└── sim/
    ├── README.md
    ├── core.R
    ├── beta.R
    ├── PRESCCO.R
    ├── conformal.R
    ├── coverage.R
    ├── visualization.R
    ├── Rdata_files/
    └── results/
```

The `sim/` directory is excluded from the installed package.

## Simulation study

The simulation scripts can be run independently from the package root:

```sh
Rscript sim/beta.R
Rscript sim/PRESCCO.R
Rscript sim/conformal.R
Rscript sim/coverage.R
Rscript sim/visualization.R
```

Simulation `.Rdata` files are stored in `sim/Rdata_files/`. The supplied files
can be used directly to reproduce the summary table and figures without
rerunning the simulations:

```sh
Rscript sim/visualization.R
```

Outputs from the visualization script are written to `sim/results/`.

See `sim/README.md` for the simulation design and reproducibility details.

## Development checks

After modifying the package source:

```r
devtools::document()
devtools::check()
```
