# These exercise the interval-construction layer. Several fit a SPARCC beta
# internally and so are slow; they are skipped on CRAN but run locally.

test_that("get_zeta_seq_r brackets the empirical residual quantile", {
  d <- make_data(n = 300, seed = 100)

  seq_out <- get_zeta_seq_r(alpha = 0.1, r = r1,
                            y_data = d$y, w_data = d$w, delta_data = d$delta,
                            beta = d$beta, alpha1_star_r = d$alpha1[1],
                            tau1 = d$tau1,
                            w_min = d$w_min, w_max = d$w_max,
                            margin = 0.1, length = 5)

  res <- vapply(seq_len(d$n), function(i) {
    r1(d$y[i], d$w[i], d$delta[i], d$beta, d$alpha1[1], d$tau1,
       NULL, NULL, d$w_min, d$w_max)
  }, numeric(1))
  q <- unname(quantile(res, 0.9))

  expect_length(seq_out, 5)
  expect_true(all(diff(seq_out) > 0))          # strictly increasing
  expect_equal(seq_out[1], q * 0.9, tolerance = 1e-8)
  expect_equal(seq_out[5], q * 1.1, tolerance = 1e-8)
  expect_lt(seq_out[1], q)
  expect_gt(seq_out[5], q)
})

test_that("get_zeta_seq_r honours the length argument", {
  d <- make_data(n = 150, seed = 101)
  for (L in c(3, 5, 9)) {
    out <- get_zeta_seq_r(0.1, r1, d$y, d$w, d$delta,
                          beta = d$beta, alpha1_star_r = d$alpha1[1],
                          tau1 = d$tau1, w_min = d$w_min, w_max = d$w_max,
                          length = L)
    expect_length(out, L)
  }
})

test_that("prediction_coverage_rate returns a proportion and is monotone in zeta", {
  d <- make_data(n = 400, seed = 110)

  cvg <- function(z) {
    prediction_coverage_rate(z, r1, d$y, d$w, d$delta, NULL, NULL,
             d$beta, d$alpha1[1], d$tau1, d$w_min, d$w_max)
  }

  grid <- c(0.1, 0.5, 1, 2, 5, 20)
  vals <- vapply(grid, cvg, numeric(1))

  expect_true(all(vals >= 0 & vals <= 1))
  expect_true(all(diff(vals) >= 0))      # non-decreasing in zeta
  expect_equal(cvg(0), 0)                # no residual is <= 0 a.s.
  expect_equal(cvg(1e6), 1)              # everything covered
})

test_that("prediction_coverage_rate agrees with a direct calculation", {
  d <- make_data(n = 200, seed = 111)
  res <- vapply(seq_len(d$n), function(i) {
    r2(d$y[i], d$w[i], d$delta[i], d$beta, d$alpha1[1], d$tau1,
       NULL, NULL, d$w_min, d$w_max)
  }, numeric(1))
  z <- unname(quantile(res, 0.8))

  expect_equal(
    prediction_coverage_rate(z, r2, d$y, d$w, d$delta, NULL, NULL,
             d$beta, d$alpha1[1], d$tau1, d$w_min, d$w_max),
    mean(res <= z)
  )
})

test_that("split conformal returns positive, finite half-lengths", {
  skip_on_cran()
  d <- make_data(n = 300, seed = 120, alpha2 = 2)

  set.seed(7)
  out <- split_conformal_prediction_interval(
    y_data = d$y, w_data = d$w, delta_data = d$delta,
    alpha1 = d$alpha1, alpha2 = d$alpha2, alpha1_star_r = d$alpha1[1],
    alpha = 0.1, split_rate = 0.5,
    tau1 = d$tau1, tau2 = d$tau2, tau1_r = d$tau1,
    w_min = d$w_min, w_max = d$w_max
  )

  expect_equal(out$method, "split_conformal")
  expect_named(out$zeta, c("r1", "r2", "r1star"))
  expect_length(out$beta, 2)
  expect_true(all(is.finite(out$zeta)))
  expect_true(all(out$zeta > 0))
  # no test set supplied, so the coverage rate is NA by design
  expect_named(out$coverage_rate, c("r1", "r2", "r1star"))
  expect_true(all(is.na(out$coverage_rate)))
})

test_that("split conformal half-lengths widen as alpha shrinks", {
  skip_on_cran()
  d <- make_data(n = 400, seed = 121, alpha2 = 2)

  fit <- function(a) {
    set.seed(7)
    split_conformal_prediction_interval(
      y_data = d$y, w_data = d$w, delta_data = d$delta,
      alpha1 = d$alpha1, alpha2 = d$alpha2, alpha1_star_r = d$alpha1[1],
      alpha = a, split_rate = 0.5,
      tau1 = d$tau1, tau2 = d$tau2, tau1_r = d$tau1,
      w_min = d$w_min, w_max = d$w_max
    )
  }

  wide   <- fit(0.05)
  narrow <- fit(0.20)
  expect_gte(wide$zeta[["r1"]], narrow$zeta[["r1"]])
  expect_gte(wide$zeta[["r2"]], narrow$zeta[["r2"]])
})

test_that("split conformal attains roughly nominal coverage on fresh data", {
  skip_on_cran()
  d <- make_data(n = 600, seed = 130, alpha2 = 2)

  set.seed(7)
  out <- split_conformal_prediction_interval(
    y_data = d$y, w_data = d$w, delta_data = d$delta,
    alpha1 = d$alpha1, alpha2 = d$alpha2, alpha1_star_r = d$alpha1[1],
    alpha = 0.1, split_rate = 0.5,
    tau1 = d$tau1, tau2 = d$tau2, tau1_r = d$tau1,
    w_min = d$w_min, w_max = d$w_max
  )

  test <- make_data(n = 3000, seed = 131, alpha2 = 2)
  cvg <- prediction_coverage_rate(out$zeta[["r2"]], r2,
                  test$y, test$w, test$delta, NULL, NULL,
                  out$beta, d$alpha1[1], d$tau1, d$w_min, d$w_max)

  # Split conformal is finite-sample valid; allow generous Monte Carlo slack.
  expect_gt(cvg, 0.80)
  expect_lt(cvg, 0.99)
})


# The remaining exports (build_b123_arrays, find_zeta_param_int,
# PRESCCO_prediction_interval, full_conformal_prediction_interval, jackknife_plus_prediction_interval) are the expensive
# end-to-end paths. These smoke tests use the smallest grids that still
# exercise the code, and are skipped unless PRESCCO_SLOW_TESTS is set:
#
#   Sys.setenv(PRESCCO_SLOW_TESTS = "true"); devtools::test()

skip_unless_slow <- function() {
  testthat::skip_if_not(
    nzchar(Sys.getenv("PRESCCO_SLOW_TESTS")),
    "set PRESCCO_SLOW_TESTS to run the end-to-end tests"
  )
}

test_that("build_b123_arrays returns arrays of the documented shape", {
  skip_unless_slow()
  d <- make_data(n = 120, seed = 140, alpha2 = 2)
  zeta_seq <- get_zeta_seq_r(0.1, r1, d$y, d$w, d$delta,
                             beta = d$beta, alpha1_star_r = d$alpha1[1],
                             tau1 = d$tau1, w_min = d$w_min, w_max = d$w_max,
                             length = 3)

  out <- build_b123_arrays(
    zeta_seq = zeta_seq, alpha = 0.1, r = r1,
    beta = d$beta, x_a = seq(d$w_min, d$w_max, length.out = 5),
    z_c_data = NULL, z_d_data = NULL,
    alpha1_star = d$alpha1, alpha2_star = d$alpha2,
    alpha1_star_r = d$alpha1[1],
    sigma = d$sigma, tau1 = d$tau1, tau2 = d$tau2,
    w_min = d$w_min, w_max = d$w_max,
    b1_tt = 5, b1_tt2 = 10, b2_tt = 5, b2_tt2 = 10,
    b3_tt = 5, b3_tt2 = 10, b3_tt3 = 5
  )

  expect_type(out, "list")
  expect_true(all(vapply(out, function(a) all(is.finite(a)), logical(1))))
})

test_that("PRESCCO_prediction_interval runs end to end and returns a positive half-length", {
  skip_unless_slow()
  d <- make_data(n = 150, seed = 141, alpha2 = 2)

  out <- PRESCCO_prediction_interval(
    y_data = d$y, w_data = d$w, delta_data = d$delta,
    residual = r1, alpha = 0.1,
    w_min = d$w_min, w_max = d$w_max,
    seq_length = 3, tt = 5, m = 5,
    # small grids keep this test fast; build_b123_arrays' defaults are much
    # larger, and the study's values larger still
    b_args = list(b1_tt = 5, b1_tt2 = 10, b2_tt = 5, b2_tt2 = 10,
                  b3_tt = 5, b3_tt2 = 10, b3_tt3 = 5)
  )

  expect_type(out, "list")
  expect_true(all(c("method", "alpha", "residual", "zeta", "coverage_rate",
                    "beta", "sigma", "zeta_seq", "zeta_list") %in% names(out)))
  expect_length(out$beta, 2)
  expect_gt(out$sigma, 0)
  expect_equal(out$method, "PRESCCO")
  expect_length(out$zeta_seq, 3)
  expect_named(out$zeta, "r1")
  expect_true(is.finite(out$zeta))
  expect_gt(unname(out$zeta), 0)
  expect_true(is.na(out$coverage_rate[["r1"]]))
})

test_that("full conformal and jackknife+ return a coverage rate in [0, 1]", {
  skip_unless_slow()
  d  <- make_data(n = 60, seed = 142, alpha2 = 2)
  te <- make_data(n = 30, seed = 143, alpha2 = 2)

  args <- list(
    alpha = 0.1,
    y_data = d$y, w_data = d$w, delta_data = d$delta,
    test_y_data = te$y, test_w_data = te$w, test_delta_data = te$delta,
    alpha1 = d$alpha1, tau1 = d$tau1,
    alpha1_star_r = d$alpha1[1], tau1_r = d$tau1,
    w_min = d$w_min, w_max = d$w_max
  )

  for (fn in list(full_conformal_prediction_interval, jackknife_plus_prediction_interval)) {
    out <- do.call(fn, args)
    expect_named(out$zeta, c("r1", "r2", "r1star"))
    expect_named(out$coverage_rate, c("r1", "r2", "r1star"))
    expect_true(all(is.finite(out$zeta)))
    expect_true(all(out$zeta > 0))
    expect_true(all(out$coverage_rate >= 0 & out$coverage_rate <= 1))
    expect_equal(out$n_test, 30)
  }
})


test_that("all four methods share one return shape", {
  skip_unless_slow()
  d  <- make_data(n = 60, seed = 150, alpha2 = 2)
  te <- make_data(n = 30, seed = 151, alpha2 = 2)

  common <- c("method", "alpha", "residual", "zeta", "coverage_rate")

  scp <- split_conformal_prediction_interval(
    y_data = d$y, w_data = d$w, delta_data = d$delta,
    alpha1 = d$alpha1, alpha2 = d$alpha2, alpha1_star_r = d$alpha1[1],
    alpha = 0.1, tau1 = d$tau1, tau2 = d$tau2, tau1_r = d$tau1,
    test_y_data = te$y, test_w_data = te$w, test_delta_data = te$delta,
    w_min = d$w_min, w_max = d$w_max
  )

  fcp_args <- list(
    alpha = 0.1,
    y_data = d$y, w_data = d$w, delta_data = d$delta,
    test_y_data = te$y, test_w_data = te$w, test_delta_data = te$delta,
    alpha1 = d$alpha1, tau1 = d$tau1,
    alpha1_star_r = d$alpha1[1], tau1_r = d$tau1,
    w_min = d$w_min, w_max = d$w_max
  )

  for (out in list(scp,
                   do.call(full_conformal_prediction_interval, fcp_args),
                   do.call(jackknife_plus_prediction_interval, fcp_args))) {
    expect_true(all(common %in% names(out)), info = out$method)
    expect_named(out$zeta, c("r1", "r2", "r1star"))
    expect_named(out$coverage_rate, c("r1", "r2", "r1star"))
    expect_true(all(out$coverage_rate >= 0 & out$coverage_rate <= 1),
                info = out$method)
  }

  # split conformal now reports a real coverage rate, not NA, when given data
  expect_false(any(is.na(scp$coverage_rate)))
})
