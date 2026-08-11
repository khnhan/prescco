test_that("find_beta_cc reproduces lm() on the complete cases", {
  d <- make_data(n = 300, seed = 40, p_c = 1, p_d = 1)
  fit <- find_beta_cc(d$y, d$w, d$delta, z_c_data = d$z_c, z_d_data = d$z_d)

  keep <- d$delta == 1
  ref <- lm(d$y[keep] ~ d$w[keep] + d$z_c[keep, 1] + d$z_d[keep, 1] +
              I(d$w[keep] * d$z_c[keep, 1]) + I(d$w[keep] * d$z_d[keep, 1]))

  expect_equal(fit$beta_cc, unname(coef(ref)), tolerance = 1e-8)
  expect_equal(fit$p_c, 1L)
  expect_equal(fit$p_d, 1L)
  expect_equal(fit$idx_cc, which(keep))
})

test_that("find_beta_cc returns the right dimensions in every covariate case", {
  for (p_c in 0:2) for (p_d in 0:1) {
    d <- make_data(n = 150, seed = 50 + p_c * 3 + p_d, p_c = p_c, p_d = p_d)
    fit <- find_beta_cc(d$y, d$w, d$delta, d$z_c, d$z_d)
    expect_length(fit$beta_cc, 2 + 2 * (p_c + p_d))
    expect_true(all(is.finite(fit$beta_cc)))
    expect_gt(fit$sigma_cc, 0)
  }
})

test_that("find_beta_cc is nearly unbiased when nothing is censored", {
  # With delta == 1 throughout, CC is just OLS on the full sample.
  d <- make_data(n = 4000, seed = 60, sigma = 1, censor = FALSE)
  fit <- find_beta_cc(d$y, d$w, d$delta)
  expect_equal(fit$beta_cc, d$beta, tolerance = 0.1)
  expect_equal(fit$sigma_cc, d$sigma, tolerance = 0.1)
})

test_that("find_beta_cc validates its inputs", {
  expect_error(find_beta_cc(1:5, 1:4, rep(1L, 5)), "same length")
  expect_error(find_beta_cc(1:5, 1:5, rep(0L, 5)), "No complete cases")
})

test_that("find_beta_sparcc recovers beta on lightly censored data", {
  skip_on_cran()
  d <- make_data(n = 400, seed = 70, alpha2 = 3, sigma = 1)  # light censoring

  fit <- find_beta_sparcc(
    y_data = d$y, w_data = d$w, delta_data = d$delta,
    z_c_data = NULL, z_d_data = NULL,
    alpha1_star = d$alpha1, alpha2_star = d$alpha2,
    tau1 = d$tau1, tau2 = d$tau2,
    w_min = d$w_min, w_max = d$w_max
  )

  expect_length(fit$beta_hat, 2)
  expect_true(all(is.finite(fit$beta_hat)))
  expect_gt(fit$sigma_hat, 0)
  expect_equal(fit$beta_hat, d$beta, tolerance = 0.5)
})

test_that("find_beta_sparcc is silent unless verbose", {
  skip_on_cran()
  d <- make_data(n = 150, seed = 80, alpha2 = 3)
  args <- list(y_data = d$y, w_data = d$w, delta_data = d$delta,
               z_c_data = NULL, z_d_data = NULL,
               alpha1_star = d$alpha1, alpha2_star = d$alpha2,
               tau1 = d$tau1, tau2 = d$tau2,
               w_min = d$w_min, w_max = d$w_max)

  # expect_no_message rather than expect_silent: nleqslv may legitimately warn
  # on a hard solve, which is not what this test is about.
  expect_no_message(do.call(find_beta_sparcc, args))
  expect_message(do.call(find_beta_sparcc, c(args, list(verbose = TRUE))),
                 "nleqslv")
})
