test_that("find_alpha1_MLE recovers the mean of X", {
  skip_on_cran()
  d <- make_data(n = 1500, seed = 90, alpha1 = 0.3, alpha2 = 2, tau1 = 1)
  est <- find_alpha1_MLE(d$w, d$delta, w_min = d$w_min, w_max = d$w_max)

  expect_true(all(is.finite(est)))
  expect_equal(unname(est[1]), 0.3, tolerance = 0.35)
})

test_that("find_alpha2_MLE recovers the mean of C", {
  skip_on_cran()
  d <- make_data(n = 1500, seed = 91, alpha1 = 2, alpha2 = 0.3, tau2 = 1)
  est <- find_alpha2_MLE(d$w, d$delta, w_min = d$w_min, w_max = d$w_max)

  expect_true(all(is.finite(est)))
  expect_equal(unname(est[1]), 0.3, tolerance = 0.35)
})

test_that("the nuisance MLEs return one coefficient per covariate plus a scale", {
  skip_on_cran()
  for (p_c in 0:1) for (p_d in 0:1) {
    d <- make_data(n = 400, seed = 92 + p_c * 2 + p_d, p_c = p_c, p_d = p_d)
    e1 <- find_alpha1_MLE(d$w, d$delta, d$z_c, d$z_d, d$w_min, d$w_max)
    e2 <- find_alpha2_MLE(d$w, d$delta, d$z_c, d$z_d, d$w_min, d$w_max)
    # (intercept, p_c, p_d coefficients) + tau
    expect_length(e1, 1 + p_c + p_d + 1)
    expect_length(e2, 1 + p_c + p_d + 1)
    expect_true(all(is.finite(e1)))
    expect_true(all(is.finite(e2)))
  }
})

test_that("the estimated scale parameters are positive", {
  skip_on_cran()
  d <- make_data(n = 500, seed = 95)
  e1 <- find_alpha1_MLE(d$w, d$delta, w_min = d$w_min, w_max = d$w_max)
  e2 <- find_alpha2_MLE(d$w, d$delta, w_min = d$w_min, w_max = d$w_max)
  expect_gt(unname(e1[length(e1)]), 0)
  expect_gt(unname(e2[length(e2)]), 0)
})

test_that("find_alpha2_MLE is silent unless verbose", {
  skip_on_cran()
  d <- make_data(n = 200, seed = 96)
  expect_silent(find_alpha2_MLE(d$w, d$delta, w_min = d$w_min, w_max = d$w_max))
  expect_message(
    find_alpha2_MLE(d$w, d$delta, w_min = d$w_min, w_max = d$w_max,
                    verbose = TRUE)
  )
})

test_that("the nuisance MLEs reject mismatched input lengths", {
  expect_error(find_alpha1_MLE(1:5, rep(1L, 4)), "length n")
  expect_error(find_alpha2_MLE(1:5, rep(1L, 4)), "length n")
})
