# The *_vec functions are the fast paths used inside the conformal routines.
# They must agree with their scalar counterparts elementwise. These tests cover
# every combination of covariate blocks, because the p_c > 0 AND p_d > 0 case
# exercises code that the no-covariate simulation study never touches.

test_that("m0_vec matches m0 elementwise, with and without covariates", {
  for (p_c in 0:1) for (p_d in 0:1) {
    d <- make_data(n = 40, seed = 10 + p_c * 2 + p_d, p_c = p_c, p_d = p_d)

    vec <- m0_vec(d$w, d$beta, z_c = d$z_c, z_d = d$z_d)
    scal <- vapply(seq_len(d$n), function(i) {
      m0(d$w[i], d$beta, row_or_null(d$z_c, i), row_or_null(d$z_d, i))
    }, numeric(1))

    expect_equal(vec, scal,
                 info = sprintf("p_c = %d, p_d = %d", p_c, p_d))
  }
})

test_that("m0_vec keeps the discrete block when a continuous block is present", {
  # Regression test: an unparenthesized `if (p_c > 0) A else 0 + if (p_d > 0) B
  # else 0` parses as A + (0 + B-or-0) only when p_c == 0, silently dropping the
  # discrete contribution whenever p_c > 0. Changing only z_d must change m0_vec.
  n <- 5
  w  <- seq(-0.8, 0.8, length.out = n)
  zc <- matrix(0.5, n, 1)
  beta <- c(0, 3, 1.5, -2.5, 0.4, 0.9)   # (1, x, z_c, z_d, x*z_c, x*z_d)

  out_zd0 <- m0_vec(w, beta, z_c = zc, z_d = matrix(0, n, 1))
  out_zd1 <- m0_vec(w, beta, z_c = zc, z_d = matrix(1, n, 1))

  expect_false(isTRUE(all.equal(out_zd0, out_zd1)))

  # and the difference is exactly the z_d main effect plus its interaction
  expect_equal(out_zd1 - out_zd0, beta[4] + beta[6] * w)
})

test_that("m1_vec matches m1 elementwise, with and without covariates", {
  for (p_c in 0:1) for (p_d in 0:1) {
    d <- make_data(n = 40, seed = 20 + p_c * 2 + p_d, p_c = p_c, p_d = p_d)

    vec <- m1_vec(d$w, d$delta, d$beta, d$alpha1, tau1 = d$tau1,
                  z_c = d$z_c, z_d = d$z_d, w_max = d$w_max)
    scal <- vapply(seq_len(d$n), function(i) {
      m1(d$w[i], d$delta[i], d$beta, d$alpha1, d$tau1,
         row_or_null(d$z_c, i), row_or_null(d$z_d, i),
         w_min = d$w_min, w_max = d$w_max)
    }, numeric(1))

    expect_equal(vec, scal, tolerance = 1e-8,
                 info = sprintf("p_c = %d, p_d = %d", p_c, p_d))
  }
})

test_that("m1_vec keeps the discrete block when a continuous block is present", {
  n <- 6
  w     <- seq(-0.8, 0.5, length.out = n)
  delta <- c(1L, 0L, 1L, 0L, 1L, 0L)
  zc    <- matrix(0.5, n, 1)
  beta  <- c(0, 3, 1.5, -2.5, 0.4, 0.9)
  a1    <- c(0, 0.3, -0.2)

  out_zd0 <- m1_vec(w, delta, beta, a1, tau1 = 1,
                    z_c = zc, z_d = matrix(0, n, 1), w_max = W_MAX)
  out_zd1 <- m1_vec(w, delta, beta, a1, tau1 = 1,
                    z_c = zc, z_d = matrix(1, n, 1), w_max = W_MAX)

  expect_false(isTRUE(all.equal(out_zd0, out_zd1)))
})

test_that("r1_vec and r2_vec match their scalar counterparts", {
  d <- make_data(n = 40, seed = 30, p_c = 1, p_d = 1)

  r1v <- r1_vec(d$y, d$w, d$delta, d$beta, d$alpha1, tau1 = d$tau1,
                z_c = d$z_c, z_d = d$z_d, w_max = d$w_max)
  r2v <- r2_vec(d$y, d$w, d$delta, d$beta, d$alpha1, tau1 = d$tau1,
                z_c = d$z_c, z_d = d$z_d, w_max = d$w_max)

  r1s <- vapply(seq_len(d$n), function(i) {
    r1(d$y[i], d$w[i], d$delta[i], d$beta, d$alpha1, d$tau1,
       row_or_null(d$z_c, i), row_or_null(d$z_d, i),
       w_min = d$w_min, w_max = d$w_max)
  }, numeric(1))
  r2s <- vapply(seq_len(d$n), function(i) {
    r2(d$y[i], d$w[i], d$delta[i], d$beta, d$alpha1, d$tau1,
       row_or_null(d$z_c, i), row_or_null(d$z_d, i),
       w_min = d$w_min, w_max = d$w_max)
  }, numeric(1))

  expect_equal(r1v, r1s, tolerance = 1e-8)
  expect_equal(r2v, r2s, tolerance = 1e-8)
  expect_true(all(r1v >= 0))
  expect_true(all(r2v >= 0))
})

test_that("vectorized functions reject mismatched lengths", {
  expect_error(r1_vec(1:3, 1:2, c(1L, 1L), c(0, 1), 0), "same length")
  expect_error(r2_vec(1:3, 1:2, c(1L, 1L), c(0, 1), 0), "same length")
  expect_error(m1_vec(1:3, c(1L, 1L), c(0, 1), 0, z_c = 1:2), "length n")
})

test_that("m1_vec checks the alpha1_star dimension", {
  expect_error(
    m1_vec(c(0.1, 0.2), c(1L, 0L), c(0, 3, 1, 0.5),
           alpha1_star = c(0, 0, 0), z_c = c(0.5, 0.5)),
    "1 \\+ p_c \\+ p_d"
  )
})
