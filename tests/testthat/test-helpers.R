test_that("phi_xz has the documented layout and length", {
  expect_equal(phi_xz(2), c(1, 2))
  expect_equal(phi_xz(2, z_c = 5), c(1, 2, 5, 10))
  expect_equal(phi_xz(2, z_d = 3), c(1, 2, 3, 6))
  expect_equal(phi_xz(2, z_c = 5, z_d = 3), c(1, 2, 5, 3, 10, 6))

  # length grows as 2 + 2 * (p_c + p_d)
  for (p_c in 0:2) for (p_d in 0:2) {
    zc <- if (p_c > 0) seq_len(p_c) else NULL
    zd <- if (p_d > 0) seq_len(p_d) else NULL
    expect_length(phi_xz(1.5, zc, zd), 2 + 2 * (p_c + p_d))
  }
})

test_that("m0 is the inner product of the basis and beta", {
  beta <- c(0.5, -1.2, 2, 3, -0.4, 0.9)
  expect_equal(m0(0.3, beta, z_c = 0.7, z_d = 1),
               sum(phi_xz(0.3, 0.7, 1) * beta))

  # no-covariate case reduces to a straight line
  expect_equal(m0(0.25, c(1, 4)), 1 + 4 * 0.25)
})

test_that("m1 equals m0(w) whenever X is observed", {
  # delta == 1 means X = W, so no integration should happen
  beta <- c(0.2, 2.5, 1.1, -0.6, 0.3, 0.8)
  for (w in c(-0.9, -0.2, 0, 0.4, 0.99)) {
    expect_equal(
      m1(w, delta = 1, beta, alpha1_star = c(0.1, 0.2, -0.3), tau1 = 1,
         z_c = 0.5, z_d = 1, w_min = W_MIN, w_max = W_MAX),
      m0(w, beta, z_c = 0.5, z_d = 1)
    )
  }
})

test_that("m1 under censoring lies beyond m0(w) in the direction of the slope", {
  # delta == 0 implies X > w, so the imputed X exceeds w and m1 moves away
  # from m0(w) in whichever direction the slope points.
  beta_up   <- c(0, 3)    # positive slope
  beta_down <- c(0, -3)   # negative slope
  w <- -0.5

  m1_up <- m1(w, delta = 0, beta_up, alpha1_star = 0, tau1 = 1,
              w_min = W_MIN, w_max = W_MAX)
  m1_dn <- m1(w, delta = 0, beta_down, alpha1_star = 0, tau1 = 1,
              w_min = W_MIN, w_max = W_MAX)

  expect_gt(m1_up, m0(w, beta_up))
  expect_lt(m1_dn, m0(w, beta_down))
})

test_that("m1 stays inside the support bounds", {
  # The imputed X is a truncated mean on (w, w_max], so m1 must sit between
  # m0(w) and m0(w_max) for a monotone outcome model.
  beta <- c(0, 3)
  w <- -0.5
  val <- m1(w, delta = 0, beta, alpha1_star = 0, tau1 = 1,
            w_min = W_MIN, w_max = W_MAX)
  expect_gte(val, m0(w, beta))
  expect_lte(val, m0(W_MAX, beta))
})

test_that("m1 falls back to m0 at or above the upper bound", {
  beta <- c(0, 3)
  expect_equal(
    m1(W_MAX, delta = 0, beta, alpha1_star = 0, tau1 = 1,
       w_min = W_MIN, w_max = W_MAX),
    m0(W_MAX, beta)
  )
})

test_that("r1 and r2 are non-negative and match their definitions", {
  beta <- c(0.2, 2.5)
  args <- list(w = -0.3, delta = 0L, beta = beta,
               alpha1_star = 0, tau1 = 1, w_min = W_MIN, w_max = W_MAX)
  y <- 1.7

  r1_val <- do.call(r1, c(list(y = y), args))
  r2_val <- do.call(r2, c(list(y = y), args))

  expect_gte(r1_val, 0)
  expect_gte(r2_val, 0)
  expect_equal(r2_val, abs(y - m0(args$w, beta)))
  expect_equal(r1_val, abs(y - do.call(m1, args)))
})

test_that("r1 and r2 agree exactly when X is observed", {
  # delta == 1 collapses m1 to m0(w), so the two residuals coincide
  beta <- c(0.2, 2.5)
  expect_equal(
    r1(1.0, 0.4, 1L, beta, alpha1_star = 0, tau1 = 1,
       w_min = W_MIN, w_max = W_MAX),
    r2(1.0, 0.4, 1L, beta, alpha1_star = 0, tau1 = 1,
       w_min = W_MIN, w_max = W_MAX)
  )
})

test_that("r1 is centered differently from r1* (misspecified center)", {
  # r1* is r1 evaluated at a wrong alpha1_star_r; on a censored observation the
  # two must differ, which is the whole point of the r1* scenarios.
  beta <- c(0, 3)
  common <- list(y = 1.0, w = -0.5, delta = 0L, beta = beta, tau1 = 1,
                 w_min = W_MIN, w_max = W_MAX)
  r1_correct <- do.call(r1, c(common, list(alpha1_star =  0)))
  r1_star    <- do.call(r1, c(common, list(alpha1_star = -2)))
  expect_false(isTRUE(all.equal(r1_correct, r1_star)))
})
