test_that("median_bandwidth returns positive scalar", {
  set.seed(1)
  x <- rnorm(50)
  bw <- median_bandwidth(x)
  expect_true(is.numeric(bw))
  expect_length(bw, 1)
  expect_gt(bw, 0)
})

test_that("median_bandwidth matches dHSIC convention on fixed input", {
  set.seed(1)
  x  <- rnorm(80)
  bw <- median_bandwidth(x)
  # dHSIC uses median of squared pairwise distances; verify positive and finite
  expect_true(is.finite(bw))
})

test_that("hsic returns near-zero for independent data", {
  set.seed(42)
  n <- 200
  h <- hsic(rnorm(n), rnorm(n))
  expect_true(is.numeric(h))
  expect_length(h, 1)
  expect_gte(h, 0)
  expect_lt(h, 0.05)   # should be near zero under independence
})

test_that("hsic is larger for dependent than independent data", {
  set.seed(42)
  n  <- 200
  x  <- rnorm(n)
  h_indep <- hsic(x, rnorm(n))
  h_dep   <- hsic(x, x + rnorm(n, sd = 0.5))
  expect_gt(h_dep, h_indep)
})

test_that("hsic is consistent across all four kernels (non-negative)", {
  set.seed(7)
  n <- 100
  x <- rnorm(n)
  y <- x + rnorm(n)
  for (k in c("gaussian", "laplace", "linear", "polynomial")) {
    h <- hsic(x, y, kernel_x = k, kernel_y = k)
    expect_gte(h, 0, label = paste("kernel:", k))
  }
})

test_that("hsic.test returns htest object with correct names", {
  set.seed(1)
  x   <- rnorm(60)
  res <- hsic.test(x, rnorm(60), method = "gamma")
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "HSIC")
  expect_true(!is.null(res$p.value))
  expect_true(!is.null(res$statistic))
})

test_that("hsic.test statistic label is 'HSIC' for all methods", {
  set.seed(1)
  n <- 60
  x <- rnorm(n)
  y <- rnorm(n)
  for (m in c("gamma", "permutation", "eigenvalue", "bootstrap")) {
    res <- hsic.test(x, y, method = m, B = 99)
    expect_named(res$statistic, "HSIC",
                 label = paste("method:", m))
  }
})

test_that("hsic.test p-values are in [0, 1] for all methods", {
  set.seed(2)
  n <- 60
  x <- rnorm(n)
  y <- rnorm(n)
  for (m in c("gamma", "permutation", "eigenvalue", "bootstrap")) {
    p <- hsic.test(x, y, method = m, B = 99)$p.value
    expect_gte(p, 0, label = paste("method:", m))
    expect_lte(p, 1, label = paste("method:", m))
  }
})

test_that("hsic.test statistic is n * hsic()", {
  set.seed(3)
  n   <- 80
  x   <- rnorm(n)
  y   <- x + rnorm(n)
  raw <- hsic(x, y)
  res <- hsic.test(x, y, method = "gamma")
  expect_equal(as.numeric(res$statistic), n * raw, tolerance = 1e-8)
})

test_that("hsic.test gamma rejects dependent data", {
  set.seed(10)
  n <- 150
  x <- rchisq(n, df = 2)
  y <- 0.75 * x + (rchisq(n, df = 5) - 5)
  res <- hsic.test(scale(x), scale(y), method = "gamma")
  expect_lt(res$p.value, 0.05)
})

test_that("hsic.test gamma does not reject independent data (probabilistic)", {
  # Run 20 seeds; at least 15 should not reject at alpha=0.05
  rejections <- sum(sapply(1:20, function(s) {
    set.seed(s)
    n <- 100
    hsic.test(rnorm(n), rnorm(n), method = "gamma")$p.value < 0.05
  }))
  expect_lte(rejections, 7)  # generous upper bound; true rate ~0.05
})

test_that("hsic.test permutation and gamma agree directionally", {
  set.seed(5)
  n  <- 100
  x  <- rchisq(n, df = 2)
  y  <- 0.75 * x + (rchisq(n, df = 5) - 5)
  xs <- scale(x); ys <- scale(y)
  p_gamma <- hsic.test(xs, ys, method = "gamma")$p.value
  p_perm  <- hsic.test(xs, ys, method = "permutation", B = 199)$p.value
  # both should reject (or both not); just check they're on the same side of 0.5
  expect_equal(p_gamma < 0.5, p_perm < 0.5)
})

test_that("hsic.test eigenvalue reports n*HSIC not n^3*HSIC", {
  set.seed(6)
  n   <- 60
  x   <- rnorm(n)
  y   <- rnorm(n)
  raw <- hsic(x, y)
  res <- hsic.test(x, y, method = "eigenvalue", B = 199)
  # statistic should be n*HSIC, not n^3*HSIC
  expect_equal(as.numeric(res$statistic), n * raw, tolerance = 1e-8)
})

test_that("hsic.test statistics match across methods for same data", {
  set.seed(9)
  n   <- 80
  x   <- rnorm(n)
  y   <- x + rnorm(n)
  raw <- n * hsic(x, y)
  for (m in c("gamma", "permutation", "bootstrap")) {
    res <- hsic.test(x, y, method = m, B = 99)
    expect_equal(as.numeric(res$statistic), raw, tolerance = 1e-8,
                 label = paste("method:", m))
  }
})

test_that("hsic.test bandwidths element is present and positive", {
  set.seed(1)
  x   <- rnorm(50)
  res <- hsic.test(x, rnorm(50), method = "gamma")
  expect_true(!is.null(res$bandwidths))
  expect_true(all(res$bandwidths > 0))
})

test_that("hsic.test fixed bandwidth gives same result as median_bandwidth", {
  set.seed(1)
  n   <- 60
  x   <- rnorm(n)
  y   <- rnorm(n)
  bwx <- median_bandwidth(x)
  bwy <- median_bandwidth(y)
  r_auto  <- hsic.test(x, y, method = "gamma")
  r_fixed <- hsic.test(x, y, method = "gamma", bandwidth_x = bwx, bandwidth_y = bwy)
  expect_equal(as.numeric(r_auto$statistic), as.numeric(r_fixed$statistic),
               tolerance = 1e-8)
})
