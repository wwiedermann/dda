set.seed(123)

n <- 500

# --- generate moderator

z <- sort(rnorm(n))
z1 <- z[z <= 0]
z2 <- z[z > 0]

# --- x -> y when m <= 0

x1 <- rchisq(length(z1), df = 4) - 4
e1 <- rchisq(length(z1), df = 3) - 3
y1 <- 0.5 * x1 + e1

# --- y -> x when m > 0

y2 <- rchisq(length(z2), df = 4) - 4
e2 <- rchisq(length(z2), df = 3) - 3
x2 <- 0.5 * y2 + e2

y <- c(y1, y2)
x <- c(x1, x2)
dat <- data.frame(x,y,z)

m <- lm(y ~ x*z, data = dat)

test_that("cdda.indep works", {
  expect_silent(test.cdda.indep <- cdda.indep(m, pred = "x", mod = "z", dcor = TRUE,
                           diff = TRUE, modval = c(0, 1),
                           data = dat, nlfun = 2, B = 200))
})

test_that("cdda.indep summary works", {
  expect_silent(summary(test.cdda.indep, hsic.diff = TRUE, dcor.diff = TRUE, mi.diff = TRUE))
})


