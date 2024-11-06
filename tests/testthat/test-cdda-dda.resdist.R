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

test_that("cdda.vardist works", {
  expect_silent(test.dda.resdist <- dda.resdist(m, pred = "x", mod = "z",
                                                modval = c(0, 1), data = dat,
                                                B = 200))
})

test_that("cdda.vardist print works", {
  expect_silent(print(test.dda.vardist))
})
