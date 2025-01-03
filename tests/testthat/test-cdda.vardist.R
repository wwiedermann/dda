# Load necessary libraries
library(testthat)

# Set seed for reproducibility
set.seed(123)

# --- Data generation ---

n <- 500
z <- sort(rnorm(n))
z1 <- z[z <= 0]
z2 <- z[z > 0]

# --- x -> y when z <= 0 ---
x1 <- rchisq(length(z1), df = 4) - 4
e1 <- rchisq(length(z1), df = 3) - 3
y1 <- 0.5 * x1 + e1

# --- y -> x when z > 0 ---
y2 <- rchisq(length(z2), df = 4) - 4
e2 <- rchisq(length(z2), df = 3) - 3
x2 <- 0.5 * y2 + e2

# Combine data
y <- c(y1, y2)
x <- c(x1, x2)
dat <- data.frame(x, y, z)

# Fit model
m <- lm(y ~ x * z, data = dat)

# --- Define test object outside test_that blocks ---
test.cdda.vardist <- NULL

# --- Tests ---

test_that("cdda.vardist works", {
  # Test the function and assign the result to test.cdda.vardist
  test.cdda.vardist <- cdda.vardist(
    formula = m,
    pred = "x",
    mod = "z",
    modval = c(0, 1),
    data = dat,
    B = 200
  )

  # Check that it runs without warnings or messages
  expect_silent(test.cdda.vardist)
})

test_that("cdda.vardist summary works", {

  # Test the function and assign the result to test.cdda.vardist
  test.cdda.vardist <- cdda.vardist(
    formula = m,
    pred = "x",
    mod = "z",
    modval = c(0, 1),
    data = dat,
    B = 200
  )
  # Verify that summary can run without issues
  expect_output(summary(test.cdda.vardist, skew = TRUE, coskew = TRUE,
                        kurt = TRUE, cokurt = TRUE))
})
