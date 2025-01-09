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

# --- Tests ---

testthat::test_that("dda.vardist works", {
  # Test the function and assign the result to test.dda.vardist
  test.dda.vardist <- dda.vardist(
    formula = m,
    pred = "x",
    data = dat,
    B = 200
  )

  # Check that it runs without warnings or messages
  expect_silent(test.dda.vardist)
})

testthat::test_that("dda.vardist print works", {
  # Test the print functionality
  test.dda.vardist <- dda.vardist(
    formula = m,
    pred = "x",
    data = dat,
    B = 200
  )

  # Verify that print can run without issues
  testthat::expect_output(print(test.dda.vardist))
})
