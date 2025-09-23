#' Robust Breusch-Pagan Test for Heteroskedasticity
#'
#' This function performs a Breusch-Pagan test for heteroskedasticity,
#' but calculates the initial model residuals using a robust regression
#' (mblm::mblm) to be resistant to outliers.
#'
#' @param formula A regression model formula of the form y ~ x1 + x2.
#' @param varformula Optional formula for the variance model. Defaults to the main model's regressors.
#' @param studentize Logical. If TRUE (default), the studentized (Koenker-Bassett) version of the test is used.
#' @param data A data frame containing the variables.
#'
#' @return An object of class "htest" containing the test results.
#' @export
#'
#' @examples
#' # Create heteroskedastic data with an outlier
#' set.seed(123)
#' x <- 1:100
#' y <- 1 + 2*x + rnorm(100, 0, x)
#' y[50] <- 250 # Add a significant outlier
#'
#' # Run the robust Breusch-Pagan test
#' bptestrobust(y ~ x)

bptestrobust <- function(formula, varformula = NULL, studentize = TRUE, data = list())
{
  # Check for mblm package dependency
  if (!requireNamespace("mblm", quietly = TRUE)) {
    stop("Package 'mblm' is required. Please install it with install.packages('mblm').", call. = FALSE)
  }

  # --- Data Setup ---
  dname <- paste(deparse(substitute(formula)))
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  X <- model.matrix(formula, data = data)
  Z <- if (is.null(varformula)) X else model.matrix(varformula, data = data)

  if (ncol(Z) < 2L) {
    stop("The auxiliary variance regression requires at least an intercept and a regressor.")
  }

  n <- NROW(X)

  # --- Robust Residuals Calculation (mblm only) ---
  # Prepare data by combining the response and predictors into one data frame
  robust_df <- as.data.frame(cbind(y, X[, colnames(X) != "(Intercept)"]))

  # Fit the robust model and extract residuals

  robust_fit <- mblm::mblm(formula = as.formula(paste("y ~", paste(colnames(X)[colnames(X) != "(Intercept)"],
                                                                   collapse = "+"))),
                           data = robust_df, repeated = TRUE)

  resi <- residuals(robust_fit)

  # --- Test Statistic Calculation (using robust residuals) ---
  # The test uses an unweighted calculation
  sigma2 <- sum(resi^2) / n

  if (studentize) {
    w <- resi^2 - sigma2
    aux <- lm.fit(Z, w)
    bp <- n * sum(aux$fitted.values^2) / sum(w^2)
    method <- "Robust residual studentized Breusch-Pagan test"
  } else {
    f <- resi^2 / sigma2 - 1
    aux <- lm.fit(Z, f)
    bp <- 0.5 * sum(aux$fitted.values^2)
    method <- "Robust residual Breusch-Pagan test"
  }

  # --- Format and Return Results ---
  names(bp) <- "BP"
  df <- c(df = aux$rank - 1)

  RVAL <- list(statistic = bp,
               parameter = df,
               method = method,
               p.value = pchisq(bp, df, lower.tail = FALSE),
               data.name = dname)

  class(RVAL) <- "htest"
  return(RVAL)
}
