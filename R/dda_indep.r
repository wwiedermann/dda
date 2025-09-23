#' Robust Breusch-Pagan Test for Heteroskedasticity
#'
#' This function performs a Breusch-Pagan test using residuals from a
#' robust linear model (mblm) to be resistant to outliers. It accepts
#' either a formula or a fitted model object as input.
#'
#' @param formula A regression model formula or a fitted model object.
#' @param varformula Optional formula for the variance model. Defaults to the main model's regressors.
#' @param studentize Logical. If TRUE (default), the studentized (Koenker-Bassett) version of the test is used.
#' @param data A data frame containing the variables.
#' @param weights Optional vector of weights. NOTE: Weights are ignored in the robust residual calculation.
#'
#' @return An object of class "htest" containing the test results.

bptestrobust <- function(formula, varformula = NULL, studentize = TRUE, data = list(),
                         weights = NULL)
{
  # 1. Check for required package
  if (!requireNamespace("mblm", quietly = TRUE)) {
    stop("Package 'mblm' is required. Please install it with install.packages('mblm').", call. = FALSE)
  }

  # 2. Original, comprehensive data extraction logic (as requested)
  dname <- paste(deparse(substitute(formula)))
  if (!inherits(formula, "formula")) {
    X <- if (is.matrix(formula$x))
      formula$x
    else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(formula$y))
      formula$y
    else model.response(model.frame(formula))
    Z <- if (is.null(varformula))
      X
    else model.matrix(varformula, data = data)
    wts <- weights(formula)
  } else {
    mf <- if (is.null(weights)) {
      model.frame(formula, data = data)
    } else {
      model.frame(formula, weights = weights, data = data)
    }
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    Z <- if (is.null(varformula))
      X
    else model.matrix(varformula, data = data)
    wts <- model.weights(mf)
  }
  if (is.null(wts))
    wts <- rep.int(1, NROW(X))
  if (!(all(c(row.names(X) %in% row.names(Z), row.names(Z) %in%
              row.names(X))))) {
    allnames <- row.names(X)[row.names(X) %in% row.names(Z)]
    X <- X[allnames, ]
    Z <- Z[allnames, ]
    y <- y[allnames]
    wts <- wts[row.names(X) %in% row.names(Z)]
  }
  if (ncol(Z) < 2L)
    stop("the auxiliary variance regression requires at least an intercept and a regressor")

  # Warn if weights were provided, as mblm is unweighted
  if (!all(wts == 1)) {
    warning("Weights were provided but are ignored in the robust residual calculation.", call. = FALSE)
  }

  n <- NROW(X)

  # 3. Robust residual calculation using mblm
  model_data <- as.data.frame(X)
  model_data$y <- y
  predictor_names <- colnames(X)[colnames(X) != "(Intercept)"]

  if (length(predictor_names) == 0) {
    formula_str <- "y ~ 1"
  } else {
    formula_str <- paste("y ~", paste(predictor_names, collapse = " + "))
  }
  robust_formula <- as.formula(formula_str)

  robust_fit <- mblm::mblm(robust_formula, data = model_data, repeated = TRUE)
  resi <- residuals(robust_fit)

  # 4. Test statistic calculation (unweighted, as residuals are from an unweighted fit)
  sigma2 <- sum(resi^2) / n

  if (studentize) {
    w <- resi^2 - sigma2
    # Use lm.fit for unweighted auxiliary regression
    aux <- lm.fit(Z, w)
    bp <- n * sum(aux$fitted.values^2) / sum(w^2)
    method <- "Robust residual studentized Breusch-Pagan test"
  } else {
    f <- resi^2 / sigma2 - 1
    aux <- lm.fit(Z, f)
    bp <- 0.5 * sum(aux$fitted.values^2)
    method <- "Robust residual Breusch-Pagan test"
  }

  # 5. Format and return results
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
