#' Standard Deviation helper function
#' @description
#' Calculates the standard deviation of a vector.
#' @param x a vector
#' @param ...   Additional arguments to be passed to the function.
#' @returns Invisibly returns the standard deviation of a vector
#' @noRd
mysd <- function(x){sqrt(sum((x-mean(x))^2)/length(x))} #dda.vardist

#' Correlation helper function
#' @description Calculates the correlation between two vectors raised to the power of i and j.
#' @param x a vector
#' @param y a vector
#' @param i a numeric value
#' @param j a numeric value
#' @param ... other arguments not used by this method.
#' @returns Invisibly returns the correlation between two vectors raised to the power of i and j.
#' @noRd
cor.ij <- function(x,y, i=1, j=1){ #dda.vardist
  n <- length(x)
  mx <- mean(x)
  my <- mean(y)
  Cov <- sum((x - mx)^i * (y - my)^j)/n
  Cov/(mysd(x)^i * mysd(y)^j)
}

#' Bootstrap Difference helper function
#' @description Calculates the bootstrap difference between two vectors.
#' @param dat a data frame containing the predictor and outcome
#' @param g a vector of indices
#' @param ... other arguments not used by this method.
#' @returns Invisibly returns the difference between the skewness and kurtosis of two vectors.
#' @noRd
boot.diff <- function(dat, g){ #dda.indep & dda.vardist
  dat <- dat[g, ]
  x <- dat[, 1]  # "purified" predictor
  y <- dat[, 2]  # "purified" outcome

  x <- as.vector(scale(x))
  y <- as.vector(scale(y))

  skew.diff <- (moments::skewness(x)^2) - (moments::skewness(y)^2)
  kurt.diff <- (moments::kurtosis(x)-3)^2 - (moments::kurtosis(y)-3)^2
  cor12.diff <- (cor.ij(x, y, i = 2, j = 1)^2) - (cor.ij(x, y, i = 1, j = 2)^2)
  cor13.diff <- ((cor.ij(x, y, i = 3, j = 1)^2) - (cor.ij(x, y, i = 1, j = 3)^2)) * sign(moments::kurtosis(x)-3)

  Rtanh <- cor(x, y) * mean(x * tanh(y) - tanh(x) * y)

  Cxy <- mean(x^3 * y) - 3*cor(x,y)*var(x)
  Cyx <- mean(x * y^3) - 3*cor(x,y)*var(y)
  RCC <- (Cxy + Cyx) * (Cxy - Cyx)

  xx <- sign(moments::skewness(x)) * x
  yy <- sign(moments::skewness(y)) * y
  RHS <- cor(xx, yy) * mean( (xx^2 * yy) - (xx * yy^2) )

  result <- c(skew.diff, kurt.diff, cor12.diff, cor13.diff, RHS, RCC, Rtanh)
  names(result) <- c("skew.diff", "kurt.diff", "cor21.diff", "cor13.diff", "RHS", "RCC", "Rtanh")
  return(result)
}

#' Entropy helper function
#' @description Calculates the entropy of a vector.
#' @param x a vector
#' @param ... other arguments not used by this method.
#' @returns Invisibly returns the entropy of a vector
#' @noRd
max.entropy <- function(x){ #dda.indep
  sdx <- sd(x)
  x  <- as.vector(scale(x))
  k1 <- 79.047
  k2 <- 7.412889
  g  <- 0.37457
  GaussE <- log(2*pi)/2+1/2
  NegE <- k1 * (mean(log(cosh(x))) - g)^2 + k2 * mean(x * exp(-x^2/2))^2
  entropy <- GaussE - NegE + log(sdx)
  return(entropy)
}

#' Non-linear correlation (\code{nlcor}) tests
#' @description Calculates the non-linear correlation between two vectors.
#' @param x a numeric vector representing the tentative predictor.
#' @param y a numeric vector representing the tentative outcome.
#' @param fun a numeric value or a function of .Primitive type used for non-linear correlation tests. When \code{fun} is numeric the value is used in a power transformation.
#' @param ... other arguments not used by this method.
#' @returns Invisibly returns the correlation coefficient, t-value, degrees of freedom, and p-value for three different tests.
#' @noRd
nlcorTest <- function(x, y, fun, fname=NULL){ #dda.indep

  varnames <- c(deparse(substitute(x)), deparse(substitute(y)))

  if (length(x) != length(y)) stop("Variables must have same length")

  n <- length(x)
  x <- as.vector(scale(x))
  y <- as.vector(scale(y))

  if (is.numeric(fun)){
    func <- as.character(fun)
    r1 <- cor(x^fun, y)
    r2 <- cor(x, y^fun)
    r3 <- cor(x^fun, y^fun)

    if( any(is.na( c(r1, r2, r3) ) ) || any( is.nan( c(r1, r2, r3) ) ) ){

      x <- x + abs( min(x) ) + 0.1
      y <- y + abs( min(y) ) + 0.1

      r1 <- cor(x^fun, y)
      r2 <- cor(x, y^fun)
      r3 <- cor(x^fun, y^fun)
    }
  } # end if
  else {

    func <- paste(substitute(fun))

    test.run <- suppressWarnings( c(fun(x), fun(y) ) )

    if( any(is.na( test.run ) ) || any( is.nan( test.run ) ) ){
      x <- x + abs( min(x) ) + 0.1
      y <- y + abs( min(y) ) + 0.1
    } # end if

    r1 <- cor(fun(x), y)
    r2 <- cor(x, fun(y))
    r3 <- cor(fun(x), fun(y))

  } # end else = not is.numeric(fun)

  tval1 <- r1 * sqrt( ( n - 2)/(1 - r1^2))
  tval2 <- r2 * sqrt( ( n - 2)/(1 - r2^2))
  tval3 <- r3 * sqrt( ( n - 2)/(1 - r3^2))

  pval1 <- pt(abs(tval1), df = n - 2, lower.tail=FALSE) * 2
  pval2 <- pt(abs(tval2), df = n - 2, lower.tail=FALSE) * 2
  pval3 <- pt(abs(tval3), df = n - 2, lower.tail=FALSE) * 2

  output <- list(t1 = c(r1, tval1, n - 2, pval1),
                 t2 = c(r2, tval2, n - 2, pval2),
                 t3 = c(r3, tval3, n - 2, pval3),
                 func = fname,
                 varnames = varnames)

}
