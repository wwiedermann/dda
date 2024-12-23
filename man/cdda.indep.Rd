% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cdda.indep.R
\name{cdda.indep}
\alias{cdda.indep}
\title{Conditional Direction Dependence Analysis: Independence Assumptions}
\usage{
cdda.indep(
  formula = NULL,
  pred = NULL,
  mod = NULL,
  modval = "mean",
  data = list(),
  hetero = TRUE,
  diff = FALSE,
  nlfun = NULL,
  hsic.method = "gamma",
  B = 200,
  boot.type = "perc",
  conf.level = 0.95,
  parallelize = FALSE,
  cores = 1,
  ...
)
}
\arguments{
\item{formula}{Symbolic formula of the model to be tested or an \code{lm} object.}

\item{pred}{A character indicating the variable name of the predictor which serves as the outcome in the alternative model.}

\item{mod}{A character indicating the variable name of the moderator.}

\item{modval}{Characters or a numeric sequence specifying the moderator values used in post-hoc probing. Possible characters include \code{c("mean", "median", "JN")}. \code{modval = "mean"} tests the interaction effect at the moderator values M – 1SD, M, and M + 1SD; \code{modval = "median"} uses Q1, Md, and Q3. The Johnson-Neyman approach is applied when \code{modval = "JN"}. When a numeric sequence is specified, the pick-a-point approach will be utilized for the selected numeric values.}

\item{data}{A required data frame containing the variables in the model.}

\item{hetero}{A logical value indicating whether separate homoscedasticity tests (i.e., standard and robust Breusch-Pagan tests) should be computed.}

\item{diff}{A logical value indicating whether differences in HSIC, dCor, and MI values should be computed. Bootstrap confidence intervals are computed using B bootstrap samples.}

\item{nlfun}{Either a numeric value or a function of .Primitive type used for non-linear correlation tests. When nlfun is numeric the value is used in a power tranformation.}

\item{hsic.method}{A character indicating the inference method for Hilbert-Schmidt Independence Criterion. Must be one of the four values \code{c("gamma", "eigenvalue", "boot", "permutation")}.\code{hsic.method = "gamma"}is the default.}

\item{B}{Number of permutations for separate dCor tests and number of resamples when \code{hsic.method = c("boot", "permutation")} or \code{diff = TRUE}}

\item{boot.type}{A vector of character strings representing the type of bootstrap confidence intervals required. Must be one of the two values \code{c("perc", "bca")}. \code{boot.type = "bca"} is the default.}

\item{conf.level}{Confidence level for bootstrap confidence intervals}

\item{parallelize}{A logical value indicating whether bootstrapping is performed on multiple cores. Only used if \code{diff = TRUE}.}

\item{cores}{A numeric value indicating the number of cores. Only used if parallelize = TRUE}
}
\value{
A list of class \code{cddaindep} containing the results of CDDA
         independence tests for pre-specific moderator values.
}
\description{
\code{cdda.indep} computes CDDA test statistics to evaluate
             asymmetries of predictor-error independence of competing
             conditional models (\code{y ~ x | m} vs. \code{x ~ y | m}
             with \code{m} being a continuous or categorical moderator).
}
\examples{
set.seed(123)
n <- 1000

## --- generate moderator
z <- sort(rnorm(n))
z1 <- z[z <= 0]
z2 <- z[z > 0]

x1 <- rchisq(length(z1), df = 4) - 4
e1 <- rchisq(length(z1), df = 3) - 3
y1 <- 0.5 * x1 + e1

## --- y -> x when m > 0
y2 <- rchisq(length(z2), df = 4) - 4
e2 <- rchisq(length(z2), df = 3) - 3
x2 <- 0.25 * y2 + e2

y <- c(y1, y2); x <- c(x1, x2)

d <- data.frame(x, y, z)

cdda.indep(y ~ x * z, pred = "x", mod = "z",
           diff = TRUE, nlfun = 2, data = d)

m <- lm(y ~ x * z, data = d)

#' result <- cdda.indep(m, pred = "x", mod = "z", B = 500,
                     diff = TRUE, nlfun = 2, data = d)
 print(results)
 summary(result, hsic.diff = TRUE)

result <- cdda.indep(m, pred = "x", mod = "z", B = 500, modval = c(-0.5, 0.5),
                    diff = TRUE, nlfun = 2, data = d)
summary(result, hsic.diff = TRUE)

print(result)

}
\references{
Wiedermann, W., & von Eye, A. (2025). Direction Dependence Analysis: Foundations and Statistical Methods. Cambridge, UK: Cambridge University Press.
}
\seealso{
\code{\link{dda.indep}} for a non-conditional version of the function.
}
