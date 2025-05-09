% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dda_resdist.r
\name{dda.resdist}
\alias{dda.resdist}
\alias{print.dda.resdist}
\title{Direction Dependence Analysis: Residual Distributions}
\usage{
dda.resdist(
  formula,
  pred = NULL,
  data = list(),
  B = 200,
  boot.type = "perc",
  prob.trans = FALSE,
  conf.level = 0.95
)

\method{print}{dda.resdist}(x, ...)
}
\arguments{
\item{formula}{Symbolic formula of the target model to be tested or a \code{lm} object.}

\item{pred}{Variable name of the predictor which serves as the outcome in the alternative model.}

\item{data}{An optional data frame containing the variables in the
model (by default variables are taken from the environment
which \code{dda.resdist} is called from).}

\item{B}{Number of bootstrap samples.}

\item{boot.type}{A vector of character strings representing the type of
bootstrap confidence intervals required. Must be one of
the two values \code{c("perc", "bca")};
\code{boot.type = "perc"} is the default.}

\item{prob.trans}{A logical value indicating whether a probability integral
transformation should be performed prior computation of
skewness and kurtosis difference tests.}

\item{conf.level}{Confidence level for bootstrap confidence intervals.}

\item{x}{An object of class \code{dda.resdist} when using \code{print}.}

\item{...}{Additional arguments to be passed to the method.}
}
\value{
An object of class \code{ddaresdist} containing the results of DDA
          tests of asymmetry patterns of error distributions obtained from
          the causally competing models.
}
\description{
\code{dda.resdist} evaluates patterns of asymmetry of error
distributions of causally competing models (\code{y ~ x} vs. \code{x ~ y}).

\code{print} returns DDA test statistics associated with \code{dda.resdist} objects.
}
\examples{
set.seed(123)
n <- 500
x <- rchisq(n, df = 4) - 4
e <- rchisq(n, df = 3) - 3
y <- 0.5 * x + e
d <- data.frame(x, y)

result <- dda.resdist(y ~ x, pred = "x", data = d,
            B = 50, conf.level = 0.90, prob.trans = TRUE)

print(result)

}
\references{
Wiedermann, W., & von Eye, A. (2025). Direction Dependence Analysis: Foundations and Statistical Methods. Cambridge, UK: Cambridge University Press.
}
