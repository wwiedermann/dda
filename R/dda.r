#' @title Direction Dependence Analysis: Independence Properties
#' @description \code{dda.indep} computes DDA test statistics to
#'              evaluate asymmetries of predictor-error independence of
#'              causally competing models (\code{y ~ x} vs. \code{x ~ y}).
#'
#' @param formula      Symbolic formula of the model to be tested or a \code{lm} object.
#' @param pred         A character indicating the variable name of the predictor which serves as the outcome in the alternative model.
#' @param data         An optional data frame containing the variables in the model (by default variables are taken from the environment which \code{dda.indep} is called from).
#' @param nlfun        Either a numeric value or a function of .Primitive type used for non-linear correlation tests. When \code{nlfun} is numeric the value is used in a power transformation.
#' @param hetero       A logical value indicating whether separate homoscedasticity tests (i.e., standard and robust Breusch-Pagan tests) should be computed.
#' @param hsic.method  A character indicating the inference method for the Hilbert-Schmidt Independence Criterion (HSIC). Must be one of the four specifications \code{c("gamma", "eigenvalue", "boot", "permutation")}. \code{hsic.method = "gamma"}is the default.
#' @param diff         A logical value indicating whether differences in HSIC, Distance Correlation (dCor), and MI values should be computed. Bootstrap confidence intervals are computed using B bootstrap samples.
#' @param B            Number of permutations for separate dCor tests and number of resamples if \code{hsic.method = c("boot", "permutation")} or \code{diff = TRUE}.
#' @param boot.type    A vector of character strings representing the type of bootstrap confidence intervals. Must be one of the two specifications \code{c("perc", "bca")}.\code{boot.type = "perc"} is the default.
#' @param conf.level   Confidence level for bootstrap confidence intervals.
#' @param parallelize  A logical value indicating whether bootstrapping is performed on multiple cores. Only used if \code{diff = TRUE}.
#' @param cores        A numeric value indicating the number of cores. Only used if \code{parallelize = TRUE}.
#' @param robust       A logical value indicating whether Siegel non-parametric estimators should be used for residual extraction. If \code{robust = TRUE} Seigel estimation is used, otherwise ordinary least squares estimation is used.
#'
#' @returns
#' An object of class \code{dda.indep} containing the results of DDA independence tests.
#' @references Wiedermann, W., & von Eye, A. (2025). \emph{Direction Dependence Analysis: Foundations and Statistical Methods}. Cambridge, UK: Cambridge University Press.
#'
#' @examples
#'
#' set.seed(123)
#' n <- 500
#' x <- rchisq(n, df = 4) - 4
#' e <- rchisq(n, df = 3) - 3
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' result <- dda.indep(y ~ x, pred = "x", data = d, parallelize = TRUE, cores = 2,
#'           nlfun = 2, B = 50, hetero = TRUE, diff = TRUE)
#'
#'
#' @seealso \code{\link{cdda.indep}} for a conditional version.
#' @export
#' @rdname dda
#'
dda <- function(
             formula,
             pred = NULL,
             data = list(),
             nlfun = NULL,
             hetero = FALSE,
             hsic.method = "gamma",
             diff = FALSE,
             B = 200,
             boot.type = "perc",
             conf.level = 0.95,
             parallelize = FALSE,
             cores = 1,
             robust = FALSE)
  {

   ### --- helper functions for independence difference statistics

    max.entropy <- function(x){
           sdx <- sd(x)
           x  <- as.vector(scale(x))
           k1 <- 79.047
           k2 <- 7.412889
           g  <- 0.37457
           GaussE <- log(2 * pi) / 2 + 1 / 2
           NegE <- k1 * (mean(log(cosh(x))) - g) ^ 2 + k2 * mean(x * exp(-x ^ 2 / 2)) ^ 2
           entropy <- GaussE - NegE + log(sdx)
           return(entropy)
    }

    boot.diff <- function(dat, g){

      dat <- dat[g, ]

  	  ry     <- dat[,1] # purified outcome
      err.xy <- dat[,2] # errors of alternative model
      rx     <- dat[,3] # purified predictor
      err.yx <- dat[,4] # errors of target model

      diff.hsic <- dHSIC::dhsic.test(err.xy, ry, method = "gamma")$statistic - dHSIC::dhsic.test(err.yx, rx, method = "gamma")$statistic
      diff.dcor <- energy::dcor.test(err.xy, ry)$statistic - energy::dcor.test(err.yx, rx)$statistic
      diff.mi <- (max.entropy(ry) + max.entropy(err.xy)) - (max.entropy(rx) + max.entropy(err.yx))
      c(diff.hsic, diff.dcor, diff.mi)
    }


   ### --- non-linear correlation test function

   nlcor.test <- function(x, y, fun, fname=NULL){

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


   ### --- start checking validity of input

	if(is.null(pred)) stop( "Tentative predictor is missing." )
	if(B <= 0) stop( "Number of resamples 'B' must be positive." )
	if(conf.level < 0 || conf.level > 1) stop("'conf.level' must be between 0 and 1")
	if( !boot.type %in% c("bca", "perc") ) stop( "Unknown argument in boot.type." )

	if( !is.null(hsic.method)){
	    if (!hsic.method %in% c(NULL, "gamma", "boot", "permutation", "eigenvalue")) stop( "Unknown argument in hsic.method.")
	}

	### --- prepare outcome, predictor, and model matrix for covariates

  if (!inherits(formula, "formula") ) {
           X <- if (is.matrix(formula$x) ) formula$x
           else model.matrix(terms(formula), model.frame(formula) )
           y <- if (is.vector(formula$y) ) formula$y
           else model.response(model.frame(formula))

           delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
           if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

		   x <- X[,  delete.pred]  # tentative predictor
		   X <- X[, -delete.pred]  # model matrix with covariates
		   if ( !is.matrix(X) ) X <- as.matrix(X)
	}
	else {
       mf <- model.frame(formula, data = data)
 		   y  <- model.response(mf)   # tentative outcome
		   X  <- model.matrix(formula, data = data)

		   delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
       if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

		   x <- X[,  delete.pred] # tentative predictor
		   X <- X[, -delete.pred] # model matrix with covariates
		   if (!is.matrix(X)) X <- as.matrix(X)
	}

  ry <- lm.fit(X, y)$residuals
	rx <- lm.fit(X, x)$residuals

	resid_df <- data.frame(ry, rx) # create data frame with residuals

	if (robust == TRUE){
	  m.yx <- mblm::mblm(ry ~ rx, repeated = TRUE)
	  m.xy <- mblm::mblm(rx ~ ry, repeated = TRUE)
	}

	else if (robust == FALSE) {
	  m.yx <- lm(ry ~ rx)
	  m.xy <- lm(rx ~ ry)
	}

	else stop("Invalid specification for robust argument. Please use TRUE or FALSE.")

	err.yx <- resid(m.yx)
	err.xy <- resid(m.xy)


	# ?What output do we want to store?
	# Options: m.yx, m.xy, resid_df, err.yx, err.xy, ry, rx

   ### --- prepare output object
  output <- c(output, list(var.names = c(response.name, pred)))

  class(output) <- "dda"
  return(output)
}
