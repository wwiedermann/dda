#' @title Direction Dependence Analysis: Independence Properties
#' @description \code{dda.indep} computes DDA test statistics to
#'              evaluate asymmetries of predictor-error independence of
#'              causally competing models (\code{y ~ x} vs. \code{x ~ y}).
#'
#' @param formula      Symbolic formula of the model to be tested or a \code{lm} object.
#' @param pred         A character indicating the variable name of the predictor which serves as the outcome in the alternative model.
#' @param data         An optional data frame containing the variables in the model (by default variables are taken from the environment which \code{dda.indep}is called from).
#' @param nlfun        Either a numeric value or a function of .Primitive type used for non-linear correlation tests. When \code{nlfun} is numeric the value is used in a power transformation.
#' @param hetero       A logical value indicating whether separate homoscedasticity tests (i.e., standard and robust Breusch-Pagan tests) should be computed.
#' @param hsic.method  A character indicating the inference method for the Hilbert-Schmidt Independence Criterion (HSIC). Must be one of the four specifications \code{c("gamma", "eigenvalue", "boot", "permutation")}. \code{hsic.method = "gamma"}is the default.
#' @param diff         A logical value indicating whether differences in HSIC, Distance Correlation (dCor), and MI values should be computed. Bootstrap confidence intervals are computed using B bootstrap samples.
#' @param B            Number of permutations for separate dCor tests and number of resamples if \code{hsic.method = c("boot", "permutation")} or \code{diff = TRUE}.
#' @param boot.type    A vector of character strings representing the type of bootstrap confidence intervals. Must be one of the two specifications \code{c("perc", "bca")}.\code{boot.type = "perc"} is the default.
#' @param conf.level   Confidence level for bootstrap confidence intervals.
#' @param parallelize  A logical value indicating whether bootstrapping is performed on multiple cores. Only used if \code{diff = TRUE}.
#' @param cores        A numeric value indicating the number of cores. Only used if \code{parallelize = TRUE}.
#' @param ...          Additional arguments to be passed to the function.
#' @returns  An object of class \code{dda.indep} containing the results of DDA independence tests.
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
#' @rdname dda.indep
dda.indep <- function(
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
             cores = 1)
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

	m.yx <- lm(ry ~ rx)
	m.xy <- lm(rx ~ ry)

	err.yx <- resid(m.yx) #check here
	err.xy <- resid(m.xy)


	### --- Separate HSIC Tests

	if(hsic.method %in% c("gamma", "eigenvalue")){

	   hsic.yx <- dHSIC::dhsic.test(rx, err.yx, method = hsic.method, kernel = "gaussian")
	   hsic.xy <- dHSIC::dhsic.test(ry, err.xy, method = hsic.method, kernel = "gaussian")

	   output <- list(hsic.yx = hsic.yx, hsic.xy = hsic.xy, hsic.method = hsic.method)

	}


	if(hsic.method %in% c("boot", "permutation")){

	   hsic.yx <- dHSIC::dhsic.test(rx, err.yx, method = hsic.method, kernel = "gaussian", B = B)
	   hsic.xy <- dHSIC::dhsic.test(ry, err.xy, method = hsic.method, kernel = "gaussian", B = B)

	   output <- list(hsic.yx = hsic.yx, hsic.xy = hsic.xy, hsic.method = c(hsic.method, as.character(B)) )

	}

	### --- separate dCor Tests

	dcor_yx <- energy::dcor.test(as.vector(rx), as.vector(err.yx), R = B) #rx & ry have an SPSS attribute?
	dcor_xy <- energy::dcor.test(as.vector(ry), as.vector(err.xy), R = B)

    output <- c(output,
	            distance_cor = list(dcor_yx = dcor_yx, dcor_xy = dcor_xy, dcor.method = as.character(B))
				)

	### --- Homoscedasticity tests

    if(hetero){

	  bp_yx <- lmtest::bptest(m.yx, studentize = FALSE)
	  bp_xy <- lmtest::bptest(m.xy, studentize = FALSE)

	  rbp_yx <- lmtest::bptest(m.yx, studentize = TRUE)
	  rbp_xy <- lmtest::bptest(m.xy, studentize = TRUE)

	  output <- c(output,
	              list(breusch_pagan = list( bp_yx, rbp_yx, bp_xy, rbp_xy ) )
	              )
	}


	### --- Non-linear correlation tests

	if(!is.null(nlfun)){

    fname = deparse(substitute(nlfun))

    nlout_yx <- unclass( nlcor.test( err.yx, rx, fun = nlfun, fname=fname ) )
	  nlout_xy <- unclass( nlcor.test( err.xy, ry, fun = nlfun, fname=fname ) )

	  output <- c(output,
	              list(nlcor.yx = nlout_yx, nlcor.xy = nlout_xy, nlfun = nlfun)
				 )
	}


    ### --- Difference statistics

    if(diff){

      dat.tmp <- data.frame(ry, err.xy, rx, err.yx)


	  if(!parallelize){

	      boot.res <- boot::boot(dat.tmp, boot.diff, R = B)

		}

	  if(parallelize){

	      boot.res <- boot::boot(dat.tmp, boot.diff, R = B, parallel = "snow", ncpus = cores)

	    }

      if( boot.type == "bca" && any( is.na( boot::empinf( boot.res ) ) ) ) stop("Acceleration constant cannot be calculated. Increase the number of resamples or use boot.type = 'perc'")

	  suppressWarnings(boot.out <- lapply(as.list(1:3), function(i, boot.res) boot::boot.ci(boot.res, conf=conf.level, type=boot.type, t0=boot.res$t0[i], t=boot.res$t[,i]), boot.res=boot.res))

	  names(boot.out) <- c("diff.hsic", "diff.dcor", "diff.mi")

      diff.hsic.ci  <- unclass(boot.out$diff.hsic)[[4]][4:5] ; names(diff.hsic.ci) <- c("lower", "upper")
      diff.dcor.ci  <- unclass(boot.out$diff.dcor)[[4]][4:5] ; names(diff.dcor.ci) <- c("lower", "upper")
      diff.mi.ci    <- unclass(boot.out$diff.mi)[[4]][4:5] ; names(diff.mi.ci) <- c("lower", "upper")

      out.diff <- cbind(boot.res$t0, rbind(diff.hsic.ci, diff.dcor.ci, diff.mi.ci))
      colnames(out.diff) <- c("estimate", "lower", "upper")
      rownames(out.diff) <- c("HSIC", "dCor", "MI")

			 output <- c(output,
			             out.diff = list(out.diff),
						 list(boot.args = c(boot.type, conf.level, B)),
						 list(boot.warning = FALSE)
						)

     }

  response.name <- all.vars(formula(formula))[1]  # get name of response variable
  output <- c(output, list(var.names = c(response.name, pred)))

  class(output) <- "dda.indep"
  return(output)
}

#' @title Print Method for \code{dda.indep} Objects
#'
#' @description Prints the test statistics and p-values associated with \code{dda.indep} objects.
#'
#' @param x An object of class \code{dda.indep}.
#' @param ... Additional arguments to be passed to the method.
#'
#' @examples
#' print(result)
#'
#' @export
#' @rdname dda.indep
#' @method print dda.indep
print.dda.indep <- function(x, ...){

   varnames <- x$var.names

	 cat("\n")
   cat("DIRECTION DEPENDENCE ANALYSIS: Independence Properties", "\n", "\n")

     # ------------------------------------------------------------------------------------------- Print Target Model:

	 cat(paste("Target Model:", varnames[2], "->", varnames[1], sep = " "), "\n", "\n")

	 cat("Omnibus Independence Tests:", "\n")

	 #if(x$hsic.method[1] == "gamma") cat("Hilbert-Schmidt Independence Criterion: Gamma-Approximation", "\n")
	 #if(x$hsic.method[1] == "boot") cat(paste("Hilbert-Schmidt Independence Criterion: Bootstrap-Approximation", " (", x$hsic.method[2], " resamples)", sep = ""), "\n")

	 cat(paste("HSIC = ", round(x$hsic.yx$statistic, 4), ", p-value = ", round(x$hsic.yx$p.value, 4), sep = ""))
	 cat("\n")
	 #cat(paste("Distance Correlation: Permutation", " (", x$distance_cor.dcor.method, " resamples)", sep = ""), "\n")
	 cat(paste("dCor = ", round(x$distance_cor.dcor_yx$statistic, 4), ", p-value = ", round(x$distance_cor.dcor_yx$p.value, 4), sep = ""))
     cat("\n", "\n")

	 if(!is.null(x$breusch_pagan)){

	 cat("Homoscedasticity Tests:", "\n")

	     sigtests1.yx <- rbind( c(x$breusch_pagan[[1]]$statistic, x$breusch_pagan[[1]]$parameter, x$breusch_pagan[[1]]$p.value),
	                            c(x$breusch_pagan[[2]]$statistic, x$breusch_pagan[[2]]$parameter, x$breusch_pagan[[2]]$p.value)
	                        )
         sigtests1.yx <- round(sigtests1.yx, 4)
         rownames(sigtests1.yx) <- c("BP-test", "Robust BP-test")
		 colnames(sigtests1.yx) <- c("X-squared", "df", "p-value")
         print.default(format( sigtests1.yx, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
	 cat("\n")

	 }


   if( !is.null(x$nlfun) ){

	     sigtests2.yx <- rbind(x$nlcor.yx$t1, x$nlcor.yx$t2, x$nlcor.yx$t3)
         sigtests2.yx <- round(sigtests2.yx, 4)

	      if(is.na(suppressWarnings(as.numeric(x$nlcor.yx$func)))){
	                cat(paste("Non-linear Correlation Tests:", x$nlcor.yx$func, "Transformation"))

			        rownames(sigtests2.yx) <- c(paste("Cor[", x$nlcor.yx$func, "(", "r_", varnames[1], "), ", varnames[2],"]", sep=""),
	                                            paste("Cor[", "r_", varnames[1], ", ", x$nlcor.yx$func, "(", varnames[2], ")]", sep=""),
	                                            paste("Cor[", x$nlcor.yx$func, "(", "r_", varnames[1], "), ", x$nlcor.yx$func, "(", varnames[2],")]", sep="")
							                   )
		    } else{
	 	            cat(paste("Non-linear Correlation Tests: Power Transformation using", x$nlcor.yx$func))

				    rownames(sigtests2.yx) <- c(paste("Cor[", "r_", varnames[1], "^", x$nlcor.yx$func, ", ", varnames[2],"]", sep=""),
	                                            paste("Cor[", "r_", varnames[1], ", ", varnames[2], "^", x$nlcor.yx$func, "]", sep=""),
	                                            paste("Cor[", "r_", varnames[1], "^", x$nlcor.yx$func, ", ", varnames[2], "^", x$nlcor.yx$func, "]", sep="")
							                   )
	        }
     cat("\n")

	     colnames(sigtests2.yx) <- c("estimate", "t-value", "df", "Pr(>|t|)")
         print.default(format( sigtests2.yx, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

	 cat("\n")
    }

	 # ------------------------------------------------------------------------------------------- Print Alternative Model:

	 cat(paste("Alternative Model:", varnames[1], "->", varnames[2], sep = " "), "\n", "\n")

	 cat("Omnibus Independence Tests:", "\n")

	 #if(x$hsic.method[1] == "gamma") cat("Hilbert-Schmidt Independence Criterion: Gamma-Approximation", "\n")
	 #if(x$hsic.method[1] == "boot") cat(paste("Hilbert-Schmidt Independence Criterion: Bootstrap-Approximation", " (", x$hsic.method[2], " resamples)", sep = ""), "\n")

	 cat(paste("HSIC = ", round(x$hsic.xy$statistic, 4), ", p-value = ", round(x$hsic.xy$p.value, 4), sep = ""))
	 cat("\n")
	 #cat(paste("Distance Correlation: Permutation", " (", x$distance_cor.dcor.method, " resamples)", sep = ""), "\n")
	 cat(paste("dCor = ", round(x$distance_cor.dcor_xy$statistic, 4), ", p-value = ", round(x$distance_cor.dcor_xy$p.value, 4), sep = ""))
     cat("\n", "\n")

     if(!is.null(x$breusch_pagan)){

	 cat("Homoscedasticity Tests:", "\n")

	     sigtests1.xy <- rbind( c(x$breusch_pagan[[3]]$statistic, x$breusch_pagan[[3]]$parameter, x$breusch_pagan[[3]]$p.value),
	                            c(x$breusch_pagan[[4]]$statistic, x$breusch_pagan[[4]]$parameter, x$breusch_pagan[[4]]$p.value)
	                        )
         sigtests1.xy <- round(sigtests1.xy, 4)
         rownames(sigtests1.xy) <- c("BP-test", "Robust BP-test")
		 colnames(sigtests1.xy) <- c("X-squared", "df", "p-value")
         print.default(format( sigtests1.xy, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
	cat("\n")

     }

    if(!is.null(x$nlfun)){

	     sigtests2.xy <- rbind(x$nlcor.xy$t1, x$nlcor.xy$t2, x$nlcor.xy$t3)
         sigtests2.xy <- round(sigtests2.xy, 4)

	      if(is.na(suppressWarnings(as.numeric(x$nlcor.xy$func)))){
	                cat(paste("Non-linear Correlation Tests:", x$nlcor.xy$func, "Transformation"))

			        rownames(sigtests2.xy) <- c(paste("Cor[", x$nlcor.xy$func, "(", "r_", varnames[2], "), ", varnames[1],"]", sep=""),
	                                            paste("Cor[", "r_", varnames[2], ", ", x$nlcor.xy$func, "(", varnames[1], ")]", sep=""),
	                                            paste("Cor[", x$nlcor.xy$func, "(", "r_", varnames[2], "), ", x$nlcor.xy$func, "(", varnames[1],")]", sep="")
							                   )
		    } else{
	 	            cat(paste("Non-linear Correlation Tests: Power Transformation using", x$nlcor.xy$func))

				    rownames(sigtests2.xy) <- c(paste("Cor[", "r_", varnames[2], "^", x$nlcor.xy$func, ", ", varnames[1],"]", sep=""),
	                                            paste("Cor[", "r_", varnames[2], ", ", varnames[1], "^", x$nlcor.xy$func, "]", sep=""),
	                                            paste("Cor[", "r_", varnames[2], "^", x$nlcor.xy$func, ", ", varnames[1], "^", x$nlcor.xy$func, "]", sep="")
							                   )
	        }
     cat("\n")

	     colnames(sigtests2.xy) <- c("estimate", "t-value", "df", "Pr(>|t|)")
         print.default(format( sigtests2.xy, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

	 cat("\n")
	}


if(!is.null(x$out.diff)){

     ci.level <- as.numeric(x$boot.args[2]) * 100
     if(x$boot.args[1] == "bca") cat(ci.level, "% ", "BCa Bootstrap CIs for Difference Statistics", " (", x$boot.args[3], " samples):", "\n", sep = "")
     if(x$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile Bootstrap CIs for Difference Statistics", " (", x$boot.args[3], " samples):", "\n", sep = "")

	 print.default(format( x$out.diff, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
	 cat("\n")
	 cat("---")
	 cat("\n")
	 cat(paste("Note: Difference statistics > 0 suggest", varnames[2], "->", varnames[1], sep = " "))
	 cat("\n")
    }
}


