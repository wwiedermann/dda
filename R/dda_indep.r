#' @title Direction Dependence Analysis: Independence Distribution
#' @description \code{dda.indep} tests the skewness and kurtosis of the variables of two competing models. It also tests the difference in skewness and kurtosis to assess independence properties of the two models. The function also provides bootstrap confidence intervals for the difference in skewness and kurtosis.
#' @name dda.indep
#'
#' @param formula:      symbolic formula of the model to be tested or a \code{lm} object.
#' @param pred:         a character indicating the variable name of the predictor which serves as the outcome in the alternative model.
#' @param data:         an optional data frame containing the variables in the model (by default variables are taken from the environment which \code{dda.indep} is called from)
#' @param nlfun:        Either a numeric value or a function of .Primitive type used for non-linear correlation tests. When \code{nlfun} is numeric the value is used in a power tranformation.
#' @param hetero:       A logical value indicating whether separate homoscedasticity tests (i.e., standard and robust Breusch-Pagan tests) should be computed.
#' @param hsic.method:  a character indicating the inference method for Hilbert-Schmidt Independence Criterion (HSIC). Must be one of the four values \code{c("gamma", "eigenvalue", "boot", "permutation")}. \code{hsic.method = "gamma"} is the default.
#' @param diff:         A logical value indicating whether differences in HSIC, Distance Correlation (dCor), and MI values should be computed. Bootstrap confidence intervals are computed using \code{B} bootstrap samples.
#' @param B:            Number of permutations for separate dCor tests and number of resamples if \code{hsic.method = c("boot", "permutation")} or \code{diff = TRUE}
#' @param boot.type:    A vector of character strings representing the type of bootstrap confidence intervals required. Must be one of the two values \code{c("perc", "bca")}. \code{boot.type = "perc"} is the default.
#' @param conf.level:   confidence level for bootstrap confidence intervals
#' @param parallelize:  A logical value indicating whether boostrapping is performed on multiple cores. Only used if \code{diff = TRUE.}
#' @param cores:        a numeric value indicating the number of cores. Only used if parallelize = TRUE
#'
#' @examples dda.car.indep <- dda.indep(mpg ~ wt + qsec, pred = "wt",
#'                                      diff = TRUE, data = mtcars)
#'           dda.car.indep
#'           #OR
#'           car.test <- lm(mpg ~ wt + qsec, data = mtcars)
#'           dda.indep(car.test, pred = "wt", diff = TRUE, data = mtcars)
#'
#' @returns An object of class \code{ddaindep} containing the results of the independence tests.
#' @export
dda.indep <- function(formula, pred = NULL, data = list(), nlfun = NULL,
                      hetero = FALSE, hsic.method = "gamma", diff = FALSE,
                      B = 200, boot.type = "perc", conf.level = 0.95,
                      parallelize = FALSE, cores = 1, ...) {
  library(dHSIC)
  library(lmtest)
  library(energy)
  library(boot)

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

	   hsic.yx <- dhsic.test(rx, err.yx, method = hsic.method, kernel = "gaussian")
	   hsic.xy <- dhsic.test(ry, err.xy, method = hsic.method, kernel = "gaussian")

	   output <- list(hsic.yx = hsic.yx, hsic.xy = hsic.xy, hsic.method = hsic.method)

	}


	if(hsic.method %in% c("boot", "permutation")){

	   hsic.yx <- dhsic.test(rx, err.yx, method = hsic.method, kernel = "gaussian", B = B)
	   hsic.xy <- dhsic.test(ry, err.xy, method = hsic.method, kernel = "gaussian", B = B)

	   output <- list(hsic.yx = hsic.yx, hsic.xy = hsic.xy, hsic.method = c(hsic.method, as.character(B)) )

	}

	### --- separate dCor Tests

	dcor_yx <- dcor.test(as.vector(rx), as.vector(err.yx), R = B) #rx & ry have an SPSS attribute?
	dcor_xy <- dcor.test(as.vector(ry), as.vector(err.xy), R = B)

    output <- c(output,
	            distance_cor = list(dcor_yx = dcor_yx, dcor_xy = dcor_xy, dcor.method = as.character(B))
				)

	### --- Homoscedasticity tests

    if(hetero){

	  bp_yx <- bptest(m.yx, studentize = FALSE)
	  bp_xy <- bptest(m.xy, studentize = FALSE)

	  rbp_yx <- bptest(m.yx, studentize = TRUE)
	  rbp_xy <- bptest(m.xy, studentize = TRUE)

	  output <- c(output,
	              list(breusch_pagan = list( bp_yx, rbp_yx, bp_xy, rbp_xy ) )
	              )
	}


	### --- Non-linear correlation tests

	if(!is.null(nlfun)){

    fname = deparse(substitute(nlfun))

    nlout_yx <- unclass( nlcorTest( err.yx, rx, fun = nlfun, fname=fname ) )
	  nlout_xy <- unclass( nlcorTest( err.xy, ry, fun = nlfun, fname=fname ) )

	  output <- c(output,
	              list(nlcor.yx = nlout_yx, nlcor.xy = nlout_xy, nlfun = nlfun)
				 )
	}


    ### --- Difference statistics

    if(diff){

      dat.tmp <- data.frame(ry, err.xy, rx, err.yx)


	  if(!parallelize){

	      boot.res <- boot(dat.tmp, boot.diff, R = B)

		}

	  if(parallelize){

	      boot.res <- boot(dat.tmp, boot.diff, R = B, parallel = "snow", ncpus = cores)

	    }

    if( boot.type == "bca" && any( is.na( empinf( boot.res ) ) ) ) stop("Acceleration constant cannot be calculated. Increase the number of resamples or use boot.type = 'perc'")

	  suppressWarnings(boot.out <- lapply(as.list(1:3), function(i, boot.res) boot.ci(boot.res, conf=conf.level, type=boot.type, t0=boot.res$t0[i], t=boot.res$t[,i]), boot.res=boot.res))

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
  #new ("ddaindep", output )

  varnames <- output$var.names

	 cat("\n")
     cat("DIRECTION DEPENDENCE ANALYSIS: Independence Properties", "\n", "\n")

     # ------------------------------------------------------------------------------------------- Print Target Model:

	 cat(paste("Target Model:", varnames[2], "->", varnames[1], sep = " "), "\n", "\n")

	 cat("Omnibus Independence Tests:", "\n")

	 #if(output$hsic.method[1] == "gamma") cat("Hilbert-Schmidt Independence Criterion: Gamma-Approximation", "\n")
	 #if(output$hsic.method[1] == "boot") cat(paste("Hilbert-Schmidt Independence Criterion: Bootstrap-Approximation", " (", output$hsic.method[2], " resamples)", sep = ""), "\n")

	 cat(paste("HSIC = ", round(output$hsic.yx$statistic, 4), ", p-value = ", round(output$hsic.yx$p.value, 4), sep = ""))
	 cat("\n")
	 #cat(paste("Distance Correlation: Permutation", " (", output$distance_cor.dcor.method, " resamples)", sep = ""), "\n")
	 cat(paste("dCor = ", round(output$distance_cor.dcor_yx$statistic, 4), ", p-value = ", round(output$distance_cor.dcor_yx$p.value, 4), sep = ""))
     cat("\n", "\n")

	 if(!is.null(output$breusch_pagan)){

	 cat("Homoscedasticity Tests:", "\n")

	     sigtests1.yx <- rbind( c(output$breusch_pagan[[1]]$statistic, output$breusch_pagan[[1]]$parameter, output$breusch_pagan[[1]]$p.value),
	                            c(output$breusch_pagan[[2]]$statistic, output$breusch_pagan[[2]]$parameter, output$breusch_pagan[[2]]$p.value)
	                        )
         sigtests1.yx <- round(sigtests1.yx, 4)
         rownames(sigtests1.yx) <- c("BP-test", "Robust BP-test")
		 colnames(sigtests1.yx) <- c("X-squared", "df", "p-value")
         print.default(format( sigtests1.yx, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
	 cat("\n")

	 }


   if( !is.null(output$nlfun) ){

	     sigtests2.yx <- rbind(output$nlcor.yx$t1, output$nlcor.yx$t2, output$nlcor.yx$t3)
         sigtests2.yx <- round(sigtests2.yx, 4)

	      if(is.na(suppressWarnings(as.numeric(output$nlcor.yx$func)))){
	                cat(paste("Non-linear Correlation Tests:", output$nlcor.yx$func, "Transformation"))

			        rownames(sigtests2.yx) <- c(paste("Cor[", output$nlcor.yx$func, "(", "r_", varnames[1], "), ", varnames[2],"]", sep=""),
	                                            paste("Cor[", "r_", varnames[1], ", ", output$nlcor.yx$func, "(", varnames[2], ")]", sep=""),
	                                            paste("Cor[", output$nlcor.yx$func, "(", "r_", varnames[1], "), ", output$nlcor.yx$func, "(", varnames[2],")]", sep="")
							                   )
		    } else{
	 	            cat(paste("Non-linear Correlation Tests: Power Transformation using", output$nlcor.yx$func))

				    rownames(sigtests2.yx) <- c(paste("Cor[", "r_", varnames[1], "^", output$nlcor.yx$func, ", ", varnames[2],"]", sep=""),
	                                            paste("Cor[", "r_", varnames[1], ", ", varnames[2], "^", output$nlcor.yx$func, "]", sep=""),
	                                            paste("Cor[", "r_", varnames[1], "^", output$nlcor.yx$func, ", ", varnames[2], "^", output$nlcor.yx$func, "]", sep="")
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

	 #if(output$hsic.method[1] == "gamma") cat("Hilbert-Schmidt Independence Criterion: Gamma-Approximation", "\n")
	 #if(output$hsic.method[1] == "boot") cat(paste("Hilbert-Schmidt Independence Criterion: Bootstrap-Approximation", " (", output$hsic.method[2], " resamples)", sep = ""), "\n")

	 cat(paste("HSIC = ", round(output$hsic.xy$statistic, 4), ", p-value = ", round(output$hsic.xy$p.value, 4), sep = ""))
	 cat("\n")
	 #cat(paste("Distance Correlation: Permutation", " (", output$distance_cor.dcor.method, " resamples)", sep = ""), "\n")
	 cat(paste("dCor = ", round(output$distance_cor.dcor_xy$statistic, 4), ", p-value = ", round(output$distance_cor.dcor_xy$p.value, 4), sep = ""))
     cat("\n", "\n")

     if(!is.null(output$breusch_pagan)){

	 cat("Homoscedasticity Tests:", "\n")

	     sigtests1.xy <- rbind( c(output$breusch_pagan[[3]]$statistic, output$breusch_pagan[[3]]$parameter, output$breusch_pagan[[3]]$p.value),
	                            c(output$breusch_pagan[[4]]$statistic, output$breusch_pagan[[4]]$parameter, output$breusch_pagan[[4]]$p.value)
	                        )
         sigtests1.xy <- round(sigtests1.xy, 4)
         rownames(sigtests1.xy) <- c("BP-test", "Robust BP-test")
		 colnames(sigtests1.xy) <- c("X-squared", "df", "p-value")
         print.default(format( sigtests1.xy, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
	cat("\n")

     }

    if(!is.null(output$nlfun)){

	     sigtests2.xy <- rbind(output$nlcor.xy$t1, output$nlcor.xy$t2, output$nlcor.xy$t3)
         sigtests2.xy <- round(sigtests2.xy, 4)

	      if(is.na(suppressWarnings(as.numeric(output$nlcor.xy$func)))){
	                cat(paste("Non-linear Correlation Tests:", output$nlcor.xy$func, "Transformation"))

			        rownames(sigtests2.xy) <- c(paste("Cor[", output$nlcor.xy$func, "(", "r_", varnames[2], "), ", varnames[1],"]", sep=""),
	                                            paste("Cor[", "r_", varnames[2], ", ", output$nlcor.xy$func, "(", varnames[1], ")]", sep=""),
	                                            paste("Cor[", output$nlcor.xy$func, "(", "r_", varnames[2], "), ", output$nlcor.xy$func, "(", varnames[1],")]", sep="")
							                   )
		    } else{
	 	            cat(paste("Non-linear Correlation Tests: Power Transformation using", output$nlcor.xy$func))

				    rownames(sigtests2.xy) <- c(paste("Cor[", "r_", varnames[2], "^", output$nlcor.xy$func, ", ", varnames[1],"]", sep=""),
	                                            paste("Cor[", "r_", varnames[2], ", ", varnames[1], "^", output$nlcor.xy$func, "]", sep=""),
	                                            paste("Cor[", "r_", varnames[2], "^", output$nlcor.xy$func, ", ", varnames[1], "^", output$nlcor.xy$func, "]", sep="")
							                   )
	        }
     cat("\n")

	     colnames(sigtests2.xy) <- c("estimate", "t-value", "df", "Pr(>|t|)")
         print.default(format( sigtests2.xy, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

	 cat("\n")
	}


if(!is.null(output$out.diff)){

     ci.level <- as.numeric(output$boot.args[2]) * 100
     if(output$boot.args[1] == "bca") cat(ci.level, "% ", "BCa Bootstrap CIs for Difference Statistics", " (", output$boot.args[3], " samples):", "\n", sep = "")
     if(output$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile Bootstrap CIs for Difference Statistics", " (", output$boot.args[3], " samples):", "\n", sep = "")

	 print.default(format( output$out.diff, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
	 cat("\n")
	 cat("---")
	 cat("\n")
	 cat(paste("Note: Difference statistics > 0 suggest", varnames[2], "->", varnames[1], sep = " "))
	 cat("\n")
    }

   class(output) <- "ddaindep"

}

#q: I need major help, why isn't this function documenting correctly and appearing as a function under dda?
#a: You need to add the @export tag to the function. This will make the function available to the user when the package is loaded.

