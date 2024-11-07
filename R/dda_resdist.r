#' @title Direction Dependence Analysis: Residual Distributions
#' @description This function tests the skewness and kurtosis of the residuals of two competing models. It also tests the difference in skewness and kurtosis between the residuals of the two models. The function also provides bootstrap confidence intervals for the difference in skewness and kurtosis.
#' @name dda.resdist
#'
#' @param formula    symbolic formula of the target model to be tested or a \code{lm} object
#' @param pred       variable name of the predictor which serves as the outcome in the alternative model
#' @param data       an optional data frame containing the variables in the model (by default variables are taken from the environment which \code{dda.vardist} is called from)
#' @param B          number of bootstrap samples
#' @param boot.type  A vector of character strings representing the type of bootstrap confidence intervals required. Must be one of the two values \code{c("perc", "bca")}; \code{boot.type = "bca" }is the default.
#' @param conf.level confidence level for boostrap confidence intervals
#'
#' @examples         dda.car.resdist <- dda.resdist(mpg ~ wt + qsec, pred = "wt",
#'                                                 boot.type = "bca", data = mtcars)
#'                   #OR
#'                   car.test <- lm(mpg ~ wt + qsec, data = mtcars)
#'                   dda.resdist(car.test, pred = "wt", boot.type = "bca", data = mtcars)
#'
#' @returns          An object of class \code{ddaresdist} containing the results of skewness and kurtosis tests, the difference in skewness and kurtosis, and bootstrap confidence intervals for the difference in skewness and kurtosis.
#' @export
dda.resdist <- function(formula, pred = NULL, data = list(), B = 100,
                        boot.type = "bca", conf.level = 0.95) {
   ### --- helper function for bootstrap CIs
  library(boot)
  library(moments)

    boot.diff <- function(dat, g){
              dat <- dat[g, ]
              x <- dat[, 1]  # "purified" predictor
              y <- dat[, 2]  # "purified" outcome

              skew.diff <- abs(skewness(x)) - abs(skewness(y))
              kurt.diff <- abs(kurtosis(x)-3) - abs(kurtosis(y)-3)

              result <- c(skew.diff, kurt.diff)
              names(result) <- c("skew.diff", "kurt.diff")
              return(result)
    }

    skew.diff.test <- function(x, y){

           agostino.zvalue <- function(x){
                     n  <- length(x)
                     s3 <- (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
                     y <- s3 * sqrt((n + 1) * (n + 3)/(6 * (n - 2)))
                     b2 <- 3 * (n * n + 27 * n - 70) * (n + 1) * (n + 3)/((n - 2) * (n + 5) * (n + 7) * (n + 9))
                     w <- -1 + sqrt(2 * (b2 - 1))
                     d <- 1/sqrt(log(sqrt(w)))
                     a <- sqrt(2/(w - 1))
                     z <- d * log(y/a + sqrt((y/a)^2 + 1))
                     return(z)
            }

    zval <- (agostino.zvalue(x) - agostino.zvalue(y))/sqrt(2 - 2*cor(x,y)^3)
		pval <- (1 - pnorm(abs(zval))) * 2 # two-sided pvalue
		return(list(z.value = abs(zval), p.value = pval))
	}

	kurt.diff.test <- function(x, y){

           anscombe.zvalue <- function(x){
                     n   <- length(x)
                     b   <- n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
                     eb2 <- 3 * (n - 1)/(n + 1)
                     vb2 <- 24 * n * (n - 2) * (n - 3)/((n + 1)^2 * (n + 3) * (n + 5))
                     m3  <- (6 * (n^2 - 5 * n + 2)/((n + 7) * (n + 9))) * sqrt((6 * (n + 3) * (n + 5))/(n * (n - 2) * (n - 3)))
                     a   <- 6 + (8/m3) * (2/m3 + sqrt(1 + 4/m3^2))
                     xx  <- (b - eb2)/sqrt(vb2)
                     z   <- (1 - 2/(9 * a) - ((1 - 2/a)/(1 + xx * sqrt(2/(a - 4))))^(1/3))/sqrt(2/(9 * a))
                     return(z)
            }

    zval <- (anscombe.zvalue(x) - anscombe.zvalue(y))/sqrt(2 - 2*cor(x,y)^4)
		pval <- (1 - pnorm(abs(zval))) * 2 # two-sided pvalue
		return(list(z.value = abs(zval), p.value = pval))
	}

  if(is.null(pred)) stop( "Tentative predictor is missing." )
	if(B < 0) stop( "Number of resamples 'B' must be positive." )
	if(conf.level < 0 || conf.level > 1) stop("'conf.level' must be between 0 and 1")
	if( !boot.type %in% c("bca", "perc") ) stop( "Unknown argument in boot.type." )

	### --- prepare outcome, predictor, and model matrix for covariates

    if (!inherits(formula, "formula") ) {
           X <- if (is.matrix(formula$x) ) formula$x
           else model.matrix(terms(formula), model.frame(formula) )
           y <- if (is.vector(formula$y) ) formula$y
           else model.response(model.frame(formula))

           delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
           if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

		   x <- X[, delete.pred]   # tentative predictor
		   X <- X[ , -delete.pred] # model matrix with covariates
		   if ( !is.matrix(X) ) X <- as.matrix(X)
	}
	else {
       mf <- model.frame(formula, data = data)
		   y <- model.response(mf)   # tentative outcome
		   X <- model.matrix(formula, data = data)

		   delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
           if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

		   x <- X[, delete.pred]   # tentative predictor
		   X <- X[ , -delete.pred] # model matrix with covariates
		   if (!is.matrix(X)) X <- as.matrix(X)
	}

  ry <- lm.fit(X, y)$residuals
	rx <- lm.fit(X, x)$residuals

	### --- estimate competing models

	tar <- lm(ry ~ rx)
	alt <- lm(rx ~ ry)
	dat <- data.frame(alternative = resid(alt), target = resid(tar))

	### --- run separate normality tests

	agostino.out <- apply(dat, 2, agostino.test)
	agostino.out <- lapply(agostino.out, unclass)
	agostino.out$alternative[3:5] <- NULL
    agostino.out$target[3:5] <- NULL

	anscombe.out <- apply(dat, 2, anscombe.test)
	anscombe.out <- lapply(anscombe.out, unclass)
	anscombe.out$alternative[3:5] <- NULL
    anscombe.out$target[3:5] <- NULL

	output <- list(agostino.out, anscombe.out)
	names(output) <- c("agostino", "anscombe")

	output$anscombe$target$statistic[1] <- output$anscombe$target$statistic[1] - 3 # change kurtosis to ex-kurtosis
	output$anscombe$alternative$statistic[1] <- output$anscombe$alternative$statistic[1] - 3 # change kurtosis to ex-kurtosis

	### --- run asymptotic difference tests

	output <- c(output,
	            list(skewdiff = unlist(skew.diff.test(dat$alternative, dat$target))),
				list(kurtdiff = unlist(kurt.diff.test(dat$alternative, dat$target)))
	)

	output$skewdiff <- c(abs(skewness(dat$alternative)) - abs(skewness(dat$target)), output$skewdiff)  # add point estimtes of skew and kurt differences to vector
	output$kurtdiff <- c(abs(kurtosis(dat$alternative)-3) - abs(kurtosis(dat$target)-3), output$kurtdiff)

    ### --- run bootstrap confidence intervals

   if(B > 0){
             boot.res <- boot(dat, boot.diff, R = B)
			 if( boot.type == "bca" && any( is.na( empinf( boot.res ) ) ) ) stop("Acceleration constant cannot be calculated. Increase the number of resamples or use boot.type = 'perc'")
             suppressWarnings(boot.out <- lapply(as.list(1:2), function(i, boot.res) boot.ci(boot.res, conf=conf.level, type=boot.type, t0=boot.res$t0[i], t=boot.res$t[,i]), boot.res=boot.res))
             names(boot.out) <- c("skew.diff", "kurt.diff")
             ci.skewdiff <- unclass(boot.out$skew.diff)[[4]][4:5] ; names(ci.skewdiff) <- c("lower", "upper")
			 ci.kurtdiff <- unclass(boot.out$kurt.diff)[[4]][4:5] ; names(ci.kurtdiff) <- c("lower", "upper")

			 output$skewdiff <- c(output$skewdiff, ci.skewdiff)
			 output$kurtdiff <-	c(output$kurtdiff, ci.kurtdiff)
			 output <- c(output, list(boot.args = c(boot.type, conf.level, B)))
	}

	response.name <- all.vars(formula(formula))[1]  # get name of response variable
	output <- c(output, list(var.names = c(response.name, pred)))

  varnames <- output$var.names

	 cat("\n")
     cat("DIRECTION DEPENDENCE ANALYSIS: Residual Distributions", "\n", "\n")
     cat("Skewness and kurtosis tests:", "\n")

	     sigtests <- rbind( c(output[[1]]$target$statistic, output[[1]]$target$p.value, output[[1]]$alternative$statistic, output[[1]]$alternative$p.value),
	                        c(output[[2]]$target$statistic, output[[2]]$target$p.value, output[[2]]$alternative$statistic, output[[2]]$alternative$p.value)
	                      )
       sigtests <- round(sigtests, 4)
       rownames(sigtests) <- c("Skewness", "Kurtosis")
		   colnames(sigtests) <- c("target", "z-value", "Pr(>|z|)", "alternative", "z-value", "Pr(>|z|)")
      print.default(format( sigtests, digits = max(3L, getOption("digits") - 3L), scientific = NA, scipen = 999), print.gap = 2L, quote = FALSE)

	 if(is.null(output$boot.args)){
		   cat("\n")
	       cat("Skewness and kurtosis difference tests:", "\n")

	       citests <- rbind(output$skewdiff, output$kurtdiff)
		     citests <- round(citests, 4)
	       rownames(citests) <- c("Skewness", "Kurtosis")
	       colnames(citests) <- c("diff", "z-value", "Pr(>|z|)")
	       print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
	       cat("\n")
        }

	 if(!is.null(output$boot.args)){
	       ci.level <- as.numeric(output$boot.args[2]) * 100
		     cat("\n")
	       cat("Skewness and kurtosis difference tests:", "\n")

	       citests <- rbind(output$skewdiff, output$kurtdiff)
                citests <- round(citests, 4)
	       rownames(citests) <- c("Skewness", "Kurtosis")
	       colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
	       print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
	       cat("\n")
	       cat(paste("Number of resamples:", output$boot.args[3]))
		     cat("\n")
		     if(output$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs reported", "\n", "\n", sep = "")
		     if(output$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs reported", "\n", "\n", sep = "")
          }

	  cat("---")
	  cat("\n")
	  cat(paste("Note: Target is", varnames[2], "->", varnames[1], sep = " "))
	  cat("\n")
	  cat(paste("      Alternative is", varnames[1], "->", varnames[2], sep = " "))
	  cat("\n")

	  #class(output) <- "ddaresdist"

}

