#' @title Non-linear correlation (nlcor) tests
#' @description The function \code{nlcor.test} computes non-linear correlation tests between two variables. The function can be used to test for non-linear relationships between two variables. The function returns the correlation coefficient, t-value, degrees of freedom, and p-value for three different tests: (1) Correlation between x^fun and y, (2) Correlation between x and y^fun, and (3) Correlation between x^fun and y^fun. The function can be used to test for non-linear relationships between two variables. The function returns the correlation coefficient, t-value, degrees of freedom, and p-value for three different tests: (1) Correlation between x^fun and y, (2) Correlation between x and y^fun, and (3) Correlation between x^fun and y^fun.
#' @name nlcor.test
#' 
#' @param x     a numeric vector representing the tentative predictor.
#' @param y     a numeric vector representing the tentative outcome.
#' @param fun   a numeric value or a function of .Primitive type used for non-linear correlation tests. When \code{fun} is numeric the value is used in a power transformation.
#'
#' @returns     An object of class \code{dda.nlcor} containing the correlation coefficient, t-value, degrees of freedom, and p-value for three different tests.
#' @export

setClass("dda.nlcor", representation("list")) 

nlcor.test <- function(x, y, fun) {
 varnames <- c(deparse(substitute(x)), deparse(substitute(y)))


 if(length(x) != length(y)) stop("Variables must have same length")
 n <- length(x)
 
 x <- as.vector(scale(x))
 y <- as.vector(scale(y))
 
 if(is.numeric(fun)){
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
	 }
 else{
     func <- paste(substitute(fun))
	 test.run <- suppressWarnings( c(fun(x), fun(y) ) )
	   
	   if( any(is.na( test.run ) ) || any( is.nan( test.run ) ) ){
	
	       x <- x + abs( min(x) ) + 0.1
	       y <- y + abs( min(y) ) + 0.1
        }

     r1 <- cor(fun(x), y) 
	 r2 <- cor(x, fun(y)) 
	 r3 <- cor(fun(x), fun(y)) 

    }

  tval1 <- r1 * sqrt( ( n - 2)/(1 - r1^2))
  tval2 <- r2 * sqrt( ( n - 2)/(1 - r2^2))
  tval3 <- r3 * sqrt( ( n - 2)/(1 - r3^2))
  
  pval1 <- pt(abs(tval1), df = n - 2, lower.tail=FALSE) * 2 
  pval2 <- pt(abs(tval2), df = n - 2, lower.tail=FALSE) * 2 
  pval3 <- pt(abs(tval3), df = n - 2, lower.tail=FALSE) * 2 
  
  output <- list(t1 = c(r1, tval1, n - 2, pval1), 
                 t2 = c(r2, tval2, n - 2, pval2),
                 t3 = c(r3, tval3, n - 2, pval3), 
			     func = func, 
				 varnames = varnames)
  
  new ("dda.nlcor", output )
  
}

setMethod("show", "dda.nlcor", function(object){

     varnames <- object$varnames 
	 
	 sigtests <- rbind( object$t1, object$t2, object$t3)  
     sigtests <- round(sigtests, 4)
	 
     cat("\n")
	 if(is.na(suppressWarnings(as.numeric(object$func)))){
	            cat(paste("Non-linear correlation tests:", object$func, "transformation"))
	    
			    rownames(sigtests) <- c(paste("Cor[", object$func, "(", varnames[1], "), ", varnames[2],"]", sep=""),
	                                    paste("Cor[", varnames[1], ", ", object$func, "(", varnames[2], ")]", sep=""),
	                                    paste("Cor[", object$func, "(", varnames[1], "), ", object$func, "(", varnames[2],")]", sep="")
							           )
		} else{
	 	        cat(paste("Non-linear correlation tests: Power transformation using", object$func))
				
				rownames(sigtests) <- c(paste("Cor[", varnames[1], "^", object$func, ", ", varnames[2],"]", sep=""),
	                                    paste("Cor[", varnames[1], ", ", varnames[2], "^", object$func, "]", sep=""),
	                                    paste("Cor[", varnames[1], "^", object$func, ", ", varnames[2], "^", object$func, "]", sep="")
							           )
	    }
	 cat("\n")
	 
	 colnames(sigtests) <- c("estimate", "t-value", "df", "Pr(>|t|)")
     print.default(format( sigtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

})


# ----- nlcor.test() function for internal use in dda.indep(...)

nlcorTest <- function(x, y, fun, fname=NULL)
{

## ------------------------------------------------------------------------------------------------------
## x ... tentative predictor
## y ... tentative outcome 
## fun ... non-linear function. When fun is a numeric value the power transformation x^fun is performed
## Note: Function is used internally in the dda.indep(...)
## ------------------------------------------------------------------------------------------------------

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

