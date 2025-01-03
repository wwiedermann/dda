#' @title Hilbert-Schmidt Independence Criterion (HSIC) Test
#' @description This function computes the Hilbert-Schmidt Independence Criterion (HSIC) test statistic for testing independence between two variables. The HSIC test is a non-parametric test that is based on the Hilbert-Schmidt norm of the cross-covariance operator. The null hypothesis is that the two variables are independent. The function also provides a p-value based on a bootstrap procedure.

#' @param model       An "lm" object representing the model to be tested.
#' @param x           A vector specifying the target predictor that is used to test independence.
#' @param hx          Numeric value specifying the bandwidth parameter for "x".
#' @param hy          Numeric value specifying the bandwidth parameter for "y".
#' @param B           Numeric value specifying the number of resamples used to compute the HSIC p-value.
#' @param parallelize A logical value indicating whether boostrapping is performed on multiple cores.
#' @param cores       A numeric value indicating the number of cores. Only used if parallelize = TRUE.
#' @param ...         Additional arguments to be passed to the function
#'
#' @examples    m <- lm(y ~ x + z)
#'              hsic.test(m, x, B = 500, parallelize = TRUE, cores = 4)
#' @noRd

## Delete all hsic.test document
hsic.test <- function(model, x = NULL, hx = 1, hy = 1, B = 1000, parallelize = FALSE, cores = 2)
{
    X <- model.matrix(model)
	b <- as.matrix(coef(model))
	n <- nobs(model)
	e <- resid(model)

	if(is.null(x)) stop("Predictor 'x' is missing.")
	if(!is.matrix(x)) x <- as.matrix(x)

	T_hat <- HSIC(x, e, 2*hx*hx, 2*hy*hy)   ## Compute the test-statistic

    # Implementing the bootstrap procedure (parallelized if specified)

	e_0 <- e - mean(e)      ## centered residuals

    if(parallelize){

      require(doParallel)
	  cl <- makeCluster(cores)
      registerDoParallel(cl)

         T_hat_B <- foreach(icount(B), .combine = c, .export = "HSIC") %dopar% {

               idx <- sample(n, n, replace=TRUE)    ## with replacement samples from the errors
        	   e_B <- e_0[idx]

        	   idx2 <- sample(n, n, replace=TRUE)   ## with replacement samples from the predictors
        	   x_B <- x[idx2,]
        	   X_B <- X[idx2,]

        	   yhat_B  <- X_B %*% b + e_B                              ## Create the bootstrap response values
        	   bhat_B  <- solve(t(X_B) %*% X_B) %*% t(X_B) %*% yhat_B  ## Get new slope estimates
        	   e_hat_B <- yhat_B - X_B %*% bhat_B                      ## Bootstrap residuals

        	   hx_B <- hx; hy_B <- hy
        	   HSIC(x_B, e_hat_B, 2*hx_B*hx_B, 2*hy_B*hy_B)
               }

	   stopCluster(cl)

	} else {

      T_hat_B <- matrix(0, B, 1)
      for(j in 1:B){
      	    idx <- sample(n, n, replace=TRUE)                ## with replacement samples from the errors
        	e_B <- e_0[idx]

        	idx2 <- sample(n, n, replace=TRUE)               ## with replacement samples from the predictors
        	x_B <- x[idx2,]
        	X_B <- X[idx2,]

        	yhat_B  <- X_B %*% b + e_B                               ## Create the bootstrap response values
        	bhat_B  <- solve(t(X_B) %*% X_B) %*% t(X_B) %*% yhat_B   ## Get new slope estimates
        	e_hat_B <- yhat_B - X_B %*% bhat_B                       ## Bootstrap residuals

        	hx_B <- hx; hy_B <- hy
        	T_hat_B[j] <- HSIC(x_B, e_hat_B, 2*hx_B*hx_B, 2*hy_B*hy_B)

 	    }
	}

	pval <- mean(T_hat_B >= T_hat)
	names(T_hat) <- "HSIC"
	#hist(T_hat_B, col = "grey"); abline(v = T_hat)
	output <- list(statistic = T_hat, p.value = pval,
	               method = "Hilbert-Schmidt Independence Test", data.name = deparse(formula(model)))
	class(output) <- "htest"
	return(output)
}

## --- Helper function to calculate HSIC value

HSIC <- function(x, y, hx, hy)
{
   n <- length(y)
   if(n != length(x)) stop("Variables must have same length.")
   H <- diag(n) - matrix(1, n, n)/n
   K <- exp(-as.matrix(dist(x, method = "euclidean", diag = TRUE, upper = TRUE))^2/hx)
   L <- exp(-as.matrix(dist(y, method = "euclidean", diag = TRUE, upper = TRUE))^2/hy)
   hsic <- sum(diag(K %*% H %*% L %*% H))/n
   return(hsic)
}

