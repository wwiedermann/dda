# ---- helpers --------------------------------------------------------

#' Build a Kernel (Gram) Matrix
#'
#' @description Constructs an \eqn{n \times n} kernel matrix for a dataset
#'   using one of four supported kernel functions: Gaussian (RBF), Laplace,
#'   linear, or polynomial.  This is an internal workhorse called by
#'   \code{\link{hsic}}, \code{\link{hsic_test}}, and
#'   \code{\link{hsic_resid_test}}.
#'
#' @param x A numeric vector (length \eqn{n}) or matrix (\eqn{n \times p}).
#' @param kernel Character string.  One of \code{"gaussian"} (default),
#'   \code{"laplace"}, \code{"linear"}, or \code{"polynomial"}.
#' @param bandwidth Scalar bandwidth parameter.
#'   \itemize{
#'     \item \strong{Gaussian}: \eqn{h = 2\sigma^2}, so the kernel is
#'       \eqn{k(x,y) = \exp(-\|x-y\|^2 / h)}.  This is the same
#'       parameterisation used in the \pkg{dHSIC} source (Peters et al., 2022)
#'       where the argument is called \code{hx} or \code{hy}.  Suzuki (2025,
#'       Example 34) uses \eqn{\sigma^2 = 1}, which corresponds to
#'       \code{bandwidth = 2} here.
#'     \item \strong{Laplace}: \eqn{\sigma}, so
#'       \eqn{k(x,y) = \exp(-\|x-y\| / \sigma)}.
#'     \item \strong{Linear / Polynomial}: ignored.
#'   }
#'   Default \code{NULL} (triggers median heuristic).
#' @param degree Integer degree \eqn{d} for the polynomial kernel
#'   \eqn{k(x,y) = (x^\top y + c)^d}.  Default \code{2}.
#' @param coef0 Constant \eqn{c} for the polynomial kernel.  Default \code{1}.
#'
#' @return An \eqn{n \times n} numeric matrix.
#'
#' @keywords internal
build_kernel_matrix <- function(x,
                                kernel    = c("gaussian", "laplace",
                                              "linear",  "polynomial"),
                                bandwidth = NULL,
                                degree    = 2,
                                coef0     = 1) {
  kernel <- match.arg(kernel)
  if (is.vector(x)) x <- matrix(x, ncol = 1)

  # Safeguard: If bandwidth isn't provided, apply median heuristic safely
  if (is.null(bandwidth) && kernel %in% c("gaussian", "laplace")) {
    bandwidth <- median_bandwidth(x)
  }

  if (bandwidth <= 0 && kernel %in% c("gaussian", "laplace"))
    stop("'bandwidth' must be strictly positive.")

  D <- as.matrix(dist(x, method = "euclidean", diag = TRUE, upper = TRUE))

  switch(kernel,
         gaussian   = exp(-D^2 / bandwidth),
         laplace    = exp(-D  / bandwidth),
         linear     = tcrossprod(x),
         polynomial = (tcrossprod(x) + coef0)^degree
  )
}


# Kernel Matrix
#'
#' @description Applies the centering transform
#'   \eqn{\widetilde{K} = H K H}
#'   where \eqn{H = I_n - n^{-1}\mathbf{1}\mathbf{1}^\top} is the
#'   \eqn{n \times n} centering matrix.
#'
#' @param K An \eqn{n \times n} kernel matrix.
#'
#' @return An \eqn{n \times n} centered kernel matrix.
#'
#' @keywords internal
center_kernel <- function(K) {
  n <- nrow(K)
  H <- diag(n) - matrix(1 / n, n, n)
  H %*% K %*% H
}

#' Compute Empirical HSIC
#'
#' @description Computes the empirical Hilbert-Schmidt Independence Criterion
#'   (HSIC) between two random variables \eqn{X} and \eqn{Y}.
#'
#' @param x A numeric vector (length \eqn{n}) or matrix (\eqn{n \times p})
#'   of observations of \eqn{X}.
#' @param y A numeric vector (length \eqn{n}) or matrix (\eqn{n \times q})
#'   of observations of \eqn{Y}.
#' @param kernel_x Character string.  Kernel for \eqn{X}.
#' @param kernel_y Character string.  Kernel for \eqn{Y}.  Defaults to \code{kernel_x}.
#' @param bandwidth_x Bandwidth for the \eqn{X} kernel. Default \code{NULL} (median heuristic).
#' @param bandwidth_y Bandwidth for the \eqn{Y} kernel. Default \code{NULL} (median heuristic).
#' @param degree Integer.  Polynomial degree.  Default \code{2}.
#' @param coef0 Numeric.  Polynomial constant.  Default \code{1}.
#'
#' @return A single non-negative numeric: \eqn{\widehat{\mathrm{HSIC}}(X, Y)}.
#'
#' @export
hsic <- function(x, y,
                 kernel_x    = "gaussian",
                 kernel_y    = kernel_x,
                 bandwidth_x = NULL,
                 bandwidth_y = NULL,
                 degree      = 2,
                 coef0       = 1) {
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  if (is.vector(y)) y <- matrix(y, ncol = 1)
  n <- nrow(x)
  if (nrow(y) != n)
    stop("'x' and 'y' must have the same number of observations.")

  # Dynamic bandwidth defaults matching dHSIC
  if (is.null(bandwidth_x)) bandwidth_x <- median_bandwidth(x)
  if (is.null(bandwidth_y)) bandwidth_y <- median_bandwidth(y)

  Kx <- build_kernel_matrix(x, kernel = kernel_x, bandwidth = bandwidth_x,
                            degree = degree, coef0 = coef0)
  Ky <- build_kernel_matrix(y, kernel = kernel_y, bandwidth = bandwidth_y,
                            degree = degree, coef0 = coef0)

  Kx_c <- center_kernel(Kx)
  Ky_c <- center_kernel(Ky)

  sum(diag(Kx_c %*% Ky_c)) / n^2
}


#' HSIC Independence Test
#'
#' @description Tests whether two random variables \eqn{X} and \eqn{Y} are
#'   independent using the Hilbert-Schmidt Independence Criterion (HSIC).
#'
#' @param x A numeric vector or matrix (\eqn{n \times p}) of observations of \eqn{X}.
#' @param y A numeric vector or matrix (\eqn{n \times q}) of observations of \eqn{Y}.
#' @param method Character string.  Null-distribution method.
#' @param kernel_x Character string.  Kernel for \eqn{X}.
#' @param kernel_y Character string.  Kernel for \eqn{Y}.
#' @param bandwidth_x Scalar bandwidth for the \eqn{X} kernel. Default \code{NULL}.
#' @param bandwidth_y Scalar bandwidth for the \eqn{Y} kernel. Default \code{NULL}.
#' @param degree Integer.  Polynomial degree. Default \code{2}.
#' @param coef0 Numeric.  Polynomial constant. Default \code{1}.
#' @param B Integer.  Number of replications. Default \code{1000}.
#' @param parallelize Logical.  If \code{TRUE}, parallelises the replication loop.
#' @param cores Integer.  Number of cores when \code{parallelize = TRUE}.
#'
#' @return An object of class \code{"htest"}
#'
#' @export
hsic_test <- function(x, y,
                      method      = c("permutation", "asymptotic", "bootstrap"),
                      kernel_x    = "gaussian",
                      kernel_y    = kernel_x,
                      bandwidth_x = NULL,
                      bandwidth_y = NULL,
                      degree      = 2,
                      coef0       = 1,
                      B           = 1000,
                      parallelize = FALSE,
                      cores       = 2) {

  method <- match.arg(method)
  dname  <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  if (is.vector(x)) x <- matrix(x, ncol = 1)
  if (is.vector(y)) y <- matrix(y, ncol = 1)
  n <- nrow(x)
  if (nrow(y) != n)
    stop("'x' and 'y' must have the same number of observations.")

  # Dynamic bandwidth defaults
  if (is.null(bandwidth_x)) bandwidth_x <- median_bandwidth(x)
  if (is.null(bandwidth_y)) bandwidth_y <- median_bandwidth(y)

  # Kernel matrices (computed once)
  Kx   <- build_kernel_matrix(x, kernel = kernel_x, bandwidth = bandwidth_x,
                              degree = degree, coef0 = coef0)
  Ky   <- build_kernel_matrix(y, kernel = kernel_y, bandwidth = bandwidth_y,
                              degree = degree, coef0 = coef0)
  Kx_c <- center_kernel(Kx)
  Ky_c <- center_kernel(Ky)

  T_obs_raw <- sum(diag(Kx_c %*% Ky_c)) / n^2
  T_obs_scaled <- T_obs_raw * n # Scale by n to mirror dHSIC::dhsic.test

  # ---------- permutation --------------------------------------------------
  if (method == "permutation") {

    perm_one <- function() {
      idx   <- sample(n)
      Kx_p  <- center_kernel(Kx[idx, idx])
      (sum(diag(Kx_p %*% Ky_c)) / n^2) * n
    }

    if (parallelize) {
      cl     <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      T_null <- foreach::foreach(
        iterators::icount(B), .combine = c,
        .export = c("center_kernel")
      ) %dopar% perm_one()
      parallel::stopCluster(cl)
    } else {
      T_null <- replicate(B, perm_one())
    }

    pval <- mean(T_null >= T_obs_scaled)
    stat <- T_obs_scaled
    names(stat) <- "HSIC (scaled by n)"
    meth <- sprintf(
      "HSIC permutation test (B = %d, %s kernel)", B, kernel_x
    )
  }

  # ---------- asymptotic (Kun Zhang / Suzuki Proposition 16) ---------------
  if (method == "asymptotic") {

    # n^3 * HSIC  =  n^2 * scaled_statistic
    stat_asymp <- T_obs_raw * n^3

    lambda_x  <- Re(eigen(Kx_c, symmetric = TRUE)$values)
    lambda_y  <- Re(eigen(Ky_c, symmetric = TRUE)$values)
    lambda_xy <- as.vector(outer(lambda_x, lambda_y, `*`))

    # Null: sum_ij lambda_Xi * lambda_Yj * chi^2_1
    r      <- length(lambda_xy)
    T_null <- replicate(B, sum(lambda_xy * rchisq(r, df = 1)))

    pval <- mean(T_null >= stat_asymp)
    stat <- stat_asymp
    names(stat) <- "n^3 * HSIC"
    meth <- sprintf(
      "HSIC asymptotic test — Zhang/Suzuki (B = %d MC draws, %s kernel)",
      B, kernel_x
    )
  }

  # ---------- bootstrap ----------------------------------------------------
  if (method == "bootstrap") {

    boot_one <- function() {
      ix   <- sample(n, replace = TRUE)
      iy   <- sample(n, replace = TRUE)
      Kxb  <- build_kernel_matrix(x[ix, , drop = FALSE], kernel = kernel_x,
                                  bandwidth = bandwidth_x, degree = degree,
                                  coef0 = coef0)
      Kyb  <- build_kernel_matrix(y[iy, , drop = FALSE], kernel = kernel_y,
                                  bandwidth = bandwidth_y, degree = degree,
                                  coef0 = coef0)
      Kxbc <- center_kernel(Kxb)
      Kybc <- center_kernel(Kyb)
      (sum(diag(Kxbc %*% Kybc)) / n^2) * n
    }

    if (parallelize) {
      cl     <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      T_null <- foreach::foreach(
        iterators::icount(B), .combine = c,
        .export = c("build_kernel_matrix", "center_kernel")
      ) %dopar% boot_one()
      parallel::stopCluster(cl)
    } else {
      T_null <- replicate(B, boot_one())
    }

    pval <- mean(T_null >= T_obs_scaled)
    stat <- T_obs_scaled
    names(stat) <- "HSIC (scaled by n)"
    meth <- sprintf(
      "HSIC bootstrap test (B = %d, %s kernel)", B, kernel_x
    )
  }

  structure(
    list(statistic = stat, p.value = pval, method = meth, data.name = dname),
    class = "htest"
  )
}


#' HSIC-Based Regression Residual Independence Test
#'
#' @description Tests the null hypothesis that the predictors \eqn{X} are
#'   independent of the residuals \eqn{\hat{e} = Y - X\hat{\beta}} from a
#'   fitted linear model.
#'
#' @export
hsic_resid_test <- function(model,
                            x           = NULL,
                            kernel_x    = "gaussian",
                            kernel_y    = kernel_x,
                            bandwidth_x = NULL,
                            bandwidth_y = NULL,
                            degree      = 2,
                            coef0       = 1,
                            B           = 1000,
                            parallelize = FALSE,
                            cores       = 2) {

  X <- model.matrix(model)
  b <- as.matrix(coef(model))
  n <- nobs(model)
  e <- resid(model)

  # Predictor matrix for HSIC (drop intercept by default)
  if (is.null(x)) {
    int_col <- which(colnames(X) == "(Intercept)")
    x <- if (length(int_col)) X[, -int_col, drop = FALSE] else X
  }
  if (!is.matrix(x)) x <- as.matrix(x)
  if (nrow(x) != n)
    stop("'x' must have the same number of rows as observations in 'model'.")

  # Calculate bandwidths ONCE here if NULL, to save repeating in bootstrap loop
  if (is.null(bandwidth_x)) bandwidth_x <- median_bandwidth(x)
  if (is.null(bandwidth_y)) bandwidth_y <- median_bandwidth(e)

  # Observed test statistic (Scaled by n)
  T_obs_raw <- hsic(x, e,
                    kernel_x    = kernel_x,    kernel_y    = kernel_y,
                    bandwidth_x = bandwidth_x, bandwidth_y = bandwidth_y,
                    degree = degree, coef0 = coef0)

  T_obs_scaled <- T_obs_raw * n

  e0 <- e - mean(e)   # centred residuals (dHSIC source convention)

  one_boot <- function() {
    idx_e  <- sample(n, replace = TRUE)
    e_B    <- e0[idx_e]
    idx_x  <- sample(n, replace = TRUE)
    x_B    <- x[idx_x, , drop = FALSE]
    X_B    <- X[idx_x, , drop = FALSE]
    Yhat_B <- X_B %*% b + e_B
    bhat_B <- tryCatch(
      solve(crossprod(X_B), crossprod(X_B, Yhat_B)),
      error = function(err) MASS::ginv(crossprod(X_B)) %*% crossprod(X_B, Yhat_B)
    )
    ehat_B <- Yhat_B - X_B %*% bhat_B

    # Scale bootstrap result by n
    hsic(x_B, ehat_B,
         kernel_x    = kernel_x,    kernel_y    = kernel_y,
         bandwidth_x = bandwidth_x, bandwidth_y = bandwidth_y,
         degree = degree, coef0 = coef0) * n
  }

  if (parallelize) {
    cl     <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    T_null <- foreach::foreach(
      iterators::icount(B), .combine = c,
      .export = c("hsic", "build_kernel_matrix", "center_kernel")
    ) %dopar% one_boot()
    parallel::stopCluster(cl)
  } else {
    T_null <- replicate(B, one_boot())
  }

  pval <- mean(T_null >= T_obs_scaled)
  stat <- T_obs_scaled
  names(stat) <- "HSIC (scaled by n)"

  structure(
    list(
      statistic = stat,
      p.value   = pval,
      method    = sprintf(
        "HSIC regression residual test — dHSIC bootstrap (B = %d, %s kernel)",
        B, kernel_x
      ),
      data.name = deparse(formula(model))
    ),
    class = "htest"
  )
}


# ---- Utilities ---------------------------------------------------------------

#' Median Heuristic Bandwidth for the Gaussian Kernel
#'
#' @description Computes a data-driven bandwidth using the median heuristic.
#'
#' @export
median_bandwidth <- function(x) {
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  d2 <- as.vector(dist(x, method = "euclidean"))^2
  med <- median(d2[d2 > 0])
  if (is.na(med) || med == 0)
    warning("Median bandwidth is zero or NA; returning 1.")
  if (is.na(med) || med == 0) return(1)
  med
}


#' Plot the HSIC Null Distribution
#'
#' @description Visualises the scaled permutation or bootstrap null distribution.
#'
#' @export
hsic_plot <- function(x, y,
                      method      = c("permutation", "bootstrap"),
                      kernel_x    = "gaussian",
                      kernel_y    = kernel_x,
                      bandwidth_x = NULL,
                      bandwidth_y = NULL,
                      degree      = 2,
                      coef0       = 1,
                      B           = 500,
                      alpha       = 0.05,
                      ...) {
  method <- match.arg(method)

  if (is.null(bandwidth_x)) bandwidth_x <- median_bandwidth(x)
  if (is.null(bandwidth_y)) bandwidth_y <- median_bandwidth(y)

  res    <- hsic_test(x, y, method = method,
                      kernel_x = kernel_x, kernel_y = kernel_y,
                      bandwidth_x = bandwidth_x, bandwidth_y = bandwidth_y,
                      degree = degree, coef0 = coef0,
                      B = B)

  # Re-generate null for plotting using scaled values matching hsic_test
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  if (is.vector(y)) y <- matrix(y, ncol = 1)
  n    <- nrow(x)
  Kx   <- build_kernel_matrix(x, kernel = kernel_x, bandwidth = bandwidth_x,
                              degree = degree, coef0 = coef0)
  Ky   <- build_kernel_matrix(y, kernel = kernel_y, bandwidth = bandwidth_y,
                              degree = degree, coef0 = coef0)
  Kx_c <- center_kernel(Kx)
  Ky_c <- center_kernel(Ky)
  T_obs <- (sum(diag(Kx_c %*% Ky_c)) / n^2) * n

  if (method == "permutation") {
    T_null <- replicate(B, {
      idx <- sample(n)
      Kxp <- center_kernel(Kx[idx, idx])
      (sum(diag(Kxp %*% Ky_c)) / n^2) * n
    })
  } else {
    T_null <- replicate(B, {
      ix  <- sample(n, replace = TRUE)
      iy  <- sample(n, replace = TRUE)
      Kxb <- center_kernel(
        build_kernel_matrix(x[ix, , drop=FALSE], kernel=kernel_x,
                            bandwidth=bandwidth_x, degree=degree, coef0=coef0))
      Kyb <- center_kernel(
        build_kernel_matrix(y[iy, , drop=FALSE], kernel=kernel_y,
                            bandwidth=bandwidth_y, degree=degree, coef0=coef0))
      (sum(diag(Kxb %*% Kyb)) / n^2) * n
    })
  }

  crit <- quantile(T_null, 1 - alpha)

  xlim <- range(c(T_null, T_obs, crit))
  plot(density(T_null), xlim = xlim,
       main = "HSIC null distribution",
       xlab = "HSIC (scaled by n)", ylab = "Density", ...)
  abline(v = crit,  col = "red",  lty = 2, lwd = 2)
  abline(v = T_obs, col = "blue", lty = 1, lwd = 2)
  legend("topright",
         legend = c(sprintf("Critical value (alpha = %.2f)", alpha),
                    "Observed HSIC"),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2, bty = "n")

  invisible(list(T_obs = T_obs, T_null = T_null, critical = crit))
}
