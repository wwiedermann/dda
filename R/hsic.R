# ==============================================================================
# hsic.R
#
# Hilbert-Schmidt Independence Criterion (HSIC)
# C++-Accelerated Implementation for R Package
#
# References:
#   Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Schölkopf, B., &
#     Smola, A. J. (2008). A kernel statistical test of independence.
#     Advances in Neural Information Processing Systems, 20.
#
#   Peters, J., Pfister, N., & Mooij, J. M. (2022). dHSIC: Independence
#     Testing via Hilbert Schmidt Independence Criterion. R package version
#     2.1. https://CRAN.R-project.org/package=dHSIC
#
#   Suzuki, T. (2025). Statistical Learning Theory: Kernel Methods,
#     Sparsity, and Related Topics. Springer. Section 4.2,
#     Propositions 11-16.
#   Zhang, K. (eigenvalue-based asymptotic null; cited in Suzuki, 2025,
#     Proposition 16 and Section 4.2).
#
# C++ code: src/hsic.cpp (compiled via Rcpp)
# ==============================================================================


# ==============================================================================
# INTERNAL HELPERS
# ==============================================================================

#' Bandwidth to a Positive Scalar
#'
#' @description Translates the flexible \code{bandwidth} argument into a
#'   single positive numeric ready for the C++ kernel builders.
#'
#' @param x A numeric matrix (coerced from vector if needed).
#' @param bw \code{NULL} or \code{"median"} triggers \code{median_bandwidth};
#'   a strictly positive scalar is used directly.
#'
#' @return A single strictly positive numeric.
#' @keywords internal
#' @noRd
.get_bandwidth <- function(x, bw) {
  if (is.null(bw) || identical(bw, "median"))
    return(median_bandwidth_cpp(x))
  if (is.numeric(bw) && length(bw) == 1L && bw > 0)
    return(bw)
  stop("'bandwidth' must be NULL, \"median\", or a strictly positive numeric scalar.",
       call. = FALSE)
}


# ==============================================================================
# median_bandwidth
# ==============================================================================

#' @title Median Heuristic Bandwidth
#'
#' @description Returns the median of all strictly positive squared
#'   pairwise Euclidean distances, interpreted as the variance
#'   parameter \eqn{\sigma^2} for the Gaussian kernel
#'   \eqn{K(x,y) = \exp(-\|x-y\|^2 / (2\sigma^2))}. This is the same
#'   convention used in the \code{dHSIC} source code (Peters et al.,
#'   2022).
#'
#' @name median_bandwidth
#'
#' @param x A numeric vector or matrix of observations (n rows).
#'
#' @return A single strictly positive numeric: the bandwidth
#'   \eqn{\sigma^2}, ready to pass as \code{bandwidth_x} or
#'   \code{bandwidth_y}. Returns 1 when all pairwise distances are zero.
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., &
#'   Smola, A. J. (2008). A kernel statistical test of independence.
#'   \emph{Advances in Neural Information Processing Systems}, 20.
#'
#' Peters, J., Pfister, N., & Mooij, J. M. (2022). \emph{dHSIC:
#'   Independence Testing via Hilbert Schmidt Independence Criterion}.
#'   R package version 2.1.
#'
#' @seealso \code{\link{hsic.test}}
#'
#' @examples
#' x <- rnorm(80)
#' median_bandwidth(x)
#'
#' @export
median_bandwidth <- function(x) {
  if (is.vector(x)) x <- matrix(x, ncol = 1L)
  bw <- median_bandwidth_cpp(x)
  if (bw == 1.0 && all(dist(x) == 0))
    warning("All pairwise distances are zero; bandwidth defaulting to 1.",
            call. = FALSE)
  bw
}


# ==============================================================================
# hsic
# ==============================================================================

#' @title Compute the Empirical HSIC
#'
#' @description \code{hsic} computes the empirical Hilbert-Schmidt
#'   Independence Criterion between two variables X and Y using the
#'   biased V-statistic
#'   \eqn{(1/n^2) \, \mathrm{tr}(\tilde{K}_X \tilde{K}_Y)}.
#'   All kernel and centering computations are performed in C++.
#'
#' @name hsic
#'
#' @param x           A numeric vector of length n or matrix (n x p).
#' @param y           A numeric vector of length n or matrix (n x q).
#' @param kernel_x    Kernel for X. One of \code{c("gaussian",
#'   "laplace", "linear", "polynomial")}. Default is
#'   \code{"gaussian"}.
#' @param kernel_y    Kernel for Y. Defaults to \code{kernel_x}.
#' @param bandwidth_x Bandwidth for the X kernel. The median heuristic
#'   (\code{\link{median_bandwidth}}) is always the default and is
#'   applied when \code{bandwidth_x = NULL} (default) or
#'   \code{bandwidth_x = "median"}. Alternatively, a strictly positive
#'   numeric value is used directly as \eqn{\sigma^2} for the Gaussian
#'   kernel \eqn{K(x,y) = \exp(-\|x-y\|^2 / (2\sigma^2))}.
#' @param bandwidth_y Bandwidth for the Y kernel. Same options as
#'   \code{bandwidth_x}; the median heuristic is the default.
#' @param degree      Integer degree for the polynomial kernel. Default
#'   \code{2}.
#' @param coef0       Constant term for the polynomial kernel. Default
#'   \code{1}.
#'
#' @details
#' The Gaussian kernel uses \eqn{K(x,y) = \exp(-\|x-y\|^2 / (2\sigma^2))}
#' where \code{bandwidth} = \eqn{\sigma^2}. The median heuristic sets
#' \eqn{\sigma^2} to the median of all strictly positive squared pairwise
#' distances, matching the convention in the \pkg{dHSIC} package.
#'
#' @return A single non-negative numeric value: the raw HSIC estimate
#'   \eqn{(1/n^2) \, \mathrm{tr}(\tilde{K}_X \tilde{K}_Y)}. A value of
#'   zero (with a characteristic kernel) implies independence.
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., &
#'   Smola, A. J. (2008). A kernel statistical test of independence.
#'   \emph{Advances in Neural Information Processing Systems}, 20.
#'
#' Peters, J., Pfister, N., & Mooij, J. M. (2022). \emph{dHSIC:
#'   Independence Testing via Hilbert Schmidt Independence Criterion}.
#'   R package version 2.1.
#'   \url{https://CRAN.R-project.org/package=dHSIC}
#'
#' @seealso \code{\link{hsic.test}}, \code{\link{median_bandwidth}}
#'
#' @examples
#' set.seed(12)
#' x <- rnorm(100)
#' hsic(x, rnorm(100))
#' hsic(x, x + rnorm(100, sd = 0.5))
#' hsic(x, rnorm(100), bandwidth_x = 1, bandwidth_y = 1)
#'
#' @export
hsic <- function(
  x,
  y,
  kernel_x    = "gaussian",
  kernel_y    = kernel_x,
  bandwidth_x = NULL,
  bandwidth_y = NULL,
  degree      = 2L,
  coef0       = 1
) {
  x <- as.matrix(x); y <- as.matrix(y)
  n <- nrow(x)
  if (nrow(y) != n)
    stop("x and y must have the same number of rows.", call. = FALSE)

  bwx <- .get_bandwidth(x, bandwidth_x)
  bwy <- .get_bandwidth(y, bandwidth_y)

  Kx <- build_kernel_matrix_cpp(x, kernel_x, bwx, degree, coef0)
  Ky <- build_kernel_matrix_cpp(y, kernel_y, bwy, degree, coef0)

  Kx_c <- center_kernel_cpp(Kx)
  Ky_c <- center_kernel_cpp(Ky)

  hsic_trace_cpp(Kx_c, Ky_c)
}


# ==============================================================================
# hsic.test
# ==============================================================================

#' @title HSIC Independence Test
#'
#' @description \code{hsic.test} tests whether X and Y are independent
#'   using the Hilbert-Schmidt Independence Criterion. The test statistic
#'   is \eqn{n \times \widehat{\mathrm{HSIC}}}, labelled \code{"HSIC"},
#'   and is consistent across all four inference methods.
#'
#'   Available null-distribution methods:
#'   \describe{
#'     \item{\code{"gamma"}}{Fits a Gamma distribution to the first two
#'       null moments of HSIC analytically (Peters, Pfister & Mooij,
#'       2022). No resampling required.}
#'     \item{\code{"permutation"}}{Permutes the row index of X to
#'       simulate the null distribution entirely in C++.}
#'     \item{\code{"eigenvalue"}}{Derives the null distribution from the
#'       eigenvalue spectrum of the centered kernel matrices (Zhang, 2011).
#'      More accurate in small samples; uses \code{B} Monte
#'       Carlo draws. Reports \eqn{n \times \widehat{\mathrm{HSIC}}} with
#'       p-value from the spectral null.}
#'     \item{\code{"bootstrap"}}{Independently resamples rows of the
#'       kernel matrices with replacement in C++ (Peters, Pfister &
#'       Mooij, 2022).}
#'   }
#'
#' @name hsic.test
#'
#' @param x           A numeric vector of length n or matrix (n x p).
#' @param y           A numeric vector of length n or matrix (n x q).
#' @param method      Inference method for the null distribution. One of
#'   \code{c("gamma", "permutation", "eigenvalue", "bootstrap")}.
#'   Default is \code{"gamma"}.
#' @param kernel_x    Kernel for X. One of \code{c("gaussian",
#'   "laplace", "linear", "polynomial")}. Default is
#'   \code{"gaussian"}.
#' @param kernel_y    Kernel for Y. Defaults to \code{kernel_x}.
#' @param bandwidth_x Bandwidth (\eqn{\sigma^2}) for the X kernel. The
#'   median heuristic is always the default and is applied when
#'   \code{bandwidth_x = NULL}. A strictly positive numeric value is
#'   used directly.
#' @param bandwidth_y Bandwidth for the Y kernel. Same options as
#'   \code{bandwidth_x}; the median heuristic is the default.
#' @param degree      Integer degree for the polynomial kernel. Default
#'   \code{2}.
#' @param coef0       Constant term for the polynomial kernel. Default
#'   \code{1}.
#' @param B           Number of permutation or bootstrap replicates, or
#'   Monte Carlo draws for \code{method = "eigenvalue"}. Ignored for
#'   \code{method = "gamma"}. Default is \code{1000}.
#'
#' @details
#' All four methods report \eqn{n \times \widehat{\mathrm{HSIC}}} as
#' the test statistic. Permutation and bootstrap p-values use the
#' Laplace correction
#' \eqn{(\#\{T_b \geq T_{\mathrm{obs}}\} + 1) / (B + 1)}.
#'
#' @return An object of class \code{"htest"} with components:
#'   \describe{
#'     \item{\code{statistic}}{The test statistic
#'       \eqn{n \times \widehat{\mathrm{HSIC}}}, labelled
#'       \code{"HSIC"}.}
#'     \item{\code{p.value}}{P-value for the test of independence.}
#'     \item{\code{bandwidths}}{Resolved bandwidths \code{c(x, y)}.}
#'   }
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., &
#'   Smola, A. J. (2008). A kernel statistical test of independence.
#'   \emph{Advances in Neural Information Processing Systems}, 20.
#'
#' Peters, J., Pfister, N., & Mooij, J. M. (2022). \emph{dHSIC:
#'   Independence Testing via Hilbert Schmidt Independence Criterion}.
#'   R package version 2.1.
#'   \url{https://CRAN.R-project.org/package=dHSIC}
#'
#' Zhang, K., Peters, J., Janzing, D., & Scholkopf, B. (2011).
#'   Kernel-based conditional independence test and application in
#'   causal discovery. In \emph{Proceedings of the Twenty-Seventh
#'   Conference on Uncertainty in Artificial Intelligence (UAI 2011)}
#'   (pp. 804-813).
#'
#' @seealso \code{\link{hsic}}, \code{\link{median_bandwidth}}
#'
#' @examples
#' set.seed(7)
#' n <- 80
#' x <- rnorm(n)
#' hsic.test(x, rnorm(n), method = "gamma")
#' hsic.test(x, rnorm(n), method = "permutation", B = 499)
#' hsic.test(x, x + rnorm(n), method = "permutation", B = 499)
#'
#' @export
hsic.test <- function(
  x,
  y,
  method      = c("gamma", "permutation", "eigenvalue", "bootstrap"),
  kernel_x    = "gaussian",
  kernel_y    = kernel_x,
  bandwidth_x = NULL,
  bandwidth_y = NULL,
  degree      = 2L,
  coef0       = 1,
  B           = 1000L
) {
  method <- match.arg(method)
  x <- as.matrix(x); y <- as.matrix(y)
  n <- nrow(x)
  if (nrow(y) != n)
    stop("x and y must have the same number of rows.", call. = FALSE)
  if (B < 1L)
    stop("B must be a positive integer.", call. = FALSE)

  bwx <- .get_bandwidth(x, bandwidth_x)
  bwy <- .get_bandwidth(y, bandwidth_y)

  Kx <- build_kernel_matrix_cpp(x, kernel_x, bwx, degree, coef0)
  Ky <- build_kernel_matrix_cpp(y, kernel_y, bwy, degree, coef0)

  Kx_c <- center_kernel_cpp(Kx)
  Ky_c <- center_kernel_cpp(Ky)

  T_obs   <- hsic_trace_cpp(Kx_c, Ky_c)
  T_obs_n <- T_obs * n

  # ---- Gamma approximation (Peters & Mooij, 2022)
  if (method == "gamma") {
    ea <- c(sum(Kx) / n^2, sum(Ky) / n^2)
    eb <- c(sum(Kx^2) / n^2, sum(Ky^2) / n^2)
    ec <- c(sum(rowSums(Kx)^2) / n^3, sum(rowSums(Ky)^2) / n^3)

    prod_a <- ea[1] * ea[2]; prod_b <- eb[1] * eb[2]; prod_c <- ec[1] * ec[2]
    exp_est <- (1 - (ea[1] + ea[2]) + prod_a) / n

    t1 <- prod_b
    t2 <- prod_a^2
    t3 <- 2 * prod_c
    t4 <- eb[1] * ea[2]^2 + eb[2] * ea[1]^2
    t5 <- -2 * (eb[1] * ec[2] + eb[2] * ec[1])
    t6 <- -2 * (ec[1] * ea[2]^2 + ec[2] * ea[1]^2)
    t7 <- 2 * ec[1] * ec[2]

    f1 <- (n - 4) * (n - 5)
    f2 <- n * (n - 1) * (n - 2) * (n - 3)
    var_est <- 2 * f1 / f2 * (t1 + t2 + t3 + t4 + t5 + t6 + t7)

    alpha_g <- exp_est^2 / var_est
    beta_g  <- n * var_est / exp_est
    pval    <- stats::pgamma(T_obs_n, shape = alpha_g, scale = beta_g, lower.tail = FALSE)

    stat <- T_obs_n; names(stat) <- "HSIC"
    meth <- sprintf("HSIC gamma test (%s kernel; Peters & Mooij 22)", kernel_x)
  }

  # ---- Permutation test
  if (method == "permutation") {
    T_null <- permute_hsic_cpp(Kx_c, Ky_c, B)
    pval   <- (sum(T_null >= T_obs) + 1L) / (B + 1L)

    stat <- T_obs_n; names(stat) <- "HSIC"
    meth <- sprintf("HSIC permutation test (B = %d, %s kernel)", B, kernel_x)
  }

  # ---- Eigenvalue / spectral (Zhang; reports n*HSIC, p-value from eigenvalue null)
  if (method == "eigenvalue") {
    stat_sc <- T_obs * n^3

    lam_x  <- Re(eigen(Kx_c, symmetric = TRUE)$values)
    lam_y  <- Re(eigen(Ky_c, symmetric = TRUE)$values)
    lam_xy <- as.vector(outer(lam_x, lam_y, `*`))

    T_null <- replicate(B, sum(lam_xy * stats::rchisq(length(lam_xy), df = 1)))
    pval   <- mean(T_null >= stat_sc)

    stat <- T_obs_n; names(stat) <- "HSIC"
    meth <- sprintf("HSIC eigenvalue test (B = %d MC draws, %s kernel; Zhang method)",
                    B, kernel_x)
  }

  # ---- Bootstrap test
  if (method == "bootstrap") {
    T_null <- bootstrap_hsic_cpp(Kx, Ky, B)
    pval   <- (sum(T_null >= T_obs) + 1L) / (B + 1L)

    stat <- T_obs_n; names(stat) <- "HSIC"
    meth <- sprintf("HSIC bootstrap test (B = %d, %s kernel)", B, kernel_x)
  }

  structure(list(
    statistic = stat,
    p.value   = pval,
    method    = meth,
    data.name = paste(deparse(substitute(x)), "and", deparse(substitute(y))),
    bandwidths = c(bandwidth_x = bwx, bandwidth_y = bwy)
  ), class = "htest")
}
