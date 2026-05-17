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
# Dependencies: Rcpp (all kernel and null-distribution work is done in C++)
# ==============================================================================


# ==============================================================================
# 1. C++ BACKEND
# ==============================================================================

.hsic_cpp_code <- "
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// median_bandwidth_cpp
//
// Returns the median of all strictly positive squared pairwise Euclidean
// distances, interpreted as the variance parameter sigma^2 for the Gaussian
// kernel  K(x,y) = exp( -||x-y||^2 / (2*sigma^2) ).
//
// This matches the dHSIC package convention (Peters et al., 2022):
//   gaussian_grammat_rcpp uses bandwidth = sigma^2, and
//   median_bandwidth_rcpp returns median(d_ij^2) as sigma^2.
//
// Uses nth_element for O(n^2) average-case cost instead of a full sort.
// Returns 1.0 when all pairwise distances are zero (identical observations).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double median_bandwidth_cpp(NumericMatrix x) {
    int n = x.nrow();
    int p = x.ncol();
    std::vector<double> dists;
    dists.reserve((size_t)(n) * (n - 1) / 2);

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d2 = 0.0;
            for (int k = 0; k < p; k++) {
                double diff = x(i, k) - x(j, k);
                d2 += diff * diff;
            }
            if (d2 > 0.0) dists.push_back(d2);
        }
    }

    if (dists.empty()) return 1.0;

    size_t mid = dists.size() / 2;
    std::nth_element(dists.begin(), dists.begin() + mid, dists.end());
    double med = dists[mid];

    // Even-length median: average two middle elements
    if (dists.size() % 2 == 0 && dists.size() > 1) {
        std::nth_element(dists.begin(), dists.begin() + mid - 1, dists.end());
        med = (med + dists[mid - 1]) / 2.0;
    }

    return (med == 0.0) ? 1.0 : med;
}

// ---------------------------------------------------------------------------
// build_kernel_matrix_cpp
//
// Constructs an n x n kernel (Gram) matrix.  Supported kernels:
//
//   'gaussian'   : K(x,y) = exp( -||x-y||^2 / (2*bandwidth) )
//                  bandwidth = sigma^2.  With bandwidth = 1: exp(-d^2/2).
//                  Matches dHSIC's gaussian_grammat_rcpp convention
//                  (Peters et al., 2022).
//
//   'laplace'    : K(x,y) = exp( -||x-y|| / bandwidth )
//                  bandwidth = sigma (length scale).
//
//   'linear'     : K(x,y) = x'y
//
//   'polynomial' : K(x,y) = (x'y + coef0)^degree
//
// Only the upper triangle is computed; result is symmetrised in-place.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix build_kernel_matrix_cpp(NumericMatrix x,
                                      std::string kernel,
                                      double bandwidth,
                                      int degree,
                                      double coef0) {
    int n = x.nrow();
    int p = x.ncol();
    NumericMatrix K(n, n);

    bool is_gauss  = (kernel == \"gaussian\");
    bool is_lap    = (kernel == \"laplace\");
    bool is_linear = (kernel == \"linear\");
    bool is_poly   = (kernel == \"polynomial\");

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double val = 0.0;

            if (is_linear || is_poly) {
                for (int k = 0; k < p; k++) val += x(i, k) * x(j, k);
                if (is_poly) val = std::pow(val + coef0, (double)degree);
            } else {
                double d2 = 0.0;
                for (int k = 0; k < p; k++) {
                    double diff = x(i, k) - x(j, k);
                    d2 += diff * diff;
                }
                if (is_gauss) {
                    val = std::exp(-d2 / (2.0 * bandwidth));   // sigma^2 convention
                } else {  // laplace
                    val = std::exp(-std::sqrt(d2) / bandwidth);
                }
            }
            K(i, j) = val;
            K(j, i) = val;
        }
    }
    return K;
}

// ---------------------------------------------------------------------------
// center_kernel_cpp
//
// Applies the centering transform  K_c = H K H,  H = I - (1/n)*11'.
// Uses the O(n^2) row-mean identity:
//   K_c(i,j) = K(i,j) - row_mean_i - row_mean_j + grand_mean
// avoiding the O(n^3) triple matrix product.
// Reference: Suzuki (2025), Equation (4.6).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix center_kernel_cpp(NumericMatrix K) {
    int n = K.nrow();
    NumericMatrix Kc(n, n);
    std::vector<double> row_means(n, 0.0);
    double grand_mean = 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) row_means[i] += K(i, j);
        row_means[i] /= n;
        grand_mean   += row_means[i];
    }
    grand_mean /= n;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Kc(i, j) = K(i, j) - row_means[i] - row_means[j] + grand_mean;

    return Kc;
}

// ---------------------------------------------------------------------------
// hsic_trace_cpp
//
// Computes (1/n^2) * trace(Kx_c * Ky_c) for symmetric centered matrices.
// Uses trace(A*B) = sum_{ij} A_{ij} B_{ij} (valid when B is symmetric),
// reducing O(n^3) matrix multiply + trace to an O(n^2) element-wise sum.
// Reference: Suzuki (2025), Proposition 15.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double hsic_trace_cpp(NumericMatrix Kx_c, NumericMatrix Ky_c) {
    int n = Kx_c.nrow();
    double s = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            s += Kx_c(i, j) * Ky_c(i, j);
    return s / ((double)n * n);
}

// ---------------------------------------------------------------------------
// permute_hsic_cpp
//
// Generates B permutation-null HSIC values.
// Since H commutes with permutation matrices, permuting the already-centered
// Kx_c is equivalent to permuting raw Kx and then centering.
// Reference: Suzuki (2025), Section 4.2.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector permute_hsic_cpp(NumericMatrix Kx_c, NumericMatrix Ky_c, int B) {
    int n = Kx_c.nrow();
    NumericVector null_dist(B);
    IntegerVector idx = seq(0, n - 1);

    for (int b = 0; b < B; b++) {
        IntegerVector p = sample(idx, n, false);
        double val = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                val += Kx_c(p[i], p[j]) * Ky_c(i, j);
        null_dist[b] = val / ((double)n * n);
    }
    return null_dist;
}

// ---------------------------------------------------------------------------
// bootstrap_hsic_cpp
//
// Generates B bootstrap-null HSIC values by independently resampling rows of
// Kx and Ky with replacement and re-centering.
// Mirrors the dHSIC bootstrap strategy (Peters et al., 2022).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector bootstrap_hsic_cpp(NumericMatrix Kx, NumericMatrix Ky, int B) {
    int n = Kx.nrow();
    NumericVector null_dist(B);
    IntegerVector seq_idx = seq(0, n - 1);

    for (int b = 0; b < B; b++) {
        IntegerVector ix = sample(seq_idx, n, true);
        IntegerVector iy = sample(seq_idx, n, true);

        NumericMatrix Kx_b(n, n), Ky_b(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                Kx_b(i, j) = Kx(ix[i], ix[j]);
                Ky_b(i, j) = Ky(iy[i], iy[j]);
            }

        // Re-center resampled submatrices inline
        std::vector<double> rx(n, 0.0), ry(n, 0.0);
        double gx = 0.0, gy = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                rx[i] += Kx_b(i, j);
                ry[i] += Ky_b(i, j);
            }
            rx[i] /= n; ry[i] /= n;
            gx += rx[i]; gy += ry[i];
        }
        gx /= n; gy /= n;

        double val = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double kxc = Kx_b(i, j) - rx[i] - rx[j] + gx;
                double kyc = Ky_b(i, j) - ry[i] - ry[j] + gy;
                val += kxc * kyc;
            }
        null_dist[b] = val / ((double)n * n);
    }
    return null_dist;
}
"

Rcpp::sourceCpp(code = .hsic_cpp_code)


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
.get_bandwidth <- function(x, bw) {
  if (is.null(bw) || identical(bw, "median"))
    return(median_bandwidth_cpp(x))
  if (is.numeric(bw) && length(bw) == 1L && bw > 0)
    return(bw)
  stop("'bandwidth' must be NULL, \"median\", or a strictly positive numeric scalar.",
       call. = FALSE)
}


# ------------------------------------------------------------------------------
# median_bandwidth
# ------------------------------------------------------------------------------

#' Median Heuristic Bandwidth
#'
#' @description Returns the median of all strictly positive squared pairwise
#'   Euclidean distances, interpreted as the variance parameter
#'   \eqn{\sigma^2} for the Gaussian kernel
#'   \eqn{K(x,y) = \exp(-\|x-y\|^2 / (2\sigma^2))}.
#'   This is the same convention used in the \pkg{dHSIC} source code
#'   (Peters et al., 2022).
#'
#' @param x A numeric vector or matrix of observations (\eqn{n} rows).
#'
#' @return A single strictly positive numeric: the bandwidth \eqn{\sigma^2},
#'   ready to pass as \code{bandwidth_x} or \code{bandwidth_y}.
#'   Returns \code{1} (with a warning) when all pairwise distances are zero.
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., &
#'   Smola, A. J. (2008). A kernel statistical test of independence.
#'   \emph{Advances in Neural Information Processing Systems}, 20.
#'
#' Peters, J., Pfister, N., & Mooij, J. M. (2022). \emph{dHSIC: Independence
#'   Testing via Hilbert Schmidt Independence Criterion}. R package version
#'   2.1.
#'
#' @keywords internal
median_bandwidth <- function(x) {
  if (is.vector(x)) x <- matrix(x, ncol = 1L)
  bw <- median_bandwidth_cpp(x)
  if (bw == 1.0 && all(dist(x) == 0))
    warning("All pairwise distances are zero; bandwidth defaulting to 1.", call. = FALSE)
  bw
}


# ------------------------------------------------------------------------------
# hsic
# ------------------------------------------------------------------------------

#' Compute the Empirical HSIC
#'
#' @description Computes the empirical Hilbert-Schmidt Independence Criterion
#'   (HSIC) between two random variables \eqn{X} and \eqn{Y} using the
#'   biased V-statistic estimator (Gretton et al., 2008):
#'   \deqn{
#'     \widehat{\mathrm{HSIC}}(X, Y)
#'     = \frac{1}{n^2}\operatorname{tr}(\widetilde{K}_X \widetilde{K}_Y).
#'   }
#'   This is algebraically identical to the \code{dHSIC} three-term formula
#'   (Peters et al., 2022).  All kernel and centering work is done in C++.
#'
#' @param x A numeric vector (length \eqn{n}) or matrix (\eqn{n \times p}).
#' @param y A numeric vector (length \eqn{n}) or matrix (\eqn{n \times q}).
#' @param kernel_x Kernel for \eqn{X}: \code{"gaussian"} (default),
#'   \code{"laplace"}, \code{"linear"}, or \code{"polynomial"}.
#' @param kernel_y Kernel for \eqn{Y}.  Defaults to \code{kernel_x}.
#' @param bandwidth_x Bandwidth for the \eqn{X} kernel.  \code{NULL} (default)
#'   or \code{"median"} uses the median heuristic (\code{\link{median_bandwidth}}).
#'   A strictly positive numeric is used directly as \eqn{\sigma^2} for the
#'   Gaussian kernel \eqn{K(x,y) = \exp(-\|x-y\|^2 / (2\sigma^2))}, matching
#'   the \pkg{dHSIC} parameterisation (Peters et al., 2022).
#' @param bandwidth_y Bandwidth for the \eqn{Y} kernel.  Same options as
#'   \code{bandwidth_x}.
#' @param degree Integer degree for the polynomial kernel.  Default \code{2}.
#' @param coef0 Constant for the polynomial kernel.  Default \code{1}.
#'
#' @return A single non-negative numeric: raw HSIC \eqn{= (1/n^2)\,\mathrm{tr}(\widetilde{K}_X \widetilde{K}_Y)}.
#'   Zero (with a characteristic kernel) implies independence (Suzuki, 2025,
#'   Proposition 12).
#'
#' @details
#' \subsection{Kernel and bandwidth convention}{
#'   The Gaussian kernel uses \eqn{K(x,y) = \exp(-\|x-y\|^2 / (2\sigma^2))}
#'   where \code{bandwidth} = \eqn{\sigma^2}.  This matches the
#'   \pkg{dHSIC} source (Peters et al., 2022).  With \code{bandwidth = 1}
#'   both functions compute \eqn{\exp(-d^2/2)}.  The median heuristic sets
#'   \eqn{\sigma^2 = \mathrm{median}\{\|x_i - x_j\|^2 : i \neq j\}},
#'   also matching \pkg{dHSIC}.
#' }
#'
#' @seealso \code{\link{hsic_test}}, \code{\link{median_bandwidth}}
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., &
#'   Smola, A. J. (2008). A kernel statistical test of independence.
#'   \emph{Advances in Neural Information Processing Systems}, 20.
#'
#' Suzuki, T. (2025). \emph{Statistical Learning Theory}. Springer.
#'   Propositions 12, 13, 15.
#'
#' Peters, J., Pfister, N., & Mooij, J. M. (2022). \emph{dHSIC: Independence
#'   Testing via Hilbert Schmidt Independence Criterion}. R package version
#'   2.1.
#'
#' @examples
#' set.seed(12)
#' x <- rnorm(100)
#' hsic(x, rnorm(100))             # near zero (independent)
#' hsic(x, x + rnorm(100, sd=0.5)) # positive (dependent)
#' hsic(x, rnorm(100), bandwidth_x = 1, bandwidth_y = 1)  # fixed bw, matches dHSIC
#'
#' @export
hsic <- function(x, y,
                 kernel_x    = "gaussian",
                 kernel_y    = kernel_x,
                 bandwidth_x = NULL,
                 bandwidth_y = NULL,
                 degree      = 2L,
                 coef0       = 1) {

  kernel_x <- match.arg(kernel_x, c("gaussian", "laplace", "linear", "polynomial"))
  kernel_y <- match.arg(kernel_y, c("gaussian", "laplace", "linear", "polynomial"))

  if (is.vector(x)) x <- matrix(x, ncol = 1L)
  if (is.vector(y)) y <- matrix(y, ncol = 1L)
  n <- nrow(x)
  if (nrow(y) != n)
    stop("'x' and 'y' must have the same number of observations.", call. = FALSE)

  bw_x <- if (kernel_x %in% c("gaussian", "laplace")) .get_bandwidth(x, bandwidth_x) else 1
  bw_y <- if (kernel_y %in% c("gaussian", "laplace")) .get_bandwidth(y, bandwidth_y) else 1

  Kx_c <- center_kernel_cpp(build_kernel_matrix_cpp(x, kernel_x, bw_x, degree, coef0))
  Ky_c <- center_kernel_cpp(build_kernel_matrix_cpp(y, kernel_y, bw_y, degree, coef0))

  hsic_trace_cpp(Kx_c, Ky_c)
}


# ------------------------------------------------------------------------------
# hsic_test
# ------------------------------------------------------------------------------

#' HSIC Independence Test
#'
#' @description Tests whether \eqn{X \perp\!\!\!\perp Y} using HSIC.
#'   The test statistic is \eqn{n \times \widehat{\mathrm{HSIC}}}, matching
#'   the \code{dHSIC::dhsic.test()} convention (\code{stat = dHSIC * len}),
#'   so printed statistics are directly comparable.
#'
#'   Four null-distribution methods are available:
#'   \describe{
#'     \item{\code{"gamma"}}{Exact Gamma approximation matching the
#'       \pkg{dHSIC} source moment formula (Peters et al., 2022).  Fits a
#'       \eqn{\Gamma(\alpha, \beta)} to the first two moments of the \eqn{n \times
#'       \widehat{\mathrm{HSIC}}} null distribution, computed from the
#'       \emph{uncentered} Gram matrices via sufficient statistics
#'       \eqn{\hat{a}_j = n^{-2}\sum K_j},
#'       \eqn{\hat{b}_j = n^{-2}\sum K_j^2},
#'       \eqn{\hat{c}_j = n^{-3}\sum_i (\sum_k K_j(i,k))^2}.
#'       Instantaneous; no replications required.}
#'     \item{\code{"permutation"}}{Permutes indices of \eqn{x} in C++ (no R
#'       loop).  Exact conditional level (Suzuki, 2025, Section 4.2).}
#'     \item{\code{"eigenvalue"}}{asymptotic spectral mixture: under
#'       \eqn{H_0}, \eqn{n^3\,\widehat{\mathrm{HSIC}} \xrightarrow{d}
#'       \sum_{i,j}\lambda_{X,i}\lambda_{Y,j}\chi^2_1} (Zhang; Suzuki, 2025,
#'       Proposition 16).  \code{B} Monte Carlo draws estimate the p-value.}
#'     \item{\code{"bootstrap"}}{Independently resamples rows of \eqn{K_X}
#'       and \eqn{K_Y} with replacement in C++ (Peters et al., 2022).}
#'   }
#'
#' @param x A numeric vector or matrix (\eqn{n \times p}).
#' @param y A numeric vector or matrix (\eqn{n \times q}).
#' @param method One of \code{"gamma"} (default), \code{"permutation"},
#'   \code{"eigenvalue"}, \code{"bootstrap"}.
#' @param kernel_x Kernel for \eqn{X}: \code{"gaussian"} (default),
#'   \code{"laplace"}, \code{"linear"}, \code{"polynomial"}.
#' @param kernel_y Kernel for \eqn{Y}.  Defaults to \code{kernel_x}.
#' @param bandwidth_x Bandwidth (\eqn{\sigma^2}) for the \eqn{X} kernel.
#'   \code{NULL} (default) uses the median heuristic.  See \code{\link{hsic}}.
#' @param bandwidth_y Bandwidth for the \eqn{Y} kernel.
#' @param degree Polynomial degree.  Default \code{2}.
#' @param coef0 Polynomial constant.  Default \code{1}.
#' @param B Permutation / bootstrap replications, or MC draws for
#'   \code{"eigenvalue"}.  Ignored for \code{"gamma"}.  Default \code{1000}.
#'
#' @return An object of class \code{"htest"}:
#'   \describe{
#'     \item{\code{statistic}}{For \code{"gamma"}, \code{"permutation"},
#'       and \code{"bootstrap"}: \eqn{n \times \widehat{\mathrm{HSIC}}},
#'       matching \code{dHSIC::dhsic.test()$statistic}.
#'       For \code{"eigenvalue"}: \eqn{n^3 \times \widehat{\mathrm{HSIC}}}.}
#'     \item{\code{p.value}}{For \code{"permutation"} and \code{"bootstrap"}:
#'       \eqn{(\#\{T_b \geq T_{\mathrm{obs}}\} + 1)/(B + 1)}, the Laplace
#'       correction used in \code{dHSIC::dhsic.test()}.}
#'     \item{\code{bandwidths}}{Resolved \code{c(x, y)} bandwidths.}
#'   }
#'
#' @seealso \code{\link{hsic}}, \code{\link{hsic_resid_test}},
#'   \code{\link{median_bandwidth}}
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., &
#'   Smola, A. J. (2007/2008). A kernel statistical test of independence.
#'   \emph{Advances in Neural Information Processing Systems}, 20.
#'
#' Suzuki, T. (2025). \emph{Statistical Learning Theory}. Springer.
#'   Section 4.2, Propositions 12, 13, 16.
#'
#' Zhang, K. (eigenvalue-based asymptotic null; cited in Suzuki, 2025,
#'   Proposition 16).
#'
#' Peters, J., Pfister, N., & Mooij, J. M. (2022). \emph{dHSIC}. R package
#'   version 2.1.
#'
#' @examples
#' set.seed(7); n <- 80; x <- rnorm(n)
#' hsic_test(x, rnorm(n), method = "gamma")
#' hsic_test(x, rnorm(n), method = "permutation", B = 499)
#' hsic_test(x, x + rnorm(n), method = "permutation", B = 499)
#'
#' @export
hsic_test <- function(x, y,
                      method      = c("gamma", "permutation", "eigenvalue", "bootstrap"),
                      kernel_x    = "gaussian",
                      kernel_y    = kernel_x,
                      bandwidth_x = NULL,
                      bandwidth_y = NULL,
                      degree      = 2L,
                      coef0       = 1,
                      B           = 1000L) {

  method   <- match.arg(method)
  kernel_x <- match.arg(kernel_x, c("gaussian", "laplace", "linear", "polynomial"))
  kernel_y <- match.arg(kernel_y, c("gaussian", "laplace", "linear", "polynomial"))
  dname    <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  if (is.vector(x)) x <- matrix(x, ncol = 1L)
  if (is.vector(y)) y <- matrix(y, ncol = 1L)
  n <- nrow(x)
  if (nrow(y) != n)
    stop("'x' and 'y' must have the same number of observations.", call. = FALSE)

  bw_x <- if (kernel_x %in% c("gaussian", "laplace")) .get_bandwidth(x, bandwidth_x) else 1
  bw_y <- if (kernel_y %in% c("gaussian", "laplace")) .get_bandwidth(y, bandwidth_y) else 1

  Kx   <- build_kernel_matrix_cpp(x, kernel_x, bw_x, degree, coef0)
  Ky   <- build_kernel_matrix_cpp(y, kernel_y, bw_y, degree, coef0)
  Kx_c <- center_kernel_cpp(Kx)
  Ky_c <- center_kernel_cpp(Ky)

  T_obs   <- hsic_trace_cpp(Kx_c, Ky_c)   # raw HSIC = (1/n^2) tr(Kx_c Ky_c)
  T_obs_n <- T_obs * n                      # n * HSIC, matching dHSIC statistic

  # ---- Gamma (exact dHSIC moment formula, d=2 case) --------------------------
  # Reference: Peters et al. (2022) dHSIC source, dhsic.test(), method="gamma"
  # ---------------------------------------------------------------------------
  if (method == "gamma") {
    # Sufficient statistics on uncentered kernels (dHSIC est.a, est.b, est.c)
    ea <- c(sum(Kx) / n^2,         sum(Ky) / n^2)          # E[K] per variable
    eb <- c(sum(Kx^2) / n^2,       sum(Ky^2) / n^2)        # E[K^2]
    ec <- c(sum(rowSums(Kx)^2) / n^3, sum(rowSums(Ky)^2) / n^3)  # E[row-sum^2]

    prod_a <- ea[1] * ea[2]
    prod_b <- eb[1] * eb[2]
    prod_c <- ec[1] * ec[2]

    # E[HSIC] under H0  (dHSIC: exp.est = (1 - sum(oneoutprod.a) + (d-1)*prod.a)/len)
    # For d=2: oneoutprod.a = c(ea[2], ea[1])
    exp_est <- (1 - (ea[1] + ea[2]) + prod_a) / n

    # Var[HSIC] under H0: seven terms from dHSIC source for d=2 ---------------
    # t1 = prod.b
    t1 <- prod_b
    # t2 = (d-1)^2 * prod.d = prod_a^2  (d=2, prod.d = prod_a^2)
    t2 <- prod_a^2
    # t3 = 2*(d-1)*prod.c = 2*prod_c
    t3 <- 2 * prod_c
    # t4: loop r=1..d-1 then add r=d term
    #   = eb[1]*ea[2]^2 + eb[2]*ea[1]^2
    t4 <- eb[1]*ea[2]^2 + eb[2]*ea[1]^2
    # t5 = -2*(eb[1]*ec[2] + eb[2]*ec[1])
    t5 <- -2 * (eb[1]*ec[2] + eb[2]*ec[1])
    # t6 = -2*(d-1)*(ec[1]*ea[2]^2 + ec[2]*ea[1]^2)  [d=2 => factor -2]
    t6 <- -2 * (ec[1]*ea[2]^2 + ec[2]*ea[1]^2)
    # t7: double loop r<s: 2*ec[1]*ec[2]*ea[2]^2/ea[2]^2 = 2*ec[1]*ec[2]
    t7 <- 2 * ec[1] * ec[2]

    # Factorial correction for d=2:
    #   factor1 = (n-4)*(n-5)
    #   factor2 = n*(n-1)*(n-2)*(n-3)
    f1 <- (n - 4) * (n - 5)
    f2 <- n * (n - 1) * (n - 2) * (n - 3)

    var_est <- 2 * f1 / f2 * (t1 + t2 + t3 + t4 + t5 + t6 + t7)

    if (var_est <= 0 || exp_est <= 0)
      stop("Gamma moment estimates are non-positive; try a different kernel or larger n.",
           call. = FALSE)

    # Gamma shape/scale for the n*HSIC statistic (dHSIC: a, b = len*var/exp)
    alpha_g <- exp_est^2 / var_est
    beta_g  <- n * var_est / exp_est
    pval    <- stats::pgamma(T_obs_n, shape = alpha_g, scale = beta_g,
                             lower.tail = FALSE)

    stat <- T_obs_n; names(stat) <- "n * HSIC"
    meth <- sprintf("HSIC gamma test (%s kernel; dHSIC formula, Peters et al., 2022)",
                    kernel_x)
  }

  # ---- Permutation (C++ loop) -----------------------------------------------
  if (method == "permutation") {
    T_null <- permute_hsic_cpp(Kx_c, Ky_c, B)  # raw HSIC null values
    # Laplace p-value correction matching dHSIC: (sum >= obs + 1) / (B + 1)
    pval   <- (sum(T_null >= T_obs) + 1L) / (B + 1L)

    stat <- T_obs_n; names(stat) <- "n * HSIC"
    meth <- sprintf("HSIC permutation test (B = %d, %s kernel)", B, kernel_x)
  }

  # ---- Asymptotic / spectral (Zhang; Suzuki 2025, Proposition 16) -----------
  if (method == "eigenvalue") {
    stat_sc <- T_obs * n^3  # n^3 * HSIC for the spectral null scale

    lam_x  <- Re(eigen(Kx_c, symmetric = TRUE)$values)
    lam_y  <- Re(eigen(Ky_c, symmetric = TRUE)$values)
    lam_xy <- as.vector(outer(lam_x, lam_y, `*`))

    T_null <- replicate(B, sum(lam_xy * stats::rchisq(length(lam_xy), df = 1)))
    pval   <- mean(T_null >= stat_sc)

    stat <- stat_sc; names(stat) <- "n^3 * HSIC"
    meth <- sprintf(
      "HSIC eigenvalue test — Zhang/Suzuki Prop. 16 (B = %d MC draws, %s kernel)",
      B, kernel_x)
  }

  # ---- Bootstrap (C++ loop) -------------------------------------------------
  if (method == "bootstrap") {
    T_null <- bootstrap_hsic_cpp(Kx, Ky, B)    # raw HSIC null values
    # Laplace p-value correction matching dHSIC
    pval   <- (sum(T_null >= T_obs) + 1L) / (B + 1L)

    stat <- T_obs_n; names(stat) <- "n * HSIC"
    meth <- sprintf("HSIC bootstrap test (B = %d, %s kernel)", B, kernel_x)
  }

  structure(
    list(statistic  = stat,
         p.value    = pval,
         method     = meth,
         data.name  = dname,
         bandwidths = c(x = bw_x, y = bw_y)),
    class = "htest"
  )
}


# ------------------------------------------------------------------------------
# hsic_resid_test
# ------------------------------------------------------------------------------

#' HSIC Regression Residual Independence Test
#'
#' @description Tests \eqn{H_0: X \perp\!\!\!\perp \hat{e}} where
#'   \eqn{\hat{e}} are OLS residuals from a fitted linear model.  Rejection
#'   indicates misspecification: a missing nonlinear term, omitted variable,
#'   or heteroscedasticity.
#'
#'   The null distribution is estimated via the pairs-bootstrap of the
#'   \pkg{dHSIC} source (Peters et al., 2022).  The test statistic is
#'   \eqn{n \times \widehat{\mathrm{HSIC}}}, consistent with
#'   \code{\link{hsic_test}}.
#'
#' @param model A fitted model (typically \code{\link[stats]{lm}}) with
#'   working \code{model.matrix}, \code{coef}, \code{nobs}, \code{resid}.
#' @param x Numeric predictor matrix for the HSIC computation.  If
#'   \code{NULL} (default), \code{model.matrix(model)} is used with the
#'   intercept column removed.
#' @param kernel_x Kernel for \eqn{X}.  Default \code{"gaussian"}.
#' @param kernel_y Kernel for residuals.  Default \code{kernel_x}.
#' @param bandwidth_x Bandwidth (\eqn{\sigma^2}) for \eqn{X}.  \code{NULL}
#'   uses the median heuristic.
#' @param bandwidth_y Bandwidth for residuals.  \code{NULL} uses the median
#'   heuristic.
#' @param degree Polynomial degree.  Default \code{2}.
#' @param coef0 Polynomial constant.  Default \code{1}.
#' @param B Bootstrap replications.  Default \code{1000}.
#'
#' @return An \code{"htest"} object with \code{statistic} = \eqn{n \times
#'   \widehat{\mathrm{HSIC}}(X, \hat{e})}, \code{p.value} from the
#'   Laplace-corrected bootstrap \eqn{(\#\{T_b \geq T_{\mathrm{obs}}\}+1)/(B+1)},
#'   and \code{bandwidths}.
#'
#' @seealso \code{\link{hsic}}, \code{\link{hsic_test}}
#'
#' @references
#' Peters, J., Pfister, N., & Mooij, J. M. (2022). \emph{dHSIC}. R package
#'   version 2.1.
#'
#' @examples
#' \dontrun{
#' set.seed(99); n <- 100; x <- rnorm(n)
#' hsic_resid_test(lm(2*x + rnorm(n) ~ x), B = 299)       # large p-value
#' hsic_resid_test(lm(x^2 + rnorm(n, sd=0.5) ~ x), B = 299) # small p-value
#' }
#'
#' @keywords internal
hsic_resid_test <- function(model,
                             x           = NULL,
                             kernel_x    = "gaussian",
                             kernel_y    = kernel_x,
                             bandwidth_x = NULL,
                             bandwidth_y = NULL,
                             degree      = 2L,
                             coef0       = 1,
                             B           = 1000L) {

  kernel_x <- match.arg(kernel_x, c("gaussian", "laplace", "linear", "polynomial"))
  kernel_y <- match.arg(kernel_y, c("gaussian", "laplace", "linear", "polynomial"))

  X <- stats::model.matrix(model)
  b <- as.matrix(stats::coef(model))
  n <- stats::nobs(model)
  e <- as.matrix(stats::resid(model))

  if (is.null(x)) {
    int_col <- which(colnames(X) == "(Intercept)")
    x <- if (length(int_col)) X[, -int_col, drop = FALSE] else X
  }
  if (!is.matrix(x)) x <- as.matrix(x)
  if (nrow(x) != n)
    stop("'x' must have the same number of rows as observations in 'model'.",
         call. = FALSE)

  bw_x <- if (kernel_x %in% c("gaussian", "laplace")) .get_bandwidth(x, bandwidth_x) else 1
  bw_y <- if (kernel_y %in% c("gaussian", "laplace")) .get_bandwidth(e, bandwidth_y) else 1

  # Observed raw HSIC(X, residuals)
  T_obs <- hsic_trace_cpp(
    center_kernel_cpp(build_kernel_matrix_cpp(x, kernel_x, bw_x, degree, coef0)),
    center_kernel_cpp(build_kernel_matrix_cpp(e, kernel_y, bw_y, degree, coef0))
  )

  e0 <- e - mean(e)  # centered residuals (dHSIC pairs-bootstrap convention)

  T_null <- replicate(B, {
    idx_e  <- sample(n, replace = TRUE)
    e_B    <- e0[idx_e, , drop = FALSE]
    idx_x  <- sample(n, replace = TRUE)
    x_B    <- x[idx_x, , drop = FALSE]
    X_B    <- X[idx_x, , drop = FALSE]
    Yhat_B <- X_B %*% b + e_B
    bhat_B <- tryCatch(
      solve(crossprod(X_B), crossprod(X_B, Yhat_B)),
      error = function(err) MASS::ginv(crossprod(X_B)) %*% crossprod(X_B, Yhat_B)
    )
    ehat_B <- Yhat_B - X_B %*% bhat_B
    hsic_trace_cpp(
      center_kernel_cpp(build_kernel_matrix_cpp(x_B,    kernel_x, bw_x, degree, coef0)),
      center_kernel_cpp(build_kernel_matrix_cpp(ehat_B, kernel_y, bw_y, degree, coef0))
    )
  })

  # Laplace-corrected p-value and n*HSIC statistic (consistent with hsic_test)
  pval <- (sum(T_null >= T_obs) + 1L) / (B + 1L)
  stat <- T_obs * n; names(stat) <- "n * HSIC"

  structure(
    list(
      statistic  = stat,
      p.value    = pval,
      method     = sprintf(
        "HSIC regression residual test — dHSIC bootstrap (B = %d, %s kernel)",
        B, kernel_x),
      data.name  = deparse(stats::formula(model)),
      bandwidths = c(x = bw_x, y = bw_y)
    ),
    class = "htest"
  )
}
