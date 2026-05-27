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

    bool is_gauss  = (kernel == "gaussian");
    bool is_lap    = (kernel == "laplace");
    bool is_linear = (kernel == "linear");
    bool is_poly   = (kernel == "polynomial");

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
