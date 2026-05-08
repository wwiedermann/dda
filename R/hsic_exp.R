library(Rcpp)

# ==============================================================================
# 1. C++ BACKEND (Rcpp)
# ==============================================================================
cpp_code <- "
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
double median_bandwidth_cpp(NumericMatrix x) {
    int n = x.nrow();
    int p = x.ncol();
    std::vector<double> dists;
    dists.reserve((n * (n - 1)) / 2);
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d2 = 0;
            for (int k = 0; k < p; k++) {
                double diff = x(i, k) - x(j, k);
                d2 += diff * diff;
            }
            if (d2 > 0) dists.push_back(d2);
        }
    }
    
    if (dists.empty()) return 1.0;
    
    // Efficient median finding
    size_t mid = dists.size() / 2;
    std::nth_element(dists.begin(), dists.begin() + mid, dists.end());
    double med = dists[mid];
    
    return (med == 0) ? 1.0 : med;
}

// [[Rcpp::export]]
NumericMatrix build_kernel_matrix_cpp(NumericMatrix x, String kernel, double bandwidth, int degree, double coef0) {
    int n = x.nrow();
    int p = x.ncol();
    NumericMatrix K(n, n);
    
    bool is_linear = (kernel == \"linear\");
    bool is_poly = (kernel == \"polynomial\");
    bool is_gauss = (kernel == \"gaussian\");
    
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double val = 0;
            if (is_linear || is_poly) {
                for (int k = 0; k < p; k++) val += x(i, k) * x(j, k);
                if (is_poly) val = std::pow(val + coef0, degree);
            } else {
                double d2 = 0;
                for (int k = 0; k < p; k++) {
                    double diff = x(i, k) - x(j, k);
                    d2 += diff * diff;
                }
                if (is_gauss) {
                    val = std::exp(-d2 / bandwidth);
                } else { // laplace
                    val = std::exp(-std::sqrt(d2) / bandwidth);
                }
            }
            K(i, j) = val;
            K(j, i) = val;
        }
    }
    return K;
}

// [[Rcpp::export]]
NumericMatrix center_kernel_cpp(NumericMatrix K) {
    int n = K.nrow();
    NumericMatrix Kc(n, n);
    std::vector<double> row_sums(n, 0.0);
    double total_sum = 0.0;
    
    // O(N^2) centering instead of O(N^3) matrix multiplication
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            row_sums[i] += K(i, j);
            total_sum += K(i, j);
        }
    }
    
    double mean_total = total_sum / (n * n);
    for (int i = 0; i < n; i++) {
        double mean_i = row_sums[i] / n;
        for (int j = 0; j < n; j++) {
            Kc(i, j) = K(i, j) - mean_i - (row_sums[j] / n) + mean_total;
        }
    }
    return Kc;
}

// [[Rcpp::export]]
double hsic_trace_cpp(NumericMatrix Kx_c, NumericMatrix Ky_c) {
    int n = Kx_c.nrow();
    double sum = 0;
    // Trace(A * B) for symmetric matrices equals sum of element-wise products
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum += Kx_c(i, j) * Ky_c(i, j);
        }
    }
    return sum / (n * n);
}

// [[Rcpp::export]]
NumericVector permute_hsic_cpp(NumericMatrix Kx_c, NumericMatrix Ky_c, int B) {
    int n = Kx_c.nrow();
    NumericVector null_dist(B);
    IntegerVector idx = seq(0, n - 1);
    
    // Since H commutes with permutations, we permute the already-centered matrix
    for (int b = 0; b < B; b++) {
        IntegerVector p = sample(idx, n, false);
        double val = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                val += Kx_c(p[i], p[j]) * Ky_c(i, j);
            }
        }
        null_dist[b] = (val / (n * n)) * n; // Scaled by n
    }
    return null_dist;
}

// [[Rcpp::export]]
NumericVector bootstrap_hsic_cpp(NumericMatrix Kx, NumericMatrix Ky, int B) {
    int n = Kx.nrow();
    NumericVector null_dist(B);
    IntegerVector seq_idx = seq(0, n - 1);
    
    for (int b = 0; b < B; b++) {
        IntegerVector ix = sample(seq_idx, n, true);
        IntegerVector iy = sample(seq_idx, n, true);
        
        NumericMatrix Kx_b(n, n);
        NumericMatrix Ky_b(n, n);
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Kx_b(i, j) = Kx(ix[i], ix[j]);
                Ky_b(i, j) = Ky(iy[i], iy[j]);
            }
        }
        
        NumericMatrix Kx_bc = center_kernel_cpp(Kx_b);
        NumericMatrix Ky_bc = center_kernel_cpp(Ky_b);
        
        double val = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                val += Kx_bc(i, j) * Ky_bc(i, j);
            }
        }
        null_dist[b] = (val / (n * n)) * n; // Scaled by n
    }
    return null_dist;
}
"
Rcpp::sourceCpp(code = cpp_code)


# ==============================================================================
# 2. R WRAPPERS (_exp functions)
# ==============================================================================

median_bandwidth_exp <- function(x) {
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  median_bandwidth_cpp(x)
}

build_kernel_matrix_exp <- function(x,
                                    kernel    = c("gaussian", "laplace", "linear", "polynomial"),
                                    bandwidth = NULL,
                                    degree    = 2,
                                    coef0     = 1) {
  kernel <- match.arg(kernel)
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  
  if (is.null(bandwidth) && kernel %in% c("gaussian", "laplace")) {
    bandwidth <- median_bandwidth_exp(x)
  }
  if (bandwidth <= 0 && kernel %in% c("gaussian", "laplace")) {
    stop("'bandwidth' must be strictly positive.")
  }
  
  build_kernel_matrix_cpp(x, kernel, bandwidth, degree, coef0)
}

center_kernel_exp <- function(K) {
  center_kernel_cpp(K)
}

hsic_exp <- function(x, y,
                     kernel_x    = "gaussian",
                     kernel_y    = kernel_x,
                     bandwidth_x = NULL,
                     bandwidth_y = NULL,
                     degree      = 2,
                     coef0       = 1) {
  
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  if (is.vector(y)) y <- matrix(y, ncol = 1)
  n <- nrow(x)
  if (nrow(y) != n) stop("'x' and 'y' must have the same number of observations.")
  
  if (is.null(bandwidth_x)) bandwidth_x <- median_bandwidth_exp(x)
  if (is.null(bandwidth_y)) bandwidth_y <- median_bandwidth_exp(y)
  
  Kx <- build_kernel_matrix_exp(x, kernel = kernel_x, bandwidth = bandwidth_x, degree = degree, coef0 = coef0)
  Ky <- build_kernel_matrix_exp(y, kernel = kernel_y, bandwidth = bandwidth_y, degree = degree, coef0 = coef0)
  
  Kx_c <- center_kernel_exp(Kx)
  Ky_c <- center_kernel_exp(Ky)
  
  hsic_trace_cpp(Kx_c, Ky_c)
}

hsic_test_exp <- function(x, y,
                          method      = c("permutation", "asymptotic", "bootstrap"),
                          kernel_x    = "gaussian",
                          kernel_y    = kernel_x,
                          bandwidth_x = NULL,
                          bandwidth_y = NULL,
                          degree      = 2,
                          coef0       = 1,
                          B           = 1000) {
  
  method <- match.arg(method)
  dname  <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  if (is.vector(y)) y <- matrix(y, ncol = 1)
  n <- nrow(x)
  if (nrow(y) != n) stop("'x' and 'y' must have the same number of observations.")
  
  if (is.null(bandwidth_x)) bandwidth_x <- median_bandwidth_exp(x)
  if (is.null(bandwidth_y)) bandwidth_y <- median_bandwidth_exp(y)
  
  Kx <- build_kernel_matrix_exp(x, kernel = kernel_x, bandwidth = bandwidth_x, degree = degree, coef0 = coef0)
  Ky <- build_kernel_matrix_exp(y, kernel = kernel_y, bandwidth = bandwidth_y, degree = degree, coef0 = coef0)
  
  Kx_c <- center_kernel_exp(Kx)
  Ky_c <- center_kernel_exp(Ky)
  
  T_obs_raw <- hsic_trace_cpp(Kx_c, Ky_c)
  T_obs_scaled <- T_obs_raw * n 
  
  if (method == "permutation") {
    T_null <- permute_hsic_cpp(Kx_c, Ky_c, B)
    pval <- mean(T_null >= T_obs_scaled)
    stat <- T_obs_scaled
    names(stat) <- "HSIC (scaled by n)"
    meth <- sprintf("HSIC permutation test (B = %d, %s kernel, C++ Optimized)", B, kernel_x)
    
  } else if (method == "asymptotic") {
    stat_asymp <- T_obs_raw * n^3
    lambda_x  <- Re(eigen(Kx_c, symmetric = TRUE)$values)
    lambda_y  <- Re(eigen(Ky_c, symmetric = TRUE)$values)
    lambda_xy <- as.vector(outer(lambda_x, lambda_y, `*`))
    
    r <- length(lambda_xy)
    T_null <- replicate(B, sum(lambda_xy * rchisq(r, df = 1)))
    
    pval <- mean(T_null >= stat_asymp)
    stat <- stat_asymp
    names(stat) <- "n^3 * HSIC"
    meth <- sprintf("HSIC asymptotic test — Zhang/Suzuki (B = %d MC draws, %s kernel)", B, kernel_x)
    
  } else if (method == "bootstrap") {
    T_null <- bootstrap_hsic_cpp(Kx, Ky, B)
    pval <- mean(T_null >= T_obs_scaled)
    stat <- T_obs_scaled
    names(stat) <- "HSIC (scaled by n)"
    meth <- sprintf("HSIC bootstrap test (B = %d, %s kernel, C++ Optimized)", B, kernel_x)
  }
  
  structure(
    list(statistic = stat, p.value = pval, method = meth, data.name = dname),
    class = "htest"
  )
}

hsic_resid_test_exp <- function(model,
                                x           = NULL,
                                kernel_x    = "gaussian",
                                kernel_y    = kernel_x,
                                bandwidth_x = NULL,
                                bandwidth_y = NULL,
                                degree      = 2,
                                coef0       = 1,
                                B           = 1000) {
  
  X <- model.matrix(model)
  b <- as.matrix(coef(model))
  n <- nobs(model)
  e <- as.matrix(resid(model))
  
  if (is.null(x)) {
    int_col <- which(colnames(X) == "(Intercept)")
    x <- if (length(int_col)) X[, -int_col, drop = FALSE] else X
  }
  if (!is.matrix(x)) x <- as.matrix(x)
  if (nrow(x) != n) stop("'x' must have the same number of rows as observations in 'model'.")
  
  if (is.null(bandwidth_x)) bandwidth_x <- median_bandwidth_exp(x)
  if (is.null(bandwidth_y)) bandwidth_y <- median_bandwidth_exp(e)
  
  T_obs_raw <- hsic_exp(x, e, kernel_x = kernel_x, kernel_y = kernel_y, 
                        bandwidth_x = bandwidth_x, bandwidth_y = bandwidth_y, 
                        degree = degree, coef0 = coef0)
  
  T_obs_scaled <- T_obs_raw * n
  e0 <- e - mean(e)
  
  # Keeping OLS in R, but calling the C++ HSIC evaluator inside the loop
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
    
    hsic_exp(x_B, ehat_B, kernel_x = kernel_x, kernel_y = kernel_y,
             bandwidth_x = bandwidth_x, bandwidth_y = bandwidth_y,
             degree = degree, coef0 = coef0) * n
  })
  
  pval <- mean(T_null >= T_obs_scaled)
  stat <- T_obs_scaled
  names(stat) <- "HSIC (scaled by n)"
  
  structure(
    list(
      statistic = stat, p.value = pval, 
      method = sprintf("HSIC regression residual test (B = %d, %s kernel, C++ Evaluator)", B, kernel_x),
      data.name = deparse(formula(model))
    ),
    class = "htest"
  )
}

hsic_plot_exp <- function(x, y,
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
  
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  if (is.vector(y)) y <- matrix(y, ncol = 1)
  n <- nrow(x)
  
  if (is.null(bandwidth_x)) bandwidth_x <- median_bandwidth_exp(x)
  if (is.null(bandwidth_y)) bandwidth_y <- median_bandwidth_exp(y)
  
  Kx <- build_kernel_matrix_exp(x, kernel = kernel_x, bandwidth = bandwidth_x, degree = degree, coef0 = coef0)
  Ky <- build_kernel_matrix_exp(y, kernel = kernel_y, bandwidth = bandwidth_y, degree = degree, coef0 = coef0)
  
  Kx_c <- center_kernel_exp(Kx)
  Ky_c <- center_kernel_exp(Ky)
  T_obs <- hsic_trace_cpp(Kx_c, Ky_c) * n
  
  if (method == "permutation") {
    T_null <- permute_hsic_cpp(Kx_c, Ky_c, B)
  } else {
    T_null <- bootstrap_hsic_cpp(Kx, Ky, B)
  }
  
  crit <- quantile(T_null, 1 - alpha)
  xlim <- range(c(T_null, T_obs, crit))
  
  plot(density(T_null), xlim = xlim,
       main = "HSIC null distribution (C++ Optimized)",
       xlab = "HSIC (scaled by n)", ylab = "Density", ...)
  abline(v = crit,  col = "red",  lty = 2, lwd = 2)
  abline(v = T_obs, col = "blue", lty = 1, lwd = 2)
  legend("topright",
         legend = c(sprintf("Critical value (alpha = %.2f)", alpha), "Observed HSIC"),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2, bty = "n")
  
  invisible(list(T_obs = T_obs, T_null = T_null, critical = crit))
}