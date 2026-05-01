# RBF (Gaussian) Kernel Helper
rbf_kernel <- function(x, y, sigma2 = 1) {
  exp(-sum((x - y)^2) / (2 * sigma2))
}

# Kernel Matrix Builder
build_kernel_matrix <- function(X, k_func, sigma2 = 1) {
  n <- if (is.vector(X)) length(X) else nrow(X)
  K <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      x_i <- if (is.vector(X)) X[i] else X[i, ]
      x_j <- if (is.vector(X)) X[j] else X[j, ]
      K[i, j] <- k_func(x_i, x_j, sigma2)
    }
  }
  return(K)
}

# Centering Matrix Helper
center_kernel_matrix <- function(K) {
  n <- nrow(K)
  H <- diag(n) - matrix(1/n, n, n)
  return(H %*% K %*% H)
}

# Main Gamma-Approximated HSIC Test
hsic_test_gamma <- function(x, y, sigma2 = 1, alpha_level = 0.05) {
  n <- if(is.vector(x)) length(x) else nrow(x)

  # 1: center the kernel matrices
  Kx <- build_kernel_matrix(x, rbf_kernel, sigma2)
  Ky <- build_kernel_matrix(y, rbf_kernel, sigma2)

  Kx_tilde <- center_kernel_matrix(Kx)
  Ky_tilde <- center_kernel_matrix(Ky)

  # 2: test statistic (T = n^3 * HSIC)
  # sum(A * B) is equiv to tr(A %*% B) for symmetric matrices but faster
  tr_Kx_Ky <- sum(Kx_tilde * Ky_tilde)
  hsic_val <- tr_Kx_Ky / (n^2)
  test_stat <- tr_Kx_Ky * n

  # 3: theoretical moments from the traces
  # E[T] = tr(Kx_tilde) * tr(Ky_tilde)
  mu_x <- sum(diag(Kx_tilde))
  mu_y <- sum(diag(Ky_tilde))
  mean_T <- mu_x * mu_y

  # Var(T) = 2 * tr(Kx_tilde^2) * tr(Ky_tilde^2)
  var_x <- sum(Kx_tilde^2) # SS elements = trace of squared matrix
  var_y <- sum(Ky_tilde^2)
  var_T <- 2 * var_x * var_y

  # 4: Map moments to Gamma dist params
  gamma_shape <- (mean_T^2) / var_T
  gamma_scale <- var_T / mean_T

  # 5: p-value and crit region using the Gamma dist
  p_value <- pgamma(test_stat, shape = gamma_shape, scale = gamma_scale, lower.tail = FALSE)
  critical_val <- qgamma(1 - alpha_level, shape = gamma_shape, scale = gamma_scale)

  #  output
  structure(list(
    method = "HSIC Test via Gamma Approximation",
    hsic_value = hsic_val,
    test_statistic = test_stat,
    critical_value = critical_val,
    p_value = p_value,
    gamma_shape = gamma_shape,
    gamma_scale = gamma_scale,
    n = n
  ), class = "htest")
}
