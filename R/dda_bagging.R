data("mtcars")

n <- 500
x <- rchisq(n, df = 4) - 4
e <- rchisq(n, df = 3) - 3
y <- 0.5 * x + e
d <- data.frame(x, y)

# 1. dda.indep: Direction Dependence Analysis - Independence
# result_indep <- dda.indep(mpg ~ wt + hp, pred = "wt", data = mtcars)
# print(result_indep)

res_ex_indep <- dda.indep(y ~ x, pred = "x", data = d, parallelize = TRUE, cores = 2,
                                     nlfun = 2, B = 2000, hetero = TRUE, diff = TRUE)
print(res_ex_indep)

init <- Sys.time()
bagged_indep <- dda_bagging(res_ex_indep, iter = 100)
Sys.time() - init #Time difference of 1.757695 hours
summary.dda_bagging_indep(bagged_indep)

# ===== DDA Bagging Summary (INDEP) =====
# Function: dda.indep
# Object Type: dda.indep
# Iterations: 100
# Completed Iterations: 100
# ----
#
#   HSIC and dCor Test Statistics & Harmonic p-values:
#        Target Stat Target p Alternative Stat Alternative p
# HSIC      0.9399        0           7.7198             0
# ---
#
#   Difference Statistics (mean estimates, lower, upper):
#      HSIC   dCor     MI
# HSIC 6.7799 4.8988 8.7575
# dCor 0.2401 0.1681 0.2788
# MI   0.2661 0.1778 0.3741
# ---
#
#   Decision proportions for hsic :
#  undecided      y->x      x->y
#  0.83      0.17      0.00
# ---


# 2. dda.resdist: Direction Dependence Analysis - Residual Distribution

#result_resdist <- dda.resdist(mpg ~ wt + hp, pred = "wt", data = mtcars)
result_resdist <- dda.resdist(y ~ x, pred = "x", data = d,
                      B = 2000, conf.level = 0.90)

print(result_resdist)

init <- Sys.time()
bagged_resdist <- dda_bagging(result_resdist, iter = 100)
Sys.time() - init # Time difference of 3.210665 mins
summary.dda_bagging_resdist(bagged_resdist)

# ===== DDA Bagging Summary =====
# Function: dda.resdist
# Object Type: dda.resdist
# Iterations: 100
# Completed Iterations: 100
# ----
#
#   DIRECTION DEPENDENCE ANALYSIS: Residual Distributions (Bagged)
#
# Skewness and kurtosis tests:
#           y         z-value   Pr(>|z|)  x         z-value   Pr(>|z|)
# Skewness   1.5826  10.7247    0.0000     0.5963   5.0876    0.0000
# Kurtosis   3.2643   6.2666    0.0000     1.7091   4.4472    0.0000
#
# Skewness and kurtosis difference tests and 90% Percentile bootstrap CIs:
#
#             diff      z-value   Pr(>|z|)  lower     upper
# Skewness   -2.1925   -3.6485    0.0062   -3.3565   -0.9905
# Kurtosis   -9.7248   -1.3542    0.3247  -27.2791    4.3254
#
# 90% Percentile bootstrap CIs for joint higher moment differences:
#                                 estimate  lower    upper
# Co-Skewness                     0.5033    0.2428   0.7980
# Hyvarinen-Smith (Co-Skewness)   0.7403    0.4641   1.0112
# Co-Kurtosis                     7.3345    1.5167  15.2955
# Hyvarinen-Smith (Co-Kurtosis)   0.6944    0.1872   1.1881
# Chen-Chan (Co-Kurtosis)         3.1095    0.3550   8.1098
#
# Number of resamples: 2000
# ---
# Note: Target is x -> y
# Alternative is y -> x
# Difference statistics > 0 suggest the model x -> y

# ===== DDA Bagging Summary =====
# Function: dda.resdist
# Object Type: dda.resdist
# Iterations: 100
# Completed Iterations: 100
# ----
#
#   DIRECTION DEPENDENCE ANALYSIS: Residual Distributions (Bagged)
#
# Skewness and kurtosis tests:
#            y         z-value   Pr(>|z|)  x         z-value   Pr(>|z|)
# Skewness   1.5826  10.7247    0.0000     0.5963   5.0876    0.0000
# Kurtosis   3.2643   6.2666    0.0000     1.7091   4.4472    0.0000
#
# Skewness and kurtosis difference tests and 90% Percentile bootstrap CIs:
#
#            diff      z-value   Pr(>|z|)  lower     upper
# Skewness   -2.1925   -3.6485    0.0062   -3.3565   -0.9905
# Kurtosis   -9.7248   -1.3542    0.3247  -27.2791    4.3254
#
# 90% Percentile bootstrap CIs for joint higher moment differences:
#                                 estimate  lower    upper
# Co-Skewness                     0.5033    0.2428   0.7980
# Hyvarinen-Smith (Co-Skewness)   0.7403    0.4641   1.0112
# Co-Kurtosis                     7.3345    1.5167  15.2955
# Hyvarinen-Smith (Co-Kurtosis)   0.6944    0.1872   1.1881
# Chen-Chan (Co-Kurtosis)         3.1095    0.3550   8.1098
#
# Number of resamples: 2000
# ---
#   Note: Target is x -> y
# Alternative is y -> x
# Difference statistics > 0 suggest the model x -> y


###############################################################################
# 3. dda.vardist: Direction Dependence Analysis - Variable Distribution
#result_vardist <- dda.vardist(mpg ~ wt + hp, pred = "wt", data = mtcars)
result_vardist <- dda.vardist(y ~ x, pred = "x", data = d, B = 2000)

print(result_vardist)

init <- Sys.time()
bagged_vardist <- dda_bagging(result_vardist, iter = 100)
Sys.time() - init #Time difference of 2.703492 mins
summary.dda_bagging_vardist(bagged_vardist)
#
# ===== DDA Bagging Summary =====
#   Function: dda.vardist
# Object Type: dda.vardist
# Iterations: 100
# Completed Iterations: 100
# ----
#
#   DIRECTION DEPENDENCE ANALYSIS: Variable Distributions (Bagged)
#
# Skewness and kurtosis tests:
#          y        z-value   Pr(>|z|)  x        z-value   Pr(>|z|)
# Skewness   1.304   9.460     0.000      1.591  10.734     0.000
# Kurtosis   2.382   5.554     0.000      3.737   6.445     0.000
#
# 95% Percentile bootstrap CIs for higher moment differences:
#   diff     lower    upper
# Skewness   0.8970  -0.5912   2.5777
# Kurtosis  12.9933  -9.5477  49.5642
#
# 95% Percentile bootstrap CIs for differences in higher-order correlations:
#   estimate  lower    upper
# Cor^2[2,1] - Cor^2[1,2]   0.5236    0.2010   0.9365
# Cor^2[3,1] - Cor^2[1,3]   9.8078   -0.9159  25.8260
#
# 95% Percentile bootstrap CIs for Likelihood Ratio approximations:
#   estimate  lower    upper
# Hyvarinen-Smith (co-skewness)   0.1882    0.0932   0.2785
# Hyvarinen-Smith (tanh)          0.0147    0.0063   0.0223
# Chen-Chan (co-kurtosis)         5.3845   -0.5292  16.9036
#
# Number of resamples: 2000
# ---
#   Note: (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests the model x -> y
#############################################################################

#' Bootstrap Aggregated DDA Analysis (harmonic mean for p-values, mean for other stats)
#'
#' @param dda_result Output from any DDA function (dda.indep, dda.vardist, dda.resdist, etc.)
#' @param iter Number of bootstrap iterations (default: 100)
#' @param progress Whether to show progress bar (default: TRUE)
#' @param save_file Optional file path to save results
#' @param alpha Significance level for decisions (default: 0.05)
#' @return A list containing bootstrap and aggregated results

dda_bagging <- function(
    dda_result,
    iter = 100,
    progress = TRUE,
    save_file = NULL,
    alpha = 0.05
) {
  get_numeric <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (is.numeric(x)) return(as.numeric(x[1]))
    if (is.list(x)) return(get_numeric(x[[1]]))
    return(NA_real_)
  }
  harmonic_p <- function(pvec) {
    pvec <- as.numeric(pvec)
    pvec <- pvec[!is.na(pvec) & pvec > 0 & pvec < 1]
    if (length(pvec) == 0) return(NA_real_)
    if (!requireNamespace("harmonicmeanp", quietly = TRUE)) {
      warning("harmonicmeanp package not found, using mean for p-values")
      return(mean(pvec))
    }
    harmonicmeanp::p.hmp(pvec, L = length(pvec))
  }
  print_decisions <- function(x) {
    levs <- c("undecided", "y->x", "x->y")
    tab <- table(factor(x, levels = levs))
    prop <- tab / sum(tab)
    names(prop) <- levs
    prop
  }

  # Extract call info and data
  if (inherits(dda_result, "dda.indep") ||
      inherits(dda_result, "dda.vardist") ||
      inherits(dda_result, "dda.resdist") ||
      inherits(dda_result, "cdda.indep") ||
      inherits(dda_result, "cdda.vardist")) {
    call_info <- dda_result$call_info
    function_name <- call_info$function_name
    all_args <- call_info$all_args
    data_name <- call_info$data_name
    original_data <- call_info$original_data
    if (is.null(original_data)) {
      tryCatch({
        original_data <- get(data_name, envir = parent.frame())
      }, error = function(e) {
        stop(paste("Could not retrieve original data:", data_name))
      })
    }
  } else if (length(dda_result) >= 5 && !is.null(dda_result[[5]])) {
    call_info <- dda_result[[5]]
    function_name <- call_info$function_name
    all_args <- call_info$all_args
    original_env <- call_info$environment
    data_name <- call_info$data_name
    original_data <- get(data_name, envir = original_env)
  } else {
    stop("Unrecognized DDA result structure. Please ensure you're using a supported DDA object.")
  }

  nobs <- nrow(original_data)
  dda_function <- get(function_name)

  bagged_results <- vector("list", iter)
  if (progress) {
    pb <- txtProgressBar(min = 0, max = iter, style = 3, width = 50, char = "=")
    cat(paste("\n", "Running", iter, "bootstrap iterations of", function_name, "\n"))
  }

  for(i in seq_len(iter)) {
    boot_indices <- sample(seq_len(nobs), nobs, replace = TRUE)
    boot_data <- original_data[boot_indices, ]
    boot_args <- all_args
    boot_args$data <- boot_data
    res <- tryCatch({
      do.call(dda_function, boot_args)
    }, error = function(e) NA)
    bagged_results[[i]] <- res
    if (progress) setTxtProgressBar(pb, i)
  }
  if (progress) { close(pb); cat("\nBootstrap iterations completed.\n") }

  valid_results <- bagged_results[!sapply(bagged_results, function(x) all(is.na(x)) || is.null(x))]
  object_type <- class(dda_result)[1]
  agg <- list()
  decisions <- list()
  n_valid <- length(valid_results)
  class_label <- "dda_bagging"

  ## --- DDA.INDEP block (now includes decisions) ---
  if (n_valid > 0 && object_type == "dda.indep") {
    hsic_yx <- sapply(valid_results, function(x) get_numeric(x$hsic.yx$statistic))
    hsic_xy <- sapply(valid_results, function(x) get_numeric(x$hsic.xy$statistic))
    hsic_yx_pval <- sapply(valid_results, function(x) get_numeric(x$hsic.yx$p.value))
    hsic_xy_pval <- sapply(valid_results, function(x) get_numeric(x$hsic.xy$p.value))
    dcor_yx <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_yx$statistic))
    dcor_xy <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_xy$statistic))
    dcor_yx_pval <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_yx$p.value))
    dcor_xy_pval <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_xy$p.value))

    agg$hsic_yx_stat <- mean(hsic_yx, na.rm = TRUE)
    agg$hsic_xy_stat <- mean(hsic_xy, na.rm = TRUE)
    agg$dcor_yx_stat <- mean(dcor_yx, na.rm = TRUE)
    agg$dcor_xy_stat <- mean(dcor_xy, na.rm = TRUE)

    agg$hsic_yx_pval <- harmonic_p(hsic_yx_pval)
    agg$hsic_xy_pval <- harmonic_p(hsic_xy_pval)
    agg$dcor_yx_pval <- harmonic_p(dcor_yx_pval)
    agg$dcor_xy_pval <- harmonic_p(dcor_xy_pval)

    # Decision proportions (robust to NAs)
    hsic_decision <- ifelse(is.na(hsic_yx_pval) | is.na(hsic_xy_pval), "undecided",
                            ifelse(hsic_yx_pval < alpha & hsic_xy_pval >= alpha, "x->y",
                                   ifelse(hsic_yx_pval >= alpha & hsic_xy_pval < alpha, "y->x",
                                          "undecided")))
    decisions$hsic <- print_decisions(hsic_decision)

    # Difference statistics if present
    if (!is.null(valid_results[[1]]$out.diff)) {
      first_diff <- valid_results[[1]]$out.diff
      diff_mat <- lapply(valid_results, function(x) {
        if (!is.null(x$out.diff)) {
          as.matrix(x$out.diff)
        } else {
          matrix(NA_real_, nrow = nrow(first_diff), ncol = ncol(first_diff),
                 dimnames = dimnames(first_diff))
        }
      })
      diff_array <- simplify2array(diff_mat)
      agg$diff_matrix <- apply(diff_array, c(1,2), function(xx) mean(xx, na.rm=TRUE))
    }

    class_label <- "dda_bagging_indep"
  }

  ## --- DDA.RESDIST block (robust to all/some missing vectors) ---
  if (n_valid > 0 && object_type == "dda.resdist") {
    agg$var.names <- if (!is.null(valid_results[[1]]$var.names)) valid_results[[1]]$var.names else c("target", "alternative")
    agg$probtrans <- if (!is.null(valid_results[[1]]$probtrans)) valid_results[[1]]$probtrans else FALSE

    agg$agostino.target.statistic <- mean(sapply(valid_results, function(x) get_numeric(x$agostino$target$statistic[1])), na.rm=TRUE)
    agg$agostino.target.z         <- mean(sapply(valid_results, function(x) get_numeric(x$agostino$target$statistic[2])), na.rm=TRUE)
    agg$agostino.target.p.value   <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$agostino$target$p.value)))
    agg$agostino.alternative.statistic <- mean(sapply(valid_results, function(x) get_numeric(x$agostino$alternative$statistic[1])), na.rm=TRUE)
    agg$agostino.alternative.z         <- mean(sapply(valid_results, function(x) get_numeric(x$agostino$alternative$statistic[2])), na.rm=TRUE)
    agg$agostino.alternative.p.value   <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$agostino$alternative$p.value)))

    agg$anscombe.target.statistic <- mean(sapply(valid_results, function(x) get_numeric(x$anscombe$target$statistic[1])), na.rm=TRUE)
    agg$anscombe.target.z         <- mean(sapply(valid_results, function(x) get_numeric(x$anscombe$target$statistic[2])), na.rm=TRUE)
    agg$anscombe.target.p.value   <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$anscombe$target$p.value)))
    agg$anscombe.alternative.statistic <- mean(sapply(valid_results, function(x) get_numeric(x$anscombe$alternative$statistic[1])), na.rm=TRUE)
    agg$anscombe.alternative.z         <- mean(sapply(valid_results, function(x) get_numeric(x$anscombe$alternative$statistic[2])), na.rm=TRUE)
    agg$anscombe.alternative.p.value   <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$anscombe$alternative$p.value)))

    mean_diff <- function(key, len) {
      out <- lapply(valid_results, function(x) {
        v <- x[[key]]
        if (is.null(v)) return(rep(NA, len))
        vnum <- suppressWarnings(as.numeric(v))
        if (length(vnum) != len) return(rep(NA, len))
        vnum
      })
      mat <- do.call(rbind, out)
      if (is.null(mat) || length(mat) == 0) return(rep(NA, len))
      colMeans(mat, na.rm = TRUE)
    }
    agg$skewdiff <- mean_diff("skewdiff", 5)
    agg$kurtdiff <- mean_diff("kurtdiff", 5)
    agg$cor12diff <- mean_diff("cor12diff", 3)
    agg$cor13diff <- mean_diff("cor13diff", 3)
    agg$RHS3      <- mean_diff("RHS3", 3)
    agg$RHS4      <- mean_diff("RHS4", 3)
    agg$RCC       <- mean_diff("RCC", 3)

    class_label <- "dda_bagging_resdist"
  }

  ## --- DDA.VARDIST block (robust; mean stats, harmonic-mean p-values) ---
  if (n_valid > 0 && object_type == "dda.vardist") {
    # Keep variable names for summary headers
    agg$var.names <- if (!is.null(valid_results[[1]]$var.names)) valid_results[[1]]$var.names else c("Outcome", "Predictor")

    # Helper: robust mean of a scalar nested path (skips NULLs/NA)
    get_mean <- function(f) {
      vals <- sapply(valid_results, function(x) {
        out <- tryCatch(f(x), error = function(e) NA_real_)
        get_numeric(out)
      })
      mean(vals, na.rm = TRUE)
    }
    # Helper: harmonic mean p for a scalar nested path (skips NULLs/NA)
    get_hmp <- function(f) {
      pvals <- sapply(valid_results, function(x) {
        out <- tryCatch(f(x), error = function(e) NA_real_)
        get_numeric(out)
      })
      harmonic_p(pvals)
    }
    # Helper: robust vector aggregator (len columns) using column means
    mean_vec <- function(key, len) {
      rows <- lapply(valid_results, function(x) {
        v <- tryCatch(x[[key]], error = function(e) NULL)
        if (is.null(v)) return(rep(NA_real_, len))
        vnum <- suppressWarnings(as.numeric(v))
        if (length(vnum) < len) {
          # pad with NAs if shorter
          vnum <- c(vnum, rep(NA_real_, len - length(vnum)))
        } else if (length(vnum) > len) {
          vnum <- vnum[seq_len(len)]
        }
        vnum
      })
      mat <- do.call(rbind, rows)
      if (is.null(mat) || length(mat) == 0) return(rep(NA_real_, len))
      as.numeric(colMeans(mat, na.rm = TRUE))
    }

    # Agostino (skewness) and Anscombe (kurtosis) tests
    # Outcome is varnames[1]; Predictor is varnames[2]
    agg$agostino.outcome.statistic.skew   <- get_mean(function(x) x$agostino$outcome$statistic[1])
    agg$agostino.outcome.statistic.z      <- get_mean(function(x) x$agostino$outcome$statistic[2])
    agg$agostino.outcome.p.value          <- get_hmp(function(x) x$agostino$outcome$p.value)

    agg$agostino.predictor.statistic.skew <- get_mean(function(x) x$agostino$predictor$statistic[1])
    agg$agostino.predictor.statistic.z    <- get_mean(function(x) x$agostino$predictor$statistic[2])
    agg$agostino.predictor.p.value        <- get_hmp(function(x) x$agostino$predictor$p.value)

    agg$anscombe.outcome.statistic.kurt   <- get_mean(function(x) x$anscombe$outcome$statistic[1])
    agg$anscombe.outcome.statistic.z      <- get_mean(function(x) x$anscombe$outcome$statistic[2])
    agg$anscombe.outcome.p.value          <- get_hmp(function(x) x$anscombe$outcome$p.value)

    agg$anscombe.predictor.statistic.kurt <- get_mean(function(x) x$anscombe$predictor$statistic[1])
    agg$anscombe.predictor.statistic.z    <- get_mean(function(x) x$anscombe$predictor$statistic[2])
    agg$anscombe.predictor.p.value        <- get_hmp(function(x) x$anscombe$predictor$p.value)

    # Higher-moment differences and CIs (mean across iterations)
    # Each expected as c(value, lower, upper)
    agg$skewdiff  <- mean_vec("skewdiff",  3)  # c(diff, lower, upper)
    agg$kurtdiff  <- mean_vec("kurtdiff",  3)  # c(diff, lower, upper)
    agg$cor12diff <- mean_vec("cor12diff", 3)  # c(estimate, lower, upper)
    agg$cor13diff <- mean_vec("cor13diff", 3)  # c(estimate, lower, upper)
    agg$RHS       <- mean_vec("RHS",       3)  # c(estimate, lower, upper)
    agg$Rtanh     <- mean_vec("Rtanh",     3)  # c(estimate, lower, upper)
    agg$RCC       <- mean_vec("RCC",       3)  # c(estimate, lower, upper)

    # Ensure the bagged object can dispatch to a vardist summary
    class_label <- "dda_bagging_vardist"
  }

  output <- list(
    bagged_results = bagged_results,
    aggregated_stats = agg,
    decision_percentages = decisions,
    n_valid_iterations = n_valid,
    original_result = dda_result,
    parameters = list(
      function_name = function_name,
      object_type = object_type,
      iter = iter,
      successful_iterations = n_valid,
      failed_iterations = iter - n_valid
    )
  )
  if (!is.null(save_file)) {
    save(output, file = save_file)
    cat(paste("Results saved to:", save_file, "\n"))
  }
  class(output) <- c(class_label, "dda_bagging", class(output))
  return(output)
}
