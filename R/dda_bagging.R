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
                                     nlfun = 2, B = 50, hetero = TRUE, diff = TRUE)
print(res_ex_indep)

bagged_indep <- dda_bagging(res_ex_indep, iter = 10)
summary.dda_bagging_indep(bagged_indep)

# ===== DDA Bagging Summary =====
# Function: dda.indep
# Object Type: dda.indep
# Iterations: 10
# Valid Iterations: 10
# ----
#
# HSIC and dCor Test Statistics & Harmonic p-values:
# Target Stat Target p Alternative Stat Alternative p
# HSIC       0.697    5e-04           6.4876             0
# ---
#
# Difference Statistics (mean estimates, lower, upper):
#          HSIC   dCor     MI
# estimate 5.7906 0.2202 0.1964
# lower    3.8147 0.1431 0.1279
# upper    7.6252 0.2491 0.2710
# ---
#
# Decision proportions for hsic :
# undecided      y->x      x->y
#       0.6       0.4       0.0
# ---


# 2. dda.resdist: Direction Dependence Analysis - Residual Distribution

# result_resdist <- dda.resdist(mpg ~ wt + hp, pred = "wt", data = mtcars)
result_resdist <- dda.resdist(y ~ x, pred = "x", data = d,
                      B = 50, conf.level = 0.90)


print(result_resdist)

bagged_resdist <- dda_bagging(result_resdist, iter = 100)
summary.dda_bagging_resdist(bagged_resdist)


# 3. dda.vardist: Direction Dependence Analysis - Variable Distribution
result_vardist <- dda.vardist(mpg ~ wt + hp, pred = "wt", data = mtcars)
print(result_vardist)

bagged_vardist <- dda_bagging(result_vardist, iter = 10)
summary.dda_bagging_vardist(bagged_vardist)


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

  bagged_results <- list()
  if (progress) {
    pb <- txtProgressBar(min = 0, max = iter, style = 3, width = 50, char = "=")
    cat(paste("\n", "Running", iter, "bootstrap iterations of", function_name, "\n"))
  }

  for(i in 1:iter) {
    boot_indices <- sample(1:nobs, nobs, replace = TRUE)
    boot_data <- original_data[boot_indices, ]
    boot_args <- all_args
    boot_args$data <- boot_data
    tryCatch({
      bagged_results[[i]] <- do.call(dda_function, boot_args)
    }, error = function(e) {
      warning(paste("DDA function failed at iteration", i, ":", e$message))
      bagged_results[[i]] <- NA
    })
    if (progress) setTxtProgressBar(pb, i)
  }
  if (progress) { close(pb); cat("\nBootstrap iterations completed.\n") }

  valid_results <- bagged_results[!sapply(bagged_results, function(x) all(is.na(x)) || is.null(x))]
  object_type <- class(dda_result)[1]
  agg <- list()
  decisions <- list()
  n_valid <- length(valid_results)


  ## --- DDA.INDEP block ---
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

    hsic_decision <- ifelse(hsic_yx_pval < alpha & hsic_xy_pval >= alpha, "x->y",
                            ifelse(hsic_yx_pval >= alpha & hsic_xy_pval < alpha, "y->x", "undecided"))
    dcor_decision <- ifelse(dcor_yx_pval < alpha & dcor_xy_pval >= alpha, "x->y",
                            ifelse(dcor_yx_pval >= alpha & dcor_xy_pval < alpha, "y->x", "undecided"))
    decisions$hsic <- print_decisions(hsic_decision)
    decisions$dcor <- print_decisions(dcor_decision)

    # Difference statistics if exist
    if (!is.null(valid_results[[1]]$out.diff)) {
      diff_estimates <- sapply(valid_results, function(x) {
        if (!is.null(x$out.diff)) as.numeric(x$out.diff[, "estimate"]) else rep(NA_real_, 3)
      })
      diff_lower <- sapply(valid_results, function(x) {
        if (!is.null(x$out.diff)) as.numeric(x$out.diff[, "lower"]) else rep(NA_real_, 3)
      })
      diff_upper <- sapply(valid_results, function(x) {
        if (!is.null(x$out.diff)) as.numeric(x$out.diff[, "upper"]) else rep(NA_real_, 3)
      })
      col_names <- rownames(valid_results[[1]]$out.diff)
      agg$diff_matrix <- matrix(NA_real_, ncol = 3, nrow = 3)
      agg$diff_matrix[1, ] <- apply(diff_estimates, 1, mean, na.rm = TRUE)
      agg$diff_matrix[2, ] <- apply(diff_lower,    1, mean, na.rm = TRUE)
      agg$diff_matrix[3, ] <- apply(diff_upper,    1, mean, na.rm = TRUE)
      colnames(agg$diff_matrix) <- col_names
      rownames(agg$diff_matrix) <- c("estimate", "lower", "upper")
    }
  }

  ## --- DDA.RESDIST block ---
  if (n_valid > 0 && object_type == "dda.resdist") {
    # Variable names and probtrans
    agg$var.names <- if (!is.null(valid_results[[1]]$var.names)) valid_results[[1]]$var.names else c("target", "alternative")
    agg$probtrans <- if (!is.null(valid_results[[1]]$probtrans)) valid_results[[1]]$probtrans else FALSE

    # Mean aggregate test statistics, harmonic mean aggregate p-values
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

    # Difference tests (classic print style)
    agg$skewdiff <- apply(sapply(valid_results, function(x) x$skewdiff), 1, mean, na.rm=TRUE)
    agg$kurtdiff <- apply(sapply(valid_results, function(x) x$kurtdiff), 1, mean, na.rm=TRUE)
    agg$cor12diff <- apply(sapply(valid_results, function(x) x$cor12diff), 1, mean, na.rm=TRUE)
    agg$cor13diff <- apply(sapply(valid_results, function(x) x$cor13diff), 1, mean, na.rm=TRUE)
    agg$RHS3      <- apply(sapply(valid_results, function(x) x$RHS3), 1, mean, na.rm=TRUE)
    agg$RHS4      <- apply(sapply(valid_results, function(x) x$RHS4), 1, mean, na.rm=TRUE)
    agg$RCC       <- apply(sapply(valid_results, function(x) x$RCC), 1, mean, na.rm=TRUE)
  }

  #############################################################################

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
  return(output)
}

