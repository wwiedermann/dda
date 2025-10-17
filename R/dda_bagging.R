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

bagged_resdist <- dda_bagging(result_resdist, iter = 10)
summary.dda_bagging_resdist(bagged_resdist)


# 3. dda.vardist: Direction Dependence Analysis - Variable Distribution
result_vardist <- dda.vardist(mpg ~ wt + hp, pred = "wt", data = mtcars)
print(result_vardist)

bagged_vardist <- dda_bagging(result_vardist, iter = 10)
summary.dda_bagging_vardist(bagged_vardist)

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

  ## --- DDA.INDEP block ---
  if (n_valid > 0 && object_type == "dda.indep") {
    agg$hsic_yx_stat <- mean(sapply(valid_results, function(x) get_numeric(x$hsic.yx$statistic)), na.rm = TRUE)
    agg$hsic_xy_stat <- mean(sapply(valid_results, function(x) get_numeric(x$hsic.xy$statistic)), na.rm = TRUE)
    agg$dcor_yx_stat <- mean(sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_yx$statistic)), na.rm = TRUE)
    agg$dcor_xy_stat <- mean(sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_xy$statistic)), na.rm = TRUE)

    agg$hsic_yx_pval <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$hsic.yx$p.value)))
    agg$hsic_xy_pval <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$hsic.xy$p.value)))
    agg$dcor_yx_pval <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_yx$p.value)))
    agg$dcor_xy_pval <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_xy$p.value)))

    # Difference statistics if present
    if (!is.null(valid_results[[1]]$out.diff)) {
      diff_mat <- lapply(valid_results, function(x) {
        if (!is.null(x$out.diff)) {
          as.matrix(x$out.diff)
        } else {
          matrix(NA_real_, nrow = nrow(valid_results[[1]]$out.diff), ncol = ncol(valid_results[[1]]$out.diff),
                 dimnames = dimnames(valid_results[[1]]$out.diff))
        }
      })
      diff_array <- simplify2array(diff_mat)
      agg$diff_matrix <- apply(diff_array, c(1,2), function(xx) mean(xx, na.rm=TRUE))
      # names are inherited from the first valid result
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
