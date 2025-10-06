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
summary.dda_bagging(bagged_indep)

# ===== DDA Bagging Summary =====
# Function: dda.indep
# Object Type: dda.indep
# Iterations: 10
# Valid Iterations: 10
# ----
#
# Aggregated Statistics:
# hsic_yx_stat : 0.6324
# hsic_xy_stat : 4.9462
# hsic_yx_pval : 0.1094
# hsic_xy_pval : 0
# diff_estimates : 4.3139, 0.1809, 0.1142
# diff_lower : 2.4641, 0.0974, 0.0614
# diff_upper : 6.1486, 0.2171, 0.1765
#
# Decision Percentages (alpha = 0.05 ):
#   hsic :
#   hsic_decision
# undecided      y->x
# 0.4       0.6

# 2. dda.resdist: Direction Dependence Analysis - Residual Distribution

# result_resdist <- dda.resdist(mpg ~ wt + hp, pred = "wt", data = mtcars)
result_resdist <- dda.resdist(y ~ x, pred = "x", data = d,
                      B = 50, conf.level = 0.90)


print(result_resdist)

bagged_resdist <- dda_bagging(result_resdist, iter = 100)
summary.dda_bagging(bagged_resdist)


# 3. dda.vardist: Direction Dependence Analysis - Variable Distribution
result_vardist <- dda.vardist(mpg ~ wt + hp, pred = "wt", data = mtcars)
print(result_vardist)

bagged_vardist <- dda_bagging(result_vardist, iter = 10)
summary.dda_bagging(bagged_vardist)



#' Bootstrap Aggregated DDA Analysis
#'
#' @param dda_result Output from any DDA function (dda.indep, dda.vardist, dda.resdist, etc.)
#' @param iter Number of bootstrap iterations (default: 100)
#' @param progress Whether to show progress bar (default: TRUE)
#' @param save_file Optional file path to save results
#' @param override_args Named list of arguments to override from original call
#' @param method Method for aggregating results ("mean", "median", "harmonic_p")
#' @param alpha Significance level for decisions (default: 0.05)
#'
#' @rdname dda.indep
#' @export
#' @return A list containing bootstrap and aggregated results
# NA-safe DDA Bagging for dda.resdist, dda.indep, etc.
# NA-safe DDA Bagging for dda.resdist, dda.indep, etc.

dda_bagging <- function(
    dda_result,
    iter = 100,
    progress = TRUE, # Show progress bar
    #save_file = NULL, # Eh?
    #override_args = NULL, # Eh?
    method = "mean",
    alpha = 0.05
) {
  get_numeric <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (is.numeric(x)) return(as.numeric(x[1]))
    if (is.list(x)) return(get_numeric(x[[1]]))
    return(NA_real_)
  }

  # Detect object type and extract call information
  if (inherits(dda_result, "dda.indep") ||
      inherits(dda_result, "dda.vardist") ||
      inherits(dda_result, "dda.resdist") ||
      inherits(dda_result, "cdda.indep") ||
      inherits(dda_result, "cdda.vardist") ){

    if (is.null(dda_result$call_info)) {
      stop("DDA result does not contain function call information. Please ensure you're using an updated DDA function that stores call information.")
    }

    call_info <- dda_result$call_info
    function_name <- call_info$function_name
    all_args <- call_info$all_args
    data_name <- call_info$data_name
    original_data <- call_info$original_data
    if (is.null(original_data)) {
      tryCatch({
        original_data <- get(data_name, envir = parent.frame())
      }, error = function(e) {
        stop(paste("Could not retrieve original data:", data_name,
                   "Please ensure the data is available in the environment."))
      })
    }
  } else if (length(dda_result) >= 5 && !is.null(dda_result[[5]])) {
    call_info <- dda_result[[5]] #[[5]] is function_call, function_name, all_args, formula, data_name, and original_data
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

  if (!is.null(override_args)) {
    for (arg_name in names(override_args)) {
      all_args[[arg_name]] <- override_args[[arg_name]]
    }
  }

  bagged_results <- list()
  if (progress) {
    pb <- txtProgressBar(min = 0, max = iter, style = 3, width = 50, char = "=")
    cat(paste("\n", "Running", iter, "bootstrap iterations of", function_name, "\n"))
  }

  for(i in 1:iter) { #argument "resamplin
    boot_indices <- sample(1:nobs, nobs, replace = TRUE) #options for subsampling
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

  aggregate_numeric <- function(vals, method) {
    vals <- as.numeric(vals)
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NA_real_)
    if (method == "median") return(median(vals))
    mean(vals)
  }

  ## --- DDA.INDEP block ---
  if (n_valid > 0 && object_type == "dda.indep") {
    # If method is mean or median
    hsic_yx <- sapply(valid_results, function(x) get_numeric(x$hsic.yx$statistic))
    hsic_xy <- sapply(valid_results, function(x) get_numeric(x$hsic.xy$statistic))
    hsic_yx_pval <- sapply(valid_results, function(x) get_numeric(x$hsic.yx$p.value))
    hsic_xy_pval <- sapply(valid_results, function(x) get_numeric(x$hsic.xy$p.value))
    dcor_yx <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_yx$statistic))
    dcor_xy <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_xy$statistic))
    dcor_yx_pval <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_yx$p.value))
    dcor_xy_pval <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_xy$p.value))

    agg$hsic_yx_stat <- aggregate_numeric(hsic_yx, method)
    agg$hsic_xy_stat <- aggregate_numeric(hsic_xy, method)
    agg$dcor_yx_stat <- aggregate_numeric(dcor_yx, method)
    agg$dcor_xy_stat <- aggregate_numeric(dcor_xy, method)

    if (method == "harmonic_p") {
      if (!requireNamespace("harmonicmeanp", quietly = TRUE)) {
        warning("harmonicmeanp package not available, using mean for p-values")
        agg$hsic_yx_pval <- aggregate_numeric(hsic_yx_pval, "mean")
        agg$hsic_xy_pval <- aggregate_numeric(hsic_xy_pval, "mean")
        agg$dcor_yx_pval <- aggregate_numeric(dcor_yx_pval, "mean")
        agg$dcor_xy_pval <- aggregate_numeric(dcor_xy_pval, "mean")
      } else {
        agg$hsic_yx_pval <- harmonicmeanp::p.hmp(hsic_yx_pval[!is.na(hsic_yx_pval)], L = sum(!is.na(hsic_yx_pval)))
        agg$hsic_xy_pval <- harmonicmeanp::p.hmp(hsic_xy_pval[!is.na(hsic_xy_pval)], L = sum(!is.na(hsic_xy_pval)))
        agg$dcor_yx_pval <- harmonicmeanp::p.hmp(dcor_yx_pval[!is.na(dcor_yx_pval)], L = sum(!is.na(dcor_yx_pval)))
        agg$dcor_xy_pval <- harmonicmeanp::p.hmp(dcor_xy_pval[!is.na(dcor_xy_pval)], L = sum(!is.na(dcor_xy_pval)))
      }
    } else {
      agg$hsic_yx_pval <- aggregate_numeric(hsic_yx_pval, method)
      agg$hsic_xy_pval <- aggregate_numeric(hsic_xy_pval, method)
      agg$dcor_yx_pval <- aggregate_numeric(dcor_yx_pval, method)
      agg$dcor_xy_pval <- aggregate_numeric(dcor_xy_pval, method)
    }

    hsic_decision <- ifelse(hsic_yx_pval < alpha & hsic_xy_pval >= alpha, "x->y",
                            ifelse(hsic_yx_pval >= alpha & hsic_xy_pval < alpha, "y->x", "undecided"))
    dcor_decision <- ifelse(dcor_yx_pval < alpha & dcor_xy_pval >= alpha, "x->y",
                            ifelse(dcor_yx_pval >= alpha & dcor_xy_pval < alpha, "y->x", "undecided"))
    decisions$hsic <- prop.table(table(hsic_decision))
    decisions$dcor <- prop.table(table(dcor_decision))

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
      agg$diff_estimates <- apply(diff_estimates, 1, aggregate_numeric, method)
      agg$diff_lower <- apply(diff_lower, 1, aggregate_numeric, method)
      agg$diff_upper <- apply(diff_upper, 1, aggregate_numeric, method)
      rownames_matrix <- rownames(valid_results[[1]]$out.diff)
      names(agg$diff_estimates) <- rownames_matrix
      names(agg$diff_lower) <- rownames_matrix
      names(agg$diff_upper) <- rownames_matrix
    }
  }

  ## --- DDA.RESDIST block ---
  if (n_valid > 0 && object_type == "dda.resdist") {
    .get_val <- function(x, key) {
      val <- tryCatch(x[[key]], error = function(e) NA_real_)
      get_numeric(val)
    }
    keys <- c("agostino.alternative.statistic.skew",
              "agostino.alternative.statistic.z",
              "agostino.alternative.p.value",
              "agostino.target.statistic.skew",
              "agostino.target.statistic.z",
              "agostino.target.p.value",
              "anscombe.alternative.statistic.kurt",
              "anscombe.alternative.statistic.z",
              "anscombe.alternative.p.value",
              "anscombe.target.statistic.kurt",
              "anscombe.target.statistic.z",
              "anscombe.target.p.value",
              "skewdiff1",
              "skewdiff.z.value",
              "skewdiff.p.value",
              "skewdiff.lower",
              "skewdiff.upper",
              "kurtdiff1",
              "kurtdiff.z.value",
              "kurtdiff.p.value",
              "kurtdiff.lower",
              "kurtdiff.upper",
              "cor12diff.cor21.diff",
              "cor12diff.lower",
              "cor12diff.upper",
              "cor13diff.cor13.diff",
              "cor13diff.lower",
              "cor13diff.upper",
              "RHS3.RHS3",
              "RHS3.lower",
              "RHS3.upper",
              "RCC.RCC",
              "RCC.lower",
              "RCC.upper",
              "RHS4.RHS4",
              "RHS4.lower",
              "RHS4.upper"
    )
    for (key in keys) {
      vals <- sapply(valid_results, function(x) .get_val(x, key))
      agg[[key]] <- aggregate_numeric(vals, method)
    }

    pval_keys <- c("agostino.alternative.p.value",
                   "agostino.target.p.value",
                   "anscombe.alternative.p.value",
                   "anscombe.target.p.value",
                   "skewdiff.p.value",
                   "kurtdiff.p.value")
    zval_keys <- c("agostino.alternative.statistic.z",
                   "agostino.target.statistic.z",
                   "anscombe.alternative.statistic.z",
                   "anscombe.target.statistic.z",
                   "skewdiff.z.value",
                   "kurtdiff.z.value")
    for (i in seq_along(pval_keys)) {
      if (!is.na(agg[[zval_keys[i]]])) {
        agg[[pval_keys[i]]] <- 2 * pnorm(abs(agg[[zval_keys[i]]]), lower.tail = FALSE)
      }
    }

    z_crit <- qnorm(1 - alpha/2)
    dec_agost <- sapply(valid_results, function(x) {
      z_alt <- .get_val(x, "agostino.alternative.statistic.z")
      z_tgt <- .get_val(x, "agostino.target.statistic.z")
      if (is.na(z_alt) || is.na(z_tgt)) {
        "undecided"
      } else if (abs(z_alt) >= z_crit && abs(z_tgt) < z_crit) {
        "x->y"
      } else if (abs(z_alt) < z_crit && abs(z_tgt) >= z_crit) {
        "y->x"
      } else {
        "undecided"
      }
    })
    dec_anscom <- sapply(valid_results, function(x) {
      z_alt <- .get_val(x, "anscombe.alternative.statistic.z")
      z_tgt <- .get_val(x, "anscombe.target.statistic.z")
      if (is.na(z_alt) || is.na(z_tgt)) {
        "undecided"
      } else if (abs(z_alt) >= z_crit && abs(z_tgt) < z_crit) {
        "x->y"
      } else if (abs(z_alt) < z_crit && abs(z_tgt) >= z_crit) {
        "y->x"
      } else {
        "undecided"
      }
    })
    dec_skewdiff <- sapply(valid_results, function(x) {
      l <- .get_val(x, "skewdiff.lower")
      u <- .get_val(x, "skewdiff.upper")
      if (is.na(l) || is.na(u)) {
        "undecided"
      } else if (l > 0 & u > 0) {
        "x->y"
      } else if (l < 0 & u < 0) {
        "y->x"
      } else {
        "undecided"
      }
    })
    dec_kurtdiff <- sapply(valid_results, function(x) {
      l <- .get_val(x, "kurtdiff.lower")
      u <- .get_val(x, "kurtdiff.upper")
      if (is.na(l) || is.na(u)) {
        "undecided"
      } else if (l > 0 & u > 0) {
        "x->y"
      } else if (l < 0 & u < 0) {
        "y->x"
      } else {
        "undecided"
      }
    })
    dec_cor12diff <- sapply(valid_results, function(x) {
      l <- .get_val(x, "cor12diff.lower")
      u <- .get_val(x, "cor12diff.upper")
      if (is.na(l) || is.na(u)) {
        "undecided"
      } else if (l > 0 & u > 0) {
        "x->y"
      } else if (l < 0 & u < 0) {
        "y->x"
      } else {
        "undecided"
      }
    })
    dec_cor13diff <- sapply(valid_results, function(x) {
      l <- .get_val(x, "cor13diff.lower")
      u <- .get_val(x, "cor13diff.upper")
      if (is.na(l) || is.na(u)) {
        "undecided"
      } else if (l > 0 & u > 0) {
        "x->y"
      } else if (l < 0 & u < 0) {
        "y->x"
      } else {
        "undecided"
      }
    })
    dec_RHS3 <- sapply(valid_results, function(x) {
      l <- .get_val(x, "RHS3.lower")
      u <- .get_val(x, "RHS3.upper")
      if (is.na(l) || is.na(u)) {
        "undecided"
      } else if (l > 0 & u > 0) {
        "x->y"
      } else if (l < 0 & u < 0) {
        "y->x"
      } else {
        "undecided"
      }
    })
    dec_RCC <- sapply(valid_results, function(x) {
      l <- .get_val(x, "RCC.lower")
      u <- .get_val(x, "RCC.upper")
      if (is.na(l) || is.na(u)) {
        "undecided"
      } else if (l > 0 & u > 0) {
        "x->y"
      } else if (l < 0 & u < 0) {
        "y->x"
      } else {
        "undecided"
      }
    })
    dec_RHS4 <- sapply(valid_results, function(x) {
      l <- .get_val(x, "RHS4.lower")
      u <- .get_val(x, "RHS4.upper")
      if (is.na(l) || is.na(u)) {
        "undecided"
      } else if (l > 0 & u > 0) {
        "x->y"
      } else if (l < 0 & u < 0) {
        "y->x"
      } else {
        "undecided"
      }
    })

    decisions$agostino <- prop.table(table(dec_agost))
    decisions$anscombe <- prop.table(table(dec_anscom))
    decisions$skewdiff <- prop.table(table(dec_skewdiff))
    decisions$kurtdiff <- prop.table(table(dec_kurtdiff))
    decisions$cor12diff <- prop.table(table(dec_cor12diff))
    decisions$cor13diff <- prop.table(table(dec_cor13diff))
    decisions$RHS3 <- prop.table(table(dec_RHS3))
    decisions$RCC <- prop.table(table(dec_RCC))
    decisions$RHS4 <- prop.table(table(dec_RHS4))
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
  return(output)
}
#' Summary for dda_bagging Output
#'
#' @param object Output from dda_bagging()
#' @param digits Number of digits for rounding (default: 4)
#' @param alpha Significance level for decisions (default: 0.05)
#'
#' @method summary dda_bagging
#' @export
#'
#' @return Prints an organized summary of aggregated statistics and decision percentages

summary.dda_bagging <- function(object, digits = 4, alpha = 0.05) {
  if (is.null(object$aggregated_stats)) {
    cat("No aggregated statistics found.\n")
    return(invisible(NULL))
  }
  stats <- object$aggregated_stats
  decisions <- object$decision_percentages

  clean_print <- function(x, digits) {
    if (is.null(x)) return(NULL)
    if (is.numeric(x)) {
      if (length(x) == 0) return(NULL)
      if (all(is.na(x))) return(NULL)
      x <- x[!is.na(x)]
      if (length(x) == 0) return(NULL)
      return(paste(round(x, digits), collapse = ", "))
    }
    if (is.matrix(x) || is.data.frame(x)) {
      x_nona <- x[!is.na(x)]
      if (length(x_nona) == 0) return(NULL)
      return(capture.output(print(round(x, digits))))
    }
    if (is.list(x)) {
      x_filt <- Filter(function(y) {
        if (is.null(y)) return(FALSE)
        if (is.numeric(y) && (length(y) == 0 || all(is.na(y)))) return(FALSE)
        if (is.matrix(y) || is.data.frame(y)) {
          y_nona <- y[!is.na(y)]
          return(length(y_nona) > 0)
        }
        TRUE
      }, x)
      if (length(x_filt) == 0) return(NULL)
      return(capture.output(str(x_filt)))
    }
    if (is.atomic(x) && length(x) == 0) return(NULL)
    if (is.atomic(x) && all(is.na(x))) return(NULL)
    return(as.character(x))
  }

  cat("\n===== DDA Bagging Summary =====\n")
  cat("Function:", object$parameters$function_name, "\n")
  cat("Object Type:", object$parameters$object_type, "\n")
  cat("Iterations:", object$parameters$iter, "\n")
  cat("Valid Iterations:", object$n_valid_iterations, "\n")
  cat("----\n")

  cat("\nAggregated Statistics:\n")
  for (stat_name in names(stats)) {
    out <- clean_print(stats[[stat_name]], digits)
    if (!is.null(out)) {
      cat(stat_name, ":", out, "\n")
    }
  }

  if (!is.null(decisions) && length(decisions) > 0) {
    cat("\nDecision Percentages (alpha =", alpha, "):\n")
    for (dec_name in names(decisions)) {
      out <- clean_print(decisions[[dec_name]], digits)
      if (!is.null(out)) {
        cat(dec_name, ":\n")
        print(round(decisions[[dec_name]], digits))
        cat("\n")
      }
    }
  }

  invisible(object)
}
