data("mtcars")
result <- dda.indep(mpg ~ wt + hp, pred = "wt", data = mtcars)
bagged_result <- dda_bagging(result, iter = 10)
summary.dda_bagging(bagged_result)

# result <- cdda.indep(mpg ~ wt * hp, pred = "wt", mod = "hp",
#                      modval = "mean", data = mtcars)
# bagged_result <- dda_bagging(result, iter = 10)
#
# result <- dda.resdist(mpg ~ wt + hp, pred = "wt", data = mtcars)
# bagged_result <- dda_bagging(result, iter = 10)
# summary.dda_bagging(bagged_result)


#' Bootstrap Aggregated DDA Analysis
#'
#' @param dda_result Output from any DDA function (dda.indep, dda.vardist, etc.) or CDDA function
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

dda_bagging <- function(
    dda_result,
    iter = 100,
    progress = TRUE,
    save_file = NULL,     ### Eh?
    override_args = NULL, ### Also unsure about this argument
    method = "mean",
    alpha = 0.05
) {

  # Helper: robustly extract first numeric value from any R object
  get_numeric <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (is.numeric(x)) return(as.numeric(x[1]))
    if (is.list(x)) return(get_numeric(x[[1]])) # recursively check first element
    return(NA_real_)
  }

  # Detect object type and extract call information
  if (inherits(dda_result, "dda.indep") ||
      inherits(dda_result, "dda.vardist") ||
      inherits(dda_result, "dda.resdist")) {

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
    call_info <- dda_result[[5]]
    function_name <- call_info$function_name
    all_args <- call_info$all_args
    original_env <- call_info$environment
    data_name <- call_info$data_name
    original_data <- get(data_name, envir = original_env)
  } else {
    stop("Unrecognized DDA result structure. Please ensure you're using a supported DDA function with call information storage.")
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
    cat(paste("Running", iter, "bootstrap iterations of", function_name, "\n"))
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

  if (n_valid > 0) {
    if (object_type == "dda.indep") {
      hsic_yx <- sapply(valid_results, function(x) get_numeric(x$hsic.yx$statistic))
      hsic_xy <- sapply(valid_results, function(x) get_numeric(x$hsic.xy$statistic))
      hsic_yx_pval <- sapply(valid_results, function(x) get_numeric(x$hsic.yx$p.value))
      hsic_xy_pval <- sapply(valid_results, function(x) get_numeric(x$hsic.xy$p.value))
      dcor_yx <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_yx$statistic))
      dcor_xy <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_xy$statistic))
      dcor_yx_pval <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_yx$p.value))
      dcor_xy_pval <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_xy$p.value))

      aggregate_numeric <- function(vals, method) {
        vals <- as.numeric(vals)
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) return(NA_real_)
        if (method == "median") return(median(vals))
        mean(vals) # default mean
      }

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

      # Decision percentages
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
    } else {
      # Generic aggregation
      extract_all_numeric <- function(result) {
        numeric_vals <- list()
        for (i in seq_along(result)) {
          if (is.numeric(result[[i]])) {
            numeric_vals <- c(numeric_vals, result[[i]])
          } else if (is.list(result[[i]])) {
            numeric_vals <- c(numeric_vals, unlist(result[[i]], use.names = TRUE))
          }
        }
        return(numeric_vals)
      }
      all_numeric <- lapply(valid_results, extract_all_numeric)
      all_names <- unique(unlist(lapply(all_numeric, names)))
      for (name in all_names) {
        vals <- sapply(all_numeric, function(x) {
          val <- x[[name]]
          get_numeric(val)
        })
        agg[[name]] <- aggregate_numeric(vals, method)
      }
    }
  } else {
    warning("No valid results to aggregate")
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

  # Helper for pretty printing non-NA, non-numeric(0), safe for all types
  clean_print <- function(x, digits) {
    if (is.null(x)) return(NULL)
    # If it's all NA or empty numeric, skip
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
      # Filter out nulls/all NA/empty numerics recursively
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
    # For atomic types and character
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
