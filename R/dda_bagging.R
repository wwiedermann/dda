#' Bootstrap Aggregated DDA Analysis
#'
#' @param dda_result Output from any DDA function (dda.indep, dda.vardist, etc.) or CDDA function
#' @param iter Number of bootstrap iterations (default: 100)
#' @param progress Whether to show progress bar (default: TRUE)
#' @param save_file Optional file path to save results
#' @param override_args Named list of arguments to override from original call
#'
#' @return A list containing aggregated results from bootstrap iterations
dda_bagging <- function(dda_result,
                        iter = 100,
                        progress = TRUE,
                        save_file = NULL,
                        override_args = NULL) {

  # Detect object type and extract call information
  if (inherits(dda_result, "dda.indep") ||
      inherits(dda_result, "dda.vardist") ||
      inherits(dda_result, "dda.resdist")) {

    # Handle DDA objects
    if (is.null(dda_result$call_info)) {
      stop("DDA result does not contain function call information. Please ensure you're using an updated DDA function that stores call information.")
    }

    call_info <- dda_result$call_info
    original_call <- call_info$function_call
    function_name <- call_info$function_name
    all_args <- call_info$all_args
    original_formula <- call_info$formula
    data_name <- call_info$data_name
    original_data <- call_info$original_data

    # If data wasn't stored, try to get it from environment
    if (is.null(original_data)) {
      tryCatch({
        original_data <- get(data_name, envir = parent.frame())
      }, error = function(e) {
        stop(paste("Could not retrieve original data:", data_name,
                   "Please ensure the data is available in the environment."))
      })
    }

  } else if (length(dda_result) >= 5 && !is.null(dda_result[[5]])) {

    # Handle CDDA objects (your original structure)
    call_info <- dda_result[[5]]
    original_call <- call_info$function_call
    function_name <- call_info$function_name
    all_args <- call_info$all_args
    original_env <- call_info$environment
    data_name <- call_info$data_name

    # Get the original data
    original_data <- get(data_name, envir = original_env)

  } else {
    stop("Unrecognized DDA result structure. Please ensure you're using a supported DDA function with call information storage.")
  }

  nobs <- nrow(original_data)

  # Get the DDA function
  dda_function <- get(function_name)

  # Override any arguments if specified
  if (!is.null(override_args)) {
    for (arg_name in names(override_args)) {
      all_args[[arg_name]] <- override_args[[arg_name]]
    }
  }

  # Initialize storage
  bagged_results <- list()

  # Progress bar setup
  if (progress) {
    pb <- txtProgressBar(min = 0, max = iter, style = 3, width = 50, char = "=")
    cat(paste("Running", iter, "bootstrap iterations of", function_name, "\n"))
  }

  # Bootstrap loop
  for(i in 1:iter) {

    # Bootstrap sample
    boot_indices <- sample(1:nobs, nobs, replace = TRUE)
    boot_data <- original_data[boot_indices, ]

    # Update the data argument for this iteration
    boot_args <- all_args
    boot_args$data <- boot_data

    # Execute the DDA function with bootstrap data
    tryCatch({
      bagged_results[[i]] <- do.call(dda_function, boot_args)
    }, error = function(e) {
      warning(paste("DDA function failed at iteration", i, ":", e$message))
      bagged_results[[i]] <- NA
    })

    # Update progress bar
    if (progress) {
      setTxtProgressBar(pb, i)
    }
  }

  # Close progress bar
  if (progress) {
    close(pb)
    cat("\nBootstrap iterations completed.\n")
  }

  # Compile output
  output <- list(
    bagged_results = bagged_results,
    original_result = dda_result,
    parameters = list(
      original_call = original_call,
      function_name = function_name,
      object_type = class(dda_result)[1],
      iter = iter,
      successful_iterations = sum(!is.na(bagged_results)),
      failed_iterations = sum(is.na(bagged_results))
    ),
    aggregated_results = NULL  # Will be filled by aggregate function
  )

  # Save results if requested
  if (!is.null(save_file)) {
    save(output, file = save_file)
    cat(paste("Results saved to:", save_file, "\n"))
  }

  return(output)
}

#' Aggregate DDA Bagging Results
#'
#' @param bagging_output Output from dda_bagging function
#' @param method Method for aggregating results ("mean", "median", "harmonic_p")
#' @param alpha Significance level for decisions (default: 0.05)
#'
#' @return List containing aggregated statistics and decision percentages
aggregate_dda_results <- function(bagging_output, method = "mean", alpha = 0.05) {

  bagged_results <- bagging_output$bagged_results
  object_type <- bagging_output$parameters$object_type

  # Remove NA results
  valid_results <- bagged_results[!is.na(bagged_results)]

  if (length(valid_results) == 0) {
    stop("No valid results to aggregate")
  }

  cat(paste("Aggregating", length(valid_results), "valid results out of",
            length(bagged_results), "total iterations\n"))

  # Handle different DDA object types
  if (object_type == "dda.indep") {
    aggregated_results <- aggregate_dda_indep(valid_results, method, alpha)
  } else if (object_type == "dda.vardist") {
    aggregated_results <- aggregate_dda_vardist(valid_results, method, alpha)
  } else if (object_type == "dda.resdist") {
    aggregated_results <- aggregate_dda_resdist(valid_results, method, alpha)
  } else {
    # Generic aggregation for other types
    aggregated_results <- aggregate_generic_dda(valid_results, method, alpha)
  }

  # Update the bagging output with aggregated results
  bagging_output$aggregated_results <- aggregated_results

  return(bagging_output)
}

#' Aggregate dda.indep results specifically
aggregate_dda_indep <- function(valid_results, method = "mean", alpha = 0.05) {

  # Extract HSIC statistics
  hsic_yx <- sapply(valid_results, function(x) x$hsic.yx$statistic)
  hsic_xy <- sapply(valid_results, function(x) x$hsic.xy$statistic)
  hsic_yx_pval <- sapply(valid_results, function(x) x$hsic.yx$p.value)
  hsic_xy_pval <- sapply(valid_results, function(x) x$hsic.xy$p.value)

  # Extract distance correlation statistics
  dcor_yx <- sapply(valid_results, function(x) x$dcor_yx$statistic)
  dcor_xy <- sapply(valid_results, function(x) x$dcor_xy$statistic)
  dcor_yx_pval <- sapply(valid_results, function(x) x$dcor_yx$p.value)
  dcor_xy_pval <- sapply(valid_results, function(x) x$dcor_xy$p.value)

  # Aggregate statistics
  aggregated_stats <- list(
    hsic_yx_stat = mean(hsic_yx, na.rm = TRUE),
    hsic_xy_stat = mean(hsic_xy, na.rm = TRUE),
    dcor_yx_stat = mean(dcor_yx, na.rm = TRUE),
    dcor_xy_stat = mean(dcor_xy, na.rm = TRUE)
  )

  # Aggregate p-values using harmonic mean if requested
  if (method == "harmonic_p") {
    require(harmonicmeanp)
    aggregated_stats$hsic_yx_pval <- p.hmp(hsic_yx_pval[!is.na(hsic_yx_pval)],
                                           L = length(hsic_yx_pval[!is.na(hsic_yx_pval)]))
    aggregated_stats$hsic_xy_pval <- p.hmp(hsic_xy_pval[!is.na(hsic_xy_pval)],
                                           L = length(hsic_xy_pval[!is.na(hsic_xy_pval)]))
    aggregated_stats$dcor_yx_pval <- p.hmp(dcor_yx_pval[!is.na(dcor_yx_pval)],
                                           L = length(dcor_yx_pval[!is.na(dcor_yx_pval)]))
    aggregated_stats$dcor_xy_pval <- p.hmp(dcor_xy_pval[!is.na(dcor_xy_pval)],
                                           L = length(dcor_xy_pval[!is.na(dcor_xy_pval)]))
  } else {
    aggregated_stats$hsic_yx_pval <- mean(hsic_yx_pval, na.rm = TRUE)
    aggregated_stats$hsic_xy_pval <- mean(hsic_xy_pval, na.rm = TRUE)
    aggregated_stats$dcor_yx_pval <- mean(dcor_yx_pval, na.rm = TRUE)
    aggregated_stats$dcor_xy_pval <- mean(dcor_xy_pval, na.rm = TRUE)
  }

  # Calculate decision percentages
  decisions <- data.frame(
    hsic_decision = ifelse(hsic_yx_pval < alpha & hsic_xy_pval >= alpha, "x->y",
                           ifelse(hsic_yx_pval >= alpha & hsic_xy_pval < alpha, "y->x", "undecided")),
    dcor_decision = ifelse(dcor_yx_pval < alpha & dcor_xy_pval >= alpha, "x->y",
                           ifelse(dcor_yx_pval >= alpha & dcor_xy_pval < alpha, "y->x", "undecided"))
  )

  decision_percentages <- list(
    hsic_decisions = prop.table(table(decisions$hsic_decision)),
    dcor_decisions = prop.table(table(decisions$dcor_decision))
  )

  # Handle difference statistics if they exist
  if (!is.null(valid_results[[1]]$out.diff)) {
    diff_estimates <- sapply(valid_results, function(x) x$out.diff[,"estimate"])
    diff_lower <- sapply(valid_results, function(x) x$out.diff[,"lower"])
    diff_upper <- sapply(valid_results, function(x) x$out.diff[,"upper"])

    aggregated_stats$diff_estimates <- apply(diff_estimates, 1, mean, na.rm = TRUE)
    aggregated_stats$diff_lower <- apply(diff_lower, 1, mean, na.rm = TRUE)
    aggregated_stats$diff_upper <- apply(diff_upper, 1, mean, na.rm = TRUE)
  }

  return(list(
    method = method,
    aggregated_stats = aggregated_stats,
    decision_percentages = decision_percentages,
    n_valid_iterations = length(valid_results)
  ))
}

#' Generic aggregation for other DDA types
aggregate_generic_dda <- function(valid_results, method = "mean", alpha = 0.05) {

  # This is a fallback for any DDA types not specifically handled
  # Extract all numeric values from each result
  extract_all_numeric <- function(result) {
    numeric_vals <- list()
    for (i in seq_along(result)) {
      if (is.numeric(result[[i]])) {
        numeric_vals <- c(numeric_vals, result[[i]])
      } else if (is.list(result[[i]])) {
        # Recursively extract from sublists
        numeric_vals <- c(numeric_vals, unlist(result[[i]], use.names = TRUE))
      }
    }
    return(numeric_vals)
  }

  all_numeric <- lapply(valid_results, extract_all_numeric)

  # Find common names across all results
  all_names <- unique(unlist(lapply(all_numeric, names)))

  # Aggregate each statistic
  aggregated_stats <- list()
  for (name in all_names) {
    values <- sapply(all_numeric, function(x) x[[name]])
    values <- values[!is.na(values)]
    if (length(values) > 0) {
      aggregated_stats[[name]] <- mean(values)
    }
  }

  return(list(
    method = method,
    aggregated_stats = aggregated_stats,
    decision_percentages = NULL,
    n_valid_iterations = length(valid_results)
  ))
}
