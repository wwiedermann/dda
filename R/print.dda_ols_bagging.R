#' Print OLS Summary from Bagged DDA
#'
#' @param object Output from dda_bagging()
#' @param agg_stat Method for aggregating test statistics. If NULL, uses the method established in dda_bagging().
#' @param trim_prob Proportion of observations to be trimmed (default: 0.10).
#' @param win_prob Proportion of observations to be Winsorized (default: 0.10).
#' @param digits Number of digits for rounding
#' @param ... Additional arguments passed to print
#' @export
print_ols_summary <- function(object,
                              agg_stat = NULL,
                              trim_prob = 0.10,
                              win_prob = 0.10,
                              digits = 4,
                              ...) {

  if (!inherits(object, "dda_bagging")) {
    stop("Object must be a bagged DDA result.")
  }

  # Ensure dynamic re-aggregation of coefficients occurs with split trim/win probs
  object <- reaggregate_bagging(object, agg_stat, trim_prob, win_prob)
  stats <- object$aggregated_stats
  raw <- object$raw_stats

  if (is.null(stats$ols_target)) {
    cat("No OLS summary available in this object.\n")
    return(invisible(NULL))
  }

  # Local helper for R-squared aggregation
  current_agg <- if (!is.null(agg_stat)) agg_stat else object$agg_stat_used

  agg_helper <- function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    switch(current_agg,
           "mean" = mean(x),
           "median" = median(x),
           "trimmed" = mean(x, trim = trim_prob),
           "winsorized" = {
             q_low <- quantile(x, probs = win_prob, na.rm = TRUE, names = FALSE)
             q_high <- quantile(x, probs = 1 - win_prob, na.rm = TRUE, names = FALSE)
             x[x < q_low] <- q_low
             x[x > q_high] <- q_high
             mean(x)
           },
           "midhinge" = {
             q <- quantile(x, probs = c(0.25, 0.75), names = FALSE)
             mean(q)
           },
           "tukey" = {
             q <- quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE)
             (q[1] + 2*q[2] + q[3]) / 4
           }
    )
  }

  # Print Aggregation Method
  cat("\nAggregation method:", current_agg, "\n\n")

  # --- Target Model Print ---
  cat("OLS Summary: Target Model\n")
  print.default(format(stats$ols_target, digits = digits), print.gap = 2L, quote = FALSE)

  if (!is.null(raw$ols_tar_rsq)) {
    r2 <- agg_helper(raw$ols_tar_rsq[, 1])
    adj_r2 <- agg_helper(raw$ols_tar_rsq[, 2])
    cat(sprintf("\nR-squared: %.*f, Adjusted R-squared: %.*f\n", digits, r2, digits, adj_r2))
  }

  cat("\n---\n\n")

  # --- Alternative Model Print ---
  cat("OLS Summary: Alternative Model\n")
  print.default(format(stats$ols_alternative, digits = digits), print.gap = 2L, quote = FALSE)

  if (!is.null(raw$ols_alt_rsq)) {
    r2_alt <- agg_helper(raw$ols_alt_rsq[, 1])
    adj_r2_alt <- agg_helper(raw$ols_alt_rsq[, 2])
    cat(sprintf("\nR-squared: %.*f, Adjusted R-squared: %.*f\n", digits, r2_alt, digits, adj_r2_alt))
  }

  invisible(object)
}
