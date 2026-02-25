#' Print OLS Summary from Bagged DDA
#'
#' @param object Output from dda_bagging()
#' @param digits Number of digits for rounding (default: 4)
#' @export
print_ols_summary <- function(object, digits = 4) {

  if (!inherits(object, "dda_bagging")) {
    stop("Object must be a bagged DDA result.")
  }

  stats <- object$aggregated_stats

  if (is.null(stats$ols_target)) {
    cat("No OLS summary available in this object.\n")
    return(invisible(NULL))
  }

  cat("OLS Summary: Target Model\n")
  print.default(format(stats$ols_target, digits = digits), print.gap = 2L, quote = FALSE)

  cat("\nOLS Summary: Alternative Model\n")
  print.default(format(stats$ols_alternative, digits = digits), print.gap = 2L, quote = FALSE)

  invisible(object)
}
