#' Summary for dda_bagging Output (INDEP)
#'
#' @param object Output from dda_bagging() for dda.indep objects (class: dda_bagging_indep)
#' @param ... Additional arguments
#' @export
#' @method summary dda_bagging_indep
summary.dda_bagging_indep <- function(object, ...) {
  decisions <- object$decision_percentages

  if (!is.null(decisions) && length(decisions) > 0) {
    for (dname in names(decisions)) {
      prop <- decisions[[dname]]

      # Prepare title
      cat(paste0("$", dname), "\n")
      cat("x\n")

      # Format into a data frame for clean printing
      # prop is named c("Undecided", "Target", "Alternative")
      df_print <- data.frame(
        Undecided   = sprintf("%.2f", prop["Undecided"]),
        Target      = sprintf("%.2f", prop["Target"]),
        Alternative = sprintf("%.2f", prop["Alternative"]),
        row.names   = "" # suppress row names
      )

      print(df_print)
      cat("\n")
    }
  } else {
    cat("No decision proportions available.\n")
  }
  invisible(object)
}

#' Summary for dda_bagging Output (VARDIST)
#'
#' @param object Output from dda_bagging() for dda.vardist objects (class: dda_bagging_vardist)
#' @param ... Additional arguments
#' @export
#' @method summary dda_bagging_vardist
summary.dda_bagging_vardist <- function(object, ...) {
  decisions <- object$decision_percentages

  if (!is.null(decisions) && length(decisions) > 0) {
    for (dname in names(decisions)) {
      prop <- decisions[[dname]]

      cat(paste0("$", dname), "\n")
      cat("x\n")

      df_print <- data.frame(
        Undecided   = sprintf("%.2f", prop["Undecided"]),
        Target      = sprintf("%.2f", prop["Target"]),
        Alternative = sprintf("%.2f", prop["Alternative"]),
        row.names   = ""
      )

      print(df_print)
      cat("\n")
    }
  } else {
    cat("No decision proportions available.\n")
  }
  invisible(object)
}

#' Summary for dda_bagging Output (RESDIST)
#'
#' @param object Output from dda_bagging() for dda.resdist objects (class: dda_bagging_resdist)
#' @param ... Additional arguments
#' @export
#' @method summary dda_bagging_resdist
summary.dda_bagging_resdist <- function(object, ...) {
  decisions <- object$decision_percentages

  if (!is.null(decisions) && length(decisions) > 0) {
    for (dname in names(decisions)) {
      prop <- decisions[[dname]]

      cat(paste0("$", dname), "\n")
      cat("x\n")

      df_print <- data.frame(
        Undecided   = sprintf("%.2f", prop["Undecided"]),
        Target      = sprintf("%.2f", prop["Target"]),
        Alternative = sprintf("%.2f", prop["Alternative"]),
        row.names   = ""
      )

      print(df_print)
      cat("\n")
    }
  } else {
    cat("No decision proportions available.\n")
  }
  invisible(object)
}
