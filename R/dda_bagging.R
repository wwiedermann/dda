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

      # Initialize with 0 for all expected categories
      out <- c("Undecided" = 0, "Target" = 0, "Alternative" = 0)

      # Map internal names to display names
      # internal: "undecided", "x->y", "y->x"
      # display:  "Undecided", "Target", "Alternative"

      if (!is.null(prop) && length(prop) > 0) {
        if ("undecided" %in% names(prop)) out["Undecided"] <- prop["undecided"]
        if ("x->y" %in% names(prop))      out["Target"]    <- prop["x->y"]
        if ("y->x" %in% names(prop))      out["Alternative"] <- prop["y->x"]
      }

      # Mapping internal stat names to display headers
      display_name <- switch(dname,
                             "hsic" = "$hsic",
                             "dcor" = "$dcor",
                             dname)

      cat(display_name, "\n")

      # Create a clean data.frame for printing
      # Using data.frame handles column alignment automatically
      df_print <- data.frame(
        Undecided = sprintf("%.2f", out["Undecided"]),
        Target    = sprintf("%.2f", out["Target"]),
        Alternative = sprintf("%.2f", out["Alternative"]),
        row.names = "x"
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

      out <- c("Undecided" = 0, "Target" = 0, "Alternative" = 0)

      if (!is.null(prop) && length(prop) > 0) {
        if ("undecided" %in% names(prop)) out["Undecided"] <- prop["undecided"]
        if ("x->y" %in% names(prop))      out["Target"]    <- prop["x->y"]
        if ("y->x" %in% names(prop))      out["Alternative"] <- prop["y->x"]
      }

      cat(paste0("$", dname), "\n")

      df_print <- data.frame(
        Undecided = sprintf("%.2f", out["Undecided"]),
        Target    = sprintf("%.2f", out["Target"]),
        Alternative = sprintf("%.2f", out["Alternative"]),
        row.names = "x"
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
      out <- c("Undecided" = 0, "Target" = 0, "Alternative" = 0)

      if (!is.null(prop) && length(prop) > 0) {
        if ("undecided" %in% names(prop)) out["Undecided"] <- prop["undecided"]
        if ("x->y" %in% names(prop))      out["Target"]    <- prop["x->y"]
        if ("y->x" %in% names(prop))      out["Alternative"] <- prop["y->x"]
      }

      cat(paste0("$", dname), "\n")

      df_print <- data.frame(
        Undecided = sprintf("%.2f", out["Undecided"]),
        Target    = sprintf("%.2f", out["Target"]),
        Alternative = sprintf("%.2f", out["Alternative"]),
        row.names = "x"
      )

      print(df_print)
      cat("\n")
    }
  } else {
    cat("No decision proportions available.\n")
  }
  invisible(object)
}
