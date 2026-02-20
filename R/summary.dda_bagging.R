# Helper: Largest Remainder Method for rounding proportions to sum exactly to 1
round_preserve_sum <- function(x, digits = 2) {
  if (sum(x, na.rm = TRUE) == 0) return(x)

  multiplier <- 10^digits
  scaled <- x * multiplier
  rounded <- floor(scaled)
  remainder <- scaled - rounded

  diff <- round(multiplier - sum(rounded))

  if (diff > 0 && diff <= length(x)) {
    # Add 1 to the elements with the largest decimal remainders
    idx <- order(remainder, decreasing = TRUE)[1:diff]
    rounded[idx] <- rounded[idx] + 1
  }

  return(rounded / multiplier)
}

# Internal helper to print decision summary
#' @noRd
print_bagging_decisions <- function(object, show = NULL, moment = NULL, type = "generic") {

  # --- Header ---
  header_type <- switch(type,
                        "indep" = "Independence Properties",
                        "resdist" = "Residual Distributions",
                        "vardist" = "Variable Distributions",
                        "Generic")

  cat("\n")
  cat(paste0("DIRECTION DEPENDENCE ANALYSIS SUMMARY: ", header_type, " (Bagged)"), "\n")
  cat(paste0("Number of Bootstrap Aggregated Results: ", object$n_valid_iterations), "\n\n")

  decisions <- object$decision_percentages
  if (length(decisions) == 0) {
    cat("No decision proportions available.\n")
    return()
  }

  # --- Mappings ---
  alias_map <- list(
    "skew"   = c("dec_agost", "dec_skewdiff"),
    "kurt"   = c("dec_anscom", "dec_kurtdiff"),
    "coskew" = c("dec_cor12diff", "dec_RHS", "dec_RHS3"),
    "cokurt" = c("dec_cor13diff", "dec_RCC", "dec_RHS4", "dec_Rtanh"),
    "hsic"   = c("hsic", "diff_hsic"),
    "dcor"   = c("dcor", "diff_dcor"),
    "mi"     = c("diff_mi"),
    "bp"     = c("dec_bp"),
    "nlcor"  = c("dec_nl.min")
  )

  keys_m3 <- c("dec_agost", "dec_skewdiff", "dec_cor12diff", "dec_RHS", "dec_RHS3")
  keys_m4 <- c("dec_anscom", "dec_kurtdiff", "dec_cor13diff", "dec_RCC", "dec_RHS4", "dec_Rtanh")

  label_map <- list(
    "hsic"          = "HSIC",
    "dcor"          = "dCor",
    "diff_hsic"     = "HSIC Difference",
    "diff_dcor"     = "dCor Difference",
    "diff_mi"       = "MI Difference",
    "dec_bp"        = "Breusch-Pagan",
    "dec_nl.min"    = "Non-linear Correlation (minimum)",
    "dec_agost"     = "Separate D'Agostino Tests",
    "dec_anscom"    = "Separate Anscombe-Glynn Tests",
    "dec_skewdiff"  = "Skewness Difference",
    "dec_kurtdiff"  = "Kurtosis Difference",
    "dec_cor12diff" = "Co-Skewness Difference",
    "dec_cor13diff" = "Co-Kurtosis Difference",
    "dec_RHS"       = "Hyvarinen-Smith (Co-Skewness)",
    "dec_RHS3"      = "Hyvarinen-Smith (Co-Skewness)",
    "dec_RHS4"      = "Hyvarinen-Smith (Co-Kurtosis)",
    "dec_RCC"       = "Chen-Chan (Co-Kurtosis)",
    "dec_Rtanh"     = "Rtanh"
  )

  all_keys <- names(decisions)
  keys_to_show <- c()

  # --- Logic for 'show' (including "all" option) and 'moment' ---
  if (!is.null(show)) {
    if ("all" %in% tolower(show)) {
      keys_to_show <- all_keys
    } else {
      expanded_show <- c()
      for (s in show) {
        if (s %in% names(alias_map)) {
          expanded_show <- c(expanded_show, alias_map[[s]])
        } else {
          expanded_show <- c(expanded_show, s)
        }
      }
      keys_to_show <- intersect(expanded_show, all_keys)
    }
  }

  moment_keys <- c()
  if (!is.null(moment)) {
    if (3 %in% moment) moment_keys <- c(moment_keys, keys_m3)
    if (4 %in% moment) moment_keys <- c(moment_keys, keys_m4)
    moment_keys <- intersect(moment_keys, all_keys)
  }

  if (is.null(show) && is.null(moment)) {
    keys_to_show <- all_keys
  } else if (is.null(show) && !is.null(moment)) {
    keys_to_show <- moment_keys
  } else if (!is.null(show) && !is.null(moment)) {
    keys_to_show <- unique(c(keys_to_show, moment_keys))
  }

  keys_to_show <- intersect(all_keys, keys_to_show)

  if (length(keys_to_show) == 0) {
    if(!is.null(show) || !is.null(moment)) {
      cat("No statistics matched the specified filters.\n")
    }
    return()
  }

  # --- Print Decision Tables ---
  for (dname in keys_to_show) {
    prop <- decisions[[dname]]

    if (all(is.nan(prop)) || all(is.na(prop))) next

    display_title <- if (!is.null(label_map[[dname]])) label_map[[dname]] else dname
    cat(display_title, "\n")

    # Extract exact values
    t_val <- if ("Target" %in% names(prop)) prop["Target"] else 0
    a_val <- if ("Alternative" %in% names(prop)) prop["Alternative"] else 0
    u_val <- if ("Undecided" %in% names(prop)) prop["Undecided"] else 0

    # Apply Largest Remainder Method rounding
    vals <- c(t_val, a_val, u_val)
    vals_rounded <- round_preserve_sum(vals, digits = 2)

    t_val_rd <- vals_rounded[1]
    a_val_rd <- vals_rounded[2]
    u_val_rd <- vals_rounded[3]

    # Output formatting
    if (type == "indep") {
      df_print <- data.frame(
        Target      = sprintf("%.2f", t_val_rd),
        Alternative = sprintf("%.2f", a_val_rd),
        Confounding = sprintf("%.2f", u_val_rd),
        check.names = FALSE
      )
    } else {
      df_print <- data.frame(
        Target      = sprintf("%.2f", t_val_rd),
        Alternative = sprintf("%.2f", a_val_rd),
        Undecided   = sprintf("%.2f", u_val_rd),
        check.names = FALSE
      )
    }

    print(df_print, row.names = FALSE)
    cat("\n")
  }

  # --- Footnote with Actual Variable Names ---
  stats <- object$aggregated_stats
  varnames <- NULL

  # Try to get variable names from aggregated stats
  if (!is.null(stats$var.names) && length(stats$var.names) == 2) {
    varnames <- stats$var.names
  } else {
    # Fallback defaults
    if (type == "indep") {
      varnames <- c("y", "x")
    } else if (type == "resdist") {
      varnames <- c("target", "alternative")
    } else {
      varnames <- c("Outcome", "Predictor")
    }
  }

  # Print footnote using actual variable names
  cat("---\n")

  if (type == "indep") {
    cat(paste("Note: Target is", varnames[2], "->", varnames[1]), "\n")
    cat(paste("      Alternative is", varnames[1], "->", varnames[2]), "\n")
    cat(paste("      Difference statistics > 0 suggest", varnames[2], "->", varnames[1]), "\n")
  } else if (type == "resdist") {
    cat(paste("Note: Target is", varnames[2], "->", varnames[1]), "\n")
    cat(paste("      Alternative is", varnames[1], "->", varnames[2]), "\n")

    probtrans <- if(!is.null(stats$probtrans)) stats$probtrans else FALSE
    if(isFALSE(probtrans)){
      cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1]), "\n")
    } else {
      cat(paste("      Under prob.trans = TRUE, skewness and kurtosis differences < 0 and", "\n",
                "      co-skewness and co-kurtosis differences > 0 suggest", varnames[2], "->", varnames[1]), "\n")
    }
  } else if (type == "vardist") {
    cat(paste("Note: (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests the model", varnames[2], "->", varnames[1]), "\n")
  }
}

#' Summary for dda_bagging Output (INDEP)
#'
#' @param object Output from dda_bagging() for dda.indep objects
#' @param show Character vector of stats to show (e.g. c("hsic", "dcor", "bp", "all"))
#' @param ... Additional arguments
#' @export
#' @method summary dda_bagging_indep
summary.dda_bagging_indep <- function(object, show = NULL, ...) {
  print_bagging_decisions(object, show = show, type = "indep")
  invisible(object)
}

#' Summary for dda_bagging Output (VARDIST)
#'
#' @param object Output from dda_bagging() for dda.vardist objects
#' @param show Character vector of stats to show (e.g. c("skew", "cokurt", "all"))
#' @param moment Numeric vector for moments to include (3, 4)
#' @param ... Additional arguments
#' @export
#' @method summary dda_bagging_vardist
summary.dda_bagging_vardist <- function(object, show = NULL, moment = NULL, ...) {
  print_bagging_decisions(object, show = show, moment = moment, type = "vardist")
  invisible(object)
}

#' Summary for dda_bagging Output (RESDIST)
#'
#' @param object Output from dda_bagging() for dda.resdist objects
#' @param show Character vector of stats to show (e.g. c("all"))
#' @param moment Numeric vector for moments to include (3, 4)
#' @param ... Additional arguments
#' @export
#' @method summary dda_bagging_resdist
summary.dda_bagging_resdist <- function(object, show = NULL, moment = NULL, ...) {
  print_bagging_decisions(object, show = show, moment = moment, type = "resdist")
  invisible(object)
}
