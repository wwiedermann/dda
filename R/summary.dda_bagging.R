#' Internal helper to print decision summary
#' @noRd
print_bagging_decisions <- function(object, show = NULL, moment = NULL, type = "generic") {
  fn <- object$parameters$function_name
  iter <- object$parameters$iter
  cat(paste0(fn, ", ", iter, " bootstrap aggregations\n\n"))

  decisions <- object$decision_percentages
  if (length(decisions) == 0) {
    cat("No decision proportions available.\n")
    return()
  }

  # --- Mappings ---
  # Map aliases to internal keys for 'show' argument
  alias_map <- list(
    "skew"   = c("dec_agost", "dec_skewdiff"),
    "kurt"   = c("dec_anscom", "dec_kurtdiff"),
    "coskew" = c("dec_cor12diff", "dec_RHS", "dec_RHS3"),
    "cokurt" = c("dec_cor13diff", "dec_RCC", "dec_RHS4", "dec_Rtanh"),
    "hsic"   = c("hsic", "diff_hsic"),
    "dcor"   = c("dcor", "diff_dcor"),
    "mi"     = c("diff_mi")
  )

  # Map Moments to internal keys
  keys_m3 <- c("dec_agost", "dec_skewdiff", "dec_cor12diff", "dec_RHS", "dec_RHS3")
  keys_m4 <- c("dec_anscom", "dec_kurtdiff", "dec_cor13diff", "dec_RCC", "dec_RHS4", "dec_Rtanh")

  # Label Map for Display Titles
  label_map <- list(
    "hsic"          = "HSIC",
    "dcor"          = "dCor",
    "diff_hsic"     = "HSIC Difference",
    "diff_dcor"     = "dCor Difference",
    "diff_mi"       = "MI Difference",
    "dec_agost"     = "Agostino",
    "dec_anscom"    = "Anscombe",
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

  # --- Logic for 'show' and 'moment' ---

  # 1. Resolve 'show' argument
  if (!is.null(show)) {
    expanded_show <- c()
    for (s in show) {
      if (s %in% names(alias_map)) {
        expanded_show <- c(expanded_show, alias_map[[s]])
      } else {
        expanded_show <- c(expanded_show, s)
      }
    }
    # Only keep keys that actually exist in the object
    keys_to_show <- intersect(expanded_show, all_keys)
  }

  # 2. Resolve 'moment' argument
  moment_keys <- c()
  if (!is.null(moment)) {
    if (3 %in% moment) moment_keys <- c(moment_keys, keys_m3)
    if (4 %in% moment) moment_keys <- c(moment_keys, keys_m4)
    moment_keys <- intersect(moment_keys, all_keys)
  }

  # 3. Combine Logic
  if (is.null(show) && is.null(moment)) {
    # Default: Show everything
    keys_to_show <- all_keys
  } else if (is.null(show) && !is.null(moment)) {
    # Only moment specified
    keys_to_show <- moment_keys
  } else if (!is.null(show) && !is.null(moment)) {
    # Union: Show items in 'show' OR items in 'moment'
    keys_to_show <- unique(c(keys_to_show, moment_keys))
  }

  # Ensure print order follows original object order
  keys_to_show <- intersect(all_keys, keys_to_show)

  if (length(keys_to_show) == 0) {
    return()
  }

  for (dname in keys_to_show) {
    prop <- decisions[[dname]]

    # SKIP if all NaNs (stat not calculated/requested in original function)
    if (all(is.nan(prop)) || all(is.na(prop))) next

    # Title
    display_title <- if (!is.null(label_map[[dname]])) label_map[[dname]] else dname
    cat(display_title, "\n")

    # Extract values (names in prop are "Undecided", "Target", "Alternative")
    u_val <- if ("Undecided" %in% names(prop)) prop["Undecided"] else 0
    t_val <- if ("Target" %in% names(prop)) prop["Target"] else 0
    a_val <- if ("Alternative" %in% names(prop)) prop["Alternative"] else 0

    if (type == "indep") {
      # Indep columns: Target, Alternative, Confounding (mapped from Undecided)
      df_print <- data.frame(
        Target      = sprintf("%.2f", t_val),
        Alternative = sprintf("%.2f", a_val),
        Confounding = sprintf("%.2f", u_val),
        check.names = FALSE
      )
    } else {
      # Var/Res columns: Undecided, Target, Alternative
      df_print <- data.frame(
        Undecided   = sprintf("%.2f", u_val),
        Target      = sprintf("%.2f", t_val),
        Alternative = sprintf("%.2f", a_val),
        check.names = FALSE
      )
    }

    print(df_print, row.names = FALSE)
    cat("\n")
  }
}

#' Summary for dda_bagging Output (INDEP)
#'
#' @param object Output from dda_bagging() for dda.indep objects
#' @param show Character vector of stats to show (e.g. c("hsic", "dcor")).
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
#' @param show Character vector of stats to show (e.g. c("skew", "cokurt")).
#' @param moment Numeric vector for moments to include (3, 4).
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
#' @param show Character vector of stats to show.
#' @param moment Numeric vector for moments to include (3, 4).
#' @param ... Additional arguments
#' @export
#' @method summary dda_bagging_resdist
summary.dda_bagging_resdist <- function(object, show = NULL, moment = NULL, ...) {
  print_bagging_decisions(object, show = show, moment = moment, type = "resdist")
  invisible(object)
}
