# Helper: Largest Remainder Method for rounding proportions to sum exactly to 1
#' @noRd
round_preserve_sum <- function(x, digits = 3) {
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

# Internal helper to dynamically re-aggregate statistics from raw vectors
#' @noRd
reaggregate_bagging <- function(object, agg_stat = NULL, trim_prob = 0.20) {
  if (is.null(agg_stat)) return(object)

  agg_stat <- match.arg(agg_stat, c("mean", "median", "trimmed", "winsorized", "midhinge", "tukey"))
  if (agg_stat == object$agg_stat_used) return(object)

  agg_helper <- function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)

    switch(agg_stat,
           "mean" = mean(x),
           "median" = median(x),
           "trimmed" = mean(x, trim = trim_prob),
           "winsorized" = {
             q_low <- quantile(x, probs = trim_prob, na.rm = TRUE, names = FALSE)
             q_high <- quantile(x, probs = 1 - trim_prob, na.rm = TRUE, names = FALSE)
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

  raw <- object$raw_stats
  agg <- object$aggregated_stats
  obj_type <- class(object)[1]

  # Re-aggregate OLS
  if (!is.null(raw$ols_tar_coefs)) agg$ols_target[,"estimate"] <- apply(raw$ols_tar_coefs, 2, agg_helper)
  if (!is.null(raw$ols_alt_coefs)) agg$ols_alternative[,"estimate"] <- apply(raw$ols_alt_coefs, 2, agg_helper)

  # Re-aggregate specific DDA statistics (p-values are deliberately ignored)
  if (obj_type == "dda_bagging_indep") {
    agg$hsic_yx_stat <- agg_helper(raw$hsic_yx_stat)
    agg$hsic_xy_stat <- agg_helper(raw$hsic_xy_stat)

    if (!is.null(raw$dcor_yx_stat)) {
      agg$dcor_yx_stat <- agg_helper(raw$dcor_yx_stat)
      agg$dcor_xy_stat <- agg_helper(raw$dcor_xy_stat)
    }
    if (!is.null(raw$bp_yx_stat)) {
      agg$breusch_pagan[[1]]$statistic <- agg_helper(raw$bp_yx_stat)
      agg$breusch_pagan[[1]]$parameter <- agg_helper(raw$bp_yx_df)
      agg$breusch_pagan[[2]]$statistic <- agg_helper(raw$rbp_yx_stat)
      agg$breusch_pagan[[2]]$parameter <- agg_helper(raw$rbp_yx_df)
      agg$breusch_pagan[[3]]$statistic <- agg_helper(raw$bp_xy_stat)
      agg$breusch_pagan[[3]]$parameter <- agg_helper(raw$bp_xy_df)
      agg$breusch_pagan[[4]]$statistic <- agg_helper(raw$rbp_xy_stat)
      agg$breusch_pagan[[4]]$parameter <- agg_helper(raw$rbp_xy_df)
    }
    if (!is.null(raw$nlcor_yx_t1)) {
      agg$nlcor.yx$t1[1:3] <- c(agg_helper(raw$nlcor_yx_t1[,1]), agg_helper(raw$nlcor_yx_t1[,2]), agg_helper(raw$nlcor_yx_t1[,3]))
      agg$nlcor.yx$t2[1:3] <- c(agg_helper(raw$nlcor_yx_t2[,1]), agg_helper(raw$nlcor_yx_t2[,2]), agg_helper(raw$nlcor_yx_t2[,3]))
      agg$nlcor.yx$t3[1:3] <- c(agg_helper(raw$nlcor_yx_t3[,1]), agg_helper(raw$nlcor_yx_t3[,2]), agg_helper(raw$nlcor_yx_t3[,3]))

      agg$nlcor.xy$t1[1:3] <- c(agg_helper(raw$nlcor_xy_t1[,1]), agg_helper(raw$nlcor_xy_t1[,2]), agg_helper(raw$nlcor_xy_t1[,3]))
      agg$nlcor.xy$t2[1:3] <- c(agg_helper(raw$nlcor_xy_t2[,1]), agg_helper(raw$nlcor_xy_t2[,2]), agg_helper(raw$nlcor_xy_t2[,3]))
      agg$nlcor.xy$t3[1:3] <- c(agg_helper(raw$nlcor_xy_t3[,1]), agg_helper(raw$nlcor_xy_t3[,2]), agg_helper(raw$nlcor_xy_t3[,3]))
    }
    if (!is.null(raw$diff_arr)) {
      agg$diff_matrix <- apply(raw$diff_arr, c(1,2), agg_helper)
    }

  } else if (obj_type == "dda_bagging_resdist") {
    agg$agostino.target.statistic <- agg_helper(raw$agost_tar_stat)
    agg$agostino.target.z <- agg_helper(raw$agost_tar_z)
    agg$agostino.alternative.statistic <- agg_helper(raw$agost_alt_stat)
    agg$agostino.alternative.z <- agg_helper(raw$agost_alt_z)
    agg$anscombe.target.statistic <- agg_helper(raw$anscom_tar_stat)
    agg$anscombe.target.z <- agg_helper(raw$anscom_tar_z)
    agg$anscombe.alternative.statistic <- agg_helper(raw$anscom_alt_stat)
    agg$anscombe.alternative.z <- agg_helper(raw$anscom_alt_z)

    for(k in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS3", "RCC", "RHS4")) {
      if (!is.null(raw[[k]])) agg[[k]] <- apply(raw[[k]], 2, agg_helper)
    }

  } else if (obj_type == "dda_bagging_vardist") {
    agg$agostino.predictor.statistic.skew <- agg_helper(raw$agost_pre_stat)
    agg$agostino.predictor.statistic.z <- agg_helper(raw$agost_pre_z)
    agg$agostino.outcome.statistic.skew <- agg_helper(raw$agost_out_stat)
    agg$agostino.outcome.statistic.z <- agg_helper(raw$agost_out_z)
    agg$anscombe.predictor.statistic.kurt <- agg_helper(raw$anscom_pre_stat)
    agg$anscombe.predictor.statistic.z <- agg_helper(raw$anscom_pre_z)
    agg$anscombe.outcome.statistic.kurt <- agg_helper(raw$anscom_out_stat)
    agg$anscombe.outcome.statistic.z <- agg_helper(raw$anscom_out_z)

    for(k in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS", "RCC", "Rtanh")) {
      if (!is.null(raw[[k]])) agg[[k]] <- apply(raw[[k]], 2, agg_helper)
    }
  }

  object$aggregated_stats <- agg
  object$agg_stat_used <- agg_stat
  return(object)
}

# Internal helper to print decision summary
#' @noRd
print_bagging_decisions <- function(object, show = NULL, moment = NULL, type = "generic", digits = 3) {

  # --- Header ---
  header_type <- switch(type,
                        "indep" = "Independence Properties",
                        "resdist" = "Residual Distributions",
                        "vardist" = "Variable Distributions",
                        "Generic")

  cat("\n")
  cat(paste0("BOOTSTRAP AGGREGATED DDA: ", header_type), "\n")
  cat(paste0("Number of Bootstrap Samples: ", object$n_valid_iterations), "\n\n")

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

  # Updated labels to match Docx track changes
  label_map <- list(
    "hsic"          = "HSIC",
    "dcor"          = "dCor",
    "diff_hsic"     = "HSIC Difference",
    "diff_dcor"     = "dCor Difference",
    "diff_mi"       = "MI Difference",
    "dec_bp"        = "Breusch-Pagan",
    "dec_nl.min"    = "Non-linear Correlation",
    "dec_agost"     = "Separate D'Agostino Tests",
    "dec_anscom"    = "Separate Anscombe-Glynn Tests",
    "dec_skewdiff"  = "Skewness Difference",
    "dec_kurtdiff"  = "Kurtosis Difference",
    "dec_cor12diff" = "Co-Skewness Difference",
    "dec_cor13diff" = "Co-Kurtosis Difference",
    "dec_RHS"       = "Hyvarinen-Smith Co-Skewness Difference",
    "dec_RHS3"      = "Hyvarinen-Smith Co-Skewness Difference",
    "dec_RHS4"      = "Hyvarinen-Smith Co-Kurtosis Difference",
    "dec_RCC"       = "Chen-Chan Co-Kurtosis Difference",
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

  # Format string for dynamic decimal places
  fmt <- paste0("%.", digits, "f")

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

    # Apply Largest Remainder Method rounding dynamically based on user digits
    vals <- c(t_val, a_val, u_val)
    vals_rounded <- round_preserve_sum(vals, digits = digits)

    t_val_rd <- vals_rounded[1]
    a_val_rd <- vals_rounded[2]
    u_val_rd <- vals_rounded[3]

    # Output formatting
    if (type == "indep") {
      df_print <- data.frame(
        Target      = sprintf(fmt, t_val_rd),
        Alternative = sprintf(fmt, a_val_rd),
        Confounding = sprintf(fmt, u_val_rd),
        check.names = FALSE
      )
    } else {
      df_print <- data.frame(
        Target      = sprintf(fmt, t_val_rd),
        Alternative = sprintf(fmt, a_val_rd),
        Undecided   = sprintf(fmt, u_val_rd),
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
  } else if (type == "resdist") {
    cat(paste("Note: Target is", varnames[2], "->", varnames[1]), "\n")
    cat(paste("      Alternative is", varnames[1], "->", varnames[2]), "\n")

    probtrans <- if(!is.null(stats$probtrans)) stats$probtrans else FALSE
    if(isTRUE(probtrans)){
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
#' @param agg_stat Method for aggregating test statistics. Options: "mean", "median", "trimmed", "winsorized", "midhinge", "tukey". If NULL, uses the aggregation method established in dda_bagging().
#' @param trim_prob Proportion of observations to be trimmed or Winsorized from each end (default: 0.20).
#' @param digits Number of decimal places to print for proportions (default: 3)
#' @param ... Additional arguments
#' @export
#' @method summary dda_bagging_indep
summary.dda_bagging_indep <- function(object, show = NULL, agg_stat = NULL, trim_prob = 0.20, digits = 3, ...) {
  object <- reaggregate_bagging(object, agg_stat, trim_prob)
  print_bagging_decisions(object, show = show, type = "indep", digits = digits)
  invisible(object)
}

#' Summary for dda_bagging Output (VARDIST)
#'
#' @param object Output from dda_bagging() for dda.vardist objects
#' @param show Character vector of stats to show (e.g. c("skew", "cokurt", "all"))
#' @param moment Numeric vector for moments to include (3, 4)
#' @param agg_stat Method for aggregating test statistics. Options: "mean", "median", "trimmed", "winsorized", "midhinge", "tukey". If NULL, uses the aggregation method established in dda_bagging().
#' @param trim_prob Proportion of observations to be trimmed or Winsorized from each end (default: 0.20).
#' @param digits Number of decimal places to print for proportions (default: 3)
#' @param ... Additional arguments
#' @export
#' @method summary dda_bagging_vardist
summary.dda_bagging_vardist <- function(object, show = NULL, moment = NULL, agg_stat = NULL, trim_prob = 0.20, digits = 3, ...) {
  object <- reaggregate_bagging(object, agg_stat, trim_prob)
  print_bagging_decisions(object, show = show, moment = moment, type = "vardist", digits = digits)
  invisible(object)
}

#' Summary for dda_bagging Output (RESDIST)
#'
#' @param object Output from dda_bagging() for dda.resdist objects
#' @param show Character vector of stats to show (e.g. c("all"))
#' @param moment Numeric vector for moments to include (3, 4)
#' @param agg_stat Method for aggregating test statistics. Options: "mean", "median", "trimmed", "winsorized", "midhinge", "tukey". If NULL, uses the aggregation method established in dda_bagging().
#' @param trim_prob Proportion of observations to be trimmed or Winsorized from each end (default: 0.20).
#' @param digits Number of decimal places to print for proportions (default: 3)
#' @param ... Additional arguments
#' @export
#' @method summary dda_bagging_resdist
summary.dda_bagging_resdist <- function(object,
                                        show = NULL,
                                        moment = NULL,
                                        agg_stat = NULL,
                                        trim_prob = 0.20,
                                        digits = 3, ...) {

  object <- reaggregate_bagging(object, agg_stat, trim_prob)
  print_bagging_decisions(object, show = show, moment = moment, type = "resdist", digits = digits)
  invisible(object)
}
