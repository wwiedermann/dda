#' Summary for dda_bagging Output (INDEP, robust S3, no type check)
#'
#' @param object Output from dda_bagging() for dda.indep objects (class: dda_bagging_indep)
#' @param digits Number of digits for rounding (default: 4)
#' @param alpha Significance level for decisions (default: 0.05)
#' @export
#' @method summary dda_bagging_indep

summary.dda_bagging_indep <- function(object, digits = 4, alpha = 0.05) {
  stats <- object$aggregated_stats
  decisions <- object$decision_percentages

  cat("\n===== DDA Bagging Summary (INDEP) =====\n")
  cat("Function:", object$parameters$function_name, "\n")
  cat("Object Type:", object$parameters$object_type, "\n")
  cat("Iterations:", object$parameters$iter, "\n")
  cat("Completed Iterations:", object$n_valid_iterations, "\n")
  cat("----\n")

  # Build stat table even if some stats are NA
  stat_tab <- matrix(NA, nrow = 2, ncol = 4)
  rownames(stat_tab) <- c("HSIC", "dCor")
  colnames(stat_tab) <- c("Target Stat", "Target p", "Alternative Stat", "Alternative p")

  stat_tab[1, ] <- c(
    if (!is.null(stats$hsic_yx_stat)) round(stats$hsic_yx_stat, digits) else NA,
    if (!is.null(stats$hsic_yx_pval)) round(stats$hsic_yx_pval, digits) else NA,
    if (!is.null(stats$hsic_xy_stat)) round(stats$hsic_xy_stat, digits) else NA,
    if (!is.null(stats$hsic_xy_pval)) round(stats$hsic_xy_pval, digits) else NA
  )
  stat_tab[2, ] <- c(
    if (!is.null(stats$dcor_yx_stat)) round(stats$dcor_yx_stat, digits) else NA,
    if (!is.null(stats$dcor_yx_pval)) round(stats$dcor_yx_pval, digits) else NA,
    if (!is.null(stats$dcor_xy_stat)) round(stats$dcor_xy_stat, digits) else NA,
    if (!is.null(stats$dcor_xy_pval)) round(stats$dcor_xy_pval, digits) else NA
  )

  # Print only rows with any non-NA value
  keep_rows <- apply(stat_tab, 1, function(row) any(!is.na(row) & !is.nan(row)))
  stat_tab_print <- stat_tab[keep_rows, , drop = FALSE]

  if (nrow(stat_tab_print) > 0) {
    cat("\nHSIC and dCor Test Statistics & Harmonic p-values:\n")
    print(stat_tab_print)
    cat("---\n")
  } else {
    cat("\nNo HSIC/dCor statistics found in bagged output.\n")
  }

  # Print difference statistics if present
  if (!is.null(stats$diff_matrix)) {
    diffmat <- stats$diff_matrix
    colnames(diffmat) <- c("HSIC", "dCor", "MI")
    keep_cols <- apply(diffmat, 2, function(col) any(!is.na(col) & !is.nan(col)))
    diffmat_print <- diffmat[, keep_cols, drop = FALSE]
    keep_rows <- apply(diffmat_print, 1, function(row) any(!is.na(row) & !is.nan(row)))
    diffmat_print <- diffmat_print[keep_rows, , drop = FALSE]

    if (nrow(diffmat_print) > 0 && ncol(diffmat_print) > 0) {
      cat("\nDifference Statistics (mean estimates, lower, upper):\n")
      print(round(diffmat_print, digits))
      cat("---\n")
    }
  }

  # Decision proportions
  for (dname in names(decisions)) {
    prop <- decisions[[dname]]
    out <- rep(0, 3)
    names(out) <- c("undecided", "y->x", "x->y")
    out[names(prop)] <- prop

    if (any(!is.na(out) & out != 0)) {
      cat(paste("\nDecision proportions for", dname, ":\n"))
      print(round(out, digits))
      cat("---\n")
    }
  }

  invisible(object)
}

###############################################################################

#' Summary for dda_bagging Output (dda.resdist, classic style, robust)
#'
#' @param object Output from dda_bagging() for dda.resdist objects
#' @param digits Number of digits for rounding (default: 4)
#' @return Prints organized summary tables, matching classic print.dda.resdist style
#' @export

summary.dda_bagging_resdist <- function(object, digits = 4) {
  stats <- object$aggregated_stats
  params <- object$parameters %||% list()
  varnames <- if (!is.null(stats$var.names) && length(stats$var.names) == 2) stats$var.names else c("target", "alternative")

  # indep-style preamble
  cat("\n===== DDA Bagging Summary =====\n")
  cat("Function:", if (!is.null(params$function_name)) params$function_name else "unknown", "\n")
  cat("Object Type:", if (!is.null(params$object_type)) params$object_type else "unknown", "\n")
  cat("Iterations:", if (!is.null(params$iter)) params$iter else NA, "\n")
  cat("Completed Iterations:", if (!is.null(object$n_valid_iterations)) object$n_valid_iterations else NA, "\n")
  cat("----\n")

  # Classic header
  cat("\nDIRECTION DEPENDENCE ANALYSIS: Residual Distributions (Bagged)\n\n")
  cat("Skewness and kurtosis tests:\n")

  # Skewness & kurtosis test table
  skew_row <- c(
    if (!is.null(stats$agostino.target.statistic)) stats$agostino.target.statistic else NA,
    if (!is.null(stats$agostino.target.z)) stats$agostino.target.z else NA,
    if (!is.null(stats$agostino.target.p.value)) stats$agostino.target.p.value else NA,
    if (!is.null(stats$agostino.alternative.statistic)) stats$agostino.alternative.statistic else NA,
    if (!is.null(stats$agostino.alternative.z)) stats$agostino.alternative.z else NA,
    if (!is.null(stats$agostino.alternative.p.value)) stats$agostino.alternative.p.value else NA
  )
  kurt_row <- c(
    if (!is.null(stats$anscombe.target.statistic)) stats$anscombe.target.statistic else NA,
    if (!is.null(stats$anscombe.target.z)) stats$anscombe.target.z else NA,
    if (!is.null(stats$anscombe.target.p.value)) stats$anscombe.target.p.value else NA,
    if (!is.null(stats$anscombe.alternative.statistic)) stats$anscombe.alternative.statistic else NA,
    if (!is.null(stats$anscombe.alternative.z)) stats$anscombe.alternative.z else NA,
    if (!is.null(stats$anscombe.alternative.p.value)) stats$anscombe.alternative.p.value else NA
  )

  sigtests <- rbind(skew_row, kurt_row)
  sigtests <- round(sigtests, digits)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  # varnames[1] is the response (y) = target; varnames[2] is pred (x) = alternative
  colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")

  keep_rows <- apply(sigtests, 1, function(row) any(!is.na(row) & !is.nan(row)))
  sigtests_print <- sigtests[keep_rows, , drop = FALSE]
  if (nrow(sigtests_print) > 0) {
    print.default(format(sigtests_print, digits = max(3L, getOption("digits") - 3L), scientific = NA, scipen = 999),
                  print.gap = 2L, quote = FALSE)
  }

  # Skewness and kurtosis difference tests with CI
  if (!is.null(stats$skewdiff) || !is.null(stats$kurtdiff)) {
    cat("\n")
    # Attempt to detect CI level/type from any bagged result
    ci.level <- NULL; boot.type <- NULL; n.resamples <- NULL
    if (!is.null(object$bagged_results) && length(object$bagged_results) > 0) {
      for (res in object$bagged_results) {
        if (!is.null(res$boot.args)) {
          boot.type <- res$boot.args[1]
          ci.level <- as.numeric(res$boot.args[2]) * 100
          n.resamples <- as.numeric(res$boot.args[3])
          break
        }
      }
    }
    if (!is.null(ci.level) && !is.null(boot.type)) {
      if (boot.type == "bca")  cat("Skewness and kurtosis difference tests and ", ci.level, "% BCa bootstrap CIs:\n\n", sep = "")
      if (boot.type == "perc") cat("Skewness and kurtosis difference tests and ", ci.level, "% Percentile bootstrap CIs:\n\n", sep = "")
    } else {
      cat("Skewness and kurtosis difference tests and bootstrap CIs:\n\n")
    }

    skewdiff <- if (!is.null(stats$skewdiff)) stats$skewdiff else rep(NA, 5)
    kurtdiff <- if (!is.null(stats$kurtdiff)) stats$kurtdiff else rep(NA, 5)
    citests <- rbind(skewdiff, kurtdiff)
    rownames(citests) <- c("Skewness", "Kurtosis")
    colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
    citests <- round(citests, digits)
    keep_rows <- apply(citests, 1, function(row) any(!is.na(row) & !is.nan(row)))
    citests_print <- citests[keep_rows, , drop = FALSE]
    if (nrow(citests_print) > 0) {
      print.default(format(citests_print, digits = max(3L, getOption("digits") - 3L)),
                    print.gap = 2L, quote = FALSE)
    }
  }

  # Joint higher moment difference CIs
  joint_names <- c("Co-Skewness", "Hyvarinen-Smith (Co-Skewness)",
                   "Co-Kurtosis", "Hyvarinen-Smith (Co-Kurtosis)", "Chen-Chan (Co-Kurtosis)")
  joint_keys <- c("cor12diff", "RHS3", "cor13diff", "RHS4", "RCC")
  joint_vals <- lapply(joint_keys, function(k) if (!is.null(stats[[k]])) stats[[k]] else rep(NA, 3))
  joint_mat <- do.call(rbind, joint_vals)
  rownames(joint_mat) <- joint_names
  colnames(joint_mat) <- c("estimate", "lower", "upper")
  joint_mat <- round(joint_mat, digits)
  keep_rows <- apply(joint_mat, 1, function(row) any(!is.na(row) & !is.nan(row)))
  joint_print <- joint_mat[keep_rows, , drop = FALSE]
  if (nrow(joint_print) > 0) {
    cat("\n")
    ci.level <- NULL; boot.type <- NULL; n.resamples <- NULL
    if (!is.null(object$bagged_results) && length(object$bagged_results) > 0) {
      for (res in object$bagged_results) {
        if (!is.null(res$boot.args)) {
          boot.type <- res$boot.args[1]
          ci.level <- as.numeric(res$boot.args[2]) * 100
          n.resamples <- as.numeric(res$boot.args[3])
          break
        }
      }
    }
    if (!is.null(ci.level) && !is.null(boot.type)) {
      if (boot.type == "bca")  cat(ci.level, "% BCa bootstrap CIs for joint higher moment differences:\n", sep = "")
      if (boot.type == "perc") cat(ci.level, "% Percentile bootstrap CIs for joint higher moment differences:\n", sep = "")
    } else {
      cat("Bootstrap CIs for joint higher moment differences:\n")
    }
    print.default(format(joint_print, digits = max(3L, getOption("digits") - 3L)),
                  print.gap = 2L, quote = FALSE)
  }

  # Number of resamples (if available from any iteration)
  if (!is.null(n.resamples)) {
    cat("\nNumber of resamples:", n.resamples, "\n")
  }

  cat("---\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], "\n"))
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], "\n"))
  if (isFALSE(stats$probtrans)) {
    cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1], "\n"))
  } else {
    cat(paste("      Under prob.trans = TRUE, skewness and kurtosis differences < 0 and\n     co-skewness and co-kurtosis differences > 0 suggest", varnames[2], "->", varnames[1], "\n"))
  }

  invisible(object)
}

# Helper to provide a default value when NULL
`%||%` <- function(lhs, rhs) if (is.null(lhs)) rhs else lhs

################################################################################
#' Summary for dda_bagging Output (dda.vardist, classic style, robust)
#'
#' @param object Output from dda_bagging() for dda.vardist objects
#' @param digits Number of digits for rounding (default: 4)
#' @return Prints organized summary tables, matching classic print.dda.vardist style
#' @export
summary.dda_bagging_vardist <- function(object, digits = 4) {
  stats <- object$aggregated_stats
  params <- object$parameters
  # Outcome first, Predictor second (matches original print.dda.vardist)
  varnames <- if (!is.null(stats$var.names) && length(stats$var.names) == 2) stats$var.names else c("Outcome", "Predictor")

  # indep-style preamble
  cat("\n===== DDA Bagging Summary =====\n")
  cat("Function:", if (!is.null(params$function_name)) params$function_name else "unknown", "\n")
  cat("Object Type:", if (!is.null(params$object_type)) params$object_type else "unknown", "\n")
  cat("Iterations:", if (!is.null(params$iter)) params$iter else NA, "\n")
  cat("Completed Iterations:", if (!is.null(object$n_valid_iterations)) object$n_valid_iterations else NA, "\n")
  cat("----\n")

  cat("\nDIRECTION DEPENDENCE ANALYSIS: Variable Distributions (Bagged)\n\n")
  cat("Skewness and kurtosis tests:\n")

  # Skewness & kurtosis test table (robust to missing)
  skew_row <- c(
    if (!is.null(stats$agostino.outcome.statistic.skew)) stats$agostino.outcome.statistic.skew else NA,
    if (!is.null(stats$agostino.outcome.statistic.z))    stats$agostino.outcome.statistic.z    else NA,
    if (!is.null(stats$agostino.outcome.p.value))        stats$agostino.outcome.p.value        else NA,
    if (!is.null(stats$agostino.predictor.statistic.skew)) stats$agostino.predictor.statistic.skew else NA,
    if (!is.null(stats$agostino.predictor.statistic.z))    stats$agostino.predictor.statistic.z    else NA,
    if (!is.null(stats$agostino.predictor.p.value))        stats$agostino.predictor.p.value        else NA
  )
  kurt_row <- c(
    if (!is.null(stats$anscombe.outcome.statistic.kurt)) stats$anscombe.outcome.statistic.kurt else NA,
    if (!is.null(stats$anscombe.outcome.statistic.z))    stats$anscombe.outcome.statistic.z    else NA,
    if (!is.null(stats$anscombe.outcome.p.value))        stats$anscombe.outcome.p.value        else NA,
    if (!is.null(stats$anscombe.predictor.statistic.kurt)) stats$anscombe.predictor.statistic.kurt else NA,
    if (!is.null(stats$anscombe.predictor.statistic.z))    stats$anscombe.predictor.statistic.z    else NA,
    if (!is.null(stats$anscombe.predictor.p.value))        stats$anscombe.predictor.p.value        else NA
  )

  sigtests <- rbind(skew_row, kurt_row)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")
  sigtests <- round(sigtests, digits)
  keep_rows <- apply(sigtests, 1, function(row) any(!is.na(row) & !is.nan(row)))
  sigtests_print <- sigtests[keep_rows, , drop = FALSE]
  if (nrow(sigtests_print) > 0) {
    print.default(format(sigtests_print, digits = max(3L, getOption("digits") - 3L), scientific = NA, scipen = 999),
                  print.gap = 2L, quote = FALSE)
  }

  # Skewness/Kurtosis differences (diff, lower, upper)
  # Try to extract CI level & type from any bagged result
  ci.level <- NULL; boot.type <- NULL; n.resamples <- NULL
  if (!is.null(object$bagged_results) && length(object$bagged_results) > 0) {
    for (res in object$bagged_results) {
      if (!is.null(res$boot.args)) {
        boot.type <- res$boot.args[1]
        ci.level <- as.numeric(res$boot.args[2]) * 100
        n.resamples <- as.numeric(res$boot.args[3])
        break
      }
    }
  }

  # Differences heading
  if (!is.null(stats$skewdiff) || !is.null(stats$kurtdiff)) {
    cat("\n")
    if (!is.null(ci.level) && !is.null(boot.type)) {
      if (boot.type == "bca")  cat(ci.level, "% BCa bootstrap CIs for higher moment differences:\n", sep = "")
      if (boot.type == "perc") cat(ci.level, "% Percentile bootstrap CIs for higher moment differences:\n", sep = "")
    } else {
      cat("Bootstrap CIs for higher moment differences:\n")
    }
    skewdiff <- if (!is.null(stats$skewdiff)) stats$skewdiff else rep(NA_real_, 3)
    kurtdiff <- if (!is.null(stats$kurtdiff)) stats$kurtdiff else rep(NA_real_, 3)
    citests <- rbind(skewdiff, kurtdiff)
    rownames(citests) <- c("Skewness", "Kurtosis")
    colnames(citests) <- c("diff", "lower", "upper")
    citests <- round(citests, digits)
    keep_rows <- apply(citests, 1, function(row) any(!is.na(row) & !is.nan(row)))
    citests_print <- citests[keep_rows, , drop = FALSE]
    if (nrow(citests_print) > 0) {
      print.default(format(citests_print, digits = max(3L, getOption("digits") - 3L)),
                    print.gap = 2L, quote = FALSE)
    }
  }

  # Higher-order correlations differences (coskew/cokurt)
  hoctests <- rbind(
    if (!is.null(stats$cor12diff)) stats$cor12diff else rep(NA_real_, 3),
    if (!is.null(stats$cor13diff)) stats$cor13diff else rep(NA_real_, 3)
  )
  rownames(hoctests) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]")
  colnames(hoctests) <- c("estimate", "lower", "upper")
  hoctests <- round(hoctests, digits)
  keep_rows <- apply(hoctests, 1, function(row) any(!is.na(row) & !is.nan(row)))
  hoctests_print <- hoctests[keep_rows, , drop = FALSE]
  if (nrow(hoctests_print) > 0) {
    cat("\n")
    if (!is.null(ci.level) && !is.null(boot.type)) {
      if (boot.type == "bca")  cat(ci.level, "% BCa bootstrap CIs for differences in higher-order correlations:\n", sep = "")
      if (boot.type == "perc") cat(ci.level, "% Percentile bootstrap CIs for differences in higher-order correlations:\n", sep = "")
    } else {
      cat("Bootstrap CIs for differences in higher-order correlations:\n")
    }
    print.default(format(hoctests_print, digits = max(3L, getOption("digits") - 3L)),
                  print.gap = 2L, quote = FALSE)
  }

  # Likelihood ratio approximations
  LRtests <- rbind(
    if (!is.null(stats$RHS))   stats$RHS   else rep(NA_real_, 3),
    if (!is.null(stats$Rtanh)) stats$Rtanh else rep(NA_real_, 3),
    if (!is.null(stats$RCC))   stats$RCC   else rep(NA_real_, 3)
  )
  rownames(LRtests) <- c("Hyvarinen-Smith (co-skewness)", "Hyvarinen-Smith (tanh)", "Chen-Chan (co-kurtosis)")
  colnames(LRtests) <- c("estimate", "lower", "upper")
  LRtests <- round(LRtests, digits)
  keep_rows <- apply(LRtests, 1, function(row) any(!is.na(row) & !is.nan(row)))
  LRtests_print <- LRtests[keep_rows, , drop = FALSE]
  if (nrow(LRtests_print) > 0) {
    cat("\n")
    if (!is.null(ci.level) && !is.null(boot.type)) {
      if (boot.type == "bca")  cat(ci.level, "% BCa bootstrap CIs for Likelihood Ratio approximations:\n", sep = "")
      if (boot.type == "perc") cat(ci.level, "% Percentile bootstrap CIs for Likelihood Ratio approximations:\n", sep = "")
    } else {
      cat("Bootstrap CIs for Likelihood Ratio approximations:\n")
    }
    print.default(format(LRtests_print, digits = max(3L, getOption("digits") - 3L)),
                  print.gap = 2L, quote = FALSE)
  }

  # Number of resamples (if available from any iteration)
  if (!is.null(n.resamples)) {
    cat("\nNumber of resamples:", n.resamples, "\n")
  }

  cat("---\n")
  cat(paste("Note: (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests the model", varnames[2], "->", varnames[1], "\n"))

  invisible(object)
}
