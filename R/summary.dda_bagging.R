#' Summary for dda_bagging Output (dda.indep only)
#'
#' @param object Output from dda_bagging() for dda.indep objects
#' @param digits Number of digits for rounding (default: 4)
#' @param alpha Significance level for decisions (default: 0.05)
#' @return Prints organized summary tables
#' @export
#' @method summary dda_bagging_indep

summary.dda_bagging_indep <- function(object, digits = 4, alpha = 0.05) {
  stats <- object$aggregated_stats
  decisions <- object$decision_percentages

  cat("\n===== DDA Bagging Summary (INDEP) =====\n")
  cat("Function:", object$parameters$function_name, "\n")
  cat("Object Type:", object$parameters$object_type, "\n")
  cat("Iterations:", object$parameters$iter, "\n")
  cat("Valid Iterations:", object$n_valid_iterations, "\n")
  cat("----\n")

  # Only run for dda.indep
  if (object$parameters$object_type == "dda.indep" && !is.null(stats$hsic_yx_stat)) {
    stat_tab <- matrix(NA, nrow = 2, ncol = 4)
    rownames(stat_tab) <- c("HSIC", "dCor")
    colnames(stat_tab) <- c("Target Stat", "Target p", "Alternative Stat", "Alternative p")

    stat_tab[1, ] <- c(round(stats$hsic_yx_stat, digits),
                       round(stats$hsic_yx_pval, digits),
                       round(stats$hsic_xy_stat, digits),
                       round(stats$hsic_xy_pval, digits))
    stat_tab[2, ] <- c(round(stats$dcor_yx_stat, digits),
                       round(stats$dcor_yx_pval, digits),
                       round(stats$dcor_xy_stat, digits),
                       round(stats$dcor_xy_pval, digits))

    keep_rows <- apply(stat_tab, 1, function(row) any(!is.na(row) & !is.nan(row)))
    stat_tab_print <- stat_tab[keep_rows, , drop = FALSE]

    if (nrow(stat_tab_print) > 0) {
      cat("\nHSIC and dCor Test Statistics & Harmonic p-values:\n")
      print(stat_tab_print)
      cat("---\n")
    }

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

      # Print only if any category is nonzero/non-NA
      if (any(!is.na(out) & out != 0)) {
        cat(paste("\nDecision proportions for", dname, ":\n"))
        print(round(out, digits))
        cat("---\n")
      }
    }
  } else {
    cat("\nThis summary method only supports dda.indep objects.\n")
  }

  invisible(object)
}




#' Summary for dda_bagging Output (dda.vardist)
#'
#' @param object Output from dda_bagging() for dda.vardist objects
#' @param digits Number of digits for rounding (default: 4)
#' @return Prints organized summary tables, matching classic print.dda.vardist style
#' @export


summary.dda_bagging_vardist <- function(object, digits = 4) {
  stats <- object$aggregated_stats
  # Try to get variable names from stats or fallback
  varnames <- if (!is.null(stats$var.names) && length(stats$var.names) == 2) stats$var.names else c("Outcome", "Predictor")

  cat("\nDIRECTION DEPENDENCE ANALYSIS: Variable Distributions (Bagged)\n\n")
  cat("Skewness and kurtosis tests:\n")

  # --- Skewness & kurtosis test table ---
  skew_row <- c(stats$agostino.outcome.statistic.skew, stats$agostino.outcome.statistic.z, stats$agostino.outcome.p.value,
                stats$agostino.predictor.statistic.skew, stats$agostino.predictor.statistic.z, stats$agostino.predictor.p.value)
  kurt_row <- c(stats$anscombe.outcome.statistic.kurt, stats$anscombe.outcome.statistic.z, stats$anscombe.outcome.p.value,
                stats$anscombe.predictor.statistic.kurt, stats$anscombe.predictor.statistic.z, stats$anscombe.predictor.p.value)
  sigtests <- rbind(skew_row, kurt_row)
  sigtests <- round(sigtests, digits)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")
  # Only print if at least one non-NA value per row
  for (i in 1:nrow(sigtests)) {
    if (any(!is.na(sigtests[i, ]))) {
      cat("\n", rownames(sigtests)[i], ":\n", sep = "")
      print.default(format(sigtests[i, , drop = FALSE], digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    }
  }

  # --- Skewness/Kurtosis differences table (Bootstrap CIs) ---
  skewdiff <- c(stats$skewdiff.skew.diff, stats$skewdiff.lower, stats$skewdiff.upper)
  kurtdiff <- c(stats$kurtdiff.kurt.diff, stats$kurtdiff.lower, stats$kurtdiff.upper)
  citests <- rbind(skewdiff, kurtdiff)
  citests <- round(citests, digits)
  rownames(citests) <- c("Skewness", "Kurtosis")
  colnames(citests) <- c("diff", "lower", "upper")
  for (i in 1:nrow(citests)) {
    if (any(!is.na(citests[i, ]))) {
      cat("\n", rownames(citests)[i], " Difference (Bootstrap CI):\n", sep = "")
      print.default(format(citests[i, , drop = FALSE], digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    }
  }

  # --- Co-Skewness / Co-Kurtosis ---
  cor12diff <- c(stats$cor12diff.cor21.diff, stats$cor12diff.lower, stats$cor12diff.upper)
  cor13diff <- c(stats$cor13diff.cor13.diff, stats$cor13diff.lower, stats$cor13diff.upper)
  hoctests <- rbind(cor12diff, cor13diff)
  hoctests <- round(hoctests, digits)
  rownames(hoctests) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]")
  colnames(hoctests) <- c("estimate", "lower", "upper")
  for (i in 1:nrow(hoctests)) {
    if (any(!is.na(hoctests[i, ]))) {
      cat("\n", rownames(hoctests)[i], " (Bootstrap CI):\n", sep = "")
      print.default(format(hoctests[i, , drop = FALSE], digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    }
  }

  # --- LR approximations ---
  RHS <- c(stats$RHS.RHS, stats$RHS.lower, stats$RHS.upper)
  Rtanh <- c(stats$Rtanh.Rtanh, stats$Rtanh.lower, stats$Rtanh.upper)
  RCC <- c(stats$RCC.RCC, stats$RCC.lower, stats$RCC.upper)
  LRtests <- rbind(RHS, Rtanh, RCC)
  LRtests <- round(LRtests, digits)
  rownames(LRtests) <- c("Hyvarinen-Smith (co-skewness)", "Hyvarinen-Smith (tanh)", "Chen-Chan (co-kurtosis)")
  colnames(LRtests) <- c("estimate", "lower", "upper")
  for (i in 1:nrow(LRtests)) {
    if (any(!is.na(LRtests[i, ]))) {
      cat("\n", rownames(LRtests)[i], " (Bootstrap CI):\n", sep = "")
      print.default(format(LRtests[i, , drop = FALSE], digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    }
  }

  cat("\n---\n")
  cat(paste("Note: (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests the model", varnames[2], "->", varnames[1], "\n"))

  invisible(object)
}




#' Summary for dda_bagging Output (dda.resdist, classic style, no decisions)
#'
#' @param object Output from dda_bagging() for dda.resdist objects
#' @param digits Number of digits for rounding (default: 4)
#' @return Prints organized summary tables, matching classic print.dda.resdist style
#' @export

summary.dda_bagging_resdist <- function(object, digits = 4) {
  stats <- object$aggregated_stats
  varnames <- if (!is.null(stats$var.names) && length(stats$var.names) == 2) stats$var.names else c("target", "alternative")

  cat("\nDIRECTION DEPENDENCE ANALYSIS: Residual Distributions (Bagged)\n\n")
  cat("Skewness and kurtosis tests:\n")

  # Skewness & kurtosis test table (classic style)
  sigtests <- rbind(
    c(stats$agostino.target.statistic, stats$agostino.target.z, stats$agostino.target.p.value,
      stats$agostino.alternative.statistic, stats$agostino.alternative.z, stats$agostino.alternative.p.value),
    c(stats$anscombe.target.statistic, stats$anscombe.target.z, stats$anscombe.target.p.value,
      stats$anscombe.alternative.statistic, stats$anscombe.alternative.z, stats$anscombe.alternative.p.value)
  )
  sigtests <- round(sigtests, digits)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")
  print.default(format(sigtests, digits = max(3L, getOption("digits") - 3L), scientific = NA, scipen = 999), print.gap = 2L, quote = FALSE)

  # Difference tests (classic style: means across bagged runs)
  cat("\nSkewness and kurtosis difference tests:\n")
  citests <- rbind(stats$skewdiff, stats$kurtdiff)
  citests <- round(citests, digits)
  rownames(citests) <- c("Skewness", "Kurtosis")
  colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
  print.default(format(citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
  cat("\n")

  cat("Co-skewness and co-kurtosis difference tests:\n")
  jointtests <- rbind(stats$cor12diff, stats$RHS3, stats$cor13diff, stats$RHS4, stats$RCC)
  jointtests <- round(jointtests, digits)
  rownames(jointtests) <- c("Co-Skewness", "Hyvarinen-Smith (Co-Skewness)", "Co-Kurtosis", "Hyvarinen-Smith (Co-Kurtosis)", "Chen-Chan (Co-Kurtosis)")
  colnames(jointtests) <- c("estimate", "lower", "upper")
  print.default(format(jointtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

  cat("\n---\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], "\n"))
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], "\n"))
  cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1], "\n"))

  invisible(object)
}
