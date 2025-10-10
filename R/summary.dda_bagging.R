#' Summary for dda_bagging Output (dda.indep only)
#'
#' @param object Output from dda_bagging() for dda.indep objects
#' @param digits Number of digits for rounding (default: 4)
#' @param alpha Significance level for decisions (default: 0.05)
#' @return Prints organized summary tables
#' @export
#' @method summary dda_bagging

summary.dda_bagging <- function(object, digits = 4, alpha = 0.05) {
  stats <- object$aggregated_stats
  decisions <- object$decision_percentages

  cat("\n===== DDA Bagging Summary =====\n")
  cat("Function:", object$parameters$function_name, "\n")
  cat("Object Type:", object$parameters$object_type, "\n")
  cat("Iterations:", object$parameters$iter, "\n")
  cat("Valid Iterations:", object$n_valid_iterations, "\n")
  cat("----\n")

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

      if (any(!is.na(out) & out != 0)) {
        cat(paste("\nDecision proportions for", dname, ":\n"))
        print(round(out, digits))
        cat("---\n")
      }
    }
  }

  invisible(object)
}


#' Summary for dda_bagging Output (dda.vardist)
#'
#' @param object Output from dda_bagging() for dda.vardist objects
#' @param digits Number of digits for rounding (default: 4)
#' @param alpha Significance level for decisions (default: 0.05)
#' @return Prints organized summary tables
#' @export

summary.dda_bagging_vardist <- function(object, digits = 4, alpha = 0.05) {
  stats <- object$aggregated_stats
  decisions <- object$decision_percentages

  cat("\n===== DDA Bagging Summary: Variable Distributions =====\n")
  cat("Function:", object$parameters$function_name, "\n")
  cat("Object Type:", object$parameters$object_type, "\n")
  cat("Iterations:", object$parameters$iter, "\n")
  cat("Valid Iterations:", object$n_valid_iterations, "\n")
  cat("----\n")

  # Skewness and Kurtosis Tests
  if (!is.null(stats$agostino.outcome.statistic.skew)) {
    test_tab <- matrix(NA, nrow = 2, ncol = 6)
    rownames(test_tab) <- c("Skewness", "Kurtosis")
    colnames(test_tab) <- c("Outcome", "z-value", "p-value", "Predictor", "z-value", "p-value")

    test_tab[1, ] <- c(round(stats$agostino.outcome.statistic.skew, digits),
                       round(stats$agostino.outcome.statistic.z, digits),
                       round(stats$agostino.outcome.p.value, digits),
                       round(stats$agostino.predictor.statistic.skew, digits),
                       round(stats$agostino.predictor.statistic.z, digits),
                       round(stats$agostino.predictor.p.value, digits))

    test_tab[2, ] <- c(round(stats$anscombe.outcome.statistic.kurt, digits),
                       round(stats$anscombe.outcome.statistic.z, digits),
                       round(stats$anscombe.outcome.p.value, digits),
                       round(stats$anscombe.predictor.statistic.kurt, digits),
                       round(stats$anscombe.predictor.statistic.z, digits),
                       round(stats$anscombe.predictor.p.value, digits))

    keep_rows <- apply(test_tab, 1, function(row) any(!is.na(row) & !is.nan(row)))
    test_tab_print <- test_tab[keep_rows, , drop = FALSE]

    if (nrow(test_tab_print) > 0) {
      cat("\nSkewness and Kurtosis Tests:\n")
      print(test_tab_print)
      cat("---\n")
    }
  }

  # Higher Moment Differences
  if (!is.null(stats$skewdiff.skew.diff)) {
    moment_tab <- matrix(NA, nrow = 2, ncol = 3)
    rownames(moment_tab) <- c("Skewness", "Kurtosis")
    colnames(moment_tab) <- c("diff", "lower", "upper")

    moment_tab[1, ] <- c(round(stats$skewdiff.skew.diff, digits),
                         round(stats$skewdiff.lower, digits),
                         round(stats$skewdiff.upper, digits))

    moment_tab[2, ] <- c(round(stats$kurtdiff.kurt.diff, digits),
                         round(stats$kurtdiff.lower, digits),
                         round(stats$kurtdiff.upper, digits))

    keep_rows <- apply(moment_tab, 1, function(row) any(!is.na(row) & !is.nan(row)))
    moment_tab_print <- moment_tab[keep_rows, , drop = FALSE]

    if (nrow(moment_tab_print) > 0) {
      cat("\nHigher Moment Differences (Bootstrap CIs):\n")
      print(moment_tab_print)
      cat("---\n")
    }
  }

  # Higher-Order Correlations
  if (!is.null(stats$cor12diff.cor21.diff)) {
    hoc_tab <- matrix(NA, nrow = 2, ncol = 3)
    rownames(hoc_tab) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]")
    colnames(hoc_tab) <- c("estimate", "lower", "upper")

    hoc_tab[1, ] <- c(round(stats$cor12diff.cor21.diff, digits),
                      round(stats$cor12diff.lower, digits),
                      round(stats$cor12diff.upper, digits))

    hoc_tab[2, ] <- c(round(stats$cor13diff.cor13.diff, digits),
                      round(stats$cor13diff.lower, digits),
                      round(stats$cor13diff.upper, digits))

    keep_rows <- apply(hoc_tab, 1, function(row) any(!is.na(row) & !is.nan(row)))
    hoc_tab_print <- hoc_tab[keep_rows, , drop = FALSE]

    if (nrow(hoc_tab_print) > 0) {
      cat("\nHigher-Order Correlations:\n")
      print(hoc_tab_print)
      cat("---\n")
    }
  }

  # Likelihood Ratio Approximations
  if (!is.null(stats$RHS.RHS)) {
    lr_tab <- matrix(NA, nrow = 3, ncol = 3)
    rownames(lr_tab) <- c("Hyvarinen-Smith (co-skewness)",
                          "Hyvarinen-Smith (tanh)",
                          "Chen-Chan (co-kurtosis)")
    colnames(lr_tab) <- c("estimate", "lower", "upper")

    lr_tab[1, ] <- c(round(stats$RHS.RHS, digits),
                     round(stats$RHS.lower, digits),
                     round(stats$RHS.upper, digits))

    lr_tab[2, ] <- c(round(stats$Rtanh.Rtanh, digits),
                     round(stats$Rtanh.lower, digits),
                     round(stats$Rtanh.upper, digits))

    lr_tab[3, ] <- c(round(stats$RCC.RCC, digits),
                     round(stats$RCC.lower, digits),
                     round(stats$RCC.upper, digits))

    keep_rows <- apply(lr_tab, 1, function(row) any(!is.na(row) & !is.nan(row)))
    lr_tab_print <- lr_tab[keep_rows, , drop = FALSE]

    if (nrow(lr_tab_print) > 0) {
      cat("\nLikelihood Ratio Approximations:\n")
      print(lr_tab_print)
      cat("---\n")
    }
  }

  # Decision Percentages
  if (!is.null(decisions) && length(decisions) > 0) {
    cat("\nDecision Percentages (alpha =", alpha, "):\n")
    for (dec_name in names(decisions)) {
      prop <- decisions[[dec_name]]
      out <- rep(0, 3)
      names(out) <- c("undecided", "y->x", "x->y")
      out[names(prop)] <- prop

      if (any(!is.na(out) & out != 0)) {
        cat(dec_name, ":\n")
        print(round(out, digits))
        cat("---\n")
      }
    }
  }

  cat("\nNote: (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests predictor -> outcome\n")

  invisible(object)
}


#' Summary for dda_bagging Output (dda.resdist)
#'
#' @param object Output from dda_bagging() for dda.resdist objects
#' @param digits Number of digits for rounding (default: 4)
#' @param alpha Significance level for decisions (default: 0.05)
#' @return Prints organized summary tables
#' @export

summary.dda_bagging_resdist <- function(object, digits = 4, alpha = 0.05) {
  stats <- object$aggregated_stats
  decisions <- object$decision_percentages

  cat("\n===== DDA Bagging Summary: Residual Distributions =====\n")
  cat("Function:", object$parameters$function_name, "\n")
  cat("Object Type:", object$parameters$object_type, "\n")
  cat("Iterations:", object$parameters$iter, "\n")
  cat("Valid Iterations:", object$n_valid_iterations, "\n")
  cat("----\n")

  # Skewness and Kurtosis Tests
  if (!is.null(stats$agostino.target.statistic.skew)) {
    test_tab <- matrix(NA, nrow = 2, ncol = 6)
    rownames(test_tab) <- c("Skewness", "Kurtosis")
    colnames(test_tab) <- c("Target", "z-value", "p-value", "Alternative", "z-value", "p-value")

    test_tab[1, ] <- c(round(stats$agostino.target.statistic.skew, digits),
                       round(stats$agostino.target.statistic.z, digits),
                       round(stats$agostino.target.p.value, digits),
                       round(stats$agostino.alternative.statistic.skew, digits),
                       round(stats$agostino.alternative.statistic.z, digits),
                       round(stats$agostino.alternative.p.value, digits))

    test_tab[2, ] <- c(round(stats$anscombe.target.statistic.kurt, digits),
                       round(stats$anscombe.target.statistic.z, digits),
                       round(stats$anscombe.target.p.value, digits),
                       round(stats$anscombe.alternative.statistic.kurt, digits),
                       round(stats$anscombe.alternative.statistic.z, digits),
                       round(stats$anscombe.alternative.p.value, digits))

    keep_rows <- apply(test_tab, 1, function(row) any(!is.na(row) & !is.nan(row)))
    test_tab_print <- test_tab[keep_rows, , drop = FALSE]

    if (nrow(test_tab_print) > 0) {
      cat("\nSkewness and Kurtosis Tests:\n")
      print(test_tab_print)
      cat("---\n")
    }
  }

  # Skewness and Kurtosis Difference Tests
  if (!is.null(stats$skewdiff1)) {
    diff_tab <- matrix(NA, nrow = 2, ncol = 5)
    rownames(diff_tab) <- c("Skewness", "Kurtosis")
    colnames(diff_tab) <- c("diff", "z-value", "p-value", "lower", "upper")

    diff_tab[1, ] <- c(round(stats$skewdiff1, digits),
                       round(stats$skewdiff.z.value, digits),
                       round(stats$skewdiff.p.value, digits),
                       round(stats$skewdiff.lower, digits),
                       round(stats$skewdiff.upper, digits))

    diff_tab[2, ] <- c(round(stats$kurtdiff1, digits),
                       round(stats$kurtdiff.z.value, digits),
                       round(stats$kurtdiff.p.value, digits),
                       round(stats$kurtdiff.lower, digits),
                       round(stats$kurtdiff.upper, digits))

    keep_rows <- apply(diff_tab, 1, function(row) any(!is.na(row) & !is.nan(row)))
    diff_tab_print <- diff_tab[keep_rows, , drop = FALSE]

    if (nrow(diff_tab_print) > 0) {
      cat("\nSkewness and Kurtosis Difference Tests:\n")
      print(diff_tab_print)
      cat("---\n")
    }
  }

  # Joint Higher Moment Differences
  if (!is.null(stats$cor12diff.cor21.diff)) {
    joint_tab <- matrix(NA, nrow = 5, ncol = 3)
    rownames(joint_tab) <- c("Co-Skewness",
                             "Hyvarinen-Smith (Co-Skewness)",
                             "Co-Kurtosis",
                             "Hyvarinen-Smith (Co-Kurtosis)",
                             "Chen-Chan (Co-Kurtosis)")
    colnames(joint_tab) <- c("estimate", "lower", "upper")

    joint_tab[1, ] <- c(round(stats$cor12diff.cor21.diff, digits),
                        round(stats$cor12diff.lower, digits),
                        round(stats$cor12diff.upper, digits))

    joint_tab[2, ] <- c(round(stats$RHS3.RHS3, digits),
                        round(stats$RHS3.lower, digits),
                        round(stats$RHS3.upper, digits))

    joint_tab[3, ] <- c(round(stats$cor13diff.cor13.diff, digits),
                        round(stats$cor13diff.lower, digits),
                        round(stats$cor13diff.upper, digits))

    joint_tab[4, ] <- c(round(stats$RHS4.RHS4, digits),
                        round(stats$RHS4.lower, digits),
                        round(stats$RHS4.upper, digits))

    joint_tab[5, ] <- c(round(stats$RCC.RCC, digits),
                        round(stats$RCC.lower, digits),
                        round(stats$RCC.upper, digits))

    keep_rows <- apply(joint_tab, 1, function(row) any(!is.na(row) & !is.nan(row)))
    joint_tab_print <- joint_tab[keep_rows, , drop = FALSE]

    if (nrow(joint_tab_print) > 0) {
      cat("\nJoint Higher Moment Differences:\n")
      print(joint_tab_print)
      cat("---\n")
    }
  }

  # Decision Percentages
  if (!is.null(decisions) && length(decisions) > 0) {
    cat("\nDecision Percentages (alpha =", alpha, "):\n")
    for (dec_name in names(decisions)) {
      prop <- decisions[[dec_name]]
      out <- rep(0, 3)
      names(out) <- c("undecided", "y->x", "x->y")
      out[names(prop)] <- prop

      if (any(!is.na(out) & out != 0)) {
        cat(dec_name, ":\n")
        print(round(out, digits))
        cat("---\n")
      }
    }
  }

  cat("\nNote: Target = predictor -> outcome; Alternative = outcome -> predictor\n")
  cat("      Difference statistics > 0 suggest the target model\n")

  invisible(object)
}
