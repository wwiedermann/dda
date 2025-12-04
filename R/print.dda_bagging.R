#' Print for dda_bagging Output (INDEP)
#'
#' @param object Output from dda_bagging() for dda. indep objects (class: dda_bagging_indep)
#' @param digits Number of digits for rounding (default: 4)
#' @param alpha Significance level for decisions (default: 0.05)
#' @export
#' @method print dda_bagging_indep
print.dda_bagging_indep <- function(object, digits = 4, alpha = 0.05) {
  stats <- object$aggregated_stats

  # Extract variable names - try multiple locations where they might be stored
  # First try aggregated_stats$var.names, then try from first bagged result
  varnames <- stats$var.names
  if (is.null(varnames) || length(varnames) != 2) {
    # Try to get from the first valid bagged result
    if (! is.null(object$bagged_results) && length(object$bagged_results) > 0) {
      varnames <- object$bagged_results[[1]]$var.names
    }
  }
  # Final fallback
  if (is.null(varnames) || length(varnames) != 2) {
    varnames <- c("y", "x")
  }

  cat("\n")
  cat("DIRECTION DEPENDENCE ANALYSIS: Independence Properties (Bagged)", "\n")
  cat("Number of Bootstrap Aggregated Results:", object$n_valid_iterations, "\n", "\n")

  # -------------------- Target Model (yx) --------------------
  cat(paste("Target Model:", varnames[2], "->", varnames[1], sep = " "), "\n", "\n")

  # Omnibus Tests
  if (! is.null(stats$hsic_yx_stat) || !is.null(stats$dcor_yx_stat)) {
    cat("Omnibus Independence Tests:", "\n")
    if (!is.null(stats$hsic_yx_stat)) {
      cat(paste0("HSIC = ", round(stats$hsic_yx_stat, digits),
                 ", p-value = ", round(stats$hsic_yx_pval, digits)), "\n")
    }
    if (! is.null(stats$dcor_yx_stat)) {
      cat(paste0("dCor = ", round(stats$dcor_yx_stat, digits),
                 ", p-value = ", round(stats$dcor_yx_pval, digits)), "\n")
    }
    cat("\n")
  }

  # Homoscedasticity Tests (Target)
  if (! is.null(stats$breusch_pagan) && length(stats$breusch_pagan) >= 2) {
    cat("Homoscedasticity Tests:", "\n")
    bp_tab <- rbind(
      unlist(stats$breusch_pagan[[1]]),
      unlist(stats$breusch_pagan[[2]])
    )
    rownames(bp_tab) <- c("BP-test", "Robust BP-test")
    if(ncol(bp_tab) >= 3) {
      bp_tab <- bp_tab[, c("statistic", "parameter", "p. value"), drop=FALSE]
      colnames(bp_tab) <- c("X-squared", "df", "p-value")
      print(round(bp_tab, digits))
    }
    cat("\n")
  }

  # Non-linear Correlation (Target)
  if (!is.null(stats$nlcor.yx)) {
    nl <- stats$nlcor.yx
    func <- nl$func
    if (is.na(suppressWarnings(as.numeric(func)))) {
      cat(paste("Non-linear Correlation Tests:", func, "Transformation"), "\n")
      labs <- c(
        paste0("Cor[", func, "(", "r_", varnames[1], "), ", varnames[2], "]"),
        paste0("Cor[", "r_", varnames[1], ", ", func, "(", varnames[2], ")]"),
        paste0("Cor[", func, "(", "r_", varnames[1], "), ", func, "(", varnames[2], ")]")
      )
    } else {
      cat(paste("Non-linear Correlation Tests: Power Transformation using", func), "\n")
      labs <- c(
        paste0("Cor[r_", varnames[1], "^", func, ", ", varnames[2], "]"),
        paste0("Cor[r_", varnames[1], ", ", varnames[2], "^", func, "]"),
        paste0("Cor[r_", varnames[1], "^", func, ", ", varnames[2], "^", func, "]")
      )
    }

    tab <- rbind(nl$t1, nl$t2, nl$t3)
    colnames(tab) <- c("estimate", "t-value", "df", "Pr(>|t|)")
    rownames(tab) <- labs
    print(round(tab, digits))
    cat("\n")
  }

  # -------------------- Alternative Model (xy) --------------------
  cat(paste("Alternative Model:", varnames[1], "->", varnames[2], sep = " "), "\n", "\n")

  # Omnibus Tests
  if (!is.null(stats$hsic_xy_stat) || !is.null(stats$dcor_xy_stat)) {
    cat("Omnibus Independence Tests:", "\n")
    if (!is.null(stats$hsic_xy_stat)) {
      cat(paste0("HSIC = ", round(stats$hsic_xy_stat, digits),
                 ", p-value = ", round(stats$hsic_xy_pval, digits)), "\n")
    }
    if (! is.null(stats$dcor_xy_stat)) {
      cat(paste0("dCor = ", round(stats$dcor_xy_stat, digits),
                 ", p-value = ", round(stats$dcor_xy_pval, digits)), "\n")
    }
    cat("\n")
  }

  # Homoscedasticity Tests (Alternative)
  if (!is.null(stats$breusch_pagan) && length(stats$breusch_pagan) >= 4) {
    cat("Homoscedasticity Tests:", "\n")
    bp_tab <- rbind(
      unlist(stats$breusch_pagan[[3]]),
      unlist(stats$breusch_pagan[[4]])
    )
    rownames(bp_tab) <- c("BP-test", "Robust BP-test")
    if(ncol(bp_tab) >= 3) {
      bp_tab <- bp_tab[, c("statistic", "parameter", "p.value"), drop=FALSE]
      colnames(bp_tab) <- c("X-squared", "df", "p-value")
      print(round(bp_tab, digits))
    }
    cat("\n")
  }

  # Non-linear Correlation (Alternative)
  if (!is.null(stats$nlcor.xy)) {
    nl <- stats$nlcor.xy
    func <- nl$func
    if (is.na(suppressWarnings(as.numeric(func)))) {
      cat(paste("Non-linear Correlation Tests:", func, "Transformation"), "\n")
      labs <- c(
        paste0("Cor[", func, "(", "r_", varnames[2], "), ", varnames[1], "]"),
        paste0("Cor[", "r_", varnames[2], ", ", func, "(", varnames[1], ")]"),
        paste0("Cor[", func, "(", "r_", varnames[2], "), ", func, "(", varnames[1], ")]")
      )
    } else {
      cat(paste("Non-linear Correlation Tests: Power Transformation using", func), "\n")
      labs <- c(
        paste0("Cor[r_", varnames[2], "^", func, ", ", varnames[1], "]"),
        paste0("Cor[r_", varnames[2], ", ", varnames[1], "^", func, "]"),
        paste0("Cor[r_", varnames[2], "^", func, ", ", varnames[1], "^", func, "]")
      )
    }

    tab <- rbind(nl$t1, nl$t2, nl$t3)
    colnames(tab) <- c("estimate", "t-value", "df", "Pr(>|t|)")
    rownames(tab) <- labs
    print(round(tab, digits))
    cat("\n")
  }

  # Difference Statistics CIs
  if (! is.null(stats$diff_matrix)) {
    boot_type <- "Percentile"
    if (! is.null(object$bagged_results[[1]]$boot.args)) {
      type_arg <- object$bagged_results[[1]]$boot.args[1]
      if (! is.null(type_arg)) {
        if(type_arg == "bca") boot_type <- "BCa"
        if(type_arg == "perc") boot_type <- "Percentile"
      }
    }
    conf_lev <- "95%"
    if (!is.null(object$bagged_results[[1]]$boot.args)) {
      lev <- as.numeric(object$bagged_results[[1]]$boot.args[2])
      if(! is.na(lev)) conf_lev <- paste0(lev*100, "%")
    }

    cat(paste0(conf_lev, " ", boot_type, " Bootstrap CIs for Difference Statistics (", object$n_valid_iterations, " samples):"), "\n")
    print(round(stats$diff_matrix, digits))
    cat("\n")
  }

  cat("---\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], sep = " "), "\n")
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], sep = " "), "\n")
  cat(paste("      Difference statistics > 0 suggest", varnames[2], "->", varnames[1], sep = " "), "\n")

  invisible(object)
}


#' Print for dda_bagging Output (RESDIST)
#' @param object Output from dda_bagging() for dda.resdist objects (class: dda_bagging_resdist)
#' @param digits Number of digits for rounding (default: 4)
#' @export
#' @method print dda_bagging_resdist
print.dda_bagging_resdist <- function(object, digits = 4) {
  stats <- object$aggregated_stats

  # Extract variable names - try multiple locations where they might be stored
  varnames <- stats$var.names
  if (is.null(varnames) || length(varnames) != 2) {
    if (! is.null(object$bagged_results) && length(object$bagged_results) > 0) {
      varnames <- object$bagged_results[[1]]$var.names
    }
  }
  if (is.null(varnames) || length(varnames) != 2) {
    varnames <- c("y", "x")
  }

  cat("\n")
  cat("DIRECTION DEPENDENCE ANALYSIS: Residual Distributions (Bagged)", "\n")
  cat("Number of Bootstrap Aggregated Results:", object$n_valid_iterations, "\n", "\n")

  cat("Skewness and kurtosis tests:", "\n")

  skew_row <- c(stats$agostino.target.statistic, stats$agostino.target.z, stats$agostino.target.p.value,
                stats$agostino.alternative.statistic, stats$agostino.alternative.z, stats$agostino.alternative.p.value)
  kurt_row <- c(stats$anscombe.target.statistic, stats$anscombe.target.z, stats$anscombe.target.p.value,
                stats$anscombe.alternative.statistic, stats$anscombe.alternative.z, stats$anscombe.alternative.p.value)

  sigtests <- rbind(skew_row, kurt_row)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c("target", "z-value", "Pr(>|z|)", "alternative", "z-value", "Pr(>|z|)")
  print(round(sigtests, digits))

  cat("\n")

  # Determine bootstrap info
  boot_type <- "Percentile"
  conf_lev <- "95%"
  if (!is.null(object$bagged_results[[1]]$boot.args)) {
    type_arg <- object$bagged_results[[1]]$boot.args[1]
    if (!is.null(type_arg)) {
      if(type_arg == "bca") boot_type <- "BCa"
      if(type_arg == "perc") boot_type <- "Percentile"
    }
    lev <- as.numeric(object$bagged_results[[1]]$boot.args[2])
    if(!is.na(lev)) conf_lev <- paste0(lev*100, "%")
  }

  cat(paste0("Skewness and kurtosis difference tests and ", conf_lev, " ", boot_type, " bootstrap CIs:", "\n", "\n"))

  citests <- rbind(stats$skewdiff, stats$kurtdiff)
  rownames(citests) <- c("Skewness", "Kurtosis")
  if (ncol(citests) == 5) colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
  if (ncol(citests) == 3) colnames(citests) <- c("diff", "lower", "upper")
  print(round(citests, digits))

  cat("\n")
  cat(paste0(conf_lev, " ", boot_type, " bootstrap CIs for joint higher moment differences:", "\n"))

  joint_mat <- rbind(stats$cor12diff, stats$RHS3, stats$cor13diff, stats$RHS4, stats$RCC)
  rownames(joint_mat) <- c("Co-Skewness", "Hyvarinen-Smith (Co-Skewness)", "Co-Kurtosis", "Hyvarinen-Smith (Co-Kurtosis)", "Chen-Chan (Co-Kurtosis)")
  colnames(joint_mat) <- c("estimate", "lower", "upper")
  print(round(joint_mat, digits))

  cat("\n")
  cat(paste("Number of resamples:", object$n_valid_iterations), "\n")
  cat("---\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], sep = " "), "\n")
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], sep = " "), "\n")
  cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1], sep = " "), "\n")

  invisible(object)
}


#' Print for dda_bagging Output (VARDIST)
#' @param object Output from dda_bagging() for dda. vardist objects (class: dda_bagging_vardist)
#' @param digits Number of digits for rounding (default: 4)
#' @export
#' @method print dda_bagging_vardist
print.dda_bagging_vardist <- function(object, digits = 4) {
  stats <- object$aggregated_stats

  # Extract variable names - try multiple locations where they might be stored
  varnames <- stats$var.names
  if (is.null(varnames) || length(varnames) != 2) {
    if (! is.null(object$bagged_results) && length(object$bagged_results) > 0) {
      varnames <- object$bagged_results[[1]]$var.names
    }
  }
  if (is.null(varnames) || length(varnames) != 2) {
    varnames <- c("y", "x")
  }

  cat("\n")
  cat("DIRECTION DEPENDENCE ANALYSIS: Variable Distributions (Bagged)", "\n")
  cat("Number of Bootstrap Aggregated Results:", object$n_valid_iterations, "\n", "\n")

  cat("Skewness and kurtosis tests:", "\n")

  skew_row <- c(stats$agostino.outcome.statistic.skew, stats$agostino.outcome.statistic.z, stats$agostino.outcome.p.value,
                stats$agostino.predictor.statistic.skew, stats$agostino.predictor.statistic.z, stats$agostino.predictor.p.value)
  kurt_row <- c(stats$anscombe.outcome.statistic.kurt, stats$anscombe.outcome.statistic.z, stats$anscombe.outcome.p.value,
                stats$anscombe.predictor.statistic.kurt, stats$anscombe.predictor.statistic.z, stats$anscombe.predictor.p.value)

  sigtests <- rbind(skew_row, kurt_row)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")
  print(round(sigtests, digits))

  cat("\n")

  # Determine bootstrap info
  boot_type <- "Percentile"
  conf_lev <- "95%"
  if (!is.null(object$bagged_results[[1]]$boot.args)) {
    type_arg <- object$bagged_results[[1]]$boot.args[1]
    if (!is.null(type_arg)) {
      if(type_arg == "bca") boot_type <- "BCa"
      if(type_arg == "perc") boot_type <- "Percentile"
    }
    lev <- as.numeric(object$bagged_results[[1]]$boot.args[2])
    if(! is.na(lev)) conf_lev <- paste0(lev*100, "%")
  }

  cat(paste0(conf_lev, " ", boot_type, " bootstrap CIs for higher moment differences:", "\n"))
  citests <- rbind(stats$skewdiff, stats$kurtdiff)
  rownames(citests) <- c("Skewness", "Kurtosis")
  colnames(citests) <- c("diff", "lower", "upper")
  print(round(citests, digits))

  cat("\n")
  cat(paste0(conf_lev, " ", boot_type, " bootstrap CIs for differences in higher-order correlations:", "\n"))
  hoctests <- rbind(stats$cor12diff, stats$cor13diff)
  rownames(hoctests) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]")
  colnames(hoctests) <- c("estimate", "lower", "upper")
  print(round(hoctests, digits))

  cat("\n")
  cat(paste0(conf_lev, " ", boot_type, " bootstrap CIs for Likelihood Ratio approximations:", "\n"))
  LRtests <- rbind(stats$RHS, stats$Rtanh, stats$RCC)
  rownames(LRtests) <- c("Hyvarinen-Smith (co-skewness)", "Hyvarinen-Smith (tanh)", "Chen-Chan (co-kurtosis)")
  colnames(LRtests) <- c("estimate", "lower", "upper")
  print(round(LRtests, digits))

  cat("\n")
  cat(paste("Number of resamples:", object$n_valid_iterations), "\n")
  cat("---\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], sep = " "), "\n")
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], sep = " "), "\n")
  cat(paste("      (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests the model", varnames[2], "->", varnames[1], sep = " "), "\n")

  invisible(object)
}
