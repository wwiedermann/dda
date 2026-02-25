#' Print for dda_bagging Output (INDEP)
#'
#' @param object Output from dda_bagging() for dda.indep objects (class: dda_bagging_indep)
#' @param digits Number of digits for rounding (default: 4)
#' @param alpha Significance level for decisions (default: 0.05)
#' @export
#' @method print dda_bagging_indep
print.dda_bagging_indep <- function(object, digits = 4, alpha = 0.05) {
  stats <- object$aggregated_stats

  # Try to find variable names (default y, x if missing)
  varnames <- NULL
  if (!is.null(stats$var.names) && length(stats$var.names) == 2) {
    varnames <- stats$var.names
  } else {
    varnames <- c("y", "x")
  }

  cat("\nBOOTSTRAP AGGREGATED DDA: Independence Properties\n")
  cat("Number of Bootstrap Samples:", object$n_valid_iterations, "\n\n")

  # Helper to check if a value should be printed (not null, not na, not nan)
  should_print <- function(val) {
    !is.null(val) && !is.na(val) && !is.nan(val)
  }

  # -------------------- Target Model (yx) --------------------
  cat(paste("Target Model:", varnames[2], "->", varnames[1], sep = " "), "\n\n")

  # Omnibus Tests
  print_hsic <- should_print(stats$hsic_yx_stat)
  print_dcor <- should_print(stats$dcor_yx_stat)

  if (print_hsic || print_dcor) {
    cat("Omnibus Independence Tests:\n")
    if (print_hsic) {
      cat(paste0("HSIC = ", signif(as.numeric(stats$hsic_yx_stat), digits),
                 ", p-value = ", signif(as.numeric(stats$hsic_yx_pval), digits)), "\n")
    }
    if (print_dcor) {
      cat(paste0("dCor = ", signif(as.numeric(stats$dcor_yx_stat), digits),
                 ", p-value = ", signif(as.numeric(stats$dcor_yx_pval), digits)), "\n")
    }
    cat("\n")
  }

  # Homoscedasticity Tests (Target)
  if (!is.null(stats$breusch_pagan) && length(stats$breusch_pagan) >= 2) {
    cat("Homoscedasticity Tests:\n")
    cat(paste0("Standard Breusch-Pagan test: BP = ", signif(as.numeric(stats$breusch_pagan[[1]]$statistic), digits),
               ", df = ", as.numeric(stats$breusch_pagan[[1]]$parameter),
               ", p-value = ", signif(as.numeric(stats$breusch_pagan[[1]]$p.value), digits)), "\n")
    cat(paste0("Robust Breusch-Pagan test:   BP = ", signif(as.numeric(stats$breusch_pagan[[2]]$statistic), digits),
               ", df = ", as.numeric(stats$breusch_pagan[[2]]$parameter),
               ", p-value = ", signif(as.numeric(stats$breusch_pagan[[2]]$p.value), digits)), "\n")
    cat("\n")
  }

  # Non-linear Correlation (Target)
  if (!is.null(stats$nlcor.yx)) {
    nl <- stats$nlcor.yx
    func <- nl$func
    # Construct labels based on variable names and function
    if (is.na(suppressWarnings(as.numeric(func)))) {
      cat(paste("Non-linear Correlation Tests:", func, "Transformation\n"))
      labs <- c(
        paste0("Cor[", func, "(", "r_", varnames[1], "), ", varnames[2], "]"),
        paste0("Cor[", "r_", varnames[1], ", ", func, "(", varnames[2], ")]"),
        paste0("Cor[", func, "(", "r_", varnames[1], "), ", func, "(", varnames[2], ")]")
      )
    } else {
      cat(paste("Non-linear Correlation Tests: Power Transformation using", func, "\n"))
      labs <- c(
        paste0("Cor[r_", varnames[1], "^", func, ", ", varnames[2], "]"),
        paste0("Cor[r_", varnames[1], ", ", varnames[2], "^", func, "]"),
        paste0("Cor[r_", varnames[1], "^", func, ", ", varnames[2], "^", func, "]")
      )
    }

    tab <- rbind(as.numeric(nl$t1), as.numeric(nl$t2), as.numeric(nl$t3))
    colnames(tab) <- c("estimate", "t-value", "df", "Pr(>|t|)")
    rownames(tab) <- labs
    print.default(format(tab, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  # -------------------- Alternative Model (xy) --------------------
  cat(paste("Alternative Model:", varnames[1], "->", varnames[2], sep = " "), "\n\n")

  # Omnibus Tests
  print_hsic_alt <- should_print(stats$hsic_xy_stat)
  print_dcor_alt <- should_print(stats$dcor_xy_stat)

  if (print_hsic_alt || print_dcor_alt) {
    cat("Omnibus Independence Tests:\n")
    if (print_hsic_alt) {
      cat(paste0("HSIC = ", signif(as.numeric(stats$hsic_xy_stat), digits),
                 ", p-value = ", signif(as.numeric(stats$hsic_xy_pval), digits)), "\n")
    }
    if (print_dcor_alt) {
      cat(paste0("dCor = ", signif(as.numeric(stats$dcor_xy_stat), digits),
                 ", p-value = ", signif(as.numeric(stats$dcor_xy_pval), digits)), "\n")
    }
    cat("\n")
  }

  # Homoscedasticity Tests (Alternative)
  if (!is.null(stats$breusch_pagan) && length(stats$breusch_pagan) >= 4) {
    cat("Homoscedasticity Tests:\n")
    cat(paste0("Standard Breusch-Pagan test: BP = ", signif(as.numeric(stats$breusch_pagan[[3]]$statistic), digits),
               ", df = ", as.numeric(stats$breusch_pagan[[3]]$parameter),
               ", p-value = ", signif(as.numeric(stats$breusch_pagan[[3]]$p.value), digits)), "\n")
    cat(paste0("Robust Breusch-Pagan test:   BP = ", signif(as.numeric(stats$breusch_pagan[[4]]$statistic), digits),
               ", df = ", as.numeric(stats$breusch_pagan[[4]]$parameter),
               ", p-value = ", signif(as.numeric(stats$breusch_pagan[[4]]$p.value), digits)), "\n")
    cat("\n")
  }

  # Non-linear Correlation (Alternative)
  if (!is.null(stats$nlcor.xy)) {
    nl <- stats$nlcor.xy
    func <- nl$func
    if (is.na(suppressWarnings(as.numeric(func)))) {
      cat(paste("Non-linear Correlation Tests:", func, "Transformation\n"))
      labs <- c(
        paste0("Cor[", func, "(", "r_", varnames[2], "), ", varnames[1], "]"),
        paste0("Cor[", "r_", varnames[2], ", ", func, "(", varnames[1], ")]"),
        paste0("Cor[", func, "(", "r_", varnames[2], "), ", func, "(", varnames[1], ")]")
      )
    } else {
      cat(paste("Non-linear Correlation Tests: Power Transformation using", func, "\n"))
      labs <- c(
        paste0("Cor[r_", varnames[2], "^", func, ", ", varnames[1], "]"),
        paste0("Cor[r_", varnames[2], ", ", varnames[1], "^", func, "]"),
        paste0("Cor[r_", varnames[2], "^", func, ", ", varnames[1], "^", func, "]")
      )
    }

    tab <- rbind(as.numeric(nl$t1), as.numeric(nl$t2), as.numeric(nl$t3))
    colnames(tab) <- c("estimate", "t-value", "df", "Pr(>|t|)")
    rownames(tab) <- labs
    print.default(format(tab, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  # Difference Statistics CIs
  if (!is.null(stats$diff_matrix)) {
    boot_type <- "Percentile"
    if (!is.null(object$bagged_results[[1]]$boot.args)) {
      type_arg <- object$bagged_results[[1]]$boot.args[1]
      if (!is.null(type_arg)) {
        if(type_arg == "bca") boot_type <- "BCa"
        if(type_arg == "perc") boot_type <- "Percentile"
      }
    }
    conf_lev <- "95%"
    if (!is.null(object$bagged_results[[1]]$boot.args)) {
      lev <- as.numeric(object$bagged_results[[1]]$boot.args[2])
      if(!is.na(lev)) conf_lev <- paste0(lev*100, "%")
    }

    cat(paste0(conf_lev, " ", boot_type, " Bootstrap CIs for Difference Statistics (", object$n_valid_iterations, " samples):\n"))

    dmat <- stats$diff_matrix
    # ensure it is strictly numeric and un-named
    dmat_clean <- matrix(as.numeric(dmat), nrow=nrow(dmat), ncol=ncol(dmat))
    rownames(dmat_clean) <- rownames(dmat)
    colnames(dmat_clean) <- colnames(dmat)

    print.default(format(dmat_clean, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  cat("---\n")
  cat(paste("Note: Difference statistics > 0 suggest", varnames[2], "->", varnames[1], "\n"))

  invisible(object)
}

#' Print for dda_bagging Output (RESDIST)
#' @export
#' @method print dda_bagging_resdist
print.dda_bagging_resdist <- function(object, digits = 4) {
  stats <- object$aggregated_stats
  varnames <- if (!is.null(stats$var.names)) stats$var.names else c("target", "alternative")

  cat("\nBOOTSTRAP AGGREGATED DDA: Residual Distributions\n")
  cat("Number of Bootstrap Samples:", object$n_valid_iterations, "\n\n")

  cat("Skewness and kurtosis tests:\n")

  skew_row <- as.numeric(c(stats$agostino.target.statistic, stats$agostino.target.z, stats$agostino.target.p.value,
                           stats$agostino.alternative.statistic, stats$agostino.alternative.z, stats$agostino.alternative.p.value))
  kurt_row <- as.numeric(c(stats$anscombe.target.statistic, stats$anscombe.target.z, stats$anscombe.target.p.value,
                           stats$anscombe.alternative.statistic, stats$anscombe.alternative.z, stats$anscombe.alternative.p.value))

  sigtests <- rbind(skew_row, kurt_row)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")
  print.default(format(sigtests, digits = digits, scientific = NA, scipen = 999), print.gap = 2L, quote = FALSE)

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

  cat(paste0("Skewness and kurtosis difference tests and ", conf_lev, " ", boot_type, " bootstrap CIs:\n\n"))

  citests <- rbind(as.numeric(stats$skewdiff), as.numeric(stats$kurtdiff))
  rownames(citests) <- c("Skewness", "Kurtosis")
  if (ncol(citests) == 5) colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
  if (ncol(citests) == 3) colnames(citests) <- c("diff", "lower", "upper")
  print.default(format(citests, digits = digits), print.gap = 2L, quote = FALSE)

  cat(paste0("\n", conf_lev, " ", boot_type, " bootstrap CIs for joint higher moment differences:\n"))
  joint_mat <- rbind(as.numeric(stats$cor12diff), as.numeric(stats$RHS3), as.numeric(stats$cor13diff), as.numeric(stats$RHS4), as.numeric(stats$RCC))
  rownames(joint_mat) <- c("Co-Skewness", "Hyvarinen-Smith (Co-Skewness)", "Co-Kurtosis", "Hyvarinen-Smith (Co-Kurtosis)", "Chen-Chan (Co-Kurtosis)")
  colnames(joint_mat) <- c("estimate", "lower", "upper")
  print.default(format(joint_mat, digits = digits), print.gap = 2L, quote = FALSE)

  cat("\n---\n")
  probtrans <- if(!is.null(stats$probtrans)) stats$probtrans else FALSE
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], "\n"))
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], "\n"))
  if(isFALSE(probtrans)){
    cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1], "\n"))
  } else {
    cat(paste("      Under prob.trans = TRUE, skewness and kurtosis differences < 0 and\n",
              "      co-skewness and co-kurtosis differences > 0 suggest", varnames[2], "->", varnames[1], "\n"))
  }

  invisible(object)
}

#' Print for dda_bagging Output (VARDIST)
#' @export
#' @method print dda_bagging_vardist
print.dda_bagging_vardist <- function(object, digits = 4) {
  stats <- object$aggregated_stats
  varnames <- if (!is.null(stats$var.names)) stats$var.names else c("Outcome", "Predictor")

  cat("\nBOOTSTRAP AGGREGATED DDA: Variable Distributions\n")
  cat("Number of Bootstrap Samples:", object$n_valid_iterations, "\n\n")

  cat("Skewness and kurtosis tests:\n")

  skew_row <- as.numeric(c(stats$agostino.outcome.statistic.skew, stats$agostino.outcome.statistic.z, stats$agostino.outcome.p.value,
                           stats$agostino.predictor.statistic.skew, stats$agostino.predictor.statistic.z, stats$agostino.predictor.p.value))
  kurt_row <- as.numeric(c(stats$anscombe.outcome.statistic.kurt, stats$anscombe.outcome.statistic.z, stats$anscombe.outcome.p.value,
                           stats$anscombe.predictor.statistic.kurt, stats$anscombe.predictor.statistic.z, stats$anscombe.predictor.p.value))

  sigtests <- rbind(skew_row, kurt_row)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")
  print.default(format(sigtests, digits = digits), print.gap = 2L, quote = FALSE)

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

  cat(paste0(conf_lev, " ", boot_type, " bootstrap CIs for higher moment differences:\n"))
  citests <- rbind(as.numeric(stats$skewdiff), as.numeric(stats$kurtdiff))
  rownames(citests) <- c("Skewness", "Kurtosis")
  colnames(citests) <- c("diff", "lower", "upper")
  print.default(format(citests, digits = digits), print.gap = 2L, quote = FALSE)

  cat(paste0("\n", conf_lev, " ", boot_type, " bootstrap CIs for differences in higher-order correlations:\n"))
  hoctests <- rbind(as.numeric(stats$cor12diff), as.numeric(stats$cor13diff))
  rownames(hoctests) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]")
  colnames(hoctests) <- c("estimate", "lower", "upper")
  print.default(format(hoctests, digits = digits), print.gap = 2L, quote = FALSE)

  cat(paste0("\n", conf_lev, " ", boot_type, " bootstrap CIs for Likelihood Ratio approximations:\n"))
  LRtests <- rbind(as.numeric(stats$RHS), as.numeric(stats$Rtanh), as.numeric(stats$RCC))
  rownames(LRtests) <- c("Hyvarinen-Smith (co-skewness)", "Hyvarinen-Smith (co-kurtosis)", "Chen-Chan (co-kurtosis)")
  colnames(LRtests) <- c("estimate", "lower", "upper")
  print.default(format(LRtests, digits = digits), print.gap = 2L, quote = FALSE)

  cat("\n---\n")
  cat(paste("Note: (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests the model", varnames[2], "->", varnames[1], "\n"))

  invisible(object)
}

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
