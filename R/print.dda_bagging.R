### TO DO:
### Get output to match classic print styles for
### dda.resdist and dda.vardist and dda.indep



#' Summary for dda_bagging Output (INDEP, robust S3, no type check)
#'
#' @param object Output from dda_bagging() for dda.indep objects (class: dda_bagging_indep)
#' @param digits Number of digits for rounding (default: 4)
#' @param alpha Significance level for decisions (default: 0.05)
#' @export
#' @method summary dda_bagging_indep
#'
print.dda_bagging_indep <- function(object, digits = 4, alpha = 0.05) {
  # helpers
  has_value <- function(x) {
    if (is.null(x)) return(FALSE)
    if (is.matrix(x) || is.data.frame(x)) return(any(!is.na(as.vector(as.matrix(x))) & !is.nan(as.vector(as.matrix(x)))))
    if (is.list(x)) return(any(sapply(x, function(z) !is.null(z) && !all(is.na(z) | is.nan(z)))))
    if (length(x) == 0) return(FALSE)
    return(!all(is.na(x) | is.nan(x)))
  }

  # Print matrix-like objects only keeping rows/cols that have values
  print_clean_matrix <- function(mat, digits = 4, header = NULL) {
    if (is.null(mat)) return(FALSE)
    matm <- as.matrix(mat)
    # drop rows and cols that are entirely NA or NaN
    keep_rows <- apply(matm, 1, function(r) any(!is.na(r) & !is.nan(r)))
    keep_cols <- apply(matm, 2, function(c) any(!is.na(c) & !is.nan(c)))
    matm <- matm[keep_rows, keep_cols, drop = FALSE]
    if (nrow(matm) == 0 || ncol(matm) == 0) return(FALSE)
    if (!is.null(header)) cat(header, "\n")
    print.default(format(round(matm, digits), digits = max(3L, getOption("digits") - 3L)),
                  print.gap = 2L, quote = FALSE)
    return(TRUE)
  }

  stats <- object$aggregated_stats
  decisions <- object$decision_percentages

  # try to find variable names
  varnames <- NULL
  if (!is.null(object$var.names)) varnames <- object$var.names
  if (is.null(varnames) && !is.null(object$parameters$var.names)) varnames <- object$parameters$var.names
  if (is.null(varnames) && !is.null(stats$var.names)) varnames <- stats$var.names
  if (is.null(varnames)) varnames <- c("y", "x")

  # attempt to find iteration / bootstrap count for subheader
  iter_label <- NULL
  if (!is.null(object$parameters$iter)) iter_label <- object$parameters$iter
  if (is.null(iter_label) && !is.null(object$parameters$B)) iter_label <- object$parameters$B
  if (is.null(iter_label) && !is.null(object$n_valid_iterations)) iter_label <- object$n_valid_iterations
  if (is.null(iter_label) && !is.null(object$parameters$iterations)) iter_label <- object$parameters$iterations

  cat("\n")
  cat("DIRECTION DEPENDENCE ANALYSIS: Independence Properties", "\n", "\n")
  if (!is.null(iter_label)) {
    # per request: no quotes around the number
    cat(paste0("bootstrapped ", iter_label, " times"), "\n", "\n")
  }

  # -------------------- Target Model (yx) --------------------
  cat(paste("Target Model:", varnames[2], "->", varnames[1], sep = " "), "\n", "\n")

  # Omnibus independence tests header: only print if at least one test present
  # HSIC and dCor stats may be stored under various names in aggregated_stats
  hsic_yx_stat <- NULL; hsic_yx_pval <- NULL
  if (!is.null(stats$hsic_yx_stat)) hsic_yx_stat <- stats$hsic_yx_stat
  if (is.null(hsic_yx_stat) && !is.null(stats$hsic_yx_statistic)) hsic_yx_stat <- stats$hsic_yx_statistic
  if (is.null(hsic_yx_stat) && !is.null(stats$hsic.yx.stat)) hsic_yx_stat <- stats$hsic.yx.stat
  if (!is.null(stats$hsic_yx_pval)) hsic_yx_pval <- stats$hsic_yx_pval
  if (is.null(hsic_yx_pval) && !is.null(stats$hsic_yx_p.value)) hsic_yx_pval <- stats$hsic_yx_p.value
  if (is.null(hsic_yx_pval) && !is.null(stats$hsic.yx.p.value)) hsic_yx_pval <- stats$hsic.yx.p.value

  dcor_yx_stat <- NULL; dcor_yx_pval <- NULL
  if (!is.null(stats$dcor_yx_stat)) dcor_yx_stat <- stats$dcor_yx_stat
  if (is.null(dcor_yx_stat) && !is.null(stats$dcor_yx_statistic)) dcor_yx_stat <- stats$dcor_yx_statistic
  if (is.null(dcor_yx_stat) && !is.null(stats$dcor.yx.stat)) dcor_yx_stat <- stats$dcor.yx.stat
  if (!is.null(stats$dcor_yx_pval)) dcor_yx_pval <- stats$dcor_yx_pval
  if (is.null(dcor_yx_pval) && !is.null(stats$dcor_yx_p.value)) dcor_yx_pval <- stats$dcor_yx_p.value
  if (is.null(dcor_yx_pval) && !is.null(stats$dcor.yx.p.value)) dcor_yx_pval <- stats$dcor.yx.p.value

  if (has_value(hsic_yx_stat) || has_value(hsic_yx_pval) || has_value(dcor_yx_stat) || has_value(dcor_yx_pval)) {
    cat("Omnibus Independence Tests:", "\n")
    if (has_value(hsic_yx_stat) || has_value(hsic_yx_pval)) {
      cat(paste0("HSIC = ",
                 ifelse(has_value(hsic_yx_stat), format(round(as.numeric(hsic_yx_stat), digits), nsmall = digits), "NA"),
                 ", p-value = ",
                 ifelse(has_value(hsic_yx_pval), format(round(as.numeric(hsic_yx_pval), digits), nsmall = digits), "NA")))
      cat("\n")
    }
    # Only print dCor if present and not NA (user requested omit if NA)
    if (has_value(dcor_yx_stat) || has_value(dcor_yx_pval)) {
      if (has_value(dcor_yx_stat) || has_value(dcor_yx_pval)) {
        cat(paste0("dCor = ",
                   ifelse(has_value(dcor_yx_stat), format(round(as.numeric(dcor_yx_stat), digits), nsmall = digits), "NA"),
                   ", p-value = ",
                   ifelse(has_value(dcor_yx_pval), format(round(as.numeric(dcor_yx_pval), digits), nsmall = digits), "NA")))
        cat("\n")
      }
    }
    cat("\n")
  } else {
    cat("No omnibus independence statistics available for Target Model.\n\n")
  }

  # Homoscedasticity Tests (target)
  # support different naming conventions in aggregated stats
  bp_yx <- NULL
  if (!is.null(stats$breusch_pagan_yx)) bp_yx <- stats$breusch_pagan_yx
  if (is.null(bp_yx) && !is.null(stats$breusch_pagan)) {
    # if aggregated had a combined list of 4 like dda.indep, take first two as yx
    if (is.list(stats$breusch_pagan) && length(stats$breusch_pagan) >= 2) {
      bp_yx <- list(stats$breusch_pagan[[1]], stats$breusch_pagan[[2]])
    }
  }
  if (has_value(bp_yx)) {
    cat("Homoscedasticity Tests:", "\n")
    # try to construct matrix as dda.indep did
    if (is.list(bp_yx) && length(bp_yx) >= 2 &&
        all(sapply(bp_yx[1:2], function(z) !is.null(z$statistic) && !is.null(z$parameter) && !is.null(z$p.value)))) {
      sigtests1.yx <- rbind(c(bp_yx[[1]]$statistic, bp_yx[[1]]$parameter, bp_yx[[1]]$p.value),
                            c(bp_yx[[2]]$statistic, bp_yx[[2]]$parameter, bp_yx[[2]]$p.value))
      sigtests1.yx <- round(sigtests1.yx, digits)
      rownames(sigtests1.yx) <- c("BP-test", "Robust BP-test")
      colnames(sigtests1.yx) <- c("X-squared", "df", "p-value")
      print.default(format(sigtests1.yx, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
      cat("\n")
    } else if (is.matrix(bp_yx) || is.data.frame(bp_yx)) {
      print_clean_matrix(bp_yx, digits = digits)
      cat("\n")
    }
  }

  # Non-linear correlation tests (target)
  nl_yx <- NULL
  if (!is.null(stats$nlcor_yx)) nl_yx <- stats$nlcor_yx
  if (is.null(nl_yx) && !is.null(stats$nlcor.yx)) nl_yx <- stats$nlcor.yx
  if (has_value(nl_yx) && !is.null(nl_yx$t1) && !is.null(nl_yx$t2) && !is.null(nl_yx$t3)) {
    sigtests2.yx <- rbind(nl_yx$t1, nl_yx$t2, nl_yx$t3)
    sigtests2.yx <- round(sigtests2.yx, digits)
    if (is.null(nl_yx$func) && !is.null(nl_yx$fname)) nl_yx$func <- nl_yx$fname
    if (is.na(suppressWarnings(as.numeric(nl_yx$func)))) {
      cat(paste("Non-linear Correlation Tests:", nl_yx$func, "Transformation"), "\n\n")
      rownames(sigtests2.yx) <- c(paste("Cor[", nl_yx$func, "(", "r_", varnames[1], "), ", varnames[2], "]", sep = ""),
                                  paste("Cor[", "r_", varnames[1], ", ", nl_yx$func, "(", varnames[2], ")]", sep = ""),
                                  paste("Cor[", nl_yx$func, "(", "r_", varnames[1], "), ", nl_yx$func, "(", varnames[2], ")]", sep = ""))
    } else {
      cat(paste("Non-linear Correlation Tests: Power Transformation using", nl_yx$func), "\n\n")
      rownames(sigtests2.yx) <- c(paste("Cor[", "r_", varnames[1], "^", nl_yx$func, ", ", varnames[2], "]", sep = ""),
                                  paste("Cor[", "r_", varnames[1], ", ", varnames[2], "^", nl_yx$func, "]", sep = ""),
                                  paste("Cor[", "r_", varnames[1], "^", nl_yx$func, ", ", varnames[2], "^", nl_yx$func, "]", sep = ""))
    }
    colnames(sigtests2.yx) <- c("estimate", "t-value", "df", "Pr(>|t|)")
    print.default(format(sigtests2.yx, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  # -------------------- Alternative Model (xy) --------------------
  cat(paste("Alternative Model:", varnames[1], "->", varnames[2], sep = " "), "\n", "\n")

  hsic_xy_stat <- NULL; hsic_xy_pval <- NULL
  if (!is.null(stats$hsic_xy_stat)) hsic_xy_stat <- stats$hsic_xy_stat
  if (is.null(hsic_xy_stat) && !is.null(stats$hsic_xy_statistic)) hsic_xy_stat <- stats$hsic_xy_statistic
  if (is.null(hsic_xy_stat) && !is.null(stats$hsic.xy.stat)) hsic_xy_stat <- stats$hsic.xy.stat
  if (!is.null(stats$hsic_xy_pval)) hsic_xy_pval <- stats$hsic_xy_pval
  if (is.null(hsic_xy_pval) && !is.null(stats$hsic_xy_p.value)) hsic_xy_pval <- stats$hsic_xy_p.value
  if (is.null(hsic_xy_pval) && !is.null(stats$hsic.xy.p.value)) hsic_xy_pval <- stats$hsic.xy.p.value

  dcor_xy_stat <- NULL; dcor_xy_pval <- NULL
  if (!is.null(stats$dcor_xy_stat)) dcor_xy_stat <- stats$dcor_xy_stat
  if (is.null(dcor_xy_stat) && !is.null(stats$dcor_xy_statistic)) dcor_xy_stat <- stats$dcor_xy_statistic
  if (is.null(dcor_xy_stat) && !is.null(stats$dcor.xy.stat)) dcor_xy_stat <- stats$dcor.xy.stat
  if (!is.null(stats$dcor_xy_pval)) dcor_xy_pval <- stats$dcor_xy_pval
  if (is.null(dcor_xy_pval) && !is.null(stats$dcor_xy_p.value)) dcor_xy_pval <- stats$dcor_xy_p.value
  if (is.null(dcor_xy_pval) && !is.null(stats$dcor.xy.p.value)) dcor_xy_pval <- stats$dcor.xy.p.value

  if (has_value(hsic_xy_stat) || has_value(hsic_xy_pval) || has_value(dcor_xy_stat) || has_value(dcor_xy_pval)) {
    cat("Omnibus Independence Tests:", "\n")
    if (has_value(hsic_xy_stat) || has_value(hsic_xy_pval)) {
      cat(paste0("HSIC = ",
                 ifelse(has_value(hsic_xy_stat), format(round(as.numeric(hsic_xy_stat), digits), nsmall = digits), "NA"),
                 ", p-value = ",
                 ifelse(has_value(hsic_xy_pval), format(round(as.numeric(hsic_xy_pval), digits), nsmall = digits), "NA")))
      cat("\n")
    }
    if (has_value(dcor_xy_stat) || has_value(dcor_xy_pval)) {
      cat(paste0("dCor = ",
                 ifelse(has_value(dcor_xy_stat), format(round(as.numeric(dcor_xy_stat), digits), nsmall = digits), "NA"),
                 ", p-value = ",
                 ifelse(has_value(dcor_xy_pval), format(round(as.numeric(dcor_xy_pval), digits), nsmall = digits), "NA")))
      cat("\n")
    }
    cat("\n")
  } else {
    cat("No omnibus independence statistics available for Alternative Model.\n\n")
  }

  # Homoscedasticity Tests (alternative)
  bp_xy <- NULL
  if (!is.null(stats$breusch_pagan_xy)) bp_xy <- stats$breusch_pagan_xy
  if (is.null(bp_xy) && !is.null(stats$breusch_pagan) && is.list(stats$breusch_pagan) && length(stats$breusch_pagan) >= 4) {
    bp_xy <- list(stats$breusch_pagan[[3]], stats$breusch_pagan[[4]])
  }
  if (has_value(bp_xy)) {
    cat("Homoscedasticity Tests:", "\n")
    if (is.list(bp_xy) && length(bp_xy) >= 2 &&
        all(sapply(bp_xy[1:2], function(z) !is.null(z$statistic) && !is.null(z$parameter) && !is.null(z$p.value)))) {
      sigtests1.xy <- rbind(c(bp_xy[[1]]$statistic, bp_xy[[1]]$parameter, bp_xy[[1]]$p.value),
                            c(bp_xy[[2]]$statistic, bp_xy[[2]]$parameter, bp_xy[[2]]$p.value))
      sigtests1.xy <- round(sigtests1.xy, digits)
      rownames(sigtests1.xy) <- c("BP-test", "Robust BP-test")
      colnames(sigtests1.xy) <- c("X-squared", "df", "p-value")
      print.default(format(sigtests1.xy, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
      cat("\n")
    } else if (is.matrix(bp_xy) || is.data.frame(bp_xy)) {
      print_clean_matrix(bp_xy, digits = digits)
      cat("\n")
    }
  }

  # Non-linear correlation tests (alternative)
  nl_xy <- NULL
  if (!is.null(stats$nlcor_xy)) nl_xy <- stats$nlcor_xy
  if (is.null(nl_xy) && !is.null(stats$nlcor.xy)) nl_xy <- stats$nlcor.xy
  if (has_value(nl_xy) && !is.null(nl_xy$t1) && !is.null(nl_xy$t2) && !is.null(nl_xy$t3)) {
    sigtests2.xy <- rbind(nl_xy$t1, nl_xy$t2, nl_xy$t3)
    sigtests2.xy <- round(sigtests2.xy, digits)
    if (is.null(nl_xy$func) && !is.null(nl_xy$fname)) nl_xy$func <- nl_xy$fname
    if (is.na(suppressWarnings(as.numeric(nl_xy$func)))) {
      cat(paste("Non-linear Correlation Tests:", nl_xy$func, "Transformation"), "\n\n")
      rownames(sigtests2.xy) <- c(paste("Cor[", nl_xy$func, "(", "r_", varnames[2], "), ", varnames[1], "]", sep = ""),
                                  paste("Cor[", "r_", varnames[2], ", ", nl_xy$func, "(", varnames[1], ")]", sep = ""),
                                  paste("Cor[", nl_xy$func, "(", "r_", varnames[2], "), ", nl_xy$func, "(", varnames[1], ")]", sep = ""))
    } else {
      cat(paste("Non-linear Correlation Tests: Power Transformation using", nl_xy$func), "\n\n")
      rownames(sigtests2.xy) <- c(paste("Cor[", "r_", varnames[2], "^", nl_xy$func, ", ", varnames[1], "]", sep = ""),
                                  paste("Cor[", "r_", varnames[2], ", ", varnames[1], "^", nl_xy$func, "]", sep = ""),
                                  paste("Cor[", "r_", varnames[2], "^", nl_xy$func, ", ", varnames[1], "^", nl_xy$func, "]", sep = ""))
    }
    colnames(sigtests2.xy) <- c("estimate", "t-value", "df", "Pr(>|t|)")
    print.default(format(sigtests2.xy, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  # Difference statistics (bootstrap CIs) - try to follow dda.indep print style
  out_diff <- NULL
  if (!is.null(stats$out.diff)) out_diff <- stats$out.diff
  if (is.null(out_diff) && !is.null(stats$diff_matrix)) out_diff <- stats$diff_matrix
  # boot.args may be in several places
  boot_args <- NULL
  if (!is.null(object$boot.args)) boot_args <- object$boot.args
  if (is.null(boot_args) && !is.null(object$parameters$boot.args)) boot_args <- object$parameters$boot.args
  if (is.null(boot_args) && !is.null(stats$boot.args)) boot_args <- stats$boot.args
  if (has_value(out_diff)) {
    ci_level <- NA
    if (!is.null(boot_args) && length(boot_args) >= 2) ci_level <- as.numeric(boot_args[2]) * 100
    # choose label for type
    boot_type <- if (!is.null(boot_args) && length(boot_args) >= 1) as.character(boot_args[1]) else NULL
    n_samples <- if (!is.null(boot_args) && length(boot_args) >= 3) boot_args[3] else NULL

    if (!is.null(ci_level) && !is.na(ci_level) && !is.null(boot_type) && !is.null(n_samples)) {
      pretty_type <- if (boot_type == "bca") "BCa" else if (boot_type == "perc") "Percentile" else boot_type
      cat(paste0(ci_level, "% ", pretty_type, " Bootstrap CIs for Difference Statistics",
                 if (!is.null(n_samples) && !is.na(n_samples)) paste0(" (", n_samples, " samples):") else ":"), "\n")
    } else {
      cat("Bootstrap CIs for Difference Statistics:\n")
    }

    # out_diff expected to be a matrix: estimate, lower, upper
    printed <- FALSE
    if (is.matrix(out_diff) || is.data.frame(out_diff)) {
      printed <- print_clean_matrix(out_diff, digits = digits)
    } else if (is.list(out_diff) && length(out_diff) > 0) {
      # try to coerce if it's a one-element list with matrix
      if (length(out_diff) == 1 && (is.matrix(out_diff[[1]]) || is.data.frame(out_diff[[1]]))) {
        printed <- print_clean_matrix(out_diff[[1]], digits = digits)
      }
    }
    if (printed) {
      cat("\n")
      cat("---\n")
      cat("\n")
      # Note about interpretation similar to dda.indep (direction based on difference > 0)
      cat(paste("Note: Difference statistics > 0 suggest", varnames[2], "->", varnames[1], sep = " "), "\n")
    } else {
      cat("No difference statistics available to print.\n")
    }
  }

  # Decision proportions (bagging-level) - only print non-empty, non-NA proportions
  if (!is.null(decisions) && length(decisions) > 0) {
    for (dname in names(decisions)) {
      prop <- decisions[[dname]]
      out <- c(undecided = 0, `y->x` = 0, `x->y` = 0)
      if (!is.null(prop) && length(prop) > 0) {
        prop_names <- names(prop)
        out[prop_names] <- prop
      }
      if (any(!is.na(out) & out != 0)) {
        cat(paste("\nDecision proportions for", dname, ":\n"))
        print(round(out, digits))
        cat("---\n")
      }
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
#' Print method for bagged residual-distribution results (resdist)
#' Robust to missing/NA components and attempts to match print.dda.resdist output.

print.dda_bagging_resdist <- function(object, digits = 4) {
  # Helpers
  has_value <- function(x) {
    if (is.null(x)) return(FALSE)
    if (is.matrix(x) || is.data.frame(x)) return(any(!is.na(as.vector(as.matrix(x))) & !is.nan(as.vector(as.matrix(x)))))
    if (is.list(x)) return(any(sapply(x, function(z) !is.null(z) && !all(is.na(z) | is.nan(z)))))
    if (length(x) == 0) return(FALSE)
    return(!all(is.na(x) | is.nan(x)))
  }
  first_nonnull <- function(...) {
    vals <- list(...)
    for (v in vals) if (!is.null(v)) return(v)
    return(NULL)
  }

  stats <- object$aggregated_stats
  params <- if (!is.null(object$parameters)) object$parameters else list()

  # varnames: try aggregated stats, then parameters, else sensible defaults
  varnames <- NULL
  if (!is.null(stats$var.names) && length(stats$var.names) == 2) varnames <- stats$var.names
  if (is.null(varnames) && !is.null(params$var.names) && length(params$var.names) == 2) varnames <- params$var.names
  if (is.null(varnames) && !is.null(object$var.names) && length(object$var.names) == 2) varnames <- object$var.names
  if (is.null(varnames)) varnames <- c("target", "alternative")

  # iteration / counts for preamble
  iter_val <- first_nonnull(params$iter, params$B, params$iterations, object$n_valid_iterations)
  completed_iter <- if (!is.null(object$n_valid_iterations)) object$n_valid_iterations else NA

  # Preamble similar to indep-style bagging summary
  # cat("\n===== DDA Bagging Summary =====\n")
  # cat("Function:", if (!is.null(params$function_name)) params$function_name else "unknown", "\n")
  # cat("Object Type:", if (!is.null(params$object_type)) params$object_type else "unknown", "\n")
  # cat("Iterations:", if (!is.null(iter_val)) iter_val else NA, "\n")
  # cat("Completed Iterations:", if (!is.null(completed_iter)) completed_iter else NA, "\n")
  # cat("----\n\n")

  # Classic resdist header
  cat("DIRECTION DEPENDENCE ANALYSIS: Residual Distributions (Bagged)\n\n")
  cat("Bootstrapped", iter_val, "times", "\n")


  cat("Skewness and kurtosis tests:\n")

  # Build skew/kurtosis rows robustly; try several possible names used in aggregated_stats
  # For skewness (agostino) and kurtosis (anscombe) we look for target/alternative values
  # Accept either aggregated naming like agostino.target.statistic or nested structures
  get_val <- function(prefix, side, what) {
    # prefix examples: "agostino", "anscombe"
    # side: "target" or "alternative"; what: "statistic","z","p.value"
    # possible keys:
    keys <- c(
      paste0(prefix, ".", side, ".", what),
      paste0(prefix, ".", side, ".", what, ifelse(what == "p.value", "", "")),
      paste0(prefix, ".", side, ".", what), # fallback duplicate but safe
      paste0(prefix, ".", side, ".", what),
      paste0(prefix, ".", side, ".", what) # keep in case of slight variations
    )
    for (k in keys) {
      if (!is.null(stats[[k]])) return(stats[[k]])
    }
    # handle case where stats$agostino is a list with $target and $alternative entries
    if (!is.null(stats[[prefix]]) && is.list(stats[[prefix]])) {
      entry <- stats[[prefix]]
      if (!is.null(entry[[side]]) && is.list(entry[[side]])) {
        if (!is.null(entry[[side]][[what]])) return(entry[[side]][[what]])
        # sometimes p-value named "p.value" or "pvalue"
        if (what == "p.value" && !is.null(entry[[side]][["pvalue"]])) return(entry[[side]][["pvalue"]])
      }
    }
    return(NA)
  }

  skew_row <- c(
    get_val("agostino", "target", "statistic"),
    get_val("agostino", "target", "z"),
    get_val("agostino", "target", "p.value"),
    get_val("agostino", "alternative", "statistic"),
    get_val("agostino", "alternative", "z"),
    get_val("agostino", "alternative", "p.value")
  )
  kurt_row <- c(
    get_val("anscombe", "target", "statistic"),
    get_val("anscombe", "target", "z"),
    get_val("anscombe", "target", "p.value"),
    get_val("anscombe", "alternative", "statistic"),
    get_val("anscombe", "alternative", "z"),
    get_val("anscombe", "alternative", "p.value")
  )

  sigtests <- rbind(skew_row, kurt_row)
  sigtests <- round(sigtests, digits)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")

  # drop rows that are entirely NA/NaN
  keep_rows <- apply(sigtests, 1, function(row) any(!is.na(row) & !is.nan(row)))
  sigtests_print <- sigtests[keep_rows, , drop = FALSE]
  if (nrow(sigtests_print) > 0) {
    print.default(format(sigtests_print, digits = max(3L, getOption("digits") - 3L), scientific = NA, scipen = 999),
                  print.gap = 2L, quote = FALSE)
  } else {
    cat("No skewness/kurtosis statistics available.\n")
  }

  # Skewness/kurtosis difference tests and bootstrap CIs detection
  # Determine boot args (type, level, n) from any bagged result element if present
  ci_level <- NULL; boot_type <- NULL; n_resamples <- NULL
  if (!is.null(object$bagged_results) && length(object$bagged_results) > 0) {
    for (res in object$bagged_results) {
      if (!is.null(res$boot.args)) {
        boot_type <- res$boot.args[1]
        ci_level <- as.numeric(res$boot.args[2]) * 100
        n_resamples <- as.numeric(res$boot.args[3])
        break
      }
    }
  }
  # fallback: maybe boot.args attached to top-level aggregated_stats or object
  if (is.null(boot_type)) {
    boot_args_top <- first_nonnull(object$boot.args, stats$boot.args, params$boot.args)
    if (!is.null(boot_args_top) && length(boot_args_top) >= 3) {
      boot_type <- boot_args_top[1]
      ci_level <- as.numeric(boot_args_top[2]) * 100
      n_resamples <- as.numeric(boot_args_top[3])
    }
  }

  # decide whether to print difference tests (presence of either skewdiff/kurtdiff)
  skewdiff <- if (!is.null(stats$skewdiff)) stats$skewdiff else NULL
  kurtdiff <- if (!is.null(stats$kurtdiff)) stats$kurtdiff else NULL
  if (has_value(skewdiff) || has_value(kurtdiff)) {
    cat("\n")
    if (!is.null(ci_level) && !is.null(boot_type)) {
      if (boot_type == "bca") {
        cat("Skewness and kurtosis difference tests and ", ci_level, "% BCa bootstrap CIs:\n\n", sep = "")
      } else if (boot_type == "perc") {
        cat("Skewness and kurtosis difference tests and ", ci_level, "% Percentile bootstrap CIs:\n\n", sep = "")
      } else {
        cat("Skewness and kurtosis difference tests and bootstrap CIs:\n\n")
      }
    } else {
      cat("Skewness and kurtosis difference tests and bootstrap CIs:\n\n")
    }

    # Prepare citests rows; if boot present include lower/upper
    # Normalize shape to (2 x 5) if lower/upper present, else (2 x 3)
    make_row <- function(x) {
      if (is.null(x)) return(rep(NA, 5))
      # if vector length 5 already, return; if length 3 return padded with NA for lower/upper
      if (length(x) == 5) return(as.numeric(x))
      if (length(x) == 3) return(c(as.numeric(x), NA, NA))
      # if named with lower/upper, try to pick elements
      xnum <- as.numeric(x)
      if (length(xnum) >= 5) return(xnum[1:5])
      # fallback pad
      out <- rep(NA, 5); out[1:length(xnum)] <- xnum; return(out)
    }

    srow <- make_row(skewdiff)
    krow <- make_row(kurtdiff)
    citests <- rbind(srow, krow)
    rownames(citests) <- c("Skewness", "Kurtosis")
    # choose colnames depending on whether lower/upper present anywhere
    if (any(!is.na(citests[, 4:5]))) {
      colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
    } else {
      colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
    }
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
  joint_vals <- lapply(joint_keys, function(k) {
    v <- stats[[k]]
    if (is.null(v)) return(rep(NA, 3))
    # expected 3-length numeric: estimate, lower, upper
    if (is.numeric(v) && length(v) >= 3) return(as.numeric(v[1:3]))
    if (is.matrix(v) || is.data.frame(v)) {
      # if matrix, take first row
      m <- as.matrix(v)
      if (nrow(m) >= 1 && ncol(m) >= 3) return(as.numeric(m[1, 1:3]))
    }
    return(rep(NA, 3))
  })
  joint_mat <- do.call(rbind, joint_vals)
  rownames(joint_mat) <- joint_names
  colnames(joint_mat) <- c("estimate", "lower", "upper")
  joint_mat <- round(joint_mat, digits)
  keep_rows <- apply(joint_mat, 1, function(row) any(!is.na(row) & !is.nan(row)))
  joint_print <- joint_mat[keep_rows, , drop = FALSE]
  if (nrow(joint_print) > 0) {
    cat("\n")
    if (!is.null(ci_level) && !is.null(boot_type)) {
      if (boot_type == "bca") cat(ci_level, "% BCa bootstrap CIs for joint higher moment differences:\n", sep = "")
      if (boot_type == "perc") cat(ci_level, "% Percentile bootstrap CIs for joint higher moment differences:\n", sep = "")
    } else {
      cat("Bootstrap CIs for joint higher moment differences:\n")
    }
    print.default(format(joint_print, digits = max(3L, getOption("digits") - 3L)),
                  print.gap = 2L, quote = FALSE)
  }

  # Number of resamples (if available)
  if (!is.null(n_resamples) && !is.na(n_resamples)) {
    cat("\nNumber of resamples:", n_resamples, "\n")
  } else {
    # also check top-level boot.args
    boot_args_top <- first_nonnull(object$boot.args, stats$boot.args, params$boot.args)
    if (!is.null(boot_args_top) && length(boot_args_top) >= 3) {
      cat("\nNumber of resamples:", boot_args_top[3], "\n")
    }
  }

  cat("---\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], "\n"))
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], "\n"))

  probtrans_flag <- if (!is.null(stats$probtrans)) stats$probtrans else if (!is.null(object$probtrans)) object$probtrans else FALSE
  if (isFALSE(probtrans_flag)) {
    cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1], "\n"))
  } else {
    cat(paste("      Under prob.trans = TRUE, skewness and kurtosis differences < 0 and", "\n",
              "     co-skewness and co-kurtosis differences > 0 suggest", varnames[2], "->", varnames[1], "\n"))
  }

  # optional boot warning (match print.dda.resdist behavior)
  boot_warning_flag <- first_nonnull(stats$boot.warning, object$boot.warning, params$boot.warning)
  if (!is.null(boot_warning_flag) && isTRUE(boot_warning_flag)) {
    cat("\nWarning: Excess-kurtosis values of residuals have unequal signs\n")
    cat("         Also compute Co-Kurtosis and Hyvarinen-Smith Co-Kurtosis for", varnames[1], "->", varnames[2], "\n")
  }

  cat("\n")
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
print.dda_bagging_vardist <- function(object, digits = 4) {
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
