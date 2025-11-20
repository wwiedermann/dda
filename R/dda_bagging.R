data("mtcars")

n <- 1000
x <- rchisq(n, df = 4) - 4
e <- rchisq(n, df = 3) - 3
y <- 0.5 * x + e
d <- data.frame(x, y)

# 1. dda.indep: Direction Dependence Analysis - Independence
# result_indep <- dda.indep(mpg ~ wt + hp, pred = "wt", data = mtcars)
# print(result_indep)

res_ex_indep <- dda.indep(y ~ x, pred = "x", data = d, parallelize = TRUE, cores = 2,
                                     nlfun = 2, B = 50, hetero = TRUE, diff = TRUE)
print(res_ex_indep)

init <- Sys.time()
bagged_indep <- dda_bagging(res_ex_indep, iter = 15)
Sys.time() - init
#print.dda_bagging_indep(bagged_indep)
summary.dda_bagging_indep(bagged_indep)



# 2. dda.resdist: Direction Dependence Analysis - Residual Distribution

#result_resdist <- dda.resdist(mpg ~ wt + hp, pred = "wt", data = mtcars)
result_resdist <- dda.resdist(y ~ x, pred = "x", data = d,
                      B = 2000, conf.level = 0.90)


bagged_resdist <- dda_bagging(result_resdist, iter = 10)
Sys.time() - init # Time difference of 3.210665 mins
print.dda_bagging_resdist(bagged_resdist)

summary.dda_bagging_resdist(bagged_resdist)


###############################################################################
# 3. dda.vardist: Direction Dependence Analysis - Variable Distribution
#result_vardist <- dda.vardist(mpg ~ wt + hp, pred = "wt", data = mtcars)
result_vardist <- dda.vardist(y ~ x, pred = "x", data = d, B = 2000)
print(result_vardist)

init <- Sys.time()
bagged_vardist <- dda_bagging(result_vardist, iter = 100)
Sys.time() - init #Time difference of 2.703492 mins
summary.dda_bagging_vardist(bagged_vardist)


#############################################################################
#' Bootstrap Aggregated DDA Analysis (harmonic mean for p-values, mean for other stats)
#'
#' @param dda_result Output from any DDA function (dda.indep, dda.vardist, dda.resdist, etc.)
#' @param iter Number of bootstrap iterations (default: 100)
#' @param progress Whether to show progress bar (default: TRUE)
#' @param save_file Optional file path to save results
#' @param alpha Significance level for decisions (default: 0.05)
#' @return A list containing bootstrap and aggregated results
#' @export
dda_bagging <- function(
    dda_result,
    iter = 100,
    progress = TRUE,
    save_file = NULL,
    alpha = 0.05
) {
  get_numeric <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (is.numeric(x)) return(as.numeric(x[1]))
    if (is.list(x)) return(get_numeric(x[[1]]))
    return(NA_real_)
  }
  harmonic_p <- function(pvec) {
    pvec <- as.numeric(pvec)
    pvec <- pvec[!is.na(pvec) & pvec > 0 & pvec < 1]
    if (length(pvec) == 0) return(NA_real_)
    if (!requireNamespace("harmonicmeanp", quietly = TRUE)) {
      warning("harmonicmeanp package not found, using mean for p-values")
      return(mean(pvec))
    }
    harmonicmeanp::p.hmp(pvec, L = length(pvec))
  }

  # Helper to calculate decision proportions with fixed labels
  calc_props <- function(dec_vec) {
    levs <- c("Undecided", "Target", "Alternative")
    # Create factor with explicit levels to ensure table always has 3 entries (even if 0 count)
    dec_fac <- factor(dec_vec, levels = levs)
    tab <- table(dec_fac)
    prop <- tab / sum(tab)
    prop
  }

  # Extract call info and data
  if (inherits(dda_result, "dda.indep") ||
      inherits(dda_result, "dda.vardist") ||
      inherits(dda_result, "dda.resdist") ||
      inherits(dda_result, "cdda.indep") ||
      inherits(dda_result, "cdda.vardist")) {
    call_info <- dda_result$call_info
    function_name <- call_info$function_name
    all_args <- call_info$all_args
    data_name <- call_info$data_name
    original_data <- call_info$original_data
    if (is.null(original_data)) {
      tryCatch({
        original_data <- get(data_name, envir = parent.frame())
      }, error = function(e) {
        stop(paste("Could not retrieve original data:", data_name))
      })
    }
  } else if (length(dda_result) >= 5 && !is.null(dda_result[[5]])) {
    call_info <- dda_result[[5]]
    function_name <- call_info$function_name
    all_args <- call_info$all_args
    original_env <- call_info$environment
    data_name <- call_info$data_name
    original_data <- get(data_name, envir = original_env)
  } else {
    stop("Unrecognized DDA result structure. Please ensure you're using a supported DDA object.")
  }

  nobs <- nrow(original_data)
  dda_function <- get(function_name)

  bagged_results <- vector("list", iter)
  if (progress) {
    pb <- txtProgressBar(min = 0, max = iter, style = 3, width = 50, char = "=")
    cat(paste("\n", "Running", iter, "bootstrap iterations of", function_name, "\n"))
  }

  for(i in seq_len(iter)) {
    boot_indices <- sample(seq_len(nobs), nobs, replace = TRUE)
    boot_data <- original_data[boot_indices, ]
    boot_args <- all_args
    boot_args$data <- boot_data
    res <- tryCatch({
      do.call(dda_function, boot_args)
    }, error = function(e) NA)
    bagged_results[[i]] <- res
    if (progress) setTxtProgressBar(pb, i)
  }
  if (progress) { close(pb); cat("\nBootstrap iterations completed.\n") }

  valid_results <- bagged_results[!sapply(bagged_results, function(x) all(is.na(x)) || is.null(x))]
  object_type <- class(dda_result)[1]
  agg <- list()
  decisions <- list()
  n_valid <- length(valid_results)
  class_label <- "dda_bagging"

  # Critical value for Z-tests (approx 1.96 for alpha=0.05)
  crit_val <- qnorm(1 - alpha/2)

  ## --- DDA.INDEP block ---
  if (n_valid > 0 && object_type == "dda.indep") {
    hsic_yx <- sapply(valid_results, function(x) get_numeric(x$hsic.yx$statistic))
    hsic_xy <- sapply(valid_results, function(x) get_numeric(x$hsic.xy$statistic))
    hsic_yx_pval <- sapply(valid_results, function(x) get_numeric(x$hsic.yx$p.value))
    hsic_xy_pval <- sapply(valid_results, function(x) get_numeric(x$hsic.xy$p.value))

    dcor_yx <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_yx$statistic))
    dcor_xy <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_xy$statistic))
    dcor_yx_pval <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_yx$p.value))
    dcor_xy_pval <- sapply(valid_results, function(x) get_numeric(x$distance_cor$dcor_xy$p.value))

    agg$hsic_yx_stat <- mean(hsic_yx, na.rm = TRUE)
    agg$hsic_xy_stat <- mean(hsic_xy, na.rm = TRUE)
    agg$dcor_yx_stat <- mean(dcor_yx, na.rm = TRUE)
    agg$dcor_xy_stat <- mean(dcor_xy, na.rm = TRUE)

    agg$hsic_yx_pval <- harmonic_p(hsic_yx_pval)
    agg$hsic_xy_pval <- harmonic_p(hsic_xy_pval)
    agg$dcor_yx_pval <- harmonic_p(dcor_yx_pval)
    agg$dcor_xy_pval <- harmonic_p(dcor_xy_pval)

    # --- Decision Proportions (INDEP) ---
    # Logic:
    #   If Y|X is indep (p >= alpha) AND X|Y is dep (p < alpha) => x->y (Target)
    #   If Y|X is dep (p < alpha) AND X|Y is indep (p >= alpha) => y->x (Alternative)

    dec_hsic <- ifelse(hsic_yx_pval >= alpha & hsic_xy_pval < alpha, "Target",
                       ifelse(hsic_yx_pval < alpha & hsic_xy_pval >= alpha, "Alternative", "Undecided"))
    decisions$hsic <- calc_props(dec_hsic)

    dec_dcor <- ifelse(dcor_yx_pval >= alpha & dcor_xy_pval < alpha, "Target",
                       ifelse(dcor_yx_pval < alpha & dcor_xy_pval >= alpha, "Alternative", "Undecided"))
    decisions$dcor <- calc_props(dec_dcor)

    # Difference statistics if present
    if (!is.null(valid_results[[1]]$out.diff)) {
      first_diff <- valid_results[[1]]$out.diff
      diff_mat <- lapply(valid_results, function(x) {
        if (!is.null(x$out.diff)) {
          as.matrix(x$out.diff)
        } else {
          matrix(NA_real_, nrow = nrow(first_diff), ncol = ncol(first_diff),
                 dimnames = dimnames(first_diff))
        }
      })
      diff_array <- simplify2array(diff_mat)
      agg$diff_matrix <- apply(diff_array, c(1,2), function(xx) mean(xx, na.rm=TRUE))

      # Calculate difference decisions (if CI > 0 => Target)
      # We need row-wise extraction. out.diff row 1=HSIC, 2=dCor, 3=MI. cols 1=est, 2=lower, 3=upper

      calc_diff_decision <- function(row_idx) {
        lowers <- sapply(valid_results, function(x) {
          mat <- x$out.diff
          if (!is.null(mat) && nrow(mat) >= row_idx) mat[row_idx, 2] else NA
        })
        uppers <- sapply(valid_results, function(x) {
          mat <- x$out.diff
          if (!is.null(mat) && nrow(mat) >= row_idx) mat[row_idx, 3] else NA
        })
        d <- ifelse(!is.na(lowers) & !is.na(uppers) & lowers > 0 & uppers > 0, "Target",
                    ifelse(!is.na(lowers) & !is.na(uppers) & lowers < 0 & uppers < 0, "Alternative", "Undecided"))
        calc_props(d)
      }

      decisions$diff_hsic <- calc_diff_decision(1)
      if(nrow(first_diff) >= 2) decisions$diff_dcor <- calc_diff_decision(2)
      if(nrow(first_diff) >= 3) decisions$diff_mi   <- calc_diff_decision(3)
    }

    class_label <- "dda_bagging_indep"
  }

  ## --- DDA.RESDIST block ---
  if (n_valid > 0 && object_type == "dda.resdist") {
    agg$var.names <- if (!is.null(valid_results[[1]]$var.names)) valid_results[[1]]$var.names else c("target", "alternative")
    agg$probtrans <- if (!is.null(valid_results[[1]]$probtrans)) valid_results[[1]]$probtrans else FALSE

    # Collect stats for aggregation
    agg$agostino.target.statistic <- mean(sapply(valid_results, function(x) get_numeric(x$agostino$target$statistic[1])), na.rm=TRUE)
    agg$agostino.target.z         <- mean(sapply(valid_results, function(x) get_numeric(x$agostino$target$statistic[2])), na.rm=TRUE)
    agg$agostino.target.p.value   <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$agostino$target$p.value)))
    agg$agostino.alternative.statistic <- mean(sapply(valid_results, function(x) get_numeric(x$agostino$alternative$statistic[1])), na.rm=TRUE)
    agg$agostino.alternative.z         <- mean(sapply(valid_results, function(x) get_numeric(x$agostino$alternative$statistic[2])), na.rm=TRUE)
    agg$agostino.alternative.p.value   <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$agostino$alternative$p.value)))

    agg$anscombe.target.statistic <- mean(sapply(valid_results, function(x) get_numeric(x$anscombe$target$statistic[1])), na.rm=TRUE)
    agg$anscombe.target.z         <- mean(sapply(valid_results, function(x) get_numeric(x$anscombe$target$statistic[2])), na.rm=TRUE)
    agg$anscombe.target.p.value   <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$anscombe$target$p.value)))
    agg$anscombe.alternative.statistic <- mean(sapply(valid_results, function(x) get_numeric(x$anscombe$alternative$statistic[1])), na.rm=TRUE)
    agg$anscombe.alternative.z         <- mean(sapply(valid_results, function(x) get_numeric(x$anscombe$alternative$statistic[2])), na.rm=TRUE)
    agg$anscombe.alternative.p.value   <- harmonic_p(sapply(valid_results, function(x) get_numeric(x$anscombe$alternative$p.value)))

    mean_diff <- function(key, len) {
      out <- lapply(valid_results, function(x) {
        v <- x[[key]]
        if (is.null(v)) return(rep(NA, len))
        vnum <- suppressWarnings(as.numeric(v))
        if (length(vnum) != len) return(rep(NA, len))
        vnum
      })
      mat <- do.call(rbind, out)
      if (is.null(mat) || length(mat) == 0) return(rep(NA, len))
      colMeans(mat, na.rm = TRUE)
    }
    agg$skewdiff <- mean_diff("skewdiff", 5)
    agg$kurtdiff <- mean_diff("kurtdiff", 5)
    agg$cor12diff <- mean_diff("cor12diff", 3)
    agg$cor13diff <- mean_diff("cor13diff", 3)
    agg$RHS3      <- mean_diff("RHS3", 3)
    agg$RHS4      <- mean_diff("RHS4", 3)
    agg$RCC       <- mean_diff("RCC", 3)

    # --- Decision Proportions (RESDIST) ---

    calc_ci_decision <- function(key) {
      # x[[key]] is c(est, lower, upper)
      lowers <- sapply(valid_results, function(x) { v <- x[[key]]; if(length(v)>=2) v[2] else NA })
      uppers <- sapply(valid_results, function(x) { v <- x[[key]]; if(length(v)>=3) v[3] else NA })

      d <- ifelse(!is.na(lowers) & !is.na(uppers) & lowers > 0 & uppers > 0, "Target",
                  ifelse(!is.na(lowers) & !is.na(uppers) & lowers < 0 & uppers < 0, "Alternative", "Undecided"))
      calc_props(d)
    }

    decisions$dec_skewdiff  <- calc_ci_decision("skewdiff")
    decisions$dec_kurtdiff  <- calc_ci_decision("kurtdiff")
    decisions$dec_cor12diff <- calc_ci_decision("cor12diff")
    decisions$dec_cor13diff <- calc_ci_decision("cor13diff")
    decisions$dec_RHS3      <- calc_ci_decision("RHS3")
    decisions$dec_RCC       <- calc_ci_decision("RCC")
    decisions$dec_RHS4      <- calc_ci_decision("RHS4")

    # Agostino/Anscombe for ResDist (Errors)
    # Target (e_yx) should be Gaussian (low stat), Alt (e_xy) should be Non-Gaussian (high stat)
    # Logic: abs(z_target) < crit & abs(z_alt) >= crit => Target (x->y)

    ago_tar_z <- sapply(valid_results, function(x) get_numeric(x$agostino$target$statistic[2]))
    ago_alt_z <- sapply(valid_results, function(x) get_numeric(x$agostino$alternative$statistic[2]))
    dec_ago <- ifelse(abs(ago_tar_z) < crit_val & abs(ago_alt_z) >= crit_val, "Target",
                      ifelse(abs(ago_tar_z) >= crit_val & abs(ago_alt_z) < crit_val, "Alternative", "Undecided"))
    decisions$dec_agost <- calc_props(dec_ago)

    ans_tar_z <- sapply(valid_results, function(x) get_numeric(x$anscombe$target$statistic[2]))
    ans_alt_z <- sapply(valid_results, function(x) get_numeric(x$anscombe$alternative$statistic[2]))
    dec_ans <- ifelse(abs(ans_tar_z) < crit_val & abs(ans_alt_z) >= crit_val, "Target",
                      ifelse(abs(ans_tar_z) >= crit_val & abs(ans_alt_z) < crit_val, "Alternative", "Undecided"))
    decisions$dec_anscom <- calc_props(dec_ans)

    class_label <- "dda_bagging_resdist"
  }

  ## --- DDA.VARDIST block ---
  if (n_valid > 0 && object_type == "dda.vardist") {
    agg$var.names <- if (!is.null(valid_results[[1]]$var.names)) valid_results[[1]]$var.names else c("Outcome", "Predictor")

    # Helpers for aggregation
    get_mean <- function(f) {
      vals <- sapply(valid_results, function(x) {
        out <- tryCatch(f(x), error = function(e) NA_real_)
        get_numeric(out)
      })
      mean(vals, na.rm = TRUE)
    }
    get_hmp <- function(f) {
      pvals <- sapply(valid_results, function(x) {
        out <- tryCatch(f(x), error = function(e) NA_real_)
        get_numeric(out)
      })
      harmonic_p(pvals)
    }
    mean_vec <- function(key, len) {
      rows <- lapply(valid_results, function(x) {
        v <- tryCatch(x[[key]], error = function(e) NULL)
        if (is.null(v)) return(rep(NA_real_, len))
        vnum <- suppressWarnings(as.numeric(v))
        if (length(vnum) < len) {
          vnum <- c(vnum, rep(NA_real_, len - length(vnum)))
        } else if (length(vnum) > len) {
          vnum <- vnum[seq_len(len)]
        }
        vnum
      })
      mat <- do.call(rbind, rows)
      if (is.null(mat) || length(mat) == 0) return(rep(NA_real_, len))
      as.numeric(colMeans(mat, na.rm = TRUE))
    }

    agg$agostino.outcome.statistic.skew   <- get_mean(function(x) x$agostino$outcome$statistic[1])
    agg$agostino.outcome.statistic.z      <- get_mean(function(x) x$agostino$outcome$statistic[2])
    agg$agostino.outcome.p.value          <- get_hmp(function(x) x$agostino$outcome$p.value)

    agg$agostino.predictor.statistic.skew <- get_mean(function(x) x$agostino$predictor$statistic[1])
    agg$agostino.predictor.statistic.z    <- get_mean(function(x) x$agostino$predictor$statistic[2])
    agg$agostino.predictor.p.value        <- get_hmp(function(x) x$agostino$predictor$p.value)

    agg$anscombe.outcome.statistic.kurt   <- get_mean(function(x) x$anscombe$outcome$statistic[1])
    agg$anscombe.outcome.statistic.z      <- get_mean(function(x) x$anscombe$outcome$statistic[2])
    agg$anscombe.outcome.p.value          <- get_hmp(function(x) x$anscombe$outcome$p.value)

    agg$anscombe.predictor.statistic.kurt <- get_mean(function(x) x$anscombe$predictor$statistic[1])
    agg$anscombe.predictor.statistic.z    <- get_mean(function(x) x$anscombe$predictor$statistic[2])
    agg$anscombe.predictor.p.value        <- get_hmp(function(x) x$anscombe$predictor$p.value)

    agg$skewdiff  <- mean_vec("skewdiff",  3)
    agg$kurtdiff  <- mean_vec("kurtdiff",  3)
    agg$cor12diff <- mean_vec("cor12diff", 3)
    agg$cor13diff <- mean_vec("cor13diff", 3)
    agg$RHS       <- mean_vec("RHS",       3)
    agg$Rtanh     <- mean_vec("Rtanh",     3)
    agg$RCC       <- mean_vec("RCC",       3)

    # --- Decision Proportions (VARDIST) ---
    # Logic: Predictor Z significant (Non-Normal), Outcome Z non-significant (Normal) => Target (x->y)

    ago_pred_z <- sapply(valid_results, function(x) get_numeric(x$agostino$predictor$statistic[2]))
    ago_out_z  <- sapply(valid_results, function(x) get_numeric(x$agostino$outcome$statistic[2]))

    dec_ago <- ifelse(abs(ago_pred_z) >= crit_val & abs(ago_out_z) < crit_val, "Target",
                      ifelse(abs(ago_pred_z) < crit_val & abs(ago_out_z) >= crit_val, "Alternative", "Undecided"))
    decisions$dec_agost <- calc_props(dec_ago)

    ans_pred_z <- sapply(valid_results, function(x) get_numeric(x$anscombe$predictor$statistic[2]))
    ans_out_z  <- sapply(valid_results, function(x) get_numeric(x$anscombe$outcome$statistic[2]))

    dec_ans <- ifelse(abs(ans_pred_z) >= crit_val & abs(ans_out_z) < crit_val, "Target",
                      ifelse(abs(ans_pred_z) < crit_val & abs(ans_out_z) >= crit_val, "Alternative", "Undecided"))
    decisions$dec_anscom <- calc_props(dec_ans)

    calc_ci_decision <- function(key) {
      lowers <- sapply(valid_results, function(x) { v <- x[[key]]; if(length(v)>=2) v[2] else NA })
      uppers <- sapply(valid_results, function(x) { v <- x[[key]]; if(length(v)>=3) v[3] else NA })

      d <- ifelse(!is.na(lowers) & !is.na(uppers) & lowers > 0 & uppers > 0, "Target",
                  ifelse(!is.na(lowers) & !is.na(uppers) & lowers < 0 & uppers < 0, "Alternative", "Undecided"))
      calc_props(d)
    }

    decisions$dec_skewdiff  <- calc_ci_decision("skewdiff")
    decisions$dec_kurtdiff  <- calc_ci_decision("kurtdiff")
    decisions$dec_cor12diff <- calc_ci_decision("cor12diff")
    decisions$dec_cor13diff <- calc_ci_decision("cor13diff")
    decisions$dec_RHS       <- calc_ci_decision("RHS")
    decisions$dec_RCC       <- calc_ci_decision("RCC")
    decisions$dec_Rtanh     <- calc_ci_decision("Rtanh")

    class_label <- "dda_bagging_vardist"
  }

  output <- list(
    bagged_results = bagged_results,
    aggregated_stats = agg,
    decision_percentages = decisions,
    n_valid_iterations = n_valid,
    original_result = dda_result,
    parameters = list(
      function_name = function_name,
      object_type = object_type,
      iter = iter,
      successful_iterations = n_valid,
      failed_iterations = iter - n_valid
    )
  )
  if (!is.null(save_file)) {
    save(output, file = save_file)
    cat(paste("Results saved to:", save_file, "\n"))
  }
  class(output) <- c(class_label, "dda_bagging", class(output))
  return(output)
}
