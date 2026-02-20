#' Bootstrap Aggregated DDA Analysis
#'
#' @param dda_result Output from any DDA function (dda.indep, dda.vardist, dda.resdist)
#' @param iter Number of bootstrap iterations (default: 100)
#' @param progress Whether to show progress bar (default: TRUE)
#' @param save_file Optional file path to save results
#' @param alpha Significance level for decisions (default: 0.05)
#' @param data Raw data frame with ALL variables (outcome, predictor, covariates)
#' @return A list containing bootstrap and aggregated results
#' @export
dda_bagging <- function(
    dda_result,
    iter = 100,
    progress = TRUE,
    save_file = NULL,
    alpha = 0.05,
    data = NULL
) {

  # --- Helper: Safe Extraction ---
  get_numeric <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (is.numeric(x)) return(as.numeric(x[1]))
    if (is.list(x)) return(get_numeric(x[[1]]))
    return(NA_real_)
  }

  # --- Helper: Harmonic Mean P-values (Wilson 2019) ---
  harmonic_p <- function(pvec) {
    pvec <- as.numeric(pvec)
    pvec <- pvec[!is.na(pvec)]
    pvec[pvec <= 0] <- 1e-300 # Prevent harmonic mean crash on exact 0
    if (length(pvec) == 0) return(NA_real_)
    if (!requireNamespace("harmonicmeanp", quietly = TRUE)) return(mean(pvec))
    harmonicmeanp::p.hmp(pvec, L = length(pvec))
  }

  # --- Helper: Decision Percentages ---
  calc_props <- function(dec_vec) {
    levs <- c("Target", "Alternative", "Undecided")
    tab <- table(factor(dec_vec, levels = levs))
    return(tab / sum(tab))
  }

  # --- Validate Input ---
  if (!inherits(dda_result, c("dda.indep", "dda.resdist", "dda.vardist"))) {
    stop("Unsupported DDA object. Must be dda.indep, dda.resdist, or dda.vardist")
  }

  call_info <- dda_result$call_info
  if (is.null(data)) {
    stop("Please provide the original raw data (with all variables) via the 'data' argument.")
  }

  original_data <- data
  nobs <- nrow(original_data)
  dda_func <- get(call_info$function_name)
  obj_type <- class(dda_result)[1]

  var_names <- dda_result$var.names
  if (is.null(var_names) || length(var_names) != 2) {
    stop("Variable names not found in DDA result.")
  }

  y_name <- var_names[1]
  x_name <- var_names[2]

  # --- Extract Formula (to avoid Dummy Var issues) ---
  original_formula <- NULL
  if (!is.null(call_info$formula)) {
    if (inherits(call_info$formula, "formula")) {
      original_formula <- call_info$formula
    } else if (inherits(call_info$formula, "lm")) {
      original_formula <- formula(call_info$formula)
    } else if (!is.null(call_info$all_args$formula)) {


      if (inherits(call_info$all_args$formula, "formula")) {
        original_formula <- call_info$all_args$formula
      } else if (inherits(call_info$all_args$formula, "lm")) {
        original_formula <- formula(call_info$all_args$formula)
      } else {
        original_formula <- tryCatch(as.formula(call_info$all_args$formula), error = function(e) NULL)
      }
    }
  }

  if (is.null(original_formula)) {
    stop("Could not extract formula from DDA result.")
  }

  all_vars <- all.vars(original_formula)
  cov_names <- setdiff(all_vars, c(y_name, x_name))
  has_covariates <- length(cov_names) > 0

  # --- Bootstrap Execution Loop ---
  bagged_results <- vector("list", iter)
  if (progress) pb <- txtProgressBar(min = 0, max = iter, style = 3)

  for(i in 1:iter) {

    # 1. Resample Data
    boot_indices <- sample(1:nobs, nobs, replace = TRUE)
    datboot <- original_data[boot_indices, ]

    # 2. Residualize Variables
    if (has_covariates) {
      cov_formula_y <- as.formula(paste(y_name, "~", paste(cov_names, collapse = " + ")))
      cov_formula_x <- as.formula(paste(x_name, "~", paste(cov_names, collapse = " + ")))

      tryCatch({
        ry <- as.vector(scale(resid(lm(cov_formula_y, data = datboot))))
        rx <- as.vector(scale(resid(lm(cov_formula_x, data = datboot))))
      }, error = function(e) {
        stop(paste("Covariate residualization failed:", e$message))
      })
    } else {
      ry <- as.vector(scale(datboot[[y_name]]))
      rx <- as.vector(scale(datboot[[x_name]]))
    }

    # 3. boot_args obj
    boot_args <- call_info$all_args
    if (length(boot_args) > 0 && names(boot_args)[1] == "") {
      boot_args[[1]] <- NULL # remove the original lm object
    }

    boot_args$formula <- as.formula(paste(y_name, "~", x_name))
    boot_args$pred <- x_name

    boot_processed <- data.frame(rx, ry)
    names(boot_processed) <- c(x_name, y_name)
    boot_args$data <- boot_processed

    # 4. Call DDA Module
    bagged_results[[i]] <- tryCatch(
      do.call(dda_func, boot_args),
      error = function(e) {
        warning(paste("Bootstrap iteration", i, "failed:", e$message))
        return(NA)
      }
    )

    if (progress) setTxtProgressBar(pb, i)
  }

  if (progress) close(pb)

  # --- Filter Valid Results ---
  valid_res <- bagged_results[!sapply(bagged_results, function(x) all(is.na(x)) || is.null(x))]
  n_valid <- length(valid_res)

  if (n_valid == 0) {
    stop("No valid bootstrap iterations. Check your data and DDA function arguments.")
  }

  agg <- list()
  decs <- list()
  crit <- qnorm(1 - alpha/2)

  ## ============================================================================
  ## INDEP Block
  ## ============================================================================
  if (obj_type == "dda.indep") {
    agg$var.names <- var_names

    agg$hsic_yx_stat <- mean(sapply(valid_res, function(x) get_numeric(x$hsic.yx$statistic)), na.rm=TRUE)
    agg$hsic_xy_stat <- mean(sapply(valid_res, function(x) get_numeric(x$hsic.xy$statistic)), na.rm=TRUE)
    agg$hsic_yx_pval <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$hsic.yx$p.value)))
    agg$hsic_xy_pval <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$hsic.xy$p.value)))

    h_p_yx <- sapply(valid_res, function(x) get_numeric(x$hsic.yx$p.value))
    h_p_xy <- sapply(valid_res, function(x) get_numeric(x$hsic.xy$p.value))
    decs$hsic <- calc_props(ifelse(h_p_yx > alpha & h_p_xy <= alpha, "Target",
                                   ifelse(h_p_xy > alpha & h_p_yx <= alpha, "Alternative", "Undecided")))

    if (!is.null(valid_res[[1]]$dcor.yx)) {
      agg$dcor_yx_stat <- mean(sapply(valid_res, function(x) get_numeric(x$dcor.yx$statistic)), na.rm=TRUE)
      agg$dcor_xy_stat <- mean(sapply(valid_res, function(x) get_numeric(x$dcor.xy$statistic)), na.rm=TRUE)
      agg$dcor_yx_pval <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$dcor.yx$p.value)))
      agg$dcor_xy_pval <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$dcor.xy$p.value)))

      d_p_yx <- sapply(valid_res, function(x) get_numeric(x$dcor.yx$p.value))
      d_p_xy <- sapply(valid_res, function(x) get_numeric(x$dcor.xy$p.value))
      decs$dcor <- calc_props(ifelse(d_p_yx > alpha & d_p_xy <= alpha, "Target",
                                     ifelse(d_p_xy > alpha & d_p_yx <= alpha, "Alternative", "Undecided")))
    }

    if (!is.null(valid_res[[1]]$breusch_pagan)) {
      bp_yx_stat <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[1]]$statistic))
      bp_yx_df <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[1]]$parameter))
      bp_yx_p <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[1]]$p.value))

      rbp_yx_stat <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[2]]$statistic))
      rbp_yx_df <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[2]]$parameter))
      rbp_yx_p <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[2]]$p.value))

      bp_xy_stat <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[3]]$statistic))
      bp_xy_df <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[3]]$parameter))
      bp_xy_p <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[3]]$p.value))

      rbp_xy_stat <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[4]]$statistic))
      rbp_xy_df <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[4]]$parameter))
      rbp_xy_p <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[4]]$p.value))

      agg$breusch_pagan <- list(
        list(statistic = mean(bp_yx_stat, na.rm=TRUE), parameter = mean(bp_yx_df, na.rm=TRUE), p.value = harmonic_p(bp_yx_p)),
        list(statistic = mean(rbp_yx_stat, na.rm=TRUE), parameter = mean(rbp_yx_df, na.rm=TRUE), p.value = harmonic_p(rbp_yx_p)),
        list(statistic = mean(bp_xy_stat, na.rm=TRUE), parameter = mean(bp_xy_df, na.rm=TRUE), p.value = harmonic_p(bp_xy_p)),
        list(statistic = mean(rbp_xy_stat, na.rm=TRUE), parameter = mean(rbp_xy_df, na.rm=TRUE), p.value = harmonic_p(rbp_xy_p))
      )
      decs$dec_bp <- calc_props(ifelse(rbp_yx_p > alpha & rbp_xy_p <= alpha, "Target",
                                       ifelse(rbp_yx_p <= alpha & rbp_xy_p > alpha, "Alternative", "Undecided")))
    }

    if (!is.null(valid_res[[1]]$nlcor.yx)) {
      nlcor_yx_t1 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.yx$t1)))
      nlcor_yx_t2 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.yx$t2)))
      nlcor_yx_t3 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.yx$t3)))
      nlcor_xy_t1 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.xy$t1)))
      nlcor_xy_t2 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.xy$t2)))
      nlcor_xy_t3 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.xy$t3)))

      agg$nlcor.yx <- list(
        t1 = c(mean(nlcor_yx_t1[,1], na.rm=TRUE), mean(nlcor_yx_t1[,2], na.rm=TRUE), mean(nlcor_yx_t1[,3], na.rm=TRUE), harmonic_p(nlcor_yx_t1[,4])),
        t2 = c(mean(nlcor_yx_t2[,1], na.rm=TRUE), mean(nlcor_yx_t2[,2], na.rm=TRUE), mean(nlcor_yx_t2[,3], na.rm=TRUE), harmonic_p(nlcor_yx_t2[,4])),
        t3 = c(mean(nlcor_yx_t3[,1], na.rm=TRUE), mean(nlcor_yx_t3[,2], na.rm=TRUE), mean(nlcor_yx_t3[,3], na.rm=TRUE), harmonic_p(nlcor_yx_t3[,4])),
        func = valid_res[[1]]$nlcor.yx$func
      )

      agg$nlcor.xy <- list(
        t1 = c(mean(nlcor_xy_t1[,1], na.rm=TRUE), mean(nlcor_xy_t1[,2], na.rm=TRUE), mean(nlcor_xy_t1[,3], na.rm=TRUE), harmonic_p(nlcor_xy_t1[,4])),
        t2 = c(mean(nlcor_xy_t2[,1], na.rm=TRUE), mean(nlcor_xy_t2[,2], na.rm=TRUE), mean(nlcor_xy_t2[,3], na.rm=TRUE), harmonic_p(nlcor_xy_t2[,4])),
        t3 = c(mean(nlcor_xy_t3[,1], na.rm=TRUE), mean(nlcor_xy_t3[,2], na.rm=TRUE), mean(nlcor_xy_t3[,3], na.rm=TRUE), harmonic_p(nlcor_xy_t3[,4])),
        func = valid_res[[1]]$nlcor.xy$func
      )

      nlcor_yx_min <- apply(cbind(nlcor_yx_t1[,4], nlcor_yx_t2[,4], nlcor_yx_t3[,4]), 1, min)
      nlcor_xy_min <- apply(cbind(nlcor_xy_t1[,4], nlcor_xy_t2[,4], nlcor_xy_t3[,4]), 1, min)
      decs$dec_nl.min <- calc_props(ifelse(nlcor_yx_min > alpha & nlcor_xy_min <= alpha, "Target",
                                           ifelse(nlcor_yx_min <= alpha & nlcor_xy_min > alpha, "Alternative", "Undecided")))
    }

    if (!is.null(valid_res[[1]]$out.diff)) {
      diff_arr <- simplify2array(lapply(valid_res, function(x) as.matrix(x$out.diff)))
      agg$diff_matrix <- apply(diff_arr, c(1,2), mean, na.rm=TRUE)

      calc_diff <- function(row) {
        nc <- ncol(valid_res[[1]]$out.diff)
        low_idx <- nc - 1
        upp_idx <- nc
        low <- sapply(valid_res, function(x) x$out.diff[row, low_idx])
        upp <- sapply(valid_res, function(x) x$out.diff[row, upp_idx])
        calc_props(ifelse(low > 0 & upp > 0, "Target",
                          ifelse(low < 0 & upp < 0, "Alternative", "Undecided")))
      }
      decs$diff_hsic <- calc_diff(1)
      decs$diff_dcor <- calc_diff(2)
      decs$diff_mi <- calc_diff(3)
    }
  }

  ## ============================================================================
  ## RESDIST Block
  ## ============================================================================
  if (obj_type == "dda.resdist") {
    agg$var.names <- var_names

    z_tar <- sapply(valid_res, function(x) get_numeric(x$agostino$target$statistic[2]))
    z_alt <- sapply(valid_res, function(x) get_numeric(x$agostino$alternative$statistic[2]))
    decs$dec_agost <- calc_props(ifelse(abs(z_alt) >= crit & abs(z_tar) < crit, "Target",
                                        ifelse(abs(z_tar) >= crit & abs(z_alt) < crit, "Alternative", "Undecided")))

    z_tar_kurt <- sapply(valid_res, function(x) get_numeric(x$anscombe$target$statistic[2]))
    z_alt_kurt <- sapply(valid_res, function(x) get_numeric(x$anscombe$alternative$statistic[2]))
    decs$dec_anscom <- calc_props(ifelse(abs(z_alt_kurt) >= crit & abs(z_tar_kurt) < crit, "Target",
                                         ifelse(abs(z_tar_kurt) >= crit & abs(z_alt_kurt) < crit, "Alternative", "Undecided")))

    agg$agostino.target.statistic <- mean(sapply(valid_res, function(x) get_numeric(x$agostino$target$statistic[1])), na.rm=TRUE)
    agg$agostino.target.z <- mean(z_tar, na.rm=TRUE)
    agg$agostino.target.p.value <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$agostino$target$p.value)))

    agg$agostino.alternative.statistic <- mean(sapply(valid_res, function(x) get_numeric(x$agostino$alternative$statistic[1])), na.rm=TRUE)
    agg$agostino.alternative.z <- mean(z_alt, na.rm=TRUE)
    agg$agostino.alternative.p.value <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$agostino$alternative$p.value)))

    agg$anscombe.target.statistic <- mean(sapply(valid_res, function(x) get_numeric(x$anscombe$target$statistic[1])), na.rm=TRUE)
    agg$anscombe.target.z <- mean(z_tar_kurt, na.rm=TRUE)
    agg$anscombe.target.p.value <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$anscombe$target$p.value)))

    agg$anscombe.alternative.statistic <- mean(sapply(valid_res, function(x) get_numeric(x$anscombe$alternative$statistic[1])), na.rm=TRUE)
    agg$anscombe.alternative.z <- mean(z_alt_kurt, na.rm=TRUE)
    agg$anscombe.alternative.p.value <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$anscombe$alternative$p.value)))

    # Joint moments - Dynamic Matrix Indexing
    for(k in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS3", "RCC", "RHS4")) {
      if (is.null(valid_res[[1]][[k]])) next
      mat <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x[[k]])))
      agg[[k]] <- colMeans(mat, na.rm=TRUE)

      nc <- ncol(mat)
      low_idx <- nc - 1
      upp_idx <- nc
      decs[[paste0("dec_", k)]] <- calc_props(ifelse(mat[,low_idx] > 0 & mat[,upp_idx] > 0, "Target",
                                                     ifelse(mat[,low_idx] < 0 & mat[,upp_idx] < 0, "Alternative", "Undecided")))
    }
  }

  ## ============================================================================
  ## VARDIST Block
  ## ============================================================================
  if (obj_type == "dda.vardist") {
    agg$var.names <- var_names

    z_pre <- sapply(valid_res, function(x) get_numeric(x$agostino$predictor$statistic[2]))
    z_out <- sapply(valid_res, function(x) get_numeric(x$agostino$outcome$statistic[2]))
    decs$dec_agost <- calc_props(ifelse(abs(z_pre) >= crit & abs(z_out) < crit, "Target",
                                        ifelse(abs(z_pre) < crit & abs(z_out) >= crit, "Alternative", "Undecided")))

    z_pre_kurt <- sapply(valid_res, function(x) get_numeric(x$anscombe$predictor$statistic[2]))
    z_out_kurt <- sapply(valid_res, function(x) get_numeric(x$anscombe$outcome$statistic[2]))
    decs$dec_anscom <- calc_props(ifelse(abs(z_pre_kurt) >= crit & abs(z_out_kurt) < crit, "Target",
                                         ifelse(abs(z_pre_kurt) < crit & abs(z_out_kurt) >= crit, "Alternative", "Undecided")))

    agg$agostino.predictor.statistic.skew <- mean(sapply(valid_res, function(x) get_numeric(x$agostino$predictor$statistic[1])), na.rm=TRUE)
    agg$agostino.predictor.statistic.z <- mean(z_pre, na.rm=TRUE)
    agg$agostino.predictor.p.value <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$agostino$predictor$p.value)))

    agg$agostino.outcome.statistic.skew <- mean(sapply(valid_res, function(x) get_numeric(x$agostino$outcome$statistic[1])), na.rm=TRUE)
    agg$agostino.outcome.statistic.z <- mean(z_out, na.rm=TRUE)
    agg$agostino.outcome.p.value <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$agostino$outcome$p.value)))

    agg$anscombe.predictor.statistic.kurt <- mean(sapply(valid_res, function(x) get_numeric(x$anscombe$predictor$statistic[1])), na.rm=TRUE)
    agg$anscombe.predictor.statistic.z <- mean(z_pre_kurt, na.rm=TRUE)
    agg$anscombe.predictor.p.value <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$anscombe$predictor$p.value)))

    agg$anscombe.outcome.statistic.kurt <- mean(sapply(valid_res, function(x) get_numeric(x$anscombe$outcome$statistic[1])), na.rm=TRUE)
    agg$anscombe.outcome.statistic.z <- mean(z_out_kurt, na.rm=TRUE)
    agg$anscombe.outcome.p.value <- harmonic_p(sapply(valid_res, function(x) get_numeric(x$anscombe$outcome$p.value)))

    # Joint moments - Use Dynamic Matrix Indexing
    for(k in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS", "RCC", "Rtanh")) {
      if (is.null(valid_res[[1]][[k]])) next
      mat <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x[[k]])))
      agg[[k]] <- colMeans(mat, na.rm=TRUE)

      nc <- ncol(mat)
      low_idx <- nc - 1
      upp_idx <- nc
      decs[[paste0("dec_", k)]] <- calc_props(ifelse(mat[,low_idx] > 0 & mat[,upp_idx] > 0, "Target",
                                                     ifelse(mat[,low_idx] < 0 & mat[,upp_idx] < 0, "Alternative", "Undecided")))
    }
  }

  if (!is.null(save_file)) {
    saveRDS(list(bagged_results = bagged_results,
                 aggregated_stats = agg,
                 decision_percentages = decs,
                 n_valid_iterations = n_valid),
            file = save_file)
  }

  out <- list(bagged_results = bagged_results,
              aggregated_stats = agg,
              decision_percentages = decs,
              n_valid_iterations = n_valid)
  class(out) <- c(paste0("dda_bagging_", gsub("dda.", "", obj_type)), "dda_bagging")
  return(out)
}
