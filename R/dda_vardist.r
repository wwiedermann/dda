#' @title Direction Dependence Analysis: Variable Distributions
#' @description \code{dda.vardist} evaluates patterns of asymmetry of variable
#'              distributions for causally competing models
#'              (\code{y ~ x} vs. \code{x ~ y}).
#' @name dda.vardist
#'
#' @param formula     Symbolic formula of the model to be tested or a \code{lm}object.
#' @param pred        Variable name of the predictor which serves as the outcome in the alternative model.
#' @param data        An optional data frame containing the variables in the
#'                    model (by default variables are taken from the environment
#'                    which \code{dda.vardist} is called from).
#' @param B           Number of bootstrap samples.
#' @param boot.type   A character indicating the type of bootstrap confidence intervals. Must be one of c("perc","bca"). Default is "perc".
#' @param conf.level  Confidence level for bootstrap confidence intervals.
#' @param ...         Additional arguments to be passed to the function.
#'
#' @returns  An object of class \code{dda.vardist} containing the results of DDA tests
#'           of asymmetry patterns of variable distributions.
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' x <- rchisq(n, df = 4) - 4
#' e <- rchisq(n, df = 3) - 3
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#' result <- dda.vardist(y ~ x, pred = "x", data = d, B = 50)
#'
#' @references Wiedermann, W., & von Eye, A. (2025). Direction Dependence Analysis: Foundations and Statistical Methods. Cambridge University Press.
#' @seealso \code{\link{cdda.vardist}} for a conditional version.
#' @export
#' @rdname dda.vardist
dda.vardist <- function(
    formula,
    pred = NULL,
    data = list(),
    B = 200,
    boot.type = "perc",
    conf.level = 0.95,
    ...
) {

  # --- helper functions for bootstrap CIs

  mysd <- function(x) sqrt(sum((x - mean(x))^2) / length(x))

  cor.ij <- function(x, y, i = 1, j = 1) {
    n  <- length(x)
    mx <- mean(x)
    my <- mean(y)
    Cov <- sum((x - mx)^i * (y - my)^j) / n
    Cov / (mysd(x)^i * mysd(y)^j)
  }

  boot.stat <- function(dat, g) {
    dat <- dat[g, ]
    x <- dat[, 1]  # "purified" predictor
    y <- dat[, 2]  # "purified" outcome

    x <- as.vector(scale(x))
    y <- as.vector(scale(y))

    skew.diff <- (moments::skewness(x)^2) - (moments::skewness(y)^2)
    kurt.diff <- (moments::kurtosis(x) - 3)^2 - (moments::kurtosis(y) - 3)^2
    cor21.diff <- (cor.ij(x, y, i = 2, j = 1)^2) - (cor.ij(x, y, i = 1, j = 2)^2)
    cor13.diff <- (cor.ij(x, y, i = 3, j = 1)^2 - cor.ij(x, y, i = 1, j = 3)^2) * sign(moments::kurtosis(x) - 3)

    Rtanh <- stats::cor(x, y) * mean(x * tanh(y) - tanh(x) * y)

    Cxy <- mean(x^3 * y) - 3 * stats::cor(x, y) * stats::var(x)
    Cyx <- mean(x * y^3) - 3 * stats::cor(x, y) * stats::var(y)
    RCC <- (Cxy + Cyx) * (Cxy - Cyx)

    xx <- sign(moments::skewness(x)) * x
    yy <- sign(moments::skewness(y)) * y
    RHS <- stats::cor(xx, yy) * mean((xx^2 * yy) - (xx * yy^2))

    out <- c(skew.diff, kurt.diff, cor21.diff, cor13.diff, RHS, RCC, Rtanh)
    names(out) <- c("skew.diff", "kurt.diff", "cor21.diff", "cor13.diff", "RHS", "RCC", "Rtanh")
    out
  }

  extract_limits <- function(ci_obj, type, t_vals, conf.level) {
    # Returns c(lower, upper) given a boot.ci object and desired type.
    # Falls back to simple percentiles if structure is missing.
    find_last_two <- function(mat) {
      mat <- as.matrix(mat)
      c(as.numeric(mat[1, ncol(mat) - 1]), as.numeric(mat[1, ncol(mat)]))
    }
    if (!inherits(ci_obj, "bootci")) {
      # direct percentile fallback
      alpha <- (1 - conf.level) / 2
      qs <- stats::quantile(t_vals, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
      return(c(qs[1], qs[2]))
    }
    if (type == "bca" && !is.null(ci_obj$bca)) {
      return(find_last_two(ci_obj$bca))
    }
    if (type == "perc" && !is.null(ci_obj$percent)) {
      return(find_last_two(ci_obj$percent))
    }
    # ultimate fallback
    alpha <- (1 - conf.level) / 2
    qs <- stats::quantile(t_vals, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
    c(qs[1], qs[2])
  }

  if (is.null(pred)) stop("Tentative predictor is missing.")
  if (B < 0) stop("Number of resamples 'B' must be non-negative.")
  if (conf.level <= 0 || conf.level >= 1) stop("'conf.level' must be between 0 and 1 (exclusive).")
  if (!boot.type %in% c("bca", "perc")) stop("Unknown argument in boot.type.")

  # --- prepare outcome, predictor, and model matrix for covariates

  if (!inherits(formula, "formula")) {
    X <- if (is.matrix(formula$x)) formula$x else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(formula$y)) formula$y else model.response(model.frame(formula))

    delete.pred <- which(colnames(X) == pred) # get position of tentative predictor
    if (length(delete.pred) == 0) stop("Specified predictor not found in the target model.")

    x <- X[, delete.pred]   # tentative predictor
    X <- X[, -delete.pred]  # model matrix with covariates
    if (!is.matrix(X)) X <- as.matrix(X)
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)   # tentative outcome
    X <- model.matrix(formula, data = data)

    delete.pred <- which(colnames(X) == pred) # get position of tentative predictor
    if (length(delete.pred) == 0) stop("Specified predictor not found in the target model.")

    x <- X[, delete.pred]   # tentative predictor
    X <- X[, -delete.pred]  # model matrix with covariates
    if (!is.matrix(X)) X <- as.matrix(X)
  }

  ry <- stats::lm.fit(X, y)$residuals
  rx <- stats::lm.fit(X, x)$residuals

  ry <- as.vector(scale(ry))
  rx <- as.vector(scale(rx))

  dat <- data.frame(predictor = rx, outcome = ry)

  # --- run separate normality tests

  agostino.out <- apply(dat, 2, moments::agostino.test)
  agostino.out <- lapply(agostino.out, unclass)
  agostino.out$predictor[3:5] <- NULL
  agostino.out$outcome[3:5] <- NULL

  anscombe.out <- apply(dat, 2, moments::anscombe.test)
  anscombe.out <- lapply(anscombe.out, unclass)
  anscombe.out$predictor[3:5] <- NULL
  anscombe.out$outcome[3:5] <- NULL

  output <- list(agostino = agostino.out, anscombe = anscombe.out)

  # Convert kurtosis to excess-kurtosis for printing consistency
  output$anscombe$outcome$statistic[1]   <- output$anscombe$outcome$statistic[1] - 3
  output$anscombe$predictor$statistic[1] <- output$anscombe$predictor$statistic[1] - 3

  # --- point estimates (always compute once)
  point <- boot.stat(dat, seq_len(nrow(dat)))

  # sign mismatch warning flag for excess-kurtosis
  kx <- moments::kurtosis(rx) - 3
  ky <- moments::kurtosis(ry) - 3
  boot.warning <- is.finite(kx) && is.finite(ky) && (sign(kx) != sign(ky)) && (kx != 0) && (ky != 0)

  # --- bootstrap confidence intervals
  boot_args <- NULL
  ci_rows <- vector("list", 7)

  if (B > 0) {
    suppressWarnings(boot.res <- boot::boot(dat, boot.stat, R = B))

    # Try BCa if requested, else go straight to percentile
    boot_out <- NULL
    effective_type <- boot.type

    if (boot.type == "bca") {
      boot_out <- tryCatch(
        {
          suppressWarnings(lapply(seq_len(7), function(i, res) {
            boot::boot.ci(res, conf = conf.level, type = "bca", t0 = res$t0[i], t = res$t[, i])
          }, res = boot.res))
        },
        error = function(e) {
          warning("Acceleration constant cannot be calculated. Falling back to percentile bootstrap method. Consider increasing the number of resamples (B) for more stable results.", call. = FALSE)
          NULL
        }
      )
      if (is.null(boot_out)) {
        effective_type <- "perc"
      }
    }

    if (is.null(boot_out) || effective_type == "perc") {
      suppressWarnings(boot_out <- lapply(seq_len(7), function(i, res) {
        boot::boot.ci(res, conf = conf.level, type = "perc", t0 = res$t0[i], t = res$t[, i])
      }, res = boot.res))
    }

    # Build CI rows (estimate, lower, upper) for each statistic
    for (i in seq_len(7)) {
      est <- as.numeric(boot.res$t0[i])
      lims <- extract_limits(boot_out[[i]], if (effective_type == "bca") "bca" else "perc", boot.res$t[, i], conf.level)
      ci_rows[[i]] <- c(est, lims[1], lims[2])
    }

    boot_args <- c(effective_type, conf.level, B)
  }

  # Map rows to the names the print method expects
  # Index order: 1=skew.diff, 2=kurt.diff, 3=cor21.diff, 4=cor13.diff, 5=RHS, 6=RCC, 7=Rtanh
  if (B > 0) {
    output$skewdiff <- setNames(ci_rows[[1]], c("estimate", "lower", "upper"))
    output$kurtdiff <- setNames(ci_rows[[2]], c("estimate", "lower", "upper"))
    # The print method expects 'cor12diff' as the name, but it's Cor^2[2,1] - Cor^2[1,2]
    output$cor12diff <- setNames(ci_rows[[3]], c("estimate", "lower", "upper"))
    output$cor13diff <- setNames(ci_rows[[4]], c("estimate", "lower", "upper"))
    output$RHS      <- setNames(ci_rows[[5]], c("estimate", "lower", "upper"))
    output$RCC      <- setNames(ci_rows[[6]], c("estimate", "lower", "upper"))
    output$Rtanh    <- setNames(ci_rows[[7]], c("estimate", "lower", "upper"))
  } else {
    # No bootstrap: provide point estimates only so print() can still show them
    output$skewdiff <- unname(point[1])
    output$kurtdiff <- unname(point[2])
    output$cor12diff <- unname(point[3])
    output$cor13diff <- unname(point[4])
    output$RHS <- unname(point[5])
    output$RCC <- unname(point[6])
    output$Rtanh <- unname(point[7])
  }

  # Attach metadata and flags
  response.name <- all.vars(stats::formula(formula))[1]
  output$var.names <- c(response.name, pred)
  output$boot.args <- boot_args
  output$boot.warning <- isTRUE(boot.warning)

  class(output) <- "dda.vardist"
  return(output)
}


#' @name print.dda.vardist
#' @title Print Method for \code{dda.vardist} Objects
#'
#' @description \code{print} returns DDA test statistics associated with \code{dda.vardist} objects.
#' @param x     An object of class \code{dda.vardist} when using \code{print}.
#' @param ...   Additional arguments to be passed to the function.
#'
#' @examples
#' # print(result)
#' @returns An object of class \code{dda.vardist}.
#'
#' @export
#' @rdname dda.vardist
#' @method print dda.vardist
print.dda.vardist <- function(x, ...) {
  varnames <- x$var.names

  cat("\n")
  cat("DIRECTION DEPENDENCE ANALYSIS: Variable Distributions", "\n", "\n")
  cat("Skewness and kurtosis tests:", "\n")

  sigtests <- rbind(
    c(x[[1]]$outcome$statistic, x[[1]]$outcome$p.value, x[[1]]$predictor$statistic, x[[1]]$predictor$p.value),
    c(x[[2]]$outcome$statistic, x[[2]]$outcome$p.value, x[[2]]$predictor$statistic, x[[2]]$predictor$p.value)
  )
  sigtests <- round(sigtests, 4)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")
  print.default(format(sigtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

  # No bootstrap available: print point estimates
  if (is.null(x$boot.args)) {
    cat("\n")
    cat("Higher moment differences and Likelihood Ratio approximations (point estimates):", "\n")

    # Gather estimates (first element if vectors)
    get_est <- function(v) if (length(v) >= 1) as.numeric(v[1]) else NA_real_

    all_ests <- c(
      get_est(x$skewdiff),
      get_est(x$kurtdiff),
      get_est(x$cor12diff),
      get_est(x$cor13diff),
      get_est(x$RHS),
      get_est(x$Rtanh),
      get_est(x$RCC)
    )
    all_ests <- matrix(as.numeric(all_ests), ncol = 1)
    rownames(all_ests) <- c(
      "Skewness Diff",
      "Kurtosis Diff",
      "Cor^2[2,1] - Cor^2[1,2]",
      "Cor^2[3,1] - Cor^2[1,3]",
      "Hyvarinen-Smith (co-skewness)",
      "Hyvarinen-Smith (tanh)",
      "Chen-Chan (co-kurtosis)"
    )
    colnames(all_ests) <- c("estimate")

    all_ests <- round(all_ests, 4)
    print.default(format(all_ests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")
    return(invisible(x))
  }

  # With bootstrap: print CI sections
  ci.level <- as.numeric(x$boot.args[2]) * 100
  cat("\n")

  if (x$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for higher moment differences:", "\n", sep = "")
  if (x$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for higher moment differences:", "\n", sep = "")

  citests <- rbind(x$skewdiff, x$kurtdiff)
  citests <- round(citests, 4)
  rownames(citests) <- c("Skewness", "Kurtosis")
  colnames(citests) <- c("diff", "lower", "upper")
  print.default(format(citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
  cat("\n")

  if (x$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for differences in higher-order correlations:", "\n", sep = "")
  if (x$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for differences in higher-order correlations:", "\n", sep = "")

  hoctests <- rbind(x$cor12diff, x$cor13diff)
  hoctests <- round(hoctests, 4)
  rownames(hoctests) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]")
  colnames(hoctests) <- c("estimate", "lower", "upper")
  print.default(format(hoctests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
  cat("\n")

  if (x$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for Likelihood Ratio approximations:", "\n", sep = "")
  if (x$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for Likelihood Ratio approximations:", "\n", sep = "")

  LRtests <- rbind(x$RHS, x$Rtanh, x$RCC)
  LRtests <- round(LRtests, 4)
  rownames(LRtests) <- c("Hyvarinen-Smith (co-skewness)", "Hyvarinen-Smith (tanh)", "Chen-Chan (co-kurtosis)")
  colnames(LRtests) <- c("estimate", "lower", "upper")
  print.default(format(LRtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

  cat("\n")
  cat(paste("Number of resamples:", x$boot.args[3]))
  cat("\n")
  cat("---")
  cat("\n")
  cat(paste("Note: (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests the model", varnames[2], "->", varnames[1], sep = " "))
  cat("\n")
  if (isTRUE(x$boot.warning)) {
    cat("Warning: Excess-kurtosis values of", varnames[2], "and", varnames[1], "have unequal signs", "\n",
        "        Cor^2[3,1] - Cor^2[1,3] should also be computed for the model", varnames[1], "->", varnames[2], "\n")
  }

  invisible(x)
}
