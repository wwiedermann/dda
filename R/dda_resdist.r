#' @title Direction Dependence Analysis: Residual Distributions
#' @description \code{dda.resdist} evaluates patterns of asymmetry of error
#' distributions of causally competing models (\code{y ~ x} vs. \code{x ~ y}).
#' @name dda.resdist
#'
#' @param formula Symbolic formula of the target model to be tested or a \code{lm} object.
#' @param pred Variable name of the predictor which serves as the outcome in the alternative model.
#' @param data Optional data frame containing the variables in the model.
#' @param B Number of bootstrap samples.
#' @param boot.type Bootstrap CI type, one of \code{c("perc", "bca")}; default "perc".
#' @param prob.trans Logical, probability integral transformation before skew/kurtosis tests.
#' @param conf.level Confidence level for bootstrap CIs.
#' @param robust Use Siegel non-parametric estimators for residual extraction.
#' @returns An object of class \code{dda.resdist}
#' @export
dda.resdist <- function(formula,
                        pred = NULL,
                        data = list(),
                        B = 200,
                        boot.type = "perc",
                        prob.trans = FALSE,
                        conf.level = 0.95,
                        robust = FALSE
) {
  # --- Helper functions ---
  mysd <- function(x) sqrt(sum((x - mean(x))^2) / length(x))
  cor.ij <- function(x, y, i = 1, j = 1) {
    n <- length(x)
    mx <- mean(x)
    my <- mean(y)
    Cov <- sum((x - mx)^i * (y - my)^j) / n
    Cov / (mysd(x)^i * mysd(y)^j)
  }

  # Probability integral transformation (safe for ties, NAs, edge cases)
  prob.int <- function(x, y) {
    # Remove NAs just in case
    x <- as.numeric(x)
    y <- as.numeric(y)
    # Empirical CDF for y, including ties
    fval <- ecdf(y)(y)
    # Map to quantiles via interpolation, always monotonic
    quant_y <- approxfun(fval, y, rule = 2)
    # Empirical CDF for x
    cum_x <- ecdf(x)
    # Map x to y's scale using y's ECDF and quantiles
    x_trans <- quant_y(cum_x(x))
    # If any NA (e.g., outside quantile range), set to min(y)
    x_trans[is.na(x_trans)] <- min(y, na.rm = TRUE)
    list(x = x_trans, y = y)
  }

  # --- Bootstrap statistic function ---
  boot.diff <- function(dat, g, prob.trans, robust) {
    # Always try-catch!
    out <- tryCatch({
      dat <- dat[g, ]
      x <- as.vector(scale(dat[, 1]))
      y <- as.vector(scale(dat[, 2]))
      # Degenerate check
      if(anyNA(x) || anyNA(y) || length(unique(x)) < 3 || length(unique(y)) < 3) return(rep(NA, 7))
      if (isTRUE(prob.trans)) {
        xboot <- dat[, 3]
        yboot <- dat[, 4]
        # Degenerate check for transformation inputs
        if(anyNA(xboot) || anyNA(yboot) || length(unique(xboot)) < 3 || length(unique(yboot)) < 3) return(rep(NA, 7))
        trans.boot <- prob.int(xboot, yboot)
        xboot.trans <- trans.boot$x
        yboot.trans <- trans.boot$y
        # Degenerate check after transformation
        if(anyNA(xboot.trans) || anyNA(yboot.trans) || length(unique(xboot.trans)) < 3 || length(unique(yboot.trans)) < 3) return(rep(NA, 7))
        # Fit models
        tar.boot.trans <- tryCatch(
          if (robust) mblm::mblm(yboot.trans ~ xboot.trans) else lm(yboot.trans ~ xboot.trans),
          error = function(e) NULL)
        alt.boot.trans <- tryCatch(
          if (robust) mblm::mblm(xboot.trans ~ yboot.trans) else lm(xboot.trans ~ yboot.trans),
          error = function(e) NULL)
        # If model fitting fails or output is degenerate
        if (is.null(tar.boot.trans) || is.null(alt.boot.trans)) return(rep(NA, 7))
        xtrans <- as.vector(scale(resid(alt.boot.trans)))
        ytrans <- as.vector(scale(resid(tar.boot.trans)))
        if(anyNA(xtrans) || anyNA(ytrans) || length(unique(xtrans)) < 3 || length(unique(ytrans)) < 3) return(rep(NA, 7))
        skew.diff <- (moments::skewness(xtrans)^2) - (moments::skewness(ytrans)^2)
        kurt.diff <- (moments::kurtosis(xtrans) - 3)^2 - (moments::kurtosis(ytrans) - 3)^2
        x <- xtrans
        y <- ytrans
      } else {
        # Already checked for degenerate x/y above
        skew.diff <- (moments::skewness(x)^2) - (moments::skewness(y)^2)
        kurt.diff <- (moments::kurtosis(x) - 3)^2 - (moments::kurtosis(y) - 3)^2
      }
      cor12.diff <- (cor.ij(y, x, 2, 1)^2) - (cor.ij(y, x, 1, 2)^2)
      cor13.diff <- ((cor.ij(y, x, 3, 1)^2) - (cor.ij(y, x, 1, 3)^2)) * sign(moments::kurtosis(y) - 3)
      xx <- sign(moments::skewness(x)) * x
      yy <- sign(moments::skewness(y)) * y
      RHS3 <- cor(xx, yy) * mean((yy^2 * xx) - (yy * xx^2))
      RHS4 <- cor(x, y) * mean((y^3 * x) - (y * x^3)) * sign(moments::kurtosis(y) - 3)
      C1 <- mean(y^3 * x) - 3 * cor(x, y) * var(y)
      C2 <- mean(y * x^3) - 3 * cor(x, y) * var(x)
      RCC <- (C1 + C2) * (C1 - C2)
      c(skew.diff, kurt.diff, cor12.diff, cor13.diff, RHS3, RCC, RHS4)
    }, error = function(e) rep(NA, 7))
    if (!is.numeric(out) || length(out) != 7) out <- rep(NA, 7)
    out
  }
  # --- Normality tests ---
  skew.diff.test <- function(x, y){
    agostino.zvalue <- function(x){
      n  <- length(x)
      s3 <- (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
      y <- s3 * sqrt((n + 1) * (n + 3)/(6 * (n - 2)))
      b2 <- 3 * (n^2 + 27*n - 70) * (n + 1) * (n + 3) / ((n - 2) * (n + 5) * (n + 7) * (n + 9))
      w <- -1 + sqrt(2 * (b2 - 1))
      d <- 1/sqrt(log(sqrt(w)))
      a <- sqrt(2/(w - 1))
      z <- d * log(y/a + sqrt((y/a)^2 + 1))
      return(z)
    }
    zval <- (agostino.zvalue(x) - agostino.zvalue(y))/sqrt(2 - 2*cor(x,y)^3)
    pval <- (1 - pnorm(abs(zval))) * 2
    return(list(z.value = zval, p.value = pval))
  }
  kurt.diff.test <- function(x, y){
    anscombe.zvalue <- function(x){
      n   <- length(x)
      b   <- n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
      eb2 <- 3 * (n - 1)/(n + 1)
      vb2 <- 24 * n * (n - 2) * (n - 3)/((n + 1)^2 * (n + 3) * (n + 5))
      m3  <- (6 * (n^2 - 5*n + 2)/((n + 7) * (n + 9))) * sqrt((6 * (n + 3) * (n + 5))/(n * (n - 2) * (n - 3)))
      a   <- 6 + (8/m3) * (2/m3 + sqrt(1 + 4/m3^2))
      xx  <- (b - eb2)/sqrt(vb2)
      cr <- sign((1 - 2/a)/(1 + xx * sqrt(2/(a - 4)))) * abs((1 - 2/a)/(1 + xx * sqrt(2/(a - 4))))^(1/3)
      z   <- (1 - 2/(9 * a) - cr ) / sqrt(2/(9 * a))
      return(z)
    }
    zval <- (anscombe.zvalue(x) - anscombe.zvalue(y))/sqrt(2 - 2*cor(x,y)^4)
    pval <- (1 - pnorm(abs(zval))) * 2
    return(list(z.value = zval, p.value = pval))
  }
  myanscombe.test <- function(x, alternative = c("two.sided", "less", "greater")){
    n   <- length(x)
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    s <- match.arg(alternative)
    alter <- switch(s, two.sided = 0, less = 1, greater = 2)
    b   <- n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    eb2 <- 3 * (n - 1)/(n + 1)
    vb2 <- 24 * n * (n - 2) * (n - 3)/((n + 1)^2 * (n + 3) * (n + 5))
    m3  <- (6 * (n^2 - 5*n + 2)/((n + 7) * (n + 9))) * sqrt((6 * (n + 3) * (n + 5))/(n * (n - 2) * (n - 3)))
    a   <- 6 + (8/m3) * (2/m3 + sqrt(1 + 4/m3^2))
    xx  <- (b - eb2)/sqrt(vb2)
    cr <- sign((1 - 2/a)/(1 + xx * sqrt(2/(a - 4)))) * abs((1 - 2/a)/(1 + xx * sqrt(2/(a - 4))))^(1/3)
    z   <- (1 - 2/(9 * a) - cr ) / sqrt(2/(9 * a))
    pval <- pnorm(z, lower.tail = FALSE)
    if (alter == 0) { pval <- 2 * pval; if (pval > 1) pval <- 2 - pval; alt <- "kurtosis is not equal to 3" }
    else if (alter == 1) { alt <- "kurtosis is greater than 3" }
    else { pval <- 1 - pval; alt <- "kurtosis is lower than 3" }
    RVAL <- list(statistic = c(kurt = b, z = z), p.value = pval,
                 alternative = alt, method = "Anscombe-Glynn kurtosis test",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }

  # --- Input checks ---
  if (is.null(pred)) stop("Tentative predictor is missing.")
  if (B < 0) stop("Number of resamples 'B' must be positive.")
  if (conf.level < 0 || conf.level > 1) stop("'conf.level' must be between 0 and 1")
  if (!boot.type %in% c("bca", "perc")) stop("Unknown argument in boot.type.")

  # --- Prepare outcome, predictor, and model matrix for covariates ---
  if (!inherits(formula, "formula")) {
    X <- if (is.matrix(formula$x)) formula$x else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(formula$y)) formula$y else model.response(model.frame(formula))
    delete.pred <- which(colnames(X) == pred)
    if (length(delete.pred) == 0) stop("Specified predictor not found in the target model.")
    x <- X[, delete.pred]
    X <- X[, -delete.pred]
    if (!is.matrix(X)) X <- as.matrix(X)
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    delete.pred <- which(colnames(X) == pred)
    if (length(delete.pred) == 0) stop("Specified predictor not found in the target model.")
    x <- X[, delete.pred]
    X <- X[, -delete.pred]
    if (!is.matrix(X)) X <- as.matrix(X)
  }

  ry <- as.vector(scale(lm.fit(X, y)$residuals))
  rx <- as.vector(scale(lm.fit(X, x)$residuals))
  resid_df <- data.frame(rx, ry)

  tar <- tryCatch(if (robust) mblm::mblm(ry ~ rx, data = resid_df, repeated = TRUE) else lm(ry ~ rx, data = resid_df), error = function(e) NULL)
  alt <- tryCatch(if (robust) mblm::mblm(rx ~ ry, data = resid_df, repeated = TRUE) else lm(rx ~ ry, data = resid_df), error = function(e) NULL)
  if (is.null(tar) || is.null(alt)) stop("Model fitting failed.")

  dat <- data.frame(alternative = as.vector(scale(resid(alt))),
                    target = as.vector(scale(resid(tar))),
                    pred.adj = rx,
                    out.adj = ry)

  if (isTRUE(prob.trans)) {
    trans <- prob.int(rx, ry)
    rx.trans <- trans$x
    ry.trans <- trans$y
    tar.trans <- tryCatch(if (robust) mblm::mblm(ry.trans ~ rx.trans) else lm(ry.trans ~ rx.trans), error = function(e) NULL)
    alt.trans <- tryCatch(if (robust) mblm::mblm(rx.trans ~ ry.trans) else lm(rx.trans ~ ry.trans), error = function(e) NULL)
    dat$alternative.trans <- if (!is.null(alt.trans)) as.vector(scale(resid(alt.trans))) else rep(NA, nrow(dat))
    dat$target.trans     <- if (!is.null(tar.trans)) as.vector(scale(resid(tar.trans))) else rep(NA, nrow(dat))
  }

  # --- Run separate normality tests ---
  if (isTRUE(prob.trans)) {
    agostino.out <- apply(dat[c("alternative.trans", "target.trans")], 2, moments::agostino.test)
    names(agostino.out) <- c("alternative", "target")
    agostino.out <- lapply(agostino.out, unclass)
    anscombe.out <- apply(dat[c("alternative.trans", "target.trans")], 2, myanscombe.test)
    names(anscombe.out) <- c("alternative", "target")
    anscombe.out <- lapply(anscombe.out, unclass)
  } else {
    agostino.out <- apply(dat[c("alternative", "target")], 2, moments::agostino.test)
    agostino.out <- lapply(agostino.out, unclass)
    anscombe.out <- apply(dat[c("alternative", "target")], 2, myanscombe.test)
    anscombe.out <- lapply(anscombe.out, unclass)
  }
  agostino.out$alternative[3:5] <- NULL
  agostino.out$target[3:5] <- NULL
  anscombe.out$alternative[3:5] <- NULL
  anscombe.out$target[3:5] <- NULL

  output <- list(agostino.out, anscombe.out)
  names(output) <- c("agostino", "anscombe")
  output$anscombe$target$statistic[1] <- output$anscombe$target$statistic[1] - 3
  output$anscombe$alternative$statistic[1] <- output$anscombe$alternative$statistic[1] - 3

  # --- Run asymptotic difference tests ---
  if (isTRUE(prob.trans)) {
    output <- c(output,
                list(skewdiff = unlist(skew.diff.test(dat$alternative.trans, dat$target.trans))),
                list(kurtdiff = unlist(kurt.diff.test(dat$alternative.trans, dat$target.trans)))
    )
    output$skewdiff <- c(moments::skewness(dat$alternative.trans)^2 - moments::skewness(dat$target.trans)^2, output$skewdiff)
    output$kurtdiff <- c((moments::kurtosis(dat$alternative.trans)-3)^2 - (moments::kurtosis(dat$target.trans)-3)^2, output$kurtdiff)
  } else {
    output <- c(output,
                list(skewdiff = unlist(skew.diff.test(dat$alternative, dat$target))),
                list(kurtdiff = unlist(kurt.diff.test(dat$alternative, dat$target)))
    )
    output$skewdiff <- c(moments::skewness(dat$alternative)^2 - moments::skewness(dat$target)^2, output$skewdiff)
    output$kurtdiff <- c((moments::kurtosis(dat$alternative)-3)^2 - (moments::kurtosis(dat$target)-3)^2, output$kurtdiff)
  }

  # --- Run bootstrap confidence intervals ---
  boot.warning <- FALSE
  if (B > 0) {
    suppressWarnings(boot.res <- boot::boot(dat, boot.diff, R = B, prob.trans = prob.trans, robust = robust))
    if (boot.type == "bca" && any(is.na(boot::empinf(boot.res)))) {
      boot.warning <- "Bootstrap failed: NA statistics. Try using fewer covariates or different data."
    } else {
      boot.out <- lapply(as.list(1:7), function(i, boot.res) boot::boot.ci(boot.res, conf=conf.level, type=boot.type, t0=boot.res$t0[i], t=boot.res$t[,i]), boot.res=boot.res)
      names(boot.out) <- c("skew.diff", "kurt.diff", "cor12.diff", "cor13.diff", "RHS3", "RCC", "RHS4")
      ci.skewdiff  <- unclass(boot.out$skew.diff)[[4]][4:5] ; names(ci.skewdiff) <- c("lower", "upper")
      ci.kurtdiff  <- unclass(boot.out$kurt.diff)[[4]][4:5] ; names(ci.kurtdiff) <- c("lower", "upper")
      ci.cor12diff <- unclass(boot.out$cor12.diff)[[4]][4:5] ; names(ci.cor12diff) <- c("lower", "upper")
      ci.cor13diff <- unclass(boot.out$cor13.diff)[[4]][4:5] ; names(ci.cor13diff) <- c("lower", "upper")
      ci.RHS3      <- unclass(boot.out$RHS3)[[4]][4:5] ; names(ci.RHS3) <- c("lower", "upper")
      ci.RCC       <- unclass(boot.out$RCC)[[4]][4:5] ; names(ci.RCC) <- c("lower", "upper")
      ci.RHS4      <- unclass(boot.out$RHS4)[[4]][4:5] ; names(ci.RHS4) <- c("lower", "upper")
      output$skewdiff <- c(output$skewdiff, ci.skewdiff)
      output$kurtdiff <- c(output$kurtdiff, ci.kurtdiff)
      output <- c(output,
                  list(cor12diff = c(boot.res$t0[3], ci.cor12diff)),
                  list(cor13diff = c(boot.res$t0[4], ci.cor13diff)),
                  list(RHS3 = c(boot.res$t0[5], ci.RHS3)),
                  list(RCC = c(boot.res$t0[6], ci.RCC)),
                  list(RHS4 = c(boot.res$t0[7], ci.RHS4)),
                  list(boot.args = c(boot.type, conf.level, B)),
                  list(boot.warning = FALSE)
      )
      if (sign(output$anscombe$alternative$statistic[1]) != sign(output$anscombe$target$statistic[1])) { output$boot.warning <- TRUE }
    }
  }

  response.name <- all.vars(formula(formula))[1]
  output <- c(output, list(var.names = c(response.name, pred), probtrans = prob.trans))
  call_info <- list(
    "function_call" = match.call(),
    "function_name" = "dda.indep",
    "all_args" = as.list(match.call())[-1],
    "formula" = formula,
    "data_name" = deparse(substitute(data)),
    "original_data" = if(missing(data) || is.null(data)) NULL else data
  )
  output <- c(output, list(call_info = call_info))
  class(output) <- "dda.resdist"
  return(output)
}

#' @name print.dda.resdist
#' @title Print Method for \code{dda.resdist} Objects
#'
#' @description \code{print} returns DDA test statistics associated with \code{dda.resdist} objects.
#'
#' @param x An object of class \code{dda.resdist} when using \code{print}.
#' @param ... Additional arguments to be passed to the method.
#'
#' @examples
#' print(result)
#'
#' @export
#' @rdname dda.resdist
#' @method print dda.resdist
print.dda.resdist <- function(x, ...){
  object <- x

  varnames <- object$var.names
  cat("\n")
  cat("DIRECTION DEPENDENCE ANALYSIS: Residual Distributions", "\n", "\n")
  cat("Skewness and kurtosis tests:", "\n")

  sigtests <- rbind( c(object[[1]]$target$statistic, object[[1]]$target$p.value, object[[1]]$alternative$statistic, object[[1]]$alternative$p.value),
                     c(object[[2]]$target$statistic, object[[2]]$target$p.value, object[[2]]$alternative$statistic, object[[2]]$alternative$p.value)
  )
  sigtests <- round(sigtests, 4)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c("target", "z-value", "Pr(>|z|)", "alternative", "z-value", "Pr(>|z|)")
  print.default(format( sigtests, digits = max(3L, getOption("digits") - 3L), scientific = NA, scipen = 999), print.gap = 2L, quote = FALSE)

  if(is.null(object$boot.args)){
    cat("\n")
    cat("Skewness and kurtosis difference tests:", "\n")

    citests <- rbind(object$skewdiff, object$kurtdiff)
    citests <- round(citests, 4)
    rownames(citests) <- c("Skewness", "Kurtosis")
    colnames(citests) <- c("diff", "z-value", "Pr(>|z|)")
    print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  if(!is.null(object$boot.args)){
    ci.level <- as.numeric(object$boot.args[2]) * 100
    cat("\n")

    if(object$boot.args[1] == "bca")  cat("Skewness and kurtosis difference tests and ", ci.level, "% ", "BCa bootstrap CIs:", "\n", "\n", sep = "")
    if(object$boot.args[1] == "perc") cat("Skewness and kurtosis difference tests and ", ci.level, "% ", "Percentile bootstrap CIs:", "\n", "\n", sep = "")

    citests <- rbind(object$skewdiff, object$kurtdiff)
    citests <- round(citests, 4)
    rownames(citests) <- c("Skewness", "Kurtosis")
    colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
    print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")

    if(object$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for joint higher moment differences:", "\n", sep = "")
    if(object$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for joint higher moment differences:", "\n", sep = "")

    jointtests <- rbind(object$cor12diff, object$RHS3, object$cor13diff, object$RHS4, object$RCC)
    jointtests <- round(jointtests, 4)
    rownames(jointtests) <- c("Co-Skewness", "Hyvarinen-Smith (Co-Skewness)", "Co-Kurtosis", "Hyvarinen-Smith (Co-Kurtosis)", "Chen-Chan (Co-Kurtosis)")
    colnames(jointtests) <- c("estimate", "lower", "upper")
    print.default(format( jointtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

    cat("\n")
    cat(paste("Number of resamples:", object$boot.args[3]))
    cat("\n")
  }

  cat("---")
  cat("\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], sep = " "))
  cat("\n")
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], sep = " "))
  cat("\n")
  if(isFALSE(object$probtrans)){
    cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1], sep = " "))
  } else {
    cat(paste("      Under prob.trans = TRUE, skewness and kurtosis differences < 0 and", "\n", "     co-skewness and co-kurtosis differences > 0 suggest", varnames[2], "->", varnames[1], sep = " "))
  }
  cat("\n")
  if(object$boot.warning) { cat("Warning: Excess-kurtosis values of residuals have unequal signs", "\n", "        Also compute Co-Kurtosis and Hyvarinen-Smith Co-Kurtosis for", varnames[1], "->", varnames[2], "\n") }
  cat("\n")
}

#' @name print.dda.resdist
#' @title Print Method for \code{dda.resdist} Objects
#'
#' @description \code{print} returns DDA test statistics associated with \code{dda.resdist} objects.
#'
#' @param x An object of class \code{dda.resdist} when using \code{print}.
#' @param ... Additional arguments to be passed to the method.
#'
#' @examples
#' print(result)
#'
#' @export
#' @rdname dda.resdist
#' @method print dda.resdist
print.dda.resdist <- function(x, ...){
  object <- x

  varnames <- object$var.names
  cat("\n")
  cat("DIRECTION DEPENDENCE ANALYSIS: Residual Distributions", "\n", "\n")
  cat("Skewness and kurtosis tests:", "\n")

  sigtests <- rbind( c(object[[1]]$target$statistic, object[[1]]$target$p.value, object[[1]]$alternative$statistic, object[[1]]$alternative$p.value),
                     c(object[[2]]$target$statistic, object[[2]]$target$p.value, object[[2]]$alternative$statistic, object[[2]]$alternative$p.value)
  )
  sigtests <- round(sigtests, 4)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c("target", "z-value", "Pr(>|z|)", "alternative", "z-value", "Pr(>|z|)")
  print.default(format( sigtests, digits = max(3L, getOption("digits") - 3L), scientific = NA, scipen = 999), print.gap = 2L, quote = FALSE)

  if(is.null(object$boot.args)){
    cat("\n")
    cat("Skewness and kurtosis difference tests:", "\n")

    citests <- rbind(object$skewdiff, object$kurtdiff)
    citests <- round(citests, 4)
    rownames(citests) <- c("Skewness", "Kurtosis")
    colnames(citests) <- c("diff", "z-value", "Pr(>|z|)")
    print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  if(!is.null(object$boot.args)){
    ci.level <- as.numeric(object$boot.args[2]) * 100
    cat("\n")

    # Print skewness and kurtosis results

    if(object$boot.args[1] == "bca")  cat("Skewness and kurtosis difference tests and ", ci.level, "% ", "BCa bootstrap CIs:", "\n", "\n", sep = "")
    if(object$boot.args[1] == "perc") cat("Skewness and kurtosis difference tests and ", ci.level, "% ", "Percentile bootstrap CIs:", "\n", "\n", sep = "")

    citests <- rbind(object$skewdiff, object$kurtdiff)
    citests <- round(citests, 4)
    rownames(citests) <- c("Skewness", "Kurtosis")
    colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
    print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")

    # Print co-skewness and co-kurtosis results

    if(object$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for joint higher moment differences:", "\n", sep = "")
    if(object$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for joint higher moment differences:", "\n", sep = "")

    jointtests <- rbind(object$cor12diff, object$RHS3, object$cor13diff, object$RHS4, object$RCC)
    jointtests <- round(jointtests, 4)
    rownames(jointtests) <- c("Co-Skewness", "Hyvarinen-Smith (Co-Skewness)", "Co-Kurtosis", "Hyvarinen-Smith (Co-Kurtosis)", "Chen-Chan (Co-Kurtosis)")
    colnames(jointtests) <- c("estimate", "lower", "upper")
    print.default(format( jointtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

    # hoctests <- rbind(object$cor12diff, object$cor13diff)
    # hoctests <- round(hoctests, 4)
    # rownames(hoctests) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]" )
    # colnames(hoctests) <- c("estimate", "lower", "upper")
    # print.default(format( hoctests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    # cat("\n")

    # # Print LR approximative measures

    # if(object$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for Likelihood Ratio approximations:", "\n", sep = "")
    # if(object$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for Likelihood Ratio approximations:", "\n", sep = "")

    # LRtests <- rbind(object$RHS3, object$RHS4, object$RCC )
    # LRtests <- round(LRtests, 4)
    # rownames(LRtests) <- c("Hyvarinen-Smith (co-skewness)", "Hyvarinen-Smith (co-kurtosis)", "Chen-Chan (co-kurtosis)")
    # colnames(LRtests) <- c("estimate", "lower", "upper")
    # print.default(format( LRtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

    cat("\n")
    cat(paste("Number of resamples:", object$boot.args[3]))
    cat("\n")
  }

  cat("---")
  cat("\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], sep = " "))
  cat("\n")
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], sep = " "))
  cat("\n")
  if(isFALSE(object$probtrans)){
    cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1], sep = " "))
  } else {
    cat(paste("      Under prob.trans = TRUE, skewness and kurtosis differences < 0 and", "\n", "     co-skewness and co-kurtosis differences > 0 suggest", varnames[2], "->", varnames[1], sep = " "))
  }
  cat("\n")
  if(object$boot.warning) { cat("Warning: Excess-kurtosis values of residuals have unequal signs", "\n", "        Also compute Co-Kurtosis and Hyvarinen-Smith Co-Kurtosis for", varnames[1], "->", varnames[2], "\n") }
	cat("\n")
  }


