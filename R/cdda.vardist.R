#' @title Conditional Directional Dependence Analysis: Variable Distributions
#' @description \code{cdda.vardist} computes DDA test statistics for observed
#'              variable distributions of competing conditional models
#'              (\code{y ~ x * m} vs.\code{x ~ y * m} with \code{m}
#'              being a continuous or categorical moderator).
#' @name cdda.vardist
#'
#' @param formula     Symbolic formula of the model to be tested or a \code{lm} object
#' @param pred        A character indicating the variable name of the predictor which serves as the outcome in the alternative model.
#' @param mod         A character indicating the variable name of the moderator.
#' @param modval      Characters or a numeric sequence specifying the moderator
#'                    values used in post-hoc probing. Possible characters
#'                    include \code{c("mean", "median", "JN")}. \code{modval = "mean"}
#'                    tests the interaction effect at the moderator values
#'                    \code{M - 1SD}, \code{M}, and \code{M + 1SD};
#'                    \code{modval = "median"} uses \code{Q1}, \code{Md},
#'                    and \code{Q3}. The Johnson-Neyman approach is applied
#'                    when \code{modval = "JN"} with conditional effects being
#'                    evaluated at the boundary values of the significance
#'                    regions. When a numeric sequence is specified, the
#'                    pick-a-point approach is used for the selected numeric values.
#' @param data        A required data frame containing the variables in the model.
#' @param conf.level  Confidence level for bootstrap confidence intervals.
#' @param B           Number of bootstrap samples.
#' @param boot.type  A character indicating the type of bootstrap confidence intervals. Must be one of the two values \code{c("perc", "bca")}. \code{boot.type = "bca"} is the default.
#'
#' @returns A list of class \code{cdda.vardist} containing the results of
#'          CDDA tests to evaluate distributional properties of observed
#'          variables for pre-specified moderator values.
#'
#' @examples
#' set.seed(321)
#' n <- 700
#'
#' ## --- generate moderator
#'
#' z <- sort(rnorm(n))
#' z1 <- z[z <= 0]
#' z2 <- z[z > 0]
#'
#' ## --- x -> y when z <= 0
#'
#' x1 <- rchisq(length(z1), df = 4) - 4
#' e1 <- rchisq(length(z1), df = 3) - 3
#' y1 <- 0.5 * x1 + e1
#'
#' ## --- y -> x when m z > 0
#'
#' y2 <- rchisq(length(z2), df = 4) - 4
#' e2 <- rchisq(length(z2), df = 3) - 3
#' x2 <- 0.25 * y2 + e2
#'
#' y <- c(y1, y2); x <- c(x1, x2)
#'
#' d <- data.frame(x, y, z)
#'
#' m <- lm(y ~ x * z, data = d)
#'
#' result <- cdda.vardist(m, pred = "x", mod = "z", B = 50,
#'                       modval = c(-1, 1), data = d)
#'
#' @references Wiedermann, W., & von Eye, A. (2025). \emph{Direction Dependence Analysis: Foundations and Statistical Methods}. Cambridge, UK: Cambridge University Press.
#' @seealso \code{\link{dda.vardist}} for an unconditional version.
#'
#' @export
#' @rdname cdda.vardist
cdda.vardist <- function(formula,
                         pred = NULL,
                         mod = NULL,
                         data = list(),
                         modval = NULL,
                         B = 200,
                         boot.type = "perc",
                         conf.level = 0.95
                      ){

  if (length(data) == 0) stop("Please specify a data frame.")

  if (!inherits(formula, "formula") ) {
    X <- if (is.matrix(formula$x) ) formula$x
    else model.matrix(terms(formula), model.frame(formula) )
    y <- if (is.vector(formula$y) ) formula$y
    else model.response(model.frame(formula))

    delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
    if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

    if(is.factor(data[,mod])){

      delete.mod <- which(colnames(X) %in% paste(mod, levels(data[,mod]), sep ="")) # get position of moderator
    }
    if(is.numeric(data[,mod])){
      delete.mod <- which(colnames(X) == mod)
    }

    if ( is.null(delete.mod) ) stop( "Specified moderator not found in the target model." )

    delete.interaction <- grep(":", colnames(X))
    if ( length(delete.interaction) == 0 ) stop( "No interaction term specified in the target model." )

    x   <- X[, delete.pred] # tentative predictor
    moderator <- X[, delete.mod]  # moderator

    X <- X[, -c(delete.pred, delete.mod, delete.interaction)] # model matrix with only covariates
    if (!is.matrix(X)) X <- as.matrix(X)

  }

  else {
    mf <- model.frame(formula, data = data)
    y  <- model.response(mf)   # tentative outcome
    X  <- model.matrix(formula, data = data)

    delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
    if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

    if(is.factor(data[,mod])){
      delete.mod <- which(colnames(X) %in% paste(mod, levels(data[,mod]), sep ="")) # get position of moderator
    }
    if(is.numeric(data[,mod])){
      delete.mod <- which(colnames(X) == mod)
    }
    if (  is.null(delete.mod) == TRUE ) stop( "Specified moderator not found in the target model." )

    delete.interaction <- grep(":", colnames(X))
    if ( length(delete.interaction) == 0 ) stop( "No interaction term specified in the target model." )

    x   <- X[,  delete.pred] # tentative predictor
    moderator <- X[, delete.mod]  # moderator

    X   <- X[, -c(delete.pred, delete.mod, delete.interaction)] # model matrix with only covariates

    if (!is.matrix(X)) X <- as.matrix(X)
  }

  # compute competing models

  target_model <- round(summary(lm(y ~ x + moderator*x + X))$coefficients, 4)

  alternate_model <- round(summary(lm(x ~ y + moderator*y + X))$coefficients, 4)

  response.name <- all.vars(formula(formula))[1]  # get name of response variable

  nam_tar <- rownames(target_model)
  nam_tar <- gsub("moderator", mod, rownames(target_model))
  nam_tar <- gsub("x", pred, nam_tar)
  rownames(target_model) <- sub("^X", "", nam_tar)

  nam_alt <- gsub("moderator", mod, rownames(alternate_model))
  nam_alt <- gsub("y", response.name, nam_alt)
  rownames(alternate_model) <- sub("^X", "", nam_alt)

  if (is.factor(data[,mod])) {

    moderator_levels <- unique(data[,mod])

    indep.temp.yx <- list()
    indep.temp.xy <- list()

    vardist.temp.yx <- list()
    vardist.temp.xy <- list()

    for (i in (levels(data[,mod]))) {

      moderator_rlvl <- relevel(data[,mod], ref = i)

      #Target model: x -> y
      design.tar <- model.matrix(~ X + moderator_rlvl + moderator_rlvl:x + x)

      design.tar <- as.data.frame(design.tar[,-1])
      design.tar$x <- NULL

      aux.yx.tar <- lm(y ~ -1 + as.matrix(design.tar))
      aux.xy.tar <- lm(x ~ -1 + as.matrix(design.tar))

      ry.aux.tar <- resid(aux.yx.tar)
      rx.aux.tar <- resid(aux.xy.tar)

      #Alternative model: y -> x
      design.alt<- model.matrix(~ X + moderator_rlvl + moderator_rlvl:y + y)
      design.alt <- as.data.frame(design.alt[,-1])
      design.alt$y <- NULL #(paste(head(design.alt)))

      aux.yx.alt <- lm(y ~ -1 + as.matrix(design.alt))
      aux.xy.alt <- lm(x ~ -1 + as.matrix(design.alt))

      ry.aux.alt <- resid(aux.yx.alt)
      rx.aux.alt <- resid(aux.xy.alt)

      #DDA Implementation
      indep.temp.yx[[i]] <- unclass(dda.vardist(ry.aux.tar ~ rx.aux.tar, pred = "rx.aux.tar",
                                                B = B, boot.type = boot.type, conf.level = conf.level
      ))
      indep.temp.xy[[i]] <- unclass(dda.vardist(rx.aux.alt ~ ry.aux.alt, pred = "ry.aux.alt",
                                                 B = B, boot.type = boot.type, conf.level = conf.level
      ))
    }

    names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", unique(data[,mod]), sep = ' ')
  }

  else {

    if (is.numeric(modval)){

      #if( is.null(modval) ) stop("No moderation values identified for pick-a-point approach.")

      modmat <- data.frame( matrix(NA, length(moderator), length(modval)) )

      for(i in seq_along(modval)){ modmat[,i] <- moderator - modval[i] }  # compute transformed moderator for each value
      colnames(modmat) <- modval
    }

    else if (modval == "mean"){

      modmat <- data.frame( (moderator - (mean(moderator) - sd(moderator)) ),
                            (moderator - mean(moderator) ),
                            (moderator - (mean(moderator) + sd(moderator)) )
      )
      colnames(modmat) <- c("M-1SD", "M", "M+1SD")
    }

    else if (modval == "median"){
      modmat <- data.frame( "Q1" = (moderator - quantile(moderator, probs = c(0.25)) ),
                            "Md" = (moderator - quantile(moderator, probs = c(0.50)) ),
                            "Q3" = (moderator - quantile(moderator, probs = c(0.75)) )
      )
    }


    else if (modval == "JN"){

      mdat <- data.frame(x, y, X, moderator)

      unwrpjn <- unclass(interactions::johnson_neyman( lm(y ~ x + moderator + moderator:x + X, data = mdat), pred = x,
                                         modx = moderator, control.fdr = FALSE ) )
      twobounds <- unwrpjn[["bounds"]]

      if(twobounds[1] < min(moderator, na.rm = TRUE) &
         twobounds[2] > max(moderator, na.rm = TRUE)) values <- NA

      if(twobounds[1] > min(moderator, na.rm = TRUE) &
         twobounds[2] > max(moderator, na.rm = TRUE)) values <- c(min(moderator), twobounds[1])

      if(twobounds[1] < min(moderator, na.rm = TRUE) &
         twobounds[2] < max(moderator, na.rm = TRUE)) values <- c(twobounds[2], max(moderator))

      if(twobounds[1] > min(moderator, na.rm = TRUE) &
         twobounds[2] < max(moderator, na.rm = TRUE)) values <- c(min(moderator), twobounds[1], twobounds[2], max(moderator))

      modmat <- data.frame( matrix(NA, length(moderator), length(values)) )
      for( i in seq_along(values) ){ modmat[,i] <- moderator - values[i] }  # compute transformed moderator for each value
      colnames(modmat) <- paste0(round(values, 2))
    }

    else stop("Invalid or no moderator value specified.")

    moderator_levels <- colnames(modmat)

    indep.temp.yx <- list()
    indep.temp.xy <- list()

    vardist.temp.yx <- list()
    vardist.temp.xy <- list()

    for ( i in 1:(ncol(modmat)) ) {

      #Target model: x -> y
      aux.yx.tar <- lm(y ~ X + modmat[,i] + modmat[,i]:x)
      aux.xy.tar <- lm(x ~ X + modmat[,i] + modmat[,i]:x)

      ry.aux.tar <- resid(aux.yx.tar)
      rx.aux.tar <- resid(aux.xy.tar)

      #Alternative model: y -> x
      aux.yx.alt <- lm(y ~ X +  modmat[,i] +  modmat[,i]:y)
      aux.xy.alt <- lm(x ~ X +  modmat[,i] +  modmat[,i]:y)

      ry.aux.alt <- resid(aux.yx.alt)
      rx.aux.alt <- resid(aux.xy.alt)

      #DDA Implementation
      indep.temp.yx[[i]] <- unclass(dda.vardist(ry.aux.tar ~ rx.aux.tar, B = B,
                                                pred = "rx.aux.tar", boot.type = boot.type))
      indep.temp.xy[[i]] <- unclass(dda.vardist(rx.aux.alt ~ ry.aux.alt, B = B,
                                                pred = "ry.aux.alt", boot.type = boot.type))

    }

    names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", colnames(modmat), sep = ' ')
  }

  cdda.var.output <- list()
  cdda.var.output[[1]] <- indep.temp.yx
  cdda.var.output[[2]] <- indep.temp.xy
  cdda.var.output[[3]] <- list("target_model" = target_model, "alternate_model" = alternate_model)
  cdda.var.output[[4]] <- list("response_name" = response.name,
                           "pred_name"= pred,
                           "mod_name" = mod,
                           "mod_levels" = moderator_levels,
                           "mod_data" = data[,mod]
  )
  names(cdda.var.output) <- c("cdda_target", "cdda_alternative", "models", "df_original")
  class(cdda.var.output) <- "cdda.vardist"
  return(cdda.var.output)

}

#' @name print.cdda.vardist
#' @title       Print Method for \code{cdda.vardist} Objects
#' @description \code{print} returns the output of standard linear model coefficients
#'              for competing target and alternative models.
#' @param x     An object of class \code{cdda.vardist} when using \code{print} or \code{plot}.
#' @param ...   Additional arguments to be passed to the function.
#'
#' @examples print(result)
#'
#' @export
#' @rdname cdda.vardist
#' @method print cdda.vardist
print.cdda.vardist <- function(x, ...){

  cdda.var.output <- x

  cat(paste("-----------------------------------------------------"))
  cat("\n")

  cat("OLS Summary: Target Model", "\n")
  cat("\n")
  print(cdda.var.output[[3]]$target_model)

  cat(paste("-----------------------------------------------------"))
  cat("\n", "\n")

  cat("OLS Summary: Alternative Model", "\n")
  cat("\n")
  print(cdda.var.output[[3]]$alternate_model)
}
