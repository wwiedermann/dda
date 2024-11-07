#' @title Conditional Direction Dependence Analysis: Independence Assumptions
#' @description \code{cdda.indep} tests the skewness and kurtosis of the variables of two competing models. It also tests the difference in skewness and kurtosis to assess independence properties of the two models. The function also provides bootstrap confidence intervals for the difference in skewness and kurtosis.
#' @name cdda.indep
#'
#' @param formula symbolic formula of the model to be tested or a \code{lm} object.
#' @param pred a character indicating the variable name of the predictor which serves as the outcome in the alternative model.
#' @param mod a character indicating the variable name of the moderator which serves as an interaction term in models.
#' @param modval For continuous moderators, character or numeric sequence specifying how the moderator is treated in the model. Possible characters include c("mean", "median", "JN"), and if a numeric sequence is specified, the 'pick a point' approach will be utilized.
#' @param data a required data frame containing the variables in the model (by default variables are taken from the environment which dda.vardist is called from)
#' @param hetero a logical value indicating whether separate homoscedasticity tests (i.e., standard and robust Breusch-Pagan tests) should be computed.
#' @param nlfun Either a numeric value or a function of .Primitive type used for non-linear correlation tests. When nlfun is numeric the value is used in a power tranformation.
#' @param hsic.method a character indicating the inference method for Hilbert-Schmidt Independence Criterion. Must be one of the two values c("gamma", "boot").
#' @param diff A logical value indicating whether differences in HSIC, dCor, and MI values should be computed. Bootstrap confidence intervals are computed using B bootstrap samples.
#' @param B number of bootstrap samples. Only used if hsic.method = "boot"
#' @param parallelize A logical value indicating whether bootstrapping is performed on multiple cores. Only used if hsic.method = "boot".
#' @param cores a numeric value indicating the number of cores. Only used if parallelize = TRUE
#'
#' @examples cdda.car.indep <- cdda.indep(mpg ~ wt * hp + qsec, pred = "wt",
#'                                    mod = "hp", diff = TRUE, data = mtcars)
#'           cdda.car.indep
#'           #OR
#'           car.test.mod <- lm(mpg ~ wt * hp + qsec, data = mtcars)
#'           cdda.indep(car.test.mod, pred = "wt", mod = "hp", diff = TRUE, data = mtcars)
#'
#' @returns A list of class \code{cddaindep} containing with information on competing models under moderation
#' @export
cdda.indep <- function(formula = NULL, pred = NULL, mod = NULL, modval = "mean",
                       data = list(), hetero = TRUE, diff = FALSE,
                       nlfun = NULL, hsic.method = "gamma",
                       B = 200, boot.type = "perc", conf.level = 0.95,
                       parallelize = FALSE, cores = 1,
                       values = NULL, JN.length = 3, ...) {
  library(boot)
  library(dHSIC)
  library(energy)
  library(lmtest)
  library(interactions)

  if (length(data) == 0) stop( "Please specify a data frame." )

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

    delete.interaction <- grep(":", colnames(X))  # which column refers to interaction term
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

    if ( is.null(delete.mod) == TRUE ) stop( "Specified moderator not found in the target model." )

    delete.interaction <- grep(":", colnames(X))  # which column refers to interaction term
    if ( length(delete.interaction) == 0 ) stop( "No interaction term specified in the target model." )

    x   <- X[,  delete.pred] # tentative predictor
    moderator <- X[, delete.mod]  # moderator

    X   <- X[, -c(delete.pred, delete.mod, delete.interaction)] # model matrix with only covariates

    if (!is.matrix(X)) X <- as.matrix(X)
  }

  #compute competing models
  target_model <- round(summary(lm(y ~ x + moderator * x + X))$coefficients, 4)
  alternate_model <- round(summary(lm(x ~ y + moderator * y + X))$coefficients, 4)

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
      design.alt$y <- NULL

      aux.yx.alt <- lm(y ~ -1 + as.matrix(design.alt))
      aux.xy.alt <- lm(x ~ -1 + as.matrix(design.alt))

      ry.aux.alt <- resid(aux.yx.alt)
      rx.aux.alt <- resid(aux.xy.alt)


      #DDA Implementation
      indep.temp.yx[[i]] <- unclass(dda.indep(ry.aux.tar ~ rx.aux.tar, pred = "rx.aux.tar",
                                              hsic.method = hsic.method, nlfun = nlfun,
                                              hsic.diff = hsic.diff, B = B, boot.type = boot.type,
                                              conf.level = conf.level, diff = diff, hetero = hetero
      ))
      indep.temp.xy[[i]] <- unclass(dda.indep(rx.aux.alt ~ ry.aux.alt, pred = "ry.aux.alt",
                                              hsic.method = hsic.method, nlfun = nlfun,
                                              hsic.diff = hsic.diff, B = B, boot.type = boot.type,
                                              conf.level = conf.level, diff = diff, hetero = hetero
      ))
    }

    names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", unique(data[,mod]), sep = ' ')
  }

  else {

    if (is.numeric(modval)){

      if( is.null(modval) ) stop("No moderation values identified for pick-a-point approach.")

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
      require(interactions)

      mdat <- data.frame(x, y, X, moderator)

      jnoutput <- unclass(johnson_neyman( lm(y ~ x + moderator + moderator:x + X, data = mdat), pred = x,
                                          modx = moderator, control.fdr = FALSE ) )
      twobounds <- jnoutput[["bounds"]]

      # Needs labeling, use values[i]  & round to 3 digits
      lwr <- max(twobounds[1], min(moderator))
      upr <- min(twobounds[2], max(moderator))
      values <- seq(from = lwr, to = upr, length.out = JN.length)

      modmat <- data.frame( matrix(NA, length(moderator), length(values)) )
      for( i in seq_along(values) ){ modmat[,i] <- moderator - values[i] }  # compute transformed moderator for each value
      colnames(modmat) <- paste0(round(values, 2))
    }

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
      indep.temp.yx[[i]] <- unclass(dda.indep(ry.aux.tar ~ rx.aux.tar, pred = "rx.aux.tar",
                                              hsic.method = hsic.method, nlfun = nlfun,
                                              hsic.diff = hsic.diff, B = B, boot.type = boot.type,
                                              conf.level = conf.level, diff = diff, hetero = hetero
      ))
      indep.temp.xy[[i]] <- unclass(dda.indep(rx.aux.alt ~ ry.aux.alt, pred = "ry.aux.alt",
                                              hsic.method = hsic.method, nlfun = nlfun,
                                              hsic.diff = hsic.diff, B = B, boot.type = boot.type,
                                              conf.level = conf.level, diff = diff, hetero = hetero
      ))

    }


    names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", colnames(modmat), sep = ' ')

  }

  response.name <- all.vars(formula(formula))[1]  # get name of response variable
  #output <- c(output, list(var.names = c(response.name, pred)))

  cdda.output <- list()
  cdda.output[[1]] <- indep.temp.yx
  cdda.output[[2]] <- indep.temp.xy
  cdda.output[[3]] <- list("target_model" = target_model, "alternate_model" = alternate_model)
  cdda.output[[4]] <- list("response_name" = response.name,
                           "pred_name"= pred,
                           "mod_name" = mod,
                           "mod_levels" = moderator_levels, # previously: modval,
                           "mod_data" = data[,mod] #raw for categ, not for contin?
                           )

  names(cdda.output) <- c("cdda_target", "cdda_alternative", "models", "df_original")
  class(cdda.output) <- "cddaindep"
  return(cdda.output)
}

#' @title Print method for "cddaindep" class
#' @description Calling `print` on a `cddaindep` object will display the output of the standard linear model coefficients for competing models.
#' @param x An object of class `cddaindep`
#'
#' @export
#'
print.cddaindep <- function(x, ...){
  cdda.output <- x
  cat("\n")
  cat(paste0("-----------------------------------------------------"))
  cat("\n")

  cat("OLS Summary: Target Model", "\n")
  cat("\n")
  print(round(cdda.output[[3]]$target_model, 4))

  cat(paste0("-----------------------------------------------------"))
  cat("\n")

  cat("OLS Summary: Alternative Model", "\n")
  cat("\n")
  print(round(cdda.output[[3]]$alternate_model, 4))
  cat("\n")
}

