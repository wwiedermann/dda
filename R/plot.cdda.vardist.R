#' @title Plot of \code{cdda.vardist} class object
#' @description \code{plot} returns test statistic results from the \code{cdda.vardist} class object.
#' @name plot.cddavardist
#'
#' @param obj    a \code{cddavardist} class object
#' @param stat   a character indicating the statistic to be plotted, default is \code{"rhs"}, with options \code{c("rhs", "cokurt", "rcc", "rtanh")}
#' @param ylim   a numeric vector of length 2 indicating the y-axis limits if \code{NULL}, the function will set the limits automatically
#' @returns      A plot of a test statistic result from the \code{cdda.vardist} class object.
#' @export
plot.cddavardist <- function(obj = NULL, stat = NULL,
                             ylim =  NULL, alpha = 0.05,
                             xlab = NULL, ylab = NULL, ...){

  if(class(obj) != "cddavardist"){
    stop("Object must be of class 'cddaindep' or 'cddavardist'.")
  }
  if(is.null(stat)){
    stop("stat argument must be specified. as 'rhs', 'cokurt', 'rcc', or 'rtanh'")
  }
  else if(stat != "rhs" & stat != "cokurt" & stat != "rcc" & stat != "rtanh"){
    stop("stat must be one of 'rhs', 'cokurt', 'rcc', or 'rtanh'.")
  }

  ## Plot Label and Assignment ##
  mod.vals <- obj[[4]][["mod_data"]] #modvalues, raw data
  mod.levels <- x.axis.labels <- obj[[4]][["mod_levels"]] #"modval levels for pick-a-point, not raw data"
  y.title <- obj[[3]][[1]] #Test statistic CI header

  tar.model.label <- paste0(obj[[4]]["response_name"], "|", obj[[4]]["mod_name"],
                            "\u2192", obj[[4]]["pred_name"]) # y ~ x | m
  alt.model.label <- paste0(obj[[4]]["pred_name"], "|", obj[[4]]["mod_name"],
                            "\u2192", obj[[4]]["response_name"]) # x ~ y | m

  ### Limits for X-axis ####

  if(is.numeric(mod.levels)){
    x.range <- c(min(mod.levels), max(mod.levels))
    plot.axis <- mod.levels
  }
  else {x.range <- c(1, length(mod.levels))
  plot.axis <- 1:length(mod.levels)
  }

  out <- matrix(NA, length(plot.axis), 6)

  if (stat == "rhs"){
    y.title <- "Hyvarinen-Smith Co-Skewness Differences (95% CI)" #make ci flexible
    if(is.null(obj[[1]][[1]]$RHS)){
      stop( "Hyvarinen-Smith Co-Skewness not found. Specify coskew = TRUE." )
      y.title <- "RHS"
    }

    for (i in 1:length(plot.axis)) {
      out[i, ] <- c(obj[[1]][[i]]$RHS, obj[[2]][[i]]$RHS)
    }
  }

  if (stat == "cokurt"){

    y.title <- "Co-Kurtosis Differences (95% CI)"
    if(is.null(obj[[1]][[1]]$cor13diff)){
      stop( "mitests not found. Specify cokurt = TRUE." )
      y.title <- "Co-Kurtosis"
    }

    for(i in 1:length(plot.axis)){
      out[i, ] <- c(obj[[1]][[i]]$cor13diff, obj[[2]][[i]]$cor13diff)
    }
  }

  if (stat == "rcc"){
    y.title <- "Chen-Chan Co-Kurtosis Differences (95% CI)"
    if(is.null(obj[[1]][[1]]$RCC)){
      stop( "Chen-Chan Co-Kurtosis Differences not found. Specify cokurt = TRUE." )
    }

    for (i in 1:length(plot.axis)) {
      out[i, ] <- c(obj[[1]][[i]]$RCC, obj[[2]][[i]]$RCC)
    }
  }

  if (stat == "rtanh"){ #check if the wrong stats are grabbed
    y.title <- "Hyperbolic Tangent Differences (95% CI)"
    if(is.null(obj[[1]][[1]]$Rtanh)){
      stop( "Hyperbolic Tangent Differences not found. Specify cokurt = TRUE." )
    }

    for (i in 1:length(plot.axis)) {
      out[i, ] <- c(obj[[1]][[i]]$Rtanh, obj[[2]][[i]]$Rtanh)
    }
  }

  ### Vardist Plotting ################################################

  ### Matrix Creation ###

  out.tar <- data.frame(condition = plot.axis, out.mean = out[,1],
                        out.low = out[,2], out.upp = out[,3])
  out.alt <- data.frame(condition = plot.axis, out.mean = out[,4],
                        out.low = out[,5], out.upp = out[,6])

  ### Set the y-axis ###

  if(is.null(ylim)) { y.range <- c(min(c(out.tar$out.low, out.alt$out.low)) - 0.1,
                                   max(c(out.tar$out.upp, out.alt$out.upp) + 0.1) )
  } else{y.range <- ylim}

  ### Start two plots ###

  par(mfrow = c(1, 2))

  plot(plot.axis, out.tar$out.mean, type = "n",
       ylim = y.range, xlim = x.range, xaxt = "n",
       xlab = tar.model.label, ylab = y.title, main = "Target Model")
  axis(1, at = 1:length(out.tar$condition), labels = x.axis.labels)
  polygon(x = c(out.tar$condition, rev(out.tar$condition)),
          y = c(out.tar$out.low, rev(out.tar$out.upp)),
          border = FALSE, col = "lightgrey")
  points(out.tar$condition, out.tar$out.mean, type = "l")
  abline(h = 0, lty = "dashed")


  plot(plot.axis, out.alt$out.mean, type = "n",
       ylim = y.range, xlim = x.range, xaxt = "n",
       xlab = alt.model.label, ylab = y.title, main = "Alternative Model")
  axis(1, at = 1:length(out.tar$condition), labels = x.axis.labels)
  polygon(x = c(out.alt$condition, rev(out.alt$condition)),
          y = c(out.alt$out.low, rev(out.alt$out.upp)),
          border = FALSE, col = "lightgrey")
  points(out.alt$condition, out.alt$out.mean, type = "l")
  abline(h = 0, lty = "dashed")
}

par(mfrow = c(1, 1))


