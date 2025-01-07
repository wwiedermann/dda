#' @title Plot of \code{cdda.indep} class object
#' @description \code{plot} returns test statistic results from the \code{cdda.indep} class object.
#' @name plot.cddaindep
#'
#' @param obj           a \code{cdda.indep} class object
#' @param stat          a character indicating the stat of statistic to be plotted with options \code{c("hsic.diff", "dcor.diff", "mi.diff")}
#' @param ylim          a numeric vector of length 2 indicating the y-axis limits if \code{NULL}, the function will set the limits automatically
#' @returns             A plot of a test statistic result from the \code{cdda.indep} class object.
#'
#' @export
#' @rdname cdda.indep
#' @method plot cdda.indep
plot.cdda.indep <- function(obj = NULL, stat = NULL, ylim =  NULL){ #alpha = 0.05,

  if(is.null(stat)){
    stop("stat argument must be specified. as 'hsic.diff', 'dcor.diff', or 'mi.diff'")
  }
  if(stat != "hsic.diff" & stat != "dcor.diff" & stat != "mi.diff"){
    stop("stat must be one of 'hsic.diff', 'dcor.diff', or 'mi.diff'.")
  }

  ## Plot Label and Assignment ##
    mod.vals <- obj[[4]][["mod_data"]] #modvalues, raw data
    mod.levels <- x.axis.labels <- obj[[4]][["mod_levels"]] #"modval levels for pick-a-point, not raw data"
    y.title <- obj[[3]][[1]] #Test statistic CI header
    ci.level <- as.numeric(obj[[1]][[1]]$boot.args[2]) * 100

    ### Limits for X-axis ####

    if(is.numeric(mod.levels)){
      x.range <- c(min(mod.levels), max(mod.levels))
      plot.axis <- mod.levels
    }
    else {x.range <- c(1, length(mod.levels))
    plot.axis <- 1:length(mod.levels)
    }

    out <- matrix(NA, length(plot.axis), 6)

    if (stat == "hsic.diff"){
      if(is.null(obj[[1]][[1]]$out.diff[1, ])) {
        stop("HSIC Differences not found, set 'diff = TRUE'.")
      }
      y.title <- paste0("HSIC Differences (", ci.level, "% CI)")

      for (i in 1:length(plot.axis)) {
        out[i, ] <- c(obj[[1]][[i]]$out.diff[1, ],
                      obj[[2]][[i]]$out.diff[1, ])
      }
    }

    if (stat == "dcor.diff"){
      if(is.null(obj[[1]][[1]]$out.diff[2, ])){
        stop( "dCor Differences not found. Specify diff = TRUE." )
      }

      y.title <- paste0("dCor Differences (", ci.level, "% CI)") #y-title CI flexible
      for (i in 1:length(plot.axis)) {
        out[i, ] <- c(obj[[1]][[i]]$out.diff[2, ],
                      obj[[2]][[i]]$out.diff[2, ])
        }
    }

    if (stat == "mi.diff"){
      if(is.null(obj[[1]][[1]]$out.diff[3, ])){
        stop( "Mutual Information differences not found. Specify diff = TRUE." )
      }

      y.title <- paste0("Mutual Information Differences (", ci.level, "% CI)")

      for (i in 1:length(plot.axis)) {
        out[i, ] <- c(obj[[1]][[i]]$out.diff[3, ],
                      obj[[2]][[i]]$out.diff[3, ])
        }
    }

    ### Indep Plotting ################################################

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
         xlab = "Moderator Values", ylab = y.title, main = "Target Model")
    axis(1, at = 1:length(out.tar$condition), labels = x.axis.labels)
    polygon(x = c(out.tar$condition, rev(out.tar$condition)),
            y = c(out.tar$out.low, rev(out.tar$out.upp)),
            border = FALSE, col = "lightgrey")
    points(out.tar$condition, out.tar$out.mean, type = "l")
    abline(h = 0, lty = "dashed")

    plot(plot.axis, out.alt$out.mean, type = "n",
         ylim = y.range, xlim = x.range, xaxt = "n",
         xlab = "Moderator Values", ylab = y.title, main = "Alternative Model")
    axis(1, at = 1:length(out.tar$condition), labels = x.axis.labels)
    polygon(x = c(out.alt$condition, rev(out.alt$condition)),
            y = c(out.alt$out.low, rev(out.alt$out.upp)),
            border = FALSE, col = "lightgrey")
    points(out.alt$condition, out.alt$out.mean, type = "l")
    abline(h = 0, lty = "dashed")

    par(mfrow = c(1, 1))
}
