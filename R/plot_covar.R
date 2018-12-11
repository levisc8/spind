#' Plot wavelet variance/covariance
#'
#' @description Plots the wavelet variance or covariance for the specified formula.
#' The scale-dependent results are graphically displayed.
#'
#' @param formula   With specified notation according to names in data frame.
#' @param data      Data frame.
#' @param coord     A matrix of 2 columns with
#'  corresponding x,y-coordinates which have to be integer.
#' @param wavelet  Type of wavelet: \code{haar}, \code{d4}, or \code{la8}.
#' @param wtrafo   Type of wavelet transform: \code{dwt} or \code{modwt}.
#' @param plot      Either \code{var} for wavelet variance analysis
#'           or \code{covar} for wavelet covariance analysis.
#' @param customize_plot Additional plotting parameters passed to \code{ggplot}.
#' NOW DEPRECATED
#'
#' @details Each variable or pair of variables in \code{formula} is passed to \code{wavevar} or
#' \code{wavecovar} internally, and the result is plotted as a function of \code{level}.
#'
#' @return    A list containing
#'
#'  1. \code{result} =  a vector of results.
#'
#'  2. \code{plot} = a \code{ggplot} object
#'
#' @author Gudrun Carl
#'
#' @examples
#' data(carlinadata)
#' coords <- carlinadata[,4:5]
#'
#' covariance <- covar.plot(carlina.horrida ~ aridity + land.use - 1,
#'                          data = carlinadata,
#'                          coord = coords,
#'                          wavelet = "d4",
#'                          wtrafo = 'modwt',
#'                          plot = 'covar')
#'
#' covariance$plot
#' covariance$result
#'
#' variance <- covar.plot(carlina.horrida ~ aridity + land.use - 1,
#'                        data = carlinadata,
#'                        coord = coords,
#'                        wavelet = "d4",
#'                        wtrafo = 'modwt',
#'                        plot = 'var')
#'
#' variance$plot
#' variance$result
#'
#' @seealso \code{\link{wavevar}}, \code{\link{wavecovar}}
#'
#' @importFrom ggplot2 theme element_blank element_text ggplot aes
#' geom_point geom_line scale_x_continuous scale_y_continuous
#' @importFrom stats model.matrix model.frame
#' @importFrom rlang quo !!
#' @export
#'


covar.plot<-function(formula,data,coord,wavelet="haar",wtrafo="dwt",
                     plot="covar", customize_plot = NULL){


  if(!is.null(customize_plot)) {
    warning('"customize_plot" argument is now soft deprecated.\n',
            'Use object_name$plot to print the ggplot2 object and for \n',
            'subsequent modification.')
  }


  x <- coord[ ,1]
  y <- coord[ ,2]
  X <- stats::model.matrix(formula, data)
  namvar <- dimnames(X)[[2]]
  if(namvar[1] != "(Intercept)") nvar1 <- 1
  if(namvar[1] == "(Intercept)") nvar1 <- 2
  nvar2 <- dim(X)[2]
  namresp <- as.character(formula[[2]])
  resp <- stats::model.frame(formula, data)[[1]]
  wvar0 <- wavevar(resp, coord, wavelet = wavelet,
                   wtrafo = wtrafo)
  nscale <- length(wvar0)

  if(plot == "var"){
    # Variance
    wvar <- matrix(NA, nvar2, nscale)
    for(kk in nvar1:nvar2){
      wvar[kk, ] <- wavevar(X[ ,kk], coord,
                            wavelet = wavelet, wtrafo = wtrafo)
    }

    VarCol <- character()
    Level <- rep(1:nscale, length(namvar[nvar1:nvar2]) + 1)

    for(i in nvar1:nvar2){
      tempdata <- rep(namvar[i], nscale)
      VarCol <- c(VarCol, tempdata)
    }

    VarCol <- c(VarCol, rep(namresp, nscale))

    Var <- as.vector(t(wvar))
    Var <- c(Var[!is.na(Var)], wvar0[1:nscale])

    PltData <- data.frame(Variable = VarCol, Level = Level,
                          Variance = Var)

    VarType <- "Variance"

  }

  if(plot == "covar"){
    # Covariance
    wcvar <- matrix(NA, nvar2, nscale)
    VarCol <- character()
    for(kk in nvar1:nvar2){
      wcvar[kk, ] <- wavecovar(resp, X[ ,kk], coord,
                               wavelet = wavelet, wtrafo = wtrafo)
      tempVar <- rep(paste(namresp, colnames(X)[kk], sep = ' - '), nscale)
      VarCol <- c(VarCol, tempVar)
    }

    Level <- rep(1:nscale, length(namvar[nvar1:nvar2]))

    Covar <- as.vector(t(wcvar))
    Covar <- Covar[!is.na(Covar)]

    PltData <- data.frame(Variable = VarCol, Level = Level,
                          Variance = Covar)

    VarType <- "Covariance"
  }




  SeqEnd <- max(PltData$Variance, na.rm = TRUE) + .1
  if(min(PltData$Variance, na.rm = TRUE) < 0){
    SeqBegin <- min(PltData$Variance, na.rm = TRUE)
  } else {
    SeqBegin <- 0
  }

  if(max(PltData$Variance > 0.5)){
    y_scale <- seq(SeqBegin, 1, .2)
    ylim <- c(SeqBegin,1)
  } else {
    y_scale <- seq(SeqBegin, SeqEnd, .1)
    ylim <- c(SeqBegin, SeqEnd)
  }

  Level <- rlang::quo(Level)
  Variable <- rlang::quo(Variable)
  Variance <- rlang::quo(Variance)

  plt.blank <-  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.background = ggplot2::element_blank(),
                               axis.line = ggplot2::element_line(colour = "black"))

  Plt <- ggplot2::ggplot(PltData,
                         ggplot2::aes(x = !! Level,
                                       y = !! Variance)) +
    plt.blank +
    ggplot2::geom_point(ggplot2::aes(colour = !! Variable,
                                      shape = !! Variable),
                        size = 3) +
    ggplot2::geom_line(ggplot2::aes(colour = !! Variable),
                       linetype = 2,
                       size = 1) +
    ggplot2::scale_x_continuous("Level", breaks = 1:nscale) +
    ggplot2::scale_y_continuous(paste("Wavelet ", VarType),
                                breaks = y_scale,
                                limits = ylim) +
    customize_plot



  if(plot == "var"){
    res <- rbind(wvar0, wvar)
    rownames(res) <- c(namresp, namvar)
  }
  if(plot == "covar"){
    res <- wcvar
    rownames(res) <- paste(namresp, namvar, sep = "-")
  }

  fit <- list(result = res,
              plot = Plt)
  return(fit)
}
