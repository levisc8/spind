#' @title Relative Variable Importance
#'
#' @description
#' Creates model selection tables, calculates and plots relative
#' variable importance based on the scale level of a given model.
#'
#' @details Calculates the relative importance of each variable
#' using multi-model inference methods in a wavelet multi-resolution regression
#' framework implemented in \code{mmiWMRR}. The scale level dependent
#' results are then graphically displayed.
#'
#'
#' @param formula   A model formula
#' @param family \code{gaussian}, \code{binomial}, and \code{poisson}
#'  are supported.
#' @param data A data frame or set of vectors of equal length.
#' @param coord X,Y coordinates for each observation. Coordinates should be
#' consecutive integers.
#' @param maxlevel   An integer for maximum scale level
#' @param detail   Remove smooth wavelets? If \code{TRUE}, only detail components are analyzed.
#' If set to \code{FALSE}, smooth and detail components are analyzed. Default is \code{TRUE}.
#' @param wavelet  Type of wavelet: \code{haar}, \code{d4}, or \code{la8}
#' @param wtrafo   Type of wavelet transform: \code{dwt} or \code{modwt}
#' @param n.eff    A numeric value of effective sample size
#' @param trace Should R print progress updates to the console? Default is FALSE
#' @param customize_plot Additional plotting parameters passed to \code{ggplot}.
#' NOW DEPRECATED.
#'
#' @return A list containing
#'
#' 1. A matrix containing the relative importance of each variable
#' in the regression at each value of the scale level.
#'
#' 2. A \code{ggplot} object containing a plot of the relative
#' variable importance
#'
#' @examples
#' data(carlinadata)
#' coords<- carlinadata[,4:5]
#'
#'\dontrun{
#'
#' wrm <- WRM(carlina.horrida ~ aridity + land.use,
#'            family = "poisson",
#'            data = carlinadata,
#'            coord = coords,
#'            level = 1,
#'            wavelet = "d4")
#'
#' mmi <- mmiWMRR(wrm, data = carlinadata, scale = 3, detail = TRUE)
#'
#'
#' # Plot scale-dependent relative variable importance
#' rvi <- rvi.plot(carlina.horrida ~ aridity + land.use,
#'                 family = "poisson",
#'                 data = carlinadata,
#'                 coord = coords,
#'                 maxlevel = 4,
#'                 detail = TRUE,
#'                 wavelet = "d4")
#'
#' rvi$plot
#' rvi$rvi
#'
#'}
#' @importFrom ggplot2 theme element_blank element_line ggplot
#' aes geom_point geom_line scale_x_continuous scale_y_continuous
#' @importFrom rlang !! quo

#' @export


rvi.plot <- function(formula, family, data, coord, maxlevel, detail = TRUE,
                     wavelet = "haar", wtrafo = "dwt",
                     n.eff = NULL, trace = FALSE, customize_plot = NULL){

  if(!is.null(customize_plot)) {
    warning('"customize_plot" argument is now soft deprecated.\n',
            'Use object_name$plot to print the ggplot2 object and for \n',
            'subsequent modification.')
  }
  if(trace){
    cat("\n","Model selection tables:","\n","\n")
  }

  wrm <- WRM(formula, family, data, coord, level = 1,
            wavelet = wavelet, wtrafo = wtrafo)

  mmi <- mmiWMRR(wrm, data, scale = 1, detail = detail)

  nrowA <- dim(mmi$result)[1]
  ncolA <- dim(mmi$result)[2]

  nvar <- dim(mmi$result)[2] - 6
  leg <- dimnames(mmi$result)[[2]][2:(nvar + 1)]

  A <- array(NA, c(nrowA, ncolA, maxlevel))
  level <- rep(NA, maxlevel)
  A[ , ,1] <- mmi$result
  level[1] <- mmi$level

  if(maxlevel >= 2){
    for (i in 2:maxlevel) {
      mmi <- mmiWMRR(wrm, data, scale = i, detail = detail, trace = trace)
      A[ , ,i] <- mmi$result
      level[i] <- mmi$level
    }
  }

  if(trace){
    cat("\n","---","\n","Relative variable importance:","\n","\n")
  }

  klimitscale <- dim(A)[3]
  ip <- dim(A)[1]

  WeightSums <- matrix(NA, nvar, klimitscale)
  for(kscale in 1:klimitscale){
    for(kvar in 2:(nvar + 1)){
      for (i in 1:ip){
        if(!is.na(A[i, kvar, kscale])){
          A[i, kvar, kscale] <- A[i, (nvar + 6), kscale]
        }
      }
    }
    B <- A[1:ip, 2:(nvar + 1), kscale]
    WeightSums[ ,kscale] <- colSums(B, na.rm = TRUE)
  } # kscale

  vec <- 1:nvar

  VarCol <- character()
  Level <- rep(1:maxlevel, length(leg))

  for(i in seq_len(length(leg))){
    tempdata <- rep(leg[i], maxlevel)
    VarCol <- c(VarCol, tempdata)
  }


  PltData <- data.frame(Variable = VarCol, Level = Level,
                        Weight = as.vector(t(WeightSums)))

  plt.blank <-  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.background = ggplot2::element_blank(),
                               axis.line = element_line(colour = "black"))


  Level <- rlang::quo(Level)
  Variable <- rlang::quo(Variable)
  Weight <- rlang::quo(Weight)

  Plt <- ggplot2::ggplot(PltData,
                         ggplot2::aes(x = !! Level,
                                       y = !! Weight)) +
         plt.blank +
    ggplot2::geom_point(ggplot2::aes(colour = !! Variable,
                                      shape = !! Variable),
                        size = 3) +
    ggplot2::geom_line(ggplot2::aes(colour = !! Variable),
                       linetype = 2,
                       size = 1) +
    ggplot2::scale_x_continuous("Level", breaks = 1:maxlevel) +
    ggplot2::scale_y_continuous("Relative Variable Importance",
                                breaks = seq(0,
                                             max(WeightSums),
                                             length.out = 6)) +
         customize_plot

  #print(Plt)

  rownames(WeightSums) <- leg
  colnames(WeightSums) <- paste("level", c(1:klimitscale), sep = "=")
  # Plot: scale-dependent relative variable importance
  if(trace){
    print(WeightSums)
  }

  fit <- list(rvi = WeightSums,
              plot = Plt)

  return(fit)

}


