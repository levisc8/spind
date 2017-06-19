#'@title Adjusted actual values
#'
#'@description Adjusts actual presence/absence data based on the autocorrelation
#'in the predictions of a model. The function will optionally plot results of
#'model predictions, un-modified actual presence/absence, and adjusted  values.
#'
#'@param data a dataframe or matrix containing actual presence/absence (binary,
#'0 or 1) values in 1st column and predicted values (numeric between 0 and 1)
#'in 2nd column.
#'@param coord a matrix of two columns of the same length providing integer,
#'consecutively numbered coordinates for each occurence and prediction in
#'\code{data}.
#'@param plot.maps A logical indicating whether maps should be plotted.
#'Default is FALSE.
#' @param color.maps A logical value. If \code{TRUE}, produces colorful maps.
#' If \code{FALSE}, produces grayscale maps. Default is grayscale.
#'
#'@return A vector of adjusted actual values.
#'
#'@author Gudrun Carl
#'
#'@examples
#'data(hook)
#'data<- hook[,1:2]
#'coord<- hook[,3:4]
#'aa<-adjusted.actuals(data,coord,plot.maps=TRUE)
#'
#'@export

adjusted.actuals<-function(data, coord, plot.maps = FALSE, color.maps = FALSE){

  x <- coord[ ,1]
  y <- coord[ ,2]
  fb <- data[ ,1]
  fa <- data[ ,2]

  if(length(x) != length(fa)) stop("coordinates and data have different dimensions")
  logic1 <- identical(as.numeric(x), round(x, 0))
  logic2 <- identical(as.numeric(y), round(y, 0))
  if(!logic1 | !logic2) stop("coordinates not integer")

  fbs <- fb
  ac01a <- acfft(coord, fa, lim1 = 0, lim2 = 1, dmax = 1)
  ac01b <- acfft(coord, fb, lim1 = 0, lim2 = 1, dmax = 1)
  if(ac01a > 0.05 & ac01b > 0.05){
    ac01 <- ac01a - ac01b
    if(ac01 > 0.02){
      alpha <- ac01
      D <- as.matrix(dist(coord))
      R <- alpha^D
      spatial.W <- R^3
      ac01s <- acfft(coord, fbs, lim1 = 0, lim2 = 1, dmax = 1)
      while(ac01a > ac01s){
        fbs <- spatial.W %*% fbs
        ac01s <- acfft(coord, fbs, lim1 = 0, lim2 = 1, dmax = 1)
      }
      fbs <- fbs - min(fbs)
      fbs <- fbs / max(fbs)
    }
  }

  if (plot.maps){
    if(color.maps){
      colours <- list(colorRampPalette(RColorBrewer::brewer.pal(10, 'Spectral'))(50))
    } else {
      colours <- list(gray((0:50)/50))
    }
    a <- lattice::levelplot(fa ~ x + y,
                            col.regions = colours[[1]],
                            colorkey = FALSE,
                            scales = list(draw = FALSE),
                            xlab = "", ylab = "",
                            main = "predictions")

    b <- lattice::levelplot(fb ~ x + y,
                            col.regions = colours[[1]],
                            colorkey = FALSE,
                            scales = list(draw = FALSE),
                            xlab = "", ylab = "",
                            main = "actuals")

    c <- lattice::levelplot(fbs ~ x + y,
                            col.regions = colours[[1]],
                            colorkey = list(space = "bottom"),
                            scales = list(draw = FALSE),
                            xlab = "", ylab = "",
                            main = "adjusted actuals")

    tp <- lattice::trellis.par.get()
    lattice::trellis.par.set(list(axis.line = list(col = "transparent")))
    print(a, position = c(0.1, 0.09, 0.92, 0.98),
          split = c(1, 1, 2, 2), more = TRUE)
    print(b, position = c(0.12, 0.09, 0.94, 0.98),
          split = c(2, 1, 2, 2), more = TRUE)
    print(c, position = c(0.12, 0, 0.94, 1.02),
          split = c(2, 2, 2, 2), more = FALSE)
  } # plot

  fbs <- as.vector(fbs)
  return(fbs)
}
