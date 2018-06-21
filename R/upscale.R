#' Upscaling of smooth components
#'
#' @description The analysis is based a wavelet multiresolution analysis
#' using only smooth wavelet components.
#' It is a 2D analysis taking the grid structure and provides scale-specific
#' results for data sampled on a contiguous geographical area. The
#' dataset is assumed to be regular gridded and the grid cells are
#' assumed to be square.
#' The scale-dependent results are graphically displayed.
#'
#'
#' @param f     A vector.
#' @param coord A matrix of two columns with corresponding cartesian
#' coordinates. Currently only supports integer coordinates.
#' @param wavelet   Name of wavelet family. \code{haar}, \code{d4}, and \code{la8}.
#' are possible. \code{haar} is the default.
#' @param wtrafo    Type of wavelet transform. Either \code{dwt} or \code{modwt}.
#' \code{dwt} is the default.
#' @param pad       A numeric value for padding the matrix
#' into a bigger square. Default is set to mean(f).
#' @param color.maps A logical value. If \code{TRUE}, produces colorful maps.
#' If \code{FALSE}, produces grayscale maps. Default is grayscale. NOW DEPRECATED,
#' color maps will not be produced in future versions.
#'
#' @return A set of plots showing the matrix image at each value for
#' \code{level}.
#'
#' @author Gudrun Carl
#'
#' @examples
#' data(carlinadata)
#' coords <- carlinadata[ ,4:5]
#'
#' # Upscaling of smooth components
#' upscale(carlinadata$land.use, coord = coords)
#'
#' @importFrom grDevices colorRampPalette gray
#' @importFrom graphics par
#' @importFrom waveslim mra.2d
#' @importFrom graphics image
#' @importFrom RColorBrewer brewer.pal
#' @export

upscale<-function(f, coord, wavelet = "haar", wtrafo = "dwt", pad = mean(f),
                  color.maps = FALSE){

  x <- coord[ ,1]
  y <- coord[ ,2]
  if(length(f) != length(x) | length(f) != length(y)) stop("error in dim")
  logic1 <- identical(as.numeric(x), round(x, 0))
  logic2 <- identical(as.numeric(y), round(y, 0))
  if(!logic1 | !logic2) stop("coordinates are not integer")
  n <- length(f)
  pdim <- max(max(y) - min(y), max(x) - min(x)) * 5 / 4
  power <- 0
  while(2^power < pdim) power <- power + 1
  xmargin <- as.integer((2^power - (max(x) - min(x))) / 2) - min(x) + 1
  ymargin <- as.integer((2^power - (max(y) - min(y))) / 2) - min(y) + 1
  Fmat <- matrix(pad, 2^power, 2^power)
  for(ii in 1:n){
    kx <- x[ii] + xmargin
    ky <- y[ii] + ymargin
    Fmat[kx,ky] <- f[ii]
  } # ii loop

  ## Plot
  graphics::par(mfrow = c(2, 2),
                mai = c(0.1, 0, 0.4, 0),
                omi = c(0, 0, 0, 0),
                pty = "s",cex.main=1)

  if(color.maps){
    colors <- list(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10,
                                                             'Spectral'))(100))
    warning('"color.maps" is now soft deprecated and will be removed in future versions')

  } else {
    colors <- list(grDevices::gray((0:50)/50))
  }
  minvec <- rep(NA, 4)
  maxvec <- rep(NA, 4)


  for (i in 1:4){
    if(i == 1){
      FTS <- Fmat
      minvec[i] <- min(FTS)
      maxvec[i] <- max(FTS)
    }
    if(i != 1){
      level <- i - 1
      FT <- waveslim::mra.2d(Fmat, wavelet, level, method = wtrafo)
      FTS <- FT[[3 * level + 1]]
      minvec[i] <- min(FTS)
      maxvec[i] <- max(FTS)
    }
  }
  minFTS <- min(minvec, na.rm = TRUE)
  maxFTS <- max(maxvec, na.rm = TRUE)
  for (i in 1:4){
    if(i == 1){
      FTS <- Fmat
    }
    if(i != 1){
      level <- i - 1
      FT <- waveslim::mra.2d(Fmat, wavelet, level, method = wtrafo)
      FTS <- FT[[3 * level + 1]]
    }
    FTS <- FTS - minFTS
    FTS <- FTS / (maxFTS - minFTS)
    graphics::image(FTS, zlim=c(0,1), axes = FALSE, col = colors[[1]],
                    main = paste("level = ", i - 1))

  } # i-loop


}
