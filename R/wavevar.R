#' Wavelet variance analysis
#'
#' @description Calculates the variance of parameters in a wavelet
#' multiresolution analysis. This is a 2D analysis taking the grid
#' structure of datasets into account and provides scale-specific
#' results for data sampled on a contiguous geographical area. The
#' dataset is assumed to be regular gridded and the grid cells are
#' assumed to be square.
#'
#' @param f    A vector
#' @param x    Corresponding x-coordinates which have to be integer
#' @param y    Corresponding y-coordinates which have to be integer
#' @param wavelet   Name of wavelet family. \code{haar}, \code{d4}, and \code{la8}.
#' are possible. \code{haar} is the default.
#' @param wtrafo    Type of wavelet transform. Either \code{dwt} or \code{modwt}.
#' \code{dwt} is the default.
#'
#' @return Wavelet variance for \code{f}.
#'
#' @seealso \pkg{waveslim}, \code{\link{WRM}}, \code{\link{covar.plot}}
#'
#' @export

wavevar<-function(f,x,y,wavelet="haar",wtrafo="dwt"){

  n<-length(f)
  pdim<- max(max(y)-min(y)+1,max(x)-min(x)+1)
  power<-0
  while(2^power<pdim) power<-power+1
  xmargin<-as.integer((2^power-(max(x)-min(x)))/2)-min(x)+1
  ymargin<-as.integer((2^power-(max(y)-min(y)))/2)-min(y)+1
  f<-scale(f) # for scaling and centering
  Fmat<-matrix(0,2^power,2^power)
  for(ii in 1:n){
    kx<-x[ii]+xmargin
    ky<-y[ii]+ymargin
    Fmat[kx,ky]<-f[ii]
  } # ii loop
  level<-power
  p<-2^power*2^power
  if(wtrafo=="dwt") F.dwt<-waveslim::dwt.2d(Fmat,wavelet,level)
  if(wtrafo=="modwt") F.dwt<-waveslim::modwt.2d(Fmat,wavelet,level)
  # wavelet variance (1/n) * sum[abs(f.dwt)^2 ]
  Var<-rep(NA,level)
  for(ik in 1:level){
    ii<-3*(ik-1)+1
    FS1 <- (1/n) * sum(abs(F.dwt[[ii]])^2)
    FS2 <- (1/n) * sum(abs(F.dwt[[ii+1]])^2)
    FS3 <- (1/n) * sum(abs(F.dwt[[ii+2]])^2)
    Var[ik]<-FS1+FS2+FS3 # all 3 components
  }
  # windows()
  # plot(Var)
  Var<-round(Var,4)
  Var
}
