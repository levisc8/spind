#' Wavelet covariance analysis
#'
#' @description The analysis is based a wavelet multiresolution analysis.
#' It is a 2D analysis taking the grid structure of datasets into
#' account, i.e. the analysis provides scale-specific
#' results for data sampled on a contiguous geographical area. The
#' dataset is assumed to be regular gridded and the grid cells are
#' assumed to be square.
#'
#' @param f1   a vectors
#' @param f2   a vectors of same length
#' @param x    corresponding x-coordinates which have to be integer
#' @param y    corresponding y-coordinates which have to be integer
#' @param wavelet  type of wavelet: "haar" or "d4" or "la8"
#' @param wtrafo   type of wavelet transform: "dwt" or "modwt"
#'
#' @return Wavelet covariance
#'
#' @seealso \pkg{waveslim}, \code{\link{WRM}}


wavecovar<-function(f1,f2,x,y,wavelet="haar",wtrafo="dwt"){

  n<-length(f1)
  pdim<- max(max(y)-min(y)+1,max(x)-min(x)+1)
  power<-0
  while(2^power<pdim) power<-power+1
  xmargin<-as.integer((2^power-(max(x)-min(x)))/2)-min(x)+1
  ymargin<-as.integer((2^power-(max(y)-min(y)))/2)-min(y)+1
  f1<-scale(f1) # for scaling and centering
  f2<-scale(f2) # for scaling and centering
  F1<-matrix(0,2^power,2^power)
  F2<-matrix(0,2^power,2^power)
  for(ii in 1:n){
    kx<-x[ii]+xmargin
    ky<-y[ii]+ymargin
    F1[kx,ky]<-f1[ii]
    F2[kx,ky]<-f2[ii]
  } # ii loop
  level<-power
  p<-2^power*2^power
  if(wtrafo=="dwt"){
    F1.dwt<-waveslim::dwt.2d(F1,wavelet,level)
    F2.dwt<-waveslim::dwt.2d(F2,wavelet,level)
  }
  if(wtrafo=="modwt"){
    F1.dwt<-waveslim::modwt.2d(F1,wavelet,level)
    F2.dwt<-waveslim::modwt.2d(F2,wavelet,level)
  }
  # wavelet covariance (1/n) * sum[abs(f1.dwt*f2.dwt) ]
  CVar<-rep(NA,level)
  for(ik in 1:level){
    ii<-3*(ik-1)+1
    FS1 <- (1/n) * sum(abs(F1.dwt[[ii]]*F2.dwt[[ii]]))
    FS2 <- (1/n) * sum(abs(F1.dwt[[ii+1]]*F2.dwt[[ii+1]]))
    FS3 <- (1/n) * sum(abs(F1.dwt[[ii+2]]*F2.dwt[[ii+2]]))
    CVar[ik]<-FS1+FS2+FS3 # all 3 components
  }
  CVar
  plot(CVar)
  iiende<-level*3+1
  CVarende <- (1/n) * sum(abs(F1.dwt[[iiende]]*F2.dwt[[iiende]]))
  CVartotal<-sum(CVar)+CVarende
  CVartotal
  # windows()
  # plot(CVar)
  CVar<-round(CVar,4)
  CVar
}