#' Upscaling of smooth components
#'
#' @description The analysis is based a wavelet multiresolution analysis
#' using only smooth wavelet components.
#' It is a 2D analysis taking the grid structure of datasets into
#' account, i.e. the analysis provides scale-specific
#' results for data sampled on a contiguous geographical area. The
#' dataset is assumed to be regular gridded and the grid cells are
#' assumed to be square.
#' The scale-dependent results are graphically displayed.
#'
#'
#' @param f     a vector
#' @param x     corresponding x-coordinates which have to be integer
#' @param y     corresponding y-coordinates which have to be integer
#' @param wavelet  type of wavelet: "haar" or "d4" or "la8"
#' @param wtrafo   type of wavelet transform: "dwt" or "modwt"
#' @param pad       a numeric value for padding into a bigger square.
#'           Mostly the mean of f values is useful.
#'           Default is set to mean(f).

upscale<-function(f,x,y,wavelet="haar",wtrafo="dwt",pad=mean(f)){

  if(length(f)!=length(x) | length(f)!=length(y)) stop("error in dim")
  logic1<-identical(as.numeric(x),round(x,0))
  logic2<-identical(as.numeric(y),round(y,0))
  if(!logic1 | !logic2) stop("coordinates are not integer")
  n<-length(f)
  pdim<- max(max(y)-min(y),max(x)-min(x))*5/4
  power<-0
  while(2^power<pdim) power<-power+1
  xmargin<-as.integer((2^power-(max(x)-min(x)))/2)-min(x)+1
  ymargin<-as.integer((2^power-(max(y)-min(y)))/2)-min(y)+1
  Fmat<-matrix(pad,2^power,2^power)
  for(ii in 1:n){
    kx<-x[ii]+xmargin
    ky<-y[ii]+ymargin
    Fmat[kx,ky]<-f[ii]
  } # ii loop

  ## Plot
  par(mfrow=c(2,2),
      mai=c(0.1,0,0.4,0),
      omi=c(0,0,0,0),
      pty="s",cex.main=1)
  colors<-gray((0:50)/50)
  minvec<-rep(NA,4)
  maxvec<-rep(NA,4)
  for (i in 1:4){
    if(i==1){
      FTS<-Fmat
      minvec[i]<-min(FTS)
      maxvec[i]<-max(FTS)
    }
    if(i!=1){
      level<-i-1
      FT<-waveslim::mra.2d(Fmat,wavelet,level,method=wtrafo)
      FTS<-FT[[3*level+1]]
      minvec[i]<-min(FTS)
      maxvec[i]<-max(FTS)
    }
  }
  minFTS<-min(minvec,na.rm = TRUE)
  maxFTS<-max(maxvec,na.rm = TRUE)
  for (i in 1:4){
    if(i==1){
      FTS<-Fmat
    }
    if(i!=1){
      level<-i-1
      FT<-waveslim::mra.2d(Fmat,wavelet,level,method=wtrafo)
      FTS<-FT[[3*level+1]]
    }
    FTS<-FTS-minFTS
    FTS<-FTS/(maxFTS-minFTS)
    image(FTS,zlim=c(0,1),axes=FALSE,col=colors,
          main=paste("level = ", i-1))
  } # i-loop
}