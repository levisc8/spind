#' Plot wavelet variance/covariance
#'
#' @description Plots the wavelet variance or covariance for the specified model.
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
#'
#' @details Each variable or pair of variables in \code{formula} is passed to \code{wavevar} or
#' \code{wavecovar} internally, and the result is plotted as a function of \code{level}.
#'
#' @return    A list containing a vector of results.
#'
#' @examples
#' data(carlinadata)
#' coords<- carlinadata[,4:5]
#'
#' covar.plot(carlina.horrida ~ aridity + land.use,
#' carlinadata,coord=coords,wavelet="d4",
#' wtrafo='dwt',plot='covar')
#'
#' covar.plot(carlina.horrida ~ aridity + land.use,
#'            carlinadata,coord=coords,wavelet="d4",
#'            wtrafo='dwt',plot='var')
#'
#' @seealso \code{\link{wavevar}}, \code{\link{wavecovar}}
#'
#' @export
#'


covar.plot<-function(formula,data,coord,wavelet="haar",wtrafo="dwt",
                     plot="covar"){

  x<-coord[,1]
  y<-coord[,2]
  X<-model.matrix(formula,data)
  namvar<-dimnames(X)[[2]]
  if(namvar[1]!="(Intercept)") nvar1<-1
  if(namvar[1]=="(Intercept)") nvar1<-2
  nvar2<-dim(X)[2]
  namresp<-as.character(formula[[2]])
  resp<-model.frame(formula,data)[[1]]
  wvar0<-wavevar(resp,x,y,wavelet=wavelet,wtrafo=wtrafo)
  nscale<-length(wvar0)

  if(plot=="var"){
    # Variance
    wvar<-matrix(NA,nvar2,nscale)
    for (kk in nvar1:nvar2){
      wvar[kk,]<-wavevar(X[,kk],x,y,wavelet=wavelet,wtrafo=wtrafo)
    }
    plot(wvar0[1:nscale],type="b",ylim=c(-.1,.9),pch=16,
         ylab="Variance", xlab="Level",
         main=paste(paste(wavelet,wtrafo)," - wavelet variance") )
    for(kk in nvar1:nvar2){
      points(wvar[kk,1:nscale],pch=kk,type="b")
    }
    leg<-c(namvar[nvar1:nvar2],namresp)
    v<-nvar1:nvar2
    v<-c(v,16)
    legend(2.5,0.8,leg,pch=v)
  }

  if(plot=="covar"){
    # Covariance
    wcvar<-matrix(NA,nvar2,nscale)
    for (kk in nvar1:nvar2){
      wcvar[kk,]<-wavecovar(resp,X[,kk],x,y,wavelet=wavelet,wtrafo=wtrafo)
    }
    plot(wcvar[nvar1,1:nscale],type="b",ylim=c(-.1,.6),pch=nvar1,
         ylab="Covariance", xlab="Level",
         main=paste(paste(wavelet,wtrafo)," - wavelet covariance") )
    for(kk in (nvar1+1):nvar2){
      points(wcvar[kk,1:nscale],pch=kk,type="b")
    }
    leg<-rep(NA,nvar2-nvar1+1)
    for (kk in nvar1:nvar2){
      leg[kk]<-paste(namresp,namvar[kk],sep=" --- ")
    }
    if(nvar1==2) leg<-leg[-1]
    v<-nvar1:nvar2
    legend(2.5,0.5,leg,pch=v)
  }

  if(plot=="var") {
    res<-rbind(wvar0,wvar)
    rownames(res)<-c(namresp,namvar)
  }
  if(plot=="covar") {
    res<-wcvar
    rownames(res)<-paste(namresp,namvar,sep="-")
  }

  fit<-list(result=res)
  fit
}
