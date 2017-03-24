#' predict.WRM
#'
#' @description
#' Model predictions for \code{WRM}
#'
#' @param object   A model object of class \code{WRM}
#' @param ... Other arguments passed to \code{predict}
#' @param newdata  A data frame containing variables used to make predictions.
#' @param sm       Logical. Should part of smooth components be included?
#' @param newcoord New coordinates corresponding to observations in \code{newdata}.
#'
#' @return   A vector of predictions based on a user-specified model
#'
#' @author Gudrun Carl,
#'         Sam Levin
#'
#' @examples
#' data(musdata)
#' coords<- musdata[,4:5]
#'
#' mwrm<-WRM(musculus ~ pollution + exposure, "poisson", musdata,
#'           coord=coords, level=0, plot=TRUE)

#' pred<-predict(mwrm,newdata=musdata)
#'
#' @export
#'


predict.WRM<-function(object,...,newdata,sm=FALSE,newcoord=NA){

  data<-newdata
  formula<-object$formula
  family<-object$family
  b<-object$b
  bsm<-object$bsm
  level<-object$level
  padzone<-object$padzone
  padform<-object$padform
  coord<-newcoord

  X<-model.matrix(formula,data)
  nvar<-dim(X)[2]
  n<-dim(data)[1]
  l<-dim(data)[2]

  if(sm) { # add part of smooth components
    if(is.na(sum(coord))) stop("coordinates are required")
    x<-coord[,1]
    y<-coord[,2]
    if(length(x)!=n) stop("error in dimension")
    logic1<-identical(as.numeric(x),round(x,0))
    logic2<-identical(as.numeric(y),round(y,0))
    if(!logic1 | !logic2) stop("coordinates not integer")

    n.level<-level
    length.s<-3*n.level+1
    s<-rep(1,length.s)
    s[length.s]<-0
    if(level==0) {s<-c(1,1,1,1) ; n.level<-1}

    pdim<- max(max(y)-min(y)+1,max(x)-min(x)+1)*padzone
    power<-0
    while(2^power<pdim) power<-power+1
    xmargin<-as.integer((2^power-(max(x)-min(x)))/2)-min(x)+1
    ymargin<-as.integer((2^power-(max(y)-min(y)))/2)-min(y)+1

    T<-array(NA,c(2^power,2^power,nvar))

    for(ii in 1:n){
      kx<-x[ii]+xmargin
      ky<-y[ii]+ymargin
      for (i3 in 1:nvar)
        T[ky,kx,i3]<-X[ii,i3]
    }  # ii loop

    P<-which(is.na(T), arr.ind = TRUE)
    if(padform==0){
      for (i3 in 1:nvar){
        i1<-P[which(P[,3]==i3),1]
        i2<-P[which(P[,3]==i3),2]
        for(i in 1:length(i1)) T[i1[i],i2[i],i3]<-0
      }
    }
    if(padform==1){
      for (i3 in 1:nvar){
        i1<-P[which(P[,3]==i3),1]
        i2<-P[which(P[,3]==i3),2]
        for(i in 1:length(i1)) T[i1[i],i2[i],i3]<-mean(T[,,i3], na.rm=TRUE)
      }
    }
    if(padform==2){
      for (i3 in 1:nvar) T[,,i3]<-padding(T[,,i3])
    }

    p<-2^power*2^power
    tt<-matrix(0,p,nvar)
    tt0<-matrix(0,p,nvar)
    for (i3 in 1:nvar){
      TT<-waveslim::mra.2d(T[,,i3],object$wavelet,n.level,method=object$wtrafo)
      TTS<-rep(0,length(TT[[1]]))
      TT0<-rep(0,length(TT[[1]]))
      for(is in 1:length(s)){
        if(s[is]==1) TTS <- TTS + TT[[is]]
        if(s[is]==0) TT0 <- TT[[is]]
      }
      tt[,i3]<-as.vector(TTS)
      if(level!=0) tt0[,i3]<-as.vector(TT0)
    }

    if(level!=0) lin<- tt%*%b  + tt0%*%bsm
    if(level==0) lin<- tt%*%b
    Fitted<-matrix(lin,2^power,2^power)
    fitted<-rep(0,n)
    for(i in 1:n) fitted[i]<-Fitted[y[i]+ymargin,x[i]+xmargin]
  } # add part of smooth components

  if(!sm) { # only part of detail components
    fitted<- X%*%b
  } # only part of detail components

  lin<-as.vector(fitted)
  if(family=="gaussian") pi<-lin
  if(family=="binomial") pi<-exp(lin)/(1+exp(lin))
  if(family=="poisson")  pi<-exp(lin)
  predict<-pi
}
