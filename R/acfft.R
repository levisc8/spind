#' @title Spatial autocorrelation diagnostics
#'
#' @description
#' A function for calculating spatial autocorrelation using Moran's I.
#'
#' @param x 	    A vector of length \emph{n} representing the x coordinates
#'         (integer, consecutively numbered cells).
#' @param y 	    A vector of length \emph{n} representing the y coordinates
#'         (integer, consecutively numbered cells)
#' @param f       A vector which is the same length as \code{x} and \code{y}
#' @param lim1    Lower bound for first bin. Default is 1
#' @param lim2    Upper bound for first bin. Default is 2
#' @param dmax    Number of distance bins to examine. Bins are formed by annuli of gradually
#' increasing radii. Default is 10.
#'
#' @return A vector of Moran's I values for each distance bin.
#'
#' @examples
#' data(musdata)
#' coords<- musdata[,4:5]
#' mglm <- glm(musculus ~ pollution + exposure, "poisson", musdata)
#'
#' ac<-acfft(coords[,1],coords[,2],resid(mglm,type="pearson"),lim1=0,lim2=1)
#' ac
#'
#' @author Gudrun Carl

#' @export
acfft<-function(x,y,f,lim1=1,lim2=2,dmax=10){



  if(length(x)!=length(f)) stop("error in dimension")
  logic1<-identical(as.numeric(x),round(x,0))
  logic2<-identical(as.numeric(y),round(y,0))
  if(!logic1 | !logic2) stop("coordinates not integer")

  reslm<-f-mean(f)
  mi<-max(x)-min(x)+1
  mk<-max(y)-min(y)+1
  n<-max(mi,mk)
  n2<-n*n
  Ares<-matrix(0,n,n)
  mask<-matrix(0,n,n)
  for(i in 1:length(x)){
    kx<-x[i]-min(x)+1
    ky<-y[i]-min(y)+1
    Ares[ky,kx]<-reslm[i]
    mask[ky,kx]<-1}
  filter1<-matrix(0,n,n)
  filter1[1,1]<-1
  leng<-length(reslm)
  ne<-convolve(convolve(Ares,filter1),Ares)[1,1]/leng
  n3<-3*n
  Ares0<-matrix(0,n3,n3)
  Ares1<-matrix(0,n3,n3)
  n1<-n+1
  nn<-n+n
  Ares0[n1:nn,n1:nn]<-Ares[1:n,1:n]
  Ares1[1:n,1:n]<-Ares[1:n,1:n]
  maske0<-matrix(0,n3,n3)
  maske1<-matrix(0,n3,n3)
  maske0[n1:nn,n1:nn]<-mask[1:n,1:n]
  maske1[1:n,1:n]<-mask[1:n,1:n]
  nx<-rep(1:n3,n3)
  ny<-as.numeric(gl(n3,n3))

  gr<-lim1
  gr1<-lim2
  h<-n*n3+n+1
  corr<-rep(0,dmax)
  kk<-0
  while(kk<dmax){
    kk<-kk+1
    filter<-matrix(0,n3,n3)
    for(i in 1:(n3*n3)) {
      d<-sqrt((nx[h]-nx[i])^2+(ny[h]-ny[i])^2)

      if(lim1!=0 & d>=gr & d<gr1) filter[ny[i],nx[i]]<-1
      if(lim1==0 & d>gr & d<=gr1) filter[ny[i],nx[i]]<-1}

    sum<-convolve(convolve(maske0,filter),maske1)[1,1]
    za<-convolve(convolve(Ares0,filter),Ares1)[1,1]/sum
    corr[kk]<-za/ne
    gr<-gr+(lim2-lim1)
    gr1<-gr1+(lim2-lim1)
  }
  corr<-as.vector(corr)
  return(corr)
}

