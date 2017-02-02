#'@title Adjusted actual values
#'
#'@description ASK GUDRUN WHAT THIS ACTUALLY DOES
#'
#'@param data a dataframe or matrix containing actual presence/absence (binary,
#'0 or 1) values in 1st column and predicted values (numeric between 0 and 1)
#'in 2nd column.
#'@param coord a matrix of two columns of the same length providing integer,
#'consecutively numbered coordinates for each occurence and prediction in
#'\code{data}.
#'@param plot.maps A logical indicating whether maps should be plotted.
#'Default is FALSE.
#'
#'@return a vector of adjusted actual values
#'
#'@examples {
#' data(hook)
#' data<- hook[,1:2]
#' coord<- hook[,3:4]
#' # plot maps
#' aa<-adjusted.actuals(data,coord,plot.maps=TRUE)}
#'
#' @export

adjusted.actuals<-function(data,coord,plot.maps=FALSE){

  x<-coord[,1]
  y<-coord[,2]
  fb<-data[,1]
  fa<-data[,2]

  if(length(x)!=length(fa)) stop("coordinates[,1] and coordinates[,2] have different dimensions")
  logic1<-identical(as.numeric(x),round(x,0))
  logic2<-identical(as.numeric(y),round(y,0))
  if(!logic1 | !logic2) stop("coordinates not integer")

  fbs<-fb
  ac01a<-acfft(x,y,fa,lim1=0,lim2=1,dmax=1)
  ac01b<-acfft(x,y,fb,lim1=0,lim2=1,dmax=1)
  if(ac01a>0.05 & ac01b>0.05){
    ac01<-ac01a-ac01b
    if(ac01>0.02){
      alpha<-ac01
      D<-as.matrix(dist(coord))
      R<-alpha^D
      spatial.W<-R^3
      ac01s<-acfft(x,y,fbs,lim1=0,lim2=1,dmax=1)
      while(ac01a>ac01s ){
        fbs<-spatial.W%*%fbs
        ac01s<-acfft(x,y,fbs,lim1=0,lim2=1,dmax=1)
      }
      fbs<-fbs-min(fbs)
      fbs<-fbs/max(fbs)
    }}

  if (plot.maps){
    colours<-rev(gray(3:32/32 -3/32 ))
    a<-lattice::levelplot(fa~x+y,
                 col.regions=colours,
                 colorkey=FALSE,
                 scales = list(draw=FALSE),
                 xlab="",ylab="",main="predictions"
    )
    b<-lattice::levelplot(fb~x+y,
                 col.regions=colours,
                 colorkey=FALSE,
                 scales = list(draw=FALSE),
                 xlab="",ylab="",main="actuals"
    )
    c<-lattice::levelplot(fbs~x+y,
                 col.regions=colours,
                 colorkey=list(space="bottom"),
                 scales = list(draw=FALSE),
                 xlab="",ylab="",main="adjusted actuals"
    )
    tp <- trellis.par.get()
    lattice::trellis.par.set(list(axis.line = list(col = "transparent")))
    print(a,position=c(0.1,0.09,0.92,0.98),split=c(1,1,2,2),more=TRUE)
    print(b,position=c(0.12,0.09,0.94,0.98),split=c(2,1,2,2),more=TRUE)
    print(c,position=c(0.12,0,0.94,1.02),split=c(2,2,2,2),more=FALSE)
  } # plot

  fbs<-as.vector(fbs)
  return(fbs)
}
