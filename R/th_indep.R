#'th.indep
#'
#'Threshold-independent metrics
#'
#'@description Calculates threshold-independent metrics for spatial models
#'(AUC,ROC,max-TSS). It is a 2-dimensional analysis taking the grid
#'structure of data set into account. Grid cells are assumed to be square.
#'
#'@param data A data frame or matrix with two columns. The first column
#'should contain actual presence/absence data (binary, 0 or 1) and the
#'second column should contain model predictions of probability of
#'occurence (numeric, between 0 and 1).
#'@param coord A data frame or matrix with two columns containing x,y
#'coordinates for each actual and predicted value. Coordinates must be
#'integer and consecutively numbered.
#'@param spatial A logical value indicating whether spatial corrected
#'indices (rather than classical indices) should be computed
#'@param plot.ROC A logical indicating whether the ROC should be plotted
#'
#'@return A list with the following components:
#'\describe{
#'  \item{\code{AUC}}{Area under curve}
#'  \item{\code{TSS}}{Maximum TSS value}
#'  \item{\code{sensitivity}}{Sensitivity}
#'  \item{\code{Specificity}}{Specificity}
#'}
#'@author
#'Gudrun Carl
#'
#'@examples
#'data(hook)
#'data<-hook[,1:2]
#'coord<-hook[,3:4]
#'si2<-th.indep(data,coord,spatial=TRUE)
#'si2$AUC
#'si2$TSS
#'
#'@import splancs
#'@export


th.indep<-function(data,coord,spatial=TRUE,plot.ROC=TRUE){

  if(dim(data)[1]!=dim(coord)[1]) stop("error in dimension")

  if(spatial){
    y<-adjusted.actuals(data,coord)
    split<-4
  }
  if(!spatial){
    y<-data[,1]
    split<-2
  }

  pi<-data[,2]
  n<-length(pi)
  o<-order(pi)
  piord<-pi[o]
  cutoff<-c(1.1,piord[n:1])

  sensitivity<-rep(NA,n)
  specificity<-rep(NA,n)
  splitlevel<-matrix(NA,split,2)
  splitlevely<-matrix(NA,split,2)

  for(i.n in 1:n){
    thresh<-cutoff[i.n]

    if(split==4) { # split=4
      # split into 4 classes
      # for split=4 -> 4x4 conf.matrix
      lower.split<-thresh/2
      upper.split<-(1+thresh)/2
      # 0 ... lower.split
      # lower.split ... thresh
      # thresh ... upper.split
      # upper.split ... 1
      splitlevel[1,]<-c(0,lower.split)
      splitlevel[2,]<-c(lower.split,thresh)
      splitlevel[3,]<-c(thresh,upper.split)
      splitlevel[4,]<-c(upper.split,1)
      splitlevely[1,]<-c(0,0.25)
      splitlevely[2,]<-c(0.25,0.5)
      splitlevely[3,]<-c(0.5,0.75)
      splitlevely[4,]<-c(0.75,1)
      pipos<-matrix(0,n,4)
      ypos<-matrix(0,n,4)
      for(k in 1:n){
        for(ksp in 1:3){
          if(splitlevel[ksp,1] <= pi[k] & pi[k] < splitlevel[ksp,2]) pipos[k,ksp]<-1
          if(splitlevely[ksp,1] <= y[k] & y[k]  < splitlevely[ksp,2]) ypos[k,ksp]<-1
        }
        if(splitlevel[4,1] <= pi[k] & pi[k] <= splitlevel[4,2]) pipos[k,4]<-1
        if(splitlevely[4,1] <= y[k] & y[k]  <= splitlevely[4,2]) ypos[k,4]<-1
      }
      cm<-matrix(0,4,4)
      for(k in 1:n){
        i<-which(ypos[k,]==1)
        j<-which(pipos[k,]==1)
        cm[i,j]<-cm[i,j]+1
      }
      # inversion in relation to the second. diagonal
      cm<-matrix(rev(as.vector(cm)),4,4, byrow = TRUE)
      # weights w
      w<-matrix(NA,4,4)
      for (i in 1:4){
        for (j in 1:4){
          w[i,j]<-ifelse(abs(i-j)<2,1,0)
        }}
      n<-sum(cm)
      sensitivity[i.n]<-sum(w[,1:2]*cm[,1:2])/sum(cm[,1:2])
      specificity[i.n]<-sum(w[,3:4]*cm[,3:4])/sum(cm[,3:4])

    } # split=4

    if(split==2) { # split=2
      # split into 2 classes
      # for split=2 -> 2x2 conf.matrix
      # 0 ... thresh
      # thresh ... 1
      splitlevel[1,]<-c(0,thresh)
      splitlevel[2,]<-c(thresh,1)
      pipos<-matrix(0,n,2)
      ypos<-matrix(0,n,2)
      for(k in 1:n){
        if(splitlevel[1,1] <= pi[k] & pi[k] < splitlevel[1,2]) pipos[k,1]<-1
        if(splitlevel[1,1] <= y[k]  & y[k]  < splitlevel[1,2]) ypos[k,1]<-1
        if(splitlevel[2,1] <= pi[k] & pi[k] <= splitlevel[2,2]) pipos[k,2]<-1
        if(splitlevel[2,1] <= y[k]  & y[k]  <= splitlevel[2,2]) ypos[k,2]<-1
      }
      cm<-matrix(0,2,2)
      for(k in 1:n){
        i<-which(ypos[k,]==1)
        j<-which(pipos[k,]==1)
        cm[i,j]<-cm[i,j]+1
      }
      # inversion in relation to the second. diagonal
      cm<-matrix(rev(as.vector(cm)),2,2, byrow = TRUE)
      # weights w
      w<-matrix(NA,2,2)
      for (i in 1:2){
        for (j in 1:2){
          if(split==2) w[i,j]<-ifelse(abs(i-j)==0,1,0)  # w = identity matrix
        }}
      n<-sum(cm)
      sensitivity[i.n]<-sum(w[,1]*cm[,1])/sum(cm[,1])
      specificity[i.n]<-sum(w[,2]*cm[,2])/sum(cm[,2])

    } # split=2

  }  # i.n loop

  # TSS
  n<-length(pi)
  sensitivity[is.na(sensitivity)]<-0
  specificity[is.na(specificity)]<-0
  TSS<-max(sensitivity+specificity)-1

  # ROC
  if(plot.ROC){
    plot(1-specificity,sensitivity,
         type="l",xlim=c(0,1),ylim=c(0,1),
         xlab="1 - specificity",ylab="sensitivity",main="ROC")
    points(c(0,1),c(0,1),type="l")
  }

  # AUC
  dat<-cbind(c(1-specificity,1,1,0), c(sensitivity, 1, 0, 0))
  AUC<-splancs::areapl(dat)

  list(AUC=AUC,TSS=TSS,sensitivity=sensitivity,specificity=specificity)

}
