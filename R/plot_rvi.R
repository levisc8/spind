#' Plot methods for scale-dependent relative variable importance
#'
#' @description The plot is based on the results provided by
#' multi-model inference for wavelet multiresolution regression (mmiWMRR).
#' Based on the results of this approach, the relative
#' importance of predictors are calculated. The scale-dependent results
#' are graphically displayed.
#'
#'
#' @param A      result (a value from mmiWMRR)
#' @param level  level (a value from mmiWMRR)

plot.rvi<-function(A,level){

  nvar<-dim(A)[2]-6
  klimitscale<-dim(A)[3]
  ip<-dim(A)[1]
  WeightSums<-matrix(NA,nvar,klimitscale)
  for (kscale in 1:klimitscale){
    for(kvar in 2:(nvar+1)){
      for (i in 1: ip){
        if(!is.na(A[i,kvar,kscale])) A[i,kvar,kscale]<-A[i,(nvar+6),kscale]
      }
    }
    B<-A[1:ip,2:(nvar+1),kscale]
    WeightSums[,kscale]<-colSums(B,na.rm=TRUE)
  } # kscale

  vec<-1:nvar

  plot(level,WeightSums[1,],type="b",ylim=c(0,2),
       xlim=c(min(level),max(level)),ylab="Relative Variable Importance",
       xlab="Level", pch=2,lty=vec[1],lwd=2)
  for (kvar in 2:nvar) {
    points(level,WeightSums[kvar,],type="b",pch=kvar+1,
           lty=vec[kvar],lwd=2)
  }

  leg<-dimnames(mmi$res)[[2]][vec+1]
  v<-2:(nvar+1)
  legend(2,2,leg,pch=v,lty=vec,lwd=2)
}
