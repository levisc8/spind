#'summary.WRM
#'
#'@description Provides \code{summary} methods for WRM model objects
#'
#'@param object A model of class \code{WRM}
#'@param ... Other arguments passed to \code{summary}
#'@return Prints results of model

#'@export

summary.WRM<-function (object,...) {
  cat("\n","Call:","\n")
  print(object$call)
  cat("\n","Pearson Residuals:","\n")
  print(summary(object$resid))
  family<-object$family
  b<-object$b
  s.e.<-object$s.e.
  z<-object$z
  p<-object$p
  it<-object$it
  n<-object$n
  n.eff<-object$n.eff
  AIC<-object$AIC
  beta<-cbind(b,s.e.,z,p)
  if(family=="gaussian")
    colnames(beta) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
  if(family=="binomial" | family=="poisson")
    colnames(beta) <- c("Estimate", "Std.Err", "z value", "Pr(>|z|)")
  cat("---","\n","Coefficients:","\n")
  printCoefmat(beta)
  cat("---","\n","Number of observations n: ",n, ",  n.eff: ",
      n.eff,",  AIC: ",AIC,"\n" )
  cat("\n","Number of iterations: ",it,"\n")
  cat("---","\n")
  ac0<-object$ac.glm
  acw<-object$ac
  cat("Autocorrelation of glm.residuals","\n")
  print(ac0)
  cat("Autocorrelation of wavelet.residuals","\n")
  print(acw)
}
