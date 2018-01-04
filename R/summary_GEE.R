#' @title Summarize the output from \code{GEE}
#' @description Returns summary of GEE parameter estimates and autocorrelations
#' of the residuals.
#' @param object An object of class \code{GEE}
#' @param ... Additional parameters to be passed \code{summary}.
#' @param printAutoCorPars A logical indicating whether to print the
#' working autocorrelation parameters
#'
#'@return Prints model details, parameter estimates, and autocorrelation values
#'for the first 10 distance bins. Additionally, if \code{printAutoCorPars} = TRUE,
#'prints working autocorrelation parameters used in the model.
#'
#' @author Sam Levin
#' @importFrom stats printCoefmat
#' @export
summary.GEE<-function(object,...,printAutoCorPars=TRUE){

   cat("\n","Call:","\n")
   print(object$call)
   family<-object$family
   b<-object$b
   s.e.<-object$s.e.
   z<-object$z
   p<-object$p
   QIC<-object$QIC
   beta<-cbind(b,s.e.,z,p)
   if(family=="gaussian")
      colnames(beta) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
   if(family=="binomial" | family=="poisson")
      colnames(beta) <- c("Estimate", "Std.Err", "z value", "Pr(>|z|)")
   cat("---","\n","Coefficients:","\n")
   stats::printCoefmat(beta)
   cat("---","\n","QIC: ",QIC,"\n" )
   cat("---","\n")
   ac0<-object$ac.glm
   acg<-object$ac.gee
   cat("Autocorrelation of GLM residuals","\n")
   print(ac0)
   cat("\n","Autocorrelation of GEE residuals","\n")
   print(acg)

   if(printAutoCorPars & object$corstr != "independence"){
     cat('---','\n','Autocorrelation parameters from ',
         object$corstr," model",'\n')
     print(object$Mat.ac)
   }
}
