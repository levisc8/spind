#' @title Summarize the output from\code{GEE}
#' @description Returns summary of GEE parameter estimates and autocorrelations
#' of the residuals.
#' @param model An object of class \code{gee}
#' @param printAutoCorPars A logical indicating whether to print the
#' autocorrelation matrix
#'
#'@return Prints model details, parameter estimates, and autocorrelation values
#'for the first 10 distance bins. Additionally, if \code{printAutoCorPars} = TRUE,
#'prints \describe{
#'
#'}
#'
#'
#' @export
summary_gee<-function(model,printAutoCorPars=FALSE){

   cat("\n","Call:","\n")
   print(model$call)
   family<-model$family
   b<-model$b
   s.e.<-model$s.e.
   z<-model$z
   p<-model$p
   QIC<-model$QIC
   beta<-cbind(b,s.e.,z,p)
   if(family=="gaussian")
      colnames(beta) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
   if(family=="binomial" | family=="poisson")
      colnames(beta) <- c("Estimate", "Std.Err", "z value", "Pr(>|z|)")
   cat("---","\n","Coefficients:","\n")
   printCoefmat(beta)
   cat("---","\n","QIC: ",QIC,"\n" )
   cat("---","\n")
   ac0<-model$ac.glm
   acg<-model$ac.gee
   cat("Autocorrelation of GLM residuals","\n")
   print(ac0)
   cat("\n","Autocorrelation of GEE residuals","\n")
   print(acg)

   if(printAutoCorPars&model$corstr!="independence"){
     cat('---','\n','Autocorrelation parameters from ',model$corstr," model",'\n')
     print(model$Mat.ac)
   }
}
