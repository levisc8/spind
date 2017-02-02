#' @export

summary_gee<-function(model,printAutoCorMat=FALSE){
#' @title Summarize the output from\code{GEE}
#' @description Returns summary of GEE parameter estimates and autocorrelations
#' of the residuals.
#' @param model An object of class \code{gee}
#' @param printAutoCorMat A logical indicating whether to print the
#' autocorrelation matrix
#'
#'@return Prints model details, parameter estimates, and autocorrelation values
#'for the first 10 distance bins. Additionally, prints the autocorrelation matrix
#'if \code{printAutoCorMat} = TRUE
#'
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

   if(printAutoCorMat) print(model$Mat.ac)
}





  # Gee function returns....
  # An object of class "gee" is a list containing the following components:
  #          call       call
  #          formula    formula
  #          family     family
  #          b          estimate of regression parameters
  #          s.e.       standard errors
  #          z          z values (or corresponding values for statistics)
#          p          probabilities
#          scale      scale parameter (dispersion parameter)
#          fitted     fitted values
#          resid      normalized Pearson residuals
#          w.ac       working autocorrelation parameters
#          Mat.ac       working autocorrelation matrix
#          QIC        quasi-information criterion
#          if plot or graph is true:
#          ac.glm     autocorrelation of glm.residuals
#          ac.gee         autocorrelation of gee.residuals
