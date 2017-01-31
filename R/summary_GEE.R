summary_gee<-function(model){
   cat("\n","Call:","\n")
   print(model$call)
   family<-model$family
   b<-model$b
   s.e.<-model$s.e.
   z<-model$z
   p<-model$p
   it<-model$it
   n<-model$n
   n.eff<-model$n.eff
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
   acg<-model$ac
   cat("Autocorrelation of GLM residuals","\n")
   print(ac0)
   cat("\n","Autocorrelation of GEE residuals","\n")
   print(acg)
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
#          W.ac       working autocorrelation matrix
#          QIC        quasi-information criterion
#          if plot or graph is true:
#          ac.glm     autocorrelation of glm.residuals
#          ac         autocorrelation of gee.residuals
