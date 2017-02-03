#' predict.GEE
#'
#' @description
#' Predicts values using a \code{GEE} model.
#'
#' @param object   A model object of class \code{GEE} to be used
#' for making predictions
#' @param newdata  A data frame containing variables with which to predict
#'
#' @return   A vector of predicted values
#' @examples
#' data(musdata)
#' coords<-musdata[,1:2]
#' mgee<-GEE(musculus ~ pollution + exposure,'poisson',musdata,
#'           coords=coords,corstr="fixed",plot=TRUE)
#' pred<-predict(mgee,musdata)
#'
predict.GEE<-function(object,newdata){


  data<-newdata
  formula<-object$formula
  family<-object$family
  b<-object$b

  x.matrix<-object.matrix(formula,data)
  fitted<-x.matrix%*%b
  fitted<-as.vector(fitted)
  if(family=="poisson") fitted<-exp(fitted)
  if(family=="binomial") fitted<-exp(fitted)/(1+exp(fitted))

  return(fitted)
}
