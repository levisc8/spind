#' predict.GEE
#'
#' @description
#' Model predictions for \code{GEE}
#'
#' @param object   A model object of class \code{GEE} to be used
#' for making predictions
#' @param ... Other arguments to be passed to \code{predict}
#' @param newdata  A data frame containing variables to base the predictions on.
#'
#' @return   A vector of predicted values
#' @examples
#' data(musdata)
#' coords<-musdata[,4:5]
#' mgee<-GEE(musculus ~ pollution + exposure,'poisson',musdata,
#'           coord=coords,corstr="fixed",plot=TRUE)
#' pred<-predict(mgee,newdata=musdata)
#'
#'@author Gudrun Carl,
#'        Sam Levin
#'
#'@export
predict.GEE<-function(object,...,newdata){


  data<-newdata
  formula<-object$formula
  family<-object$family
  b<-object$b

  x.matrix<-model.matrix(formula,data)
  fitted<-x.matrix%*%b
  fitted<-as.vector(fitted)
  if(family=="poisson") fitted<-exp(fitted)
  if(family=="binomial") fitted<-exp(fitted)/(1+exp(fitted))

  return(fitted)
}
