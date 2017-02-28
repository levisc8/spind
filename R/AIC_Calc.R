#'Aikake Information Criterion with correction for sample size
#'
#' @description Calculates AIC and AICc
#'
#' @param formula A model formula
#' @param family \code{gaussian}, \code{binomial}, or \code{poisson}
#' @param data A data frame
#' @param mu Fitted values from a model
#' @param n.eff Effective number of observations. Default is NULL
#'
#' @return A list with the following components
#' \describe{
#'   \item{\code{loglik}}{Log likelihood of the model}
#'   \item{\code{df}}{Degrees of freedom}
#'   \item{\code{AIC}}{AIC score for the specified model}
#'   \item{\code{AICc}}{AIC score corrected for small sample sizes}
#' }
#'
#' @author Gudrun Carl
#'
#' @examples
#' data(musdata)
#' coords<- musdata[,4:5]
#' mglm <- glm(musculus ~ pollution + exposure, "poisson", musdata)
#'
#' aic<-aic.calc(musculus ~ pollution + exposure, "poisson", musdata,
#'              mglm$fitted)
#' aic$AIC
#'
#'
#' @export

aic.calc<-function(formula,family,data,mu,n.eff=NULL){

  X<-model.matrix(formula,data)
  if(is.vector(model.frame(formula,data)[[1]])){
    y<-model.frame(formula,data)[[1]]
    ntr<-1
  }
  if(family=="binomial" & is.matrix(model.frame(formula,data)[[1]])){
    y<-model.frame(formula,data)[[1]][,1]
    ntr<-model.frame(formula,data)[[1]][,1] +
      model.frame(formula,data)[[1]][,2]
  }

  n <- dim(X)[1]
  nvar<-dim(X)[2]
  if(family=="gaussian"){
    sos<-sum((y-mu)^2)                    # sum of squares
    sigma2<- sos/n                        # variance
    #loglik<- -n/2*log(2*pi) - n*log(sigma2^(1/2)) - 1/(2*sigma2)*sos
    #loglik<- -n/2*log(2*pi) - n*log(sigma2^(1/2)) - n/2
    #loglik<- -n/2*log(2*pi) - n/2*log(sigma2) - n/2
    #loglik<- -n/2*(log(2*pi) + log(sigma2) +1)
    loglik<- -n/2*(log(2*pi*sigma2) +1)   # log likelihood
    if(!is.null(n.eff)) loglik<- -n.eff/2*(log(2*pi*sigma2) +1)
    K<-nvar+1            # nvar (number of pred.+interc.) +1 (for variance)
  }
  if(family=="binomial"){                         # choose= Binomialkoeff.
    loglik<- sum(y*log(mu/(1-mu)) + ntr*log(1-mu) + log(choose(ntr,y)) )
    K<-nvar              # without variance (since var. is function of mean)
  }
  if(family=="poisson"){
    # loglik<- sum(y*log(mu)-mu)  # useful for delta in multimodel inference
    loglik<- sum(y*log(mu)-mu) - sum(log(factorial(y))) # factorial=Fakult?t
    K<-nvar              # without variance (since var. is function of mean)

  }

  AIC  <- -2* loglik + 2*K              # AIC
  AICc <- -2* loglik + 2*K + 2*K*(K+1)/(n-K-1)    # AICc


  return(list(loglik=loglik,df=K,AIC=AIC,AICc=AICc))
}
