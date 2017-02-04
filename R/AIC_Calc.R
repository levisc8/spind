#'
#'
#'
#'
#'

aic.calc<-function(formula,family,data,mu,n.eff=NULL){
  ###############################################################################
  # Description
  # A function for calculating log likelihood and AIC values.
  # Arguments:
  # formula  with specified notation according to names in data frame
  # family   "gaussian", "binomial"(binary) or "poisson"
  # data     data frame
  # mu       fitted values
  # n.eff    effective number of observations
  #
  # value:   loglik     log likelihood
  #          df         degrees of freedom
  #          AIC        Akaike information criterion
  #          AICc       AIC with a correction for finite sample sizes
  ###############################################################################

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
  AIC  <- round(AIC,2)
  AICc  <- round(AICc,2)

  return(list(loglik=loglik,df=K,AIC=AIC,AICc=AICc))
}
