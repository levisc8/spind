#' Stepwise model selection for WRMs and GEEs
#'
#' @description A function that provides stepwise model selection using
#' AIC for Wavelet revised models (WRMs) and QIC for Generalized estimating
#' equations (GEEs).
#'
#' @param object A model fit of class \code{WRM} or \code{GEE}.
#' @param data The data used to fit the specified model.
#' @param scope A list with components \code{upper} and \code{lower}
#' that specify the upper and lower boundaries of the stepwise
#' model selection procedure. Default is NULL, which will step backwards through
#' all predictor variables.
#' @param direction Only \code{backward} is supported right now. The operation
#' will begin with the parameter \code{upper} and work backwards until it reaches
#' parameter \code{lower} (specified in \code{scope}). \code{forward} or \code{both} may be added later.
#'
#' @return A data frame containing model formulae. WRM models will have AIC
#' Log-likelihood, and delta-AIC information. GEE models will have QIC,
#' Quasi-likelihood, and delta-QIC information.
#'
#' @seealso \code{\link{aic.calc}},\code{\link{qic.calc}}
#'
#' @author Sam Levin
#'@examples
#' data(musdata)
#' coords<- musdata[,4:5]
#'
#' xgee<-GEE(musculus ~ pollution + exposure, "poisson", musdata,
#'       coord=coords, corstr="fixed", plot=TRUE,scale.fix=FALSE)
#'
#' mwrm<-WRM(musculus ~ pollution + exposure, "poisson", musdata,
#' coord=coords, level=1, plot=TRUE)
#'
#' mod.select.gee<-step.spind(xgee,musdata)
#' mod.select.wrm<-step.spind(mwrm,musdata)
#'
#' @export
#'

step.spind<-function(object,data,scope=NULL,direction="backward"){

  model<-class(object)
  family<-object$family
  formula<-terms(object$formula)
  vars<-all.vars(object$formula)
  coord<-object$coord
  scale.fix<-object$scale.fix
  upper<-scope$upper
  lower<-scope$lower
  if(!is.null(scope)){
    upper.int<-match(upper,vars)
    lower.int<-match(lower,vars)
  }else{
   upper.int<-length(vars)
   lower.int<-2
  }



  ic<-list(formula=rep(NA,upper.int),
           inf.crit=rep(NA,upper.int),
           log.lik=rep(NA,upper.int))

  if(model=="WRM"){
    for(i in upper.int:lower.int){

      newmodel<-drop.terms(formula,i:upper.int,keep.response = TRUE)
      mwrm<-WRM(newmodel,family,data,coord,object$level,
                object$wavelet,object$wtrafo,pad=list(object$padform,
                                                      object$padzone))

      aic<-aic.calc(mwrm$formula,mwrm$family,data,mwrm$fitted,mwrm$n.eff)

      ic$formula[i]<-deparse(newmodel)
      ic$inf.crit[i]<-aic$AIC
      ic$log.lik[i]<-aic$loglik

    }
    ic$delta<-ic$inf.crit-min(ic$inf.crit,na.rm=TRUE)
    names(ic)<-c("Formula","AIC","logLik","Delta.AIC")
  }
  if(model=="GEE"){
    if(!scale.fix){
      scale.fix<-TRUE
      warning("Warning: scale parameter is now fixed")
    }
    for(i in upper.int:lower.int){
      newmodel<-drop.terms(formula,i:upper.int,keep.response = TRUE)
      mgee<-GEE(newmodel,family,data,coord,object$corstr,object$cluster)

      ic$formula[i]<-deparse(newmodel)
      ic$inf.crit[i]<-mgee$QIC
      ic$log.lik[i]<-mgee$QLik
    }
    ic$delta<-ic$inf.crit-min(ic$inf.crit,na.rm=T)
    names(ic)<-c("Formula","QIC","Quasi.Lik","Delta")
  }

  ic<-data.frame(ic,stringsAsFactors = FALSE)
  ic<-ic[complete.cases(ic),]
  ic<-ic[order(ic$Delta),]
  return(ic)
}



