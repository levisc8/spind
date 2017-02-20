#' @title Multi-model inference for wavelet multiresolution regression
#'
#'
#' @description
#' mmiWMRR is a multimodel inference approach evaluating the relative
#' importance of predictors used in scaleWMRR. It performs automatically
#' generated model selection and creates a model
#' selection table according to the approach of multi-model inference
#' (Burnham & Anderson, 2002). The analysis is carried out for scale-specific
#' regressions, i.e. where scaleWMRR as
#' modelling function is supported. AIC is used to obtain model
#' selection weights and to rank the models.
#'
#'
#' @param formula  with specified notation according to names in data frame
#' @param family   "gaussian", "binomial"(binary) or "poisson"
#' @param data     data frame
#' @param coord    corresponding coordinates which have to be integer
#' @param scale    0 or higher integers possible (limit depends on sample size)
#' @param detail   by detail components only
#' @param wavelet  type of wavelet: "haar" or "d4" or "la8"
#' @param wtrafo   type of wavelet transform: "dwt" or "modwt"
#' @param n.eff    a numeric value of effective sample size
#'
#' @return  mmiWMRR returns a list containing the following elements
#' \describe{
#'   \item{\code{result}}{a matrix containing slopes, degrees of freedom, Likelihood,
#'          AIC, delta and weight values for the set of candidate models.
#'          The models are ranked by Akaike weights.}
#'   \item{\code{level}}{an integer corresponding to scale}
#'
#'}
#' @seealso  AIC, \pkg{MuMIn}.
#'
#' @examples
#'
#' data(carlinadata)
#' coord<- carlinadata[,4:5]
#'
#' # scale-specific regressions for detail components
#' # ranked by multimodel inference
#' A<-array(NA,c(4,8,4))
#' level<-rep(NA,4)
#' for (i in 1:4) {
#'   mmi<- mmiWMRR(carlina.horrida ~ aridity + land.use,"poisson",
#'                 carlinadata,coord,scale=i,detail=TRUE,wavelet="d4")
#'   A[,,i]<-mmi$result
#'   level[i]<-mmi$level
#' }
#'
#' # Plot: scale-dependent relative variable importance
#' plot.rvi(A,level)
#'
#' @references
#' Burnham, K.P. & Anderson, D.R. (2002) Model selection and
#' multimodel inference. Springer, New York.
#'
#' Carl G, Doktor D, Schweiger O, Kuhn I (2016)
#' Assessing relative variable importance across different spatial
#' scales: a two-dimensional wavelet analysis.
#' Journal of Biogeography 43: 2502-2512.
#'
#' @export

mmiWMRR<-function(formula,family,data,coord,scale,detail=TRUE,
                  wavelet="haar",wtrafo="dwt",n.eff=NULL){

  # Parameter: varnames, p
  X<-model.matrix(formula,data)
  if(dimnames(X)[[2]][1]!="(Intercept)") {
    formula<-update(formula, ~ . + 1)
    X<-model.matrix(formula,data)
  }
  nvar<-dim(X)[2]
  varnames<-dimnames(X)[[2]][-1]
  p<-dim(X)[2]-1 # nvar-1 (without intercept)

  pset<-rje::powerSetMat(p)
  ip<-dim(pset)[1]
  t<-terms(formula)
  # Run every model and calculate AIC (multimodel inference)
  coef.vec<-matrix(NA,ip,nvar)
  df<-rep(NA,ip)
  loglik<-rep(NA,ip)
  AIC<-rep(NA,ip)
  for (i in 1:ip) {
    if(sum(pset[i,])!=0 & sum(pset[i,])!=p){
      t1<-drop.terms(t,which(pset[i,]==0), keep.response = TRUE)
      formula1<-reformulate(attr(t1, "term.labels"), formula[[2]])
      formulae<-formula1
    }
    if(sum(pset[i,])==p) formulae<-formula
    if(sum(pset[i,])==0) formulae<-as.formula(paste(formula[[2]],"~1"))
    m0<-scaleWMRR(formulae,family,data,coord,scale=scale,
                  detail=detail,wavelet=wavelet,wtrafo=wtrafo,
                  pad=list(padzone=1.1))
    if(!m0$converged)
      m0<-scaleWMRR(formulae,family,data,coord,scale=scale,
                    detail=detail,wavelet=wavelet,wtrafo=wtrafo,
                    b.ini=glm(formulae,family,data)$coef,
                    pad=list(padzone=1.1))
    kv<-c(1,which(pset[i,]==1)+1)
    coef.vec[i,kv]<-m0$b
    mu<-m0$fitted
    if(is.null(n.eff)) aic<-aic.calc(formulae,family,data,mu)
    if(!is.null(n.eff)) aic<-aic.calc(formulae,family,data,mu,n.eff)
    df[i]<-aic$df
    loglik[i]<-aic$loglik
    AIC[i]<-aic$AIC
  }
  result<-cbind(round(coef.vec,5),df,round(loglik,3),round(AIC,1))
  # Calculate delta and weight (multimodel inference)
  delta<-AIC-min(AIC) # = delta
  weight<-exp(-delta/2)/sum(exp(-delta/2)) # = weight
  # Print results
  result<-cbind(result,round(delta,2),round(weight,3))
  ord<-order(delta)
  res<-result[ord,]
  dimnames(res)[[1]]<-ord
  dimnames(res)[[2]]<-c("(Int)",varnames,
                        "df","logLik","AIC","delta","weight")
  if(!detail & scale>=1) scale<-scale-1
  cat("---","\n","Level = ",scale,"\n")
  print(res,na.print = "")
  fit<-list(result=res,level=scale)
  fit
}