#'
#' @title Multi-model inference for GEE models
#'
#' @description
#' mmiGEE is a multimodel inference approach evaluating the relative
#' importance of predictors used in \code{\link{GEE}}.
#'
#'  @details It performs automatically
#' generated model selection and creates a model
#' selection table according to the approach of multi-model inference
#' (Burnham & Anderson, 2002). QIC is used to obtain model
#' selection weights and to rank the models. Moreover, mmiGEE calculates relative
#' variable importance of a given model.
#' Finally, this function requires that \strong{all predictor variables
#' be continuous}.
#'
#'
#' @details Calculates the relative importance of each variable
#' using multi-model inference methods in a Generalized Estimating Equations
#' framework implemented in \code{GEE}.
#'
#' @param object A model of class \code{GEE}.
#' @param data A data frame or set of vectors of equal length.
#' @param trace A logical indicating whether or not to print results to console.
#'
#' @return  \code{mmiGEE} returns a list containing the following elements
#' \describe{
#'   \item{\code{result}}{A matrix containing slopes, degrees of freedom, quasilikelihood,
#'          QIC, delta, and weight values for the set of candidate models.
#'          The models are ranked by QIC.}
#'   \item{\code{rvi}}{A vector containing the relative importance of each variable
#'          in the regression.}
#'
#'}
#'
#' @references
#' Burnham, K.P. & Anderson, D.R. (2002) Model selection and
#' multimodel inference. Springer, New York.
#'
#' Carl G & Kuehn I, 2007. Analyzing Spatial Autocorrelation in Species
#' Distributions using Gaussian and Logit Models, Ecol. Model. 207, 159 - 170
#'
#' @seealso \code{\link{GEE}}, \code{\link{qic.calc}}, \pkg{MuMIn}
#'
#' @author Gudrun Carl, Sam Levin
#'
#' @examples
#' # data (for demonstration only)
#' library(MASS)
#' data(birthwt)
#' # impose an artificial (not fully appropriate) grid structure
#' x<-rep(1:14,14)
#' y<-as.integer(gl(14,14))
#' coords<-cbind(x[-(190:196)],y[-(190:196)])
#'
#'\dontrun{
#'
#' formula<-formula(low ~ race + smoke +  bwt)
#'
#' mgee<-GEE(formula, family = "gaussian", data = birthwt,
#'          coord=coords, corstr="fixed",scale.fix=TRUE)
#'
#' mmi<-mmiGEE(mgee,birthwt)
#'}
#' @importFrom stats update reformulate drop.terms model.matrix as.formula terms
#' @importFrom rje powerSetMat
#' @export
#'
#'

mmiGEE<-function(object,data, trace = FALSE){


  if(!inherits(object, "GEE")) stop("Input model is not of class 'GEE'")

  family<-object$family
  formula<-object$formula

  coord<-object$coord
  scale.fix<-object$scale.fix
  corstr<-object$corstr
  cluster<-object$cluster
  moran.params<-object$moran.params

  if(!scale.fix){
    scale.fix<-TRUE
    message("Scale parameter is now fixed")
  }

  # Parameter: varnames, p
  X<-stats::model.matrix(formula,data)
  if(dimnames(X)[[2]][1]!="(Intercept)") {
    formula <- stats::update(formula, ~ . + 1)
    X<-stats::model.matrix(formula,data)
  }
  nvar<-dim(X)[2]
  varnames<-dimnames(X)[[2]][-1]
  p<-dim(X)[2]-1 # nvar-1 (without intercept)

  pset <- rje::powerSetMat(p)
  ip <- dim(pset)[1]
  t <- stats::terms(formula)
  # Run every model and calculate QIC (multimodel inference)
  coef.vec<-matrix(NA,ip,nvar)
  df<-rep(NA,ip)
  Qlik<-rep(NA,ip)
  QIC<-rep(NA,ip)
  for (i in 1:ip) {
    if(sum(pset[i,])!=0 & sum(pset[i,])!=p){
      t1 <- stats::drop.terms(t,which(pset[i,]==0), keep.response = TRUE)
      formula1<- stats::reformulate(attr(t1, "term.labels"), formula[[2]])
      formulae<-formula1
    }
    if(sum(pset[i,])==p) formulae<-formula
    if(sum(pset[i,])==0) formulae<-stats::as.formula(paste(formula[[2]],"~1"))

    m0<-suppressWarnings({
      GEE(formulae,family,data,coord,corstr=corstr,
          cluster=cluster,moran.params=moran.params,
          scale.fix=scale.fix)
    })

    kv<-c(1,which(pset[i,]==1)+1)
    coef.vec[i,kv]<-m0$b
    Xe<-stats::model.matrix(formulae,data)
    nvare<-dim(Xe)[2]
    K<-nvare
    if(family=="gaussian") K<-K+1
    df[i]<-K
    Qlik[i]<-m0$QLik
    QIC[i]<-m0$QIC
  }
  result<-cbind(round(coef.vec,5),df,round(Qlik,3),round(QIC,1))

  # Calculate delta and weight (multimodel inference)
  delta<-QIC-min(QIC) # = delta
  weight<-exp(-delta/2)/sum(exp(-delta/2)) # = weight

  # Print results
  result<-cbind(result,round(delta,2),round(weight,3))
  ord<-order(delta)
  res<-result[ord,]
  dimnames(res)[[1]]<-ord
  dimnames(res)[[2]]<-c("(Int)",varnames,
                        "df","QLik","QIC","delta","weight")

  if(trace) {
    cat("\n","Model selection table:","\n","\n")
    print(res,na.print = "")
  }
  nrowA<-dim(res)[1]
  ncolA<-dim(res)[2]

  nvar<-dim(res)[2]-6
  leg<-dimnames(res)[[2]][2:(nvar+1)]

  A<-matrix(NA,nrowA,ncolA)
  A<-res

  # Plot: relative variable importance

  ip<-dim(A)[1]

  WeightSums<-rep(NA,nvar)
  for(kvar in 2:(nvar+1)){
    for (i in 1: ip){
      if(!is.na(A[i,kvar])) A[i,kvar]<-A[i,(nvar+6)]
    }
  }
  B<-A[1:ip,2:(nvar+1)]
  WeightSums<-colSums(B,na.rm=TRUE)

  names(WeightSums)<-leg
  if(trace){
    cat("\n","---","\n","Relative variable importance:","\n","\n")
    print(WeightSums)
  }
  fit<-list(result=res,rvi=WeightSums)
  fit
}


