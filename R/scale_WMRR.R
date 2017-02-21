#'
#' @title Scaling by wavelet multiresolution regression (WMRR)
#'
#' @description
#' scaleWMRR performs a scale-specific regression based on a
#' wavelet multiresolution analysis. It  fits
#' generalized linear models while taking the two-dimensional grid structure of
#' datasets into account. The following error distributions (in
#' conjunction with appropriate link functions) are allowed: \code{gaussian},
#' \code{binomial}, or \code{poisson}. The model provides scale-specific
#' results for data sampled on a contiguous geographical area. The
#' dataset is assumed to be regular gridded and the grid cells are
#' assumed to be square. For the wavelet transforms, a function
#' of package 'waveslim' is used (Whitcher, 2005).
#'
#'
#' @param formula  With specified notation according to names in data frame.
#' @param family   \code{gaussian}, \code{binomial}, or \code{poisson}.
#' @param data     Data frame.
#' @param coord    Corresponding coordinates which have to be integer.
#' @param scale    0 (which is equivalent to GLM) or
#'          higher integers possible (limit depends on sample size).
#' @param detail   Remove smooth wavelets? If \code{TRUE}, only detail components are analyzed.
#' If set to \code{FALSE}, smooth and detail components are analyzed. Default is \code{TRUE}.
#' @param wavelet  Type of wavelet: \code{haar} or \code{d4} or \code{la8}
#' @param wtrafo   Type of wavelet transform: \code{dwt} or \code{modwt}.
#' @param b.ini    Initial parameter values. Default is \code{NULL}.
#' @param pad      A list of parameters for padding wavelet coefficients.
#' \itemize{
#'    \item{padform} - 0, 1, and 2 are possible.
#'     \code{padform} is automatically set to
#'     0 when either \code{level}=0 or
#'     the \code{formula} includes an intercept and has a non-\code{gaussian}
#'      \code{family}.
#'    \itemize{
#'      \item{0} - Padding with 0s.
#'      \item{1} - Padding with mean values.
#'      \item{2} - Padding with mirror values.
#'  }
#'    \item{padzone} - Factor for expanding the padding zone
#'}
#' @param control 	A list of parameters for controlling the fitting process.
#'    \itemize{
#'       \item{\code{eps}} - Positive convergence tolerance. Smaller values of
#'       \code{eps} provide better parameter estimates, but also reduce the probability
#'       of the iterations converging. In case of issues with convergence, test larger
#'       values of \code{eps}. Default is 10^-5.
#'       \item{\code{denom.eps}} - Default is 10^-20.
#'       \item{\code{itmax}} - Integer giving the maximum number of iterations.
#'       Default is 200.
#'}
#' @param moran.params    A list of parameters for calculating Moran's I.
#'   \itemize{
#'     \item\code{lim1} - Lower limit for first bin. Default is 0.
#'     \item\code{increment} - Step size for calculating Moran's I. Default is 1.
#'   }
#' @param plot     A logical value indicating whether to print parameter estimates
#' to the console
#'
#' @return  scaleWMRR returns a list containing the following elements
#' \describe{
#'   \item{\code{call}}{Model call}
#'   \item{\code{b}}{Estimates of regression parameters}
#'   \item{\code{s.e.}}{Standard errors of the parameter estimates}
#'   \item{\code{z}}{Z values (or corresponding values for statistics)}
#'   \item{\code{p}}{p-values for each parameter estimate}
#'   \item{\code{df}}{Degrees of freedom}
#'   \item{\code{fitted}}{Fitted values}
#'   \item{\code{resid}}{Pearson residuals}
#'   \item{\code{converged}}{Logical value whether the procedure converged}
#'   \item{\code{plot}}{Logical. If TRUE:}
#'
#'     \itemize{
#'       \item\code{ac.glm} Autocorrelation of glm.residuals
#'
#'       \item\code{ac} Autocorrelation of wavelet.residuals
#'   }
#'
#' }
#'
#' @seealso \pkg{waveslim},\code{\link[waveslim]{mra.2d}}
#' @examples
#' data(carlinadata)
#' coords<- carlinadata[,4:5]
#'
#' # scaleWMRR at scale=0 is equivalent to GLM
#' ms0<-scaleWMRR(carlina.horrida ~ aridity + land.use,"poisson",
#'                carlinadata,coords,scale=0,plot=TRUE)

#' # scale-specific regressions for detail components
#' ms1<-scaleWMRR(carlina.horrida ~ aridity + land.use,"poisson",
#'                carlinadata,coords,scale=1,plot=TRUE)
#' ms2<-scaleWMRR(carlina.horrida ~ aridity + land.use,"poisson",
#'                carlinadata,coords,scale=2,plot=TRUE)
#' ms3<-scaleWMRR(carlina.horrida ~ aridity + land.use,"poisson",
#'                carlinadata,coords,scale=3,plot=TRUE)
#'
#' @references
#' Carl G, Doktor D, Schweiger O, Kuhn I (2016)
#' Assessing relative variable importance across different spatial
#' scales: a two-dimensional wavelet analysis.
#' Journal of Biogeography 43: 2502?2512.
#'
#' Whitcher, B. (2005) Waveslim: basic wavelet routines for one-, two-
#' and three-dimensional signal processing. R package version 1.5.
#'
#'@export



scaleWMRR<-function(formula,family,data,coord,
                    scale=1,detail=TRUE,wavelet="haar",wtrafo="dwt",
                    b.ini=NULL,pad=list(),control=list(),moran.params=list(),
                    plot=FALSE){

  n<-dim(data)[1]
  l<-dim(data)[2]
  x<-coord[,1]
  y<-coord[,2]
  if(length(x)!=n) stop("error in dimension")
  logic1<-identical(as.numeric(x),round(x,0))
  logic2<-identical(as.numeric(y),round(y,0))
  if(!logic1 | !logic2) stop("coordinates not integer")
  converged.logi<-TRUE

  X<-model.matrix(formula,data)
  nvar<-dim(X)[2]

  if(is.vector(model.frame(formula,data)[[1]])){
    yold<-model.frame(formula,data)[[1]]
    ntr<-1
  }
  if(family=="binomial" & is.matrix(model.frame(formula,data)[[1]])){
    yold<-model.frame(formula,data)[[1]][,1]
    ntr<-model.frame(formula,data)[[1]][,1] +
      model.frame(formula,data)[[1]][,2]
  }

  n.level<-scale
  length.s<-3*n.level+1
  s<-rep(0,length.s)
  s[(length.s-3):(length.s-1)]<-1
  if(!detail) s[length.s]<-1
  if(scale==0) {
    s<-c(1,1,1,1)
    n.level<-1
  }
  beta<-matrix(NA,4,nvar)
  resi<-matrix(NA,4,n)
  ac<-matrix(NA,4,10)
  se<-matrix(NA,4,nvar)

  pad<-do.call("wrm.pad",pad)
  padform<-pad$padform
  if(scale==0) padform<-0
  if(family!="gaussian" & dimnames(X)[[2]][1]=="(Intercept)") padform<-0
  padzone<-pad$padzone
  control<-do.call("wrm.control",control)
  moran<-do.call("wrm.moran",moran.params)
  lim1<-moran$lim1
  lim2<-lim1 + moran$increment

  pdim<- max(max(y)-min(y)+1,max(x)-min(x)+1)*padzone
  power<-0
  while(2^power<pdim) power<-power+1
  xmargin0<-as.integer((2^power-(max(x)-min(x)+1))/2)-min(x)+1
  ymargin0<-as.integer((2^power-(max(y)-min(y)+1))/2)-min(y)+1

  if(power<n.level) stop("scale is too high")
  if(log2(max(max(y)-min(y)+1,max(x)-min(x)+1))==power & padzone==1)
    print("warning: insufficient padding zone")

  i4<-1
  while(i4<5) {
    if(i4==1){
      xmargin<-xmargin0
      ymargin<-ymargin0
    }
    if(i4==2){
      xmargin<-xmargin0+1
      ymargin<-ymargin0
    }
    if(i4==3){
      xmargin<-xmargin0+1
      ymargin<-ymargin0+1
    }
    if(i4==4){
      xmargin<-xmargin0
      ymargin<-ymargin0+1
    }

    # GLM for comparison
    m0<-glm(formula,family,data)
    res0<-resid(m0,type="pearson")
    beta0<-m0$coeff

    if(is.null(b.ini)){
      betaw<-rep(0,nvar)
    }else{
      betaw<-b.ini
    }

    lin<-X%*%betaw
    if(family=="gaussian") pi<-lin
    if(family=="binomial") pi<-exp(lin)/(1+exp(lin))
    if(family=="poisson") pi<-exp(lin)

    #if(family=="gaussian") pi<-rep(0,n)
    #if(family=="binomial") pi<-rep(.5,n)
    #if(family=="poisson") pi<-rep(1,n)

    it<-0

    repeat{
      it<-it+1
      if(family=="gaussian") var<-rep(1,n)
      if(family=="binomial") var<-ntr*pi*(1-pi)
      if(family=="poisson") var<-pi
      sigma<-as.vector(sqrt(var))
      W12<-diag(sigma)
      Am12<-diag(1/sigma)
      Xnew<-W12%*%X
      ynew<-W12%*%X%*%betaw+Am12%*%(yold-ntr*pi)

      FMat<-matrix(0,2^power,2^power)
      TArray<-array(0,c(2^power,2^power,nvar))
      for(ii in 1:n){
        kx<-x[ii]+xmargin
        ky<-y[ii]+ymargin
        FMat[ky,kx]<-ynew[ii]
        for (i3 in 1:nvar){
          TArray[ky,kx,i3]<-Xnew[ii,i3]
        }
      } # ii loop

      P<-which(is.na(TArray), arr.ind = TRUE)
      if(padform==0){
        FMat[is.na(FMat)]<-0
        for (i3 in 1:nvar){
          i1<-P[which(P[,3]==i3),1]
          i2<-P[which(P[,3]==i3),2]
          for(i in 1:length(i1)){
            TArray[i1[i],i2[i],i3]<-0
          }
        }
      }

      if(padform==1){
        FMat[is.na(FMat)]<-mean(FMat, na.rm=TRUE)
        for (i3 in 1:nvar){
          i1<-P[which(P[,3]==i3),1]
          i2<-P[which(P[,3]==i3),2]
          for(i in 1:length(i1)){
            TArray[i1[i],i2[i],i3]<-mean(TArray[,,i3], na.rm=TRUE)
          }
        }
      }

      if(padform==2){
        FMat<-padding(FMat)
        for (i3 in 1:nvar) TArray[,,i3]<-padding(TArray[,,i3])
      }

      p<-2^power*2^power
      tt<-matrix(0,p,nvar)
      if(is.na(max(abs(FMat)))){
        mdwt$coeff<-rep(NA,nvar)
        break
      }
      if(is.infinite(max(abs(FMat)))){
        mdwt$coeff<-rep(NA,nvar)
        break
      }

      FT<-waveslim::mra.2d(FMat,wavelet,n.level,method=wtrafo)
      FTS<-rep(0,length(FT[[1]]))
      for(is in 1:length(s)){
        if(s[is]==1) FTS <- FTS + FT[[is]]
      }

      ft<-as.vector(FTS)
      for (i3 in 1:nvar){
        TT<-waveslim::mra.2d(TArray[,,i3],wavelet,n.level,method=wtrafo)
        TTS<-rep(0,length(TT[[1]]))
        for(is in 1:length(s)){
          if(s[is]==1){
            TTS <- TTS + TT[[is]]
          }
        }
        tt[,i3]<-as.vector(TTS)
      }

      xnam<-paste("tt[,",1:nvar,"]",sep="")
      formula.dwt<-as.formula(paste("ft~",paste(xnam,collapse="+"),"-1"))
      mdwt<-lm(formula.dwt)
      if(sum(abs(tt[,1]))==0) mdwt$coeff[1]<-beta0[1]
      if(max(abs(mdwt$coeff),na.rm=TRUE)>1e+10 ) {
        mdwt$coeff<-rep(NA,nvar);break}

      lin<-X%*%mdwt$coeff
      if(family=="gaussian") pi<-lin
      if(family=="binomial") pi<-exp(lin)/(1+exp(lin))
      if(family=="poisson") {pi<-exp(lin)
        if(min(pi)<1e-10 |max(pi)>1e+10){
          mdwt$coeff<-rep(NA,nvar)
          break
        }
      }

      if (max(abs(mdwt$coeff-betaw)/
             (abs(betaw)+control$denom.eps)) <= control$eps ){
        i4<-4
        break
      }

      if (i4==4 & it > control$itmax){
        converged.logi<-FALSE
        break
      }
      if (it > control$itmax) break

      betaw<-mdwt$coeff
    } # next step of iteration

    if(sum(abs(tt[,1]))!=0 & !is.na(mdwt$coeff[1])){
      Resmdwt<-matrix(resid(mdwt),2^power,2^power)
      resmdwt<-rep(0,n)
      for(i in 1:n) resmdwt[i]<-Resmdwt[y[i]+ymargin,x[i]+xmargin]
      if(plot) {
        acw<-acfft(x,y,resmdwt,lim1,lim2)
      }
      if(!plot){
        acw<-NA
        acpw<-NA
      }
      if(scale==0) df<-n-nvar
      if(scale!=0) df<-round(n/4^(n.level-1)) -nvar

      if(family=="binomial" | family=="poisson"){
        if(scale==0) var.b<-ginv(t(tt)%*%tt)
        if(scale!=0) var.b<-ginv(t(tt)%*%tt) * (4^(n.level-1))
      }

      if(family=="gaussian"){
        var.b<-ginv(t(tt)%*%tt)
        sigma2<-sum(resmdwt^2)/df
        var.b<-sigma2*var.b
      }
      s.e.<-rep(NA,nvar)
      for(i in 1:nvar){
        s.e.[i]<-sqrt(var.b[i,i])
      }
    }

    if(sum(abs(tt[,1]))==0 | is.na(mdwt$coeff[1])) {
      acw<-NA; resmdwt<-NA; s.e.<-NA
    }

    beta[i4,1:nvar]<-mdwt$coeff[1:nvar]
    resi[i4,1:n ]<-resmdwt[1:n]
    ac[i4,1:10]<-acw[1:10]
    se[i4,1:nvar]<-s.e.[1:nvar]

    i4<-i4+1
  } # i4 loop #.........................................................

  glm.beta<-beta0
  wavelet.beta<-apply(beta,2,mean,na.rm=TRUE)
  if(plot) ac0<-acfft(x,y,res0,lim1,lim2)
  if(!plot) ac0<-NA
  acw<-apply(ac,2,mean,na.rm=TRUE)
  resw<-apply(resi,2,mean,na.rm=TRUE)
  s.e.<-apply(se,2,mean,na.rm=TRUE)
  lin<-X%*%wavelet.beta
  if(family=="gaussian") pi<-lin
  if(family=="binomial") pi<-exp(lin)/(1+exp(lin))
  if(family=="poisson") pi<-exp(lin)

  # test statistics

  z.value<-rep(NA,nvar)
  pr<-rep(NA,nvar)
  if(sum(abs(tt[,1]))!=0 & !is.na(wavelet.beta[1])) {
    z.value<-wavelet.beta/s.e.
    for(i in 1:nvar){
      if(family=="gaussian"){
        if(z.value[i]<=0) pr[i]<-2*pt(z.value[i],df)
        if(z.value[i]>0) pr[i]<-2*(1-pt(z.value[i],df))
      }
      if(family=="binomial" | family=="poisson"){
        if(z.value[i]<=0) pr[i]<-2*pnorm(z.value[i])
        if(z.value[i]>0) pr[i]<-2*(1-pnorm(z.value[i]))
      }
    }
  }

  coef<-as.vector(wavelet.beta)
  names(coef)<-dimnames(X)[[2]]
  call<-match.call()
  fit<-list(call=call,b=coef,s.e.=s.e.,z=z.value,p=pr,df=df,
            fitted=pi,resid=resw,converged=converged.logi,ac.glm=ac0,ac=acw)

  if(plot){
    cat("\n","Call:","\n")
    print(call)
    mra.comp<- matrix(s,1,length(s))
    FTnames<-names(FT)
    colnames(mra.comp)<-c(FTnames[-(length(FTnames))],
                          paste("LL",n.level,sep=""))
    rownames(mra.comp)<-"included"
    cat("---","\n","Selection of multiresolution components: ","\n")
    print(mra.comp)
    cat("---","\n","Pearson Residuals:","\n")
    print(summary(resw))
    beta<-cbind(glm.beta,wavelet.beta,s.e.,z.value,pr)
    beta<-beta[,2:5]
    if(family=="gaussian")
      colnames(beta) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
    if(family=="binomial" | family=="poisson")
      colnames(beta) <- c("Estimate", "Std.Err", "z value", "Pr(>|z|)")
    cat("---","\n","Coefficients:","\n")
    printCoefmat(beta)
  }

  fit
}
