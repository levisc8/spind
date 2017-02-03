#' Wavelet-revised models
#'
#' Description
#' A wavelet-based method to remove spatial autocorrelation
#' in multiple linear regressions
#' Details
#' WRM can be used to fit linear models for response vectors of different
#' distributions: "gaussian", "binomial"(binary) or "poisson". As spatial model,
#' it is a generalized linear model in which the residuals may be autocorrelated.
#' It accounts for spatial (2-dimensional) residual autocorrelation
#' in case of regular gridded datasets. The grid cells are assumed to be square.
#'
#' Arguments:
#' formula   with specified notation according to names in data frame
#' family    "gaussian", "binomial"(binary) or "poisson"
#' data      a data frame
#' coord     a matrix of two columns with corresponding cartesian coordinates,
#'           which have to be integer (grid cell numbering)
#' level     an integer specifying the depth of wavelet decomposition
#'           0 without autocorrelation removal ( = GLM)
#'           1 for best autocorrelation removal
#'           higher integers possible (limit depends on sample size)
#' wavelet   name of wavelet family: "haar" or "d4" or "la8"
#' wtrafo    type of wavelet transform: either "dwt" or "modwt"
#' b.ini     initial values for the parameters b
#' pad       a list of parameters for padding wavelet coefficients.
#'           This is passed to wrm.pad.
#' #########Args for wrm.pad( padform   0 for padding with zeros,
#'                             1 for padding with mean values,
#                              2 for padding with mirror values.
#                              Padform is automatically set to
#                              zero in case of either level=0 or
#                              a formula including an intercept
#                              and a non-gaussian family)
#                             padzone   factor for expanding the padding zone
#
#' control 	a list of parameters for controlling the fitting process.
#'           This is passed to wrm.control.
#' moran     a list of parameters for calculating Moran's I.
#'           This is passed to wrm.moran.
#' plot      a logical value indicating whether results should be plotted
#' graph     a logical value indicating whether results should be displayed
#'
#' Value:
#' An object of class "wrm" is a list containing the following components:
#'          call       call
#'          formula    formula
#'          family     family
#'          b          estimate of regression parameters
#'          s.e.       standard errors
#'          z          z values (or corresponding values for statistic)
#'          p          probabilities
#'          fitted     fitted values
#'          resid      Pearson residuals
#'          b.sm       estimate of regression parameters of neglected smooth part
#'          fitted.sm  fitted values of neglected smooth part
#'          level      level
#'          wavelet    wavelet
#'          wtrafo     wtrafo
#'          padzone    padzone
#'          padform    padform
#'          n.eff      effective number of observations
#'          AIC        Akaike information criterion
#'          if plot or graph is true:
#'          ac.glm     autocorrelation of glm.residuals
#'          ac         autocorrelation of wavelet.residuals
#' References
#' Carl, G., Kuhn, I. (2010): A wavelet-based extension of generalized
#' linear models to remove the effect of spatial autocorrelation.
#' Geogr. Anal. 42 (3), 323 - 337


WRM<-function(formula,family,data,coord,
              level=1,wavelet="haar",wtrafo="dwt",
              b.ini=NULL,pad=list(),control=list(),moran=list(),
              plot=FALSE,graph=FALSE){

  n<-dim(data)[1]
  l<-dim(data)[2]
  x<-coord[,1]
  y<-coord[,2]
  if(length(x)!=n) stop("error in dimension")
  logic1<-identical(as.numeric(x),round(x,0))
  logic2<-identical(as.numeric(y),round(y,0))
  if(!logic1 | !logic2) stop("coordinates not integer")

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

  n.level<-level
  length.s<-3*n.level+1
  s<-rep(1,length.s)
  s[length.s]<-0
  if(level==0) {s<-c(1,1,1,1) ; n.level<-1}
  beta<-matrix(NA,4,nvar)
  beta.smooth<-matrix(NA,4,nvar)
  resi<-matrix(NA,4,n)
  ac<-matrix(NA,4,10)
  acp<-matrix(NA,4,10)
  se<-matrix(NA,4,nvar)

  pad<-do.call("wrm.pad",pad)
  padform<-pad$padform
  if(level==0) padform<-0
  if(family!="gaussian" & dimnames(X)[[2]][1]=="(Intercept)") padform<-0
  padzone<-pad$padzone
  control<-do.call("wrm.control",control)
  moran<-do.call("wrm.moran",moran)
  lim1<-moran$lim1
  lim2<-lim1 + moran$increment

  pdim<- max(max(y)-min(y)+1,max(x)-min(x)+1)*padzone
  power<-0
  while(2^power<pdim) power<-power+1
  xmargin0<-as.integer((2^power-(max(x)-min(x)+1))/2)-min(x)+1
  ymargin0<-as.integer((2^power-(max(y)-min(y)+1))/2)-min(y)+1

  if(power<n.level) stop("level is too high")
  if(log2(max(max(y)-min(y)+1,max(x)-min(x)+1))==power & padzone==1)
    message("warning: insufficient padding zone")

  i4<-1
  while(i4<5){
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


    if(is.null(b.ini)) {
      betaw<-rep(0,nvar)
    } else {
      betaw<-b.ini
    }
    lin<-X%*%betaw
    if(family=="gaussian") pi<-lin
    if(family=="binomial") pi<-exp(lin)/(1+exp(lin))
    if(family=="poisson")  pi<-exp(lin)

    #if(family=="gaussian") pi<-rep(0,n)
    #if(family=="binomial") pi<-rep(.5,n)
    #if(family=="poisson") pi<-rep(1,n)

    it<-0

    repeat{
      it<-it+1
      if(family=="gaussian") var<-rep(1,n)
      if(family=="binomial") var<-ntr*pi*(1-pi)
      if(family=="poisson")  var<-pi
      sigma<-as.vector(sqrt(var))
      W12<-diag(sigma)
      Am12<-diag(1/sigma)
      Xnew<-W12%*%X
      ynew<-W12%*%X%*%betaw+Am12%*%(yold-ntr*pi)

      F.mat<-matrix(NA,2^power,2^power)
      T.array<-array(NA,c(2^power,2^power,nvar))
      for(ii in 1:n){
        kx<-x[ii]+xmargin
        ky<-y[ii]+ymargin
        F.mat[ky,kx]<-ynew[ii]
        for (i3 in 1:nvar)
          T.array[ky,kx,i3]<-Xnew[ii,i3]
      }  # ii loop

      P<-which(is.na(T.array), arr.ind = TRUE)
      if(padform==0){
        F.mat[is.na(F.mat)]<-0
        for (i3 in 1:nvar){
          i1<-P[which(P[,3]==i3),1]
          i2<-P[which(P[,3]==i3),2]
          for(i in 1:length(i1)) T.array[i1[i],i2[i],i3]<-0
        }
      }
      if(padform==1){
        F.mat[is.na(F.mat)]<-mean(F.mat, na.rm=TRUE)
        for (i3 in 1:nvar){
          i1<-P[which(P[,3]==i3),1]
          i2<-P[which(P[,3]==i3),2]
          for(i in 1:length(i1)) T.array[i1[i],i2[i],i3]<-mean(T.array[,,i3], na.rm=TRUE)
        }
      }
      if(padform==2){
        F.mat<-padding(F.mat)
        for (i3 in 1:nvar) T.array[,,i3]<-padding(T.array[,,i3])
      }

      p<-2^power*2^power
      tt<-matrix(0,p,nvar)
      tt0<-matrix(0,p,nvar)
      if(is.na(max(abs(F.mat)))) {
        mdwt$coeff<-rep(NA,nvar)
        break
      }
      if(is.infinite(max(abs(F.mat)))) {
        mdwt$coeff<-rep(NA,nvar)
        break
      }

      FT<-waveslim::mra.2d(F.mat,wavelet,n.level,method=wtrafo)
      FTS<-rep(0,length(FT[[1]]))
      for(is in 1:length(s)){
        if(s[is]==1) FTS <- FTS + FT[[is]]
        if(s[is]==0) FT0 <- FT[[is]]
      }
      ft<-as.vector(FTS)
      if(level!=0) ft0<-as.vector(FT0)
      for (i3 in 1:nvar){
        TT<-waveslim::mra.2d(T.array[,,i3],wavelet,n.level,method=wtrafo)
        TTS<-rep(0,length(TT[[1]]))
        TT0<-rep(0,length(TT[[1]]))
        for(is in 1:length(s)){
          if(s[is]==1) TTS <- TTS + TT[[is]]
          if(s[is]==0) TT0 <- TT[[is]]
        }
        tt[,i3]<-as.vector(TTS)
        if(level!=0) tt0[,i3]<-as.vector(TT0)
      }

      xnam<-paste("tt[,",1:nvar,"]",sep="")
      formula.dwt<-as.formula(paste("ft~",paste(xnam,collapse="+"),"-1"))
      mdwt<-lm(formula.dwt)
      if(sum(abs(tt[,1]))==0) mdwt$coeff[1]<- beta0[1]
      if(is.na(max(abs(mdwt$coeff))) ) {mdwt$coeff<-rep(NA,nvar);break}
      if(max(abs(mdwt$coeff))>1e+10 ) {mdwt$coeff<-rep(NA,nvar);break}
      if(level!=0){
        xnam0<-paste("tt0[,",1:nvar,"]",sep="")
        formula.dwt0<-as.formula(paste("ft0~",paste(xnam0,collapse="+"),"-1"))
        mdwt0<-lm(formula.dwt0)
        if(sum(abs(tt0[,1]))==0) mdwt0$coeff[1]<- beta0[1]
      }

      lin<-X%*%mdwt$coeff
      if(family=="gaussian") pi<-lin
      if(family=="binomial") pi<-exp(lin)/(1+exp(lin))
      if(family=="poisson") {
        pi<-exp(lin)
        if(min(pi)<1e-10 |max(pi)>1e+10) {
          mdwt$coeff<-rep(NA,nvar)
          break
        }
      }

      if (max(abs(mdwt$coeff-betaw)/(abs(betaw)+control$denom.eps))
          <= control$eps ) {i4<-4 ; break}
      if (i4==4 & it > control$itmax)  stop("too many iterations")
      if (it > control$itmax) break
      betaw<-mdwt$coeff
    }  # next step of iteration

    if(level!=0) Lin<-tt%*%mdwt$coeff + tt0%*%mdwt0$coeff
    if(level==0) Lin<-tt%*%mdwt$coeff
    Lin.Mat<-matrix(Lin,2^power,2^power)
    lin<-rep(0,n)
    for(i in 1:n) lin[i]<-Lin.Mat[y[i]+ymargin,x[i]+xmargin]
    if(family=="gaussian") fitted.sm<-lin
    if(family=="binomial") fitted.sm<-exp(lin)/(1+exp(lin))
    if(family=="poisson")  fitted.sm<-exp(lin)

    if(!is.na(mdwt$coeff[1])){
      Resmdwt<-matrix(resid(mdwt),2^power,2^power)
      resmdwt<-rep(0,n)
      for(i in 1:n) resmdwt[i]<-Resmdwt[y[i]+ymargin,x[i]+xmargin]
      if(plot | graph) {
        acw<-acfft(x,y,resmdwt,lim1,lim2)
      }

      if(!plot & !graph) {
        acw<-NA
        acpw<-NA
      }

      try.solve<-try(ginv(t(tt)%*%tt),silent = TRUE)
      if (inherits(try.solve, "try-error"))  var.b<-matrix(NA,nvar,nvar)
      if (!inherits(try.solve, "try-error")) var.b<-try.solve
      if(family=="gaussian"){
        df<-n-nvar
        sigma2<-sum(resmdwt^2)/df
        var.b<-sigma2*var.b
      }
      s.e.<-rep(NA,nvar)

      for(i in 1:nvar){
        s.e.[i]<-sqrt(var.b[i,i])
      }
    }

    if(is.na(mdwt$coeff[1])) {acw<-NA; acpw<-NA; resmdwt<-NA; s.e.<-NA}
    beta[i4,1:nvar]<-mdwt$coeff[1:nvar]
    if(level!=0) beta.smooth[i4,1:nvar]<-mdwt0$coeff[1:nvar]
    resi[i4,1:n ]<-resmdwt[1:n]
    ac[i4,1:10]<-acw[1:10]
    se[i4,1:nvar]<-s.e.[1:nvar]

    i4<-i4+1
  } # i4 loop #..................................................................

  glm.beta<-beta0
  wavelet.beta<-apply(beta,2,mean,na.rm=TRUE)
  if(level!=0) wavelet.beta.smooth<-apply(beta.smooth,2,mean,na.rm=TRUE)
  if(plot | graph) ac0<-acfft(x,y,res0,lim1,lim2)
  if(!plot & !graph) ac0<-NA
  acw<-apply(ac,2,mean,na.rm=TRUE)
  resw<-apply(resi,2,mean,na.rm=TRUE)
  s.e.<-apply(se,2,mean,na.rm=TRUE)
  lin<-X%*%wavelet.beta
  if(family=="gaussian") pi<-lin
  if(family=="binomial") pi<-exp(lin)/(1+exp(lin))
  if(family=="poisson")  pi<-exp(lin)
  if(level==0) n.eff<-n
  if(level!=0) n.eff<-round(n*(1- 1/(2^level*2^level)  ))
  aic<-aic.calc(formula,family,data,mu=pi,n.eff=n.eff)
  AIC<-round(aic$AIC,1)


  z.value<-rep(NA,nvar)
  pr<-rep(NA,nvar)
  if(!is.na(wavelet.beta[1]) & !is.na(s.e.[1]) ) {
    z.value<-wavelet.beta/s.e.
    for(i in 1:nvar){
      if(family=="gaussian"){
        if(z.value[i]<=0) pr[i]<-2*pt(z.value[i],df)
        if(z.value[i]>0)  pr[i]<-2*(1-pt(z.value[i],df))
      }
      if(family=="binomial" | family=="poisson"){
        if(z.value[i]<=0) pr[i]<-2*pnorm(z.value[i])
        if(z.value[i]>0)  pr[i]<-2*(1-pnorm(z.value[i]))
      }
    }
  }

  if(plot){
    beta<-cbind(glm.beta,wavelet.beta,s.e.,z.value,pr)
    beta<-beta[,2:5]
    if(family=="gaussian")
      colnames(beta) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
    if(family=="binomial" | family=="poisson")
      colnames(beta) <- c("Estimate", "Std.Err", "z value", "Pr(>|z|)")
    cat("---","\n","Coefficients:","\n")
    printCoefmat(beta)
    cat("---","\n","n.eff: ",n.eff,",  AIC: ",AIC,"\n" )
  }

  if(plot){
    cat("---","\n")
    cat("Autocorrelation of glm.residuals","\n")
    print(ac0)
    cat("Autocorrelation of wavelet.residuals","\n")
    print(acw)
  }

  if(graph & !is.na(acw[1])){
    y1<-min(min(ac0),min(acw))-.1
    y2<-max(max(ac0),max(acw))+.1
    plot(ac0,type="b",ylim=c(y1,y2),
         ylab="Autocorrelation of residuals", xlab="Lag distance",
         main=paste("Autocorrelation for level = ", level))
    points(acw,pch=2,type="b")
    v<-1:2
    leg<-c("glm.residuals","wavelet.residuals")
    legend(6,y2-.1,leg,pch=v)
  }

  glm.beta<-beta0
  coef<-as.vector(wavelet.beta)
  names(coef)<-dimnames(X)[[2]]
  if(level!=0) coefsm<-as.vector(wavelet.beta.smooth)
  if(level==0) coefsm<-rep(NA,length(coef))
  names(coefsm)<-dimnames(X)[[2]]

  call<-match.call()
  fit<-list(call=call,
            formula=formula,
            family=family,
            b=coef,
            s.e.=s.e.,
            z=z.value,
            p=pr,
            fitted=pi,
            resid=resw,
            b.sm=coefsm,
            fitted.sm=fitted.sm,
            level=level,
            wavelet=wavelet,
            wtrafo=wtrafo,
            padzone=padzone,
            padform=padform,
            it=it,
            n=n,
            n.eff=n.eff,
            AIC=AIC,
            ac.glm=ac0,
            ac=acw)

  class(fit)<-"WRM"
  return(fit)
}
