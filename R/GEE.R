GEE<-function(formula,family,data,coord,
              corstr="fixed",cluster=3,moran.params=list(),
              plot=FALSE,graph=FALSE,silent=TRUE){
  ###############################################################################
  # GEE (Generalized Estimating Equations)
  #
  # Description
  # A GEE-based method to account for spatial autocorrelation
  # in multiple linear regressions
  # Details
  # GEE can be used to fit linear models for response vectors of different
  # distributions: "gaussian", "binomial"(binary) or "poisson". As spatial model,
  # it is a generalized linear model in which the residuals may be autocorrelated.
  # It accounts for spatial (2-dimensional) residual autocorrelation
  # in case of regular gridded datasets. The grid cells are assumed to be square.
  #
  # Arguments:
  # formula  with specified notation according to names in data frame
  # family   "gaussian", "binomial"(binary) or "poisson"
  # data     data frame
  # coord    a matrix of two columns with corresponding cartesian coordinates,
  #          which have to be integer (grid cell numbering)
  # corstr   autocorrelation structure: "independence", "fixed",
  #          "exchangeable", "quadratic"  are possible.
  #          "independence" = GLM, i.e. correlation matrix = identity matrix;
  #          "fixed" for best autocorrelation removal by means of an adapted
  #          isotropic power function specifying all correlation coefficients;
  #          "exchangeable", "quadratic" for clustering, i.e. the correlation
  #          matrix has a block diagonal form:
  #          "exchangeable": all intra-block correlation coefficients are equal;
  #          "quadratic": intra-block correlation coefficients for different
  #          distances can be different.
  # cluster  cluster size for cluster models: "exchangeable", "quadratic".
  #          cluster= 2,3,4 are allowed for size 2*2,3*3,4*4, respectively
  # moran    a list of parameters for calculating Moran's I.
  #          This is passed to wrm.moran.
  # plot     a logical value indicating whether results should be plotted
  # graph    a logical value indicating whether results should be displayed
  #
  # Value:
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
  # References
  # Carl & Kühn (2007): Analyzing Spatial Autocorrelation in Species
  # Distributions using Gaussian and Logit Models, Ecol. Model. 207, 159 - 170
  ###############################################################################

  at<-intersect(names(data),all.vars(formula))
  if(length(at)==0) stop("formula: specified notation is missing")
  nn<-nrow(data)
  x<-coord[,1]
  y<-coord[,2]
  if(length(x)!=nn) stop("error in dimension")
  logic1<-identical(as.numeric(x),round(x,0))
  logic2<-identical(as.numeric(y),round(y,0))
  if(!logic1 | !logic2) stop("coordinates not integer")

  moran<-do.call("wrm.moran",moran.params)
  lim1<-moran$lim1
  lim2<-lim1 + moran$increment

  m0<-glm(formula,family,data)
  res0<-resid(m0,type="pearson")
  id<-rep(1,nn)
  dato<-data.frame(data,id)
  mgee<-gee(formula=formula,family=family,data=dato,id=id,
            corstr="independence",silent=silent)
  var.indep.naive<-mgee$naive.variance

  if(corstr=="independence"){
    ashort<-0
    A<-0
    fitted<-fitted(m0)
    resid<-resid(m0,type="pearson")
    b<-summary(m0)$coefficients[,1]
    s.e.<-summary(m0)$coefficients[,2]
    z<-summary(m0)$coefficients[,3]
    p<-summary(m0)$coefficients[,4]
    scale<-summary(m0)$dispersion
    Icrit<-qic.culc(formula,data,family,fitted,var.indep.naive,var.indep.naive)
    QIC<-Icrit$QIC
  }

  if(corstr=="fixed"){
    ac01<-acfft(x,y,res0,1,1.1,dmax=1)
    ac05<-acfft(x,y,res0,5,5.1,dmax=1)
    if(ac05<=0) v<-1
    if(ac05>0) v<-log(log(ac05)/log(ac01))/log(5)
    alpha<-ac01
    para0<-paste("n",nn,sep="=")
    para1<-paste(", alpha",round(alpha,3),sep="=")
    para2<-paste(", v",round(v,3),sep="=")
    A0<-paste(para0,para1,para2)
    id<-rep(1,nn)
    coord<-cbind(x,y)
    D<-as.matrix(dist(coord))
    R<-alpha^(D^v)
    data<-data.frame(data,id)
    mgee<-gee(formula=formula,family=family,data=data,id=id,R=R,corstr="fixed")
    var.naive<-mgee$naive.variance
    para3<-"a=alpha^(d^v) "
    ashort<-c(alpha,v)
    A<-paste(para3,para1,para2)
    b<-mgee$coeff
    res<-res.gee(formula,family,data,nn,b=mgee$coeff,R=R)
    fitted<-res$fitted
    resid<-res$resid
    s.e.=summary(mgee)$coefficients[,2]
    z<-summary(mgee)$coefficients[,3]
    p<-rep(NA,nrow(summary(mgee)$coefficients))
    for(ii in 1:nrow(summary(mgee)$coefficients)){
      if(z[ii]>=0) p[ii]<-2*(1-pnorm(z[ii]))
      if(z[ii]<0) p[ii]<-2*(pnorm(z[ii]))
    }
    scale<-summary(mgee)[[9]]
    Icrit<-qic.culc(formula,data,family,fitted,var.naive,var.indep.naive)
    QIC<-Icrit$QIC
   }


  if(corstr=="exchangeable") {
    dato<-dat.nn(data,coord,cluster)
    l<-dim(dato)[2]
    o<-dato[,l-2]
    id<-dato[,l-1]
    waves<-dato[,l]

    clusz<-clus.sz(id)
    zcor<-genZcor(clusz=clusz,waves=waves,"unstr")
    mgee<-gee(formula=formula,family=family,data=dato,id=id,corstr="exchangeable")
    var.robust<-mgee$robust.variance
    ashort<-mgee$w[1,2]
    a<-a.gee(mgee$w,cluster,type="gee",corstr="exchangeable")
    A<-mgee$w
    b<-mgee$coeff
    res<-res.gee(formula,family,dato,cluster,clusz,zcor,a,b)
    fitted<-res$fitted[order(o)]
    resid<-res$resid[order(o)]
    s.e.=summary(mgee)$coefficients[,4]
    z<-summary(mgee)$coefficients[,5]
    p<-rep(NA,nrow(summary(mgee)$coefficients))
    for(ii in 1:nrow(summary(mgee)$coefficients)){
      if(z[ii]>=0) p[ii]<-2*(1-pnorm(z[ii]))
      if(z[ii]<0) p[ii]<-2*(pnorm(z[ii]))
    }

    scale<-summary(mgee)[[9]]
    Icrit<-qic.culc(formula,data,family,fitted,var.robust,var.indep.naive)
    QIC<-Icrit$QIC

   }


  if(corstr=="quadratic"){
    dato<-dat.nn(data,coord,cluster)
    l<-dim(dato)[2]
    o<-dato[,l-2]
    id<-dato[,l-1]
    waves<-dato[,l]

    clusz<-clus.sz(id)
    zcor<-genZcor(clusz=clusz,waves=waves,"unstr")
    zcorq<-zcor.quad(zcor,cluster,quad=T)
    mgeese<-geese(formula=formula,family=family,data=dato,id=id,corstr=
                    "userdefined",zcor=zcorq)
    var.robust<-mgeese$vbeta
    ashort<-mgeese$a
    a<-a.gee(mgeese$a,cluster,type="geese",corstr="userdefined",quad=T)
    A<-cor.mat(cluster,a)
    b<-mgeese$b
    res<-res.gee(formula,family,dato,cluster,clusz,zcor,a,b)
    fitted<-res$fitted[order(o)]
    resid<-res$resid[order(o)]
    s.e.=summary(mgeese)$mean[,2]
    z<-summary(mgeese)$mean[,3]
    p<-summary(mgeese)$mean[,4]

    scale<-as.numeric(summary(mgeese)$scale[1])
    Icrit<-qic.culc(formula,data,family,fitted,var.robust,var.indep.naive)
    QIC<-Icrit$QIC
  }

  if(!plot & !graph) {ac0<-NA; ac<-NA}
  if(plot | graph) {
    x<-coord[,1]
    y<-coord[,2]
    ac0<-acfft(x,y,res0,lim1,lim2)
    ac<-acfft(x,y,resid,lim1,lim2)
  }

  if(graph){
    y1<-min(min(ac0),min(ac))-.1
    y2<-max(max(ac0),max(ac))+.1
    plot(ac0,type="b",ylim=c(y1,y2),
         ylab="Autocorrelation of residuals", xlab="Distance between points",
         main=paste("Autocorrelation for correlation structure: ", corstr))
    points(ac,pch=2,type="b")
    v<-1:2
    leg<-c("GLM Residuals","GEE Residuals")
    legend(6,y2-.1,leg,pch=v)
  }

  call<-match.call()
  fit<-list(call=call,formula=formula,family=family,b=b,s.e.=s.e.,z=z,p=p,
            scale=scale,fitted=fitted,resid=resid,w.ac=ashort,W.ac=A,QIC=QIC,
            ac.glm=ac0,ac=ac)
  class(fit)<-"gee"
  return(fit)

}

qic.culc<-function(formula,data,family,mu,var.robust,var.indep.naive){
  ###############################################################################
  # Description
  # A function for calculating quasi-likelihood and Quasi-Information Criterion
  # values.
  # QIC as defined in Hardin & Hilbe (2003).
  # Reference: Barnett et al. Methods in Ecology & Evolution 2010, 1, 1524.
  # Arguments:
  # formula  with specified notation according to names in data frame
  # family   "gaussian", "binomial"(binary) or "poisson"
  # data     data frame
  # mu       fitted values
  # var.robust        variance of b values
  # var.indep.naive   naive variance of b values under the "independence" model
  #
  # value:   QIC        quasi-information criterion
  #          loglik     quasi-likelihood
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
    loglik<- -n/2*(log(2*pi*sigma2) +1)   # log likelihood
  }
  if(family=="binomial"){                         # choose= Binomialkoeff.
    loglik<- sum(y*log(mu/(1-mu)) + ntr*log(1-mu) + log(choose(ntr,y)) )
  }
  if(family=="poisson"){
    # loglik<- sum(y*log(mu)-mu)  # useful for delta in multimodel inference
    loglik<- sum(y*log(mu)-mu) - sum(log(factorial(y))) # factorial=Fakultät
  }

  trace<-sum(diag(ginv(var.indep.naive)%*% var.robust))

  QIC<- -2*loglik + 2*trace

  return(list(QIC=QIC, loglik=loglik))

}


wrm.moran<-function(lim1=0,increment=1){
  ###############################################################################
  # Description
  # Auxiliary function for specifying the Moran's I calculation.
  # Typically only used internally by WRM.
  # Arguments:
  # lim1      lower limit for first bin
  # increment increment (=1 and lim1=0 is conform to correlog{ncf})
  #
  # Value:     A list with components named as the arguments.
  ###############################################################################
  moran.list<-list(lim1=lim1,increment=increment)
  return(moran.list)
}


dat.nn<-function(data,coord,n){
  ###############################################################################
  # Description
  # A function to generate clusters and order variables
  # Arguments
  # data      a data frame
  # coord     a matrix of two columns with corresponding coordinates
  # n         for maximal cluster size  n*n
  #
  # Value: a new data frame containing rearranged data, rearranged coordinates,
  # and 3 new parameters:
  #    o      order parameter
  #    id     parameter identifying clusters
  #    waves  parameter identifying members of clusters
  #
  ###############################################################################

  l<-dim(data)[2]
  OST<-coord[,1]
  NORD<-coord[,2]
  ko<-OST-min(OST)
  idx<-(ko-(ko%%(n)))/n+1
  ks<-NORD-min(NORD)
  idy<-(ks-(ks%%(n)))/n+1
  ie<-(idy-1)*max(idx)+idx
  idwx<-ko%%(n)+1
  idwy<-ks%%(n)+1
  wav<-(idwy-1)*n+idwx
  data<-as.matrix(data)
  o<-order(ie,wav)
  x<-OST[o]
  y<-NORD[o]
  id<-ie[o]
  waves<-wav[o]
  dat.new1<-data[o,]
  dat.new2<-cbind(dat.new1,x,y,o,id,waves)
  dat.new<-as.data.frame(dat.new2)
}


clus.sz<-function(id){
  ###############################################################################
  # Description
  # A function to calculate sizes of clusters
  # Argument
  # id     a vector that identifies clusters
  # Value:
  # A vector of numbers of cluster sizes
  ###############################################################################

  clus<-rep(0,length(id))
  k0<-0
  k1<-1
  for(i in 2:length(id)) { i1<-i-1
  if(id[i]==id[i1]) {k1<-k1+1
  if(i==length(id)) {k0<-k0+1
  clus[k0]<-k1}}
  if(id[i]!=id[i1]) {k0<-k0+1
  clus[k0]<-k1
  k1<-1
  if(i==length(id)) {k0<-k0+1
  clus[k0]<-k1 }}}
  clusz<-clus[clus>0]
}


zcor.quad<-function(zcor,n,quad=TRUE) {
  ###############################################################################
  # Description
  # A function to create a quadratic correlation structure
  # Arguments:
  # zcor    an object of class "genZcor" (see: geepack)
  # n       for maximal cluster size n*n
  # quad    by default quadratic correlation structure
  # Value:
  # A matrix describing the quadratic correlation structure
  ###############################################################################

  if(quad) {
    if(n==2)  {
      zcorn<-matrix(0,dim(zcor)[1],2)
      zcorn[,1]<-zcor[,1]+zcor[,2]+zcor[,5]+zcor[,6]
      zcorn[,2]<-zcor[,3]+zcor[,4]
    }
    if(n==3)  {
      zcorn<-matrix(0,dim(zcor)[1],5)
      zcorn[,1]<-zcor[,1]+zcor[,3]+zcor[,9]+zcor[,11]+zcor[,18]+zcor[,22]+
        zcor[,24]+zcor[,27]+zcor[,29]+zcor[,33]+zcor[,34]+zcor[,36]
      zcorn[,2]<-zcor[,2]+zcor[,6]+zcor[,14]+zcor[,21]+zcor[,23]+zcor[,35]
      zcorn[,3]<-zcor[,4]+zcor[,10]+zcor[,12]+zcor[,17]+zcor[,25]+zcor[,28]+
        zcor[,30]+zcor[,32]
      zcorn[,4]<-zcor[,5]+zcor[,7]+zcor[,13]+zcor[,15]+zcor[,16]+zcor[,20]+
        zcor[,26]+zcor[,31]
      zcorn[,5]<-zcor[,8]+zcor[,19]
    }
    if(n==4)  {
      zcorn<-matrix(0,dim(zcor)[1],9)
      zcorn[,1]<-zcor[,1]+zcor[,4]+zcor[,16]+zcor[,19]+zcor[,30]+zcor[,33]+
        zcor[,46]+zcor[,55]+zcor[,58]+zcor[,66]+zcor[,69]+zcor[,76]+
        zcor[,79]+zcor[,88]+zcor[,93]+zcor[,96]+zcor[,100]+zcor[,103]+
        zcor[,106]+zcor[,109]+zcor[,114]+zcor[,115]+zcor[,118]+zcor[,120]
      zcorn[,2]<-zcor[,2]+zcor[,8]+zcor[,17]+zcor[,23]+zcor[,37]+zcor[,50]+
        zcor[,56]+zcor[,62]+zcor[,67]+zcor[,73]+zcor[,83]+zcor[,92]+
        zcor[,94]+zcor[,101]+zcor[,116]+zcor[,119]
      zcorn[,3]<-zcor[,3]+zcor[,12]+zcor[,27]+zcor[,41]+zcor[,54]+zcor[,57]+
        zcor[,95]+zcor[,117]
      zcorn[,4]<-zcor[,5]+zcor[,18]+zcor[,20]+zcor[,32]+zcor[,34]+zcor[,45]+
        zcor[,59]+zcor[,68]+zcor[,70]+zcor[,78]+zcor[,80]+zcor[,87]+
        zcor[,97]+zcor[,102]+zcor[,104]+zcor[,108]+zcor[,110]+zcor[,113]
      zcorn[,5]<-zcor[,6]+zcor[,9]+zcor[,21]+zcor[,22]+zcor[,24]+zcor[,31]+
        zcor[,36]+zcor[,38]+zcor[,44]+zcor[,49]+zcor[,60]+zcor[,63]+
        zcor[,71]+zcor[,72]+zcor[,74]+zcor[,77]+zcor[,82]+zcor[,84]+
        zcor[,86]+zcor[,91]+zcor[,98]+zcor[,105]+zcor[,107]+zcor[,112]
      zcorn[,6]<-zcor[,7]+zcor[,13]+zcor[,26]+zcor[,28]+zcor[,40]+zcor[,42]+
        zcor[,43]+zcor[,53]+zcor[,61]+zcor[,85]+zcor[,99]+zcor[,111]
      zcorn[,7]<-zcor[,10]+zcor[,25]+zcor[,35]+zcor[,48]+zcor[,64]+zcor[,75]+
        zcor[,81]+zcor[,90]
      zcorn[,8]<-zcor[,11]+zcor[,14]+zcor[,29]+zcor[,39]+zcor[,47]+zcor[,52]+
        zcor[,65]+zcor[,89]
      zcorn[,9]<-zcor[,15]+zcor[,51]
    }
  }
  if(!quad) zcorn<-zcor
  zcorn<-as.matrix(zcorn)
}


a.gee<-function(mgee,n,type="glm",corstr="independence",quad=T) {
  ###############################################################################
  # Description
  # A function to order correlation parameters of Generalized Estimating
  # Equation Models
  # Arguments
  # mgee       matrix or vector of correlation parameters according to model
  # n          for maximal cluster size n*n
  # type       type of model:
  #            "glm", "gee", "geese" are allowed
  # corstr     correlation structure:
  #            "independence", "exchangeable", "userdefined" are allowed
  # quad       by default quadratic correlation structure
  #            for model "geese" and "userdefined" correlation only
  # Value:
  # A vector of correlation parameters
  ###############################################################################

  if(n==2)n3<-6
  if(n==3)n3<-36
  if(n==4)n3<-120
  a<-rep(0,n3)
  if(type=="glm") a<-a
  if(type=="gee"){
    if(corstr=="exchangeable") a[c(1:n3)]<-mgee[1,2]
    if(corstr=="independence") a<-a
  }
  a<-as.vector(a)

  if(type=="geese") {
    if(corstr=="userdefined"){
      if(quad) {
        if(n==2)  {
          a<-rep(0,6)
          a[c(1,2,5,6)]<-mgee[1]
          a[c(3,4)]<-mgee[2]
        }
        if(n==3)  {
          a<-rep(0,36)
          a[c(1,3,9,11,18,22,24,27,29,33,34,36)]<-mgee[1]
          a[c(2,6,14,21,23,35)]<-mgee[2]
          a[c(4,10,12,17,25,28,30,32)]<-mgee[3]
          a[c(5,7,13,15,16,20,26,31)]<-mgee[4]
          a[c(8,19)]<-mgee[5]
        }
        if(n==4)  {
          a<-rep(0,120)
          a[c(1,4,16,19,30,33,46,55,58,66,69,76,79,88,93,96,100,103,106,109,
              114,115,118,120)]<-mgee[1]
          a[c(2,8,17,23,37,50,56,62,67,73,83,92,94,101,116,119)]<-mgee[2]
          a[c(3,12,27,41,54,57,95,117)]<-mgee[3]
          a[c(5,18,20,32,34,45,59,68,70,78,80,87,97,102,104,108,110,113)]<-mgee[4]
          a[c(6,9,21,22,24,31,36,38,44,49,60,63,71,72,74,77,82,84,86,91,98,
              105,107,112)]<-mgee[5]
          a[c(7,13,26,28,40,42,43,53,61,85,99,111)]<-mgee[6]
          a[c(10,25,35,48,64,75,81,90)]<-mgee[7]
          a[c(11,14,29,39,47,52,65,89)]<-mgee[8]
          a[c(15,51)]<-mgee[9]
        }}
      if(!quad) a<-mgee
    }
    if(corstr=="exchangeable") a[c(1:n3)]<-mgee
    if(corstr=="independence") a<-a
  }
  a<-as.vector(a)
}


cor.mat<-function(cluster,a) {
  ###############################################################################
  # Description
  # A function to create a block of the correlation matrix
  # Arguments:
  # cluster     cluster size
  # a           cluster parameter
  # Value:
  # A matrix representing a block of the correlation matrix
  ###############################################################################
  n<-cluster
  n2<-cluster*cluster
  z2<-a
  v1<-matrix(0,n2,n2)
  if(n==2)
  {  v1[1,2:4]<-z2[1:3]
  v1[2,3:4]<-z2[4:5]
  v1[3,4]<-z2[6]  }
  if(n==3)
  {  v1[1,2:9]<-z2[1:8]
  v1[2,3:9]<-z2[9:15]
  v1[3,4:9]<-z2[16:21]
  v1[4,5:9]<-z2[22:26]
  v1[5,6:9]<-z2[27:30]
  v1[6,7:9]<-z2[31:33]
  v1[7,8:9]<-z2[34:35]
  v1[8,9]<-z2[36]  }
  if(n==4)
  {  v1[1,2:16]<-z2[1:15]
  v1[2,3:16]<-z2[16:29]
  v1[3,4:16]<-z2[30:42]
  v1[4,5:16]<-z2[43:54]
  v1[5,6:16]<-z2[55:65]
  v1[6,7:16]<-z2[66:75]
  v1[7,8:16]<-z2[76:84]
  v1[8,9:16]<-z2[85:92]
  v1[9,10:16]<-z2[93:99]
  v1[10,11:16]<-z2[100:105]
  v1[11,12:16]<-z2[106:110]
  v1[12,13:16]<-z2[111:114]
  v1[13,14:16]<-z2[115:117]
  v1[14,15:16]<-z2[118:119]
  v1[15,16]<-z2[120]   }
  for(i in 1:n2) v1[i,i]<-0.5
  v<-v1+t(v1)
  v
}


acfft<-function(x,y,reslm,lim1=1,lim2=2,dmax=10){
  ###############################################################################
  # Description
  # A function for calculating spatial autocorrelation (i.e. Moran's I values).
  # Arguments:
  # x 	  a vector of length n representing the x coordinates
  # y 	  a vector of length n representing the y coordinates
  # reslm a vector of length n representing the residuals at each location.
  # lim1  lower bound for first bin
  # lim2  upper bound for first bin
  # dmax  number of bins, i.e. uniformly distributed distance classes
  #
  # Value: acfft returns a vector of Moran's I values
  ###############################################################################
  reslm<-reslm-mean(reslm)
  mi<-max(x)-min(x)+1
  mk<-max(y)-min(y)+1
  n<-max(mi,mk)
  n2<-n*n
  Ares<-matrix(0,n,n)
  mask<-matrix(0,n,n)
  for(i in 1:length(x)){
    kx<-x[i]-min(x)+1
    ky<-y[i]-min(y)+1
    Ares[ky,kx]<-reslm[i]
    mask[ky,kx]<-1}
  filter1<-matrix(0,n,n)
  filter1[1,1]<-1
  leng<-length(reslm)
  ne<-convolve(convolve(Ares,filter1),Ares)[1,1]/leng
  n3<-3*n
  Ares0<-matrix(0,n3,n3)
  Ares1<-matrix(0,n3,n3)
  n1<-n+1
  nn<-n+n
  Ares0[n1:nn,n1:nn]<-Ares[1:n,1:n]
  Ares1[1:n,1:n]<-Ares[1:n,1:n]
  maske0<-matrix(0,n3,n3)
  maske1<-matrix(0,n3,n3)
  maske0[n1:nn,n1:nn]<-mask[1:n,1:n]
  maske1[1:n,1:n]<-mask[1:n,1:n]
  nx<-rep(1:n3,n3)
  ny<-as.numeric(gl(n3,n3))

  gr<-lim1
  gr1<-lim2
  h<-n*n3+n+1
  corr<-rep(0,dmax)
  kk<-0
  while(kk<dmax){
    kk<-kk+1
    filter<-matrix(0,n3,n3)
    for(i in 1:(n3*n3)) {
      d<-sqrt((nx[h]-nx[i])^2+(ny[h]-ny[i])^2)

      if(lim1!=0 & d>=gr & d<gr1) filter[ny[i],nx[i]]<-1
      if(lim1==0 & d>gr & d<=gr1) filter[ny[i],nx[i]]<-1}

    sum<-convolve(convolve(maske0,filter),maske1)[1,1]
    za<-convolve(convolve(Ares0,filter),Ares1)[1,1]/sum
    corr[kk]<-za/ne
    gr<-gr+(lim2-lim1)
    gr1<-gr1+(lim2-lim1)
  }
  corr<-as.vector(corr)
}



res.gee<-function(formula,family=gaussian,data,n,clusz=NA,zcor=NA,a=NA,b,
                  R=NA)  {
  ###############################################################################
  # Description
  # A function to calculate fitted values and residuals
  # for Generalized Estimating Equation Models
  # for gaussian or binary data (with logit link) or Poisson data (log link)
  # Arguments
  # formula     a formula expression
  # family      "gaussian", "binomial", "poisson" are allowed
  #             "binomial" = binary
  # data        a data frame
  # n           for maximal cluster size n*n
  # clusz       an object of function "clus.sz"
  # zcor        an object of function "genZcor"
  # a           a vector of correlation parameters
  #             for clusters only
  #             as an object of class "a.gee"
  # b           a vector of regression parameters beta
  # R           a square matrix of correlation parameters
  #             for full dimension (=number of observations)  only
  #
  # Value:     A list with components
  #            fitted    fitted values
  #            resid     normalized Pearson residuals
  ###############################################################################

  l<-dim(data)[2]
  ieo<-data[,l-1]
  if(n!=dim(data)[1]) {
    n2<-n*n
    n3<-n2*(n2-1)/2
    n4<-n2-1
    n5<-n2-2
    for(i in 1:dim(zcor)[1]){
      for(k in 1:n3){
        if(zcor[i,k]==1) zcor[i,k]<-a[k]  }}
    lc<-length(clusz)
    z2<-matrix(0,lc,n3)
    for( j in 1:n3) {
      k3<-0
      k2<-0
      for(i in 1:lc) {
        if(clusz[i]!=1) {
          k2<-k3+1
          k3<-clusz[i]*(clusz[i]-1)/2+k3
          for(k in k2:k3)
            z2[i,j]<-zcor[k,j]+z2[i,j] }}}
    if(n==2)
      iod<-c(1,1,2)
    if(n==3)
      iod<-c(1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,6,6,7)
    if(n==4)
      iod<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,
             3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,
             6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,9,9,9,9,9,9,10,10,10,10,10,
             11,11,11,11,12,12,12,13,13,14)
    cs<-0
    v<-matrix(0,length(ieo),length(ieo))
    vgl<-rep(0,n2)
    for(i in 1:lc) {clu<-clusz[i]
    if(clu!=1) {
      v1<-matrix(0,n2,n2)
      if(n==2)
      {  v1[1,2:4]<-z2[i,1:3]
      v1[2,3:4]<-z2[i,4:5]
      v1[3,4]<-z2[i,6]  }
      if(n==3)
      {  v1[1,2:9]<-z2[i,1:8]
      v1[2,3:9]<-z2[i,9:15]
      v1[3,4:9]<-z2[i,16:21]
      v1[4,5:9]<-z2[i,22:26]
      v1[5,6:9]<-z2[i,27:30]
      v1[6,7:9]<-z2[i,31:33]
      v1[7,8:9]<-z2[i,34:35]
      v1[8,9]<-z2[i,36]  }
      if(n==4)
      {  v1[1,2:16]<-z2[i,1:15]
      v1[2,3:16]<-z2[i,16:29]
      v1[3,4:16]<-z2[i,30:42]
      v1[4,5:16]<-z2[i,43:54]
      v1[5,6:16]<-z2[i,55:65]
      v1[6,7:16]<-z2[i,66:75]
      v1[7,8:16]<-z2[i,76:84]
      v1[8,9:16]<-z2[i,85:92]
      v1[9,10:16]<-z2[i,93:99]
      v1[10,11:16]<-z2[i,100:105]
      v1[11,12:16]<-z2[i,106:110]
      v1[12,13:16]<-z2[i,111:114]
      v1[13,14:16]<-z2[i,115:117]
      v1[14,15:16]<-z2[i,118:119]
      v1[15,16]<-z2[i,120]   }
      for(i1 in 1:length(iod)) {
        i2<-iod[i1]
        if(var(v1[i2,1:n2])==0) {for(k in i2:n5) {k1<-k+1
        v1[k,]<-v1[k1,]
        v1[k1,]<-vgl[]}}}
      for(i1 in 1:length(iod)){
        i3<-iod[i1]+1
        if(var(v1[1:n2,i3])==0) {for(k in i3:n4) {k1<-k+1
        v1[,k]<-v1[,k1]
        v1[,k1]<-vgl[]}}}

      clu1<-clu-1
      for(k in 1:clu1) {csk<-cs+k
      f1<-2
      for(k1 in f1:clu) {k2<-cs+f1
      v[csk,k2]<-v1[k,k1]
      f1<-f1+1 }}
      for(k in 1:clu) {csk<-cs+k
      v[csk,csk]<- 0.5 } }
    if(clu==1) {cs1<-cs+1
    v[cs1,cs1]<-0.5 }
    cs<- cumsum(clusz)[i]  }
    v<-v+t(v)
  }
  if(n==dim(data)[1]) v<-R
  ww<-solve(v)

  s.geese<-svd(ww)
  d.geese<-diag(sqrt(s.geese$d))
  w<-s.geese$u%*%d.geese%*%t(s.geese$u)

  x.matrix<-model.matrix(formula,data)
  fitted<-x.matrix%*%b
  fitted<-fitted[1:length(ieo)]
  if(family=="poisson") fitted<-exp(fitted)
  if(family=="binomial") fitted<-exp(fitted)/(1+exp(fitted))

  if(family=="gaussian") rgeese<- model.frame(formula,data)[[1]]-fitted
  if(family=="poisson")
    rgeese<-( model.frame(formula,data)[[1]]-fitted)/sqrt(fitted)
  if(family=="binomial")
    rgeese<-( model.frame(formula,data)[[1]]-fitted)/sqrt(fitted*(1-fitted))
  rsgeese<-w%*%rgeese
  resgeeseo<-rsgeese[1:length(ieo)]

  list(fitted=fitted,resid=resgeeseo)
}
