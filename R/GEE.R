#' @title  GEE (Generalized Estimating Equations)
#' @description
#' \code{GEE} provides GEE-based methods from the packages \pkg{gee} and \pkg{geepack}
#' to account for spatial autocorrelation in multiple linear regressions
#' @details
#' GEE can be used to fit linear models for response variables with
#' different distributions: \code{gaussian}, \code{binomial}, or \code{poisson}.
#' As a spatial model, it is a generalized linear model in which the residuals
#' may be autocorrelated. It accounts for spatial (2-dimensional)
#' autocorrelation of the residuals in cases of regular gridded datasets
#' and returns corrected parameter estimates. The grid cells are assumed to be square.
#' Futhermore, this function requires that \strong{all predictor variables
#' be continuous}.
#'
#' @param formula  Model formula. Variable names must match variables in \code{data}.
#' @param family \code{gaussian}, \code{binomial}, or \code{poisson} are supported.
#' Called using a quoted character string (i.e. \code{family} = "gaussian").
#' @param data  A data frame with variable names that match the variables
#' specified in \code{formula}.
#' @param coord    A matrix of two columns with corresponding cartesian
#' coordinates. Currently only supports integer coordinates.
#' @param  corstr   Expected autocorrelation structure: \code{independence}, \code{fixed},
#' \code{exchangeable}, and \code{quadratic}  are possible.
#'
#'  \itemize{
#'    \item\code{independence} - This is the same as a GLM, i.e. correlation matrix = identity matrix.
#'
#'    \item\code{fixed} - Uses an adapted isotropic power function specifying all correlation
#'    coefficients.
#'
#'    \item\code{exchangeable} and \code{quadratic} for clustering, i.e.
#'    the correlation matrix has a block diagonal form:
#'
#'    \itemize{
#'       \item\code{exchangeable} - All intra-block correlation coefficients are equal.
#'
#'       \item\code{quadratic} - Intra-block correlation coefficients for different
#'          distances can be different.
#'          }
#'        }
#' @param cluster  Cluster size for cluster models \code{exchangeable}
#'  and \code{quadratic}. Values of 2, 3, and 4 are allowed.
#'  \itemize{
#'    \item 2 - a 2*2 cluster
#'
#'    \item 3 - a 3*3 cluster
#'
#'    \item 4 - a 4*4 cluster
#' }
#'
#' @param moran.params    A list of parameters for calculating Moran's I.
#'   \itemize{
#'     \item\code{lim1} Lower limit for first bin. Default is 0.
#'     \item\code{increment} Step size for calculating I. Default is 1.
#'   }
#'
#' @param plot    A logical value indicating whether autocorrelation of
#' residuals should be plotted. NOW DEPRECATED in favor of \code{plot.GEE} method.
#'
#' @param scale.fix A logical indicating whether or not the scale parameter should
#' be fixed. The default is \code{FALSE}. Use \code{TRUE} when planning to use
#' stepwise model selection procedures in \code{step.spind}.
#' @param customize_plot Additional plotting parameters passed to \code{ggplot}.
#' NOW DEPRECATED in favor \code{plot.GEE} method.
#'
#'
#' @return An object of class \code{GEE}. This consists of a list with the
#' following elements:
#' \describe{
#'       \item{\code{call}}{Call}
#'       \item{\code{formula}}{Model formula}
#'       \item{\code{family}}{Family}
#'       \item{\code{coord}}{Coordinates used for the model}
#'       \item{\code{corstr}}{User-selected correlation structure}
#'       \item{\code{b}}{Estimate of regression parameters}
#'       \item{\code{s.e.}}{Standard errors of the estimates}
#'       \item{\code{z}}{Depending on the \code{family}, either a \emph{z} or \emph{t} value}
#'       \item{\code{p}}{\emph{p}-values for each parameter estimate}
#'       \item{\code{scale}}{Scale parameter (dispersion parameter) of the distribution's variance}
#'       \item{\code{scale.fix}}{Logical indicating whether \code{scale} has fixed value}
#'       \item{\code{cluster}}{User-specified cluster size for clustered models}
#'       \item{\code{fitted}}{Fitted values from the model}
#'       \item{\code{resid}}{Normalized Pearson residuals}
#'       \item{\code{w.ac}}{Working autocorrelation parameters}
#'       \item{\code{Mat.ac}}{Working autocorrelation matrix}
#'       \item{\code{QIC}}{Quasi Information Criterion. See \code{\link{qic.calc}}
#'        for further details}
#'       \item{\code{QLik}}{Quasi-likelihood. See \code{\link{qic.calc}}
#'        for further details}
#'       \item{\code{plot}}{Logical value indicating whether autocorrelation should
#'       be plotted}
#'       \item{\code{moran.params}}{Parameters for calculating Moran's I}
#'       \item{\code{v2}}{Parameter variance of the \code{GEE} model}
#'       \item{\code{var.naive}}{Paramter variance of the \code{independence} model}
#'       \item{\code{ac.glm}}{Autocorrelation of GLM residuals}
#'       \item{\code{ac.gee}}{Autocorrelation of GEE residuals}
#'       \item{\code{plot}}{An object of class \code{ggplot} containing information
#'       on the autocorrelation of residuals from the fitted \code{GEE} and a
#'       \code{GLM}}
#' }
#'
#' Elements can be viewed using the \code{\link{summary.GEE}} methods included in
#' the package.
#'
#' @note When using \code{corstr = "fixed"} on large data sets, the function
#' may return an error, as the resulting variance-covariance matrix is too
#' large for R to handle. If this happens, one will have to use one of the
#' cluster models (i.e \code{quadratic, exchangeable}).
#'
#' @seealso \code{\link{qic.calc}}, \code{\link{summary.GEE}}, \code{\link[gee]{gee}}
#'
#' @author Gudrun Carl, Sam Levin
#'
#'@examples
#' data(musdata)
#' coords<- musdata[,4:5]
#'
#' \dontrun{
#' mgee<-GEE(musculus ~ pollution + exposure, "poisson", musdata,
#'       coord=coords, corstr="fixed",scale.fix=FALSE)
#'
#' summary(mgee,printAutoCorPars=TRUE)
#'}
#' @references
#' Carl G & Kuehn I, 2007. Analyzing Spatial Autocorrelation in Species
#' Distributions using Gaussian and Logit Models, Ecol. Model. 207, 159 - 170
#'
#' Carey, V. J., 2006. Ported to R by Thomas Lumley (versions 3.13,
#' 4.4, version 4.13)., B. R. gee: Generalized Estimation Equation
#' solver. R package version 4.13-11.
#'
#' Yan, J., 2004. geepack: Generalized Estimating Equation Package.
#' R package version 0.2.10.
#'
#' @importFrom ggplot2 theme element_blank element_line element_text
#' ggplot aes_ geom_line geom_point scale_color_manual
#' scale_x_continuous scale_y_continuous
#' @importFrom gee gee
#' @importFrom geepack genZcor geese
#' @importFrom stats glm resid fitted dist pnorm
#' @importFrom utils capture.output
#' @importFrom rlang quo !!
#' @export
#'
#'
GEE <- function(formula,family,data,coord,
              corstr="fixed",cluster=3,moran.params=list(),
              plot=FALSE,scale.fix=FALSE, customize_plot = NULL){

  if(!is.null(customize_plot) | plot) {
    warning('"customize_plot" and "plot = TRUE" arguments are now soft deprecated.\n',
            'Use plot.GEE method and access the ggplot2 object using object_name$plot\n',
            'subsequent modification.')
  }

  at <- intersect(names(data), all.vars(formula))
  if(length(at) == 0) stop("formula: specified notation is missing")
  nn <- nrow(data)
  x <- coord[ ,1]
  y <- coord[ ,2]
  if(length(x) != nn) {
    stop("length of data does not match length of coordinates")
  }
  logic1 <- identical(as.numeric(x), round(x, 0))
  logic2 <- identical(as.numeric(y), round(y, 0))
  if(!logic1 | !logic2) stop("coordinates not integer")

  moran <- do.call("wrm.moran", moran.params)
  lim1 <- moran$lim1
  lim2 <- lim1 + moran$increment

  m0 <- stats::glm(formula, family, data)
  res0 <- stats::resid(m0, type = "pearson")
  id <- rep(1, nn)
  dato <- data.frame(data, id)
  suppressMessages(suppressWarnings(utils::capture.output({
    mgee <- gee::gee(formula = formula, family = family,
                     data = dato, id = id,
                     corstr = "independence", scale.fix = scale.fix)
    })))
  var.indep.naive <- mgee$naive.variance

  if(corstr == "independence"){
    ashort <- 0
    A <- 0
    fitted <- stats::fitted(m0)
    resid <- stats::resid(m0, type = "pearson")
    b <- summary(m0)$coefficients[ ,1]
    s.e. <- summary(m0)$coefficients[ ,2]
    z <- summary(m0)$coefficients[ ,3]
    p <- summary(m0)$coefficients[ ,4]
    scale <- summary(m0)$dispersion
    Icrit <- qic.calc(formula, family = family, data = data,
                      fitted, var.indep.naive, var.indep.naive)
    QIC <- Icrit$QIC
    logLik <- Icrit$loglik
    v2 <- var.indep.naive
  }

  if(corstr == "fixed"){
    ac01 <- acfft(coord, res0, 1, 1.1, dmax = 1)
    ac05 <- acfft(coord, res0, 5, 5.1, dmax = 1)
    if(ac05 <= 0) v <- 1
    if(ac05 > 0) v <- log(log(ac05) / log(ac01)) / log(5)
    alpha <- ac01
    para0 <- paste("n", nn,sep = "=")
    para1 <- paste(", alpha", round(alpha, 3), sep = "=")
    para2 <- paste(", v", round(v, 3), sep = "=")
    A0 <- paste(para0, para1, para2)
    id <- rep(1, nn)
    coord <- cbind(x, y)
    D <- as.matrix(stats::dist(coord))
    R <- alpha ^ (D^v)
    data <- data.frame(data, id)
  suppressMessages(suppressWarnings(utils::capture.output({
      mgee <- gee::gee(formula = formula, family = family,
                       data = data, id = id,R = R,corstr = "fixed",
                       scale.fix = scale.fix)
      })))
    var.naive <- mgee$naive.variance
    para3 <- "a=alpha^(d^v) "
    ashort <- c(alpha, v)
    A <- paste(para3, para1, para2)
    b <- mgee$coeff
    res <- res.gee(formula, family, data, nn, b = mgee$coeff, R = R)
    fitted <- res$fitted
    resid <- res$resid
    s.e. <- summary(mgee)$coefficients[ ,2]
    z <- summary(mgee)$coefficients[ ,3]
    p <- rep(NA, nrow(summary(mgee)$coefficients))
    for(ii in seq_len(nrow(summary(mgee)$coefficients))){
      if(z[ii] >= 0) p[ii] <- 2 * (1 - stats::pnorm(z[ii]))
      if(z[ii] < 0) p[ii] <- 2 * (stats::pnorm(z[ii]))
    }
    scale <- summary(mgee)[[9]]
    Icrit <- qic.calc(formula, family = family, data = data, fitted,
                      var.naive, var.indep.naive)
    QIC <- Icrit$QIC
    logLik <- Icrit$loglik
    v2 <- var.naive
   }


  if(corstr == "exchangeable") {
    dato <- dat.nn(data, coord, cluster)
    l <- dim(dato)[2]
    o <- dato[ ,l - 2]
    id <- dato[ ,l - 1]
    waves <- dato[ ,l]

    clusz <- clus.sz(id)
    zcor <- geepack::genZcor(clusz = clusz,
                             waves = waves, "unstr")
    suppressMessages(suppressWarnings(utils::capture.output({
      mgee <- gee::gee(formula = formula, family = family,
                     data = dato, id = id, corstr = "exchangeable",
                     scale.fix = scale.fix)
      })))
    var.robust <- mgee$robust.variance
    ashort <- mgee$w[1, 2]
    a <- a.gee(mgee$w, cluster, type = "gee", corstr = "exchangeable")
    A <- mgee$w
    b <- mgee$coeff
    res <- res.gee(formula, family, dato, cluster, clusz, zcor, a, b)
    fitted <- res$fitted[order(o)]
    resid <- res$resid[order(o)]
    s.e. <- summary(mgee)$coefficients[ ,4]
    z <- summary(mgee)$coefficients[ ,5]
    p <- rep(NA,nrow(summary(mgee)$coefficients))
    for(ii in seq_len(nrow(summary(mgee)$coefficients))) {
      if(z[ii] >= 0) p[ii] <- 2 * (1 - stats::pnorm(z[ii]))
      if(z[ii] < 0) p[ii] <- 2 * (stats::pnorm(z[ii]))
    }

    scale <- summary(mgee)[[9]]
    Icrit <- qic.calc(formula, family = family, data = data,
                      fitted, var.robust, var.indep.naive)
    QIC <- Icrit$QIC
    logLik <- Icrit$loglik
    v2 <- var.robust
   }


  if(corstr == "quadratic"){
    dato <- dat.nn(data, coord, cluster)
    l <- dim(dato)[2]
    o <- dato[ ,l - 2]
    id <- dato[ ,l - 1]
    waves <- dato[ ,l]

    clusz <- clus.sz(id)
    zcor <- geepack::genZcor(clusz = clusz, waves = waves, "unstr")
    zcorq <- zcor.quad(zcor, cluster, quad = TRUE)
    mgeese <- geepack::geese(formula = formula, family = family,
                             data = dato, id = id,
                             corstr = "userdefined", zcor = zcorq,
                             scale.fix = scale.fix)
    var.robust <- mgeese$vbeta
    ashort <- mgeese$a
    a <- a.gee(mgeese$a, cluster, type = "geese",
               corstr = "userdefined", quad = TRUE)
    A <- cor.mat(cluster, a)
    b <- mgeese$b
    res <- res.gee(formula, family, dato, cluster, clusz, zcor, a, b)
    fitted <- res$fitted[order(o)]
    resid <- res$resid[order(o)]
    s.e. <- summary(mgeese)$mean[ ,2]
    z <- summary(mgeese)$mean[ ,3]
    p <- summary(mgeese)$mean[ ,4]

    scale <- as.numeric(summary(mgeese)$scale[1])
    Icrit <- qic.calc(formula, family = family, data = data, fitted,
                      var.robust, var.indep.naive)
    QIC <- Icrit$QIC
    logLik <- Icrit$loglik
    v2 <- var.robust
  }


  ac0 <- acfft(coord, res0, lim1, lim2)
  ac <- acfft(coord, resid, lim1, lim2)

  plt.blank <-  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.background = ggplot2::element_blank(),
                               axis.line = ggplot2::element_line(colour = "black"),
                               legend.title = ggplot2::element_text(size = 9))

  plt.data <- data.frame(val = seq_len(length(ac)),
                         ac.gee = ac,
                         ac.glm = ac0)

  y.breaks <- round(seq(min(plt.data[ ,2:3])-.02,
                        max(plt.data[ ,2:3]) + .02,
                        length.out = 6), 2)

  val <- rlang::quo(val)
  ac.gee <- rlang::quo(ac.gee)
  ac.glm <- rlang::quo(ac.glm)

  plt <- ggplot2::ggplot(data = plt.data,
                         ggplot2::aes(x = !! val)) +
    plt.blank +
    ggplot2::geom_line(ggplot2::aes(y = !! ac.gee,
                                    color = "GEE Residuals"),
                       size = 0.9) +
    ggplot2::geom_line(ggplot2::aes(y = !! ac.glm,
                                    color = "GLM Residuals"),
                       size = 0.9) +
    ggplot2::geom_point(ggplot2::aes(y = !! ac.gee,
                                     color = "GEE Residuals"),
                        size = 2) +
    ggplot2::geom_point(ggplot2::aes(y = !! ac.glm,
                                     color = "GLM Residuals"),
                        size = 2) +
    ggplot2::scale_color_manual(paste("Correlation structure: ",
                                      corstr),
                                breaks = c('GEE Residuals','GLM Residuals'),
                                values = c('blue', 'red')) +
    ggplot2::scale_x_continuous('Lag Distance', breaks = 1:10) +
    ggplot2::scale_y_continuous("Autocorrelation of residuals",
                                breaks = y.breaks,
                                limits = c(min(plt.data[ ,2:3]) - .02,
                                           max(plt.data[ ,2:3]) + .02)) +
    customize_plot

  if(plot){
    print(plt)
  }

  call <- match.call()
  fit <- list(call = call,
              formula = formula,
              family = family,
              coords = coord,
              corstr = corstr,
              b = b,
              s.e. = s.e.,
              z = z, p = p,
              scale = scale,
              scale.fix = scale.fix,
              cluster = cluster,
              fitted = fitted,
              resid = resid,
              w.ac = ashort,
              Mat.ac = A,
              QIC = QIC,
              QLik = logLik,
              ac.glm = ac0,
              ac.gee = ac,
              var.gee = v2,
              var.naive = var.indep.naive,
              moran.params = moran.params,
              plot = plt)

  class(fit) <- "GEE"
  return(fit)


}


#' @name plot.GEE
#' @rdname GEE
#'
#'
#' @param x An object of class \code{GEE} or \code{WRM}
#' @param ... Not used.
#'
#'@examples
#' \dontrun{
#'library(ggplot2)
#' plot(mgee)
#'
#' my_gee_plot <- mgee$plot
#'
#' # move the legend to a new position
#' print(my_gee_plot + ggplot2::theme(legend.position = 'top'))
#'
#'}
#'
#' @export

plot.GEE <- function(x, ...) {

  print(x$plot)

  invisible(x)

}

#' @name predict.GEE
#' @rdname GEE
#'
#' @inheritParams summary.GEE
#' @param newdata  A data frame containing variables to base the predictions on.
#' @examples
#' data(musdata)
#' coords<-musdata[,4:5]
#' mgee<-GEE(musculus ~ pollution + exposure,'poisson',musdata,
#'           coord=coords,corstr="fixed")
#'
#' pred<-predict(mgee,newdata=musdata)
#'
#'
#'@importFrom stats model.matrix
#'@export
predict.GEE<-function(object,newdata,...){


  data<-newdata
  formula<-object$formula
  family<-object$family
  b<-object$b

  x.matrix<-stats::model.matrix(formula,data)
  fitted<-x.matrix%*%b
  fitted<-as.vector(fitted)
  if(family=="poisson") fitted<-exp(fitted)
  if(family=="binomial") fitted<-exp(fitted)/(1+exp(fitted))

  return(fitted)
}

#' @name summary.GEE
#' @rdname GEE
#'
#' @param object An object of class \code{GEE}.
#' @param printAutoCorPars A logical indicating whether to print the
#' working autocorrelation parameters
#' @inheritParams plot.GEE
#'
#' @importFrom stats printCoefmat
#' @export

summary.GEE<-function(object,...,printAutoCorPars=TRUE){

  cat("\n","Call:","\n")
  print(object$call)
  family<-object$family
  b<-object$b
  s.e.<-object$s.e.
  z<-object$z
  p<-object$p
  QIC<-object$QIC
  beta<-cbind(b,s.e.,z,p)
  if(family=="gaussian")
    colnames(beta) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
  if(family=="binomial" | family=="poisson")
    colnames(beta) <- c("Estimate", "Std.Err", "z value", "Pr(>|z|)")
  cat("---","\n","Coefficients:","\n")
  stats::printCoefmat(beta)
  cat("---","\n","QIC: ",QIC,"\n" )
  cat("---","\n")
  ac0<-object$ac.glm
  acg<-object$ac.gee
  cat("Autocorrelation of GLM residuals","\n")
  print(ac0)
  cat("\n","Autocorrelation of GEE residuals","\n")
  print(acg)

  if(printAutoCorPars & object$corstr != "independence"){
    cat('---','\n','Autocorrelation parameters from ',
        object$corstr," model",'\n')
    print(object$Mat.ac)
  }
}





#' @title Quasi-Information Criterion for Generalized Estimating
#' Equations
#'
#' @description
#' A function for calculating quasi-likelihood and Quasi-Information
#' Criterion values based on the method of Hardin & Hilbe (2003).
#' @param  formula  a model formula
#' @param  family   \code{gaussian}, \code{binomial}, or \code{poisson}
#' @param  data     a data frame
#' @param  mu       fitted values from a model
#' @param  var.robust        variance of model parameters
#' @param  var.indep.naive   naive variance of model parameters under the
#' \code{independence} model
#'
#' @return  A list with the following components:
#'  \describe{
#'    \item{\code{QIC}}{quasi-information criterion}
#'    \item{\code{loglik}}{quasi-likelihood}
#'  }
#' @references
#' Hardin, J.W. & Hilbe, J.M. (2003) Generalized Estimating Equations. Chapman and Hall, New York.
#'
#' Barnett et al. Methods in Ecology & Evolution 2010, 1, 15-24.
#' @importFrom stats model.matrix model.frame
#' @export
qic.calc <- function(formula, family, data, mu, var.robust, var.indep.naive){

  X <- stats::model.matrix(formula, data)
  if(is.vector(stats::model.frame(formula, data)[[1]])){
    y <- stats::model.frame(formula, data)[[1]]
    ntr <- 1
  }
  if(family == "binomial" & is.matrix(stats::model.frame(formula, data)[[1]])){
    y <- stats::model.frame(formula, data)[[1]][ ,1]
    ntr <- stats::model.frame(formula, data)[[1]][ ,1] +
           stats::model.frame(formula, data)[[1]][ ,2]
  }
  n  <-  dim(X)[1]
  nvar <- dim(X)[2]

  if(family == "gaussian"){
    sos <- sum((y - mu)^2)                    # sum of squares
    sigma2 <-  sos / n                        # variance
    loglik <-  -n / 2 * (log(2 * pi * sigma2) + 1)   # log likelihood
  }
  if(family == "binomial"){                         # choose= Binomialkoeff.
    loglik <-  sum(y * log(mu / (1 - mu)) +
                     ntr * log(1 - mu) + log(choose(ntr, y)))
  }
  if(family == "poisson"){
    # loglik <-  sum(y*log(mu)-mu)  # useful for delta in multimodel inference
    loglik <-  sum(y * log(mu) - mu) - sum(log(factorial(y)))
  }

  trace <- sum(diag(MASS::ginv(var.indep.naive) %*% var.robust))

  QIC <-  -2 * loglik + 2 * trace

  return(list(QIC = QIC, loglik = loglik))

}

dat.nn <- function(data,coord,n){
  ###########################################################################
  # Description
  # A function to generate clusters and order variables
  # Arguments
  # data      a data frame
  # coord     a matrix of two columns with corresponding coordinates
  # n         for maximal cluster size  n*n
  #
  # Value: a new data frame containing rearranged data, rearranged coordinates
  # and 3 new parameters:
  #    o      order parameter
  #    id     parameter identifying clusters
  #    waves  parameter identifying members of clusters
  #
  ############################################################################

  l <- dim(data)[2]
  OST <- coord[,1]
  NORD <- coord[,2]
  ko <- OST-min(OST)
  idx <- (ko-(ko%%(n)))/n+1
  ks <- NORD-min(NORD)
  idy <- (ks-(ks%%(n)))/n+1
  ie <- (idy-1)*max(idx)+idx
  idwx <- ko%%(n)+1
  idwy <- ks%%(n)+1
  wav <- (idwy-1)*n+idwx
  data <- as.matrix(data)
  o <- order(ie,wav)
  x <- OST[o]
  y <- NORD[o]
  id <- ie[o]
  waves <- wav[o]
  dat.new1 <- data[o,]
  dat.new2 <- cbind(dat.new1,x,y,o,id,waves)
  dat.new <- as.data.frame(dat.new2)
}


clus.sz <- function(id){
  ########################################################################
  # Description
  # A function to calculate sizes of clusters
  # Argument
  # id     a vector that identifies clusters
  # Value:
  # A vector of numbers of cluster sizes
  ########################################################################

  clus <- rep(0,length(id))
  k0 <- 0
  k1 <- 1
  for(i in 2:length(id)) {
    i1 <- i-1
    if(id[i]==id[i1]) {
      k1 <- k1+1
      if(i==length(id)) {
        k0 <- k0+1
        clus[k0] <- k1
        }
      }
   if(id[i]!=id[i1]) {
      k0 <- k0+1
      clus[k0] <- k1
      k1 <- 1
      if(i==length(id)) {
        k0 <- k0+1
        clus[k0] <- k1
      }
    }
  }
  clusz <- clus[clus>0]
}


zcor.quad <- function(zcor,n,quad=TRUE) {
  #########################################################################
  # Description
  # A function to create a quadratic correlation structure
  # Arguments:
  # zcor    an object of class "genZcor" (see: geepack)
  # n       for maximal cluster size n*n
  # quad    by default quadratic correlation structure
  # Value:
  # A matrix describing the quadratic correlation structure
  #########################################################################

  if(quad) {
    if(n==2)  {
      zcorn <- matrix(0,dim(zcor)[1],2)
      zcorn[,1] <- zcor[,1]+zcor[,2]+zcor[,5]+zcor[,6]
      zcorn[,2] <- zcor[,3]+zcor[,4]
    }
    if(n==3)  {
      zcorn <- matrix(0,dim(zcor)[1],5)
      zcorn[,1] <- zcor[,1]+zcor[,3]+zcor[,9]+zcor[,11]+zcor[,18]+zcor[,22]+
        zcor[,24]+zcor[,27]+zcor[,29]+zcor[,33]+zcor[,34]+zcor[,36]
      zcorn[,2] <- zcor[,2]+zcor[,6]+zcor[,14]+zcor[,21]+zcor[,23]+zcor[,35]
      zcorn[,3] <- zcor[,4]+zcor[,10]+zcor[,12]+zcor[,17]+zcor[,25]+zcor[,28]+
        zcor[,30]+zcor[,32]
      zcorn[,4] <- zcor[,5]+zcor[,7]+zcor[,13]+zcor[,15]+zcor[,16]+zcor[,20]+
        zcor[,26]+zcor[,31]
      zcorn[,5] <- zcor[,8]+zcor[,19]
    }
    if(n==4)  {
      zcorn <- matrix(0,dim(zcor)[1],9)
      zcorn[,1] <- zcor[,1]+zcor[,4]+zcor[,16]+zcor[,19]+zcor[,30]+zcor[,33]+
        zcor[,46]+zcor[,55]+zcor[,58]+zcor[,66]+zcor[,69]+zcor[,76]+
        zcor[,79]+zcor[,88]+zcor[,93]+zcor[,96]+zcor[,100]+zcor[,103]+
        zcor[,106]+zcor[,109]+zcor[,114]+zcor[,115]+zcor[,118]+zcor[,120]
      zcorn[,2] <- zcor[,2]+zcor[,8]+zcor[,17]+zcor[,23]+zcor[,37]+zcor[,50]+
        zcor[,56]+zcor[,62]+zcor[,67]+zcor[,73]+zcor[,83]+zcor[,92]+
        zcor[,94]+zcor[,101]+zcor[,116]+zcor[,119]
      zcorn[,3] <- zcor[,3]+zcor[,12]+zcor[,27]+zcor[,41]+zcor[,54]+zcor[,57]+
        zcor[,95]+zcor[,117]
      zcorn[,4] <- zcor[,5]+zcor[,18]+zcor[,20]+zcor[,32]+zcor[,34]+zcor[,45]+
        zcor[,59]+zcor[,68]+zcor[,70]+zcor[,78]+zcor[,80]+zcor[,87]+
        zcor[,97]+zcor[,102]+zcor[,104]+zcor[,108]+zcor[,110]+zcor[,113]
      zcorn[,5] <- zcor[,6]+zcor[,9]+zcor[,21]+zcor[,22]+zcor[,24]+zcor[,31]+
        zcor[,36]+zcor[,38]+zcor[,44]+zcor[,49]+zcor[,60]+zcor[,63]+
        zcor[,71]+zcor[,72]+zcor[,74]+zcor[,77]+zcor[,82]+zcor[,84]+
        zcor[,86]+zcor[,91]+zcor[,98]+zcor[,105]+zcor[,107]+zcor[,112]
      zcorn[,6] <- zcor[,7]+zcor[,13]+zcor[,26]+zcor[,28]+zcor[,40]+zcor[,42]+
        zcor[,43]+zcor[,53]+zcor[,61]+zcor[,85]+zcor[,99]+zcor[,111]
      zcorn[,7] <- zcor[,10]+zcor[,25]+zcor[,35]+zcor[,48]+zcor[,64]+zcor[,75]+
        zcor[,81]+zcor[,90]
      zcorn[,8] <- zcor[,11]+zcor[,14]+zcor[,29]+zcor[,39]+zcor[,47]+zcor[,52]+
        zcor[,65]+zcor[,89]
      zcorn[,9] <- zcor[,15]+zcor[,51]
    }
  }
  if(!quad) zcorn <- zcor
  zcorn <- as.matrix(zcorn)
}


a.gee <- function(mgee,n,type="glm",corstr="independence",quad=TRUE) {
  ########################################################################
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
  #########################################################################

  if(n==2)n3 <- 6
  if(n==3)n3 <- 36
  if(n==4)n3 <- 120
  a <- rep(0,n3)
  if(type=="glm") a <- a
  if(type=="gee"){
    if(corstr=="exchangeable") a[c(1:n3)] <- mgee[1,2]
    if(corstr=="independence") a <- a
  }
  a <- as.vector(a)

  if(type=="geese") {
    if(corstr=="userdefined"){
      if(quad) {
        if(n==2)  {
          a <- rep(0,6)
          a[c(1,2,5,6)] <- mgee[1]
          a[c(3,4)] <- mgee[2]
        }
        if(n==3)  {
          a <- rep(0,36)
          a[c(1,3,9,11,18,22,24,27,29,33,34,36)] <- mgee[1]
          a[c(2,6,14,21,23,35)] <- mgee[2]
          a[c(4,10,12,17,25,28,30,32)] <- mgee[3]
          a[c(5,7,13,15,16,20,26,31)] <- mgee[4]
          a[c(8,19)] <- mgee[5]
        }
        if(n==4)  {
          a <- rep(0,120)
          a[c(1,4,16,19,30,33,46,55,58,66,69,76,79,88,93,96,100,103,106,109,
              114,115,118,120)] <- mgee[1]
          a[c(2,8,17,23,37,50,56,62,67,73,83,92,94,101,116,119)] <- mgee[2]
          a[c(3,12,27,41,54,57,95,117)] <- mgee[3]
          a[c(5,18,20,32,34,45,59,68,70,78,
              80,87,97,102,104,108,110,113)] <- mgee[4]
          a[c(6,9,21,22,24,31,36,38,44,49,60,63,71,72,74,77,82,84,86,91,98,
              105,107,112)] <- mgee[5]
          a[c(7,13,26,28,40,42,43,53,61,85,99,111)] <- mgee[6]
          a[c(10,25,35,48,64,75,81,90)] <- mgee[7]
          a[c(11,14,29,39,47,52,65,89)] <- mgee[8]
          a[c(15,51)] <- mgee[9]
        }}
      if(!quad) a <- mgee
    }
    if(corstr=="exchangeable") a[c(1:n3)] <- mgee
    if(corstr=="independence") a <- a
  }
  a <- as.vector(a)
}


cor.mat <- function(cluster,a) {
  ######################################################################
  # Description
  # A function to create a block of the correlation matrix
  # Arguments:
  # cluster     cluster size
  # a           cluster parameter
  # Value:
  # A matrix representing a block of the correlation matrix
  ######################################################################
  n <- cluster
  n2 <- cluster*cluster
  z2 <- a
  v1 <- matrix(0,n2,n2)
  if(n==2)
  {  v1[1,2:4] <- z2[1:3]
  v1[2,3:4] <- z2[4:5]
  v1[3,4] <- z2[6]  }
  if(n==3)
  {  v1[1,2:9] <- z2[1:8]
  v1[2,3:9] <- z2[9:15]
  v1[3,4:9] <- z2[16:21]
  v1[4,5:9] <- z2[22:26]
  v1[5,6:9] <- z2[27:30]
  v1[6,7:9] <- z2[31:33]
  v1[7,8:9] <- z2[34:35]
  v1[8,9] <- z2[36]  }
  if(n==4)
  {  v1[1,2:16] <- z2[1:15]
  v1[2,3:16] <- z2[16:29]
  v1[3,4:16] <- z2[30:42]
  v1[4,5:16] <- z2[43:54]
  v1[5,6:16] <- z2[55:65]
  v1[6,7:16] <- z2[66:75]
  v1[7,8:16] <- z2[76:84]
  v1[8,9:16] <- z2[85:92]
  v1[9,10:16] <- z2[93:99]
  v1[10,11:16] <- z2[100:105]
  v1[11,12:16] <- z2[106:110]
  v1[12,13:16] <- z2[111:114]
  v1[13,14:16] <- z2[115:117]
  v1[14,15:16] <- z2[118:119]
  v1[15,16] <- z2[120]   }
  for(i in 1:n2) v1[i,i] <- 0.5
  v <- v1+t(v1)
  v
}


#' @importFrom stats model.matrix model.frame
#' @importFrom stats var
res.gee <- function(formula,family='gaussian',data,n,clusz=NA,zcor=NA,a=NA,b,
                  R=NA)  {
  #######################################################################
  # Description
  # A function to calculate fitted values and residuals
  # for Generalized Estimating Equation Models
  # for gaussian or binary data (with logit link) or Poisson data (log link)
  # Arguments
  # formula     a formula expression
  # family      "gaussian", "binomial", "poisson" are allowed
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
  #########################################################################

  l <- dim(data)[2]
  ieo <- data[,l-1]
  if(n!=dim(data)[1]) {
    n2 <- n*n
    n3 <- n2*(n2-1)/2
    n4 <- n2-1
    n5 <- n2-2
    for(i in 1:dim(zcor)[1]){
      for(k in 1:n3){
        if(zcor[i,k]==1) zcor[i,k] <- a[k]  }}
    lc <- length(clusz)
    z2 <- matrix(0,lc,n3)
    for( j in 1:n3) {
      k3 <- 0
      k2 <- 0
      for(i in 1:lc) {
        if(clusz[i]!=1) {
          k2 <- k3+1
          k3 <- clusz[i]*(clusz[i]-1)/2+k3
          for(k in k2:k3)
            z2[i,j] <- zcor[k,j]+z2[i,j] }}}
    if(n==2)
      iod <- c(1,1,2)
    if(n==3)
      iod <- c(1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,6,6,7)
    if(n==4)
      iod <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,2,2,2,2,2,2,
               3,3,3,3,3,3,3,3,3,3,3,3,
               4,4,4,4,4,4,4,4,4,4,4,
               5,5,5,5,5,5,5,5,5,5,
               6,6,6,6,6,6,6,6,6,
               7,7,7,7,7,7,7,7,
               8,8,8,8,8,8,8,
               9,9,9,9,9,9,
               10,10,10,10,10,
               11,11,11,11,
               12,12,12,
               13,13,
               14)
    cs <- 0
    v <- matrix(0,length(ieo),length(ieo))
    vgl <- rep(0,n2)
    for(i in 1:lc) {
      clu <- clusz[i]
      if(clu!=1) {
        v1 <- matrix(0,n2,n2)
        if(n==2) {
          v1[1,2:4] <- z2[i,1:3]
          v1[2,3:4] <- z2[i,4:5]
          v1[3,4] <- z2[i,6]
        }
        if(n==3){
          v1[1,2:9] <- z2[i,1:8]
          v1[2,3:9] <- z2[i,9:15]
          v1[3,4:9] <- z2[i,16:21]
          v1[4,5:9] <- z2[i,22:26]
          v1[5,6:9] <- z2[i,27:30]
          v1[6,7:9] <- z2[i,31:33]
          v1[7,8:9] <- z2[i,34:35]
          v1[8,9] <- z2[i,36]
        }
        if(n==4){
          v1[1,2:16] <- z2[i,1:15]
          v1[2,3:16] <- z2[i,16:29]
          v1[3,4:16] <- z2[i,30:42]
          v1[4,5:16] <- z2[i,43:54]
          v1[5,6:16] <- z2[i,55:65]
          v1[6,7:16] <- z2[i,66:75]
          v1[7,8:16] <- z2[i,76:84]
          v1[8,9:16] <- z2[i,85:92]
          v1[9,10:16] <- z2[i,93:99]
          v1[10,11:16] <- z2[i,100:105]
          v1[11,12:16] <- z2[i,106:110]
          v1[12,13:16] <- z2[i,111:114]
          v1[13,14:16] <- z2[i,115:117]
          v1[14,15:16] <- z2[i,118:119]
          v1[15,16] <- z2[i,120]
        }
        for(i1 in seq_len(length(iod))) {
          i2 <- iod[i1]
          if(stats::var(v1[i2,1:n2])==0) {
            for(k in i2:n5) {
              k1 <- k+1
              v1[k,] <- v1[k1,]
              v1[k1,] <- vgl[]
            }
          }
        }
        for(i1 in seq_len(length(iod))) {
          i3 <- iod[i1]+1
          if(stats::var(v1[1:n2,i3])==0) {
            for(k in i3:n4) {
              k1 <- k+1
              v1[,k] <- v1[,k1]
              v1[,k1] <- vgl[]
            }
          }
        }

        clu1 <- clu-1
        for(k in seq_len(clu1)) {
          csk <- cs+k
          f1 <- 2
          for(k1 in f1:clu) {
            k2 <- cs+f1
            v[csk,k2] <- v1[k,k1]
            f1 <- f1+1
          }
        }
        for(k in seq_len(clu)) {
          csk <- cs+k
          v[csk,csk] <-  0.5
        }
      }
      if(clu==1) {cs1 <- cs+1
      v[cs1,cs1] <- 0.5 }
      cs <-  cumsum(clusz)[i]  }
    v <- v+t(v)
  }
  if(n == dim(data)[1]) v <- R
  ww <- solve(v)

  s.geese <- svd(ww)
  d.geese <- diag(sqrt(s.geese$d))
  w <- s.geese$u %*% d.geese %*% t(s.geese$u)

  x.matrix <- stats::model.matrix(formula,data)
  fitted <- x.matrix %*% b
  fitted <- fitted[seq_len(length(ieo))]
  if(family=="poisson") fitted <- exp(fitted)
  if(family=="binomial") fitted <- exp(fitted)/(1+exp(fitted))

  if(family=="gaussian") {
    rgeese <-  stats::model.frame(formula,data)[[1]]-fitted
  }
  if(family=="poisson")
    rgeese <- ( stats::model.frame(formula,
                                   data)[[1]]-fitted)/sqrt(fitted)
  if(family=="binomial"){
    rgeese <- ( stats::model.frame(formula,
                                   data)[[1]]-fitted)/sqrt(fitted*(1-fitted))
  }
  rsgeese <- w %*% rgeese
  resgeeseo <- rsgeese[seq_len(length(ieo))]

  list(fitted=fitted,resid=resgeeseo)
}

