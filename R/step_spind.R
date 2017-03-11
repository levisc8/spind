#' Stepwise AIC, AICc, and QIC for GEEs and WRMs
#'
#' @param object A model of class \code{WRM} or \code{GEE}.
#' @param data The data used to fit that model.
#' @param steps Number of iterations the procedure should
#' go through before concluding. The default is to use the number of
#' variables as the number of iterations.
#' @param trace Should R print progress updates and the final, best model found
#' to the console? Default is \code{TRUE}.
#' @param AICc Logical. In the case of model selection with \code{WRM}s,
#' should AICc be used to determine which model is best rather than AIC?
#' This argument is ignored for \code{GEE}s. Default is \code{FALSE}.
#'
#' @return A list with components \code{model} and \code{table}.
#' \code{model} is always formula for the best model found by the procedure.
#' \code{table} is always a data frame, but the content varies for each type of
#' model.
#' For \code{WRM}'s, the columns
#' returned are
#' \itemize{
#'   \item\code{Deleted.Vars} Variables retained from the previous iteration
#'   which were tested in the current iteration.
#'   \item\code{LogLik} Log-likelihood of the model.
#'   \item\code{AIC} AIC score for the model.
#'   \item\code{AICc} AICc score for the model.
#' }
#'
#' For \code{GEE}s:
#' \itemize{
#'   \item\code{Deleted.Vars} Variables retained from the previous iteration
#'   which were tested in the current iteration.
#'   \item\code{QIC} Quasi-information criterion of the model.
#'   \item\code{Quasi.Lik} Quasi-likelihood of the model.
#' }
#'
#' @details This function performs stepwise variable elimination
#' for model comparison. Each iteration will try to find the best
#' combination of predictors for a given number of variables based
#' on AIC, AICc, or QIC, and then use that as the base model
#' for the next iteration until there are no more variables to eliminate.
#' Alternatively, it will terminate when reducing the number of variables
#' while respecting the model hierarchy no longer produces lower
#' information criterion values.
#'
#' @note Currently, the function only supports backwards model selection
#' (i.e. one must start with a full model and subtract variables).
#' Forward and both directions options may be added later.
#'
#' @references
#' Hardin, J.W. & Hilbe, J.M. (2003) Generalized Estimating Equations. Chapman and Hall, New York.
#'
#' @seealso \code{\link{qic.calc}}, \code{\link{aic.calc}}, \code{\link[stats]{drop1}},
#' \code{\link[stats]{step}}, \code{\link[MASS]{stepAIC}}
#'
#' @author Sam Levin
#'
#'
#' @examples
#' library(MASS)
#' data(birthwt)
#'
#'
#' x<-rep(1:14,14)
#' y<-as.integer(gl(14,14))
#' coords<-cbind(x[-(190:196)],y[-(190:196)])
#'
#' formula<-formula(low ~ age+ lwt+ race+ smoke+ ftv+  bwt)
#'
#' mgee<-GEE(formula, family = "gaussian", data = birthwt,
#'           coord=coords, corstr="fixed",scale.fix=TRUE)
#'
#' ss<-step.spind(mgee,birthwt)
#'
#' best.mgee<-GEE(ss$model, family = "gaussian", data = birthwt,
#'            coord=coords, corstr="fixed",scale.fix=TRUE)
#'
#' summary(best.mgee,printAutoCorPars=FALSE)
#'
#'
#' @export
#'


step.spind<-function (object,data,steps=NULL,trace=TRUE,AICc=FALSE){

  scope <- attr(terms(object$formula),
                "term.labels")
  model<-class(object)
  family<-object$family


  coord<-object$coord
  scale.fix<-object$scale.fix
  corstr<-object$corstr
  cluster<-object$cluster

  level<-object$level
  wavelet<-object$wavelet
  wtrafo<-object$wtrafo
  b.ini<-object$b.ini
  pad<-object$pad
  control<-object$control

  moran.params<-object$moran.params


  it<-1
  ns <- length(scope)
  polynomial.pattern<-'([I(])'
  interaction.pattern<-'([:])'
  poly.terms<-scope[stringr::str_detect(scope,polynomial.pattern)]
  inter.terms<-scope[stringr::str_detect(scope,interaction.pattern)]
  base.terms<-setdiff(scope,c(inter.terms,poly.terms))

  if(length(poly.terms)>0 & length(inter.terms)>0){
    term.col<-3
  } else if((length(poly.terms)|length(inter.terms))>0){
    term.col<-2
  } else if(length(poly.terms)==0 & length(inter.terms)==0){
    term.col<-1
  }

  terms<-matrix(nrow=length(base.terms),ncol=term.col)

  for(i in 1:length(base.terms)){
    fill<-scope[stringr::str_detect(scope,stringr::fixed(base.terms[i]))]
    terms[i,1:length(fill)]<-fill
  }

  if(model=="WRM"){
    ans <- matrix(nrow = ns + 1L, ncol = 3L,
                  dimnames = list(c("<none>", scope),
                                  c("loglik","inf.crit1","inf.crit2")))

    ans[1, ] <- c(object$LogLik,object$AIC,object$AICc)

    for (i in seq_len(ns)) {
      tt <- scope[i]

      nfit <- update.formula(object$formula, as.formula(paste("~ . -", tt)))
      newmod<-WRM(nfit,family,data,coord,level=level,
                  wavelet=wavelet,wtrafo=wtrafo,b.ini=b.ini,
                  pad=pad,control=control,moran.params=moran.params)
      ans[i + 1, ] <-c(newmod$LogLik,newmod$AIC,newmod$AICc)
    }
    aod <- data.frame(Deleted.Vars=rownames(ans),
                      LogLik=ans[, 1],
                      AIC = ans[, 2],
                      AICc=ans[, 3],
                      stringsAsFactors = F)

    rownames(aod)<-1:dim(aod)[1]
    if(AICc){
      best.mod<-aod$Deleted.Vars[which(aod$AICc==min(aod$AICc))]
    }else{
      best.mod<-aod$Deleted.Vars[which(aod$AIC==min(aod$AIC))]
    }
  }

  if(model=="GEE"){
    if(!scale.fix){
      scale.fix<-TRUE
      message("Scale parameter is now fixed")
    }
    ans <- matrix(nrow = ns + 1L, ncol = 2L,
                  dimnames = list(c("<none>", scope),
                                  c("inf.crit1","qlik")))

    ans[1, ] <- c(object$QIC,object$QLik)

    for (i in seq_len(ns)) {
      tt <- scope[i]

      nfit <- update.formula(object$formula, as.formula(paste("~ . -", tt)))
      newmod<-suppressWarnings({
        GEE(nfit,family,data,coord,corstr=corstr,
                  cluster=cluster,moran.params=moran.params,
                  scale.fix=scale.fix)
      })
      ans[i + 1, ] <-c(newmod$QIC,newmod$QLik)
    }
    aod <- data.frame(Deleted.Vars=rownames(ans),
                      QIC=ans[, 1],
                      Quasi.Lik = ans[, 2],
                      stringsAsFactors = F)

    rownames(aod)<-1:dim(aod)[1]
    best.mod<-aod$Deleted.Vars[which(aod$QIC==min(aod$QIC))]
  }
  if(trace){
    cat('Iteration: ',it,'\n','Single term deletions\n','Deleted Term: ',best.mod,
    '\n -------------------- \n')
    print(aod)
    cat('\n')
  }
  newvars<-setdiff(scope,best.mod)
  use.formula<-object$formula
  aod1<-aod
  newstart<-object$formula
  if(!is.null(steps)){
    steps<-steps
  } else{
    steps<-length(scope)
  }
  if(best.mod!='<none>'){
    while(it<=steps){
      usevars<-newvars
      aod1<-aod
      it<-it+1
      newstart<-update.formula(use.formula, as.formula(paste("~ . -", best.mod)))
      vars<-setdiff(usevars,best.mod)
      if(dim(terms)[2]==3){
        for(j in 1:dim(terms)[1]){
          if((terms[j,] %in% vars)[1]==FALSE &
             ((terms[j,] %in% vars)[2]|
              (terms[j,] %in% vars)[3])==TRUE){

            if(model=="GEE") aod1<-aod1[order(aod1$QIC),]
            if(model=="WRM" & !AICc) aod1<-aod1[order(aod1$AIC),]
            if(model=="WRM" & AICc) aod1<-aod1[order(aod1$AICc),]

            aod1<-aod1[-c(1),]

            new.best.mod<-aod1[1,"Deleted.Vars"]
            best.mod<-new.best.mod
            if(new.best.mod=='<none>') break

            newstart<-update.formula(use.formula,paste("~ . -",new.best.mod))
            cat('-----\nModel hierarchy violated by last removal\nNew Deleted Term: ',
            new.best.mod,'\nPreviously deleted term added back into model\n','-----\n')
          }
        }
      }

      if(dim(terms)[2]==2){
        for(j in 1:dim(terms)[1]){
          if((terms[j,] %in% vars)[1]==FALSE &
             ((terms[j,] %in% vars)[2])==TRUE){

            if(model=="GEE") aod1<-aod1[order(aod1$QIC),]
            if(model=="WRM" & !AICc) aod1<-aod1[order(aod1$AIC),]
            if(model=="WRM" & AICc) aod1<-aod1[order(aod1$AICc),]

            aod1<-aod1[-c(1),]

            new.best.mod<-aod1[1,"Deleted.Vars"]
            best.mod<-new.best.mod
            if(new.best.mod=='<none>') break

            newstart<-update.formula(use.formula,paste("~ . -",new.best.mod))
            cat('-----\nModel hierarchy violated by last removal\nNew deleted term: ',
                new.best.mod,'\nPreviously deleted term added back into model\n','-----\n')

          }
        }
      }



      ns<-length(attr(terms(newstart),'term.labels'))

      if(model=="WRM"){

        ans <- matrix(nrow = ns + 1L, ncol = 3L,
                      dimnames = list(c("<none>", attr(terms(newstart),'term.labels')),
                                      c("loglik","inf.crit1","inf.crit2")))

        newwrm<-WRM(newstart,family,data,coord,level=level,
                    wavelet=wavelet,wtrafo=wtrafo,b.ini=b.ini,
                    pad=pad,control=control,moran.params=moran.params)



        ans[1, ] <- c(newwrm$LogLik,newwrm$AIC,newwrm$AICc)


        for (i in seq_len(ns)) {
          tt <- attr(terms(newstart),'term.labels')[i]

          nfit <- update.formula(newstart, as.formula(paste("~ . -", tt)))
          newmod<-WRM(nfit,family,data,coord,level=level,
                      wavelet=wavelet,wtrafo=wtrafo,b.ini=b.ini,
                      pad=pad,control=control,moran.params=moran.params)
          ans[i + 1, ] <-c(newmod$LogLik,newmod$AIC,newmod$AICc)
        }
        aod <- data.frame(Deleted.Vars=rownames(ans),
                          LogLik=ans[, 1],
                          AIC = ans[, 2],
                          AICc=ans[, 3],
                          stringsAsFactors = F)

        rownames(aod)<-1:dim(aod)[1]
        if(AICc){
          best.mod<-aod$Deleted.Vars[which(aod$AICc==min(aod$AICc))]
        }else{
          best.mod<-aod$Deleted.Vars[which(aod$AIC==min(aod$AIC))]
        }
      }

      if(model=="GEE"){
        if(!scale.fix){
          scale.fix<-TRUE
        }
        ans <- matrix(nrow = ns + 1L, ncol = 2L,
                      dimnames = list(c("<none>", attr(terms(newstart),'term.labels')),
                                      c("inf.crit1","qlik")))

        newGEE<-suppressWarnings({
          GEE(newstart,family,data,coord,corstr=corstr,
                    cluster=cluster,moran.params=moran.params,
                    scale.fix=scale.fix)
        })
        ans[1, ] <- c(newGEE$QIC,newGEE$QLik)

        for (i in seq_len(ns)) {
          tt <- attr(terms(newstart),'term.labels')[i]

          nfit <- update.formula(newstart, as.formula(paste("~ . -", tt)))
          newmod<-suppressWarnings({
            GEE(nfit,family,data,coord,corstr=corstr,
                      cluster=cluster,moran.params=moran.params,
                      scale.fix=scale.fix)
          })
          ans[i + 1, ] <-c(newmod$QIC,newmod$QLik)
        }
        aod <- data.frame(Deleted.Vars=rownames(ans),
                          QIC=ans[, 1],
                          Quasi.Lik = ans[, 2],
                          stringsAsFactors = F)

        rownames(aod)<-1:dim(aod)[1]
        best.mod<-aod$Deleted.Vars[which(aod$QIC==min(aod$QIC))]

      }

      if(length(best.mod)>1){
        warning('Multiple equally parsimonious models')
      }


      aod1<-aod
      newvars<-setdiff(vars,best.mod)

      if(trace){
        cat('Iteration: ',it,'\n','Single term deletions\n','Deleted Term: ',best.mod,
            '\n -------------------- \n')
        print(aod)
        cat('\n')
      }

      suppressWarnings(if(best.mod!='<none>'){
        use.formula<-newstart
      })

      suppressWarnings(if(best.mod=='<none>') break)



    }
  }
  if(trace ){
    cat('\n---------------\nBest model found:\n')
    print(newstart)
  }

  return(list(model=newstart,
                table=aod))
}
