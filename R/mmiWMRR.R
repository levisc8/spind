#' @title Multi-model inference for wavelet multiresolution regression
#'
#'
#' @description
#' mmiWMRR is a multimodel inference approach evaluating the relative
#' importance of predictors used in \code{\link{scaleWMRR}}.
#' @details It performs automatically
#' generated model selection and creates a model
#' selection table according to the approach of multi-model inference
#' (Burnham & Anderson, 2002). The analysis is carried out for scale-specific
#' regressions (i.e. where \code{\link{scaleWMRR}} can be used). AIC is
#' used to obtain model
#' selection weights and to rank the models.
#' Futhermore, this function requires that \strong{all predictor variables
#' be continuous}.
#'
#'
#' @param object A model of class \code{WRM}.
#' @param data     Data frame.
#' @param scale    0 or higher integers possible (limit depends on sample size).
#' \code{scale}=1 is equivalent to \code{WRM} with \code{level}=1.
#' @param detail   Remove smooth wavelets? If \code{TRUE}, only detail
#' components are analyzed. If set to \code{FALSE}, smooth and detail
#' components are analyzed. Default is \code{TRUE}.
#' @param trace Logical value indicating whether to print results to console.
#'
#' @return  \code{mmiWMRR} returns a list containing the following elements
#' \describe{
#'   \item{\code{result}}{A matrix containing slopes, degrees of freedom, likelihood,
#'          AIC, delta, and weight values for the set of candidate models.
#'          The models are ranked by AIC.}
#'   \item{\code{level}}{An integer corresponding to scale}
#'
#'}
#' @seealso  \code{\link{aic.calc}}, \code{\link{rvi.plot}},
#' \pkg{MuMIn}, \code{\link{WRM}}
#'
#' @author Gudrun Carl
#'
#' @examples
#'
#' data(carlinadata)
#' coords <- carlinadata[ ,4:5]
#' \dontrun{
#'
#'
#' wrm<- WRM(carlina.horrida ~ aridity + land.use, family = "poisson",
#'           data = carlinadata, coord = coords, level = 1,
#'           wavelet = "d4")
#'
#' mmi<- mmiWMRR(wrm, data = carlinadata, scale = 3,
#'               detail = TRUE, trace = FALSE)
#'
#'}
#' @references
#' Burnham, K.P. & Anderson, D.R. (2002) Model selection and
#' multimodel inference. Springer, New York.
#'
#' Carl G, Doktor D, Schweiger O, Kuehn I (2016)
#' Assessing relative variable importance across different spatial
#' scales: a two-dimensional wavelet analysis.
#' Journal of Biogeography 43: 2502-2512.
#' @importFrom stats update reformulate drop.terms model.matrix as.formula terms
#' @importFrom rje powerSetMat
#' @export

mmiWMRR <- function(object, data, scale, detail = TRUE, trace = FALSE){

  family <- object$family
  formula <- object$formula
  coord <- object$coord
  wavelet <- object$wavelet
  wtrafo <- object$wtrafo
  n.eff <- object$n.eff


  # Parameter: varnames, p
  X <- stats::model.matrix(formula, data)
  if(dimnames(X)[[2]][1] != "(Intercept)") {
    formula <- update(formula, ~ . + 1)
    X <- stats::model.matrix(formula, data)
  }
  nvar <- dim(X)[2]
  varnames <- dimnames(X)[[2]][-1]
  p <- dim(X)[2]-1 # nvar-1 (without intercept)

  pset <- rje::powerSetMat(p)
  ip <- dim(pset)[1]
  t <- stats::terms(formula)
  # Run every model and calculate AIC (multimodel inference)
  coef.vec <- matrix(NA,ip,nvar)
  df <- rep(NA, ip)
  loglik <- rep(NA, ip)
  AIC <- rep(NA, ip)
  for (i in 1:ip) {
    if(sum(pset[i, ]) != 0 & sum(pset[i, ]) != p){
      t1 <- stats::drop.terms(t, which(pset[i, ] == 0), keep.response = TRUE)
      formula1 <- stats::reformulate(attr(t1, "term.labels"), formula[[2]])
      formulae <- formula1
    }
    if(sum(pset[i, ]) == p) formulae <- formula
    if(sum(pset[i, ]) == 0) formulae <- stats::as.formula(paste(formula[[2]],
                                                                "~1"))
    m0 <- scaleWMRR(formulae, family, data, coord, scale = scale,
                    detail = detail, wavelet = wavelet, wtrafo = wtrafo,
                    pad = list(padzone = 1.1))
    if(!m0$converged)
      m0 <- scaleWMRR(formulae, family, data, coord, scale = scale,
                    detail = detail, wavelet = wavelet, wtrafo = wtrafo,
                    b.ini = glm(formulae, family, data)$coef,
                    pad = list(padzone = 1.1))
    kv <- c(1, which(pset[i, ] == 1) + 1)
    coef.vec[i, kv] <- m0$b
    mu <- m0$fitted
    if(is.null(n.eff)) aic <-aic.calc(formulae, family, data, mu)
    if(!is.null(n.eff)) aic <- aic.calc(formulae, family, data, mu, n.eff)
    df[i] <- aic$df
    loglik[i] <- aic$loglik
    AIC[i] <- aic$AIC
  }
  result <- cbind(round(coef.vec, 5), df, round(loglik, 3), round(AIC, 1))
  # Calculate delta and weight (multimodel inference)
  delta <- AIC - min(AIC) # = delta
  weight <- exp(-delta / 2) / sum(exp(-delta / 2)) # = weight
  # Print results
  result <- cbind(result, round(delta, 2), round(weight, 3))
  ord <- order(delta)
  res <- result[ord, ]
  dimnames(res)[[1]] <- ord
  dimnames(res)[[2]] <- c("(Int)", varnames,
                        "df", "logLik", "AIC", "delta", "weight")
  if(!detail & scale >= 1) scale <- scale - 1
  if(trace){
    cat("---","\n", "Level = ", scale, "\n")
    print(res, na.print = "")
  }
  fit <- list(result = res,
              level = scale)
  fit
}
