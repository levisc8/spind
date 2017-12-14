#'@title Spatial threshold-dependent accuracy measures
#'
#'@description Calculates spatially corrected, threshold-dependent metrics for
#'an observational data set and model predictions (Kappa and confusion matrix)
#'
#'@param data A data frame or matrix with two columns. The first column
#'should contain actual presence/absence data (binary, 0 or 1) and the
#'second column should contain model predictions of probability of
#'occurence (numeric, between 0 and 1).
#'@param coord A data frame or matrix with two columns containing x,y
#'coordinates for each actual and predicted value. Coordinates must be
#'integer and consecutively numbered.
#'@param thresh A cutoff value for classifying predictions as modeled
#'presences or modeled absences. Default is 0.5.
#'@param spatial A logical indicating whether spatially corrected indices
#'(rather than classical indices) should be computed.
#'
#'@return A list with the following components:
#'\describe{
#'  \item{\code{kappa}}{Kappa statistic}
#'  \item{\code{cm}}{Confusion matrix}
#'  \item{\code{sensitivity}}{Sensitivity}
#'  \item{\code{specificity}}{Specificity}
#'  \item{\code{actuals}}{Actual occurence data or adjusted actual occurence data}
#'  \item{\code{splitlevel.pred}}{Level splitting of predicted values}
#'  \item{\code{splitlevel.act}}{Level splitting of actuals or adjusted actuals}
#'  \item{\code{splitposition.pred}}{Position splitting of predicted values}
#'  \item{\code{splitposition.act}}{Position splitting of actuals or adjusted actuals}
#'}
#'
#'@author Gudrun Carl
#'
#'@seealso \code{\link{th.indep}}
#'
#'@references Carl G, Kuehn I (2017) Spind: a package for computing
#' spatially corrected accuracy measures.
#' Ecography 40: 675-682. doi: 10.1111/ecog.02593
#'
#'@examples
#'data(hook)
#'data <- hook[ ,1:2]
#'coord <- hook[ ,3:4]
#'si1 <- th.dep(data, coord, spatial = TRUE)
#'si1$kappa
#'si1$cm
#'
#'@export


th.dep<-function(data,coord,thresh=0.5,spatial=TRUE){

  if(dim(data)[1] != dim(coord)[1]) stop("error in dimension")

  if(spatial){
    y <- adjusted.actuals(data, coord)
    split <- 4
  }
  if(!spatial){
    y <- data[ ,1]
    split <- 2
  }

  pi <- data[ ,2]
  n <- length(pi)

  splitlevel <- matrix(NA, split, 2)
  splitlevely <- matrix(NA, split, 2)

  if(split==4) {
    lower.split <- thresh / 2
    upper.split <- (1 + thresh) / 2

    splitlevel[1, ] <- c(0, lower.split)
    splitlevel[2, ] <- c(lower.split, thresh)
    splitlevel[3, ] <- c(thresh, upper.split)
    splitlevel[4, ] <- c(upper.split, 1)
    splitlevely[1, ] <- c(0, 0.25)
    splitlevely[2, ] <- c(0.25, 0.5)
    splitlevely[3, ] <- c(0.5, 0.75)
    splitlevely[4, ] <- c(0.75, 1)
    pipos <- matrix(0, n, 4)
    ypos <- matrix(0, n, 4)
    for(k in 1:n){
      for(ksp in 1:3){
        if(splitlevel[ksp,1] <= pi[k] &
           pi[k] < splitlevel[ksp,2]) pipos[k,ksp] <- 1

        if(splitlevely[ksp,1] <= y[k] &
           y[k]  < splitlevely[ksp,2]) ypos[k,ksp] <- 1
      }
      if(splitlevel[4, 1] <= pi[k] &
         pi[k] <= splitlevel[4, 2]) pipos[k,4] <- 1

      if(splitlevely[4, 1] <= y[k] &
         y[k]  <= splitlevely[4, 2]) ypos[k,4] <- 1
    }
    cm <- matrix(0, 4, 4)
    for(k in 1:n){
      i <- which(ypos[k, ] == 1)
      j <- which(pipos[k, ] == 1)
      cm[i, j] <- cm[i, j] + 1
    }

    cm <- matrix(rev(as.vector(cm)), 4, 4, byrow = TRUE)
    # weights w
    w <- matrix(NA, 4, 4)
    for (i in 1:4){
      for (j in 1:4){
        w[i, j] <- ifelse(abs(i - j) < 2, 1, 0)
      }}
    sensitivity <- sum(w[ ,1:2] * cm[ ,1:2]) / sum(cm[ ,1:2])
    specificity <- sum(w[ ,3:4] * cm[ ,3:4]) / sum(cm[ ,3:4])

  }

  if(split == 2) {
    splitlevel[1, ] <- c(0, thresh)
    splitlevel[2, ] <- c(thresh, 1)
    pipos <- matrix(0, n, 2)
    ypos <- matrix(0, n, 2)
    for(k in 1:n){
      if(splitlevel[1, 1] <= pi[k] &
         pi[k] < splitlevel[1, 2]) pipos[k, 1] <- 1
      if(splitlevel[1, 1] <= y[k]  &
         y[k]  < splitlevel[1, 2]) ypos[k, 1] <- 1
      if(splitlevel[2, 1] <= pi[k] &
         pi[k] <= splitlevel[2, 2]) pipos[k, 2] <- 1
      if(splitlevel[2, 1] <= y[k]  &
         y[k]  <= splitlevel[2, 2]) ypos[k, 2] <- 1
    }
    cm <- matrix(0, 2, 2)
    for(k in 1:n){
      i <- which(ypos[k, ] == 1)
      j <- which(pipos[k, ] == 1)
      cm[i, j] <- cm[i, j] + 1
    }

    cm <- matrix(rev(as.vector(cm)), 2, 2, byrow = TRUE)
    w<-matrix(NA, 2, 2)
    for (i in 1:2){
      for (j in 1:2){
        # w = identity matrix
        if(split == 2) w[i, j] <- ifelse(abs(i - j) == 0, 1, 0)
      }}
    sensitivity <- sum(w[ ,1] * cm[ ,1]) / sum(cm[ ,1])
    specificity <- sum(w[ ,2] * cm[ ,2]) / sum(cm[ ,2])

  }

  n <- sum(cm)
  pobserved <- sum(w * cm) / n
  pexpected <- 0
  for (i in 1:split){
    for (j in 1:split){
      pexpected <- pexpected + w[i, j] * sum(cm[i, ]) %*% sum(cm[ ,j]) / (n^2)
    }}
  kappa <- (pobserved - pexpected) / (1 - pexpected)
  kappa <- as.numeric(kappa)

  return(list(kappa = kappa, cm = cm, sensitivity = sensitivity,
              specificity = specificity,
              actuals = y, splitlevel.pred = splitlevel,
              splitlevel.act = splitlevely, splitposition.pred = pipos,
              splitposition.act = ypos))
}
