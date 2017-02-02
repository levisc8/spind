#' Mus musculus data set
#'
#' A data frame containing presence absence records of Mus musculus (a house mouse)
#'
#'@format A data frame with 400 rows and 5 columns
#'  \describe{
#'    \item{musculus}{integer - Presence/absence records of Mus musculus}
#'    \item{pollution}{numberic - describes degree of pollution in corresponding grid cell}
#'    \item{exposure}{numeric - describes degree of exposure for each grid cell}
#'    \item{x}{integer - x-coordinates for each grid cell}
#'    \item{y}{numeric - y-coordinates for each grid cell}
#'  }
#'
"musdata"

#'hook data set
#'
#'A data frame containing actual presence absence data and predicted probability
#'of occurrence
#'
#'@format A data frame with 100 rows and 4 columns
#'  \describe{
#'    \item{actuals}{integer - Presence/absence records}
#'    \item{predictions}{numeric - predicted probabilities of occurence}
#'    \item{x}{integer - x-coordinates for each grid cell}
#'    \item{y}{numeric - y-coordinates for each grid cell}
#'
#'  }
"hook"
