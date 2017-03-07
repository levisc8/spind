#' Mus musculus data set
#'
#' A data frame containing simulated count data of a house mouse.
#'
#'@format A data frame with 400 rows and 5 columns
#'  \describe{
#'    \item{musculus}{integer - Simulated count data for Mus musculus}
#'    \item{pollution}{numeric - Simulated variable that describes degree
#'     of pollution in corresponding grid cell}
#'    \item{exposure}{numeric - Simulated variable that describes degree of
#'     exposure for each grid cell}
#'    \item{x}{integer - x-coordinates for each grid cell}
#'    \item{y}{integer - y-coordinates for each grid cell}
#'  }
#'
"musdata"

#'Hook data set
#'
#'A data frame containing actual presence absence data and predicted probability
#'of occurrence values.
#'
#'@format A data frame with 100 rows and 4 columns
#'  \describe{
#'    \item{actuals}{integer - Presence/absence records}
#'    \item{predictions}{numeric - predicted probabilities of occurence}
#'    \item{x}{integer - x-coordinates for each grid cell}
#'    \item{y}{integer - y-coordinates for each grid cell}
#'
#'  }
"hook"

#' Carlina data set
#'
#' A data frame containing simulated count data for the thistle, Carlina horrida.
#'
#' @format A data frame with 961 rows and 5 columns
#' \describe{
#'   \item{carlina.horrida}{integer - Simulated count data}
#'   \item{aridity}{numeric - Simulated aridity index values. This variable has high spatial autocorrelation
#'   values.}
#'   \item{land.use}{numeric - Simulated land use intensity. This variable has no spatial autocorrelation.}
#'   \item{x}{integer - x-coordinates for each grid cell}
#'   \item{y}{integer - y-coordinates for each grid cell}
#' }
"carlinadata"
