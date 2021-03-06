% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/th_dep.R
\name{th.dep}
\alias{th.dep}
\title{Spatial threshold-dependent accuracy measures}
\usage{
th.dep(data, coord, thresh = 0.5, spatial = TRUE)
}
\arguments{
\item{data}{A data frame or matrix with two columns. The first column
should contain actual presence/absence data (binary, 0 or 1) and the
second column should contain model predictions of probability of
occurrence (numeric, between 0 and 1).}

\item{coord}{A data frame or matrix with two columns containing x,y
coordinates for each actual and predicted value. Coordinates must be
integer and consecutively numbered.}

\item{thresh}{A cutoff value for classifying predictions as modeled
presences or modeled absences. Default is 0.5.}

\item{spatial}{A logical indicating whether spatially corrected indices
(rather than classical indices) should be computed.}
}
\value{
A list with the following components:
\describe{
 \item{\code{kappa}}{Kappa statistic}
 \item{\code{cm}}{Confusion matrix}
 \item{\code{sensitivity}}{Sensitivity}
 \item{\code{specificity}}{Specificity}
 \item{\code{actuals}}{Actual occurrence data or adjusted actual occurrence data}
 \item{\code{splitlevel.pred}}{Level splitting of predicted values}
 \item{\code{splitlevel.act}}{Level splitting of actuals or adjusted actuals}
 \item{\code{splitposition.pred}}{Position splitting of predicted values}
 \item{\code{splitposition.act}}{Position splitting of actuals or adjusted actuals}
}
}
\description{
Calculates spatially corrected, threshold-dependent metrics for
an observational data set and model predictions (Kappa and confusion matrix)
}
\examples{
data(hook)
data <- hook[ ,1:2]
coord <- hook[ ,3:4]
si1 <- th.dep(data, coord, spatial = TRUE)
si1$kappa
si1$cm

}
\references{
Carl G, Kuehn I (2017) Spind: a package for computing
spatially corrected accuracy measures.
Ecography 40: 675-682. DOI: 10.1111/ecog.02593
}
\seealso{
\code{\link{th.indep}}
}
\author{
Gudrun Carl
}
