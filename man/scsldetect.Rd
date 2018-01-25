\name{scsldetect}
\alias{scsldetect}
\docType{data}
\title{
  The small carnivorous species logbook data
}
\description{
  This data.frame contains the animals detected under the small
  carnivorous species logbook program of the "Office national de la
  chasse et de la faune sauvage" (French wildlife management organization).
}
\usage{data(scsldetect)}
\format{
  A data frame with 31811 observations on the following 7 variables.
  \describe{
    \item{\code{Department}}{a character vector containing the
      department label corresponding to each detection}
    \item{\code{Year}}{a numeric vector containing the year
      corresponding to each detection}
    \item{\code{Species}}{a character vector containing the species name 
      corresponding to each detection}
    \item{\code{Status}}{a character vector containing the status label
      corresponding to each detection}
    \item{\code{SAR}}{a character vector containing the small
      agricultural region label corresponding to each detection}
    \item{\code{Car}}{a numeric vector containing the car label
      corresponding to each detection}
  }
}
\source{
  Office national de la chasse et de la faune sauvage. 
  Lieu-dit Monfort. 01330 Birieux. France
}
\examples{
data(scsldetect)
head(scsldetect)

}
\keyword{datasets}
