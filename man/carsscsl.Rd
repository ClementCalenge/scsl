\name{carsscsl}
\alias{carsscsl}
\docType{data}
\title{
  Data frame with the information related to the cars
}
\description{
  This dataset is a dataframe containing the information related to the
  cars of the Office national de la chasse et de la faune sauvage.
}
\usage{data(carsscsl)}
\format{
  A data frame with 2979 observations on the following 4 variables.
  \describe{
    \item{\code{Car}}{a numeric vector containing the car label
      corresponding to each detection}
    \item{\code{Year}}{a numeric vector containing one year during which
      the car was used (note that a given car is often used several years).}
    \item{\code{Department}}{a character vector containing the
      department label where the car was used}
    \item{\code{Km}}{a numeric vector containing the number of
      kilometers traveled by the car during the year considered.}
  }
}
\source{
  Office national de la chasse et de la faune sauvage. 
  Lieu-dit Monfort. 01330 Birieux. France
}
\examples{
data(carsscsl)
head(carsscsl)

}
\keyword{datasets}
