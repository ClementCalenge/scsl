\name{SAR}
\alias{SAR}
\docType{data}
\title{
  Map of the small agricultural regions
}
\description{
  This dataset is a map of the French small agricultural region (SAR)
}
\usage{data(SAR)}
\format{
  This is an object of class \code{'SpatialPolygonsDataFrame'} (package
  "sp"). For each SAR, the following variables are measured:
  \describe{
    \item{\code{SP_ID}}{The ID of the SAR}
    \item{\code{codeSAR}}{A factor containing the official code of the
      SAR (the first two digits correspond to the code of the
      department)}
    \item{\code{Name}}{The name of the SAR}
    \item{\code{Area}}{A numeric vector containing the area (squared
      kilometers) of each SAR}
    \item{\code{Department}}{A factor containing the official code of the
      department where the SAR is found}}
}
\source{
  Institut Geographique National, GeoFla (r)
  (\url{http://professionnels.ign.fr/}).\cr
  Ministere de l'agriculture, de l'agroalimentaire et de la foret
  (\url{http://www.agreste.agriculture.gouv.fr/definitions/zonages/})
}
\examples{
data(SAR)
plot(SAR)
}
\keyword{datasets}
