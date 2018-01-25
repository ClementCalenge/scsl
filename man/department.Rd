\name{department}
\alias{department}
\docType{data}
\title{
  Map of the French departments
}
\description{
  This dataset is a map of the French departments (French administrative
  units)
}
\usage{data(department)}
\format{
  This is an object of class \code{'SpatialPolygonsDataFrame'} (package
  "sp"). For each department, the following variables are measured:
  \describe{
    \item{\code{SP_ID}}{The ID of the department}
    \item{\code{Name}}{A factor containing the name of the department}
    \item{\code{Code}}{A factor containing the official code of the
      department}
    \item{\code{Area}}{The area (squared kilometers) of each department}
  }
}
\source{
  Institut Geographique National, GeoFla (r)
  \url{http://professionnels.ign.fr/}.
}
\examples{
data(department)
plot(department)
}
\keyword{datasets}
