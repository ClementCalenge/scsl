\name{scsl-package}
\alias{scsl-package}
\alias{scsl}
\docType{package}
\title{
  scsl -- The small carnivorous species logbook program
}
\description{
  This package is a convenient way to bundle the R code, the C code and
  the data used in the paper:\cr

  Calenge et al. (2014) The spatial distribution of Mustelidae
  species in France.
  
}
\details{
  \tabular{ll}{
    Package: \tab scsl\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2014-02-27\cr
    License: \tab GPL (>= 2)\cr
  }
  The model used in the paper is described on the help page of
  \code{penalizedmodel}.  The function \code{crossvalidation} can be
  used to identify the optimal value of the penalty.  The function
  \code{bootstrapcisstat} can be used to estimate the standard errors of
  the estimates with a bootstrap approach.\cr

  The data collected under the small carnivorous species logbook program
  are available in the following datasets: (i) the dataframe \code{scsl}
  contains the data related to the detected animals, (ii) the dataframe
  \code{carsscsl} contains the data related to the cars working in this
  program, (iii) the dataset \code{department} contains a map of the
  French department, and (iv) the dataset \code{SAR} contains a map of
  the small agricultural regions in the department.
}
\author{
  Clement Calenge (Office national de la chasse et de la faune sauvage)\cr
  Joel Chadoeuf (INRA Avignon, France)\cr
  Christophe Giraud (Universite Orsay, France)\cr
  Sylvie Huet (INRA, Jouy-en-Josas, France)\cr
  Romain Julliard (Museum national d'histoire naturelle)\cr
  Pascal Monestiez (INRA Avignon, France)\cr
  Jeremy Piffady (IRSTEA, France)\cr
  David Pinaud (CNRS Chize, France)\cr
  Sandrine Ruette (Office national de la chasse et de la faune sauvage)
}
\keyword{ package }
