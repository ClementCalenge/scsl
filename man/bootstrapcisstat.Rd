\name{bootstrapcisstat}
\alias{bootstrapcisstat}
\alias{print.bootstrapCisstat}
\alias{sdaij}
\alias{sdejk}
\alias{sdpik}
\title{
  Calculate the standard error of abundances estimates with a bootstrap approach.
}
\description{
  The function \code{bootstrapcisstat} calculates the standard error of
  abundances estimates with a bootstrap approach.  The functions
  \code{sdaij}, \code{sdeij}, and \code{sdpik} allow to extract the
  standard errors of the relative densities, relative effort and of the
  parameters pik from the resulting object.
}
\usage{
bootstrapcisstat(Sp, SU, Sta, dfVj, nu, eta = NULL, proximities = NULL,
control = NULL, trace = TRUE, B = 100)

\method{print}{bootstrapCisstat}(x, \dots)

sdaij(x, type=c("link", "response"))
sdejk(x, type=c("link", "response"))
sdpik(x)

}
\arguments{
  \item{Sp}{
    a factor containing the species label for each detected animal.
  }
  \item{SU}{
    a factor containing the spatial unit label for each detected animal.    
  }
  \item{Sta}{
    a factor containing the status label for each detected animal.    
  }
  \item{dfVj}{
    a data.frame with three columns and J rows (where J is the number of
    spatial units) containing: (i) a factor containing the spatial unit
    label (with the same levels as \code{SU}, (ii) a numeric vector
    containing the sampling effort of the dead in these units, and (iii)
    a numeric vector containing the area of the SAR.
  }
  \item{nu}{
    value of the penalty parameter used in the regularization.
  }
  \item{eta}{
    numeric vector of length I (where I is the number of species)
    allowing to introduce differing penalties for the different
    species.  Thus, the difference of relative density between spatial
    units j and m for the species i with:
    \code{nu*eta[i]*proximities[j,l]} (this parameter was introduced for
    exploratory purposes, but is not actually used in the paper)
  }
  \item{proximities}{
    J x J matrix giving the proximity between the J spatial units.
  }
  \item{control}{
    a list with named elements controlling the fit (see
    \code{?penalizedmodel} for a description of the available options. 
  }
  \item{trace}{
    Whether information about the progression of the crossvalidation
    should be printed.
  }
  \item{B}{
    The number of bootstrap samples
  }
  \item{x}{
    An object of class \code{"bootstrapCisstat"}    
  }
  \item{type}{
    Whether the resulting standard errors should be expressed on the
    log-scale (\code{type="link"}) or on the response scale
    (\code{type="response"}). 
  }
  \item{\dots}{
    Arguments passed to and from other methods. 
  }
}
\value{
  A list of class \code{"bootstrapCisstat"} containing one fitted model
  (class \code{"cisstats"}, see \code{penalizedmodel}) per bootstrap
  sample.
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
\seealso{
  \code{\link{crossvalidation}}, \code{\link{penalizedmodel}}
}
\examples{
## See the appendix of the paper
## for examples

}
\keyword{models}
