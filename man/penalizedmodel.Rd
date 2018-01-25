\name{penalizedmodel}
\alias{penalizedmodel}
\alias{print.cisstats}
\alias{predict.cisstats}
\alias{residuals.cisstats}
\alias{residualdf}
\alias{deviance.cisstats}
\title{
  Fit the regularized model described in the paper
}
\description{
  The function \code{penalizedmodel} fits the penalized model presented
  in the paper.  The functions \code{predict} and \code{residuals} can
  be used to extract the predictions and residuals of the model
  respectively.  The functions \code{deviance} and \code{residualdf} can
  be used to extract the deviance and the residual number of degrees of
  freedom of the model respectively.
}
\usage{
penalizedmodel(Sp, SU, Sta, Y, Vj, Area = rep(1, length(SU)), nu = 1,
eta = NULL, proximities = NULL, control = NULL)

\method{print}{cisstats}(x, \dots)

\method{predict}{cisstats}(object, type = c("link", "response"), \dots)

\method{residuals}{cisstats}(object, type = c("pearson", "response"), \dots)

\method{deviance}{cisstats}(object, \dots)

residualdf(object)
}
\arguments{
  \item{Sp}{
    a factor containing the species label for
    each combination of the factors Species/Spatial unit/Status.
  }
  \item{SU}{
    a factor containing the spatial unit label for
    each combination of the factors Species/Spatial unit/Status.
  }
  \item{Sta}{
    a factor containing the status label for
    each combination of the factors Species/Spatial unit/Status
  }
  \item{Y}{
    A numeric vector containing the number of detected animals
    for each combination of the factors Species/Spatial unit/Status
  }
  \item{Vj}{
    A numeric vector containing  the sampling effort of the
    dead animals for each combination of
    the factors Species/Spatial unit/Status (note that \code{Vj} is ignored if
    the status corresponds to the living animals)
  }
  \item{Area}{
    A vector containing the area of the spatial unit
    for each combination of the factors Species/Spatial unit/Status
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
    J x J matrix giving the proximity between the J spatial units
  }
  \item{control}{
    A list with named elements containing the options controlling the
    fit (see details).
  }
  \item{x,object}{
    An object of class \code{cisstats} returned by the function \code{penalizedmodel}.
  }
  \item{type}{
    For \code{predict.cisstats} a character value indicating the type of
    prediction required (either \code{"link"} which corresponds to the
    log-density, or \code{"response"} which corresponds to the actual
    density).  For \code{residuals.cisstats} a character value
    indicating the type of residuals required (either Pearson or
    response residuals).
  }
  \item{\dots}{
    Arguments passed to and from other methods.
  }
}
\details{
  The parameter \code{options} of the function \code{penalizedmodel}
  should be a list with named elements.  The following options are
  available:
  \describe{
    \item{maxIter}{maximum number of iterations allowed for the
      algorithm. Default to 10000}
    \item{h}{for the mixed and steepest descent approach, the
      step size. Default to 10}
    \item{stopCrit}{value of the threshold below which we judge
      that the algorithm has converged. Default to \code{10^-8}}
    \item{typealgo}{The type of algorithm used for the fit:
      value 1 = steepest descent; value 2 = mixed approach; value 3 =
      BFGS approach; value 4: approach consisting to use a steepest
      descent approach until a convergence evaluated with a threshold
      equal to \code{10^-6}, then a BFGS to reach convergence for the
      threshold defined with \code{stopCrit}. Default to 3.}
    \item{verbose}{logical value indicating if the evolution of
      the likelihood should be displayed in the console. Default to FALSE.}
    \item{StartingValues}{a list containing the starting values
      of the parameters, with three named elements: \code{aij} is a
      matrix or a data.frame (J spatial units x I species) corresponding
      to the parameters aij (abundances of species i in site j).
      \code{ej2} is a vector of length J containing the relative effort
      for living animals and \code{pi2} is a vector of length I containing
      the parameters pi2 for each species. If this element is not defined
      by the user, the program calculates the estimations with the
      non-penalized maximum likelihood  (which is very fast), and use them
      as starting values}
  }

}
\value{
The function \code{penalizedmodel} returns a list of class
\code{"cisstats"} containing the following elements:
\item{PenLik}{ The value of the penalized log-likelihood}
\item{aij}{ The estimated relative densities}
\item{ejk}{ The estimated relative effort for the living animals (as
  well as the known effort for the dead animals) }
\item{pi2}{ The estimated value of the parameters pi2 (complicated
  functions of the detection probability)}
\item{Penalty}{ The value of the penalty in the penalized likelihood at
  the solution}
\item{PenaltyMaxLik}{ The value of the penalty in the penalized likelihood at
  the maximum unpenalized likelihood}
\item{NiterConvergence}{ The number of iterations of the algorithm}
\item{h, eta, nu}{ The value of the parameters h, eta and nu}
\item{convergence}{ logical: whether the model has converged}
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
\code{\link{crossvalidation}}, \code{\link{bootstrapcisstat}}
}
\examples{
## See the appendix of the paper
## for examples

}
\keyword{models}
