\name{crossvalidation}
\alias{crossvalidation}
\alias{global.crossvalCisstat}
\alias{print.crossvalCisstat}
\alias{print.crossvalCisstatgl}
\alias{summary.crossvalCisstatgl}
\alias{plot.crossvalCisstatgl}
\alias{plot.crossvalCisstat}
\title{
  Cross-validation approach for the choice of the penalty parameter
}
\description{
  The function \code{crossvalidation} carries out the cross-validation
  approach to select the penalty parameter, and stores the information
  in an object of class \code{crossvalCisstat}.  The function
  \code{global.crossvalCisstat} generates an object of class
  \code{crossvalCisstatgl}, that agglomerates the results over species.
}
\usage{
crossvalidation(Sp, SU, Sta, group, dfVj, vecnu, eta = NULL,
proximities = NULL, control = NULL, trace = TRUE)

global.crossvalCisstat(x)

\method{print}{crossvalCisstat}(x, \dots)
\method{print}{crossvalCisstatgl}(x, ...)
\method{summary}{crossvalCisstatgl}(object, ...)
\method{plot}{crossvalCisstatgl}(x, criterion = c(1:4), xlim = NULL,
                                 withgroups = TRUE,
                                 groups = 1:(ncol(x[[1]][[1]])-1), ...)

\method{plot}{crossvalCisstat}(x, criterion = 1, xlim = NULL,
                               withgroups = TRUE,
                               groups = 1:(ncol(x[[1]][[1]][[1]])-1), ...)

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
  \item{group}{
    a factor containing the group label for each detected animal (these
    groups are used for the K-fold cross-validation -- i.e. the years in
    the paper).
  }
  \item{dfVj}{
    a data.frame with three columns and J rows (where J is the number of
    spatial units) containing: (i) a factor containing the spatial unit
    label (with the same levels as \code{SU}, (ii) a numeric vector
    containing the sampling effort of the dead in these units, and (iii)
    a numeric vector containing the area of the SAR.
  }
  \item{vecnu}{
    a numeric vector containing the various values of nu to be tested.
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
  \item{x}{depending on the function, an object of class
    \code{"crossvalCisstat"} or \code{"crossvalCisstatgl"}.}
  \item{object}{an object of class \code{"crossvalCisstatgl"}.}
  \item{criterion}{for \code{plot.crossvalCisstatsgl} an integer vector,
    and for \code{plot.crossvalCisstats} an integer value. This argument
    indicates which measure of the prediction error should be examinated
    (see appendix of the paper).}
  \item{xlim}{a numeric vector of length 2 storing the limits of the
    x-axis on the graph}
  \item{withgroups}{logical value. Whether to also plot the results for the
    groups used in the K-fold cross-validation (i.e., the years). }
  \item{groups}{If \code{withgroups=TRUE}, the groups that should be
    displayed on the graphs}
  \item{\dots}{Arguments passed to and from other methods.}
}
\value{
  The function \code{crossvalidation} returns a list of class
  \code{"crossvalCisstat"}.  This list contains two elements: (i) a list
  with four elements containing the results for the four possible
  measures of prediction error (criteria), each element being itself a
  list with I sublists (where I is the number of species), each sublist
  being a data.frame containing the value of the criterion for each
  year (column) and each possible value of \code{nu} (row), and (ii) the
  vector containing the tested values of \code{nu}.\cr

  The function \code{global.crossvalCisstat} returns a list of class
  \code{"crossvalCisstatgl"}.  This list contains two elements: (i) a
  list with four data.frames containing the values of the criteria
  (summed over species) in each group (column) and for each value of
  \code{nu} (row), and (ii) the cetour containing the tested values of
  \code{nu}.  
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
\code{\link{penalizedmodel}}, \code{\link{bootstrapcisstat}}
}
\examples{
## See the appendix of the paper
## for examples

}
\keyword{models}
