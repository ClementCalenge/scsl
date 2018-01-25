penalizedmodel <-
function(Sp, SU, Sta, Y, Vj, Area=rep(1, length(SU)), nu=1, eta=NULL,
                           proximities=NULL, control=NULL)
{

############ Checks

    ## checks it is a factor
    if (!is.factor(Sp))
        stop("Sp should be a factor")
    if (!is.factor(SU))
        stop("SU should be a factor")
    if (!is.factor(Sta))
        stop("Sta should be a factor")

#    Area <- rep(1,length(SU))
    nrep <- 1

    ## check the length of the vectors
    le <- sapply(list(Sp, SU, Sta, Y, Vj, Area), length)
    if (length(unique(le))>1)
        stop("the vectors Sp, SU, Sta, Y, Vj and Area are not of equal length")

    ## check that all levels are represented in the factors
    if (nlevels(Sp)!=nlevels(factor(Sp)))
        stop("some levels of Sp are not present in the data")
    if (nlevels(SU)!=nlevels(factor(SU)))
        stop("some levels of SU are not present in the data")
    if (nlevels(Sta)!=nlevels(factor(Sta)))
        stop("some levels of Sta are not present in the data")

    ## check that Y >=0, nrep Vj and area >0
    if (any(Y<0))
        stop("Negative counts are present in Y")
    if (!all(Vj>0))
        stop("The effort is negative or equal to 0 for some spatial units")
    if (!all(Area>0))
        stop("The area is negative or equal to 0 for some spatial units")
    if (nrep<1)
        stop("The number of repetitions is negative or equal to 0")

    ## check that eta has a correct length
    if (!is.null(eta)&(length(eta)!=nlevels(Sp)))
        stop("the length of eta is not equal to the number of species")
    if (nu < 0)
        stop("negative values are not allowed for nu")


    ## check the proximity matrices
    if (!is.null(proximities)){
        if (nrow(proximities)!=ncol(proximities))
            stop("The number of rows and columns of proximities are not equal")
        if (nrow(proximities)!=nlevels(SU))
            stop("The dimension of the proximity matrix does not correspond\nto the number of spatial units")
    }

    ## If eta is not specified, gives a warning
    if (is.null(eta))
        warning("No eta specified: the program generates a vector of identical eta")
    if ((is.null(proximities))&(!is.null(eta)))
        stop("No specified proximity matrices, but lambda specified: why?")

    ## If proximities are not specified, gives a warning
    if (is.null(proximities))
        warning("No matrix proximities specified. An identity matrix be generated")

    ## If control not specified
    controlb <- list(maxIter=10000, h=10, stopCrit=1e-8, typealgo=3, verbose=FALSE,
                     StartingValues = list(0))
    if (is.null(control)) {
        control <- controlb
    } else {
        control <- lapply(c("maxIter", "h", "stopCrit", "typealgo",
                            "verbose", "StartingValues"),
                          function(x) {
                              if (length(control[x][[1]])==0)
                                  return(controlb[x][[1]])

                              if (x=="StartingValues") {
                                  con <- control[x][[1]]
                                  if (!all(names(con)%in%c("aij","ej2","pi2")))
                                      stop("StartingValues should contain elements\n aij, ej2 and pi2")
                                  con2 <- list(0, as.vector(t(as.matrix(con$aij))),
                                               con$ej2, con$pi2)
                                  return(con2)
                              }

                              return(control[x][[1]])
                          })
        names(control) <- names(controlb)
    }


############ Prepare the data

    I=length(levels(Sp))
    J=length(levels(SU))

    ## checks that the data have the correct order
    dford <- data.frame(Sp, SU, Sta, Y, Vj, Area)
    dford <- dford[order(Sp, Sta, SU),]
    Sp <- dford[,1]
    SU <- dford[,2]
    Sta <- dford[,3]
    Y <- dford[,4]
    Vj <- dford[,5]
    Area <- dford[,6]

    ## If the proximities are not provided, generate fake prox
    if (is.null(proximities))
        proximities <- diag(rep(1, J))

    ## If eta is not provided, values of eta are generated (spatial regularization only)
    if (is.null(eta))
        eta <- rep(1, I)

    ### calcul de nuijl
    nuijl <- array(0,c(I,J,J))
    for (i in 1:I) {
        nuijl[i,,] <- eta[i]*proximities
    }

    ## Estimation
    res <- .Call("ModeliseBien",as.integer(Sp), as.integer(SU), as.integer(Sta), as.double(Y),
                 as.double(Vj), as.integer(I), as.integer(J), as.integer(2),
                 as.double(Area), as.integer(nrep), as.integer(control$maxIter),
                 as.double(control$stopCrit), as.double(nuijl), as.double(nu),
                 as.double(control$h), as.integer(control$typealgo), as.integer(control$verbose),
                 control$StartingValues, PACKAGE = "scsl")

    names(res) <- c("PenLik", "aij", "ejk", "pi2", "Penalty", "PenaltyMaxLik",
                    "NiterConvergence", "h")
    res$aij <- matrix(res$aij, nrow=J, byrow=TRUE)

    if ((res$NiterConvergence)==control$maxIter) {
        warning("The algorithm did not converge. Try to increase maxIter")
        res$convergence <- FALSE
    } else {
        res$convergence <- TRUE
    }

    res$ejk <- data.frame(sapply(levels(SU), function(x) mean(log(Vj[SU==x]))),
                                 res$ejk)

    ## les noms de colonnes et de lignes
    row.names(res$aij) <- row.names(res$ejk) <- levels(SU)
    names(res$pi2) <- levels(Sp)
    names(res$ejk) <- levels(Sta)
    res$eta <- eta
    names(res$eta) <- levels(Sp)
    res$nu <- nu

    ## données originales
    attr(res, "X") <- data.frame(Sp, SU, Sta)
    attr(res, "Y") <- Y
    Yhat <- Y
    SU <- as.numeric(SU)
    Sp <- as.numeric(Sp)
    Sta <- as.numeric(Sta)
    Yhat <- sapply(1:length(Sp), function(i) {
        Area[i]*exp(res$aij[SU[i],Sp[i]] +
                    res$ejk[SU[i],Sta[i]] +
                    ifelse(Sta[i]==2, res$pi2[Sp[i]],0))
    })

    attr(res, "Yhat") <- Yhat

    ## Classe
    class(res) <- c("cisstats")
    return(res)
}
