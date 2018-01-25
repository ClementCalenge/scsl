bootstrapcisstat <-
function(Sp, SU, Sta, dfVj, nu, eta=NULL,
                             proximities=NULL,
                             control=NULL, trace=TRUE, B=100)
{
    ## check that the vectors have the correct length
    le <- c(length(Sp),length(SU), length(Sta))
    if (!(all(le==le[1])))
        stop("vectors Sp, SU, Sta of different lengths")
    if (!is.factor(SU))
        stop("SU should be a factor")

    ## check that all levels are represented in the factors
    if (nlevels(Sp)!=nlevels(factor(Sp)))
        stop("some levels of Sp are not present in the data")
    if (nlevels(Sta)!=nlevels(factor(Sta)))
        stop("some levels of Sta are not present in the data")

    ## check that eta has a correct length
    if (!is.null(eta)&(length(eta)!=nlevels(Sp)))
        stop("the length of eta is not equal to the number of species")
    if (nu < 0)
        stop("negative values are not allowed for nu")

    ## checks that dfVj has the correct format
    if (!all(dfVj[,2]>0))
        stop("The effort is negative or equal to 0 for some spatial units")
    if (ncol(dfVj)!=3)
        stop("dfVj should have 3 columns")
    lev1 <- levels(SU)
    lev2 <- as.character(dfVj[,1])
    lev1 <- sort(lev1)
    lev2 <- sort(lev2)
    if (length(lev1)!=length(lev2))
        stop("the number of levels of SU does not correspond to the number of rows of dfVj")
    if (!all(lev1==lev2))
        stop("The levels of SU do not correspond to the number of levels of the first column of dfVj")
    dfVj <- dfVj[order(as.character(dfVj[,1])),]
    SU <- factor(SU, levels=lev2)
    OrderSU <- unique(SU)

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
                     StartingValues=list(0))
    if (is.null(control)) {
        control <- controlb
    } else {
        control <- lapply(c("maxIter", "h", "stopCrit", "typealgo",
                            "verbose", "StartingValues"), function(x) {
                                if (length(control[x][[1]])==0)
                                    return(controlb[x][[1]])

                                if (x=="StartingValues") {
                                    con <- control[x][[1]]
                                    if (!all(names(con)%in%c("aij","ej2","pi2")))
                                        stop("StartingValues should contain elements\n aij, ej2 and pi2")
                                    con2 <- list(0, as.vector(t(con$aij)), con$ej2, con$pi2)
                                    return(con2)
                                }

                                return(control[x][[1]])
                            })
        names(control) <- names(controlb)
    }


    ## prepares the data
    I=length(levels(Sp))
    J=length(levels(SU))

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

    resu <- list()
    for (b in 1:B) {
        if (b == 1) {
            sv <- control$StartingValues
        } else {
            sv <- list(0, as.vector(resu[[b-1]][[2]]), as.vector(resu[[b-1]][[3]]),
                       as.vector(resu[[b-1]][[4]]))
        }
        cat("Iteration", b, "\n")
        resu[[b]] <- .Call("BootstrapModele", as.integer(Sp), as.integer(SU), as.integer(Sta),
                           as.double(dfVj[,2]), as.integer(I), as.integer(J), as.integer(2),
                           as.integer(1), as.double(nu), as.integer(B), as.double(nuijl),
                           as.integer(control$maxIter),  as.double(dfVj[,3]),
                           as.double(control$stopCrit),
                           as.double(control$h), as.integer(control$typealgo),
                           as.integer(control$verbose), as.integer(trace),
                           control$StartingValues, PACKAGE = "scsl")
    }
    resu <- lapply(resu, function(res) {
        names(res) <- c("PenLik", "aij", "ejk", "pi2", "Penalty", "PenaltyMaxLik",
                        "NiterConvergence", "h")
        res$aij <- as.data.frame(matrix(res$aij, nrow=J, byrow=TRUE))
        if ((res$NiterConvergence)==control$maxIter) {
            warning("The algorithm did not converge. Try to increase maxIter")
            res$convergence <- FALSE
        } else {
            res$convergence <- TRUE
        }
        res$ejk <- as.data.frame(cbind(res$ejk, dfVj[,2]))
        ## les noms de colonnes et de lignes
        row.names(res$aij) <- row.names(res$ejk) <- levels(SU)
        names(res$aij) <- levels(Sp)
        names(res$pi2) <- levels(Sp)
        names(res$ejk) <- levels(Sta)
        res$eta <- eta
        names(res$eta) <- levels(Sp)
        res$nu <- nu
        ## Classe
        class(res) <- c("cisstats")
        return(res)
    })

    class(resu) <- "bootstrapCisstat"
    return(resu)
}
