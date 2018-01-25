deviance.cisstats <-
function(object, ...)
{
    if (!inherits(object, "cisstats"))
        stop("object should be of class cisstats")
    Y <- attr(object, "Y")
    tmp <- (Y-Y*log(Y))[Y>0]
    devmodcompl <- 2*sum(tmp)
    devcurmod <- 2*(object$PenLik - object$nu*object$Penalty)
    return(devcurmod-devmodcompl)
}
