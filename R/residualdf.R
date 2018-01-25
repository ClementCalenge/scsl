residualdf <-
function(object)
{
    if (!inherits(object, "cisstats"))
        stop("object should be of class cisstats")
    I <- ncol(object$aij)
    J <- nrow(object$aij)
    K <- 2
    return((I*J*K)-(I*J)-J-(I-1))
}
