sdpik <-
function(x)
{
    if (!inherits(x, "bootstrapCisstat"))
        stop("x should be of class bootstrapCisstat")
    sdd <- x[[1]]$pi2
    for (i in length(sdd)) {
        sdd[i] <- sd(unlist(lapply(x, function(y) y$pi2[i])))
    }
    names(sdd) <- names(x[[1]]$aij)
    return(sdd)
}
