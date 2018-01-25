sdejk <-
function(x, type=c("link", "response"))
{
    if (!inherits(x, "bootstrapCisstat"))
        stop("x should be of class bootstrapCisstat")
    type <- match.arg(type)
    sdd <- x[[1]]$ejk
    for (j in 1:length(sdd)) {
        if (type=="link")
            sdd[j] <- sd(unlist(lapply(x, function(y) y$ejk[j])))
        if (type=="response")
            sdd[j] <- sd(unlist(lapply(x, function(y) exp(y$ejk[j]))))
    }
    return(sdd)
}
