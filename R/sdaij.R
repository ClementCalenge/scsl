sdaij <-
function(x, type=c("link", "response"))
{
    if (!inherits(x, "bootstrapCisstat"))
        stop("x should be of class bootstrapCisstat")
    type <- match.arg(type)
    sdd <- as.data.frame(x[[1]]$aij)
    for (i in 1:ncol(sdd)) {
        for (j in 1:nrow(sdd)) {
            if (type=="link")
                sdd[j,i] <- sd(unlist(lapply(x, function(y) y$aij[j,i])))
            if (type=="response")
                sdd[j,i] <- sd(unlist(lapply(x, function(y) exp(y$aij[j,i]))))
        }
    }
    return(sdd)
}
