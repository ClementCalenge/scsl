global.crossvalCisstat <-
function(x)
{
    y <- x[[1]]
    y <- lapply(y, function(listecrit) {
        dfesp <- as.matrix(listecrit[[1]])
        dfesp <- dfesp[,-ncol(dfesp)]
        for (i in 2:length(listecrit))
            dfesp <- dfesp+listecrit[[i]][,-(ncol(dfesp)+1)]
        dfesp <- dfesp
        dfesp <- data.frame(as.data.frame(dfesp), mean=apply(dfesp,1,mean))
    })
    res <- list(y, x[[2]])
    class(res) <- "crossvalCisstatgl"
    return(res)
}
