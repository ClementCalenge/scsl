summary.crossvalCisstatgl <-
function(object, ...)
{
    if (!inherits(object, "crossvalCisstatgl"))
        stop("object should be of class crossvalCisstatgl")
    df <- data.frame(object[[2]],
                     as.data.frame(do.call("cbind",lapply(object[[1]],
                                                          function(x) x[,ncol(x)]))))
    names(df) <- c("nu", paste("Q", 1:4, sep=""))
    return(df)
}
