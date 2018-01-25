plot.crossvalCisstat <-
function(x, criterion=1, xlim=NULL,
                                 withgroups=TRUE, groups=1:(ncol(x[[1]][[1]][[1]])-1),...)
{
    if (!inherits(x, "crossvalCisstat"))
        stop("x should be of class crossvalCisstat")

    if (is.null(xlim))
        xlim <- range(x[[2]])

    opar <- par(mfrow=n2mfrow(length(x[[1]][[criterion]])))
    y <- x[[1]][[criterion]]
    nu <- x[[2]]

    tmp <- lapply(1:length(y), function(i) {
        dfcrit <- y[[i]]
        if (withgroups) {
            ylim <- range(unlist(dfcrit[,-ncol(dfcrit)]))
        } else {
            ylim <- range(dfcrit[nu>=xlim[1]&nu<=xlim[2],ncol(dfcrit)])
        }
        plot(nu, dfcrit[,ncol(dfcrit)], xlab="nu", ylab="Criterion",
             main=names(y)[i], ylim=ylim, ty="l", lwd=2, col="blue", xlim=xlim,...)

        points(nu[which.min(dfcrit[,ncol(dfcrit)])],
               dfcrit[which.min(dfcrit[,ncol(dfcrit)]),ncol(dfcrit)],
               pch=21, col="black",
               bg="yellow",cex=1.7)

        if (withgroups) {
            tmp <- lapply(groups, function(ngr) {
                val <- dfcrit[,ngr]
                lines(nu,val, col="grey")
                points(nu[which.min(val)], val[which.min(val)], pch=16, col="red")
            })
        }
    })
    return(invisible(NULL))
}
