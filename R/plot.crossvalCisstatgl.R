plot.crossvalCisstatgl <-
function(x, criterion=c(1:4), xlim=NULL,
                                   withgroups=TRUE, groups=1:(ncol(x[[1]][[1]])-1),...)
{
    if (!inherits(x, "crossvalCisstatgl"))
        stop("x should be of class crossvalCisstatgl")
    y <- summary(x)
    if (is.null(xlim))
        xlim <- range(y[,1])
    opar <- par(mfrow=n2mfrow(length(criterion)))

    y <- y[y[,1]>=xlim[1]&y[,1]<=xlim[2],]

    tmp <- lapply(criterion, function(crit) {
        dfcrit <- x[[1]][[crit]]
        dfcrit <- dfcrit[x[[2]]>=xlim[1]&x[[2]]<=xlim[2],]
        if (withgroups) {
            ylim <- range(unlist(dfcrit[,-ncol(dfcrit)]))
        } else {
            ylim <- range(y[,crit+1])
        }
        plot(y[,1], y[,crit+1], xlab="nu", ylab="Q",
             main=names(y)[crit+1], ylim=ylim, ty="l", lwd=2, col="blue",...)

        points(y[which.min(y[,crit+1]),1], y[which.min(y[,crit+1]),crit+1],
               pch=21, col="black",
               bg="yellow",cex=1.7)
        if (withgroups) {
            tmp <- lapply(groups, function(ngr) {
                val <- dfcrit[,ngr]
                lines(y[,1],val, col="grey")
                points(y[which.min(val),1], val[which.min(val)], pch=16, col="red")
            })
        }
    })
    return(invisible(NULL))
}
