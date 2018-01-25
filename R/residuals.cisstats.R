residuals.cisstats <-
function(object, type=c("pearson", "response"), ...)
{
    if (!inherits(object, "cisstats"))
        stop("object should be of class cisstats")
    type <- match.arg(type)
    res <- attr(object, "Y")-attr(object, "Yhat")
    if (type=="response")
        return(res)
    if (type=="pearson")
        return(res/sqrt(attr(object, "Yhat")))
}
