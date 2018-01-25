predict.cisstats <-
function(object, type=c("link","response"),...)
{
    if (!inherits(object, "cisstats"))
        stop("object should be of class cisstats")
    type <- match.arg(type)
    if (type=="link")
        return(log(attr(object, "Yhat")))
    if (type=="response")
        return(attr(object, "Yhat"))
}
