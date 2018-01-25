print.bootstrapCisstat <-
function(x, ...)
{
    if (!inherits(x, "bootstrapCisstat"))
        stop("x should be of class bootstrapCisstat")
    cat("#######################################\n##\n## Object of class bootstrapCisstat\n\n")
    cat("This object is a list with one element per bootstrap sample (",length(x), " bootstrap samples)\n")
    cat("Each element is an object of class cisstats corresponding to one model\n")
    cat("Use sdaij(object) to calculate the standard deviation associated with\n each relative abundance\n")
    cat("Use sdejk(object) to calculate the standard deviation associated with\n each effort for status 2\n")
    cat("Use sdpik(object) to calculate the standard deviation associated with\n each detection probability for status 2\n")
}
