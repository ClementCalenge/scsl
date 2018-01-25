print.crossvalCisstat <-
function(x, ...)
{
    if (!inherits(x, "crossvalCisstat"))
        stop("x should be of class crossvalCisstat")
    cat("######################################\n##\n## Object of class crossvalCisstat\n\n")
}
