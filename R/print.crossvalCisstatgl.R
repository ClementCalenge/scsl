print.crossvalCisstatgl <-
function(x, ...)
{
    if (!inherits(x, "crossvalCisstatgl"))
        stop("x should be of class crossvalCisstatgl")
    cat("######################################\n##\n## Object of class crossvalCisstatgl\n\n")
}
