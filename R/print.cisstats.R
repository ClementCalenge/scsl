print.cisstats <-
function(x, ...)
{
    if (!inherits(x, "cisstats"))
        stop("x should be of class cisstats")
    cat("########################################\n##\n## ")
    cat("Penalized model of abundance for", length(x$pi2) ,"species in",
        nrow(x$aij), "spatial units:\n\n")
    cat("Penalty coefficient (nu):", x$nu, "\n")
    cat("Penalized likelihood for the solution:", x$PenLik, "\n")
    cat("Penalty:", x$Penalty,"\n")
    cat("Penalized likelihood for the maximum likelihood solution:", x$PenaltyMaxLik, "\n\n")
}
