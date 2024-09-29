#' @title hurdle.boot
#' 
#' @description
#' Bootstrap-based arc strength calculation based on Hurdle model.
#' 
#' @param data data (row: sample, column: gene)
#' @param R replicate number
#' @param m sampling number
#' @param score scoring function after identification of skeleton.
#' @param removeAllZero remove all zero genes in replicates
#' @param debug passed to maximize function
#' @param skeleton if skeleton is provided, the bootstrap is based on
#' these skeleton. This corresponds to hurdle fitted object.
#' @param maximizeFun maximization function (default to HC)
#' @param verbose logging in scstruc, for bnlearn debugging, change the 
#' argument `debug`
#' @param return.all return all the network in bootstrapping
#' @export
#' @return bn.strength 
hurdle.boot <- function(data, R=200, m=nrow(data), score=NULL, removeAllZero=FALSE,
    debug=FALSE, skeleton=NULL, maximizeFun=hc, verbose=FALSE, return.all=FALSE) {
    nodes = names(data)
    perRun <- list()
    for (r in seq_len(R)) {
        if (verbose) {
            cat_subtle("R: ", r, "\n")
        }
        resampling = sample(nrow(data), m, replace = TRUE)
        # generate the r-th bootstrap sample.
        replicate = data[resampling, , drop = FALSE]
        run <- .Hurdle(replicate, score=score,
            removeAllZero=removeAllZero, debug=debug, skeleton=skeleton,
            maximizeFun=maximizeFun)
        perRun[[r]] <- run
    }
    if (return.all) {
        return(perRun)
    }
	st <- custom.strength(perRun, nodes)
    return(st)
}
