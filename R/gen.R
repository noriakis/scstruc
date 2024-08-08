#' @title gen
#' @description generate random.graph based on gene node name and returns the 
#' fitted network based on the user-provided data
#' @param n sample number
#' @param data data to be fitted
#' @param method arguments in bnlearn::random.graph().
#' ic-dag, melancon or ordered. default to melancon.
#' @param params pass to random.graph
#' @export
#' @return bn object and bn.fit object
gen <- function(n, data, method="melancon", params=list()) {
    nodename <- sample(colnames(data), n, replace=FALSE)

    params[["method"]] <- method
    params[["nodes"]] <- nodename
    
    rdbn <- do.call(bnlearn::random.graph, params)
    fitted <- bnlearn::bn.fit(rdbn, data[, nodename])
    return(list("net"=rdbn, "fitted"=fitted))
}