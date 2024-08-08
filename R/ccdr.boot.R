#' @title ccdr.boot
#' 
#' @description
#' Bootstrap-based arc strength calculation based on CCDr algorithm.
#' The function uses ccdrAlgorithm::ccdr.run for the sampled data, and returns the 
#' bootstrap strength for the arcs per lambda. The same lambdas are used for all the inference.
#' The networks with the same lambda are used to compute bootstapped network.
#' 
#' @param data data (row: sample, column: gene)
#' @param R replicate number
#' @param m sampling number
#' @param lambdas.length passed to generate.lambdas
#' @param alpha passed to ccdr.run, default to 2
#' @param gamma passed to ccdr.run, default to 2
#' @export
#' @return list of the bn.strength per lambda
#' 
ccdr.boot <- function(data, R=200, m=nrow(data), lambdas.length=20, alpha=2, gamma=2) {
    nodes = names(data)
    ## Use same lambda for all replicate (probably not needed)
    ## Set alpha to be sufficiently large
    lambdas <- sparsebnUtils::generate.lambdas(lambda.max=sqrt(nrow(data)), lambdas.ratio=1e-2,
        lambdas.length=lambdas.length, scale="log")
    perRun <- list()
    for (r in seq_len(R)) {
        resampling = sample(nrow(data), m, replace = TRUE)
        # generate the r-th bootstrap sample.
        replicate = data[resampling, , drop = FALSE]
        dat <- sparsebnUtils::sparsebnData(replicate, type = "continuous")
        run <- ccdrAlgorithm::ccdr.run(data = dat, lambdas=lambdas,
            alpha=alpha, gamma=gamma, verbose=FALSE)
        perRun[[r]] <- run
    }
    all_lambdas <- Reduce(intersect, lapply(perRun, function(pr) {unlist(lapply(pr, function(x) x$lambda))}))
    
    
    res <- lapply(which(lambdas %in% all_lambdas), function(l) {
    	net_list <- lapply(perRun, function(pr) {sparsebnUtils::to_bn(pr[[l]])$edges})
    	st <- custom.strength(net_list, nodes)
    	st
    })
    

    names(res) <- as.character(lambdas[lambdas %in% all_lambdas])
    res
}
