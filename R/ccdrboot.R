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
    
    
    ## No need to calculate, just use bnlearn::custom.strength()
    # probs <- lapply(which(lambdas %in% all_lambdas), function(l) {
    #     prob = matrix(0, ncol = ncol(data), nrow = ncol(data))
    #     for (pr in perRun) {
    #         tmp <- pr[[l]]
    #         tmpbn <- sparsebnUtils::to_igraph(tmp$edges)
    #         tmpadj <- igraph::get.adjacency(tmpbn, sparse=FALSE)
    #         for (i in seq_len(ncol(data))) {
    #             for (j in seq_len(ncol(data))) {
    #                 if (tmpadj[i,j] == 1) {
    #                     prob[i,j] <- prob[i,j] + 1
    #                 }
    #             }
    #         }
    #     }
    #     prob / R
    # })
    # names(probs) <- as.character(lambdas[lambdas %in% all_lambdas])
    # ## Strength and direction computing
    # res <- lapply(names(probs), function(p) {
    #     amat <- probs[[p]]
    #     ret <- NULL
    #     for (i in seq_len(ncol(data))) {
    #         for (j in seq_len(ncol(data))) {
    #             if (i==j) {
    #                 next
    #             }
    #             s <- amat[i, j] + amat[j, i]
    #             if (s==0) {d <- 0} else {
    #                 d <- amat[i,j] / s
    #             }
    #             if (d < 0) {d <- 0}
    #             if (d > 1) {d <- 1}
    #             if (s < 0) {s <- 0}
    #             if (s > 1) {s <- 1}
    #             ret <- rbind(ret, c(nodes[i], nodes[j], s, d))
    #         }
    #     }
    #     ret <- data.frame(ret) %>% `colnames<-`(c("from","to","strength","direction"))
    #     ret$strength <- as.numeric(ret$strength)
    #     ret$direction <- as.numeric(ret$direction)
    #     ret
    # })
    names(res) <- as.character(lambdas[lambdas %in% all_lambdas])
    res
}
