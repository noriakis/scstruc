#' interVal
#' @param ss sample size
interVal <- function(data, algos=c("hc","mmhc"), algorithm.args=list(list(),list()),
    mode="shd", r=10, ss=100, verbose=FALSE, returnA0=FALSE) {
    if (length(algos)!=length(algorithm.args)) {
        stop("length of algos and algorithm.args must match")
    }
	bns <- lapply(seq_along(algos), function(al) {
		.getStruc(data, algos[al], algorithm.args=algorithm.args[[al]], verbose=verbose)
	})
	A0 <- Reduce(intersection, lapply(bns, function(x) {
    	as.igraph(x)
	}))
	A0.bn <- as.bn(A0)
	subsample.res <- do.call(cbind, lapply(seq_len(r), function(cr) {
        replicates <- sample(seq_len(nrow(data)), ss)
        tmp <- data[replicates, ]
        inf.nets <- lapply(seq_along(algos), function(al) {
            .getStruc(data, algos[al], algorithm.args=algorithm.args[[al]], verbose=verbose)
        })
        lapply(inf.nets, function(tmpnet) {
            if (mode=="sid") {
    	          sid.bn(A0.bn, tmpnet) %>% unlist() %>% return()
            } else if (mode=="shd") {
                bnlearn::shd(A0.bn, tmpnet) %>% unlist() %>% return()
            } else {
                stop("sid or shd should be specified")
            }
        })
    }))
    sum.stat <- apply(subsample.res, 1, function(x) mean(as.numeric(x)))
    if (returnA0) {
        return(list("A0"=A0.bn, "stat"=sum.stat))
    } else {
        return(sum.stat)
    }
}


#' sid.bn
#' @noRd
sid.bn <- function(true.bn, est.bn) {
  res <- SID::structIntervDist(as_adjacency_matrix(as.igraph(true.bn)),
                               as_adjacency_matrix(as.igraph(est.bn)))
  res$sid
}