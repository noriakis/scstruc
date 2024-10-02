#' @title interVal
#' @description Performs `Intersection-Validation` approach to the data.
#' @param data data (row names as genes and columns as samples)
#' @param algos algorithms to be evaluated.
#' @param algorithm.args list of lists to be passed to algos
#' @param mode shd or sid
#' @param r number of iteration
#' @param ss sample size
#' @param verbose control verbosity
#' @param returnA0 returns the A0 (intersection of inferred networks), default to TRUE
#' @param symmetrize use symmetrized version of SID
#' @return summarized statistics
#' @export
#' @examples
#' data(gaussian.test)
#' interVal(head(gaussian.test, n=50), ss=10)
interVal <- function(data, algos=c("hc","mmhc"), algorithm.args=list(list(),list()),
    mode="shd", r=10, ss=100, verbose=FALSE, returnA0=TRUE, symmetrize=FALSE) {
    if (length(algos)!=length(algorithm.args)) {
        stop("length of algos and algorithm.args must match")
    }
    if (verbose) {
        cat_subtle("Total of ", length(algos), " algorithms\n")
    }
	bns <- lapply(seq_along(algos), function(al) {
        if (verbose) {
            cat_subtle("  Algorithm: ", algos[al], "\n")
        }
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
            .getStruc(tmp, algos[al], algorithm.args=algorithm.args[[al]], verbose=verbose)
        })
        lapply(inf.nets, function(tmpnet) {
            if (mode=="sid") {
                if (symmetrize) {
                    bnlearn.sid.sym(tmpnet, A0.bn) %>% return()
                } else {
                    bnlearn::sid(tmpnet, A0.bn) %>% unlist() %>% return()
                }
            } else if (mode=="shd") {
                bnlearn::shd(tmpnet, A0.bn) %>% unlist() %>% return()
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

