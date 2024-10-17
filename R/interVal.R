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
#' @param SID.cran use SID in CRAN package SID
#' @return summarized statistics, bn object of A0, and statistics for all the subsamples
#' @export
#' @examples
#' data(gaussian.test)
#' interVal(head(gaussian.test, n=50), ss=10)
interVal <- function(data, algos=c("hc","mmhc"), algorithm.args=NULL,
    mode="shd", r=10, ss=100, verbose=FALSE, returnA0=TRUE, symmetrize=FALSE, SID.cran=FALSE) {
    if (is.null(algorithm.args)) {
        algorithm.args <- lapply(seq_len(length(algos)), function(x) return(NULL))
    } else {
        if (length(algos)!=length(algorithm.args)) {
            stop("length of algos and algorithm.args must match")
        }        
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
        if (verbose) {
            cat("Inferred network:\n")
            print(x)
        }
    	as.igraph(x)
	}))
	A0.bn <- as.bn(A0)
    if (verbose) {
        cat("Inferred A0:\n")
        print(A0.bn)
    }
	subsample.res <- do.call(cbind, lapply(seq_len(r), function(cr) {
        replicates <- sample(seq_len(nrow(data)), ss)
        tmp <- data[replicates, ]
        inf.nets <- lapply(seq_along(algos), function(al) {
            .getStruc(tmp, algos[al], algorithm.args=algorithm.args[[al]], verbose=verbose)
        })
        lapply(inf.nets, function(tmpnet) {
            if (mode=="sid") {
                if (SID.cran) {
                    SID.sid(A0.bn, tmpnet, sym=symmetrize)
                } else {
                    if (symmetrize) {
                        bnlearn.sid.sym(tmpnet, A0.bn) %>% return()
                    } else {
                        bnlearn::sid(tmpnet, A0.bn) %>% unlist() %>% return()
                    }                    
                }
            } else if (mode=="shd") {
                bnlearn::shd(tmpnet, A0.bn) %>% unlist() %>% return()
            } else {
                stop("sid or shd should be specified.")
            }
        })
    }))
    sum.stat <- apply(subsample.res, 1, function(x) mean(as.numeric(x)))
    if (returnA0) {
        return(list("A0"=A0.bn, "stat"=sum.stat, "raw.stat"=subsample.res))
    } else {
        return(list("stat"=sum.stat, "raw.stat"=subsample.res))
    }
}

