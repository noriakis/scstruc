#' @title interVal
#' @description Performs `Intersection-Validation` approach to the data.
#' @param data data (row names as genes and columns as samples)
#' @param algos algorithms to be evaluated.
#' @param algorithm.args list of lists to be passed to algos
#' @param r number of iteration
#' @param ss sample size
#' @param verbose control verbosity
#' @param returnA0 returns the A0 (intersection of inferred networks), default to TRUE
#' @param symmetrize use symmetrized version of SID
#' @param SID.cran use SID in CRAN package SID
#' @param output output relevant data to the directory
#' @return summarized statistics, bn object of A0, and statistics for all the subsamples
#' @export
#' @importFrom igraph intersection
#' @importFrom bnlearn as.bn
#' @examples
#' data(gaussian.test)
#' interVal(head(gaussian.test, n=50), ss=10)
interVal <- function(data, algos=c("hc","mmhc"), algorithm.args=NULL,
    r=10, ss=100, verbose=FALSE, returnA0=TRUE, symmetrize=FALSE, SID.cran=FALSE,
    output=NULL) {
    if (!is.null(output)) {
        write.table(data, file=paste0(output, "/raw-data.txt"), sep="\t")
    }
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
	A0 <- Reduce(igraph::intersection, lapply(bns, function(x) {
        if (verbose) {
            cat("Inferred network:\n")
            print(x)
        }
    	as.igraph(x)
	}))
	A0.bn <- bnlearn::as.bn(A0)
    if (!is.null(output)) {
        save(A0.bn, file=paste0(output, "/A0.rda"))
        save(bns, file=paste0(output, "/all-bns.rda"))
    }
    if (verbose) {
        cat("Inferred A0:\n")
        print(A0.bn)
    }
    CNP <- dim(A0.bn$arcs)[1]
    if (CNP < 15) {
        message("Connected node pairs below 15, the results might not be reliable.")
    }
	subsample.res <- lapply(seq_len(r), function(cr) {
        if (verbose) {cat_subtle("Subsample: ", cr, ", Sample size: ", ss, "\n")}
        replicates <- sample(seq_len(nrow(data)), ss)
        tmp <- data[replicates, ]

        inf.nets <- lapply(seq_along(algos), function(al) {
            .getStruc(tmp, algos[al], algorithm.args=algorithm.args[[al]], verbose=verbose)
        })

        if (!is.null(output)) {
            write.table(tmp, file=paste0(output, "/subsample-data-",cr,".txt"), sep="\t")
            save(tmp, file=paste0(output, "/subsample-inf-net-",cr,".rda"))
        }
        lapply(seq_along(inf.nets), function(tmpnetnum) {
            tmpnet <- inf.nets[[tmpnetnum]]
            if (SID.cran) {
                val <- SID.sid(A0.bn, tmpnet, sym=symmetrize)
            } else {
                if (symmetrize) {
                    val <- bnlearn.sid.sym(tmpnet, A0.bn)
                } else {
                    val <- bnlearn::sid(tmpnet, A0.bn) %>% unlist()
                }                    
            }

            shd.val <- bnlearn::shd(tmpnet, A0.bn) %>% unlist()

            return(list("algonum"=tmpnetnum,
                "edgenum"=dim(tmpnet$arcs)[1], "shd"=shd.val, "sid"=val))
        })
    })

    long.res <- do.call(rbind, lapply(seq_along(subsample.res), function(x) {
        do.call(rbind, lapply(subsample.res[[x]], function(xx) {
            c(x, xx$algonum, xx$shd, xx$sid, xx$edgenum)
        }))
    })) %>% data.frame() %>% `colnames<-`(c("R","AlgoNum","SHD","SID","EdgeNumber"))
    if (!is.null(output)) {
        write.table(long.res, file=paste0(output, "/raw-results.txt"), sep="\t")
    }
    sum.stat <- long.res %>% group_by(.data$AlgoNum) %>% 
        summarise(SHD.stat=mean(.data$SHD), SID.stat=mean(.data$SID),
        en=mean(.data$EdgeNumber))
    # sum.stat <- apply(subsample.res, 1, function(x) mean(as.numeric(x)))
    if (returnA0) {
        return(list("A0"=A0.bn, "stat"=sum.stat, "raw.stat"=long.res))
    } else {
        return(list("stat"=sum.stat, "raw.stat"=long.res))
    }
}

