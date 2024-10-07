#' ges.pcalg
#' Perform Greedy Equivalent Search implemented in pcalg
#' @param data data
#' @param args arguments to pcalg::ges
#' @export
#' @importFrom pcalg ges
ges.pcalg <- function(data, args=NULL) {
	if (is.null(args)) {
		args <- list()
	}
	score <- new(pcalg::.__C__GaussL0penObsScore, data=data)
	args[["score"]] <- score
	args[["iterate"]] <- TRUE
	ges.fit <- do.call(pcalg::ges, args)

	g <- ges.fit$repr$weight.mat()
	ig <- igraph::graph_from_adjacency_matrix(g, mode="directed",weighted = TRUE, diag=TRUE)
	bn <- bnlearn::as.bn(ig)
	return(bn)	
}

#' lingam.pcalg
#' Perform LiNGAM implemented in pcalg
#' @param data data
#' @param args arguments to pcalg::ges
#' @export
#' @importFrom pcalg ges
lingam.pcalg <- function(data, args=NULL) {
    if (is.null(args)) {
        args <- list()
    }
    args[["X"]] <- data
    fit <- do.call(pcalg::lingam, args)

    g <- fit$Bpruned

    row.names(g) <- colnames(data)
    colnames(g) <- colnames(data)

    ## Bpruned : a p \times p matrix B of linear coefficients,
    ## where B_{i,j} corresponds to a directed edge from j to i.
    ig <- igraph::graph_from_adjacency_matrix(g %>% t(),
        mode="directed",weighted = TRUE, diag=TRUE)
    bn <- bnlearn::as.bn(ig)
    return(bn)  
}


#' @title pcalg.boot
#' 
#' @description
#' Bootstrap-based arc strength calculation based on LiNGAM and GES.
#' 
#' @param data data (row: sample, column: gene)
#' @param R replicate number
#' @param m sampling number
#' @param debug passed to maximize function
#' @param return.all return all the network in bootstrapping
#' @param verbose control verbosity
#' @param removeAllZero remove all zero genes per replicate
#' @export
#' @return bn.strength 
pcalg.boot <- function(data, R=200, m=nrow(data), removeAllZero=FALSE,
    verbose=FALSE, return.all=FALSE, fun="ges", args=NULL) {
	if (fun=="ges") {
		fun <- ges.pcalg
	} else {
		fun <- lingam.pcalg
	}
    nodes = names(data)
    perRun <- list()
    for (r in seq_len(R)) {
        if (verbose) {
            cat_subtle("R: ", r, "\n")
        }
        resampling = sample(nrow(data), m, replace = TRUE)
        # generate the r-th bootstrap sample.
        replicate = data[resampling, , drop = FALSE]
        if (removeAllZero) {
	        replicate = replicate[, apply(replicate, 2, function(x) sum(x))!=0, drop=FALSE]
        }
        run <- fun(replicate, args)
        perRun[[r]] <- run
    }
    if (return.all) {
        return(perRun)
    }
	st <- custom.strength(perRun, nodes)
    return(st)
}
