#' @noRd
#' @importFrom data.table :=
.Hurdle <- function(data, score=NULL, debug=FALSE,
	skeleton=NULL, hurdle.args=list(), removeAllZero=FALSE,
	noMessages=TRUE) {

	# cat("Using score", as.character(score), "\n")

	# remove whole zero genes
	if (removeAllZero) {
		data <- data[, !apply(data==0, 2, function(x) sum(x)==nrow(data))]
	}

	retl <- list()
	if (is.null(skeleton)) {
		## Fit model
		hurdle.args[["samp"]] <- data
		suppressMessages(hurdle <- do.call(fitHurdle, hurdle.args))
	} else {
		hurdle <- skeleton
	}

	## Best model selection
	BIC_etc <- hurdle$BIC_etc
	BIC_etc[BIC<=min(BIC),isMin := TRUE][,idx:=.I]

	ind <- which.min(BIC_etc$BIC)
	g <- graph.adjacency(abs(hurdle$adjMat[[ind]])>0, mode='max')
	retl[["hurdle"]] <- hurdle
	retl[["undir"]] <- g


	## Assign the direction ...
	arcs <- igraph::as_edgelist(g)
	arcs <- rbind(
	  arcs,
	  cbind(arcs[,2], arcs[,1])  
	)

	inc_node_undir <- V(g)$name

	bl <- bnlearn:::arcs.to.be.added(arcs, V(g)$name)

	if (is.null(score)) {
		net <- hc(data[,inc_node_undir],
			debug=debug, blacklist=bl)
	} else {
		net <- hc(data[,inc_node_undir],
	         blacklist=bl, debug=debug,
	         score="custom-score", fun=score)		
	}
	retl[["bn"]] <- net
	retl[["data"]] <- data[, inc_node_undir]
	return(retl)
}