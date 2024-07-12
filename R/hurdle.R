#' @noRd
#' @importFrom data.table :=
.Hurdle <- function(data, score=NULL, debug=FALSE) {
	retl <- list()

	## Fit model
	hurdle <- fitHurdle(data)

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
		net <- hc(mat[,inc_node_undir],
			debug=debug, blacklist=bl)
	} else {
		net <- hc(mat[,inc_node_undir],
	         blacklist=bl, debug = debug,
	         score="custom", fun=score)		
	}
	retl[["bn"]] <- net
	return(retl)
}