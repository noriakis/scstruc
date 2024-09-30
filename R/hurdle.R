#' @noRd
#' @importFrom data.table :=
#' @importFrom bnlearn hc
#' @importFrom igraph V
#' @param cdrAdjustment if TRUE, perform cellular detection rate adjustment based on the count of zero per gene.
#' This will apply to graphical model inference in `fitHurdle` and scoring function in HC phase.
.Hurdle <- function(data, score=NULL, debug=FALSE,
	skeleton=NULL, hurdle.args=list(), removeAllZero=FALSE,
	noMessages=TRUE, cdrAdjustment=FALSE, maximizeFun=bnlearn::hc,
	onlyBn=TRUE) {

	# cat("Using score", as.character(score), "\n")

	# remove whole zero genes
	if (removeAllZero) {
		data <- data[, !apply(data==0, 2, function(x) sum(x)==nrow(data))]
	}
	if (cdrAdjustment) {
		cdr <- rowSums(data>0)
		covariates <- cbind(1, scale(cdr))		
	}

	retl <- list()
	if (is.null(skeleton)) {
		## Fit model
		hurdle.args[["samp"]] <- data
		if (is.null(hurdle.args[["fixed"]]) & cdrAdjustment) {
			hurdle.args[["fixed"]] <- covariates
		}
		suppressMessages(hurdle <- do.call(HurdleNormal::fitHurdle, hurdle.args))
	} else {
		hurdle <- skeleton
	}

	## Best model selection
	BIC_etc <- hurdle$BIC_etc
	# BIC_etc[BIC<=min(BIC),isMin := TRUE][,idx:=.I]

	ind <- which.min(BIC_etc$BIC)
	g <- igraph::graph.adjacency(abs(hurdle$adjMat[[ind]])>0, mode='max')
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
		net <- maximizeFun(data[,inc_node_undir],
			debug=debug, blacklist=bl)
	} else {
		net <- maximizeFun(data[,inc_node_undir],
	         blacklist=bl, debug=debug,
	         score="custom-score", fun=score,
	         args=list("cdrAdjustment"=cdrAdjustment))
	}
	retl[["bn"]] <- net
	retl[["data"]] <- data[, inc_node_undir]
	if (onlyBn) {
		return(net)
	} else {
		return(retl)
	}
}