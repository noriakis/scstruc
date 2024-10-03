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
		suppressMessages(hurdle <- do.call(fitHurdle, hurdle.args))
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

	## Assign the direction
	arcs <- igraph::as_edgelist(g)
	arcs <- rbind(
	  arcs,
	  cbind(arcs[,2], arcs[,1])  
	)

	inc_node_undir <- V(g)$name

	# bl <- bnlearn:::arcs.to.be.added(arcs, V(g)$name)
	bl <- arcs.to.be.added.2(arcs, V(g)$name)

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

#' @noRd
#' @importFrom data.table :=
#' @importFrom bnlearn hc
#' @importFrom igraph V
#' @param cdrAdjustment if TRUE, perform cellular detection rate adjustment based on the count of zero per gene.
#' This will apply to graphical model inference in `fitHurdle` and scoring function in HC phase.
.Hurdle.2 <- function(data, score=NULL, debug=FALSE,
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
		suppressMessages(hurdle <- do.call(fitHurdle, hurdle.args))
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

	## Assign the direction
	nodes <- names(V(g))

	mb <- lapply(nodes, function(node) {
	    nn <- neighborhood(g, 1, nodes=node)[[1]] %>% names()
	    nn <- nn[nn!=node]
	    nn
	}) %>% setNames(nodes)

    ## Neighbors
    # mb <- lapply(nodes, function(nn) {
    #     candidate.neighbours <- mb[[nn]]
    #     if (length(candidate.neighbours) == 0) {
    #         return(list(mb = character(0), nbr = character(0)))
    #     }
    #     return(list(mb = mb[[nn]], nbr = candidate.neighbours))
    # })
    # names(mb) <- nodes
    # for (node in nodes)
    #   mb[[node]]$mb <- bnlearn:::fake.markov.blanket(mb, node)
    # mb <- bnlearn:::bn.recovery(mb, nodes = nodes)

    # arcs <- bnlearn:::nbr2arcs(mb)


    # res <- list(learning = list("method"="scstruc",
    #     "blacklist"=list(blacklist)),
    #     nodes = bnlearn:::cache.structure(names(mb), arcs = arcs),
    #     arcs = arcs, blacklist=blacklist)
    # constraints <- bnlearn:::arcs.to.be.added(res$arcs,
    #     nodes, whitelist = res$learning$blacklist)
    mb <- lapply(names(mb), function(node) {
        if (length(mb[[node]])!=0) {
            list(nbr=mb[[node]], mb=mb[[node]])
        } else {
            list(nbr=character(0), mb=character(0))
        }
    }) %>% setNames(nodes)

    for (node in nodes) {
        ## Corresponding to fake.markov.blanket
        dif <- setdiff(unique(lapply(mb[[node]]$nbr, function(current) {mb[[current]]$nbr}) %>% unlist(),
            mb[[node]]$nbr), node)
        mb[[node]]$mb <- dif
    }

    mb <- bn.recovery.2(mb)

    arcs <- do.call(rbind, lapply(names(mb), function(node) {
        candidate <- mb[[node]]$nbr
        cbind(rep(node, length(candidate)), candidate)
    }))
    arcs <- arcs[!duplicated(arcs),]
    constraints <- arcs.to.be.added.2(arcs, nodes)

    start <- empty.graph(nodes = nodes)

	if (is.null(score)) {
		net <- maximizeFun(data[,nodes],
			debug=debug, blacklist=constraints)
	} else {
		net <- maximizeFun(data[,nodes],
	         blacklist=constraints, debug=debug,
	         score="custom-score", fun=score,
	         args=list("cdrAdjustment"=cdrAdjustment))
	}
	retl[["bn"]] <- net
	retl[["data"]] <- data[, nodes]
	if (onlyBn) {
		return(net)
	} else {
		return(retl)
	}
}