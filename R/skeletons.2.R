#' skeleton.reg
#' no non-exported functions
#' @noRd
skeleton.reg <- function(data, algorithm="glmnet_CV", whitelist=NULL, blacklist=NULL,
    nFolds=5, verbose=FALSE, s="lambda.min", maximize="hc", maximize.args=list(),
    returnArcs=FALSE, lambda=0.1) {
    if (verbose) {
        cat_subtle("Algorithm: ", algorithm, "\n")
        cat_subtle("Input for structure learning: n ", dim(data)[1], " p ", dim(data)[2], "\n")
    }
    nodes <- colnames(data)
    penFun <- dplyr::case_when(
        algorithm=="glmnet_BIC" ~ ".glmnetBIC",
        algorithm=="plasso" ~ ".plasso",
        algorithm=="glmnet_CV" ~ ".glmnetCV",
        algorithm=="MCP_CV" ~ ".ncvregCV",
        algorithm=="SCAD_CV" ~ ".ncvregCV",
        algorithm=="L0L2_CV" ~ ".L0LXCV",
        algorithm=="L0L1_CV" ~ ".L0LXCV",
        algorithm=="L0_CV" ~ ".L0CV"
    )

    mb <- sapply(nodes, function(nn) {
        if (verbose) {cat_subtle(" ", nn)}
        X <- data[, setdiff(nodes, nn)] %>% as.matrix()
        y <- data[, nn]
        if (length(unique(y))==1) {return(NULL)}
        if (algorithm=="MAST") {
            sign <- .MASTHurdle(nodes, nn, data)
            if (length(sign)==0) {return(NULL)}
            if (verbose) {cat_subtle(" candidate: ", length(sign), "\n")}
            return(sign)
        } else if (algorithm == "plasso") {
        	###
        	# For precision lasso, we need to specify a fixed lambda
        	# value in inference
        	###
            included_var_index <- do.call(penFun, list("X"=X, "y"=y,
                "nFolds"=nFolds, "algorithm"=algorithm, "s"=s, "lambda"=lambda))            
            if (verbose) {cat_subtle(" candidate: ", length(included_var_index), "\n")}
            if (length(included_var_index)==0) {return(NULL)}
            return(colnames(X)[included_var_index])
        }
            included_var_index <- do.call(penFun, list("X"=X, "y"=y,
                "nFolds"=nFolds, "algorithm"=algorithm, "s"=s))            
            if (verbose) {cat_subtle(" candidate: ", length(included_var_index), "\n")}
            if (length(included_var_index)==0) {return(NULL)}
            return(colnames(X)[included_var_index])
        }
    })

    names(mb) <- nodes
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

    if (returnArcs) {
        return(igraph::graph_from_edgelist(arcs))
    }

    ## No blacklisting in this version
    constraints <- arcs.to.be.added.2(arcs, nodes)
    start <- bnlearn::empty.graph(nodes = nodes)
    maximize.args[["x"]] <- data
    maximize.args[["start"]] <- start
    maximize.args[["blacklist"]] <- constraints
    struc_res <- do.call(maximize, maximize.args)
    # hc_res <- hc(data, start=start, blacklist=constraints)
    return(struc_res)
}


#' @noRd
skeleton.from.ig <- function(net, data, maximize="hc",
    maximize.args=list(), returnArcs=FALSE) {
    nodes <- names(V(net))

    mb <- lapply(names(V(net)), function(nn) {
       tmp <- names(igraph::neighborhood(net, 1, nn)[[1]])
       tmp[tmp!=nn]
    })
    names(mb) <- nodes

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

    if (returnArcs) {
        return(igraph::graph_from_edgelist(arcs))
    }

    ## No blacklisting in this version
    constraints <- arcs.to.be.added.2(arcs, nodes)
    start <- bnlearn::empty.graph(nodes = nodes)
    maximize.args[["x"]] <- data
    maximize.args[["start"]] <- start
    maximize.args[["blacklist"]] <- constraints
    struc_res <- do.call(maximize, maximize.args)
    # hc_res <- hc(data, start=start, blacklist=constraints)
    return(struc_res)
}