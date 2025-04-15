#' ges.pcalg
#' Perform Greedy Equivalent Search implemented in pcalg
#' @param data data
#' @param sf score name
#' @param args arguments to pcalg::ges
#' @export
#' @importFrom pcalg ges
ges.pcalg <- function(data, sf="GaussL0penObsScore", args=NULL) {
	if (is.null(args)) {
		args <- list()
	}
	score <- new(sf, data=data)
	args[["score"]] <- score
	# args[["iterate"]] <- TRUE
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
#' @param return.all return all the network in bootstrapping
#' @param verbose control verbosity
#' @param removeAllZero remove all zero genes per replicate
#' @param fun pcalg functions, lingam or ges
#' @param args args to be passed to the main function of pcalg
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
        args[["data"]] <- replicate
        run <- do.call(fun, args)
        perRun[[r]] <- run
    }
    if (return.all) {
        return(perRun)
    }
	st <- custom.strength(perRun, nodes)
    return(st)
}


#' @title pcalg.boot.future
#' 
#' @description
#' Bootstrap-based arc strength calculation based on LiNGAM and GES.
#' The future.lapply will be used in this function.
#' 
#' @param data data (row: sample, column: gene)
#' @param R replicate number
#' @param m sampling number
#' @param return.all return all the network in bootstrapping
#' @param verbose control verbosity
#' @param removeAllZero remove all zero genes per replicate
#' @param fun pcalg functions, lingam or ges
#' @param args args to be passed to the main function of pcalg
#' @export
#' @return bn.strength 
pcalg.boot.future <- function(data, R=200, m=nrow(data), removeAllZero=FALSE,
    verbose=FALSE, return.all=FALSE, fun="ges", args=NULL) {
	if (fun=="ges") {
		fun <- ges.pcalg
	} else {
		fun <- lingam.pcalg
	}
    nodes = names(data)
    
    perRun <- future_lapply(seq_len(R), function(r) {
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
        return(run)
    })
    if (return.all) {
        return(perRun)
    }
	st <- custom.strength(perRun, nodes)
    return(st)
}



#' @title HurdleScore
#' uses DataScore class to return sum of BIC of 
#' continous and discrete part of hurdle model.
#' Intercept is TRUE by default.
#' Only the `glm` can be used.
#' @noRd
setRefClass("HurdleScore", contains = "Score",
        methods = list(
            initialize = function(data = matrix(1, 1, 1),
                                  targets = list(integer(0)),
                                  target.index = rep(as.integer(1), nrow(data)),
                                  nodes = colnames(data),
                                  ...) {
                
                p <- ncol(data)
                .pardag.class <<- "GaussParDAG"
                pp.dat$intercept <<- TRUE
                
                pp.dat$targets <<- targets
                pp.dat$target.index <<- target.index
                pp.dat$data <<- data
                pp.dat$vertex.count <<- ncol(data)
                
                .nodes <<- nodes
                ## Store list of index vectors of "non-interventions": for each vertex k,
                ## store the indices of the data points for which k has NOT been intervened
                A <- !targets2mat(pp.dat$vertex.count, pp.dat$targets, pp.dat$target.index)
                pp.dat$non.int <<- lapply(seq_len(ncol(A)), function(i) which(A[, i]))
                # apply() cannot be used since we need a list, not a matrix.
                pp.dat$data.count <<- as.integer(colSums(A))
                pp.dat$total.data.count <<- as.integer(nrow(data))
                
                ## Declare scores as not decomposable "by default"
                decomp <<- TRUE
                
                ## No C++ scoring object by default
                c.fcn <<- "none"
                
                ## R function objects
                pp.dat$local.score <<- function(vertex, parents) local.score(vertex, parents)
                pp.dat$global.score <<- function(edges) global.score(vertex, parents)
                pp.dat$local.fit <<- function(vertex, parents) local.fit(vertex, parents)
                pp.dat$global.fit <<- function(edges) global.fit(vertex, parents)
            },

            getNodes = function() {
                .nodes
            },
            

            node.count = function() {
                length(.nodes)
            },
            

            validate.vertex = function(vertex) {
                if (length(vertex) > 0) {
                    stopifnot(all(is.whole(vertex)))
                    min.max <- range(vertex)
                    stopifnot(1 <= min.max[1] && min.max[2] <= node.count())
                }
            },
            

            validate.parents = function(parents) {
                validate.vertex(parents)
                stopifnot(anyDuplicated(parents) == 0L)
            },
            

            getTargets = function() {
                pp.dat$targets
            },
            
            setTargets = function(targets) {
                pp.dat$targets <<- lapply(targets, sort)
            },
            

            c.fcn.options = function(DEBUG.LEVEL = 0) {
                list(DEBUG.LEVEL = DEBUG.LEVEL)
            },
            local.score = function(vertex, parents, ...) {
                
                ## Check validity of arguments
                validate.vertex(vertex)
                validate.parents(parents)
                yn <- colnames(pp.dat$data)[vertex]
                
                # cat(yn,"\n")
                
                if (length(parents) != 0) {
                    xn <- colnames(pp.dat$data)[parents]
                    fm <- as.formula(paste0(yn ,"~",paste0(xn, collapse="+")))
                    dfm <- as.formula(paste0(yn ,">0~",paste0(xn, collapse="+")))

                } else {
                    fm <- as.formula(paste0(yn ,"~ 1"))
                    dfm <- as.formula(paste0(yn ,">0~1"))
                }
                glf <- glm(dfm, data=pp.dat$data, family="binomial")
                lf <- lm(fm, data=pp.dat$data,subset=yn>0)
                sc <- -1 * (BIC(lf)+BIC(glf))
                
                ## Return local score
                return(sc)
                },
            
            local.fit = function(vertex, parents, ...) {
                validate.vertex(vertex)
                validate.parents(parents)
                yn <- colnames(pp.dat$data)[vertex]
                if (length(parents) != 0) {
                    xn <- colnames(pp.dat$data)[parents]
                    fm <- as.formula(paste0(yn ,"~",paste0(xn, collapse="+")))
                    dfm <- as.formula(paste0(yn ,">0~",paste0(xn, collapse="+")))
                    glf <- glm(dfm, data=pp.dat$data, family="binomial")
                    lf <- lm(fm, data=pp.dat$data,subset=yn>0)
                } else {
                    fm <- as.formula(paste0(yn ,"~ 1"))
                    dfm <- as.formula(paste0(yn ,">0~1"))
                    glf <- glm(dfm, data=pp.dat$data, family="binomial")
                    lf <- lm(fm, data=pp.dat$data,subset=yn>0)
                }
                sigma2 <- sum(lf$model[[yn]]^2)
                beta <-lf$coefficients
                
                ## Calculate regression coefficients
                if (length(parents) + pp.dat$intercept != 0) {
                    qrZ <- lf$qr
                    sigma2 <- sigma2 - sum((lf$model[[yn]] %*% qr.Q(qrZ))^2)
                }
                return(c(sigma2/pp.dat$data.count[vertex], beta))
            }, ## local.fit()
            global.fit = function(dag, ...) {
                in.edge <- dag$.in.edges
                lapply(1:pp.dat$vertex.count,
                       function(i) local.fit(i, in.edge[[i]], ...))
            }
        )
)
