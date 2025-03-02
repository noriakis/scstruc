#' In `ccdr`, one network with best BIC will be returned (if algorithm.args[["bestScore"]] is NULL)
#' In `ccdr.boot`, we cannot see whether the resulting averaged network is DAG, so returning all the network for lambdas.
#' @noRd
#' @importFrom bnlearn custom.strength boot.strength
.getStruc <- function(input, algorithm, algorithm.args, verbose, boot=FALSE, R=200, m=NULL) {
    pens <- c("glmnet_CV","glmnet_BIC","MCP_CV","SCAD_CV","L0_CV","L0L1_CV","L0L2_CV","plasso")
    if (verbose) {
        cat_subtle("Algorithm: ", algorithm, " selected\n")        
    }
    if (isTRUE(boot)) {
        cat_subtle("Bootstrapping specified\n")
        if (is.null(m)) {
            m = nrow(input)
        }
        algorithm.args[["data"]] <- input
        algorithm.args[["R"]] <- R
        algorithm.args[["m"]] <- m

        if (algorithm=="ccdr") {
            cat_subtle("Returning the list of bn.strength\n")
            return(do.call(scstruc::ccdr.boot, algorithm.args))
        } else if (algorithm=="Hurdle") {
            return(do.call(hurdle.boot, algorithm.args))
        } else if (algorithm %in% c("lingam","ges")) {
            algorithm.args[["fun"]] <- algorithm
            return(do.call(pcalg.boot, algorithm.args))
        } else if (algorithm %in% pens) {
            nodes <- names(input)
            perRun <- list()
            algorithm.args[["m"]] <- NULL
            algorithm.args[["R"]] <- NULL

            for (r in seq_len(R)) {
                if (verbose) {cat_subtle(r)}
                resampling = sample(nrow(input), m, replace = TRUE)
                replicate = input[resampling, , drop = FALSE]
                algorithm.args[["data"]] <- replicate
                algorithm.args[["algorithm"]] <- algorithm
                repnet <- do.call(skeleton.reg, algorithm.args)
                perRun[[r]] <- repnet
            }
            if (verbose) {cat_subtle("\n")}
            return(custom.strength(perRun, nodes))
        } else if (algorithm=="pidc") {
            nodes <- names(input)
            algorithm.args[["m"]] <- NULL
            algorithm.args[["R"]] <- NULL
            if (is.null(algorithm.args[["thresholds"]])) {
                algorithm.args[["thresholds"]] <- seq(0.1, 0.4, 0.1)
            }
            perRun <- list()
            for (r in seq_len(R)) {
                if (verbose) {cat_subtle(r)}
                resampling = sample(nrow(input), m, replace = TRUE)
                replicate = input[resampling, , drop = FALSE]
                algorithm.args[["data"]] <- replicate
                algorithm.args[["bestBIC"]] <- FALSE
                repnet <- do.call(pidc.using.julia, algorithm.args)
                perRun[[r]] <- repnet
            }
            if (verbose) {cat_subtle("\n")} 
            ## Average per threshold
            customs <- lapply(algorithm.args[["thresholds"]], function(tmpth) {
                tmpboot <- lapply(seq_len(R), function(tmpr) {
                    perRun[[tmpr]][[as.character(tmpth)]]
                })
                custom.strength(tmpboot, nodes)
            })
            names(customs) <- as.character(algorithm.args[["thresholds"]])
            return(customs)
            # return(custom.strength(perRun, nodes))
        } else {
            algorithm.args[["algorithm"]] <- algorithm
            return(do.call(boot.strength, algorithm.args))
        }
    } else {
        if (algorithm=="ccdr") {
            return(ccdr(input, algorithm, algorithm.args, verbose))         
        } else if (algorithm == "Hurdle") {
            algorithm.args[["data"]] <- input
            return(do.call(.Hurdle, algorithm.args))
        } else if (algorithm %in% c("lingam","ges")) {
            if (algorithm=="lingam") {
                algorithm.args[["data"]] <- input
                return(do.call(lingam.pcalg, algorithm.args))
            } else {
                algorithm.args[["data"]] <- input
                return(do.call(ges.pcalg, algorithm.args))                
            }
        } else if (algorithm %in% pens) {
            algorithm.args[["verbose"]] <- verbose
            algorithm.args[["data"]] <- input
            # algorithm.args[["restrict"]] <- "mmpc"
            # algorithm.args[["maximize"]] <- "hc"
            algorithm.args[["algorithm"]] <- algorithm
            net <- do.call(skeleton.reg, algorithm.args)
            return(net)    
        } else if (algorithm=="pidc") {
            algorithm.args[["data"]] <- input
            net <- do.call(pidc.using.julia, algorithm.args)                
        } else if (algorithm=="genie3") {
            ###
            # Bootstrap is not supported
            # Returning multiple thresholded networks
            ###
            nets <- list()
        	genie.threshold <- seq(0,1,0.1)
	        genie.input <- input %>% t()
	        gra <- GENIE3::GENIE3(genie.input)
	        nets <- lapply(genie.threshold, function(th) {
	            tmp <- gra
	            tmp[tmp < th] <- 0
	            tmp.net <- igraph::graph_from_adjacency_matrix(tmp, weighted = TRUE)
	            if (igraph::is.dag(tmp.net)) {
	                return(bnlearn::as.bn(tmp.net))
	            } else {
	                return(NA)
	            }
	        })
	        names(nets) <- genie.threshold
    	    if (sum(unlist(lapply(nets, function(net) {is.na(net)}))) == length(genie.threshold)) {
        	    cat_subtle("GENIE3 All threshold results in non-DAG\n")
        	} else {
	            for (nn in names(nets)) {
	                genie.net <- nets[[nn]]
	                if (is(genie.net, "bn")) {
	                    nets[[paste0("GENIE3_",nn)]] <- genie.net
	                }
	            }
	        }
            nets[["original"]] <- gra
	        return(nets)
        } else {
            if (verbose) {
                cat_subtle("Using default bnlearn algorithm: ", algorithm, " \n")
            }
            algorithm.args[["x"]] <- input
            net <- do.call(algorithm, algorithm.args)
            return(net)
        }
    }
}