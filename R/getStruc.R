#' In `ccdr`, one network with best BIC will be returned (if algorithm.args[["bestScore"]] is NULL)
#' In `ccdr.boot`, we cannot see whether the resulting averaged network is DAG, so returning all the network for lambdas.
#' @noRd
.getStruc <- function(input, algorithm, algorithm.args, verbose, boot=FALSE, R=200, m=NULL) {
    pens <- c("glmnet_CV","glmnet_BIC","MCP_CV","SCAD_CV","L0_CV","L0L1_CV","L0L2_CV","plasso")
    if (verbose) {
        cat_subtle("Algorithm: ", algorithm, " selected\n")        
    }
    if (boot) {
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
            ###
            # Bootstrap is not supported
            ###
            algorithm.args[["data"]] <- input
            net <- do.call(pidc.using.julia, algorithm.args)
        } else {
            ## Constraint-based algos are omitted currently
            if (verbose) {
                cat_subtle("Using default bnlearn algorithm: ", algorithm," \n")
            }
            algorithm.args[["x"]] <- input
            net <- do.call(algorithm, algorithm.args)
            return(net)
        }
    }
}