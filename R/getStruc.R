#' In `ccdr`, one network with best BIC will be returned (if algorithm.args[["bestScore"]] is NULL)
#' In `ccdr.boot`, we cannot see whether the resulting averaged network is DAG, so returning all the network for lambdas.
#' @noRd
.getStruc <- function(input, algorithm, algorithm.args, verbose) {
    if (verbose) {
        cat_subtle("Algorithm: ", algorithm, " selected\n")        
    }
    if (algorithm=="ccdr") {
        return(ccdr(input, algorithm, algorithm.args, verbose))
    } else if (algorithm=="ccdr.boot") {
        cat_subtle("Returning the list of bn.strength\n")
        algorithm.args[["data"]] <- input
        return(do.call(scstruc::ccdr.boot, algorithm.args))         
    } else if (algorithm == "Hurdle") {
        algorithm.args[["data"]] <- input
        return(do.call(.Hurdle, algorithm.args))
    } else if (algorithm == "Hurdle.boot") {
        algorithm.args[["data"]] <- input
        return(do.call(hurdle.boot, algorithm.args))
    ## Constraint-based algos are omitted currently
    } else if (algorithm %in% c("mmhc","rsmax2","h2pc","hc","tabu")) {
        cat_subtle("Using default bnlearn algorithm: ", algorithm," \n")
        algorithm.args[["x"]] <- input
        net <- do.call(algorithm, algorithm.args)
        return(net)
    } else {
        if (grepl("boot", algorithm)) {
            cat_subtle("Bootstrapping specified\n")
            algorithm <- strsplit(algorithm, "\\.") %>% sapply("[", 1)
            nodes <- names(input)
            perRun <- list()
            m <- nrow(input) ## m determined
            if (is.null(algorithm.args[["R"]])) {
                cat_subtle("R not specified in bootstrapping, default to 100\n")
                algorithm.args[["R"]] <- 100
            }
            R <- algorithm.args[["R"]]
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
            net <- custom.strength(perRun, nodes)                
        } else {
            algorithm.args[["verbose"]] <- verbose
            algorithm.args[["data"]] <- input
            # algorithm.args[["restrict"]] <- "mmpc"
            # algorithm.args[["maximize"]] <- "hc"
            algorithm.args[["algorithm"]] <- algorithm
            net <- do.call(skeleton.reg, algorithm.args)                
        }
        return(net)
    }
}