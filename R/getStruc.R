
#' @noRd
.getStruc <- function(input, algorithm, algorithm.args, verbose) {
    if (verbose) {
        cat("Algorithm:", algorithm, "selected\n")        
    }
    if (algorithm=="ccdr") {
        cat("Returning the bn per lambda from result of ccdr.run\n")
        dat <- sparsebnUtils::sparsebnData(input %>% as.matrix(), type = "continuous")
        algorithm.args[["data"]] <- dat
        ccdr.res <- do.call(ccdrAlgorithm::ccdr.run, algorithm.args)
        bn.res <- lapply(seq_along(ccdr.res), function(x) {
            sparsebnUtils::to_bn(ccdr.res[[x]]$edges)
        })
        names(bn.res) <- lapply(seq_along(ccdr.res), function(x) {
            ccdr.res[[x]]$lambda
        }) %>% unlist()
        return(bn.res)
    } else if (algorithm=="ccdr.boot") {
        cat("Returning the list of bn.strength\n")
        algorithm.args[["data"]] <- input
        return(do.call(scstruc::ccdr.boot, algorithm.args))         
    } else if (algorithm == "Hurdle") {
        algorithm.args[["data"]] <- input
        return(do.call(.Hurdle, algorithm.args))
    } else if (algorithm %in% c("mmhc","rsmax2","h2pc","hc","tabu")) {
        cat("Using default bnlearn algorithm\n")
        algorithm.args[["x"]] <- input
        net <- do.call(algorithm, algorithm.args)
        return(net)
    } else {
        if (grepl("boot", algorithm)) {
            cat("Bootstrapping specified\n")
            algorithm <- strsplit(algorithm, "\\.") %>% sapply("[", 1)
            nodes <- names(input)
            perRun <- list()
            m <- nrow(input) ## m determined
            if (is.null(algorithm.args[["R"]])) {
                cat("R not specified in bootstrapping, default to 100\n")
                algorithm.args[["R"]] <- 100
            }
            R <- algorithm.args[["R"]]
            algorithm.args[["R"]] <- NULL
            for (r in seq_len(R)) {
                if (verbose) {cat(r)}
                resampling = sample(nrow(input), m, replace = TRUE)
                replicate = input[resampling, , drop = FALSE]
                algorithm.args[["data"]] <- replicate
                algorithm.args[["algorithm"]] <- algorithm
                repnet <- do.call(skeleton.reg, algorithm.args)
                perRun[[r]] <- repnet
            }
            if (verbose) {cat("\n")}
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