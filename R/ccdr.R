#' @noRd
ccdr <- function(input, algorithm, algorithm.args, verbose) {
    if (!requireNamespace("ccdrAlgorithm")) {
        stop("Needs installation of ccdrAlgorithm")
    } else {
        requireNamespace("ccdrAlgorithm")
    }
    if (!requireNamespace("sparsebnUtils")) {
        stop("Needs installation of sparsebnUtils")
    } else {
        requireNamespace("sparsebnUtils")
    }
    if (is.null(algorithm.args[["bestScore"]])) {
        pickBest <- TRUE
    } else {
        if (isTRUE(algorithm.args[["bestScore"]])) {
            pickBest <- TRUE
        } else {
            pickBest <- FALSE
        }
    }
    algorithm.args[["bestScore"]] <- NULL
    dat <- sparsebnUtils::sparsebnData(input %>% as.matrix(), type = "continuous")
    algorithm.args[["data"]] <- dat
    if (is.null(algorithm.args[["lambdas.length"]]) & is.null(algorithm.args[["lambdas"]])) {
        cat_subtle("Setting `lambdas.length` to 10\n")
        algorithm.args[["lambdas.length"]] = 10
    }
    ccdr.res <- do.call(ccdrAlgorithm::ccdr.run, algorithm.args)
    bn.res <- lapply(seq_along(ccdr.res), function(x) {
        sparsebnUtils::to_bn(ccdr.res[[x]]$edges)
    })
    names(bn.res) <- lapply(seq_along(ccdr.res), function(x) {
        ccdr.res[[x]]$lambda
    }) %>% unlist()

    if (pickBest) {
        cat_subtle("Returning best scored network only, specify algorithm.args=list(bestScore=FALSE) to return all the network\n")
        bestind <- names(which.max(lapply(bn.res, function(tt) {
            bnlearn::score(tt,  data.frame(input, check.names=FALSE)[names(tt$nodes), ])
        })))
        if (length(bestind)==1) {
            return(bn.res[[bestind]])
        } else {
            return(bn.res[bestind])
        }
    } else {
        cat_subtle("Returning the bn per lambda from result of ccdr.run\n")
        return(bn.res)                
    }
}