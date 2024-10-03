#' @noRd
ccdr <- function(input, algorithm, algorithm.args, verbose) {
    if (is.null(algorithm.args[["bestScore"]])) {
        pickBest <- TRUE
    } else {
        if (!algorithm.args[["bestScore"]]) {
            pickBest <- FALSE
        } else {
            pickBest <- TRUE
        }
    }
    algorithm.args[["bestScore"]] <- NULL
    dat <- sparsebnData(input %>% as.matrix(), type = "continuous")
    algorithm.args[["data"]] <- dat
    if (is.null(algorithm.args[["lambdas.length"]]) & is.null(algorithm.args[["lambdas"]])) {
        cat_subtle("Setting `lambdas.length` to 10\n")
        algorithm.args[["lambdas.length"]] = 10
    }
    ccdr.res <- do.call(ccdr.run, algorithm.args)
    bn.res <- lapply(seq_along(ccdr.res), function(x) {
        to_bn(ccdr.res[[x]]$edges)
    })
    names(bn.res) <- lapply(seq_along(ccdr.res), function(x) {
        ccdr.res[[x]]$lambda
    }) %>% unlist()

    if (pickBest) {
        cat_subtle("Returning best scored network only, specify algorithm.args=list(bestScore=FALSE) to return all the network\n")
        bestind <- names(which.max(lapply(bn.res, function(tt) {
            score(tt,  data.frame(input, check.names=FALSE)[names(tt$nodes), ])
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