#' @title hurdle.boot
#' 
#' @description
#' Bootstrap-based arc strength calculation based on Hurdle model.
#' 
#' @param data data (row: sample, column: gene)
#' @param R replicate number
#' @param m sampling number
#' @param score scoring function after identification of skeleton.
#' @export
#' @return bn.strength 
hurdle.boot <- function(data, R=200, m=nrow(data), score=NULL) {
    nodes = names(data)
    perRun <- list()
    for (r in seq_len(R)) {
        resampling = sample(nrow(data), m, replace = TRUE)
        # generate the r-th bootstrap sample.
        replicate = data[resampling, , drop = FALSE]
        run <- .Hurdle(replicate, score=score)$bn
        perRun[[r]] <- run
    }
	st <- custom.strength(perRun, nodes)
    return(st)
}
