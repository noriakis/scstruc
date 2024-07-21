#' globalStruc
#' @param algorithm inference algorithm of BN
#' If you want to use CCDr algorithm, specify `ccdr` or `ccdr.boot`.
#' If you use CCDr algorithm, the {bn} object is returned.
#' @param cluster_labels (experimental) aggregate to speed up computation
#' @param use_assay use assay, default to logcounts
#' @param input default to NULL
#' @export
globalStruc <- function(spe, candidate_genes, label=NULL, algorithm="mmhc",
	algorithm.args=list(), return_bn=TRUE, return_data=TRUE,
    cluster_label=NULL,verbose=FALSE, use_assay="logcounts",
    barcode_column="row", change_symbol=TRUE, symbol_column="Symbol", input=NULL,
    nonzero=1) {
    ## .getInput for global label
    x <- NULL
    if (is.null(input)) {
        input <- .getInput(spe, candidate_genes, label, x, use_assay, barcode_column, cluster_label, verbose,
            change_symbol, symbol_column, nonzero)
    }
    global_graph <- .getStruc(input, algorithm, algorithm.args, verbose)

    if (algorithm %in% c("ccdr","ccdr.boot","Hurdle")) {
        if (return_data) {
            return(list(global_graph, input))
        } else {
            return(global_graph)
        }
    }
    if (return_bn) {
        if (return_data) {
            return(list(global_graph, input))
        } else {
            return(global_graph)
        }
    }
    ## When not returning bn object, return tbl_graph
    global_tbl_graph <- global_graph |> 
        bnlearn::as.igraph() |>
        tidygraph::as_tbl_graph()
    return(global_tbl_graph)
}