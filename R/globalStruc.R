#' globalStruc
#' @param algorithm inference algorithm of BN, ignored if reg=TRUE
#' @param cluster_labels (experimental) aggregate to speed up computation
#' @param penalty if {reg}=TRUE, the penalty used can be specified, default to {glmnet} (lasso).
#' If you want to use CCDr algorithm, specify `ccdr` or `ccdr.boot`.
#' If you use CCDr algorithm, the {bn} object is returned.
#' @param use_assay use assay
#' @param input default to NULL
#' @export
globalStruc <- function(spe, candidate_genes, label, algorithm="mmhc",
	reg=TRUE, algorithm.args=list(), return_bn=FALSE, return_data=FALSE,
    cluster_label=NULL, penalty="glmnet",verbose=FALSE, use_assay="logcounts",
    barcode_column="row", change_symbol=TRUE, symbol_column="Symbol", input=NULL,
    nonzero=1) {
    ## .getInput for global label
    x <- NULL
    if (is.null(input)) {
        input <- .getInput(spe, candidate_genes, label, x, use_assay, barcode_column, cluster_label, verbose,
            change_symbol, symbol_column, nonzero)
    }
    global_graph <- .getStruc(input, algorithm, reg, penalty, algorithm.args, verbose)

    if (penalty %in% c("ccdr.run","ccdr.boot")) {
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
    global_tbl_graph <- global_graph |> 
        bnlearn::as.igraph() |>
        tidygraph::as_tbl_graph()
    return(global_tbl_graph)
}